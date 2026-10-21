//! Read-centric calling: the same calls as htslib's `bam_plp` column pileup,
//! computed one read at a time.
//!
//! A read's calls depend on other reads in only two ways, both decided when
//! the read is pushed (htslib `bam_plp_push`, sam.c):
//! - **depth cap** (`--max-depth`, `bam_plp_set_maxcnt`): a read is dropped
//!   when it starts at the pileup's current position (the start of the last
//!   accepted read, since every push is drained) and the iterator holds more
//!   than `maxcnt` nodes: the live reads plus one sentinel.
//! - **mate overlap** (`bam_plp_init_overlaps`): `overlap_push` pairs mates
//!   by QNAME and `tweak_overlap_quality` (htslib 1.21, `csrc/olap.c`)
//!   lowers qualities in the overlap of both mates when the second arrives.
//!   A read leaves the QNAME table when htslib frees it (the pileup has
//!   passed its end) or when a read with its QNAME is dropped by the cap.
//!
//! `ReadPileup` replays exactly that bookkeeping (live reads by end position,
//! the QNAME table, the cap) without building columns. A read is called when
//! htslib would free it: by then its mate has either been pushed (and both
//! were tweaked) or can no longer pair with it. The tweak only changes bases
//! at or after the later mate's start, which htslib had not yet emitted when
//! the later mate was pushed, so calling the whole read at that point gives
//! the calls the column pileup gives.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use rust_htslib::bam::record::Record;
use rust_htslib::htslib;
use rustc_hash::FxHashMap;

use crate::accumulate::{ensure_cell, BarcodeIntern, CellMaps};
use crate::calling::{call_mctot, classify_strand, Strand};
use crate::context::{classify_trinucleotide, Context, TriContext};

extern "C" {
    fn taps_tweak_overlap_quality(a: *mut htslib::bam1_t, b: *mut htslib::bam1_t) -> std::os::raw::c_int;
}

/// Reference length of the alignment (`bam_cigar2rlen`: M, D, N, =, X).
#[inline]
fn cigar2rlen(rec: &Record) -> i64 {
    rec.raw_cigar()
        .iter()
        .filter(|&&c| matches!(c & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|&c| i64::from(c >> 4))
        .sum()
}

/// htslib's pileup bookkeeping for one window (fresh per window, like the
/// `bam_plp_t` it replaces).
pub struct ReadPileup {
    /// Copies of live reads (`bam_copy1`, as htslib keeps them); slots reused.
    slab: Vec<Record>,
    free: Vec<u32>,
    /// Live (linked) reads by exclusive end position.
    live: BinaryHeap<Reverse<(i64, u32)>>,
    /// QNAME -> slot of a read waiting for its mate (htslib `olap_hash`).
    olap: FxHashMap<Box<[u8]>, u32>,
    overlaps: bool,
    max_depth: usize,
    /// htslib `iter->pos` / `iter->tid` after the last drain.
    pos: i64,
    first: bool,
    /// Reads dropped by the depth cap.
    pub dropped: u64,
}

impl ReadPileup {
    pub fn new(max_depth: i32, overlaps: bool) -> Self {
        ReadPileup {
            slab: Vec::new(),
            free: Vec::new(),
            live: BinaryHeap::new(),
            olap: FxHashMap::default(),
            overlaps,
            max_depth: max_depth.max(0) as usize,
            pos: 0,
            first: true,
            dropped: 0,
        }
    }

    /// `bam_plp_push(rec)` followed by draining every ready column. Reads
    /// htslib frees during that drain are passed to `call`.
    pub fn push(&mut self, rec: &Record, call: &mut impl FnMut(&Record)) -> anyhow::Result<()> {
        // Reads freed by the drain after the previous push (columns < pos).
        self.retire(self.pos, call);
        let pos = rec.pos();
        let tid = rec.tid();
        // iter->tid starts at 0 and follows the reads after the first drain.
        let same_tid = !self.first || tid == 0;
        if same_tid && pos == self.pos && self.live.len() + 1 > self.max_depth {
            self.dropped += 1;
            if self.overlaps {
                self.olap.remove(rec.qname());
            }
            return Ok(());
        }
        let end = pos + cigar2rlen(rec);
        let linked = end > self.pos || (self.first && tid > 0);
        // The drain moves the pileup to this read's start.
        self.pos = pos;
        self.first = false;
        if !linked {
            // No reference span: htslib reuses its node and never calls it.
            return Ok(());
        }
        let slot = match self.free.pop() {
            Some(s) => s,
            None => {
                self.slab.push(Record::new());
                (self.slab.len() - 1) as u32
            }
        };
        if unsafe { htslib::bam_copy1(self.slab[slot as usize].inner_mut(), rec.inner()) }.is_null() {
            anyhow::bail!("bam_copy1 failed");
        }
        if self.overlaps {
            self.overlap_push(slot, end)?;
        }
        self.live.push(Reverse((end, slot)));
        Ok(())
    }

    /// End of input: every remaining read is called.
    pub fn finish(&mut self, call: &mut impl FnMut(&Record)) {
        self.retire(i64::MAX, call);
        self.olap.clear();
    }

    /// Free (and call) reads that end before `pos` (htslib: `end <= pos - 1`
    /// once column `pos - 1` has been emitted).
    fn retire(&mut self, pos: i64, call: &mut impl FnMut(&Record)) {
        while let Some(&Reverse((end, slot))) = self.live.peek() {
            if end >= pos {
                break;
            }
            self.live.pop();
            let rec = &self.slab[slot as usize];
            call(rec);
            if self.overlaps {
                // overlap_remove: by QNAME, whichever read the entry holds.
                self.olap.remove(rec.qname());
            }
            self.free.push(slot);
        }
    }

    /// htslib `overlap_push`.
    fn overlap_push(&mut self, slot: u32, end: i64) -> anyhow::Result<()> {
        let b = &self.slab[slot as usize];
        let flags = b.flags() as u32;
        if flags & htslib::BAM_FMUNMAP != 0 || flags & htslib::BAM_FPROPER_PAIR == 0 {
            return Ok(());
        }
        let core = &b.inner().core;
        if (core.mtid >= 0 && core.tid != core.mtid)
            || (core.isize_.unsigned_abs() as i64 >= 2 * i64::from(core.l_qseq) && core.mpos >= end)
        {
            return Ok(());
        }
        match self.olap.get(b.qname()) {
            None => {
                if core.mpos >= core.pos || (flags & htslib::BAM_FPAIRED != 0 && core.mpos == -1) {
                    let key: Box<[u8]> = b.qname().into();
                    self.olap.insert(key, slot);
                }
            }
            Some(&a) => {
                let pa = self.slab[a as usize].inner_mut() as *mut htslib::bam1_t;
                let pb = self.slab[slot as usize].inner_mut() as *mut htslib::bam1_t;
                let err = unsafe { taps_tweak_overlap_quality(pa, pb) };
                let qname: Box<[u8]> = self.slab[slot as usize].qname().into();
                self.olap.remove(&qname);
                if err < 0 {
                    anyhow::bail!("overlap quality adjustment failed (bad CIGAR?)");
                }
            }
        }
        Ok(())
    }
}

/// Window-wide inputs for `call_read`.
pub struct CallWindow<'a> {
    pub ref_seq: &'a [u8],
    pub ref_start: u64,
    pub win_start: u32,
    pub win_end: u32,
    pub params: &'a crate::extract::ExtractParams,
    pub whitelist: Option<&'a std::collections::HashSet<String>>,
}

/// Per-window outputs and scratch for `call_read`.
#[derive(Default)]
pub struct CallSink {
    pub stats: crate::extract::CallStats,
    pub intern: BarcodeIntern,
    pub cells: Vec<CellMaps>,
    cg: Vec<(u32, u8)>,
    ch: Vec<(u32, u8)>,
}

/// Every call one read makes inside the window: the per-base rules of the
/// column callback (`extract::process_window`), applied along the CIGAR.
pub fn call_read(w: &CallWindow, rec: &Record, out: &mut CallSink) {
    let flags = rec.flags();
    let Some(strand) = classify_strand(flags) else {
        return;
    };
    // OT reads call at reference C, OB reads at reference G.
    let want = match strand {
        Strand::Ot => b'C',
        Strand::Ob => b'G',
    };
    let p = w.params;
    let seq = rec.seq();
    let seq_len = seq.len();
    let qual = rec.qual();
    // 5' trim (`keep_after_5prime_trim`) as a qpos range.
    let trim = if flags & htslib::BAM_FREAD2 as u16 != 0 { p.trim_r2 } else { p.trim_r1 } as usize;
    let (qlo, qhi) = if trim == 0 {
        (0, seq_len)
    } else if trim >= seq_len {
        (0, 0)
    } else if flags & htslib::BAM_FREVERSE as u16 != 0 {
        (0, seq_len - trim)
    } else {
        (trim, seq_len)
    };
    let (win_lo, win_hi) = (i64::from(w.win_start), i64::from(w.win_end));
    let mut cell: Option<u32> = None;
    out.cg.clear();
    out.ch.clear();
    let mut rpos = rec.pos();
    let mut q = 0usize;
    for &c in rec.raw_cigar() {
        let l = (c >> 4) as usize;
        match c & 0xf {
            // M, =, X: one base per reference position.
            0 | 7 | 8 => {
                let lo = rpos.max(win_lo);
                let hi = (rpos + l as i64).min(win_hi);
                let mut r = lo;
                while r < hi {
                    let qp = q + (r - rpos) as usize;
                    let rel = (r as u64).wrapping_sub(w.ref_start) as usize;
                    r += 1;
                    let Some(&rb) = w.ref_seq.get(rel) else {
                        continue;
                    };
                    if rb != want || qp >= seq_len || qp < qlo || qp >= qhi || qual[qp] < p.min_baseq {
                        continue;
                    }
                    let Some(call) = call_mctot(strand, rb, seq[qp]) else {
                        continue;
                    };
                    let Some(tri) = classify_trinucleotide(w.ref_seq, rel) else {
                        continue;
                    };
                    let ctx = match tri {
                        TriContext::Cpg => Context::Cg,
                        _ => Context::Ch,
                    };
                    if !p.context.wants(ctx) {
                        continue;
                    }
                    if cell.is_none() {
                        let bc = crate::barcode::barcode_from_qname_bytes(rec.qname());
                        let valid = match w.whitelist {
                            Some(wl) => wl.contains(std::str::from_utf8(bc).unwrap_or("")),
                            None => true,
                        };
                        if !valid {
                            return;
                        }
                        cell = Some(out.intern.intern_bytes(bc));
                    }
                    out.stats.add(tri, ctx, call);
                    let site = (r as u32, call); // r is already pos + 1 (1-based)
                    match ctx {
                        Context::Cg => out.cg.push(site),
                        Context::Ch => out.ch.push(site),
                    }
                }
                rpos += l as i64;
                q += l;
            }
            // I, S: query only.
            1 | 4 => q += l,
            // D, N: reference only.
            2 | 3 => rpos += l as i64,
            _ => {}
        }
    }
    if let (Some(id), true) = (cell, p.accumulate) {
        let cm = ensure_cell(&mut out.cells, id);
        cm.add_run(Context::Cg, &out.cg);
        cm.add_run(Context::Ch, &out.ch);
    }
}
