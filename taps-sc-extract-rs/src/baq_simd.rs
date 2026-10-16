//! Batched SIMD BAQ (`csrc/baq_simd.c`).
//!
//! `taps_baq_batch` gives each read exactly what htslib's
//! `sam_prob_realn(b, ref, len, BAQ_APPLY|BAQ_EXTEND)` gives it, but runs the
//! profile HMM on 4 reads at once (AVX2) or 8 (AVX-512F). Reads whose HMM
//! matrices have the same shape share a vector; the rest of a batch is padded.
//! Both kernels are always compiled; the CPU is checked at run time.
//! `selftest_*` compare it with htslib read by read (`baq-selftest`, run by
//! `install.sh`).

use anyhow::Result;
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{Read, Record};
use rust_htslib::htslib;
use std::os::raw::{c_char, c_int, c_void};
use std::path::Path;

#[repr(C)]
pub struct BaqJob {
    pub b: *mut htslib::bam1_t,
    pub reference: *const c_char,
    pub ref_len: i64,
}

extern "C" {
    fn taps_baq_cpu_width() -> c_int;
    fn taps_baq_width_ok(width: c_int) -> c_int;
    fn taps_baq_ws_new() -> *mut c_void;
    fn taps_baq_ws_free(ws: *mut c_void);
    fn taps_baq_batch(ws: *mut c_void, jobs: *const BaqJob, n: c_int, width: c_int) -> c_int;
}

/// `--baq-simd`: which BAQ kernel to use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum BaqSimd {
    /// Widest kernel this CPU supports.
    Auto,
    /// 8 reads per vector (AVX-512F).
    Avx512,
    /// 4 reads per vector (AVX2).
    Avx2,
    /// htslib's scalar `sam_prob_realn`, one read at a time.
    Off,
}

/// Kernel the C side runs (its width code).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Width(pub i32);

impl Width {
    pub const HTSLIB: Width = Width(0);
    pub const AVX2: Width = Width(4);
    pub const AVX512: Width = Width(8);
    /// 8 lanes on AVX2 instructions: checks the 8-lane logic on AVX2-only hosts.
    pub const AVX2_X8: Width = Width(108);

    pub fn label(self) -> &'static str {
        match self.0 {
            4 => "avx2 (4 reads/vector)",
            8 => "avx512 (8 reads/vector)",
            108 => "avx2x8 (8 lanes on AVX2, test)",
            _ => "htslib (scalar)",
        }
    }

    pub fn supported(self) -> bool {
        unsafe { taps_baq_width_ok(self.0) != 0 }
    }

    /// Widest kernel this CPU runs.
    pub fn cpu_best() -> Width {
        Width(unsafe { taps_baq_cpu_width() })
    }

    pub fn resolve(choice: BaqSimd) -> Result<Width> {
        let w = match choice {
            BaqSimd::Auto => return Ok(Width::cpu_best()),
            BaqSimd::Avx512 => Width::AVX512,
            BaqSimd::Avx2 => Width::AVX2,
            BaqSimd::Off => Width::HTSLIB,
        };
        anyhow::ensure!(w.supported(), "--baq-simd {choice:?}: this CPU does not support it");
        Ok(w)
    }
}

/// Per-thread scratch for `taps_baq_batch` (grows to the largest batch seen).
pub struct Workspace(*mut c_void);

// The workspace is plain heap memory owned by one thread at a time.
unsafe impl Send for Workspace {}

impl Workspace {
    pub fn new() -> Self {
        let p = unsafe { taps_baq_ws_new() };
        assert!(!p.is_null(), "taps_baq_ws_new: out of memory");
        Self(p)
    }

    /// BAQ every job's read in place. Returns how many went through SIMD.
    pub fn run(&mut self, jobs: &[BaqJob], width: Width) -> usize {
        if jobs.is_empty() {
            return 0;
        }
        let n = unsafe { taps_baq_batch(self.0, jobs.as_ptr(), jobs.len() as c_int, width.0) };
        n.max(0) as usize
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        unsafe { taps_baq_ws_free(self.0) }
    }
}

/// Reads compared and reads whose qualities differ.
#[derive(Debug, Default, Clone, Copy)]
pub struct SelftestResult {
    pub reads: usize,
    pub simd_reads: usize,
    pub mismatched: usize,
}

/// BAQ `recs` (with references `refs[i]`) through htslib and through `width`,
/// on copies, and count reads whose QUAL differs.
fn compare(recs: &[Record], refs: &[&[u8]], width: Width, ws: &mut Workspace) -> SelftestResult {
    let mut a: Vec<Record> = recs.to_vec();
    let mut b: Vec<Record> = recs.to_vec();
    for (r, s) in a.iter_mut().zip(refs) {
        unsafe {
            htslib::sam_prob_realn(r.inner_mut(), s.as_ptr() as *const c_char, s.len() as i64, 3);
        }
    }
    let jobs: Vec<BaqJob> = b
        .iter_mut()
        .zip(refs)
        .map(|(r, s)| BaqJob { b: r.inner_mut(), reference: s.as_ptr() as *const c_char, ref_len: s.len() as i64 })
        .collect();
    let simd_reads = ws.run(&jobs, width);
    let mismatched = a.iter().zip(&b).filter(|(x, y)| x.qual() != y.qual()).count();
    SelftestResult { reads: recs.len(), simd_reads, mismatched }
}

/// Small xorshift generator (deterministic, no dependency).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

/// Random reads against a random reference: varied lengths, mismatches, N,
/// clips, insertions and deletions, low and high qualities, and reads running
/// off the reference end. Many reads share lengths so lanes fill.
pub fn selftest_synthetic(width: Width, n_reads: usize, seed: u64) -> Result<SelftestResult> {
    let mut rng = Rng(seed | 1);
    let ref_len = 50_000usize;
    let reference: Vec<u8> = (0..ref_len)
        .map(|_| match rng.below(100) {
            0 => b'N',
            x => b"ACGT"[(x % 4) as usize],
        })
        .collect();
    let mut recs = Vec::with_capacity(n_reads);
    for i in 0..n_reads {
        let len = match rng.below(10) {
            0 => 20 + rng.below(131) as usize,
            1..=3 => 99,
            4..=6 => 89,
            _ => 80 + rng.below(20) as usize,
        };
        let mut ops: Vec<Cigar> = Vec::new();
        let mut used = 0u32;
        let left = len as u32;
        if rng.below(10) == 0 {
            let s = 1 + rng.below(8) as u32;
            ops.push(Cigar::SoftClip(s));
            used += s;
        }
        while used < left {
            let rest = left - used;
            let m = (1 + rng.below(60) as u32).min(rest);
            ops.push(Cigar::Match(m));
            used += m;
            if used < left && rng.below(4) == 0 {
                match rng.below(3) {
                    0 => {
                        let l = (1 + rng.below(4) as u32).min(left - used);
                        ops.push(Cigar::Ins(l));
                        used += l;
                    }
                    1 => ops.push(Cigar::Del(1 + rng.below(12) as u32)),
                    _ => {
                        let l = left - used;
                        ops.push(Cigar::SoftClip(l));
                        used += l;
                    }
                }
            }
        }
        // Merge adjacent soft clips / ensure the last op is not a lone deletion.
        let cigar = CigarString(ops);
        let span = cigar.iter().map(|c| match c {
            Cigar::Match(l) | Cigar::Del(l) => *l as usize,
            _ => 0,
        }).sum::<usize>();
        let pos = if rng.below(50) == 0 {
            ref_len.saturating_sub(span / 2)
        } else {
            rng.below((ref_len - span - 1) as u64) as usize
        };
        let mut seq = Vec::with_capacity(len);
        let mut rpos = pos;
        for c in cigar.iter() {
            match c {
                Cigar::Match(l) => {
                    for _ in 0..*l {
                        let base = *reference.get(rpos).unwrap_or(&b'A');
                        let b = match rng.below(40) {
                            0 => b"ACGT"[rng.below(4) as usize],
                            1 if base == b'C' => b'T',
                            2 => b'N',
                            _ => if base == b'N' { b'A' } else { base },
                        };
                        seq.push(b);
                        rpos += 1;
                    }
                }
                Cigar::Ins(l) | Cigar::SoftClip(l) => {
                    for _ in 0..*l {
                        seq.push(b"ACGT"[rng.below(4) as usize]);
                    }
                }
                Cigar::Del(l) => rpos += *l as usize,
                _ => {}
            }
        }
        let qual: Vec<u8> = (0..len)
            .map(|_| match rng.below(10) {
                0 => rng.below(10) as u8,
                _ => 20 + rng.below(22) as u8,
            })
            .collect();
        let mut r = Record::new();
        r.set(format!("r{i}").as_bytes(), Some(&cigar), &seq, &qual);
        r.set_pos(pos as i64);
        r.set_tid(0);
        r.set_flags(if rng.below(2) == 0 { 99 } else { 83 });
        r.set_mapq(60);
        recs.push(r);
    }
    let refs: Vec<&[u8]> = vec![reference.as_slice(); recs.len()];
    let mut ws = Workspace::new();
    let mut total = SelftestResult::default();
    for (chunk, rchunk) in recs.chunks(4096).zip(refs.chunks(4096)) {
        let r = compare(chunk, rchunk, width, &mut ws);
        total.reads += r.reads;
        total.simd_reads += r.simd_reads;
        total.mismatched += r.mismatched;
    }
    Ok(total)
}

/// Every mapped primary read of `region` (contig name), up to `max_reads`,
/// against the TAPS-masked reference for its strand.
pub fn selftest_bam(bam_path: &Path, fasta: &Path, contig: &str, max_reads: usize, width: Width) -> Result<SelftestResult> {
    let fa = crate::fasta::FastFaiReader::open(fasta)?;
    let len = fa
        .reference_length(contig)
        .ok_or_else(|| anyhow::anyhow!("contig {contig} not in FASTA"))?;
    let baq_ref = crate::extract::BaqRef::new(fa.fetch(contig, 0, len)?, crate::extract::BaqMode::TapsFull);
    let mut bam = crate::input::open_indexed(bam_path)?;
    bam.fetch(contig).map_err(|e| anyhow::anyhow!("fetch {contig}: {e}"))?;
    let mut ws = Workspace::new();
    let mut total = SelftestResult::default();
    let mut batch: Vec<Record> = Vec::with_capacity(4096);
    let mut rec = Record::new();
    let mut flush = |batch: &mut Vec<Record>, total: &mut SelftestResult| {
        let refs: Vec<&[u8]> = batch
            .iter()
            .map(|r| if crate::extract::baq_is_ob(r.flags()) { &*baq_ref.ob } else { &*baq_ref.ot })
            .collect();
        let r = compare(batch, &refs, width, &mut ws);
        total.reads += r.reads;
        total.simd_reads += r.simd_reads;
        total.mismatched += r.mismatched;
        batch.clear();
    };
    while total.reads + batch.len() < max_reads {
        match bam.read(&mut rec) {
            None => break,
            Some(Err(e)) => anyhow::bail!("BAM read: {e}"),
            Some(Ok(())) => {
                if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
                    continue;
                }
                batch.push(rec.clone());
                if batch.len() == 4096 {
                    flush(&mut batch, &mut total);
                }
            }
        }
    }
    if !batch.is_empty() {
        flush(&mut batch, &mut total);
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simd_baq_matches_htslib_on_synthetic_reads() {
        let mut widths = vec![Width::HTSLIB];
        for w in [Width::AVX2, Width::AVX2_X8, Width::AVX512] {
            if w.supported() {
                widths.push(w);
            }
        }
        for w in widths {
            let r = selftest_synthetic(w, 20_000, 7).unwrap();
            assert_eq!(r.mismatched, 0, "{}: {r:?}", w.label());
            if w != Width::HTSLIB {
                assert!(r.simd_reads > r.reads / 2, "{}: {r:?}", w.label());
            }
        }
    }
}
