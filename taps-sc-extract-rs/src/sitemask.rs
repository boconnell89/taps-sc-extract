//! Per-position call filters applied after a window is called: positions whose
//! read depth exceeds `--site-depth-limit`, and `--blacklist` BED regions.
//! The pileup is never touched (depth caps, overlaps and BAQ stay the same);
//! a masked position loses its calls in every cell.
//!
//! Depth at a position = reads of kept cells that entered the pileup (after
//! the flag, MAPQ and orphan filters) whose reference span covers it,
//! deletions included. Reads of filtered-out cells are not counted, so the
//! depth matches the cells the auto limit is computed from.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::path::Path;

use anyhow::{Context as _, Result};

use crate::accumulate::CellMaps;
use crate::context::{classify_trinucleotide, Context, TriContext};
use crate::extract::CallStats;

/// Intervals (0-based, half-open) where depth exceeds `limit`, from read spans
/// pushed in start order.
pub struct DepthSweep {
    limit: u32,
    ends: BinaryHeap<Reverse<u64>>,
    cur: u64,
    pub masked: Vec<(u64, u64)>,
}

impl DepthSweep {
    pub fn new(limit: u32) -> Self {
        Self { limit, ends: BinaryHeap::new(), cur: 0, masked: Vec::new() }
    }

    fn mark(&mut self, lo: u64, hi: u64) {
        if lo >= hi || (self.ends.len() as u64) <= u64::from(self.limit) {
            return;
        }
        match self.masked.last_mut() {
            Some(last) if last.1 >= lo => last.1 = last.1.max(hi),
            _ => self.masked.push((lo, hi)),
        }
    }

    /// Depth is constant between events; close the stretch up to `to`.
    fn advance(&mut self, to: u64) {
        while let Some(&Reverse(e)) = self.ends.peek() {
            if e > to {
                break;
            }
            self.mark(self.cur, e);
            self.cur = self.cur.max(e);
            self.ends.pop();
        }
        self.mark(self.cur, to);
        self.cur = self.cur.max(to);
    }

    /// One read's reference span `[start, end)`; starts must not decrease.
    pub fn push(&mut self, start: u64, end: u64) {
        if end <= start {
            return;
        }
        self.advance(start);
        self.ends.push(Reverse(end));
    }

    pub fn finish(mut self) -> Vec<(u64, u64)> {
        let last = self.ends.iter().map(|r| r.0).max().unwrap_or(self.cur);
        self.advance(last);
        self.masked
    }
}

/// `--blacklist`: BED regions per contig (0-based, half-open), sorted and merged.
#[derive(Debug, Default)]
pub struct Blacklist {
    by_contig: HashMap<String, Vec<(u64, u64)>>,
}

impl Blacklist {
    /// BED: contig, start, end (further columns ignored); `#`, `track` and
    /// `browser` lines skipped. Gzipped files are read too.
    pub fn load(path: &Path) -> Result<Self> {
        let mut by_contig: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let text = crate::barcode::open_text(path)?;
        use std::io::BufRead as _;
        for (i, line) in text.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() < 3 {
                anyhow::bail!("{} line {}: expected contig, start, end", path.display(), i + 1);
            }
            let s: u64 = f[1].parse().with_context(|| format!("{} line {}: start", path.display(), i + 1))?;
            let e: u64 = f[2].parse().with_context(|| format!("{} line {}: end", path.display(), i + 1))?;
            if e > s {
                by_contig.entry(f[0].to_string()).or_default().push((s, e));
            }
        }
        for v in by_contig.values_mut() {
            *v = merge(std::mem::take(v));
        }
        Ok(Self { by_contig })
    }

    pub fn n_regions(&self) -> usize {
        self.by_contig.values().map(|v| v.len()).sum()
    }

    pub fn total_bp(&self) -> u64 {
        self.by_contig.values().flatten().map(|(s, e)| e - s).sum()
    }

    /// Regions overlapping `[start, end)` of `contig`.
    pub fn overlapping(&self, contig: &str, start: u64, end: u64) -> Vec<(u64, u64)> {
        let Some(v) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        let i = v.partition_point(|r| r.1 <= start);
        v[i..].iter().take_while(|r| r.0 < end).copied().collect()
    }
}

/// Sort and merge overlapping or touching intervals.
pub fn merge(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    v.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(v.len());
    for (s, e) in v {
        match out.last_mut() {
            Some(last) if last.1 >= s => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// Calls removed from one window.
#[derive(Clone, Copy, Debug, Default)]
pub struct Removed {
    pub calls: u64,
    pub sites: u64,
}

/// Remove every call at a position in `masked` (0-based, merged) from `cells`
/// (1-based positions) and from `stats`. `ref_seq` starts at `ref_start` and
/// classifies the removed sites' trinucleotide context for the stats.
pub fn apply(
    cells: &mut [CellMaps],
    stats: &mut CallStats,
    masked: &[(u64, u64)],
    ref_seq: &[u8],
    ref_start: u64,
) -> Removed {
    let mut removed = Removed::default();
    if masked.is_empty() {
        return removed;
    }
    for cell in cells.iter_mut() {
        for (ctx, v) in [(Context::Cg, &mut cell.cg), (Context::Ch, &mut cell.ch)] {
            let mut k = 0usize;
            v.retain(|&(p1, t, c)| {
                let p = u64::from(p1) - 1;
                while k < masked.len() && masked[k].1 <= p {
                    k += 1;
                }
                let hit = k < masked.len() && masked[k].0 <= p;
                if hit {
                    removed.calls += u64::from(t + c);
                    removed.sites += 1;
                    let tri = p
                        .checked_sub(ref_start)
                        .and_then(|rel| classify_trinucleotide(ref_seq, rel as usize))
                        .unwrap_or(match ctx {
                            Context::Cg => TriContext::Cpg,
                            Context::Ch => TriContext::Chh,
                        });
                    stats.remove(tri, ctx, u64::from(t), u64::from(c));
                }
                !hit
            });
        }
    }
    removed
}

/// Observations one cell keeps at a site under `cap` (both counts scaled to
/// `cap`, keeping their ratio; methylated rounded half up). Same rule as
/// `sitemask.cap_counts` in Python.
pub fn cap_counts(t: u32, c: u32, cap: u32) -> (u32, u32) {
    let n = u64::from(t) + u64::from(c);
    if cap == 0 || n <= u64::from(cap) {
        return (t, c);
    }
    let cap64 = u64::from(cap);
    let c2 = ((2 * cap64 * u64::from(c) + n) / (2 * n)) as u32;
    (cap - c2, c2)
}

/// `--cell-site-cap`: no cell keeps more than `cap` observations at a site
/// (`cap_counts`). Calls taken off come out of `stats` as well.
pub fn cap_cells(cells: &mut [CellMaps], stats: &mut CallStats, cap: u32, ref_seq: &[u8], ref_start: u64) -> Removed {
    let mut removed = Removed::default();
    if cap == 0 {
        return removed;
    }
    for cell in cells.iter_mut() {
        for (ctx, v) in [(Context::Cg, &mut cell.cg), (Context::Ch, &mut cell.ch)] {
            for e in v.iter_mut() {
                let (t, c) = cap_counts(e.1, e.2, cap);
                if (t, c) == (e.1, e.2) {
                    continue;
                }
                removed.calls += u64::from(e.1 + e.2 - t - c);
                removed.sites += 1;
                let tri = (u64::from(e.0) - 1)
                    .checked_sub(ref_start)
                    .and_then(|rel| classify_trinucleotide(ref_seq, rel as usize))
                    .unwrap_or(match ctx {
                        Context::Cg => TriContext::Cpg,
                        Context::Ch => TriContext::Chh,
                    });
                stats.remove(tri, ctx, u64::from(e.1 - t), u64::from(e.2 - c));
                e.1 = t;
                e.2 = c;
            }
        }
    }
    removed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sweep_marks_only_where_depth_exceeds_the_limit() {
        // Depth over [10, 20): 3 reads; limit 2 -> masked exactly there.
        let mut s = DepthSweep::new(2);
        s.push(0, 20);
        s.push(5, 25);
        s.push(10, 30);
        assert_eq!(s.finish(), vec![(10, 20)]);
        // Never over the limit.
        let mut s = DepthSweep::new(3);
        s.push(0, 20);
        s.push(5, 25);
        s.push(10, 30);
        assert!(s.finish().is_empty());
        // Two separate deep stretches, a touching end merges.
        let mut s = DepthSweep::new(1);
        s.push(0, 10);
        s.push(5, 10);
        s.push(10, 12);
        s.push(10, 15);
        s.push(100, 110);
        s.push(105, 120);
        assert_eq!(s.finish(), vec![(5, 12), (105, 110)]);
    }

    #[test]
    fn sweep_matches_brute_force() {
        let mut seed = 7u64;
        let mut rnd = |n: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % n
        };
        for limit in [0u32, 1, 2, 5] {
            let mut reads: Vec<(u64, u64)> = (0..300).map(|_| {
                let s = rnd(1000);
                (s, s + 1 + rnd(60))
            }).collect();
            reads.sort();
            let mut depth = vec![0u32; 1100];
            let mut sw = DepthSweep::new(limit);
            for &(s, e) in &reads {
                sw.push(s, e);
                for d in &mut depth[s as usize..e as usize] {
                    *d += 1;
                }
            }
            let got = sw.finish();
            for (p, &d) in depth.iter().enumerate() {
                let inside = got.iter().any(|&(s, e)| s <= p as u64 && (p as u64) < e);
                assert_eq!(inside, d > limit, "limit {limit} pos {p} depth {d}");
            }
        }
    }

    #[test]
    fn apply_removes_masked_calls_and_stats() {
        let mut cell = CellMaps::default();
        cell.cg = vec![(11, 1, 0), (21, 0, 2), (31, 1, 1)];
        cell.ch = vec![(15, 3, 0)];
        let mut cells = vec![cell];
        let mut stats = CallStats::default();
        stats.cg_t = 2;
        stats.cg_c = 3;
        stats.cpg_t = 2;
        stats.cpg_c = 3;
        stats.ch_t = 3;
        stats.chh_t = 3;
        // 0-based [20, 21) masks 1-based 21; [14, 15) masks 15.
        let ref_seq = b"NNNNNNNNNNACNNTANNNNACGNNNNNNNACGTNNNNNN";
        let r = apply(&mut cells, &mut stats, &[(14, 15), (20, 21)], ref_seq, 0);
        assert_eq!((r.calls, r.sites), (5, 2));
        assert_eq!(cells[0].cg, vec![(11, 1, 0), (31, 1, 1)]);
        assert!(cells[0].ch.is_empty());
        assert_eq!((stats.cg_t, stats.cg_c, stats.ch_t), (2, 1, 0));
    }

    #[test]
    fn cap_keeps_ratio_and_total() {
        assert_eq!(cap_counts(1, 1, 2), (1, 1));
        assert_eq!(cap_counts(0, 3, 2), (0, 2));
        assert_eq!(cap_counts(3, 0, 2), (2, 0));
        assert_eq!(cap_counts(1, 2, 2), (1, 1)); // 2 x 2/3 = 1.33 -> 1
        assert_eq!(cap_counts(1, 3, 2), (0, 2)); // 1.5 rounds up
        assert_eq!(cap_counts(10, 30, 4), (1, 3));
        assert_eq!(cap_counts(5, 5, 0), (5, 5));
        let mut cell = CellMaps::default();
        cell.cg = vec![(11, 1, 2), (21, 0, 1)];
        let mut cells = vec![cell];
        let mut stats = CallStats { cg_t: 1, cg_c: 3, cpg_t: 1, cpg_c: 3, ..Default::default() };
        let r = cap_cells(&mut cells, &mut stats, 2, b"NNNNNNNNNNCGNNNNNNNNCG", 0);
        assert_eq!((r.calls, r.sites), (1, 1));
        assert_eq!(cells[0].cg, vec![(11, 1, 1), (21, 0, 1)]);
        assert_eq!((stats.cg_t, stats.cg_c, stats.cpg_c), (1, 2, 2));
    }

    #[test]
    fn blacklist_overlap_and_merge() {
        let dir = std::env::temp_dir().join(format!("taps_bl_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("bl.bed");
        std::fs::write(&p, "track name=x\n#c\nchr2\t100\t200\tsat\nchr2\t150\t250\nchr2\t400\t500\nchr9\t0\t10\n").unwrap();
        let b = Blacklist::load(&p).unwrap();
        assert_eq!(b.n_regions(), 3);
        assert_eq!(b.total_bp(), 150 + 100 + 10);
        assert_eq!(b.overlapping("chr2", 0, 1000), vec![(100, 250), (400, 500)]);
        assert_eq!(b.overlapping("chr2", 250, 400), Vec::<(u64, u64)>::new());
        assert!(b.overlapping("chr1", 0, 1000).is_empty());
        std::fs::remove_dir_all(&dir).ok();
    }
}
