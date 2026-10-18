//! Per-window barcode intern and sparse (pos, t, c) maps.
//!
//! Intern ids are window-local. Merge across windows by barcode string, not id.
//! Positions are **1-based** (Amethyst).

use crate::barcode::barcode_from_qname_bytes;
use crate::context::Context;
use rustc_hash::FxHashMap;

#[derive(Clone, Debug, Default)]
pub struct BarcodeIntern {
    to_id: FxHashMap<Vec<u8>, u32>,
    names: Vec<String>,
}

impl BarcodeIntern {
    pub fn intern_bytes(&mut self, raw: &[u8]) -> u32 {
        if let Some(&id) = self.to_id.get(raw) {
            return id;
        }
        let id = self.names.len() as u32;
        self.to_id.insert(raw.to_vec(), id);
        self.names.push(String::from_utf8_lossy(raw).into_owned());
        id
    }

    pub fn intern_str(&mut self, s: &str) -> u32 {
        self.intern_bytes(s.as_bytes())
    }

    pub fn intern_qname(&mut self, qname: &[u8]) -> u32 {
        self.intern_bytes(barcode_from_qname_bytes(qname))
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn name(&self, id: u32) -> &str {
        &self.names[id as usize]
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }
}

/// Per-cell CG/CH sites: `(1-based pos, t unmethylated, c methylated)`, sorted
/// by position with one entry per position.
///
/// Pileup columns arrive in position order, so a call either lands on the last
/// entry (same position, another read) or appends a new one: no hashing and no
/// sort. A position arriving out of order (not expected from `bam_plp`) marks
/// the cell for a sort-and-merge in `normalize`.
#[derive(Clone, Debug, Default)]
pub struct CellMaps {
    pub cg: Vec<(u32, u32, u32)>,
    pub ch: Vec<(u32, u32, u32)>,
    unsorted: bool,
}

#[inline]
fn add_site(v: &mut Vec<(u32, u32, u32)>, pos: u32, t: u32, c: u32) -> bool {
    match v.last_mut() {
        Some(last) if last.0 == pos => {
            last.1 += t;
            last.2 += c;
            false
        }
        Some(last) if last.0 > pos => {
            v.push((pos, t, c));
            true
        }
        _ => {
            v.push((pos, t, c));
            false
        }
    }
}

fn sort_merge(v: &mut Vec<(u32, u32, u32)>) {
    v.sort_unstable_by_key(|x| x.0);
    let mut out: Vec<(u32, u32, u32)> = Vec::with_capacity(v.len());
    for &(p, t, c) in v.iter() {
        match out.last_mut() {
            Some(last) if last.0 == p => {
                last.1 += t;
                last.2 += c;
            }
            _ => out.push((p, t, c)),
        }
    }
    *v = out;
}

impl CellMaps {
    #[inline]
    pub fn add(&mut self, ctx: Context, pos_1based: u32, meth: u8) {
        let v = match ctx {
            Context::Cg => &mut self.cg,
            Context::Ch => &mut self.ch,
        };
        let (t, c) = if meth == 1 { (0, 1) } else { (1, 0) };
        if add_site(v, pos_1based, t, c) {
            self.unsorted = true;
        }
    }

    /// Add one read's calls in `ctx` (`(1-based pos, meth)`, positions
    /// increasing). Read-centric calling visits a cell's reads in start
    /// order, so a run usually appends; a mate that overlaps the previous
    /// read lands inside the tail, which is merged in place.
    pub fn add_run(&mut self, ctx: Context, run: &[(u32, u8)]) {
        let Some(&(first, _)) = run.first() else {
            return;
        };
        let v = match ctx {
            Context::Cg => &mut self.cg,
            Context::Ch => &mut self.ch,
        };
        let site = |&(p, m): &(u32, u8)| if m == 1 { (p, 0, 1) } else { (p, 1, 0) };
        if v.last().is_none_or(|l| l.0 < first) {
            v.extend(run.iter().map(site));
            return;
        }
        // Entries at or after the run's first position (sorted, unique) plus
        // the run (sorted, unique): sort the two together, then merge equals.
        let s = v.partition_point(|e| e.0 < first);
        v.extend(run.iter().map(site));
        v[s..].sort_by_key(|e| e.0);
        let mut w = s;
        for i in s..v.len() {
            if w > s && v[w - 1].0 == v[i].0 {
                v[w - 1].1 += v[i].1;
                v[w - 1].2 += v[i].2;
            } else {
                v[w] = v[i];
                w += 1;
            }
        }
        v.truncate(w);
    }

    /// Restore the sorted, one-entry-per-position invariant if `add` saw an
    /// out-of-order position. No-op otherwise.
    pub fn normalize(&mut self) {
        if self.unsorted {
            sort_merge(&mut self.cg);
            sort_merge(&mut self.ch);
            self.unsorted = false;
        }
    }

    pub fn merge_from(&mut self, other: &CellMaps) {
        self.cg.extend_from_slice(&other.cg);
        self.ch.extend_from_slice(&other.ch);
        sort_merge(&mut self.cg);
        sort_merge(&mut self.ch);
        self.unsorted = false;
    }

    pub fn is_empty(&self) -> bool {
        self.cg.is_empty() && self.ch.is_empty()
    }
}

pub fn ensure_cell(cells: &mut Vec<CellMaps>, id: u32) -> &mut CellMaps {
    let i = id as usize;
    if cells.len() <= i {
        cells.resize_with(i + 1, CellMaps::default);
    }
    &mut cells[i]
}

pub fn write_sites_tsv(
    path: &std::path::Path,
    intern: &BarcodeIntern,
    cells: &[CellMaps],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut rows: Vec<(&str, &str, u32, u32, u32)> = Vec::new();
    for (i, name) in intern.names().iter().enumerate() {
        if i >= cells.len() {
            continue;
        }
        let cell = &cells[i];
        for &(pos, t, c) in &cell.cg {
            rows.push((name.as_str(), "CG", pos, t, c));
        }
        for &(pos, t, c) in &cell.ch {
            rows.push((name.as_str(), "CH", pos, t, c));
        }
    }
    rows.sort_unstable();
    let mut f = std::fs::File::create(path)?;
    writeln!(f, "barcode\tctx\tpos\tt\tc")?;
    for (bc, ctx, pos, t, c) in rows {
        writeln!(f, "{bc}\t{ctx}\t{pos}\t{t}\t{c}")?;
    }
    Ok(())
}

/// Remap `src` cell maps onto `dst` intern ids (by barcode string).
pub fn merge_window_cells(
    dst_intern: &mut BarcodeIntern,
    dst_cells: &mut Vec<CellMaps>,
    src_intern: &BarcodeIntern,
    src_cells: &[CellMaps],
) {
    for (i, name) in src_intern.names().iter().enumerate() {
        let id = dst_intern.intern_str(name);
        if i >= src_cells.len() {
            continue;
        }
        ensure_cell(dst_cells, id).merge_from(&src_cells[i]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_run_appends_and_merges_overlapping_mates() {
        let mut m = CellMaps::default();
        m.add_run(Context::Ch, &[(10, 1), (12, 0), (15, 1)]);
        // A mate overlapping the tail: equal positions sum, others interleave.
        m.add_run(Context::Ch, &[(11, 0), (12, 1), (16, 0)]);
        m.add_run(Context::Ch, &[(20, 1)]);
        m.add_run(Context::Cg, &[]);
        assert_eq!(m.ch, vec![(10, 0, 1), (11, 1, 0), (12, 1, 1), (15, 0, 1), (16, 1, 0), (20, 0, 1)]);
        assert!(m.cg.is_empty());
    }

    #[test]
    fn intern_reuses_ids() {
        let mut i = BarcodeIntern::default();
        let a = i.intern_qname(b"CELL1:rest");
        let b = i.intern_qname(b"CELL1:other");
        let c = i.intern_qname(b"CELL2:x");
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_eq!(i.len(), 2);
        assert_eq!(i.name(a), "CELL1");
    }

    #[test]
    fn maps_add_and_merge() {
        let mut a = CellMaps::default();
        a.add(Context::Cg, 100, 1);
        a.add(Context::Cg, 100, 0);
        a.add(Context::Ch, 200, 0);
        assert_eq!(a.cg, vec![(100, 1, 1)]);
        let mut b = CellMaps::default();
        b.add(Context::Cg, 100, 1);
        b.add(Context::Cg, 50, 0);
        a.merge_from(&b);
        assert_eq!(a.cg, vec![(50, 1, 0), (100, 1, 2)]);
    }

    #[test]
    fn out_of_order_add_normalizes() {
        let mut a = CellMaps::default();
        a.add(Context::Ch, 10, 0);
        a.add(Context::Ch, 20, 1);
        a.add(Context::Ch, 15, 0);
        a.add(Context::Ch, 20, 0);
        a.add(Context::Ch, 10, 1);
        a.normalize();
        assert_eq!(a.ch, vec![(10, 1, 1), (15, 1, 0), (20, 1, 1)]);
    }
}
