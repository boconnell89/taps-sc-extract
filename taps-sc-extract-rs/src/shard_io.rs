//! Compact per-(shard, window) temp files. One file per shard that has data.
//!
//! Layout (`temp_dir/shard_XXX/chunk_YYYYYY.bin`), little-endian:
//!
//! ```text
//! magic      : b"TAPSCK01"
//! chunk_id   : u32
//! contig_len : u16
//! contig     : [u8; contig_len]
//! n_cells    : u32
//! for each cell (barcode lexicographic):
//!   name_len : u16
//!   name     : [u8; name_len]
//!   n_cg     : u32
//!   n_cg × (pos u32, t u32, c u32)  // 1-based pos, sorted
//!   n_ch     : u32
//!   n_ch × (pos u32, t u32, c u32)
//! ```
//!
//! `pct` is not stored; the HDF5 writer computes `100*c/(c+t)` like Python.

use crate::accumulate::{BarcodeIntern, CellMaps};
use crate::barcode::barcode_to_shard;
use crate::window::Window;
use anyhow::{bail, Context, Result};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub const MAGIC: &[u8; 8] = b"TAPSCK01";

fn write_u16<W: Write>(w: &mut W, v: u16) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn write_u32<W: Write>(w: &mut W, v: u32) -> std::io::Result<()> {
    w.write_all(&v.to_le_bytes())
}
fn read_u16<R: Read>(r: &mut R) -> std::io::Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}
fn read_u32<R: Read>(r: &mut R) -> std::io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

use rustc_hash::FxHashMap;

fn write_sites<W: Write>(w: &mut W, map: &FxHashMap<u32, (u32, u32)>) -> std::io::Result<()> {
    let mut items: Vec<(u32, u32, u32)> = map.iter().map(|(&p, &(t, c))| (p, t, c)).collect();
    items.sort_unstable_by_key(|x| x.0);
    write_u32(w, items.len() as u32)?;
    for (p, t, c) in items {
        write_u32(w, p)?;
        write_u32(w, t)?;
        write_u32(w, c)?;
    }
    Ok(())
}

fn read_sites<R: Read>(r: &mut R) -> std::io::Result<Vec<(u32, u32, u32)>> {
    let n = read_u32(r)? as usize;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        out.push((read_u32(r)?, read_u32(r)?, read_u32(r)?));
    }
    Ok(out)
}

pub fn shard_dir(temp_dir: &Path, shard: usize) -> PathBuf {
    temp_dir.join(format!("shard_{shard:03}"))
}

/// Create `temp_dir` and `shard_XXX/` on the main thread (avoids mkdir races).
pub fn prepare_shard_dirs(temp_dir: &Path, n_shards: usize) -> Result<()> {
    fs::create_dir_all(temp_dir).with_context(|| format!("mkdir {}", temp_dir.display()))?;
    for s in 0..n_shards.max(1) {
        fs::create_dir_all(shard_dir(temp_dir, s))
            .with_context(|| format!("mkdir {}", shard_dir(temp_dir, s).display()))?;
    }
    Ok(())
}

pub fn chunk_path(temp_dir: &Path, shard: usize, chunk_id: u32) -> PathBuf {
    shard_dir(temp_dir, shard).join(format!("chunk_{chunk_id:06}.bin"))
}

pub fn list_chunk_files(shard_dir: &Path) -> Result<Vec<PathBuf>> {
    let mut v = Vec::new();
    if !shard_dir.is_dir() {
        return Ok(v);
    }
    for ent in fs::read_dir(shard_dir).with_context(|| format!("read {}", shard_dir.display()))? {
        let p = ent?.path();
        if p.extension().map(|e| e == "bin").unwrap_or(false) {
            v.push(p);
        }
    }
    v.sort();
    Ok(v)
}

/// Bucket this window into per-shard in-memory chunk payloads (chunk_id order later).
pub fn window_shard_payloads(
    window: &Window,
    intern: &BarcodeIntern,
    cells: &[CellMaps],
    n_shards: usize,
) -> Vec<(usize, ChunkFile)> {
    let n_shards = n_shards.max(1);
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_shards];
    for (i, name) in intern.names().iter().enumerate() {
        if i >= cells.len() || cells[i].is_empty() {
            continue;
        }
        buckets[barcode_to_shard(name, n_shards)].push(i as u32);
    }
    let mut out = Vec::new();
    for (s, ids) in buckets.into_iter().enumerate() {
        if ids.is_empty() {
            continue;
        }
        let mut ordered = ids;
        ordered.sort_unstable_by(|&a, &b| intern.name(a).cmp(intern.name(b)));
        let mut payload = ChunkFile {
            chunk_id: window.chunk_id,
            contig: window.contig.clone(),
            cells: Vec::with_capacity(ordered.len()),
        };
        for id in ordered {
            let cell = &cells[id as usize];
            let mut cg: Vec<(u32, u32, u32)> = cell.cg.iter().map(|(&p, &(t, c))| (p, t, c)).collect();
            let mut ch: Vec<(u32, u32, u32)> = cell.ch.iter().map(|(&p, &(t, c))| (p, t, c)).collect();
            cg.sort_unstable_by_key(|x| x.0);
            ch.sort_unstable_by_key(|x| x.0);
            payload.cells.push(ChunkCell {
                barcode: intern.name(id).to_string(),
                cg,
                ch,
            });
        }
        out.push((s, payload));
    }
    out
}

/// Write this window's maps into `temp_dir/shard_XXX/chunk_YYYYYY.bin`.
/// Empty shards are skipped (same as Python).
pub fn write_window_shards(
    temp_dir: &Path,
    window: &Window,
    intern: &BarcodeIntern,
    cells: &[CellMaps],
    n_shards: usize,
) -> Result<usize> {
    let n_shards = n_shards.max(1);
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_shards];
    for (i, name) in intern.names().iter().enumerate() {
        if i >= cells.len() || cells[i].is_empty() {
            continue;
        }
        let s = barcode_to_shard(name, n_shards);
        buckets[s].push(i as u32);
    }
    let mut n_files = 0usize;
    for (s, ids) in buckets.iter().enumerate() {
        if ids.is_empty() {
            continue;
        }
        let dir = shard_dir(temp_dir, s);
        fs::create_dir_all(&dir).with_context(|| format!("mkdir {}", dir.display()))?;
        let path = chunk_path(temp_dir, s, window.chunk_id);
        let mut w = BufWriter::new(
            File::create(&path).with_context(|| format!("create {}", path.display()))?,
        );
        w.write_all(MAGIC)?;
        write_u32(&mut w, window.chunk_id)?;
        let cb = window.contig.as_bytes();
        if cb.len() > u16::MAX as usize {
            bail!("contig name too long");
        }
        write_u16(&mut w, cb.len() as u16)?;
        w.write_all(cb)?;
        write_u32(&mut w, ids.len() as u32)?;
        let mut ordered = ids.clone();
        ordered.sort_unstable_by(|&a, &b| intern.name(a).cmp(intern.name(b)));
        for id in ordered {
            let name = intern.name(id).as_bytes();
            if name.len() > u16::MAX as usize {
                bail!("barcode too long");
            }
            write_u16(&mut w, name.len() as u16)?;
            w.write_all(name)?;
            let cell = &cells[id as usize];
            write_sites(&mut w, &cell.cg)?;
            write_sites(&mut w, &cell.ch)?;
        }
        w.flush()?;
        n_files += 1;
    }
    Ok(n_files)
}

#[derive(Debug, Clone)]
pub struct ChunkCell {
    pub barcode: String,
    pub cg: Vec<(u32, u32, u32)>,
    pub ch: Vec<(u32, u32, u32)>,
}

#[derive(Debug, Clone)]
pub struct ChunkFile {
    pub chunk_id: u32,
    pub contig: String,
    pub cells: Vec<ChunkCell>,
}

pub fn read_chunk_file(path: &Path) -> Result<ChunkFile> {
    let mut r = BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut mag = [0u8; 8];
    r.read_exact(&mut mag)?;
    if &mag != MAGIC {
        bail!("bad magic in {}", path.display());
    }
    let chunk_id = read_u32(&mut r)?;
    let clen = read_u16(&mut r)? as usize;
    let mut cbuf = vec![0u8; clen];
    r.read_exact(&mut cbuf)?;
    let contig = String::from_utf8_lossy(&cbuf).into_owned();
    let n_cells = read_u32(&mut r)? as usize;
    let mut cells = Vec::with_capacity(n_cells);
    for _ in 0..n_cells {
        let nlen = read_u16(&mut r)? as usize;
        let mut nbuf = vec![0u8; nlen];
        r.read_exact(&mut nbuf)?;
        let barcode = String::from_utf8_lossy(&nbuf).into_owned();
        let cg = read_sites(&mut r)?;
        let ch = read_sites(&mut r)?;
        cells.push(ChunkCell { barcode, cg, ch });
    }
    Ok(ChunkFile {
        chunk_id,
        contig,
        cells,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::accumulate::{ensure_cell, BarcodeIntern};
    use crate::context::Context;
    use crate::window::Window;

    #[test]
    fn roundtrip_two_shards() {
        let dir = std::env::temp_dir().join("taps_rs_chunk_io");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let mut intern = BarcodeIntern::default();
        let mut cells = Vec::new();
        let a = intern.intern_str("cellA"); // shard 7/16
        let x = intern.intern_str("x"); // shard 6/8 → 6/16?
        ensure_cell(&mut cells, a).add(Context::Cg, 100, 1);
        ensure_cell(&mut cells, a).add(Context::Ch, 101, 0);
        ensure_cell(&mut cells, x).add(Context::Cg, 200, 0);
        let w = Window {
            chunk_id: 3,
            contig: "chr19".into(),
            start: 0,
            end: 10,
        };
        let n = write_window_shards(&dir, &w, &intern, &cells, 16).unwrap();
        assert!(n >= 1);
        let pa = chunk_path(&dir, crate::barcode::barcode_to_shard("cellA", 16), 3);
        let fa = read_chunk_file(&pa).unwrap();
        assert_eq!(fa.chunk_id, 3);
        assert_eq!(fa.contig, "chr19");
        assert!(fa.cells.iter().any(|c| c.barcode == "cellA" && c.cg == vec![(100, 0, 1)]));
        let _ = fs::remove_dir_all(&dir);
    }
}
