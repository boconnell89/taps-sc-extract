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
use std::io::{BufWriter, Read, Write};
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


/// `n` then `n` × (pos, t, c) little-endian u32s, read in one piece.
/// Sites are 12 bytes each (pos, t, c as little-endian u32). Moved through a
/// small stack buffer so each list costs one allocation.
const SITE_BLOCK: usize = 4096;

fn read_sites<R: Read>(r: &mut R) -> std::io::Result<Vec<(u32, u32, u32)>> {
    let n = read_u32(r)? as usize;
    let mut out = Vec::with_capacity(n);
    let mut buf = [0u8; 12 * SITE_BLOCK];
    let le = |b: &[u8]| u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    let mut left = n;
    while left > 0 {
        let k = left.min(SITE_BLOCK);
        let raw = &mut buf[..k * 12];
        r.read_exact(raw)?;
        out.extend(raw.chunks_exact(12).map(|s| (le(&s[0..4]), le(&s[4..8]), le(&s[8..12]))));
        left -= k;
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
/// Split a window's cells into per-bucket payloads. The cells' site vectors
/// are moved out (`cells` is left with empty maps).
pub fn window_shard_payloads(
    window: &Window,
    intern: &BarcodeIntern,
    cells: &mut [CellMaps],
    n_shards: usize,
    move_sites: bool,
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
            let cell = &mut cells[id as usize];
            // Moved when the payload is sent on at once (streamed writers);
            // cloned (exact size, no spare capacity) when memory mode keeps
            // it to the end of extraction.
            let (cg, ch) = if move_sites {
                (std::mem::take(&mut cell.cg), std::mem::take(&mut cell.ch))
            } else {
                (cell.cg.clone(), cell.ch.clone())
            };
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
        let mut w = BufWriter::with_capacity(
            1 << 20,
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
            write_sorted_sites(&mut w, &cell.cg)?;
            write_sorted_sites(&mut w, &cell.ch)?;
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
    // One read of the whole file; chunk files are a few MB.
    let bytes = fs::read(path).with_context(|| format!("open {}", path.display()))?;
    read_chunk_from(&mut bytes.as_slice())
        .with_context(|| format!("read {}", path.display()))?
        .ok_or_else(|| anyhow::anyhow!("empty chunk file {}", path.display()))
}

/// Next `TAPSCK01` record from a stream, or `None` at a clean end of stream.
pub fn read_chunk_from<R: Read>(r: &mut R) -> Result<Option<ChunkFile>> {
    let mut mag = [0u8; 8];
    let mut got = 0;
    while got < mag.len() {
        let n = r.read(&mut mag[got..])?;
        if n == 0 {
            if got == 0 {
                return Ok(None);
            }
            bail!("truncated chunk header");
        }
        got += n;
    }
    if &mag != MAGIC {
        bail!("bad chunk magic");
    }
    let chunk_id = read_u32(r)?;
    let clen = read_u16(r)? as usize;
    let mut cbuf = vec![0u8; clen];
    r.read_exact(&mut cbuf)?;
    let contig = String::from_utf8_lossy(&cbuf).into_owned();
    let n_cells = read_u32(r)? as usize;
    let mut cells = Vec::with_capacity(n_cells);
    for _ in 0..n_cells {
        let nlen = read_u16(r)? as usize;
        let mut nbuf = vec![0u8; nlen];
        r.read_exact(&mut nbuf)?;
        let barcode = String::from_utf8_lossy(&nbuf).into_owned();
        let cg = read_sites(r)?;
        let ch = read_sites(r)?;
        cells.push(ChunkCell { barcode, cg, ch });
    }
    Ok(Some(ChunkFile {
        chunk_id,
        contig,
        cells,
    }))
}

/// Read every `TAPSCK01` record from a stream until end of stream.
pub fn read_all_chunks_from<R: Read>(r: &mut R) -> Result<Vec<ChunkFile>> {
    let mut out = Vec::new();
    while let Some(ch) = read_chunk_from(r)? {
        out.push(ch);
    }
    Ok(out)
}

fn write_sorted_sites<W: Write>(w: &mut W, sites: &[(u32, u32, u32)]) -> std::io::Result<()> {
    write_u32(w, sites.len() as u32)?;
    let mut buf = [0u8; 12 * SITE_BLOCK];
    for block in sites.chunks(SITE_BLOCK) {
        for (o, &(p, t, c)) in buf.chunks_exact_mut(12).zip(block) {
            o[0..4].copy_from_slice(&p.to_le_bytes());
            o[4..8].copy_from_slice(&t.to_le_bytes());
            o[8..12].copy_from_slice(&c.to_le_bytes());
        }
        w.write_all(&buf[..block.len() * 12])?;
    }
    Ok(())
}

/// Serialize an in-memory chunk in the same `TAPSCK01` layout as the temp files.
pub fn write_chunk_to<W: Write>(w: &mut W, ch: &ChunkFile) -> Result<()> {
    w.write_all(MAGIC)?;
    write_u32(w, ch.chunk_id)?;
    let cb = ch.contig.as_bytes();
    if cb.len() > u16::MAX as usize {
        bail!("contig name too long");
    }
    write_u16(w, cb.len() as u16)?;
    w.write_all(cb)?;
    write_u32(w, ch.cells.len() as u32)?;
    for cell in &ch.cells {
        let name = cell.barcode.as_bytes();
        if name.len() > u16::MAX as usize {
            bail!("barcode too long");
        }
        write_u16(w, name.len() as u16)?;
        w.write_all(name)?;
        write_sorted_sites(w, &cell.cg)?;
        write_sorted_sites(w, &cell.ch)?;
    }
    Ok(())
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

        let mut buf = Vec::new();
        write_chunk_to(&mut buf, &fa).unwrap();
        write_chunk_to(&mut buf, &fa).unwrap();
        let back = read_all_chunks_from(&mut buf.as_slice()).unwrap();
        assert_eq!(back.len(), 2);
        assert_eq!(back[1].contig, fa.contig);
        assert_eq!(back[1].cells.len(), fa.cells.len());
        assert_eq!(back[1].cells[0].cg, fa.cells[0].cg);
        let _ = fs::remove_dir_all(&dir);
    }
}
