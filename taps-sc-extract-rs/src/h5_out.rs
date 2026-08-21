//! Amethyst HDF5 writer: `/CG|CH/<barcode>/1` compound datasets and `master.h5`
//! relative ExternalLinks.
//!
//! Each writer thread owns a private file handle. Concurrent writers are
//! capped (default 6). Chunk size is `min(len, 65536)`.

use crate::autotune::MAX_SHARD_WRITERS;
use crate::shard_io::{list_chunk_files, read_chunk_file, shard_dir, ChunkFile};
use anyhow::{Context, Result};
use hdf5_metno::types::FixedAscii;
use hdf5_metno::{File, H5Type};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::path::Path;

pub const AMETHYST_VERSION: &str = "amethyst2.0.0";

#[derive(H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct MethRec {
    pub chr: FixedAscii<10>,
    pub pos: i64,
    pub pct: f64,
    pub t: i64,
    pub c: i64,
}

impl MethRec {
    pub fn new(contig: &str, pos: u32, t: u32, c: u32) -> Self {
        let tot = t + c;
        let pct = if tot == 0 {
            0.0
        } else {
            100.0 * f64::from(c) / f64::from(tot)
        };
        Self {
            chr: contig_s10(contig),
            pos: i64::from(pos),
            pct,
            t: i64::from(t),
            c: i64::from(c),
        }
    }
}

fn contig_s10(s: &str) -> FixedAscii<10> {
    let mut buf = [0u8; 10];
    let b = s.as_bytes();
    let n = b.len().min(10);
    buf[..n].copy_from_slice(&b[..n]);
    FixedAscii::<10>::from_ascii(&buf).unwrap_or_else(|_| {
        FixedAscii::<10>::from_ascii(b"chrN\0\0\0\0\0\0").expect("literal")
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum H5Compression {
    Gzip(u8),
    GzipShuffle(u8),
    None,
}

impl H5Compression {
    pub fn parse(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "gzip" | "gzip-fast" | "gzip1" => Ok(Self::Gzip(1)),
            "gzip-shuffle" | "shuffle" => Ok(Self::GzipShuffle(1)),
            "gzip6" => Ok(Self::Gzip(6)),
            "none" => Ok(Self::None),
            other => anyhow::bail!("unsupported --compression {other} (options: gzip, gzip-shuffle, gzip6, none)"),
        }
    }
}

fn write_version(file: &File) -> Result<()> {
    let meta = file.create_group("metadata")?;
    let v = AMETHYST_VERSION.as_bytes();
    meta.new_dataset::<u8>()
        .shape(v.len())
        .create("version")?
        .write(v)
        .map_err(|e| anyhow::anyhow!("write metadata/version: {e}"))?;
    Ok(())
}

fn write_cell_dataset(
    file: &File,
    context: &str,
    barcode: &str,
    recs: &[MethRec],
    compression: H5Compression,
) -> Result<()> {
    if recs.is_empty() {
        return Ok(());
    }
    let g = if file.link_exists(context) {
        file.group(context)?
    } else {
        file.create_group(context)?
    };
    let bc = g.create_group(barcode)?;
    let n = recs.len();
    let chunk = n.min(65536);
    let mut b = bc.new_dataset::<MethRec>().shape(n).chunk(chunk);
    b = match compression {
        H5Compression::Gzip(level) => b.deflate(level),
        H5Compression::GzipShuffle(level) => b.shuffle().deflate(level),
        H5Compression::None => b,
    };
    b.create("1")?
        .write(recs)
        .map_err(|e| anyhow::anyhow!("write {context}/{barcode}/1: {e}"))?;
    Ok(())
}

/// Concatenate window chunks (already in chunk_id order) into Amethyst records.
pub fn records_from_chunks(chunks: &[ChunkFile]) -> FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)> {
    let mut out: FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)> = FxHashMap::default();
    for ch in chunks {
        for cell in &ch.cells {
            let e = out.entry(cell.barcode.clone()).or_default();
            e.0.extend(
                cell.cg
                    .iter()
                    .map(|&(pos, t, c)| MethRec::new(&ch.contig, pos, t, c)),
            );
            e.1.extend(
                cell.ch
                    .iter()
                    .map(|&(pos, t, c)| MethRec::new(&ch.contig, pos, t, c)),
            );
        }
    }
    out
}

pub fn write_shard_file(
    path: &Path,
    cells: &FxHashMap<String, (Vec<MethRec>, Vec<MethRec>)>,
    compression: H5Compression,
) -> Result<usize> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).ok();
    }
    let file = File::create(path).map_err(|e| anyhow::anyhow!("create {}: {e}", path.display()))?;
    write_version(&file)?;
    let mut bcs: Vec<&String> = cells.keys().collect();
    bcs.sort();
    for bc in bcs {
        let (cg, ch) = &cells[bc];
        write_cell_dataset(&file, "CG", bc, cg, compression)?;
        write_cell_dataset(&file, "CH", bc, ch, compression)?;
    }
    file.flush().ok();
    Ok(cells.len())
}

pub fn write_master_h5(out_dir: &Path, n_shards: usize) -> Result<std::path::PathBuf> {
    let master_path = out_dir.join("master.h5");
    let file = File::create(&master_path)
        .map_err(|e| anyhow::anyhow!("create {}: {e}", master_path.display()))?;
    write_version(&file)?;
    let cg = file.create_group("CG")?;
    let ch = file.create_group("CH")?;
    for s in 0..n_shards.max(1) {
        let shard_fn = format!("shard_{s:03}.h5");
        let shard_path = out_dir.join(&shard_fn);
        if !shard_path.exists() {
            continue;
        }
        let shard = File::open(&shard_path)
            .map_err(|e| anyhow::anyhow!("open {}: {e}", shard_path.display()))?;
        if let Ok(g) = shard.group("CG") {
            for bc in g.member_names()? {
                cg.link_external(&shard_fn, &format!("CG/{bc}"), &bc)
                    .map_err(|e| anyhow::anyhow!("CG external link {bc}: {e}"))?;
            }
        }
        if let Ok(g) = shard.group("CH") {
            for bc in g.member_names()? {
                ch.link_external(&shard_fn, &format!("CH/{bc}"), &bc)
                    .map_err(|e| anyhow::anyhow!("CH external link {bc}: {e}"))?;
            }
        }
    }
    file.flush().ok();
    Ok(master_path)
}

fn assemble_one_shard_from_temp(
    temp_dir: &Path,
    shard: usize,
    shard_h5: &Path,
    compression: H5Compression,
) -> Result<usize> {
    let dir = shard_dir(temp_dir, shard);
    let mut chunks = Vec::new();
    if dir.is_dir() {
        for p in list_chunk_files(&dir)? {
            chunks.push(read_chunk_file(&p)?);
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
    chunks.sort_by_key(|c| c.chunk_id);
    let cells = records_from_chunks(&chunks);
    if cells.is_empty() {
        // still create an empty shard so master layout is stable? Python writes empty shard files.
        write_shard_file(shard_h5, &cells, compression)?;
        return Ok(0);
    }
    write_shard_file(shard_h5, &cells, compression)
}

/// Write `shard_XXX.h5` (+ `master.h5` if `n_shards > 1`) from stream temp chunks.
pub fn assemble_hdf5_from_temp(
    temp_dir: &Path,
    out: &Path,
    n_shards: usize,
    compression: H5Compression,
    max_writers: usize,
) -> Result<std::path::PathBuf> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = n_shards.max(1);
    let n_writers = n_shards.min(max_writers.max(1)).min(MAX_SHARD_WRITERS);
    if n_shards == 1 {
        let path = if out.extension().map(|e| e == "h5").unwrap_or(false) {
            out.to_path_buf()
        } else {
            std::fs::create_dir_all(out).ok();
            out.join("output.h5")
        };
        assemble_one_shard_from_temp(temp_dir, 0, &path, compression)?;
        return Ok(path);
    }
    std::fs::create_dir_all(out)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_writers)
        .thread_name(|i| format!("taps-h5-{i}"))
        .build()
        .context("HDF5 writer pool")?;
    let out = out.to_path_buf();
    let temp_dir = temp_dir.to_path_buf();
    pool.install(|| {
        (0..n_shards).into_par_iter().try_for_each(|s| {
            let p = out.join(format!("shard_{s:03}.h5"));
            let n = assemble_one_shard_from_temp(&temp_dir, s, &p, compression)?;
            eprintln!("HDF5 shard {s:03} cells={n} -> {}", p.display());
            Ok::<(), anyhow::Error>(())
        })
    })?;
    write_master_h5(&out, n_shards)
}

pub fn assemble_hdf5_from_memory(
    per_shard: &[Vec<ChunkFile>],
    out: &Path,
    compression: H5Compression,
    max_writers: usize,
) -> Result<std::path::PathBuf> {
    std::env::set_var("HDF5_USE_FILE_LOCKING", "FALSE");
    let n_shards = per_shard.len().max(1);
    let n_writers = n_shards.min(max_writers.max(1)).min(MAX_SHARD_WRITERS);
    if n_shards == 1 {
        let path = if out.extension().map(|e| e == "h5").unwrap_or(false) {
            out.to_path_buf()
        } else {
            std::fs::create_dir_all(out).ok();
            out.join("output.h5")
        };
        let mut chunks = per_shard.first().cloned().unwrap_or_default();
        chunks.sort_by_key(|c| c.chunk_id);
        let cells = records_from_chunks(&chunks);
        write_shard_file(&path, &cells, compression)?;
        return Ok(path);
    }
    std::fs::create_dir_all(out)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_writers)
        .thread_name(|i| format!("taps-h5-{i}"))
        .build()
        .context("HDF5 writer pool")?;
    let out_b = out.to_path_buf();
    pool.install(|| {
        per_shard.par_iter().enumerate().try_for_each(|(s, v)| {
            let mut chunks = v.clone();
            chunks.sort_by_key(|c| c.chunk_id);
            let cells = records_from_chunks(&chunks);
            let p = out_b.join(format!("shard_{s:03}.h5"));
            write_shard_file(&p, &cells, compression)?;
            eprintln!("HDF5 shard {s:03} cells={} -> {}", cells.len(), p.display());
            Ok::<(), anyhow::Error>(())
        })
    })?;
    write_master_h5(out, n_shards)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shard_io::ChunkCell;

    #[test]
    fn s10_and_pct() {
        let r = MethRec::new("chr19", 100, 1, 1);
        assert_eq!(r.pos, 100);
        assert!((r.pct - 50.0).abs() < 1e-9);
        assert_eq!(r.t, 1);
        assert_eq!(r.c, 1);
        let s = r.chr.as_str();
        assert!(s.starts_with("chr19"));
    }

    #[test]
    fn roundtrip_single_file() {
        let dir = std::env::temp_dir().join("taps_rs_h5_rt");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let chunks = vec![ChunkFile {
            chunk_id: 0,
            contig: "chr19".into(),
            cells: vec![ChunkCell {
                barcode: "cellA".into(),
                cg: vec![(100, 1, 1), (200, 0, 2)],
                ch: vec![(150, 3, 0)],
            }],
        }];
        let cells = records_from_chunks(&chunks);
        let p = dir.join("out.h5");
        write_shard_file(&p, &cells, H5Compression::None).unwrap();
        let f = File::open(&p).unwrap();
        let ds = f.dataset("CG/cellA/1").unwrap();
        let recs: Vec<MethRec> = ds.read_raw().unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].pos, 100);
        assert_eq!(recs[1].pos, 200);
        assert_eq!(recs[0].t, 1);
        assert_eq!(recs[0].c, 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
