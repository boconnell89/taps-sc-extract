//! Noodles-based column-wise window extract (no htslib pileup).
//!
//! Walks CIGAR on records overlapping `[start, end)` and applies the same
//! TAPS calling as the rust-htslib path. No BAQ (htslib `sam_prob_realn` is
//! not available here). Overlap clip is per-QNAME higher-quality base, which
//! approximates but is not identical to `bam_mplp_init_overlaps`.

use crate::accumulate::{ensure_cell, merge_window_cells, BarcodeIntern};
use crate::barcode::barcode_from_qname_bytes;
use crate::calling::{call_mctot, classify_strand};
use crate::context::{classify_trinucleotide, Context, TriContext};
use crate::extract::{CallStats, ExtractParams, WindowResult};
use crate::fasta::FastFaiReader;
use crate::window::Window;
use anyhow::{Context as AnyhowContext, Result};
use noodles::bam;
use noodles::core::Region;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::collections::{HashMap, HashSet};
use std::path::Path;

struct Hit {
    qname: Vec<u8>,
    qual: u8,
    read_base: u8,
    flag: u16,
}

/// Query records overlapping a 0-based half-open window and pile up by CIGAR.
pub fn process_window_noodles(
    reader: &mut bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>,
    header: &noodles::sam::Header,
    fasta: &FastFaiReader,
    window: &Window,
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    let pad = 2u64;
    let ref_len = fasta
        .reference_length(&window.contig)
        .ok_or_else(|| anyhow::anyhow!("contig {} missing in FASTA", window.contig))?;
    let ref_start = window.start.saturating_sub(pad);
    let ref_end = (window.end + pad).min(ref_len);
    let ref_seq = fasta
        .fetch(&window.contig, ref_start, ref_end)
        .with_context(|| "FASTA fetch")?;

    // noodles Region is 1-based inclusive.
    let region: Region = format!(
        "{}:{}-{}",
        window.contig,
        window.start + 1,
        window.end.max(window.start + 1)
    )
    .parse()
    .map_err(|e| anyhow::anyhow!("region parse: {e}"))?;

    let query = reader
        .query(header, &region)
        .map_err(|e| anyhow::anyhow!("noodles query: {e}"))?;

    let win_start = window.start;
    let win_end = window.end;
    let mut columns: HashMap<u32, Vec<Hit>> = HashMap::new();

    for rec in query.records() {
        let rec = rec.map_err(|e| anyhow::anyhow!("noodles record: {e}"))?;
        let flags = rec.flags();
        let flag_bits = flags.bits();
        if classify_strand(flag_bits).is_none() {
            continue;
        }
        let mapq = rec.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < params.min_mapq {
            continue;
        }
        let Some(Ok(aln_start)) = rec.alignment_start() else {
            continue;
        };
        let mut ref_pos = (aln_start.get() as u64).saturating_sub(1);
        let mut qpos: usize = 0;
        let seq = rec.sequence();
        let qual = rec.quality_scores();
        let qual_slice = qual.as_ref();
        let qname = rec
            .name()
            .map(|n| <[u8]>::to_vec(n))
            .unwrap_or_default();

        for op_res in rec.cigar().iter() {
            let op = op_res.map_err(|e| anyhow::anyhow!("cigar: {e}"))?;
            let len = op.len();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..len {
                        if ref_pos >= win_start && ref_pos < win_end {
                            if qpos < seq.len() {
                                let qb = qual_slice.get(qpos).copied().unwrap_or(0);
                                if qb >= params.min_baseq {
                                    if let Some(base) = seq.get(qpos) {
                                        columns.entry(ref_pos as u32).or_default().push(Hit {
                                            qname: qname.clone(),
                                            qual: qb,
                                            read_base: base,
                                            flag: flag_bits,
                                        });
                                    }
                                }
                            }
                        }
                        ref_pos += 1;
                        qpos += 1;
                    }
                }
                Kind::Insertion | Kind::SoftClip => {
                    qpos += len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += len as u64;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    let mut stats = CallStats::default();
    let mut intern = BarcodeIntern::default();
    let mut cells = Vec::new();

    let mut positions: Vec<u32> = columns.keys().copied().collect();
    positions.sort_unstable();
    for pos in positions {
        let mut hits = columns.remove(&pos).unwrap_or_default();
        if hits.len() > params.max_depth as usize {
            hits.truncate(params.max_depth as usize);
        }
        if params.ignore_overlaps && hits.len() > 1 {
            let mut best: HashMap<Vec<u8>, usize> = HashMap::new();
            for (i, h) in hits.iter().enumerate() {
                match best.get(&h.qname) {
                    None => {
                        best.insert(h.qname.clone(), i);
                    }
                    Some(&j) => {
                        if h.qual > hits[j].qual {
                            best.insert(h.qname.clone(), i);
                        }
                    }
                }
            }
            let keep: HashSet<usize> = best.into_values().collect();
            hits = hits
                .into_iter()
                .enumerate()
                .filter_map(|(i, h)| keep.contains(&i).then_some(h))
                .collect();
        }

        let rel = (pos as u64).saturating_sub(ref_start) as usize;
        if rel >= ref_seq.len() {
            continue;
        }
        let ref_b = ref_seq[rel];
        let Some(tri) = classify_trinucleotide(&ref_seq, rel) else {
            continue;
        };
        let ctx = match tri {
            TriContext::Cpg => Context::Cg,
            _ => Context::Ch,
        };
        for h in hits {
            let Some(strand) = classify_strand(h.flag) else {
                continue;
            };
            let Some(call) = call_mctot(strand, ref_b, h.read_base) else {
                continue;
            };
            let bc_bytes = barcode_from_qname_bytes(&h.qname);
            if let Some(wl) = whitelist {
                let bc = std::str::from_utf8(bc_bytes).unwrap_or("");
                if !wl.contains(bc) {
                    continue;
                }
            }
            let id = intern.intern_bytes(bc_bytes);
            if params.accumulate {
                ensure_cell(&mut cells, id).add(ctx, pos + 1, call);
            }
            stats.add(tri, ctx, call);
        }
    }

    Ok(WindowResult {
        stats,
        intern,
        cells,
    })
}

pub fn extract_regions_stats_noodles(
    bam_path: &Path,
    fasta_path: &Path,
    windows: &[Window],
    whitelist: Option<&HashSet<String>>,
    params: &ExtractParams,
) -> Result<WindowResult> {
    let fasta = FastFaiReader::open(fasta_path)?;
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| anyhow::anyhow!("noodles open {}: {e}", bam_path.display()))?;
    let header = reader
        .read_header()
        .map_err(|e| anyhow::anyhow!("noodles header: {e}"))?;

    let mut total = WindowResult {
        stats: CallStats::default(),
        intern: BarcodeIntern::default(),
        cells: Vec::new(),
    };
    for (i, w) in windows.iter().enumerate() {
        let r = process_window_noodles(&mut reader, &header, &fasta, w, whitelist, params)?;
        total.stats.merge(&r.stats);
        if params.accumulate {
            merge_window_cells(&mut total.intern, &mut total.cells, &r.intern, &r.cells);
        } else {
            for name in r.intern.names() {
                total.intern.intern_str(name);
            }
        }
        if (i + 1) % 5 == 0 || i + 1 == windows.len() {
            eprintln!(
                "Rust noodles: [{}/{}] cells={} CG={:.2}% CH={:.3}%",
                i + 1,
                windows.len(),
                total.n_cells(),
                total.stats.cg_pct(),
                total.stats.ch_pct()
            );
        }
    }
    Ok(total)
}
