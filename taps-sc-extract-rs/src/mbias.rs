//! M-bias from a sample of reads, and the automatic trim rule (`--trim auto`).
//!
//! A port of `taps_sc_extract/mbias.py`: the same windows, count table, TSV
//! and rule, so both engines choose identical trims (tested). Design:
//! code/bisulfite-extraction-plan.md §5.1. Positions count from the read's
//! sequenced end in BAM SEQ coordinates (soft clips included, as the trims
//! are), on reads at least `2 * P` long; the plateau counts bases at least `P`
//! from both ends of reads at least `PLATEAU_MIN_LEN` long.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::Result;
use rust_htslib::bam::{self, Read};

use crate::calling::Strand;
use crate::chemistry::Chemistry;
use crate::fasta::FastFaiReader;

pub const P: usize = 20;
pub const PLATEAU_MIN_LEN: usize = 60;
pub const LENGTH_BINS: [(usize, usize); 4] = [(40, 59), (60, 99), (100, 149), (150, 1_000_000)];
pub const N_WINDOWS: usize = 64;
pub const TARGET_READS: u64 = 300_000;
pub const MIN_WINDOW: u64 = 10_000;
pub const Z: f64 = 3.0;
pub const STABLE_RUN: usize = 3;
pub const MAX_AUTO_TRIM: u32 = 15;
pub const FALLBACK_HEAD: usize = 5;
/// Reference padding past a sampling window for reads that start inside it.
const REF_PAD: u64 = 2_000;

fn delta(ctx: &str) -> f64 {
    if ctx == "CG" {
        0.02
    } else {
        0.01
    }
}

fn candidate_3p(chem: Chemistry) -> u32 {
    match chem {
        Chemistry::Taps => 0,
        Chemistry::Bisulfite => 10,
    }
}

fn length_bin(n: usize) -> Option<usize> {
    LENGTH_BINS.iter().position(|&(lo, hi)| lo <= n && n <= hi)
}

pub fn bin_label(i: usize) -> String {
    let (lo, hi) = LENGTH_BINS[i];
    if hi < 1_000_000 {
        format!("{lo}-{hi}")
    } else {
        format!("{lo}+")
    }
}

/// (context, mate, end, length bin, position) -> [methylated, unmethylated].
pub type Table = BTreeMap<(&'static str, u8, u8, usize, usize), [u64; 2]>;

pub struct Counts {
    pub table: Table,
    /// CG and CH plateau [methylated, unmethylated].
    pub plateau: [[u64; 2]; 2],
    pub reads: u64,
}

/// Windows to sample: every read if the selection holds at most
/// `TARGET_READS` mapped reads, else `N_WINDOWS` windows dealt to contigs by
/// mapped reads (largest remainder) and spaced evenly. Same formula as
/// `mbias.sample_windows` (Python).
pub fn sample_windows(contigs: &[(String, u64, u64)], start: Option<u64>, end: Option<u64>) -> Vec<(String, u64, u64)> {
    if start.is_some() || end.is_some() {
        let (name, len, _) = &contigs[0];
        return vec![(name.clone(), start.unwrap_or(0), end.unwrap_or(*len))];
    }
    let total_reads: u64 = contigs.iter().map(|c| c.2).sum();
    if total_reads <= TARGET_READS {
        return contigs.iter().map(|(n, l, _)| (n.clone(), 0, *l)).collect();
    }
    let total_len: u64 = contigs.iter().filter(|c| c.2 > 0).map(|c| c.1).sum();
    let width = ((TARGET_READS as f64 / total_reads as f64 * total_len as f64 / N_WINDOWS as f64) as u64).max(MIN_WINDOW);
    let quotas: Vec<f64> = contigs.iter().map(|c| N_WINDOWS as f64 * c.2 as f64 / total_reads as f64).collect();
    let mut counts: Vec<usize> = quotas.iter().map(|q| *q as usize).collect();
    let mut order: Vec<usize> = (0..contigs.len()).collect();
    order.sort_by(|&a, &b| {
        let ra = quotas[a] - counts[a] as f64;
        let rb = quotas[b] - counts[b] as f64;
        rb.partial_cmp(&ra).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b))
    });
    let left = N_WINDOWS.saturating_sub(counts.iter().sum());
    for &i in order.iter().take(left) {
        counts[i] += 1;
    }
    let mut out = Vec::new();
    for ((name, len, _), &k) in contigs.iter().zip(&counts) {
        for j in 0..k as u64 {
            let centre = (2 * j + 1) * len / (2 * k as u64);
            let ws = centre.saturating_sub(width / 2);
            out.push((name.clone(), ws, (ws + width).min(*len)));
        }
    }
    out
}

/// Count methylation by position from each sequenced end over `windows`.
pub fn count(
    bam_path: &Path,
    fasta: &FastFaiReader,
    windows: &[(String, u64, u64)],
    chem: Chemistry,
    min_mapq: u8,
    min_baseq: u8,
) -> Result<Counts> {
    let mut r = crate::input::open_indexed(bam_path)?;
    let skip: u16 = 0x704 | if chem == Chemistry::Bisulfite { 0x800 } else { 0 };
    let mut c = Counts { table: Table::new(), plateau: [[0; 2]; 2], reads: 0 };
    let mut rec = bam::Record::new();
    for (contig, ws, we) in windows {
        let clen = fasta.reference_length(contig).unwrap_or(0);
        // One reference fetch per window, padded for reads that run past its end.
        let (wr0, wr1) = (ws.saturating_sub(2), (we + REF_PAD).min(clen));
        let wref = fasta.fetch(contig, wr0, wr1)?;
        r.fetch((contig.as_str(), *ws as i64, *we as i64)).map_err(|e| anyhow::anyhow!("M-bias fetch {contig}: {e}"))?;
        while let Some(res) = r.read(&mut rec) {
            res.map_err(|e| anyhow::anyhow!("M-bias read: {e}"))?;
            let pos = rec.pos();
            if pos < *ws as i64 || pos >= *we as i64 {
                continue; // each read once: in the window holding its start
            }
            let flag = rec.flags();
            if flag & skip != 0 || rec.mapq() < min_mapq {
                continue;
            }
            let Some(strand) = chem.strand(flag) else { continue };
            let seq = rec.seq();
            let len = seq.len();
            if len == 0 {
                continue;
            }
            c.reads += 1;
            let qual = rec.qual();
            let reverse = flag & 16 != 0;
            let mate: u8 = if flag & 128 != 0 { 2 } else { 1 };
            let lb = if len >= 2 * P { length_bin(len) } else { None };
            let rlen: i64 = rec.raw_cigar().iter().filter(|&&x| matches!(x & 0xf, 0 | 2 | 3 | 7 | 8)).map(|&x| (x >> 4) as i64).sum();
            let r0 = (pos - 2).max(0) as u64;
            let r1 = ((pos + rlen + 2) as u64).min(clen);
            let own;
            let ref_seq: &[u8] = if r0 >= wr0 && r1 <= wr1 {
                &wref[(r0 - wr0) as usize..(r1 - wr0) as usize]
            } else {
                own = fasta.fetch(contig, r0, r1)?;
                &own
            };
            let (want, mut q, mut rp) = (if strand == Strand::Ot { b'C' } else { b'G' }, 0usize, pos);
            for &op in rec.raw_cigar() {
                let n = (op >> 4) as usize;
                match op & 0xf {
                    0 | 7 | 8 => {
                        for k in 0..n {
                            let qp = q + k;
                            let i = (rp + k as i64) as u64 - r0;
                            let i = i as usize;
                            let Some(&rb) = ref_seq.get(i) else { continue };
                            if rb != want || qual[qp] < min_baseq {
                                continue;
                            }
                            let Some(call) = chem.call(strand, rb, seq[qp]) else { continue };
                            let ctx = if rb == b'C' {
                                if ref_seq.get(i + 1) == Some(&b'G') { "CG" } else { "CH" }
                            } else if i >= 1 && ref_seq[i - 1] == b'C' {
                                "CG"
                            } else {
                                "CH"
                            };
                            let (d5, d3) = if reverse { (len - 1 - qp, qp) } else { (qp, len - 1 - qp) };
                            let idx = if call == 1 { 0 } else { 1 };
                            if let Some(lb) = lb {
                                if d5 < P {
                                    c.table.entry((ctx, mate, 5, lb, d5)).or_insert([0; 2])[idx] += 1;
                                }
                                if d3 < P {
                                    c.table.entry((ctx, mate, 3, lb, d3)).or_insert([0; 2])[idx] += 1;
                                }
                            }
                            if len >= PLATEAU_MIN_LEN && d5 >= P && d3 >= P {
                                c.plateau[if ctx == "CG" { 0 } else { 1 }][idx] += 1;
                            }
                        }
                        q += n;
                        rp += n as i64;
                    }
                    1 | 4 => q += n,
                    2 | 3 => rp += n as i64,
                    _ => {}
                }
            }
        }
    }
    Ok(c)
}

fn curve(t: &Table, ctx: &'static str, mate: u8, end: u8) -> Vec<(u64, u64)> {
    (0..P)
        .map(|p| {
            let mut s = (0u64, 0u64);
            for lb in 0..LENGTH_BINS.len() {
                if let Some(v) = t.get(&(ctx, mate, end, lb, p)) {
                    s.0 += v[0];
                    s.1 += v[1];
                }
            }
            s
        })
        .collect()
}

/// Trim for one end (see `mbias.walk`): (trim, unknown positions near the end).
fn walk(cur: &[(u64, u64)], plateau_rate: f64, delta: f64, min_calls: u64) -> (u32, bool) {
    let (mut trim, mut run, mut unknown_head) = (0u32, 0usize, false);
    for (p, &(m, u)) in cur.iter().enumerate() {
        let n = m + u;
        if n < min_calls {
            if p < FALLBACK_HEAD {
                unknown_head = true;
            }
            run += 1;
        } else {
            let rate = m as f64 / n as f64;
            let se = (rate * (1.0 - rate) / n as f64).sqrt();
            if (rate - plateau_rate).abs() > delta.max(Z * se) {
                trim = p as u32 + 1;
                run = 0;
            } else {
                run += 1;
            }
        }
        if run >= STABLE_RUN {
            break;
        }
    }
    (trim.min(MAX_AUTO_TRIM), unknown_head)
}

pub struct Chosen {
    pub trim_r1: u32,
    pub trim_r2: u32,
    pub trim_3p: u32,
    pub notes: Vec<String>,
}

/// Calls a read-end position needs before it can set or stop the trim
/// (fixed per position and context; `mbias.min_calls_default`).
pub const MIN_CALLS: u64 = 500;

pub fn min_calls(_plateau_n: u64) -> u64 {
    std::env::var("TAPS_MBIAS_MIN_CALLS").ok().and_then(|v| v.parse().ok()).unwrap_or(MIN_CALLS)
}

/// The rule: per end, the CG walk; CH only if CG lacks data (else a warning
/// when CH wants more); the chemistry candidate if the read end has too few
/// calls. Same as `mbias.choose_trims`.
pub fn choose(c: &Counts, chem: Chemistry, mates: &[u8]) -> Chosen {
    let pn = [c.plateau[0][0] + c.plateau[0][1], c.plateau[1][0] + c.plateau[1][1]];
    let mc = min_calls(pn[0] + pn[1]);
    let mut out = Chosen { trim_r1: 0, trim_r2: 0, trim_3p: 0, notes: Vec::new() };
    let mut three = Vec::new();
    for &mate in mates {
        for end in [5u8, 3] {
            let mut per: [Option<(u32, bool)>; 2] = [None, None];
            for (ci, ctx) in ["CG", "CH"].into_iter().enumerate() {
                if pn[ci] < mc {
                    continue;
                }
                let rate = c.plateau[ci][0] as f64 / pn[ci] as f64;
                per[ci] = Some(walk(&curve(&c.table, ctx, mate, end), rate, delta(ctx), mc));
            }
            let label = if mates.len() > 1 { format!("R{mate} ") } else { String::new() };
            let t = match per[0] {
                Some((t, false)) => {
                    if let Some((tch, false)) = per[1] {
                        if tch > t {
                            out.notes.push(format!("{label}{end}': CH suggests {tch} bp; using CG {t} bp"));
                        }
                    }
                    t
                }
                _ => match per[1] {
                    Some((tch, false)) => {
                        out.notes.push(format!("{label}{end}': too few CG calls; CH sets the trim ({tch} bp)"));
                        tch
                    }
                    _ => {
                        let f = if end == 3 { candidate_3p(chem) } else { 0 };
                        out.notes.push(format!("{label}{end}': too few calls near the read end; using the default {f} bp"));
                        f
                    }
                },
            };
            if t >= MAX_AUTO_TRIM {
                out.notes.push(format!("{label}{end}': M-bias reaches the {MAX_AUTO_TRIM} bp cap; check the library"));
            }
            if end == 5 {
                if mate == 2 {
                    out.trim_r2 = t;
                } else {
                    out.trim_r1 = t;
                }
            } else {
                three.push(t);
            }
        }
    }
    out.trim_3p = three.into_iter().max().unwrap_or(0);
    out
}

/// Same format as `mbias.write_tsv`.
pub fn write_tsv(path: &Path, c: &Counts, trim_r1: u32, trim_r2: u32, trim_3p: u32) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "# reads\t{}", c.reads)?;
    writeln!(f, "# plateau\tCG\t{}\t{}", c.plateau[0][0], c.plateau[0][1])?;
    writeln!(f, "# plateau\tCH\t{}\t{}", c.plateau[1][0], c.plateau[1][1])?;
    writeln!(f, "# trims\t5p={trim_r1}\tr2_5p={trim_r2}\t3p={trim_3p}")?;
    writeln!(f, "context\tmate\tend\tlength_bin\tposition\tmethylated\tunmethylated")?;
    for (&(ctx, mate, end, lb, p), v) in &c.table {
        writeln!(f, "{ctx}\t{mate}\t{end}\t{}\t{p}\t{}\t{}", bin_label(lb), v[0], v[1])?;
    }
    Ok(())
}

/// `<dir>/<bam name without extension>`: the shard directory, or beside a single file.
pub fn output_prefix(out: &Path, n_shards: usize, bam: &Path) -> PathBuf {
    let is_h5 = out.extension().is_some_and(|e| e == "h5" || e == "hdf5");
    let dir = if n_shards > 1 || out.is_dir() || !is_h5 {
        if is_h5 && n_shards > 1 {
            out.with_extension("")
        } else {
            out.to_path_buf()
        }
    } else {
        out.parent().filter(|p| !p.as_os_str().is_empty()).map(Path::to_path_buf).unwrap_or_else(|| PathBuf::from("."))
    };
    let name = bam.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let stem = [".bam", ".cram", ".sam"].iter().find_map(|e| name.strip_suffix(e)).unwrap_or(&name).to_string();
    dir.join(stem)
}

/// Draw `<prefix>.mbias.png` / `.pdf` with the Python plotter (seaborn).
pub fn plot(tsv: &Path, title: &str) {
    let ok = std::process::Command::new("python3")
        .args(["-m", "taps_sc_extract.mbias", "--plot"])
        .arg(tsv)
        .args(["--title", title])
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !ok {
        eprintln!(
            "note: M-bias plot not drawn (needs the taps_sc_extract Python package with seaborn); run: \
             python3 -m taps_sc_extract.mbias --plot {}",
            tsv.display()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A table whose CG and CH curves are flat at the plateau, with the
    /// listed (end, position) overrides of CG methylation (as a fraction).
    fn table(cg_override: &[(u8, usize, f64)], ch_override: &[(u8, usize, f64)], n: u64) -> Counts {
        let mut t = Table::new();
        for (ctx, base, ov) in [("CG", 0.8, cg_override), ("CH", 0.05, ch_override)] {
            for end in [5u8, 3] {
                for p in 0..P {
                    let rate = ov.iter().find(|o| o.0 == end && o.1 == p).map(|o| o.2).unwrap_or(base);
                    let m = (rate * n as f64).round() as u64;
                    t.insert((ctx, 1, end, 1, p), [m, n - m]);
                }
            }
        }
        Counts { table: t, plateau: [[80_000, 20_000], [5_000, 95_000]], reads: 1 }
    }

    #[test]
    fn rule_picks_the_biased_end_only() {
        let flat = choose(&table(&[], &[], 5000), Chemistry::Bisulfite, &[1]);
        assert_eq!((flat.trim_r1, flat.trim_3p), (0, 0));
        let step: Vec<(u8, usize, f64)> = (0..10).map(|p| (3u8, p, 0.55)).collect();
        let c = choose(&table(&step, &[], 5000), Chemistry::Bisulfite, &[1]);
        assert_eq!((c.trim_r1, c.trim_3p), (0, 10));
        // An interior blip after recovery does not move the trim.
        let mut blip = step.clone();
        blip.push((3, 16, 0.6));
        assert_eq!(choose(&table(&blip, &[], 5000), Chemistry::Bisulfite, &[1]).trim_3p, 10);
        // CH wanting more only warns.
        let ch: Vec<(u8, usize, f64)> = (0..12).map(|p| (3u8, p, 0.02)).collect();
        let c = choose(&table(&step, &ch, 5000), Chemistry::Bisulfite, &[1]);
        assert_eq!(c.trim_3p, 10);
        assert!(c.notes.iter().any(|n| n.contains("CH suggests 12")));
        // Too few calls near the end: chemistry default (3' 10 for bisulfite).
        let sparse = table(&[], &[], 10);
        let c = choose(&sparse, Chemistry::Bisulfite, &[1]);
        assert_eq!(c.trim_3p, 10);
        // The cap.
        let long: Vec<(u8, usize, f64)> = (0..P).map(|p| (5u8, p, 0.5)).collect();
        let c = choose(&table(&long, &[], 5000), Chemistry::Bisulfite, &[1]);
        assert_eq!(c.trim_r1, MAX_AUTO_TRIM);
    }

    #[test]
    fn windows_are_deterministic_and_cover_small_inputs() {
        let small = vec![("chr1".to_string(), 1_000_000, 1_000), ("chr2".to_string(), 500_000, 500)];
        assert_eq!(sample_windows(&small, None, None), vec![("chr1".into(), 0, 1_000_000), ("chr2".into(), 0, 500_000)]);
        let big = vec![("chr1".to_string(), 200_000_000, 3_000_000), ("chr2".to_string(), 100_000_000, 1_000_000)];
        let w = sample_windows(&big, None, None);
        assert_eq!(w.len(), N_WINDOWS);
        assert_eq!(w.iter().filter(|x| x.0 == "chr1").count(), 48);
        assert_eq!(w, sample_windows(&big, None, None));
    }
}
