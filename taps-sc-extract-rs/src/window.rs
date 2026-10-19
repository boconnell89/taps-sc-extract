//! Genomic window planner matching Python `plan_genomic_chunks`.

use crate::fasta::FastFaiReader;

/// Rank of a primary chromosome name, or `None` for anything else.
///
/// Primary: an optional `chr`/`Chr`/`CHR` prefix followed by a number
/// (autosomes) or one of X, Y, Z, W (sex chromosomes). That covers UCSC
/// (`chr1`, `chrX`) and Ensembl/NCBI (`1`, `X`) naming for mammals and birds,
/// and leaves out mitochondria (chrM, MT), unplaced and unlocalised scaffolds
/// (`chrUn_...`, `chr1_..._random`), alt haplotypes, decoys and spike-ins.
/// Autosomes sort numerically, then X, Y, Z, W.
fn primary_rank(name: &str) -> Option<(u8, u32)> {
    let core = name
        .strip_prefix("chr")
        .or_else(|| name.strip_prefix("Chr"))
        .or_else(|| name.strip_prefix("CHR"))
        .unwrap_or(name);
    if !core.is_empty() && core.len() <= 3 && core.bytes().all(|b| b.is_ascii_digit()) && !core.starts_with('0') {
        return Some((0, core.parse().ok()?));
    }
    match core {
        "X" => Some((1, 0)),
        "Y" => Some((1, 1)),
        "Z" => Some((1, 2)),
        "W" => Some((1, 3)),
        _ => None,
    }
}

/// Default contigs for a BAM: its primary chromosomes (`primary_rank`), in
/// natural order (1, 2, ..., 10, ..., X, Y), whatever order the header uses.
/// mm10 gives chr1..chr19, chrX, chrY; hg38 gives chr1..chr22, chrX, chrY.
pub fn primary_contigs(names: &[String]) -> Vec<String> {
    let mut v: Vec<((u8, u32), &String)> = names.iter().filter_map(|n| primary_rank(n).map(|r| (r, n))).collect();
    v.sort_by_key(|&(r, _)| r);
    v.dedup_by_key(|&mut (r, _)| r);
    v.into_iter().map(|(_, n)| n.clone()).collect()
}

/// One-line summary of `primary_contigs`, e.g. "22 autosomes + chrX, chrY".
pub fn describe_primary(contigs: &[String]) -> String {
    let auto = contigs.iter().filter(|c| matches!(primary_rank(c), Some((0, _)))).count();
    let sex: Vec<&str> = contigs.iter().filter(|c| matches!(primary_rank(c), Some((1, _)))).map(|c| c.as_str()).collect();
    if sex.is_empty() {
        format!("{auto} autosomes")
    } else {
        format!("{auto} autosomes + {}", sex.join(", "))
    }
}

/// Primary chromosomes from a BAM/CRAM header (see `primary_contigs`).
pub fn bam_primary_contigs(bam: &std::path::Path) -> anyhow::Result<Vec<String>> {
    use rust_htslib::bam::Read;
    let r = crate::input::open_indexed(bam)?;
    let names: Vec<String> =
        r.header().target_names().iter().map(|n| String::from_utf8_lossy(n).into_owned()).collect();
    let v = primary_contigs(&names);
    anyhow::ensure!(
        !v.is_empty(),
        "no primary chromosomes (1, 2, ..., X, Y with or without a chr prefix) in the header of {}; pass --chroms",
        bam.display()
    );
    Ok(v)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Window {
    pub chunk_id: u32,
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

pub fn plan_windows(
    fai: &FastFaiReader,
    contigs: &[&str],
    chunk_size_bp: u64,
) -> Vec<Window> {
    let chunk_size_bp = chunk_size_bp.max(1);
    let mut out = Vec::new();
    let mut chunk_id = 0u32;
    for &contig in contigs {
        let Some(len) = fai.reference_length(contig) else {
            continue;
        };
        let mut start = 0u64;
        while start < len {
            let end = (start + chunk_size_bp).min(len);
            out.push(Window {
                chunk_id,
                contig: contig.to_string(),
                start,
                end,
            });
            chunk_id += 1;
            start = end;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn primary_contigs_mm10_hg38_ensembl_birds() {
        // mm10 header in lexicographic order, with scaffolds and chrM.
        let mut mm10: Vec<String> = (1..=19).map(|i| format!("chr{i}")).collect();
        mm10.sort();
        mm10.extend(names(&["chrM", "chrX", "chrY", "chr1_GL456210_random", "chrUn_GL456239", "chr4_JH584295_random"]));
        let want: Vec<String> = (1..=19).map(|i| format!("chr{i}")).chain(names(&["chrX", "chrY"])).collect();
        assert_eq!(primary_contigs(&mm10), want);
        assert_eq!(describe_primary(&want), "19 autosomes + chrX, chrY");
        // hg38 with alts, decoys and EBV.
        let mut hg38: Vec<String> = (1..=22).map(|i| format!("chr{i}")).collect();
        hg38.extend(names(&["chrX", "chrY", "chrM", "chr1_KI270706v1_random", "chr6_GL000250v2_alt", "chrEBV", "chrUn_KI270302v1", "HLA-A*01:01:01:01"]));
        let got = primary_contigs(&hg38);
        assert_eq!(got.len(), 24);
        assert_eq!(got[21], "chr22");
        assert_eq!(describe_primary(&got), "22 autosomes + chrX, chrY");
        // Ensembl / NCBI naming without the prefix.
        let ens = names(&["1", "2", "10", "X", "Y", "MT", "KI270728.1", "GL000009.2"]);
        assert_eq!(primary_contigs(&ens), names(&["1", "2", "10", "X", "Y"]));
        // Birds (ZW).
        assert_eq!(primary_contigs(&names(&["chrZ", "chr2", "chrW", "chr1", "chrM"])), names(&["chr1", "chr2", "chrZ", "chrW"]));
        assert!(primary_contigs(&names(&["chrM", "lambda", "pUC19"])).is_empty());
    }

    #[test]
    fn plans_non_overlapping_windows() {
        use crate::fasta::FastFaiReader;
        use std::io::Write;
        let dir = std::env::temp_dir();
        let fa = dir.join("taps_rs_win.fa");
        let fai = dir.join("taps_rs_win.fa.fai");
        let mut f = std::fs::File::create(&fa).unwrap();
        write!(f, ">chr1\n").unwrap();
        for _ in 0..5 {
            writeln!(f, "AAAAAAAAAA").unwrap();
        }
        let mut i = std::fs::File::create(&fai).unwrap();
        writeln!(i, "chr1\t50\t6\t10\t11").unwrap();
        let r = FastFaiReader::open(&fa).unwrap();
        let w = plan_windows(&r, &["chr1"], 20);
        assert_eq!(w.len(), 3);
        assert_eq!(w[0], Window { chunk_id: 0, contig: "chr1".into(), start: 0, end: 20 });
        assert_eq!(w[2].end, 50);
        let _ = std::fs::remove_file(&fa);
        let _ = std::fs::remove_file(&fai);
    }
}
