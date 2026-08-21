//! Genomic window planner matching Python `plan_genomic_chunks`.

use crate::fasta::FastFaiReader;

pub const CANONICAL_CONTIGS: &[&str] = &[
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY",
];

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

    #[test]
    fn canonical_count() {
        assert_eq!(CANONICAL_CONTIGS.len(), 21);
        assert!(!CANONICAL_CONTIGS.contains(&"chrM"));
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
