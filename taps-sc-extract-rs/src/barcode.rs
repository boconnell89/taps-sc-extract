//! QNAME barcode extraction, whitelist parsing, and MD5 shard assignment.

use md5::{Digest, Md5};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// First field of QNAME as bytes: split on `:`, else `/`, else the whole name.
pub fn barcode_from_qname_bytes(qname: &[u8]) -> &[u8] {
    if let Some(i) = qname.iter().position(|&b| b == b':') {
        &qname[..i]
    } else if let Some(i) = qname.iter().position(|&b| b == b'/') {
        &qname[..i]
    } else {
        qname
    }
}

/// First field of QNAME: split on `:`, else `/`, else the whole name.
pub fn extract_barcode(qname: &str) -> &str {
    let b = barcode_from_qname_bytes(qname.as_bytes());
    std::str::from_utf8(b).unwrap_or(qname)
}

/// Lines of a plain or gzipped (`.gz`) text file.
pub fn open_text(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("gz")) {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// First whitespace-separated column; skips empty and `#` comments.
pub fn parse_annot(path: &Path) -> std::io::Result<HashSet<String>> {
    let mut barcodes = HashSet::new();
    for line in open_text(path)?.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(bc) = line.split_whitespace().next() {
            barcodes.insert(bc.to_string());
        }
    }
    Ok(barcodes)
}

/// Barcodes from a complexity file whose unique-read count is within
/// `[min, max]` (inclusive; `None` = unbounded).
///
/// Format (whitespace-separated, as written by the sciMET complexity step):
/// `rank  barcode  total_reads  unique_reads  percent_unique`. Column 2 is the
/// barcode and column 4 the unique-read count. Empty lines, `#` comments and a
/// header line whose 4th column is not a number are skipped. Returns the
/// passing barcodes and the number of barcodes in the file.
pub fn parse_complexity(
    path: &Path,
    min: Option<f64>,
    max: Option<f64>,
) -> anyhow::Result<(HashSet<String>, usize)> {
    let mut pass = HashSet::new();
    let mut total = 0usize;
    for (i, line) in open_text(path)?.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            anyhow::bail!(
                "{} line {}: expected at least 4 columns (rank, barcode, total, unique), got {}",
                path.display(),
                i + 1,
                f.len()
            );
        }
        let Ok(unique) = f[3].parse::<f64>() else {
            if total == 0 {
                continue; // header
            }
            anyhow::bail!("{} line {}: unique-read count {:?} is not a number", path.display(), i + 1, f[3]);
        };
        total += 1;
        if min.is_some_and(|m| unique < m) || max.is_some_and(|m| unique > m) {
            continue;
        }
        pass.insert(f[1].to_string());
    }
    Ok((pass, total))
}

/// The set of barcodes to keep, from `--whitelist` and/or a complexity file
/// with `--min-complexity` / `--max-complexity` (intersected when both are
/// given). `None` keeps every barcode.
pub fn load_cell_filter(
    whitelist: Option<&Path>,
    complexity: Option<&Path>,
    min: Option<f64>,
    max: Option<f64>,
) -> anyhow::Result<Option<HashSet<String>>> {
    if complexity.is_none() && (min.is_some() || max.is_some()) {
        anyhow::bail!("--min-complexity / --max-complexity need --complexity FILE");
    }
    let wl = match whitelist {
        Some(p) => {
            let w = parse_annot(p).map_err(|e| anyhow::anyhow!("read whitelist {}: {e}", p.display()))?;
            eprintln!("Whitelist: {} barcodes from {}", w.len(), p.display());
            Some(w)
        }
        None => None,
    };
    let Some(cp) = complexity else {
        return Ok(wl);
    };
    let (pass, total) = parse_complexity(cp, min, max)?;
    let fmt = |v: Option<f64>| v.map_or_else(|| "-".to_string(), |x| x.to_string());
    eprintln!(
        "Complexity filter: {} of {} barcodes with unique reads in [{}, {}] ({})",
        pass.len(),
        total,
        fmt(min),
        fmt(max),
        cp.display()
    );
    Ok(Some(match wl {
        Some(w) => {
            let both: HashSet<String> = pass.intersection(&w).cloned().collect();
            eprintln!("Whitelist and complexity filter together: {} barcodes", both.len());
            both
        }
        None => pass,
    }))
}

/// Same mapping as Python `_barcode_to_shard`: `md5(ascii) % n_shards`.
pub fn barcode_to_shard(bc: &str, n_shards: usize) -> usize {
    if n_shards <= 1 {
        return 0;
    }
    let mut hasher = Md5::new();
    hasher.update(bc.as_bytes());
    let digest = hasher.finalize();
    let mut n: u128 = 0;
    for b in digest {
        n = (n << 8) | u128::from(b);
    }
    (n % n_shards as u128) as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_tmp(name: &str, body: &[u8]) -> std::path::PathBuf {
        let dir = tempfile::tempdir().unwrap().keep();
        let p = dir.join(name);
        std::fs::write(&p, body).unwrap();
        p
    }

    const COMPLEXITY: &str = "rank\tbarcode\ttotal\tunique\tpct\n\
        1\tAAA_I1\t1000\t900\t90.0\n\
        2\tCCC_I1\t500\t400\t80.0\n\
        # comment\n\
        \n\
        3\tGGG_I1\t100\t50\t50.0\n";

    #[test]
    fn complexity_bounds_inclusive_and_header_skipped() {
        let p = write_tmp("c.txt", COMPLEXITY.as_bytes());
        let (all, n) = parse_complexity(&p, None, None).unwrap();
        assert_eq!(n, 3);
        assert_eq!(all.len(), 3);
        let (mid, _) = parse_complexity(&p, Some(400.0), Some(900.0)).unwrap();
        assert!(mid.contains("AAA_I1") && mid.contains("CCC_I1") && !mid.contains("GGG_I1"));
        let (lo, _) = parse_complexity(&p, None, Some(399.0)).unwrap();
        assert_eq!(lo.into_iter().collect::<Vec<_>>(), vec!["GGG_I1".to_string()]);
    }

    #[test]
    fn complexity_rejects_short_rows_and_bad_numbers() {
        let p = write_tmp("short.txt", b"1\tAAA\t10\n");
        assert!(parse_complexity(&p, None, None).is_err());
        let p = write_tmp("bad.txt", b"1\tAAA\t10\t5\t50\n2\tCCC\t10\tx\t50\n");
        assert!(parse_complexity(&p, None, None).is_err());
    }

    #[test]
    fn complexity_gz_and_whitelist_intersection() {
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(COMPLEXITY.as_bytes()).unwrap();
        let cp = write_tmp("c.txt.gz", &gz.finish().unwrap());
        let wl = write_tmp("wl.txt", b"CCC_I1 extra\nGGG_I1\nTTT_I1\n");
        let f = load_cell_filter(Some(&wl), Some(&cp), Some(100.0), None).unwrap().unwrap();
        assert_eq!(f.into_iter().collect::<Vec<_>>(), vec!["CCC_I1".to_string()]);
        assert!(load_cell_filter(None, None, Some(1.0), None).is_err());
        assert!(load_cell_filter(None, None, None, None).unwrap().is_none());
    }

    #[test]
    fn qname_formats() {
        assert_eq!(
            extract_barcode("ATTGGCTCAATGATCCTGGCGCTCCATTCG:77600300#0"),
            "ATTGGCTCAATGATCCTGGCGCTCCATTCG"
        );
        assert_eq!(extract_barcode("ACGT:123#0"), "ACGT");
        assert_eq!(extract_barcode("BC123/1"), "BC123");
        assert_eq!(extract_barcode("BC123"), "BC123");
    }

    #[test]
    fn annot_first_column() {
        let dir = std::env::temp_dir();
        let p = dir.join("taps_rs_annot_test.txt");
        let mut f = File::create(&p).unwrap();
        writeln!(f, "BC1\tsampleA\n# comment\n\nBC2\tsampleB\nBC3").unwrap();
        let set = parse_annot(&p).unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains("BC1") && set.contains("BC3"));
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn shard_deterministic() {
        assert_eq!(barcode_to_shard("anything", 1), 0);
        let a = barcode_to_shard("cellA", 16);
        let b = barcode_to_shard("cellA", 16);
        assert_eq!(a, b);
        assert!(a < 16);
        // Matches Python hashlib.md5(...).hexdigest() int % n
        assert_eq!(barcode_to_shard("cellA", 16), 7);
        assert_eq!(barcode_to_shard("x", 8), 6);
    }
}
