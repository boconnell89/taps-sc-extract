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

/// First whitespace-separated column; skips empty and `#` comments.
pub fn parse_annot(path: &Path) -> std::io::Result<HashSet<String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut barcodes = HashSet::new();
    for line in reader.lines() {
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
