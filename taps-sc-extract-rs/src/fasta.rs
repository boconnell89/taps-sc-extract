//! Process-safe `.fai`-indexed FASTA reader (byte seek, no shared C handles).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
struct FaiRec {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

#[derive(Clone, Debug)]
pub struct FastFaiReader {
    fa_path: PathBuf,
    index: HashMap<String, FaiRec>,
}

impl FastFaiReader {
    pub fn open(fa_path: impl AsRef<Path>) -> std::io::Result<Self> {
        let fa_path = fa_path.as_ref().to_path_buf();
        let fai_path = {
            let mut p = fa_path.clone().into_os_string();
            p.push(".fai");
            PathBuf::from(p)
        };
        let file = File::open(&fai_path)?;
        let reader = BufReader::new(file);
        let mut index = HashMap::new();
        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 5 {
                continue;
            }
            index.insert(
                parts[0].to_string(),
                FaiRec {
                    length: parts[1].parse().unwrap_or(0),
                    offset: parts[2].parse().unwrap_or(0),
                    line_bases: parts[3].parse().unwrap_or(0),
                    line_width: parts[4].parse().unwrap_or(0),
                },
            );
        }
        Ok(Self { fa_path, index })
    }

    pub fn reference_length(&self, contig: &str) -> Option<u64> {
        self.index.get(contig).map(|r| r.length)
    }

    #[allow(dead_code)]
    pub fn contigs(&self) -> impl Iterator<Item = &str> {
        self.index.keys().map(|s| s.as_str())
    }

    /// Fetch `[start, end)` as uppercase ASCII. Coordinates are 0-based half-open.
    pub fn fetch(&self, contig: &str, start: u64, end: u64) -> std::io::Result<Vec<u8>> {
        let rec = self
            .index
            .get(contig)
            .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::NotFound, contig))?;
        let end = end.min(rec.length);
        if start >= end {
            return Ok(Vec::new());
        }
        if rec.line_bases == 0 || rec.line_width == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "invalid .fai line_bases/line_width",
            ));
        }
        let start_line = start / rec.line_bases;
        let start_col = start % rec.line_bases;
        let start_byte = rec.offset + start_line * rec.line_width + start_col;
        let end_pos = end - 1;
        let end_line = end_pos / rec.line_bases;
        let end_col = end_pos % rec.line_bases;
        let end_byte = rec.offset + end_line * rec.line_width + end_col + 1;

        let mut f = File::open(&self.fa_path)?;
        f.seek(SeekFrom::Start(start_byte))?;
        let mut raw = vec![0u8; (end_byte - start_byte) as usize];
        f.read_exact(&mut raw)?;
        let mut out = Vec::with_capacity((end - start) as usize);
        for &b in &raw {
            if b != b'\n' && b != b'\r' {
                out.push(b.to_ascii_uppercase());
            }
        }
        let want = (end - start) as usize;
        if out.len() > want {
            out.truncate(want);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn fai_roundtrip_small() {
        let dir = std::env::temp_dir();
        let fa = dir.join("taps_rs_tiny.fa");
        let fai = dir.join("taps_rs_tiny.fa.fai");
        // 10 bases per line, line width 11 including newline
        let mut f = File::create(&fa).unwrap();
        write!(f, ">chrT\nACGTACGTAC\nGGGGAAAAAT\n").unwrap();
        let mut i = File::create(&fai).unwrap();
        // name, length, offset, linebases, linewidth
        // header ">chrT\n" is 6 bytes; sequence starts at 6
        writeln!(i, "chrT\t20\t6\t10\t11").unwrap();

        let r = FastFaiReader::open(&fa).unwrap();
        assert_eq!(r.reference_length("chrT"), Some(20));
        let s = r.fetch("chrT", 0, 4).unwrap();
        assert_eq!(s, b"ACGT");
        let s = r.fetch("chrT", 8, 12).unwrap();
        assert_eq!(s, b"ACGG");
        let _ = std::fs::remove_file(&fa);
        let _ = std::fs::remove_file(&fai);
    }
}
