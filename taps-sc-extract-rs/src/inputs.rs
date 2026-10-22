//! Several input BAMs in one run: `-b a.bam b.bam`, `-b @list.txt` (one path
//! per line) or `-b @sheet.tsv` (a header row with a `bam` column, and
//! optional per-file `complexity`, `min_complexity`, `max_complexity` and
//! `whitelist`). A barcode belongs to exactly one file. Same rules as Python
//! `taps_sc_extract/inputs.py::parse_bam_args`.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{Context as _, Result};

/// Per-file cell filter from a sample sheet row (blank = the global flag).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RowFilter {
    pub complexity: Option<PathBuf>,
    pub min_complexity: Option<f64>,
    pub max_complexity: Option<f64>,
    pub whitelist: Option<PathBuf>,
}

/// The run's input BAMs, in the order given, with their per-file filters.
#[derive(Clone, Debug)]
pub struct InputSpec {
    pub bams: Vec<PathBuf>,
    pub filters: Vec<RowFilter>,
    /// The list or sheet file, when inputs came from one (names the M-bias output).
    pub list: Option<PathBuf>,
}

impl InputSpec {
    pub fn is_multi(&self) -> bool {
        self.bams.len() > 1
    }

    /// Stem for outputs named after the input: the list file, else the first BAM.
    pub fn name_path(&self) -> &Path {
        self.list.as_deref().unwrap_or(&self.bams[0])
    }
}

const SHEET_COLUMNS: [&str; 5] = ["bam", "complexity", "min_complexity", "max_complexity", "whitelist"];

fn resolve(base: &Path, p: &str) -> PathBuf {
    let p = Path::new(p);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        base.join(p)
    }
}

/// `-b` values: paths, or a single `@file` naming a list or a sheet.
pub fn parse_bam_args(args: &[PathBuf]) -> Result<InputSpec> {
    let at: Vec<&PathBuf> = args.iter().filter(|a| a.to_string_lossy().starts_with('@')).collect();
    if at.is_empty() {
        if args.is_empty() {
            anyhow::bail!("no input BAM (-b)");
        }
        return Ok(InputSpec { bams: args.to_vec(), filters: vec![RowFilter::default(); args.len()], list: None });
    }
    if args.len() > 1 {
        anyhow::bail!("-b @FILE takes one list or sheet; do not mix it with other -b paths");
    }
    let list = PathBuf::from(&at[0].to_string_lossy()[1..]);
    if list.extension().is_some_and(|e| e.eq_ignore_ascii_case("xlsx") || e.eq_ignore_ascii_case("xls")) {
        anyhow::bail!("{}: Excel sheets are read by the taps-sc-extract wrapper; save as TSV for taps-sc-extract-rs", list.display());
    }
    let text = std::fs::read_to_string(&list).with_context(|| format!("read input list {}", list.display()))?;
    let base = list.parent().map(Path::to_path_buf).unwrap_or_default();
    let lines: Vec<(usize, &str)> = text
        .lines()
        .enumerate()
        .map(|(i, l)| (i + 1, l.trim_end_matches('\r')))
        .filter(|(_, l)| !l.trim().is_empty() && !l.trim_start().starts_with('#'))
        .collect();
    let Some(&(_, first)) = lines.first() else {
        anyhow::bail!("{}: no input BAMs listed", list.display());
    };
    let header: Vec<String> = first.split('\t').map(|c| c.trim().to_ascii_lowercase()).collect();
    let (bams, filters) = if header.iter().any(|c| c == "bam") {
        parse_sheet(&list, &base, &header, &lines[1..])?
    } else {
        let bams: Vec<PathBuf> = lines.iter().map(|(_, l)| resolve(&base, l.trim())).collect();
        let n = bams.len();
        (bams, vec![RowFilter::default(); n])
    };
    if bams.is_empty() {
        anyhow::bail!("{}: no input BAMs listed", list.display());
    }
    let mut seen = HashSet::new();
    for b in &bams {
        if !seen.insert(b.clone()) {
            anyhow::bail!("{}: {} is listed twice", list.display(), b.display());
        }
    }
    Ok(InputSpec { bams, filters, list: Some(list) })
}

fn parse_sheet(
    list: &Path,
    base: &Path,
    header: &[String],
    rows: &[(usize, &str)],
) -> Result<(Vec<PathBuf>, Vec<RowFilter>)> {
    for c in header {
        if !c.is_empty() && !SHEET_COLUMNS.contains(&c.as_str()) {
            anyhow::bail!("{}: unknown column {c:?} (columns: {})", list.display(), SHEET_COLUMNS.join(", "));
        }
    }
    let col = |name: &str| header.iter().position(|c| c == name);
    let (cb, cc, cmin, cmax, cw) =
        (col("bam").expect("checked"), col("complexity"), col("min_complexity"), col("max_complexity"), col("whitelist"));
    let mut bams = Vec::new();
    let mut filters = Vec::new();
    for &(line_no, line) in rows {
        let f: Vec<&str> = line.split('\t').map(str::trim).collect();
        let get = |i: Option<usize>| i.and_then(|i| f.get(i).copied()).filter(|v| !v.is_empty());
        let num = |i: Option<usize>, name: &str| -> Result<Option<f64>> {
            get(i)
                .map(|v| v.parse::<f64>().with_context(|| format!("{} line {line_no}: {name} {v:?} is not a number", list.display())))
                .transpose()
        };
        let Some(bam) = get(Some(cb)) else {
            anyhow::bail!("{} line {line_no}: empty bam", list.display());
        };
        let row = RowFilter {
            complexity: get(cc).map(|p| resolve(base, p)),
            min_complexity: num(cmin, "min_complexity")?,
            max_complexity: num(cmax, "max_complexity")?,
            whitelist: get(cw).map(|p| resolve(base, p)),
        };
        bams.push(resolve(base, bam));
        filters.push(row);
    }
    Ok((bams, filters))
}

/// Kept cells for the run: the union of every file's own filter (its row's
/// values, else the global flags). `None` keeps every barcode. Errors when
/// only some files have a filter, or when a barcode is in two files'
/// complexity files.
pub struct CellFilter {
    pub keep: Option<HashSet<String>>,
    /// Per file: (kept cells, barcodes in its complexity file), for the config printout.
    pub per_file: Vec<(Option<usize>, Option<usize>)>,
    /// Per file: every barcode in its complexity file (for the read sample check).
    pub listed: Vec<Option<HashSet<String>>>,
}

pub fn load_cell_filters(
    spec: &InputSpec,
    whitelist: Option<&Path>,
    complexity: Option<&Path>,
    min: Option<f64>,
    max: Option<f64>,
) -> Result<CellFilter> {
    if !spec.is_multi() && spec.filters[0] == RowFilter::default() {
        let keep = crate::barcode::load_cell_filter(whitelist, complexity, min, max)?;
        let listed = match complexity {
            Some(p) => Some(crate::barcode::parse_complexity(p, None, None)?.0),
            None => None,
        };
        return Ok(CellFilter {
            per_file: vec![(keep.as_ref().map(|k| k.len()), listed.as_ref().map(|l| l.len()))],
            keep,
            listed: vec![listed],
        });
    }
    let mut union: HashSet<String> = HashSet::new();
    let mut per_file = Vec::new();
    let mut listed: Vec<Option<HashSet<String>>> = Vec::new();
    let mut any_none = false;
    let mut any_some = false;
    let mut owner: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for (i, row) in spec.filters.iter().enumerate() {
        let cx = row.complexity.as_deref().or(complexity);
        let wl = row.whitelist.as_deref().or(whitelist);
        let lo = row.min_complexity.or(min);
        let hi = row.max_complexity.or(max);
        let keep = crate::barcode::load_cell_filter(wl, cx, lo, hi)
            .with_context(|| format!("cell filter for {}", spec.bams[i].display()))?;
        let all = match cx {
            Some(p) => Some(crate::barcode::parse_complexity(p, None, None)?.0),
            None => None,
        };
        // A shared global complexity file lists every file's barcodes; only
        // per-file complexity files are checked for overlap.
        if let (Some(all), Some(_)) = (&all, &row.complexity) {
            for bc in all {
                if let Some(&j) = owner.get(bc) {
                    anyhow::bail!(
                        "barcode {bc} is in the complexity files of both {} and {}; a barcode must belong to one input",
                        spec.bams[j].display(),
                        spec.bams[i].display()
                    );
                }
                owner.insert(bc.clone(), i);
            }
        }
        per_file.push((keep.as_ref().map(|k| k.len()), all.as_ref().map(|a| a.len())));
        match keep {
            Some(k) => {
                any_some = true;
                union.extend(k);
            }
            None => any_none = true,
        }
        listed.push(if row.complexity.is_some() { all } else { None });
    }
    if any_some && any_none {
        anyhow::bail!("some inputs have a cell filter (complexity/whitelist) and some do not; give every input one, or none");
    }
    Ok(CellFilter { keep: any_some.then_some(union), per_file, listed })
}

/// Every input must have the same contigs, in the same order, with the same lengths.
pub fn check_same_contigs(bams: &[PathBuf]) -> Result<()> {
    use rust_htslib::bam::Read as _;
    let sq = |p: &Path| -> Result<Vec<(Vec<u8>, u64)>> {
        let r = rust_htslib::bam::Reader::from_path(p).map_err(|e| anyhow::anyhow!("open {}: {e}", p.display()))?;
        let h = r.header();
        Ok((0..h.target_count()).map(|t| (h.tid2name(t).to_vec(), h.target_len(t).unwrap_or(0))).collect())
    };
    let first = sq(&bams[0])?;
    for b in &bams[1..] {
        if sq(b)? != first {
            anyhow::bail!(
                "{} and {} have different @SQ contigs or lengths; inputs must share one reference",
                bams[0].display(),
                b.display()
            );
        }
    }
    Ok(())
}

/// Reads sampled per input for `check_barcodes_match`.
const BARCODE_SAMPLE: usize = 5_000;

/// Catch a sheet row that points at another input's complexity file: in the
/// first `BARCODE_SAMPLE` mapped primary reads of each input with its own
/// complexity file, most barcodes must be listed in it.
pub fn check_barcodes_match(bams: &[PathBuf], listed: &[Option<HashSet<String>>]) -> Result<()> {
    use rust_htslib::bam::Read as _;
    if bams.len() < 2 {
        return Ok(());
    }
    for (b, set) in bams.iter().zip(listed) {
        let Some(set) = set else { continue };
        let mut r = rust_htslib::bam::Reader::from_path(b).map_err(|e| anyhow::anyhow!("open {}: {e}", b.display()))?;
        let (mut n, mut hit) = (0usize, 0usize);
        let mut rec = rust_htslib::bam::Record::new();
        while n < BARCODE_SAMPLE {
            match r.read(&mut rec) {
                Some(Ok(())) => {}
                Some(Err(e)) => anyhow::bail!("read {}: {e}", b.display()),
                None => break,
            }
            if rec.flags() & 0x904 != 0 {
                continue;
            }
            n += 1;
            let bc = std::str::from_utf8(crate::barcode::barcode_from_qname_bytes(rec.qname())).unwrap_or("");
            hit += usize::from(set.contains(bc));
        }
        if n > 0 && hit * 2 < n {
            anyhow::bail!(
                "only {hit} of {n} sampled reads in {} have a barcode in its complexity file; is the sheet row pointing at the right file?",
                b.display()
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("taps_inputs_{}_{name}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn plain_paths_list_and_sheet() {
        let s = parse_bam_args(&[PathBuf::from("a.bam"), PathBuf::from("b.bam")]).unwrap();
        assert_eq!(s.bams.len(), 2);
        assert!(s.list.is_none());

        let d = tmp("list");
        std::fs::write(d.join("l.txt"), "# inputs\nx.bam\n\n/abs/y.bam\n").unwrap();
        let s = parse_bam_args(&[PathBuf::from(format!("@{}", d.join("l.txt").display()))]).unwrap();
        assert_eq!(s.bams, vec![d.join("x.bam"), PathBuf::from("/abs/y.bam")]);

        std::fs::write(
            d.join("s.tsv"),
            "bam\tcomplexity\tmin_complexity\nx.bam\tx.cx.txt\t10000\ny.bam\t\t\n",
        )
        .unwrap();
        let s = parse_bam_args(&[PathBuf::from(format!("@{}", d.join("s.tsv").display()))]).unwrap();
        assert_eq!(s.bams, vec![d.join("x.bam"), d.join("y.bam")]);
        assert_eq!(s.filters[0].complexity, Some(d.join("x.cx.txt")));
        assert_eq!(s.filters[0].min_complexity, Some(10000.0));
        assert_eq!(s.filters[1], RowFilter::default());

        std::fs::write(d.join("bad.tsv"), "bam\tcutoff\nx.bam\t5\n").unwrap();
        assert!(parse_bam_args(&[PathBuf::from(format!("@{}", d.join("bad.tsv").display()))]).is_err());
        std::fs::write(d.join("dup.txt"), "x.bam\nx.bam\n").unwrap();
        assert!(parse_bam_args(&[PathBuf::from(format!("@{}", d.join("dup.txt").display()))]).is_err());
        assert!(parse_bam_args(&[PathBuf::from("@l.txt"), PathBuf::from("a.bam")]).is_err());
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn per_file_filters_union_and_overlap() {
        let d = tmp("filters");
        std::fs::write(d.join("a.txt"), "1\tA1-x\t100\t50\t50\n2\tA2-x\t100\t5\t5\n").unwrap();
        std::fs::write(d.join("b.txt"), "1\tB1-y\t100\t50\t50\n2\tB2-y\t100\t20\t20\n").unwrap();
        std::fs::write(
            d.join("s.tsv"),
            "bam\tcomplexity\tmin_complexity\na.bam\ta.txt\t10\nb.bam\tb.txt\t30\n",
        )
        .unwrap();
        let s = parse_bam_args(&[PathBuf::from(format!("@{}", d.join("s.tsv").display()))]).unwrap();
        let f = load_cell_filters(&s, None, None, None, None).unwrap();
        let keep = f.keep.unwrap();
        assert_eq!(keep, ["A1-x", "B1-y"].iter().map(|s| s.to_string()).collect());
        assert_eq!(f.per_file, vec![(Some(1), Some(2)), (Some(1), Some(2))]);

        // A barcode in two files' complexity files.
        std::fs::write(d.join("c.txt"), "1\tA1-x\t100\t50\t50\n").unwrap();
        std::fs::write(d.join("s2.tsv"), "bam\tcomplexity\na.bam\ta.txt\nc.bam\tc.txt\n").unwrap();
        let s = parse_bam_args(&[PathBuf::from(format!("@{}", d.join("s2.tsv").display()))]).unwrap();
        assert!(load_cell_filters(&s, None, None, None, None).is_err());

        // Only one of two files filtered.
        std::fs::write(d.join("s3.tsv"), "bam\tcomplexity\na.bam\ta.txt\nb.bam\t\n").unwrap();
        let s = parse_bam_args(&[PathBuf::from(format!("@{}", d.join("s3.tsv").display()))]).unwrap();
        assert!(load_cell_filters(&s, None, None, None, None).is_err());
        std::fs::remove_dir_all(&d).ok();
    }
}
