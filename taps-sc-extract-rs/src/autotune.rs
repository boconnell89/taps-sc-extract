//! Conservative parameter selection. Explicit CLI values always win.

/// Auto window sizing: about this many windows per worker, so rayon work
/// stealing has enough units to keep every worker busy to the end.
pub const AUTO_WINDOWS_PER_WORKER: u64 = 6;
pub const MIN_AUTO_CHUNK_BP: u64 = 2_000_000;
pub const MAX_AUTO_CHUNK_BP: u64 = 100_000_000;
/// In-process shard writers share hdf5-metno's global lock, so extra threads only
/// add resident shard data. Measured 2026-10-14 (mm10, 16 shards, LZF): 6 writers
/// wrote in 36 s at 33 GB peak RSS, 16 writers in 86 s at 47 GB.
pub const MAX_AUTO_WRITERS: usize = 6;
pub const DEFAULT_ESTIMATED_CELLS: usize = 10_000;

/// Read /proc/meminfo MemAvailable in GB on Linux, or fall back to 32 GB.
pub fn get_available_memory_gb() -> f64 {
    if let Ok(content) = std::fs::read_to_string("/proc/meminfo") {
        for line in content.lines() {
            if line.starts_with("MemAvailable:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    if let Ok(kb) = parts[1].parse::<f64>() {
                        return kb / (1024.0 * 1024.0);
                    }
                }
            }
        }
    }
    32.0
}

/// Default share of MemAvailable the run may use (`--max-memory-fraction`).
pub const DEFAULT_MEMORY_FRACTION: f64 = 0.6;

/// Memory budget for the whole process tree (extraction plus writer children):
/// `--max-memory-gb` if given, else `fraction` (default 0.6) of MemAvailable.
pub fn system_memory_budget_gb(max_memory_gb: Option<f64>, fraction: Option<f64>) -> (f64, String) {
    if let Some(gb) = max_memory_gb {
        (gb.max(1.0), "explicit".to_string())
    } else {
        let f = fraction.unwrap_or(DEFAULT_MEMORY_FRACTION).clamp(0.05, 1.0);
        let avail = get_available_memory_gb();
        ((avail * f).max(2.0), format!("{f:.2} x MemAvailable {avail:.1} GB"))
    }
}

/// BGZF decompression threads per worker when the user does not pass `--decomp-threads`.
/// More threads when fewer extraction workers. At 8 or more workers the cores are
/// already busy with pileup, and 0 (inflate on the worker thread) measured fastest.
pub fn default_decomp_threads(n_workers: usize) -> usize {
    if n_workers <= 4 {
        4
    } else if n_workers <= 7 {
        2
    } else {
        0
    }
}

/// Window size when `--chunk-size-mb` is 0: `selected_bp / (6 × workers)`,
/// clamped to [2, 100] Mb. mm10 canonical gives ~28 Mb at 16 workers and ~4.7 Mb at 96.
pub fn auto_chunk_size_bp(selected_bp: u64, n_workers: usize) -> u64 {
    let target = (n_workers.max(1) as u64) * AUTO_WINDOWS_PER_WORKER;
    selected_bp
        .div_ceil(target)
        .clamp(MIN_AUTO_CHUNK_BP, MAX_AUTO_CHUNK_BP)
}

fn proc_status_kb(key: &str) -> Option<f64> {
    let content = std::fs::read_to_string("/proc/self/status").ok()?;
    for line in content.lines() {
        if let Some(rest) = line.strip_prefix(key) {
            let kb = rest.split_whitespace().next()?;
            return kb.parse().ok();
        }
    }
    None
}

/// Current resident set of this process, in megabytes (`VmRSS`).
pub fn process_rss_mb() -> f64 {
    proc_status_kb("VmRSS:").unwrap_or(0.0) / 1024.0
}

/// Peak resident set of this process, in megabytes (`VmHWM`).
pub fn process_peak_rss_mb() -> f64 {
    match proc_status_kb("VmHWM:") {
        Some(kb) => kb / 1024.0,
        None => process_rss_mb(),
    }
}

/// Extraction-process memory per worker in stream mode (window maps, BAM
/// buffers); measured ~4.9 GB at 32 workers on mm10 (2026-10-15).
pub const STREAM_GB_PER_WORKER: f64 = 0.15;
pub const STREAM_BASE_GB: f64 = 1.0;

/// Window workers when `--workers` is omitted: min(nproc, workers that fit the
/// stream-mode extraction footprint in the budget).
pub fn default_workers(budget_gb: f64) -> usize {
    let nproc = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mem_workers = ((budget_gb - STREAM_BASE_GB).max(STREAM_GB_PER_WORKER) / STREAM_GB_PER_WORKER).floor() as usize;
    nproc.min(mem_workers).max(1)
}

/// Resolve worker count (0 means auto).
pub fn resolve_workers_with_budget(cli: usize, budget_gb: f64) -> (usize, &'static str) {
    if cli == 0 {
        (default_workers(budget_gb), "auto")
    } else {
        (cli, "explicit")
    }
}

/// Deprecated resolve_workers without budget for backwards compatibility.
pub fn resolve_workers(cli: usize) -> usize {
    if cli == 0 {
        default_workers(get_available_memory_gb() * 0.6)
    } else {
        cli
    }
}

/// Resolve shard count (0 means auto based on expected cell count).
pub fn resolve_shards(cli_shards: usize, cell_count: usize) -> (usize, &'static str) {
    if cli_shards > 0 {
        (cli_shards, "explicit")
    } else {
        let s = if cell_count <= 10_000 {
            1
        } else if cell_count <= 50_000 {
            8
        } else if cell_count <= 100_000 {
            16
        } else {
            32
        };
        (s, "auto")
    }
}

/// Memory model, measured 2026-10-15 on mm10 (53.5M reads, 7,351 cells, 702M sites):
/// 13.1 sites per mapped read; memory-mode extraction peak ~16.4 B/site; the 16
/// sliced `write-shard` children together ~14.5 B/site (10.2 GB). Margins added.
pub const SITES_PER_READ: f64 = 14.0;
pub const MEMORY_MODE_BYTES_PER_SITE: f64 = 18.0;
pub const WRITER_BYTES_PER_SITE: f64 = 20.0;
pub const WRITER_BASE_GB: f64 = 0.1;

/// Memory mode and writer count, chosen together so the whole process tree
/// (extraction process + concurrent writer children) fits the budget.
#[derive(Clone, Debug)]
pub struct MemoryPlan {
    pub mode: String,
    pub mode_src: &'static str,
    pub writers: usize,
    pub writers_src: &'static str,
    pub est_sites: f64,
    /// Extraction-process peak (GB) in the chosen mode.
    pub extract_gb: f64,
    /// Extraction process while writers run (memory mode still holds payloads).
    pub parent_during_write_gb: f64,
    pub writer_each_gb: f64,
    pub peak_gb: f64,
}

pub struct MemoryInputs<'a> {
    pub mode: &'a str,
    pub budget_gb: f64,
    /// Mapped reads on the selected contigs (BAM index); 0 if unknown.
    pub mapped_reads: u64,
    /// Share of sites in the extracted contexts (`ContextSel::site_share`).
    pub site_share: f64,
    /// Fallback when `mapped_reads` is 0.
    pub cell_count: usize,
    pub genome_mb: u64,
    pub n_shards: usize,
    pub n_workers: usize,
    /// `--max-writer-threads` (0 = auto).
    pub cli_writers: usize,
    /// Default from a saved machine profile, still bounded by memory.
    pub profile_writers: Option<usize>,
    pub in_process: bool,
}

pub fn plan_memory(i: &MemoryInputs) -> MemoryPlan {
    let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let sites = if i.mapped_reads > 0 {
        i.mapped_reads as f64 * SITES_PER_READ * i.site_share
    } else {
        // No index stats: the old per-cell guess.
        i.cell_count as f64 * (i.genome_mb * 500).min(200_000) as f64
    };
    let n_shards = i.n_shards.max(1);
    let stream_gb = STREAM_BASE_GB + STREAM_GB_PER_WORKER * i.n_workers as f64;
    let payload_gb = sites * MEMORY_MODE_BYTES_PER_SITE / 1e9;
    let memory_gb = stream_gb + payload_gb;
    let writer_each_gb = sites / n_shards as f64 * WRITER_BYTES_PER_SITE / 1e9 + WRITER_BASE_GB;
    let writer_cap = n_shards
        .min(nproc)
        .min(if i.in_process { MAX_AUTO_WRITERS } else { usize::MAX })
        .min(i.profile_writers.unwrap_or(usize::MAX))
        .max(1);
    // In-process writers share the extraction process: its footprint is the
    // parent's, already counted, plus one shard each.
    let fit = |parent_gb: f64| -> usize {
        let room = i.budget_gb - parent_gb;
        ((room / writer_each_gb).floor().max(1.0) as usize).min(writer_cap)
    };
    let (mode, mode_src) = match i.mode.to_lowercase().as_str() {
        "stream" => ("stream".to_string(), "explicit"),
        "memory" => ("memory".to_string(), "explicit"),
        _ => {
            // Memory mode only if it fits and does not cost write parallelism.
            let fits = memory_gb + writer_each_gb <= i.budget_gb;
            if fits && fit(payload_gb + STREAM_BASE_GB) >= fit(STREAM_BASE_GB) {
                ("memory".to_string(), "auto (fits in RAM budget)")
            } else {
                ("stream".to_string(), "auto (spill to temp disk)")
            }
        }
    };
    let (extract_gb, parent_during_write_gb) = if mode == "memory" {
        (memory_gb, payload_gb + STREAM_BASE_GB)
    } else {
        (stream_gb, STREAM_BASE_GB)
    };
    let (writers, writers_src) = if i.cli_writers > 0 {
        (i.cli_writers, "explicit")
    } else if i.profile_writers.is_some() {
        (fit(parent_during_write_gb), "machine profile, memory-bounded")
    } else {
        (fit(parent_during_write_gb), "auto, memory-bounded")
    };
    let peak_gb = extract_gb.max(parent_during_write_gb + writers as f64 * writer_each_gb);
    MemoryPlan {
        mode,
        mode_src,
        writers,
        writers_src,
        est_sites: sites,
        extract_gb,
        parent_during_write_gb,
        writer_each_gb,
        peak_gb,
    }
}

/// Cell-count source for auto-tune.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CellCountSource {
    Explicit(usize),
    Whitelist(usize),
    Estimate(usize),
}

impl CellCountSource {
    pub fn count(&self) -> usize {
        match *self {
            Self::Explicit(n) | Self::Whitelist(n) | Self::Estimate(n) => n,
        }
    }

    pub fn label(&self) -> &'static str {
        match *self {
            Self::Explicit(_) => "explicit (--expected-cells)",
            Self::Whitelist(_) => "whitelist cardinality",
            Self::Estimate(_) => "estimated default",
        }
    }

    /// Precedence: explicit `--expected-cells` > whitelist cardinality > estimate.
    pub fn resolve(
        expected_cells: Option<usize>,
        whitelist_len: Option<usize>,
        estimate: usize,
    ) -> Self {
        if let Some(n) = expected_cells {
            return Self::Explicit(n);
        }
        if let Some(n) = whitelist_len {
            return Self::Whitelist(n);
        }
        Self::Estimate(estimate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decomp_threads_inverse_to_workers() {
        assert_eq!(default_decomp_threads(1), 4);
        assert_eq!(default_decomp_threads(4), 4);
        assert_eq!(default_decomp_threads(5), 2);
        assert_eq!(default_decomp_threads(8), 0);
        assert_eq!(default_decomp_threads(24), 0);
    }

    #[test]
    fn auto_chunk_size_scales_with_workers() {
        let mm10 = 2_725_537_669u64;
        assert_eq!(auto_chunk_size_bp(mm10, 96), mm10.div_ceil(576));
        assert!(auto_chunk_size_bp(mm10, 16) > auto_chunk_size_bp(mm10, 32));
        assert_eq!(auto_chunk_size_bp(mm10, 1), MAX_AUTO_CHUNK_BP);
        assert_eq!(auto_chunk_size_bp(58_617_616, 96), MIN_AUTO_CHUNK_BP);
    }

    #[test]
    fn rss_reads_proc_status() {
        // Read RSS first: other tests allocate in parallel, so RSS read after
        // the peak could exceed it.
        let rss = process_rss_mb();
        assert!(rss > 0.0);
        assert!(process_peak_rss_mb() >= rss);
    }

    #[test]
    fn default_workers_nonzero() {
        assert!(default_workers(32.0) >= 1);
        assert_eq!(resolve_workers_with_budget(8, 32.0).0, 8);
        assert_eq!(resolve_workers_with_budget(0, 32.0).0, default_workers(32.0));
    }

    fn inputs(mode: &str, budget_gb: f64, mapped_reads: u64, n_shards: usize) -> MemoryInputs<'_> {
        MemoryInputs {
            mode,
            budget_gb,
            mapped_reads,
            site_share: 1.0,
            cell_count: 7_351,
            genome_mb: 2_725,
            n_shards,
            n_workers: 32,
            cli_writers: 0,
            profile_writers: None,
            in_process: false,
        }
    }

    #[test]
    fn memory_plan_fits_budget() {
        let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        // The 2026-10-15 test BAM on a 34 GB budget: memory mode fits.
        let p = plan_memory(&inputs("auto", 34.0, 53_461_257, 16));
        assert_eq!(p.mode, "memory");
        assert!(p.peak_gb <= 34.0 + p.writer_each_gb, "{p:?}");
        assert!(p.writers >= 1 && p.writers <= 16.min(nproc));
        // Same BAM on a 32 GB workstation (19 GB budget): stream, fewer writers.
        let p = plan_memory(&inputs("auto", 19.0, 53_461_257, 16));
        assert_eq!(p.mode, "stream");
        assert!(p.peak_gb <= 19.0 + p.writer_each_gb, "{p:?}");
        // 1B reads on 400 GB: still bounded.
        let p = plan_memory(&inputs("auto", 400.0, 1_000_000_000, 32));
        assert!(p.peak_gb <= 400.0 + p.writer_each_gb, "{p:?}");
    }

    #[test]
    fn memory_plan_explicit_and_in_process() {
        let mut i = inputs("memory", 8.0, 53_461_257, 16);
        i.cli_writers = 16;
        let p = plan_memory(&i);
        assert_eq!((p.mode.as_str(), p.writers), ("memory", 16));
        let mut i = inputs("stream", 500.0, 1_000_000, 16);
        i.in_process = true;
        assert!(plan_memory(&i).writers <= MAX_AUTO_WRITERS);
    }

    #[test]
    fn memory_fraction_scales_budget() {
        let (a, _) = system_memory_budget_gb(None, Some(0.3));
        let (b, _) = system_memory_budget_gb(None, Some(0.6));
        assert!(b >= a);
        assert_eq!(system_memory_budget_gb(Some(12.0), Some(0.1)).0, 12.0);
    }

    #[test]
    fn cell_count_whitelist_then_explicit() {
        assert_eq!(
            CellCountSource::resolve(None, Some(7355), 100).count(),
            7355
        );
        assert_eq!(
            CellCountSource::resolve(Some(10), Some(7355), 100),
            CellCountSource::Explicit(10)
        );
    }

    #[test]
    fn shards_resolution_auto() {
        assert_eq!(resolve_shards(0, 5_000).0, 1);
        assert_eq!(resolve_shards(0, 25_000).0, 8);
        assert_eq!(resolve_shards(0, 75_000).0, 16);
        assert_eq!(resolve_shards(0, 150_000).0, 32);
        assert_eq!(resolve_shards(4, 150_000).0, 4);
    }

}
