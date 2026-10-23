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
pub const DEFAULT_MEMORY_FRACTION: f64 = 0.85;

/// Memory budget for the whole process tree (extraction plus writer children):
/// `--max-memory-gb` if given, else `fraction` (default 0.85) of MemAvailable.
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
        default_workers(get_available_memory_gb() * DEFAULT_MEMORY_FRACTION)
    } else {
        cli
    }
}

/// Output sites per shard when shards are auto. The TAPS reference run
/// (2026-10-18) wrote 702M sites as 16 shards, ~44M each; a bisulfite library
/// (2026-10-21, 1.15B sites, 304 cells) took 33 s after extraction as one shard
/// and 5-10 s per shard as 16.
pub const AUTO_SITES_PER_SHARD: f64 = 48e6;

/// Resolve shard count (0 means auto). Auto takes the larger of a cell-count
/// rule and estimated sites / `AUTO_SITES_PER_SHARD` (few deep cells still
/// get parallel writers), capped at max(cell rule, nproc / 2) and the cell
/// count. `-o` naming one `.h5` file keeps one shard.
pub fn resolve_shards(
    cli_shards: usize,
    cell_count: usize,
    est_sites: f64,
    nproc: usize,
    single_file: bool,
) -> (usize, &'static str) {
    if cli_shards > 0 {
        return (cli_shards, "explicit");
    }
    if single_file {
        return (1, "auto (-o names one .h5 file)");
    }
    let by_cells = if cell_count <= 10_000 {
        1
    } else if cell_count <= 50_000 {
        8
    } else if cell_count <= 100_000 {
        16
    } else {
        32
    };
    let by_sites = (est_sites / AUTO_SITES_PER_SHARD).ceil().max(1.0) as usize;
    let cap = by_cells.max(nproc / 2).max(1);
    let s = by_sites.max(by_cells).min(cap).min(cell_count.max(1));
    (s, if s > by_cells { "auto, by estimated sites" } else { "auto" })
}

/// Output sites expected from the BAM index (or a per-cell guess without one).
pub fn estimate_sites(mapped_reads: u64, site_share: f64, cell_count: usize, genome_mb: u64) -> f64 {
    if mapped_reads > 0 {
        mapped_reads as f64 * SITES_PER_READ * site_share
    } else {
        cell_count as f64 * (genome_mb * 500).min(200_000) as f64
    }
}

/// Memory model, measured 2026-10-15 on mm10 (53.5M reads, 7,351 cells, 702M sites):
/// 13.1 sites per mapped read; memory-mode extraction peak ~16.4 B/site; the 16
/// sliced `write-shard` children together ~14.5 B/site (10.2 GB). Margins added.
pub const SITES_PER_READ: f64 = 14.0;
pub const MEMORY_MODE_BYTES_PER_SITE: f64 = 18.0;
/// A writer holds its whole shard: measured 12.6-12.7 B/site on 32 shards of a
/// 54.7B-site bisulfite run and ~14.5 B/site on the TAPS test BAM (2026-10-22).
/// Used at startup only; the write phase is re-planned from measured shard
/// sizes (`plan_write`).
pub const WRITER_BYTES_PER_SITE: f64 = 15.0;
pub const WRITER_BASE_GB: f64 = 0.1;
/// Writers fed during extraction (memory mode, one process per shard): the 16
/// children together measured 10.4 GB for 702M sites (14.8 B/site, 2026-10-18).
pub const STREAMED_WRITER_BYTES_PER_SITE: f64 = 15.0;
/// Share of the payloads the extraction process holds before buckets are dealt
/// to shards (`ShardStreams`: after 1/5 of the windows).
pub const STREAM_HOLD_FRACTION: f64 = 0.2;

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
    /// Writer children can be fed during extraction (process writers, shard
    /// directory output, `TAPS_STREAM_WRITERS` not 0).
    pub streamable: bool,
}

pub fn plan_memory(i: &MemoryInputs) -> MemoryPlan {
    let nproc = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let sites = estimate_sites(i.mapped_reads, i.site_share, i.cell_count, i.genome_mb);
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
    // Memory mode with one writer per shard streams payloads to the writers
    // during extraction: the parent holds only what arrives before buckets are
    // dealt, and the children hold the rest. Preferred whenever it fits.
    let mode_lc = i.mode.to_lowercase();
    if i.streamable
        && n_shards > 1
        && mode_lc != "stream"
        && n_shards <= writer_cap
        && (i.cli_writers == 0 || i.cli_writers >= n_shards)
    {
        let parent_gb = stream_gb + STREAM_HOLD_FRACTION * payload_gb;
        let writer_each = sites / n_shards as f64 * STREAMED_WRITER_BYTES_PER_SITE / 1e9 + WRITER_BASE_GB;
        let peak_gb = parent_gb + n_shards as f64 * writer_each;
        if peak_gb <= i.budget_gb || mode_lc == "memory" {
            let writers = if i.cli_writers > 0 { i.cli_writers } else { n_shards };
            return MemoryPlan {
                mode: "memory".to_string(),
                mode_src: if mode_lc == "memory" { "explicit" } else { "auto (writers fed during extraction)" },
                writers,
                writers_src: if i.cli_writers > 0 { "explicit" } else { "auto, one per shard, fed during extraction" },
                est_sites: sites,
                extract_gb: parent_gb,
                parent_during_write_gb: STREAM_BASE_GB,
                writer_each_gb: writer_each,
                peak_gb,
            };
        }
    }
    let (mode, mode_src) = match mode_lc.as_str() {
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

/// Writer RSS per byte of shard payload (temp-file bytes, or 12 B per site in
/// memory: the TAPSCK01 site record and the in-memory one are both 12 bytes),
/// plus a base per process. Measured 2026-10-22: RSS ~= 1.0 x payload + 0.1 GB
/// (20.6 GB for 20.5 GB on the server; 0.73 GB for 0.62 GB locally). 15% and
/// 0.15 GB of margin: the 32-writer crash came within 1% of MemAvailable.
pub const WRITER_RSS_PER_BYTE: f64 = 1.15;
pub const WRITER_RSS_BASE_BYTES: f64 = 0.25e9;

/// Most output files the write phase may split a run into (per CPU).
pub const MAX_FILES_PER_CPU: usize = 2;

/// How the write phase runs, decided after extraction from measured sizes.
#[derive(Clone, Debug)]
pub struct WritePlan {
    /// Buckets per output file (`groups.len()` files).
    pub groups: Vec<Vec<usize>>,
    /// Writer processes at once.
    pub writers: usize,
    /// Estimated RSS of the largest file's writer (bytes).
    pub writer_rss: f64,
    pub notes: Vec<String>,
}

fn writer_rss(bytes: u64) -> f64 {
    bytes as f64 * WRITER_RSS_PER_BYTE + WRITER_RSS_BASE_BYTES
}

/// Choose the number of output files and concurrent writers from the
/// per-bucket payload bytes, so the writers that run at once fit `budget_bytes`
/// (memory available for writers now). With automatic shards, files range from
/// `n_shards` up to the bucket count (whole buckets only; a cell never splits),
/// at most `MAX_FILES_PER_CPU` per CPU; the count with the shortest estimated
/// write (rounds of writers x largest file) wins, fewer files on ties. `cap`
/// bounds writers (`--max-writer-threads`, else the machine profile, else
/// CPUs). Explicit settings are followed as given: an explicit `--shards`
/// keeps its file count and an explicit `--max-writer-threads` runs that many
/// writers; when they do not fit the budget a warning is added instead.
pub fn plan_write(
    bucket_bytes: &[u64],
    n_shards: usize,
    shards_explicit: bool,
    cap: usize,
    cap_explicit: bool,
    budget_bytes: f64,
    nproc: usize,
) -> WritePlan {
    let n_shards = n_shards.max(1);
    let cap = cap.max(1);
    let most = if n_shards <= 1 {
        1
    } else if shards_explicit {
        n_shards
    } else {
        bucket_bytes.len().clamp(n_shards, n_shards.max(MAX_FILES_PER_CPU * nproc.max(1)))
    };
    let mut best: Option<(f64, WritePlan)> = None;
    for k in n_shards..=most {
        let groups = crate::h5_out::pack_buckets(bucket_bytes, k);
        let largest = groups.iter().map(|g| g.iter().map(|&b| bucket_bytes[b]).sum::<u64>()).max().unwrap_or(0);
        let rss = writer_rss(largest);
        let fit = ((budget_bytes / rss).floor().max(1.0)) as usize;
        let writers = if cap_explicit { k.min(cap) } else { k.min(cap).min(fit).max(1) };
        let rounds = k.div_ceil(writers) as f64;
        let cost = rounds * largest.max(1) as f64;
        if best.as_ref().is_none_or(|(c, _)| cost < *c * 0.999) {
            best = Some((cost, WritePlan { groups, writers, writer_rss: rss, notes: Vec::new() }));
        }
    }
    let mut plan = best.expect("at least one file count").1;
    let files = plan.groups.len();
    let gb = |b: f64| b / 1e9;
    if files > n_shards {
        plan.notes.push(format!(
            "writing {files} files instead of {n_shards} so that {} writers fit in {:.0} GB (~{:.1} GB each)",
            plan.writers,
            gb(budget_bytes),
            gb(plan.writer_rss)
        ));
    }
    let fit = ((budget_bytes / plan.writer_rss).floor().max(1.0)) as usize;
    if plan.writers > fit {
        plan.notes.push(format!(
            "WARNING: {} writers of ~{:.1} GB each need ~{:.0} GB, but only {:.0} GB is free; the write may run out of memory. \
             About {fit} writers fit{}",
            plan.writers,
            gb(plan.writer_rss),
            gb(plan.writers as f64 * plan.writer_rss),
            gb(budget_bytes),
            if shards_explicit { "; more --shards would make each writer smaller" } else { "" }
        ));
    } else if plan.writer_rss > budget_bytes {
        plan.notes.push(format!(
            "WARNING: one writer needs ~{:.1} GB but only {:.0} GB is free; the write may run out of memory",
            gb(plan.writer_rss),
            gb(budget_bytes)
        ));
    }
    plan
}

/// Memory the write phase may use now: `--max-memory-gb` minus what this
/// process holds, else `fraction` of MemAvailable read now.
pub fn write_budget_bytes(max_memory_gb: Option<f64>, fraction: Option<f64>) -> f64 {
    let f = fraction.unwrap_or(DEFAULT_MEMORY_FRACTION).clamp(0.05, 1.0);
    let avail = get_available_memory_gb() * f;
    let gb = match max_memory_gb {
        Some(m) => (m - process_rss_mb() / 1024.0).min(avail),
        None => avail,
    };
    gb.max(0.5) * 1e9
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
            streamable: false,
        }
    }

    #[test]
    fn memory_plan_streams_writers_when_it_fits() {
        if std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) < 16 {
            return; // 16 writers need 16 CPUs (writer_cap)
        }
        // A bisulfite library (2026-10-21): 91.6M reads, 16 shards, 34 GB budget.
        // Without streaming, memory mode does not fit (it spills to disk);
        // streamed, memory mode with all 16 writers does.
        let mut i = inputs("auto", 34.0, 91_642_549, 16);
        assert_eq!(plan_memory(&i).mode, "stream");
        i.streamable = true;
        let p = plan_memory(&i);
        assert_eq!((p.mode.as_str(), p.writers), ("memory", 16), "{p:?}");
        assert!(p.peak_gb <= 34.0, "{p:?}");
        // Too small a budget: falls back to the unstreamed plan.
        i.budget_gb = 19.0;
        let p = plan_memory(&i);
        assert!(p.peak_gb <= 19.0 + p.writer_each_gb, "{p:?}");
        assert!(p.writers < 16);
        // --memory-mode stream is never overridden.
        let mut i = inputs("stream", 400.0, 91_642_549, 16);
        i.streamable = true;
        assert_eq!(plan_memory(&i).mode, "stream");
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
        let mut i = inputs("memory", 400.0, 53_461_257, 16);
        i.cli_writers = 16;
        let p = plan_memory(&i);
        assert_eq!((p.mode.as_str(), p.writers), ("memory", 16));
        // Explicit memory mode is followed even when it cannot fit; the peak
        // estimate shows it (main.rs warns).
        let mut i = inputs("memory", 8.0, 53_461_257, 16);
        i.cli_writers = 16;
        let p = plan_memory(&i);
        assert_eq!(p.mode, "memory");
        assert!(p.peak_gb > 8.0);
        let mut i = inputs("stream", 500.0, 1_000_000, 16);
        i.in_process = true;
        assert!(plan_memory(&i).writers <= MAX_AUTO_WRITERS);
    }

    #[test]
    fn write_plan_fits_the_server_run() {
        // 2026-10-22 server run: 54.7B sites in 32 shards (128 buckets), 12 B
        // per site of payload, 0.85 x ~650 GB free after extraction, 96 CPUs.
        let site = 54.7e9 / 128.0;
        let buckets: Vec<u64> = (0..128).map(|i| (site * (0.95 + 0.1 * (i % 3) as f64 / 2.0) * 12.0) as u64).collect();
        let budget = 0.85 * 650e9;
        // Asked for 32 writers on 32 files: 32 x ~25 GB cannot fit, so the
        // files are split further until 32 run at once, never 32 x 25 GB.
        // Automatic shards and writers: the files are split further until the
        // writers that run at once fit, never 32 x 25 GB.
        let p = plan_write(&buckets, 32, false, 96, false, budget, 96);
        assert!(p.writers as f64 * p.writer_rss <= budget, "{p:?}");
        assert!(p.groups.len() > 32 && p.writers >= 24, "{p:?}");
        assert_eq!(p.groups.iter().map(|g| g.len()).sum::<usize>(), 128);
        assert!(!p.notes.iter().any(|n| n.contains("WARNING")));
        // Explicit --shards 32 --max-writer-threads 32 (the run that crashed):
        // followed as given, with a warning.
        let p = plan_write(&buckets, 32, true, 32, true, budget, 96);
        assert_eq!((p.groups.len(), p.writers), (32, 32));
        assert!(p.notes.iter().any(|n| n.contains("WARNING") && n.contains("may run out of memory")), "{p:?}");
        // Explicit --shards 32, automatic writers: as many as fit, no split.
        let p = plan_write(&buckets, 32, true, 96, false, budget, 96);
        assert_eq!(p.groups.len(), 32);
        assert!(p.writers < 32 && p.writers as f64 * p.writer_rss <= budget, "{p:?}");
        // One file stays one file.
        let one = plan_write(&buckets[..1], 1, false, 32, false, budget, 96);
        assert_eq!((one.groups.len(), one.writers), (1, 1));
        // Plenty of memory: the requested shards, all writers at once, no notes.
        let p = plan_write(&buckets, 32, false, 32, false, 5e12, 96);
        assert_eq!((p.groups.len(), p.writers), (32, 32));
        assert!(p.notes.is_empty());
        // Tiny budget: one writer, warned.
        let p = plan_write(&buckets, 32, false, 32, false, 1e9, 4);
        assert_eq!(p.writers, 1);
        assert!(p.notes.iter().any(|n| n.contains("may run out of memory")));
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
        // Cell rule (few sites).
        assert_eq!(resolve_shards(0, 5_000, 1e6, 32, false).0, 1);
        assert_eq!(resolve_shards(0, 25_000, 1e6, 32, false).0, 8);
        assert_eq!(resolve_shards(0, 75_000, 1e6, 32, false).0, 16);
        assert_eq!(resolve_shards(0, 150_000, 1e6, 32, false).0, 32);
        assert_eq!(resolve_shards(4, 150_000, 1e9, 32, false).0, 4);
        // Sites: TAPS test BAM (53.5M reads x 14) and a bisulfite library (91.6M reads, 304 cells).
        assert_eq!(resolve_shards(0, 7_351, 53_461_257.0 * SITES_PER_READ, 32, false).0, 16);
        assert_eq!(resolve_shards(0, 304, 91_642_549.0 * SITES_PER_READ, 32, false).0, 16);
        assert_eq!(resolve_shards(0, 304, 91_642_549.0 * SITES_PER_READ, 16, false).0, 8);
        assert_eq!(resolve_shards(0, 3, 1e10, 96, false).0, 3);
        // One .h5 file asked for: one shard.
        assert_eq!(resolve_shards(0, 304, 1e10, 32, true).0, 1);
    }

}
