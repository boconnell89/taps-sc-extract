//! Conservative parameter selection. Explicit CLI values always win.

pub const DEFAULT_CHUNK_SIZE_MB: u32 = 10;
pub const MAX_SHARD_WRITERS: usize = 6;
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

/// Compute conservative memory budget (explicit --max-memory-gb or 0.6 * MemAvailable).
pub fn system_memory_budget_gb(max_memory_gb: Option<f64>) -> (f64, &'static str) {
    if let Some(gb) = max_memory_gb {
        (gb.max(1.0), "explicit")
    } else {
        let avail = get_available_memory_gb();
        ((avail * 0.6).max(2.0), "auto (0.6 * MemAvailable)")
    }
}

/// BGZF decompression threads per worker when the user does not pass `--decomp-threads`.
/// More threads when fewer extraction workers.
pub fn default_decomp_threads(n_workers: usize) -> usize {
    if n_workers <= 4 {
        4
    } else if n_workers <= 7 {
        2
    } else {
        1
    }
}

pub fn default_chunk_size_mb() -> u32 {
    DEFAULT_CHUNK_SIZE_MB
}

/// Window workers when `--workers` is omitted. Cap at min(nproc, floor((budget-2)/0.7), 64).
pub fn default_workers(budget_gb: f64) -> usize {
    let nproc = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mem_workers = ((budget_gb - 2.0).max(0.7) / 0.7).floor() as usize;
    nproc.min(mem_workers).clamp(1, 64)
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

/// Resolve memory mode: stream vs memory.
/// "auto" chooses memory mode if expected_cells * genome_sites * 40B * 2.5 < 0.4 * budget.
pub fn resolve_memory_mode(
    mode: &str,
    budget_gb: f64,
    cell_count: usize,
    total_genome_mb: u64,
) -> (String, &'static str) {
    match mode.to_lowercase().as_str() {
        "stream" => ("stream".to_string(), "explicit"),
        "memory" => ("memory".to_string(), "explicit"),
        _ => {
            // Estimated CpG + CH sites observed per cell across genome window
            let estimated_sites_per_cell = (total_genome_mb * 500).min(200_000) as f64;
            // 40 bytes per record * 2.5 heap expansion/wrapper multiplier
            let estimated_data_bytes = cell_count as f64 * estimated_sites_per_cell * 40.0 * 2.5;
            let estimated_data_gb = estimated_data_bytes / 1e9;
            if estimated_data_gb < (0.4 * budget_gb) {
                ("memory".to_string(), "auto (fits in RAM budget)")
            } else {
                ("stream".to_string(), "auto (spill to temp disk)")
            }
        }
    }
}

pub fn shard_writer_concurrency(n_shards: usize) -> usize {
    n_shards.max(1).min(MAX_SHARD_WRITERS)
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
        assert_eq!(default_decomp_threads(8), 1);
        assert_eq!(default_decomp_threads(24), 1);
    }

    #[test]
    fn chunk_size_default_is_ten() {
        assert_eq!(default_chunk_size_mb(), 10);
    }

    #[test]
    fn default_workers_nonzero() {
        assert!(default_workers(32.0) >= 1);
        assert_eq!(resolve_workers_with_budget(8, 32.0).0, 8);
        assert_eq!(resolve_workers_with_budget(0, 32.0).0, default_workers(32.0));
    }

    #[test]
    fn writer_cap_six() {
        assert_eq!(shard_writer_concurrency(16), 6);
        assert_eq!(shard_writer_concurrency(2), 2);
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

    #[test]
    fn memory_mode_auto_heuristic() {
        // High budget + few cells => memory mode
        let (mode1, _) = resolve_memory_mode("auto", 64.0, 500, 2800);
        assert_eq!(mode1, "memory");

        // Small budget + many cells => stream mode
        let (mode2, _) = resolve_memory_mode("auto", 8.0, 50_000, 2800);
        assert_eq!(mode2, "stream");
    }
}
