//! `--profile-json`: per-phase CPU use, per-stage and per-window timings, and
//! shard-writer reports for one `extract` run. Recorded from inside the process
//! (no external sampler), so the numbers line up with the run summary.

use crate::extract::{WindowResult, WindowTiming};
use crate::h5_out::WriterReport;
use crate::window::Window;
use anyhow::{Context, Result};
use std::fmt::Write as _;
use std::path::Path;

/// Kernel clock ticks per second for `/proc/<pid>/stat` times (USER_HZ, 100 on Linux).
const CLK_TCK: f64 = 100.0;

/// CPU seconds from `/proc/self/stat`: (self user+sys, waited-for children user+sys).
#[derive(Clone, Copy, Debug, Default)]
pub struct CpuTimes {
    pub own_s: f64,
    pub children_s: f64,
}

pub fn cpu_times() -> CpuTimes {
    let Ok(stat) = std::fs::read_to_string("/proc/self/stat") else {
        return CpuTimes::default();
    };
    // Fields after the ")" that closes comm: index 0 is field 3 (state).
    let Some(rest) = stat.rsplit_once(')').map(|(_, r)| r) else {
        return CpuTimes::default();
    };
    let f: Vec<f64> = rest
        .split_whitespace()
        .map(|x| x.parse().unwrap_or(0.0))
        .collect();
    // utime=14, stime=15, cutime=16, cstime=17 (1-based) -> indices 11..=14 here.
    if f.len() < 15 {
        return CpuTimes::default();
    }
    CpuTimes {
        own_s: (f[11] + f[12]) / CLK_TCK,
        children_s: (f[13] + f[14]) / CLK_TCK,
    }
}

pub struct RunProfile<'a> {
    pub config: Vec<(&'static str, String)>,
    pub process_s: f64,
    pub write_s: f64,
    pub cpu_start: CpuTimes,
    pub cpu_after_process: CpuTimes,
    pub cpu_end: CpuTimes,
    pub workers: usize,
    pub writers: usize,
    pub result: &'a WindowResult,
    pub windows: &'a [Window],
    pub writer_reports: &'a [WriterReport],
    pub master_s: f64,
    pub peak_rss_mb: f64,
    pub writer_peak_mb: f64,
}

fn esc(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn timing_json(t: &WindowTiming) -> String {
    format!(
        "\"records\": {}, \"pushed\": {}, \"total_s\": {:.4}, \"setup_s\": {:.4}, \"read_s\": {:.4}, \"baq_s\": {:.4}, \"pileup_s\": {:.4}, \"call_s\": {:.4}",
        t.records, t.pushed, t.total_s, t.setup_s, t.read_s, t.baq_s, t.pileup_s, t.call_s
    )
}

pub fn write_profile_json(path: &Path, p: &RunProfile) -> Result<()> {
    let extract_cpu = p.cpu_after_process.own_s - p.cpu_start.own_s;
    let write_cpu_own = p.cpu_end.own_s - p.cpu_after_process.own_s;
    let write_cpu_children = p.cpu_end.children_s - p.cpu_after_process.children_s;
    let util = |cpu: f64, wall: f64, threads: usize| {
        if wall > 0.0 && threads > 0 {
            cpu / (wall * threads as f64)
        } else {
            0.0
        }
    };
    let stages = &p.result.timing;
    let mut o = String::new();
    o.push_str("{\n  \"config\": {");
    for (i, (k, v)) in p.config.iter().enumerate() {
        if i > 0 {
            o.push_str(", ");
        }
        let _ = write!(o, "\"{k}\": \"{}\"", esc(v));
    }
    o.push_str("},\n");
    let _ = writeln!(
        o,
        "  \"phases\": {{\"process_s\": {:.3}, \"write_s\": {:.3}, \"total_s\": {:.3}}},",
        p.process_s,
        p.write_s,
        p.process_s + p.write_s
    );
    let _ = writeln!(
        o,
        "  \"cpu\": {{\"extract_cpu_s\": {extract_cpu:.2}, \"extract_util\": {:.3}, \"write_cpu_s\": {:.2}, \"write_children_cpu_s\": {write_cpu_children:.2}, \"write_util\": {:.3}}},",
        util(extract_cpu, p.process_s, p.workers),
        write_cpu_own,
        util(write_cpu_own + write_cpu_children, p.write_s, p.writers),
    );
    let _ = writeln!(
        o,
        "  \"memory\": {{\"parent_peak_rss_mb\": {:.1}, \"writer_peak_rss_mb\": {:.1}}},",
        p.peak_rss_mb, p.writer_peak_mb
    );
    let _ = writeln!(o, "  \"stages\": {{{}}},", timing_json(stages));
    o.push_str("  \"windows\": [\n");
    let mut wt: Vec<&(u32, WindowTiming)> = p.result.window_timings.iter().collect();
    wt.sort_by_key(|(id, _)| *id);
    for (i, (id, t)) in wt.iter().enumerate() {
        let (contig, start, end) = p
            .windows
            .iter()
            .find(|w| w.chunk_id == *id)
            .map(|w| (w.contig.as_str(), w.start, w.end))
            .unwrap_or(("?", 0, 0));
        let _ = write!(
            o,
            "    {{\"chunk_id\": {id}, \"contig\": \"{}\", \"start\": {start}, \"end\": {end}, {}}}",
            esc(contig),
            timing_json(t)
        );
        o.push_str(if i + 1 < wt.len() { ",\n" } else { "\n" });
    }
    let _ = write!(o, "  ],\n  \"master_h5_s\": {:.4},\n  \"writers\": [\n", p.master_s);
    for (i, r) in p.writer_reports.iter().enumerate() {
        let _ = write!(
            o,
            "    {{\"shard\": {}, \"child\": {}, \"cells\": {}, \"sites\": {}, \"read_s\": {:.4}, \"write_s\": {:.4}, \"create_s\": {:.4}, \"fill_s\": {:.4}, \"h5_write_s\": {:.4}, \"bytes\": {}, \"peak_rss_mb\": {:.1}}}",
            r.shard, r.child, r.cells, r.sites, r.read_s, r.write_s, r.create_s, r.fill_s, r.h5_write_s, r.bytes, r.peak_rss_mb
        );
        o.push_str(if i + 1 < p.writer_reports.len() { ",\n" } else { "\n" });
    }
    o.push_str("  ]\n}\n");
    std::fs::write(path, o).with_context(|| format!("write {}", path.display()))
}
