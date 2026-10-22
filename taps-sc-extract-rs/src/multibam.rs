//! Reading one window from several BAMs as one coordinate-sorted stream.
//!
//! Each input's BAI is loaded once and shared by every worker thread (an
//! `hts_idx_t` is read-only after loading); each thread opens its own file
//! handle and iterator. Opened per thread, a BAI cost ~13.7 MB per reader on
//! the 2026-10-21 test files, ~17 GB at 40 files x 32 workers. CRAM indexes
//! belong to their open file, so CRAM inputs use a per-thread `IndexedReader`.
//!
//! `Merged` k-way merges the inputs by (tid, pos), ties to the earlier input,
//! with one lookahead record per input. With one input it passes reads
//! straight through, so single-BAM runs read exactly as before.

use std::ffi::CString;
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context as _, Result};
use rust_htslib::bam::{IndexedReader, Record};
use rust_htslib::htslib;

/// A window's reads from one input or several.
pub trait WindowSource {
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<()>;
    /// Next read in coordinate order: `None` at the end of the window.
    fn read(&mut self, rec: &mut Record) -> Option<Result<()>>;
}

impl WindowSource for IndexedReader {
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<()> {
        rust_htslib::bam::IndexedReader::fetch(self, (contig, start, end)).map_err(|e| anyhow::anyhow!("BAM fetch: {e}"))
    }

    fn read(&mut self, rec: &mut Record) -> Option<Result<()>> {
        rust_htslib::bam::Read::read(self, rec).map(|r| r.map_err(|e| anyhow::anyhow!("BAM read: {e}")))
    }
}

/// A BAM index loaded once and shared read-only by every thread.
pub struct SharedIndex(*mut htslib::hts_idx_t);

// SAFETY: htslib only reads an `hts_idx_t` when building iterators
// (`sam_itr_queryi` copies what it needs into the iterator); nothing mutates it
// after loading, and it is destroyed only when the last Arc drops.
unsafe impl Send for SharedIndex {}
unsafe impl Sync for SharedIndex {}

impl Drop for SharedIndex {
    fn drop(&mut self) {
        unsafe { htslib::hts_idx_destroy(self.0) };
    }
}

fn cpath(p: &Path) -> Result<CString> {
    CString::new(p.as_os_str().as_bytes()).with_context(|| format!("path {} has a NUL byte", p.display()))
}

impl SharedIndex {
    pub fn load(bam: &Path) -> Result<Self> {
        let fname = cpath(bam)?;
        let mode = CString::new("r").expect("static");
        let fp = unsafe { htslib::hts_open(fname.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            anyhow::bail!("open {}", bam.display());
        }
        let fnidx = crate::input::index_path(bam).map(|p| cpath(&p)).transpose()?;
        let idx = unsafe {
            htslib::sam_index_load2(fp, fname.as_ptr(), fnidx.as_ref().map_or(std::ptr::null(), |c| c.as_ptr()))
        };
        unsafe { htslib::hts_close(fp) };
        if idx.is_null() {
            anyhow::bail!("load the index of {}", bam.display());
        }
        Ok(Self(idx))
    }
}

/// One thread's handle on a BAM whose index is shared.
pub struct RawReader {
    fp: *mut htslib::htsFile,
    hdr: *mut htslib::sam_hdr_t,
    idx: Arc<SharedIndex>,
    itr: *mut htslib::hts_itr_t,
    path: PathBuf,
}

impl RawReader {
    pub fn open(path: &Path, idx: Arc<SharedIndex>, decomp_threads: usize) -> Result<Self> {
        let fname = cpath(path)?;
        let mode = CString::new("r").expect("static");
        let fp = unsafe { htslib::hts_open(fname.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            anyhow::bail!("open {}", path.display());
        }
        let hdr = unsafe { htslib::sam_hdr_read(fp) };
        if hdr.is_null() {
            unsafe { htslib::hts_close(fp) };
            anyhow::bail!("read the header of {}", path.display());
        }
        if decomp_threads > 0 && unsafe { htslib::hts_set_threads(fp, decomp_threads as i32) } != 0 {
            unsafe {
                htslib::sam_hdr_destroy(hdr);
                htslib::hts_close(fp);
            }
            anyhow::bail!("set {decomp_threads} decompression threads on {}", path.display());
        }
        Ok(Self { fp, hdr, idx, itr: std::ptr::null_mut(), path: path.to_path_buf() })
    }
}

impl Drop for RawReader {
    fn drop(&mut self) {
        unsafe {
            if !self.itr.is_null() {
                htslib::hts_itr_destroy(self.itr);
            }
            htslib::sam_hdr_destroy(self.hdr);
            htslib::hts_close(self.fp);
        }
    }
}

impl WindowSource for RawReader {
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<()> {
        let name = CString::new(contig).context("contig name")?;
        let tid = unsafe { htslib::sam_hdr_name2tid(self.hdr, name.as_ptr()) };
        if tid < 0 {
            anyhow::bail!("contig {contig} not in the header of {}", self.path.display());
        }
        unsafe {
            if !self.itr.is_null() {
                htslib::hts_itr_destroy(self.itr);
            }
            self.itr = htslib::sam_itr_queryi(self.idx.0, tid, start as i64, end as i64);
        }
        if self.itr.is_null() {
            anyhow::bail!("BAM fetch {contig}:{start}-{end} in {}", self.path.display());
        }
        Ok(())
    }

    fn read(&mut self, rec: &mut Record) -> Option<Result<()>> {
        if self.itr.is_null() {
            return None;
        }
        // Same call rust-htslib's IndexedReader::read makes.
        let r = unsafe {
            htslib::hts_itr_next((*self.fp).fp.bgzf, self.itr, rec.inner_mut() as *mut htslib::bam1_t as *mut _, self.fp as *mut _)
        };
        match r {
            -1 => None,
            r if r < -1 => Some(Err(anyhow::anyhow!("BAM read {}: htslib error {r}", self.path.display()))),
            _ => Some(Ok(())),
        }
    }
}

/// The run's inputs: paths plus shared indexes (None for CRAM), built once.
pub struct Sources {
    /// Unique per load, so thread-local readers from an earlier run are reopened.
    pub id: u64,
    pub paths: Vec<PathBuf>,
    indexes: Vec<Option<Arc<SharedIndex>>>,
}

impl Sources {
    pub fn load(paths: &[PathBuf]) -> Result<Self> {
        let indexes = paths
            .iter()
            .map(|p| if crate::input::is_cram_path(p) { Ok(None) } else { SharedIndex::load(p).map(|i| Some(Arc::new(i))) })
            .collect::<Result<Vec<_>>>()?;
        static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);
        Ok(Self { id: NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed), paths: paths.to_vec(), indexes })
    }

    /// A thread's readers over every input.
    pub fn open(&self, decomp_threads: usize) -> Result<Merged> {
        let mut srcs: Vec<Box<dyn WindowSource>> = Vec::with_capacity(self.paths.len());
        for (p, idx) in self.paths.iter().zip(&self.indexes) {
            srcs.push(match idx {
                Some(i) => Box::new(RawReader::open(p, Arc::clone(i), decomp_threads)?),
                None => Box::new(crate::extract::open_bam(p, decomp_threads)?),
            });
        }
        Ok(Merged::new(srcs))
    }
}

/// Raise the soft open-file limit toward the hard one: a run keeps
/// workers x inputs files open. Errors if even the hard limit is too low.
pub fn ensure_open_files(needed: u64) -> Result<()> {
    let mut lim = libc::rlimit { rlim_cur: 0, rlim_max: 0 };
    if unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut lim) } != 0 {
        return Ok(());
    }
    let want = needed + 256;
    if lim.rlim_cur >= want {
        return Ok(());
    }
    if lim.rlim_max < want {
        anyhow::bail!(
            "{needed} open files are needed (workers x inputs) but the hard limit is {}; lower -t or raise `ulimit -Hn`",
            lim.rlim_max
        );
    }
    lim.rlim_cur = want;
    if unsafe { libc::setrlimit(libc::RLIMIT_NOFILE, &lim) } != 0 {
        anyhow::bail!("raise the open-file limit to {want}");
    }
    Ok(())
}

/// k-way merge of a window's reads over several inputs.
pub struct Merged {
    srcs: Vec<Box<dyn WindowSource>>,
    heads: Vec<Record>,
    live: Vec<bool>,
}

impl Merged {
    pub fn new(srcs: Vec<Box<dyn WindowSource>>) -> Self {
        let n = srcs.len();
        Self { srcs, heads: (0..n).map(|_| Record::new()).collect(), live: vec![false; n] }
    }

    pub fn len(&self) -> usize {
        self.srcs.len()
    }
}

impl WindowSource for Merged {
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<()> {
        for s in &mut self.srcs {
            s.fetch(contig, start, end)?;
        }
        if self.srcs.len() > 1 {
            for i in 0..self.srcs.len() {
                self.live[i] = match self.srcs[i].read(&mut self.heads[i]) {
                    Some(r) => {
                        r?;
                        true
                    }
                    None => false,
                };
            }
        }
        Ok(())
    }

    fn read(&mut self, rec: &mut Record) -> Option<Result<()>> {
        if self.srcs.len() == 1 {
            return self.srcs[0].read(rec);
        }
        let mut best: Option<usize> = None;
        for i in 0..self.srcs.len() {
            if !self.live[i] {
                continue;
            }
            let key = (self.heads[i].tid(), self.heads[i].pos());
            if best.is_none_or(|b| key < (self.heads[b].tid(), self.heads[b].pos())) {
                best = Some(i);
            }
        }
        let i = best?;
        std::mem::swap(rec, &mut self.heads[i]);
        match self.srcs[i].read(&mut self.heads[i]) {
            Some(Ok(())) => {}
            Some(Err(e)) => return Some(Err(e)),
            None => self.live[i] = false,
        }
        Some(Ok(()))
    }
}
