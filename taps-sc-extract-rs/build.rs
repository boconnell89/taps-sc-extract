// Builds the C parts against the vendored htslib headers:
// - csrc/baq_simd.c: batched SIMD BAQ, plus the CPU half of GPU BAQ;
// - csrc/gpu/baq_opencl.c: OpenCL GPU backend (always; libOpenCL is loaded at
//   run time, so no OpenCL is needed to build or run);
// - csrc/gpu/baq_cuda.cu: CUDA backend when `nvcc` is found, and the same file
//   as a HIP backend when `hipcc` is found (cfg taps_cuda / taps_hip).
// TAPS_NO_CUDA / TAPS_NO_HIP skip a backend; NVCC / HIPCC pick the compiler;
// TAPS_CUDA_ARCH (e.g. "86") limits CUDA code to one GPU generation.
use std::path::{Path, PathBuf};
use std::process::Command;

fn which(prog: &str) -> Option<PathBuf> {
    let out = Command::new("sh").arg("-c").arg(format!("command -v {prog}")).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let p = String::from_utf8_lossy(&out.stdout).trim().to_string();
    (!p.is_empty()).then(|| PathBuf::from(p))
}

fn main() {
    let include = std::env::var("DEP_HTS_INCLUDE").expect("hts-sys exports its include dir");
    for f in [
        "csrc/baq_simd.c",
        "csrc/olap.c",
        "csrc/baq_kernel.h",
        "csrc/gpu/baq_hmm_body.h",
        "csrc/gpu/baq_opencl.c",
        "csrc/gpu/baq_cuda.cu",
        "csrc/gpu/taps_gpu.h",
    ] {
        println!("cargo:rerun-if-changed={f}");
    }
    for v in ["TAPS_BAQ_NO_SIMD", "TAPS_NO_CUDA", "TAPS_NO_HIP", "NVCC", "HIPCC", "TAPS_CUDA_ARCH"] {
        println!("cargo:rerun-if-env-changed={v}");
    }
    println!("cargo::rustc-check-cfg=cfg(taps_cuda)");
    println!("cargo::rustc-check-cfg=cfg(taps_hip)");

    let mut b = cc::Build::new();
    // TAPS_BAQ_NO_SIMD=1 builds the scalar-only path (what non-x86 hosts get).
    if std::env::var_os("TAPS_BAQ_NO_SIMD").is_some() {
        b.define("TAPS_BAQ_NO_SIMD", None);
    }
    check_olap_verbatim();
    b.file("csrc/baq_simd.c")
        .file("csrc/olap.c")
        .file("csrc/gpu/baq_opencl.c")
        .include(&include)
        .opt_level(3)
        // Same IEEE operations as htslib's probaln.c (also built without contraction).
        .flag("-ffp-contract=off")
        .flag_if_supported("-fno-fast-math")
        // 8-lane helpers in the AVX2 region are always inlined; no ABI crossing.
        .flag_if_supported("-Wno-psabi")
        // The wrapper keeps htslib realn.c's int/uint32 comparisons verbatim.
        .flag_if_supported("-Wno-sign-compare")
        .warnings(true)
        .compile("taps_baq_simd");
    println!("cargo:rustc-link-lib=dl");

    // CUDA backend.
    let nvcc = std::env::var_os("NVCC").map(PathBuf::from).or_else(|| which("nvcc"));
    if let (Some(nvcc), None) = (nvcc, std::env::var_os("TAPS_NO_CUDA")) {
        build_cuda(&nvcc);
    }
    // HIP backend (AMD): the same source compiled by hipcc.
    let hipcc = std::env::var_os("HIPCC").map(PathBuf::from).or_else(|| which("hipcc"));
    if let (Some(hipcc), None) = (hipcc, std::env::var_os("TAPS_NO_HIP")) {
        build_hip(&hipcc);
    }
}

/// csrc/olap.c must hold htslib's overlap functions unchanged: each function
/// in it must appear verbatim in the (patched) sam.c that hts-sys compiled.
fn check_olap_verbatim() {
    let root = std::env::var("DEP_HTS_ROOT").unwrap_or_default();
    let Ok(sam) = std::fs::read_to_string(Path::new(&root).join("htslib/sam.c")) else {
        println!("cargo:warning=htslib sam.c not found under DEP_HTS_ROOT; csrc/olap.c not checked");
        return;
    };
    let olap = std::fs::read_to_string("csrc/olap.c").expect("read csrc/olap.c");
    for sig in [
        "static inline int cigar_iref2iseq_set(",
        "static inline int cigar_iref2iseq_next(",
        "static int tweak_overlap_quality(bam1_t *a, bam1_t *b)",
    ] {
        let i = olap.find(sig).expect("function in csrc/olap.c");
        let j = i + olap[i..].find("\n}\n").expect("function end") + 3;
        assert!(sam.contains(&olap[i..j]), "csrc/olap.c: {sig} differs from htslib sam.c; copy it again");
    }
}

fn build_cuda(nvcc: &Path) {
    // nvcc is run directly (not through `cc`) so conda's CFLAGS/CXXFLAGS, which
    // nvcc does not understand, never reach it.
    let out = PathBuf::from(std::env::var("OUT_DIR").expect("OUT_DIR"));
    let obj = out.join("baq_cuda.o");
    let mut cmd = Command::new(nvcc);
    cmd.args(["-c", "csrc/gpu/baq_cuda.cu", "-o"])
        .arg(&obj)
        .args(["-O3", "-Xcompiler", "-fPIC", "-cudart=static"])
        // Exact IEEE double arithmetic (matches htslib and the CPU kernels).
        .args(["-fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false"]);
    match std::env::var("TAPS_CUDA_ARCH") {
        Ok(a) if !a.trim().is_empty() => {
            let a = a.trim();
            cmd.arg(format!("-gencode=arch=compute_{a},code=sm_{a}"));
        }
        _ => {
            for a in ["70", "75", "80", "86", "89", "90"] {
                cmd.arg(format!("-gencode=arch=compute_{a},code=sm_{a}"));
            }
            // PTX for newer GPUs (JIT-compiled by the driver).
            cmd.arg("-gencode=arch=compute_90,code=compute_90");
        }
    }
    let status = cmd.status().expect("run nvcc");
    assert!(status.success(), "nvcc failed on csrc/gpu/baq_cuda.cu");
    let lib = out.join("libtaps_baq_cuda.a");
    let _ = std::fs::remove_file(&lib);
    let ar = Command::new("ar").arg("crs").arg(&lib).arg(&obj).status().expect("run ar");
    assert!(ar.success(), "ar failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=taps_baq_cuda");
    // Static CUDA runtime from the toolkit next to nvcc; it loads the driver (libcuda) at run time.
    let root = nvcc.parent().and_then(|p| p.parent()).map(Path::to_path_buf).unwrap_or_default();
    for d in ["lib64", "lib", "targets/x86_64-linux/lib"] {
        let p = root.join(d);
        if p.join("libcudart_static.a").exists() {
            println!("cargo:rustc-link-search=native={}", p.display());
            break;
        }
    }
    println!("cargo:rustc-link-lib=static=cudart_static");
    println!("cargo:rustc-link-lib=rt");
    println!("cargo:rustc-link-lib=pthread");
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-cfg=taps_cuda");
}

fn build_hip(hipcc: &Path) {
    // Same kernel and host code as CUDA, compiled by hipcc (run directly, so
    // conda's compiler flags never reach it). Not tested on AMD hardware yet.
    let out = PathBuf::from(std::env::var("OUT_DIR").expect("OUT_DIR"));
    let obj = out.join("baq_hip.o");
    let status = Command::new(hipcc)
        .args(["-c", "-x", "hip", "csrc/gpu/baq_cuda.cu", "-o"])
        .arg(&obj)
        .args(["-DTAPS_HIP", "-O3", "-fPIC"])
        // Exact IEEE double arithmetic.
        .args(["-ffp-contract=off", "-fno-fast-math"])
        .status()
        .expect("run hipcc");
    assert!(status.success(), "hipcc failed on csrc/gpu/baq_cuda.cu");
    let lib = out.join("libtaps_baq_hip.a");
    let _ = std::fs::remove_file(&lib);
    let ar = Command::new("ar").arg("crs").arg(&lib).arg(&obj).status().expect("run ar");
    assert!(ar.success(), "ar failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=taps_baq_hip");
    let root = hipcc.parent().and_then(|p| p.parent()).map(Path::to_path_buf).unwrap_or_default();
    println!("cargo:rustc-link-search=native={}", root.join("lib").display());
    println!("cargo:rustc-link-lib=amdhip64");
    println!("cargo:rustc-link-lib=stdc++");
    println!("cargo:rustc-cfg=taps_hip");
}
