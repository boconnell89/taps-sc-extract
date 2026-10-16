// Compiles csrc/baq_simd.c (batched SIMD BAQ) against the vendored htslib headers.
fn main() {
    let include = std::env::var("DEP_HTS_INCLUDE").expect("hts-sys exports its include dir");
    println!("cargo:rerun-if-changed=csrc/baq_simd.c");
    println!("cargo:rerun-if-changed=csrc/baq_kernel.h");
    cc::Build::new()
        .file("csrc/baq_simd.c")
        .include(include)
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
}
