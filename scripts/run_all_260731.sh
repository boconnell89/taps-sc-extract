#!/usr/bin/env bash
set -euo pipefail

REF="/mnt/e/refs/mm10/mm10.fa"
BIN="/mnt/e/sciMET_TAPS/taps-sc-extract-rs/target/release/taps-sc-extract-rs"
DATA_DIR="/mnt/e/sciMET_TAPS/260731"

BAMS=(
  "taps_scTn5_sp2.srt.bam"
  "5base_scTn5_sp2.srt.bam"
  "5base_ucbTn5_sp2.srt.bam"
  "taps_ucbTn5_sp2.srt.bam"
)

echo "======================================================================"
echo "Starting batch extraction for ${#BAMS[@]} BAM files in ${DATA_DIR}"
echo "Reference: ${REF}"
echo "Start time: $(date)"
echo "======================================================================"

for bam_file in "${BAMS[@]}"; do
  sample_name="${bam_file%.srt.bam}"
  bam_path="${DATA_DIR}/${bam_file}"
  out_dir="${DATA_DIR}/${sample_name}"
  log_path="${DATA_DIR}/${sample_name}.extract.log"

  echo ""
  echo "----------------------------------------------------------------------"
  echo "Processing sample: ${sample_name}"
  echo "Input BAM:         ${bam_path}"
  echo "Output Directory:  ${out_dir}"
  echo "Log:               ${log_path}"
  echo "Start:             $(date)"
  echo "----------------------------------------------------------------------"

  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"

  /usr/bin/time -v "${BIN}" extract \
    -b "${bam_path}" \
    -f "${REF}" \
    -o "${out_dir}" \
    -t 24 \
    --shards 16 \
    --decomp-threads 0 \
    --memory-mode stream \
    --compression gzip 2>&1 | tee "${log_path}"

  echo "Finished ${sample_name} at $(date)"
  echo "Summary of ${sample_name}:"
  tail -n 25 "${log_path}"
done

echo ""
echo "======================================================================"
echo "All BAM extractions completed successfully at $(date)!"
echo "======================================================================"
