#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
date_tag="$(date +%F)"

nvfortran_bin="${NVFORTRAN_BIN:-nvfortran}"

cmake_gen="${CMAKE_GENERATOR:-Ninja}"

log_dir="/tmp"

if [[ -d /opt/cuda ]]; then
  export NVHPC_CUDA_HOME="${NVHPC_CUDA_HOME:-/opt/cuda}"
  export NVCOMPILER_CUDA_HOME="${NVCOMPILER_CUDA_HOME:-${NVHPC_CUDA_HOME}}"
fi

gomp_debug="${LIBNEO_BENCH_GOMP_DEBUG:-0}"
fast="${LIBNEO_BENCH_FAST:-1}"

bench_env_common=()
bench_env_dim=()
bench_env_set_for_dim() {
  local dim="$1"
  bench_env_common=()
  bench_env_dim=()

  if [[ "${fast}" == "1" ]]; then
    bench_env_common+=(
      "LIBNEO_BENCH_NQ=${LIBNEO_BENCH_NQ:-2}"
      "LIBNEO_BENCH_NBUILD=${LIBNEO_BENCH_NBUILD:-1}"
      "LIBNEO_BENCH_NBUILD_REPEAT=${LIBNEO_BENCH_NBUILD_REPEAT:-1}"
      "LIBNEO_BENCH_NBUILD_RESIDENT=${LIBNEO_BENCH_NBUILD_RESIDENT:-1}"
      "LIBNEO_BENCH_NEVAL_REPEAT=${LIBNEO_BENCH_NEVAL_REPEAT:-1}"
      "LIBNEO_BENCH_NITER=${LIBNEO_BENCH_NITER:-3}"
      "LIBNEO_BENCH_OLD_NPTS=${LIBNEO_BENCH_OLD_NPTS:-5000}"
    )
    case "${dim}" in
      1d)
        bench_env_dim+=(
          "LIBNEO_BENCH_NUM_POINTS=${LIBNEO_BENCH_NUM_POINTS:-512}"
          "LIBNEO_BENCH_NPTS=${LIBNEO_BENCH_NPTS_1D:-200000}"
        )
        ;;
      2d)
        bench_env_dim+=(
          "LIBNEO_BENCH_N1=${LIBNEO_BENCH_N1_2D:-64}"
          "LIBNEO_BENCH_N2=${LIBNEO_BENCH_N2_2D:-64}"
          "LIBNEO_BENCH_NPTS=${LIBNEO_BENCH_NPTS_2D:-50000}"
        )
        ;;
      3d)
        bench_env_dim+=(
          "LIBNEO_BENCH_N1=${LIBNEO_BENCH_N1_3D:-24}"
          "LIBNEO_BENCH_N2=${LIBNEO_BENCH_N2_3D:-16}"
          "LIBNEO_BENCH_N3=${LIBNEO_BENCH_N3_3D:-16}"
          "LIBNEO_BENCH_NPTS=${LIBNEO_BENCH_NPTS_3D:-30000}"
        )
        ;;
      *)
        echo "bench_env_set_for_dim: unknown dim '${dim}'" >&2
        exit 2
        ;;
    esac
  fi
}

build_dir_for() {
  local compiler_tag="$1"
  local variant="$2"
  echo "/tmp/libneo_gpu_spline_${compiler_tag}_${variant}_build"
}

run_one() {
  local exe_path="$1"
  local log_path="$2"
  shift 2
  (
    echo "==== $(date -Is) ===="
    echo "exe: ${exe_path}"
    echo "cmd: $*"
    "$@"
  ) 2>&1 | tee "${log_path}"
}

cmake_configure_build() {
  local compiler="$1"
  local build_dir="$2"
  shift 2
  cmake -S "${root_dir}" -B "${build_dir}" -G "${cmake_gen}" \
    -DCMAKE_Fortran_COMPILER="${compiler}" \
    "$@"
  cmake --build "${build_dir}" -j
}

extract_pts_per_s() {
  local log_path="$1"
  local kind="$2" # cpu|openacc|public
  awk -v kind="${kind}" '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == kind) state = 1;
        else if (state == 1 && $i == "pts_per_s") state = 2;
        else if (state == 2) { print $i; exit; }
      }
    }
  ' "${log_path}"
}

extract_warmup_s() {
  local log_path="$1"
  local kind="$2" # public
  awk -v kind="${kind}" '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == kind) state = 1;
        else if (state == 1 && $i == "warmup_s") state = 2;
        else if (state == 2) { print $(i+0); exit; }
      }
    }
  ' "${log_path}"
}

extract_build_grid_pts_per_s() {
  local log_path="$1"
  awk '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == "build") state = 1;
        else if (state == 1 && $i == "grid_pts_per_s") state = 2;
        else if (state == 2) { print $i; exit; }
      }
    }
  ' "${log_path}"
}

extract_build_lines_grid_pts_per_s() {
  local log_path="$1"
  awk '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == "build_lines") state = 1;
        else if (state == 1 && $i == "grid_pts_per_s") state = 2;
        else if (state == 2) { print $i; exit; }
      }
    }
  ' "${log_path}"
}

extract_build_resident_grid_pts_per_s() {
  local log_path="$1"
  awk '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == "build_resident") state = 1;
        else if (state == 1 && $i == "grid_pts_per_s") state = 2;
        else if (state == 2) { print $i; exit; }
      }
    }
  ' "${log_path}"
}

extract_build_device_grid_pts_per_s() {
  local log_path="$1"
  local tag="$2" # host|gpu
  awk -v tag="${tag}" '
    {
      want = "build_device_" tag
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == want) state = 1;
        else if (state == 1 && $i == "grid_pts_per_s") state = 2;
        else if (state == 2) { print $i; exit; }
      }
    }
  ' "${log_path}"
}

extract_setup_s() {
  local log_path="$1"
  local kind="$2" # openacc
  awk -v kind="${kind}" '
    {
      for (i = 1; i <= NF; i++) {
        if (state == 0 && $i == kind) state = 1;
        else if (state == 1 && $i == "setup_s") state = 2;
        else if (state == 2) { print $(i+0); exit; }
      }
    }
  ' "${log_path}"
}

echo "Building + running benchmarks (1D/2D/3D) into ${log_dir} (date tag ${date_tag})"

declare -A results

gfortran_bin="${GFORTRAN_BIN:-gfortran}"

# Derive GCC lib path from compiler location for LD_LIBRARY_PATH
gfortran_libdir=""
if [[ -n "${GFORTRAN_BIN:-}" ]]; then
  gfortran_prefix="$(dirname "$(dirname "${GFORTRAN_BIN}")")"
  if [[ -d "${gfortran_prefix}/lib64" ]]; then
    gfortran_libdir="${gfortran_prefix}/lib64"
  fi
fi

### gfortran
if command -v "${gfortran_bin}" >/dev/null 2>&1; then
  gf_tag="gfortran"

  # OpenACC build for gfortran
  gf_acc_build="$(build_dir_for "${gf_tag}" openacc)"
  if [[ "${fast}" != "1" ]]; then
    rm -rf "${gf_acc_build}"
  fi
  cmake_configure_build "${gfortran_bin}" "${gf_acc_build}" \
    -DBENCH_ENABLE_OPENACC=ON

  # Build gfortran LD_LIBRARY_PATH env array
  gf_ld_env=()
  if [[ -n "${gfortran_libdir}" ]]; then
    gf_ld_env+=("LD_LIBRARY_PATH=${gfortran_libdir}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}")
  fi

  for dim in 1d 2d 3d; do
    bench_env_set_for_dim "${dim}"
    exe="${gf_acc_build}/bench_spline${dim}_many"

    log_host="${log_dir}/libneo_gpu_spline_${gf_tag}_openacc_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_host}" env "${gf_ld_env[@]}" "${bench_env_common[@]}" "${bench_env_dim[@]}" ACC_DEVICE_TYPE=host "${exe}"
    results["${gf_tag},build,${dim}"]="$(extract_build_grid_pts_per_s "${log_host}")"
    results["${gf_tag},build_lines,${dim}"]="$(extract_build_lines_grid_pts_per_s "${log_host}")"
    results["${gf_tag},build_resident_host,${dim}"]="$(extract_build_resident_grid_pts_per_s "${log_host}")"
    results["${gf_tag},build_device_host,${dim}"]="$(extract_build_device_grid_pts_per_s "${log_host}" host)"
    results["${gf_tag},openacc_setup_host,${dim}"]="$(extract_setup_s "${log_host}" openacc)"
    results["${gf_tag},openacc_host,${dim}"]="$(extract_pts_per_s "${log_host}" openacc)"
    results["${gf_tag},old_public_host,${dim}"]="$(extract_pts_per_s "${log_host}" old_public)"
    results["${gf_tag},public_warmup_host,${dim}"]="$(extract_warmup_s "${log_host}" public)"
    results["${gf_tag},public_host,${dim}"]="$(extract_pts_per_s "${log_host}" public)"
    results["${gf_tag},cpu,${dim}"]="$(extract_pts_per_s "${log_host}" cpu)"

    log_gpu="${log_dir}/libneo_gpu_spline_${gf_tag}_openacc_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_gpu}" env "${gf_ld_env[@]}" "${bench_env_common[@]}" "${bench_env_dim[@]}" ACC_DEVICE_TYPE=nvidia "${exe}"
    results["${gf_tag},build_resident_gpu,${dim}"]="$(extract_build_resident_grid_pts_per_s "${log_gpu}")"
    results["${gf_tag},build_device_gpu,${dim}"]="$(extract_build_device_grid_pts_per_s "${log_gpu}" gpu)"
    results["${gf_tag},openacc_setup_gpu,${dim}"]="$(extract_setup_s "${log_gpu}" openacc)"
    results["${gf_tag},openacc_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" openacc)"
    results["${gf_tag},old_public_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" old_public)"
    results["${gf_tag},public_warmup_gpu,${dim}"]="$(extract_warmup_s "${log_gpu}" public)"
    results["${gf_tag},public_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" public)"
  done
else
  echo "gfortran not found (set GFORTRAN_BIN=... to override), skipping." >&2
fi

### nvfortran
if command -v "${nvfortran_bin}" >/dev/null 2>&1; then
  nv_tag="nvfortran"

  # One OpenACC build; runtime selects host vs GPU via ACC_DEVICE_TYPE.
  nv_acc_build="$(build_dir_for "${nv_tag}" openacc)"
  if [[ "${fast}" != "1" ]]; then
    rm -rf "${nv_acc_build}"
  fi
  cmake_configure_build "${nvfortran_bin}" "${nv_acc_build}" \
    -DBENCH_ENABLE_OPENACC=ON

  for dim in 1d 2d 3d; do
    bench_env_set_for_dim "${dim}"
    exe="${nv_acc_build}/bench_spline${dim}_many"

    log_host="${log_dir}/libneo_gpu_spline_${nv_tag}_openacc_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_host}" env "${bench_env_common[@]}" "${bench_env_dim[@]}" ACC_DEVICE_TYPE=host "${exe}"
    results["${nv_tag},build,${dim}"]="$(extract_build_grid_pts_per_s "${log_host}")"
    results["${nv_tag},build_lines,${dim}"]="$(extract_build_lines_grid_pts_per_s "${log_host}")"
    results["${nv_tag},build_resident_host,${dim}"]="$(extract_build_resident_grid_pts_per_s "${log_host}")"
    results["${nv_tag},build_device_host,${dim}"]="$(extract_build_device_grid_pts_per_s "${log_host}" host)"
    results["${nv_tag},openacc_setup_host,${dim}"]="$(extract_setup_s "${log_host}" openacc)"
    results["${nv_tag},openacc_host,${dim}"]="$(extract_pts_per_s "${log_host}" openacc)"
    results["${nv_tag},old_public_host,${dim}"]="$(extract_pts_per_s "${log_host}" old_public)"
    results["${nv_tag},public_warmup_host,${dim}"]="$(extract_warmup_s "${log_host}" public)"
    results["${nv_tag},public_host,${dim}"]="$(extract_pts_per_s "${log_host}" public)"
    results["${nv_tag},cpu,${dim}"]="$(extract_pts_per_s "${log_host}" cpu)"

    log_gpu="${log_dir}/libneo_gpu_spline_${nv_tag}_openacc_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_gpu}" env "${bench_env_common[@]}" "${bench_env_dim[@]}" ACC_DEVICE_TYPE=nvidia "${exe}"
    results["${nv_tag},build,${dim}"]="$(extract_build_grid_pts_per_s "${log_gpu}")"
    results["${nv_tag},build_lines,${dim}"]="$(extract_build_lines_grid_pts_per_s "${log_gpu}")"
    results["${nv_tag},build_resident_gpu,${dim}"]="$(extract_build_resident_grid_pts_per_s "${log_gpu}")"
    results["${nv_tag},build_device_gpu,${dim}"]="$(extract_build_device_grid_pts_per_s "${log_gpu}" gpu)"
    results["${nv_tag},openacc_setup_gpu,${dim}"]="$(extract_setup_s "${log_gpu}" openacc)"
    results["${nv_tag},openacc_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" openacc)"
    results["${nv_tag},old_public_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" old_public)"
    results["${nv_tag},public_warmup_gpu,${dim}"]="$(extract_warmup_s "${log_gpu}" public)"
    results["${nv_tag},public_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" public)"
  done

else
  echo "nvfortran not found (set NVFORTRAN_BIN=... to override)." >&2
  exit 1
fi

echo
echo "Summary (grid_pts_per_s, pts_per_s, setup_s)"
printf "%-10s %-22s %-4s %s\n" "compiler" "metric" "dim" "value"
for dim in 1d 2d 3d; do
  for compiler in gfortran nvfortran; do
    for metric in build build_lines build_resident_host build_resident_gpu \
      build_device_host build_device_gpu cpu \
      old_public_host old_public_gpu \
      public_warmup_host public_host public_warmup_gpu public_gpu \
      openacc_setup_host openacc_host openacc_setup_gpu openacc_gpu; do
      key="${compiler},${metric},${dim}"
      if [[ -n "${results[$key]:-}" ]]; then
        printf "%-10s %-22s %-4s %s\n" "${compiler}" "${metric}" "${dim}" "${results[$key]}"
      fi
    done
  done
done
