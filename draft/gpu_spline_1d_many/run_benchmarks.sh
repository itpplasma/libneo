#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
date_tag="$(date +%F)"

nvfortran_bin="${NVFORTRAN_BIN:-nvfortran}"
gfortran_bin="${GFORTRAN_BIN:-/opt/gcc16/bin/gfortran}"

cmake_gen="${CMAKE_GENERATOR:-Ninja}"

log_dir="/tmp"

if [[ -d /opt/cuda ]]; then
  export NVHPC_CUDA_HOME="${NVHPC_CUDA_HOME:-/opt/cuda}"
  export NVCOMPILER_CUDA_HOME="${NVCOMPILER_CUDA_HOME:-${NVHPC_CUDA_HOME}}"
fi

gomp_debug="${LIBNEO_DRAFT_GOMP_DEBUG:-0}"

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
  local kind="$2" # cpu|openacc|openmp
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

extract_setup_s() {
  local log_path="$1"
  local kind="$2" # openacc|openmp
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

### nvfortran
if command -v "${nvfortran_bin}" >/dev/null 2>&1; then
  nv_tag="nvfortran"

  # CPU-only (no OpenACC)
  nv_cpu_build="$(build_dir_for "${nv_tag}" cpu)"
  rm -rf "${nv_cpu_build}"
  cmake_configure_build "${nvfortran_bin}" "${nv_cpu_build}" \
    -DDRAFT_ENABLE_OPENACC=OFF

  for dim in 1d 2d 3d; do
    exe="${nv_cpu_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${nv_tag}_cpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" "${exe}"
    results["${nv_tag},build,${dim}"]="$(extract_build_grid_pts_per_s "${log}")"
    results["${nv_tag},cpu,${dim}"]="$(extract_pts_per_s "${log}" cpu)"
  done

  # OpenACC build (same binary; runtime selects device)
  nv_acc_build="$(build_dir_for "${nv_tag}" openacc)"
  rm -rf "${nv_acc_build}"
  cmake_configure_build "${nvfortran_bin}" "${nv_acc_build}" \
    -DDRAFT_ENABLE_OPENACC=ON

  for dim in 1d 2d 3d; do
    exe="${nv_acc_build}/bench_spline${dim}_many"

    log_host="${log_dir}/libneo_gpu_spline_${nv_tag}_openacc_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_host}" env ACC_DEVICE_TYPE=host "${exe}"
    results["${nv_tag},openacc_setup_host,${dim}"]="$(extract_setup_s "${log_host}" openacc)"
    results["${nv_tag},openacc_host,${dim}"]="$(extract_pts_per_s "${log_host}" openacc)"

    log_gpu="${log_dir}/libneo_gpu_spline_${nv_tag}_openacc_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_gpu}" env ACC_DEVICE_TYPE=nvidia "${exe}"
    results["${nv_tag},openacc_setup_gpu,${dim}"]="$(extract_setup_s "${log_gpu}" openacc)"
    results["${nv_tag},openacc_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" openacc)"
  done

  # OpenMP target build (same binary; runtime selects device)
  nv_omp_build="$(build_dir_for "${nv_tag}" openmp)"
  rm -rf "${nv_omp_build}"
  cmake_configure_build "${nvfortran_bin}" "${nv_omp_build}" \
    -DDRAFT_ENABLE_OPENACC=OFF \
    -DDRAFT_ENABLE_OPENMP=ON

  for dim in 1d 2d 3d; do
    exe="${nv_omp_build}/bench_spline${dim}_many"

    log_host="${log_dir}/libneo_gpu_spline_${nv_tag}_openmp_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_host}" env OMP_TARGET_OFFLOAD=DISABLED OMP_TEAMS_THREAD_LIMIT=256 "${exe}"
    results["${nv_tag},openmp_setup_host,${dim}"]="$(extract_setup_s "${log_host}" openmp)"
    results["${nv_tag},openmp_host,${dim}"]="$(extract_pts_per_s "${log_host}" openmp)"

    log_gpu="${log_dir}/libneo_gpu_spline_${nv_tag}_openmp_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log_gpu}" env OMP_TARGET_OFFLOAD=MANDATORY OMP_TEAMS_THREAD_LIMIT=256 "${exe}"
    results["${nv_tag},openmp_setup_gpu,${dim}"]="$(extract_setup_s "${log_gpu}" openmp)"
    results["${nv_tag},openmp_gpu,${dim}"]="$(extract_pts_per_s "${log_gpu}" openmp)"
  done
else
  echo "nvfortran not found (set NVFORTRAN_BIN=... to override), skipping NVHPC runs." >&2
fi

### gfortran (GCC16)
if [[ -x "${gfortran_bin}" ]]; then
  g_tag="gfortran"

  # CPU-only
  g_cpu_build="$(build_dir_for "${g_tag}" cpu)"
  rm -rf "${g_cpu_build}"
  cmake_configure_build "${gfortran_bin}" "${g_cpu_build}" \
    -DDRAFT_ENABLE_OPENACC=OFF

  for dim in 1d 2d 3d; do
    exe="${g_cpu_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${g_tag}_cpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" "${exe}"
    results["${g_tag},build,${dim}"]="$(extract_build_grid_pts_per_s "${log}")"
    results["${g_tag},cpu,${dim}"]="$(extract_pts_per_s "${log}" cpu)"
  done

  # OpenACC host (no offload flag)
  g_acc_host_build="$(build_dir_for "${g_tag}" openacc_host)"
  rm -rf "${g_acc_host_build}"
  cmake_configure_build "${gfortran_bin}" "${g_acc_host_build}" \
    -DDRAFT_ENABLE_OPENACC=ON \
    -DDRAFT_OPENACC_OFFLOAD=none

  for dim in 1d 2d 3d; do
    exe="${g_acc_host_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${g_tag}_openacc_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" env LD_LIBRARY_PATH=/opt/gcc16/lib64 ACC_DEVICE_TYPE=host "${exe}"
    results["${g_tag},openacc_setup_host,${dim}"]="$(extract_setup_s "${log}" openacc)"
    results["${g_tag},openacc_host,${dim}"]="$(extract_pts_per_s "${log}" openacc)"
  done

  # OpenACC GPU (nvptx offload)
  g_acc_gpu_build="$(build_dir_for "${g_tag}" openacc_gpu)"
  rm -rf "${g_acc_gpu_build}"
  cmake_configure_build "${gfortran_bin}" "${g_acc_gpu_build}" \
    -DDRAFT_ENABLE_OPENACC=ON \
    -DDRAFT_OPENACC_OFFLOAD=nvptx

  for dim in 1d 2d 3d; do
    exe="${g_acc_gpu_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${g_tag}_openacc_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" env LD_LIBRARY_PATH=/opt/gcc16/lib64 ACC_DEVICE_TYPE=nvidia GOMP_DEBUG="${gomp_debug}" "${exe}"
    results["${g_tag},openacc_setup_gpu,${dim}"]="$(extract_setup_s "${log}" openacc)"
    results["${g_tag},openacc_gpu,${dim}"]="$(extract_pts_per_s "${log}" openacc)"
  done

  # OpenMP host (no offload flag)
  g_omp_host_build="$(build_dir_for "${g_tag}" openmp_host)"
  rm -rf "${g_omp_host_build}"
  cmake_configure_build "${gfortran_bin}" "${g_omp_host_build}" \
    -DDRAFT_ENABLE_OPENACC=OFF \
    -DDRAFT_ENABLE_OPENMP=ON \
    -DDRAFT_OPENMP_OFFLOAD=none

  for dim in 1d 2d 3d; do
    exe="${g_omp_host_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${g_tag}_openmp_host_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" env LD_LIBRARY_PATH=/opt/gcc16/lib64 OMP_TARGET_OFFLOAD=DISABLED OMP_TEAMS_THREAD_LIMIT=256 "${exe}"
    results["${g_tag},openmp_setup_host,${dim}"]="$(extract_setup_s "${log}" openmp)"
    results["${g_tag},openmp_host,${dim}"]="$(extract_pts_per_s "${log}" openmp)"
  done

  # OpenMP GPU (nvptx offload)
  g_omp_gpu_build="$(build_dir_for "${g_tag}" openmp_gpu)"
  rm -rf "${g_omp_gpu_build}"
  cmake_configure_build "${gfortran_bin}" "${g_omp_gpu_build}" \
    -DDRAFT_ENABLE_OPENACC=OFF \
    -DDRAFT_ENABLE_OPENMP=ON \
    -DDRAFT_OPENMP_OFFLOAD=nvptx

  for dim in 1d 2d 3d; do
    exe="${g_omp_gpu_build}/bench_spline${dim}_many"
    log="${log_dir}/libneo_gpu_spline_${g_tag}_openmp_gpu_${dim}_${date_tag}.log"
    run_one "${exe}" "${log}" env LD_LIBRARY_PATH=/opt/gcc16/lib64 OMP_TARGET_OFFLOAD=MANDATORY OMP_TEAMS_THREAD_LIMIT=256 GOMP_DEBUG="${gomp_debug}" "${exe}"
    results["${g_tag},openmp_setup_gpu,${dim}"]="$(extract_setup_s "${log}" openmp)"
    results["${g_tag},openmp_gpu,${dim}"]="$(extract_pts_per_s "${log}" openmp)"
  done
else
  echo "gfortran not found at ${gfortran_bin} (set GFORTRAN_BIN=... to override), skipping GNU runs." >&2
fi

echo
echo "Summary (grid_pts_per_s, pts_per_s, setup_s)"
printf "%-10s %-13s %-4s %s\n" "compiler" "metric" "dim" "value"
for dim in 1d 2d 3d; do
  for compiler in nvfortran gfortran; do
    for metric in build cpu \
      openacc_setup_host openacc_host openacc_setup_gpu openacc_gpu \
      openmp_setup_host openmp_host openmp_setup_gpu openmp_gpu; do
      key="${compiler},${metric},${dim}"
      if [[ -n "${results[$key]:-}" ]]; then
        printf "%-10s %-13s %-4s %s\n" "${compiler}" "${metric}" "${dim}" "${results[$key]}"
      fi
    done
  done
done
