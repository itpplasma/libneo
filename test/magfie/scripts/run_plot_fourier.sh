#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: run_plot_fourier.sh <repo_root> <build_dir> <case>" >&2
  exit 2
fi

repo_root=$1
build_dir=$2
case_name=$3

vacfield_exe="${build_dir}/vacfield.x"
plot_script="${repo_root}/python/scripts/plot_biotsavart_fourier_cli.py"
data_dir="${repo_root}/test/magfie/test_data"
output_dir="${build_dir}/test/magfie/plot_fourier/${case_name}"
python_exec="${PYTHON_EXECUTABLE:-python3}"

mkdir -p "${output_dir}"

case "${case_name}" in
  aug_lowres)
    coil_files=("${data_dir}/aug_bl.dat" "${data_dir}/aug_bu.dat")
    grid_file="${data_dir}/vacfield_AUG_lowres.in"
    currents_file="${data_dir}/aug_currents.txt"
    ntor=2
    ;;
  axisymmetric)
    coil_files=("${data_dir}/axisymmetric_coil.dat")
    grid_file="${data_dir}/vacfield_axisymmetric.in"
    currents_file="${data_dir}/axisymmetric_currents.txt"
    ntor=0
    ;;
  single_coil)
    coil_files=("${data_dir}/single_coil.dat")
    grid_file="${data_dir}/vacfield_single_coil.in"
    currents_file="${data_dir}/single_coil_currents.txt"
    ntor=0
    ;;
  *)
    echo "Unknown case ${case_name}" >&2
    exit 3
    ;;
 esac

ref_file="${output_dir}/${case_name}_reference.h5"
vec_file="${output_dir}/${case_name}_vector.nc"
log_file="${output_dir}/${case_name}_plot.log"
sum_png="${output_dir}/${case_name}_sum.png"
per_coil_png="${output_dir}/${case_name}_per_coil.png"
deriv_png="${output_dir}/${case_name}_deriv_diff.png"

# Generate reference Fourier data
"${vacfield_exe}" GPEC "${#coil_files[@]}" "${coil_files[@]}" Fourier "${grid_file}" "${ref_file}"

# Generate vector potential data
"${vacfield_exe}" GPEC "${#coil_files[@]}" "${coil_files[@]}" vector_potential "${grid_file}" "${vec_file}"

# Run plotting script non-interactively
export MPLBACKEND=Agg
export PYTHONPATH="${repo_root}/python:${build_dir}"
export LD_LIBRARY_PATH="${build_dir}:${LD_LIBRARY_PATH-}"
"${python_exec}" "${plot_script}" \
    "${ref_file}" "${vec_file}" \
    --currents "${currents_file}" \
    --coil-files "${coil_files[@]}" \
    --ntor "${ntor}" \
    --per-coil-output "${per_coil_png}" \
    --sum-output "${sum_png}" \
    --deriv-diff-output "${deriv_png}" \
    --prefactor 0.1 \
    --dpi 150 &> "${log_file}"
