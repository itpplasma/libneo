# Draft: GPU spline 1D many-point benchmark

This folder is an isolated draft benchmark for the GPU-ready API work discussed in issue
`itpplasma/libneo#199`. It does **not** modify libneo production code; it reuses the
existing 1D batch spline coefficient layout (`BatchSplineData1D`) and benchmarks
four alternative GPU stacks for evaluating many points per call.

## What is benchmarked

Given one `BatchSplineData1D` holding `nq` quantities on a shared 1D grid, evaluate
all quantities at `npts` points:

- Inputs:
  - `coeff(nq, 0:order, nx)` from `BatchSplineData1D`
  - `x(1:npts)`
- Output layout (flat, equivalent to `y(nq, npts)`):
  - `y((ipt-1)*nq + iq)` for `ipt=1..npts`, `iq=1..nq`

## Implementations (alternatives)

- CPU reference (Fortran)
- OpenACC (NVHPC `nvfortran -acc`)
- OpenMP target offload (NVHPC `nvfortran -mp=gpu`)
- CUDA Fortran (NVHPC `nvfortran -cuda`)
- CUDA C kernel (NVCC) + Fortran wrapper (`iso_c_binding`)

Sources live in `draft/gpu_spline_1d_many/src/`.

## Build and run

This draft project expects NVHPC compilers and a CUDA toolkit available at `/opt/cuda`.
NVHPC needs the CUDA path:

```bash
export NVHPC_CUDA_HOME=/opt/cuda
```

Configure and build (out-of-tree; example uses `/tmp`):

```bash
cmake -S draft/gpu_spline_1d_many -B /tmp/libneo_gpu_spline_1d_many_build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvfortran \
  -DCMAKE_C_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvc
cmake --build /tmp/libneo_gpu_spline_1d_many_build -j
```

Run:

```bash
/tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
```

## One captured run (evidence)

On this machine (RTX 5060 Ti, driver 590.48.01, CUDA 13.1), running:

- `order=5`, `num_points=2048`, `num_quantities=8`
- `npts=2000000`, `niter=20`, `periodic=T`

produced (from `/tmp/libneo_gpu_spline1d_many_bench_4.log`):

- CPU: `best_s 0.226729` → `8.82e6 pts/s`
- OpenACC: `best_s 0.012501` → `1.60e8 pts/s`
- OpenMP target: `best_s 0.028722` → `6.96e7 pts/s`
- CUDA Fortran: `best_s 0.002480` → `8.06e8 pts/s`
- CUDA C: `best_s 0.002441` → `8.19e8 pts/s`

All variants reported `max_abs_diff 0.0` versus the CPU reference for this run.
