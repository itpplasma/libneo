# Draft: GPU spline many-point benchmarks (1D/2D/3D)

This folder is an isolated draft benchmark for the GPU-ready API work discussed in issue
`itpplasma/libneo#199`. It does **not** modify libneo production code; it reuses the
existing batch spline coefficient layouts (`BatchSplineData{1D,2D,3D}`) and benchmarks
four alternative GPU stacks for evaluating many points per call.

This is intentionally targeting the missing API shape in libneo today:
- existing: batch-over-quantities at a single point (`evaluate_batch_splines_*`)
- needed: batch-over-points (many points per call): `y(nq, npts)`

## What is benchmarked

Given one `BatchSplineData*D` holding `nq` quantities on a shared grid, evaluate all
quantities at `npts` points:

- Inputs:
  - 1D: `coeff(nq, 0:order, nx)`, `x(1:npts)`
  - 2D: `coeff(nq, 0:order1, 0:order2, nx, ny)`, `x(2, npts)`
  - 3D: `coeff(nq, 0:order1, 0:order2, 0:order3, nx, ny, nz)`, `x(3, npts)`
- Output layout (flat, equivalent to `y(nq, npts)`):
  - `y((ipt-1)*nq + iq)` for `ipt=1..npts`, `iq=1..nq`

All implementations are in double precision (`real(dp)` / CUDA `double`).

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
/tmp/libneo_gpu_spline_1d_many_build/bench_spline2d_many
/tmp/libneo_gpu_spline_1d_many_build/bench_spline3d_many
```

## Captured runs (evidence)

On this machine (RTX 5060 Ti, driver 590.48.01, CUDA 13.1), running:

- 1D: `order=5`, `num_points=2048`, `num_quantities=8`, `npts=2000000`, `niter=20`, `periodic=T`
- 2D: `order=[5,5]`, `num_points=[256,256]`, `num_quantities=8`, `npts=500000`, `niter=10`,
  `periodic=[T,T]`
- 3D: `order=[5,3,3]`, `num_points=[48,32,32]`, `num_quantities=8`, `npts=200000`, `niter=6`,
  `periodic=[T,T,T]`

produced the following best times (Fortran compiled with `-O3` in this draft CMake project):

From `/tmp/libneo_gpu_spline_many_bench_1d_omp_numteams_2025-12-21.log`:
- CPU: `best_s 0.041226` → `4.85e7 pts/s`
- OpenACC: `best_s 0.002578` → `7.76e8 pts/s`
- OpenMP target: `best_s 0.028726` → `6.96e7 pts/s`
- CUDA Fortran: `best_s 0.002572` → `7.78e8 pts/s`
- CUDA C: `best_s 0.002542` → `7.87e8 pts/s`

From `/tmp/libneo_gpu_spline_many_bench_2d_omp_numteams_2025-12-21.log`:
- CPU: `best_s 0.056710` → `8.82e6 pts/s`
- OpenACC: `best_s 0.000974` → `5.13e8 pts/s`
- OpenMP target: `best_s 0.034965` → `1.43e7 pts/s`
- CUDA Fortran: `best_s 0.000952` → `5.25e8 pts/s`
- CUDA C: `best_s 0.000976` → `5.12e8 pts/s`

From `/tmp/libneo_gpu_spline_many_bench_3d_omp_numteams_2025-12-21.log`:
- CPU: `best_s 0.064090` → `3.12e6 pts/s`
- OpenACC: `best_s 0.000979` → `2.04e8 pts/s`
- OpenMP target: `best_s 0.058352` → `3.43e6 pts/s`
- CUDA Fortran: `best_s 0.000960` → `2.08e8 pts/s`
- CUDA C: `best_s 0.000967` → `2.07e8 pts/s`

All variants reported `max_abs_diff 0.0` versus the CPU reference in these runs.

## Notes on OpenMP target performance

On this machine, OpenMP target offload via NVHPC (`-mp=gpu`) is consistently far slower than
OpenACC / CUDA for these kernels, and is only slightly faster than the CPU for 3D.

This draft includes a `num_teams(...)` hint in the OpenMP target kernel launch. It provides
only marginal improvement here; the remaining gap appears to be dominated by NVHPC OpenMP
GPU code generation/runtime overhead rather than simple launch configuration.
