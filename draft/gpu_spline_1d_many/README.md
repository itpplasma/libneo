# Draft: GPU spline many-point benchmarks (1D/2D/3D)

This folder is an isolated draft benchmark for the GPU-ready API work discussed in issue
`itpplasma/libneo#199`. It does **not** modify libneo production code; it reuses the
existing batch spline coefficient layouts (`BatchSplineData{1D,2D,3D}`) and benchmarks
an OpenACC implementation for evaluating many points per call.

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

All implementations are in double precision (`real(dp)`).

## Implementations (alternatives)

- One code path that runs on CPU and can be accelerated with OpenACC:
  - CPU: plain Fortran loop (same routine)
  - GPU: OpenACC offload via NVHPC `nvfortran -acc`

Sources live in `draft/gpu_spline_1d_many/src/`.

## Build and run

This draft project expects NVHPC `nvfortran`. For GPU offload (`-acc`), NVHPC needs a CUDA
toolkit; on this machine it is provided at `/opt/cuda`:

```bash
export NVHPC_CUDA_HOME=/opt/cuda
```

Configure and build (out-of-tree; example uses `/tmp`):

```bash
cmake -S draft/gpu_spline_1d_many -B /tmp/libneo_gpu_spline_1d_many_build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvfortran
cmake --build /tmp/libneo_gpu_spline_1d_many_build -j
```

Run:

```bash
/tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
/tmp/libneo_gpu_spline_1d_many_build/bench_spline2d_many
/tmp/libneo_gpu_spline_1d_many_build/bench_spline3d_many
```

## Captured runs (evidence)

On this machine (RTX 5060 Ti, driver 590.48.01), running:

- 1D: `order=5`, `num_points=2048`, `num_quantities=8`, `npts=2000000`, `niter=20`, `periodic=T`
- 2D: `order=[5,5]`, `num_points=[256,256]`, `num_quantities=8`, `npts=500000`, `niter=10`,
  `periodic=[T,T]`
- 3D: `order=[5,3,3]`, `num_points=[48,32,32]`, `num_quantities=8`, `npts=200000`, `niter=6`,
  `periodic=[T,T,T]`

produced the following best times (Fortran compiled with `-O3` in this draft CMake project):

From `/tmp/libneo_gpu_spline_openacc_only_1d_2025-12-21.log`:
- CPU: `best_s 0.041794` → `4.79e7 pts/s`
- OpenACC: `best_s 0.002580` → `7.75e8 pts/s`

From `/tmp/libneo_gpu_spline_openacc_only_2d_2025-12-21.log`:
- CPU: `best_s 0.060554` → `8.26e6 pts/s`
- OpenACC: `best_s 0.000974` → `5.13e8 pts/s`

From `/tmp/libneo_gpu_spline_openacc_only_3d_2025-12-21.log`:
- CPU: `best_s 0.063864` → `3.13e6 pts/s`
- OpenACC: `best_s 0.000980` → `2.04e8 pts/s`

All variants reported `max_abs_diff 0.0` versus the CPU reference in these runs.
