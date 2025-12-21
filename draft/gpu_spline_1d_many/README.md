# Draft: GPU spline many-point benchmarks (1D/2D/3D)

This folder is an isolated draft benchmark for the GPU-ready API work discussed in issue
`itpplasma/libneo#199`. It does **not** modify libneo production code; it reuses the
existing batch spline coefficient layouts (`BatchSplineData{1D,2D,3D}`) and benchmarks
one shared Fortran kernel compiled in different modes (CPU / OpenACC / OpenMP target).

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

- One shared Fortran kernel compiled in different modes:
  - CPU: plain Fortran loop (same routine)
  - OpenACC: offload via NVHPC `nvfortran -acc` or GNU `gfortran -fopenacc`
  - OpenMP target: offload via NVHPC `nvfortran -mp=gpu` or GNU `gfortran -fopenmp`

Sources live in `draft/gpu_spline_1d_many/src/`.

## Build and run

This draft project supports NVHPC `nvfortran` and GNU `gfortran`.

For NVHPC GPU offload (`-acc`), NVHPC needs a CUDA toolkit; on this machine it is provided
at `/opt/cuda`:

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

### NVHPC device selection

```bash
ACC_DEVICE_TYPE=host   /tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
ACC_DEVICE_TYPE=nvidia /tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
```

### OpenMP target device selection

```bash
OMP_TARGET_OFFLOAD=DISABLED  /tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
OMP_TARGET_OFFLOAD=MANDATORY /tmp/libneo_gpu_spline_1d_many_build/bench_spline1d_many
```

## Automated benchmark matrix

Run the full matrix (CPU, OpenACC host/GPU, OpenMP host/GPU; nvfortran vs gfortran; 1D/2D/3D):

```bash
./draft/gpu_spline_1d_many/run_benchmarks.sh
```

This writes one log per run under `/tmp` with the prefix `libneo_gpu_spline_`.

## Captured runs (evidence)

On this machine (RTX 5060 Ti, driver 590.48.01), running:

- 1D: `order=5`, `num_points=2048`, `num_quantities=8`, `npts=2000000`, `niter=20`, `periodic=T`
- 2D: `order=[5,5]`, `num_points=[256,256]`, `num_quantities=8`, `npts=500000`, `niter=10`,
  `periodic=[T,T]`
- 3D: `order=[5,3,3]`, `num_points=[48,32,32]`, `num_quantities=8`, `npts=200000`, `niter=6`,
  `periodic=[T,T,T]`

produced the following best times (Fortran compiled with `-O3` in this draft CMake project):

From `/tmp/libneo_gpu_spline_matrix_with_openmp_2025-12-21_v6.log` (matrix run; `pts_per_s`):

```
compiler   mode          dim  pts_per_s
nvfortran  cpu           1d   48381634.33160761
nvfortran  openacc_host  1d   47614512.90353308
nvfortran  openacc_gpu   1d   807102502.0177627
nvfortran  openmp_host   1d   179018976.0114579
nvfortran  openmp_gpu    1d   803535556.4483846
gfortran   cpu           1d   46511627.905371234
gfortran   openacc_host  1d   44444444.446169116
gfortran   openacc_gpu   1d   1000000001.6152626
gfortran   openmp_host   1d   44444444.446169116
gfortran   openmp_gpu    1d   400000000.20954758
nvfortran  cpu           2d   8561790.441617154
nvfortran  openacc_host  2d   8388558.006878630
nvfortran  openacc_gpu   2d   538793103.4482527
nvfortran  openmp_host   2d   75999392.00486395
nvfortran  openmp_gpu    2d   531914893.6170547
gfortran   cpu           2d   10000000.000145519
gfortran   openacc_host  2d   10000000.000145519
gfortran   openacc_gpu   2d   499999999.89813662
gfortran   openmp_host   2d   33333333.334626839
gfortran   openmp_gpu    2d   250000000.40381566
nvfortran  cpu           3d   3123584.625716472
nvfortran  openacc_host  3d   3158509.815069252
nvfortran  openacc_gpu   3d   214592274.6781271
nvfortran  openmp_host   3d   33200531.20849917
nvfortran  openmp_gpu    3d   210526315.7894723
gfortran   cpu           3d   4166666.6667651953
gfortran   openacc_host  3d   4081632.6531388024
gfortran   openacc_gpu   3d   199999999.95925465
gfortran   openmp_host   3d   22222222.226677623
gfortran   openmp_gpu    3d   199999999.95925465
```

All variants reported `max_abs_diff` around machine epsilon versus the CPU reference.

These are historical single runs; prefer `run_benchmarks.sh` for the current matrix and keep the
`/tmp/libneo_gpu_spline_*_<date>.log` logs as the source of truth.
