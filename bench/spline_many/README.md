# GPU spline many-point benchmarks (1D/2D/3D)

This folder is an isolated benchmark build that reuses libneo interpolate sources and
benchmarks one shared Fortran kernel compiled in different modes (CPU / OpenACC).

This is intentionally benchmarking the many-point API shape in libneo:
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

The benchmark binaries also time spline *construction* (building coefficients) in three
variants:
- `build`: legacy CPU construction (`construct_batch_splines_*`)
- `build_lines`: line-based CPU construction (`construct_batch_splines_*_lines`)
- `build_device_*`: device-resident construction (`construct_batch_splines_*_resident_device`)

## Implementations (alternatives)

- One shared Fortran kernel compiled in different modes:
  - CPU: plain Fortran loop (same routine)
  - OpenACC: offload via NVHPC `nvfortran -acc`

Sources live in `bench/spline_many/src/`.

## Build and run

For NVHPC GPU offload (`-acc`), NVHPC needs a CUDA toolkit; on this machine it is provided
at `/opt/cuda`:

```bash
export NVHPC_CUDA_HOME=/opt/cuda
```

Configure and build (out-of-tree; example uses `/tmp`):

```bash
cmake -S bench/spline_many -B /tmp/libneo_bench_spline_many_build -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin/nvfortran
cmake --build /tmp/libneo_bench_spline_many_build -j
```

Run:

```bash
/tmp/libneo_bench_spline_many_build/bench_spline1d_many
/tmp/libneo_bench_spline_many_build/bench_spline2d_many
/tmp/libneo_bench_spline_many_build/bench_spline3d_many
```

### NVHPC device selection

```bash
ACC_DEVICE_TYPE=host   /tmp/libneo_bench_spline_many_build/bench_spline1d_many
ACC_DEVICE_TYPE=nvidia /tmp/libneo_bench_spline_many_build/bench_spline1d_many
```

## Automated benchmark matrix

Run the matrix (CPU, OpenACC host/GPU; nvfortran; 1D/2D/3D):

```bash
./bench/spline_many/run_benchmarks.sh
```

This writes one log per run under `/tmp` with the prefix `libneo_gpu_spline_`.

## Captured runs (evidence)

On this machine (RTX 5060 Ti, driver 590.48.01), running:

- 1D: `order=5`, `num_points=2048`, `num_quantities=8`, `npts=2000000`, `niter=20`, `periodic=T`
- 2D: `order=[5,5]`, `num_points=[256,256]`, `num_quantities=8`, `npts=500000`, `niter=10`,
  `periodic=[T,T]`
- 3D: `order=[5,3,3]`, `num_points=[48,32,32]`, `num_quantities=8`, `npts=200000`, `niter=6`,
  `periodic=[T,T,T]`

These are historical single runs; prefer `run_benchmarks.sh` for the current matrix and keep the
`/tmp/libneo_gpu_spline_*_<date>.log` logs as the source of truth.
