# Spline Many-Point Evaluation Benchmark Results

Benchmark results for GPU-accelerated spline construction and evaluation.

## System Configuration

- **OS**: Linux 6.18.2-2-cachyos (x86_64)
- **GPU**: NVIDIA GeForce RTX 5060 Ti
- **Driver**: 590.48.01
- **CUDA**: 13.0
- **Compilers**:
  - gfortran 16.0.0 (GCC trunk) with nvptx offloading
  - nvfortran 25.1 (NVIDIA HPC SDK)

## Benchmark Parameters

- Grid sizes: N1=48, N2=32, N3=32, NQ=8 quantities
- Evaluation points: 100,000
- Spline order: 5 (all dimensions)
- Boundary: periodic

## Evaluation Performance (100K points)

### GCC16 (gfortran 16.0.0 + nvptx OpenACC)

| Dimension | CPU (pts/s) | GPU (pts/s) | Speedup |
|-----------|-------------|-------------|---------|
| 1D | 47M | 683M | **14.5x** |
| 2D | 9.6M | 394M | **41x** |
| 3D | 1.8M | 86M | **47x** |

### nvfortran (25.1 + native GPU)

| Dimension | CPU (pts/s) | GPU (pts/s) | Speedup |
|-----------|-------------|-------------|---------|
| 1D | 52M | 712M | **13.7x** |
| 2D | 8.5M | 469M | **55x** |
| 3D | 1.4M | 88M | **64x** |

## Construction Performance (grid pts/s)

### GCC16

| Dimension | Host Build | Device Build (GPU) | Speedup |
|-----------|------------|-------------------|---------|
| 1D | 8.5M | 2.0M | 0.24x (slower) |
| 2D | 3.1M | 4.8M | **1.5x** |
| 3D | 304K | 1.04M | **3.4x** |

### nvfortran

| Dimension | Host Build | Device Build (GPU) | Speedup |
|-----------|------------|-------------------|---------|
| 1D | 9.3M | 1.5M | 0.16x (slower) |
| 2D | 2.7M | 4.3M | **1.6x** |
| 3D | 296K | 2.6M | **8.8x** |

## Numerical Accuracy

All results are at machine epsilon precision:

| Compiler | 1D | 2D | 3D |
|----------|----|----|-----|
| GCC16 | 2.8e-16 | 1.1e-16 | 5.6e-17 |
| nvfortran | 0 | 0 | 0 |

## Key Findings

1. **Evaluation speedup**: GPU is 14-64x faster than CPU for point evaluation
2. **3D construction speedup**: GPU is 3-9x faster than CPU for spline construction
3. **1D construction overhead**: GPU is slower for 1D due to kernel launch overhead dominating small workloads
4. **Compiler comparison**: nvfortran shows slightly better GPU performance, especially for 3D construction (8.8x vs 3.4x speedup)

## GCC OpenACC Workaround

The 3D GPU construction code includes a workaround for a GCC nvptx OpenACC bug
(order-dependent struct mapping issue). The workaround keeps allocations persistent
across repeated construction calls to avoid triggering the bug. See
`src/interpolate/batch_interpolate_3d.f90` for details.

## Running the Benchmarks

```bash
# GCC16 build
cmake -S bench/spline_many -B build_bench_gcc16 -G Ninja \
  -DCMAKE_Fortran_COMPILER=/opt/gcc16/bin/gfortran
cmake --build build_bench_gcc16 -j
LD_LIBRARY_PATH=/opt/gcc16/lib64 ./build_bench_gcc16/bench_spline3d_many

# nvfortran build
cmake -S bench/spline_many -B build_bench_nvfortran -G Ninja \
  -DCMAKE_Fortran_COMPILER=nvfortran
cmake --build build_bench_nvfortran -j
./build_bench_nvfortran/bench_spline3d_many
```

Environment variables for configuration:
- `LIBNEO_BENCH_N1`, `LIBNEO_BENCH_N2`, `LIBNEO_BENCH_N3`: Grid dimensions
- `LIBNEO_BENCH_NQ`: Number of quantities
- `LIBNEO_BENCH_NPTS`: Number of evaluation points
- `LIBNEO_BENCH_NITER`: Number of timing iterations
