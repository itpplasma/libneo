#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>

struct spline2d_handle {
    int order1;
    int order2;
    int n1;
    int n2;
    int nq;
    int periodic1;
    int periodic2;
    double x_min1;
    double x_min2;
    double h1;
    double h2;
    double* coeff_dev;
    double* x_dev;
    double* y_dev;
    int capacity_npts;
};

static __device__ __forceinline__ double wrap_periodic(double x, double x_min, double period) {
    double t = x - x_min;
    int k_wrap = static_cast<int>(t / period);
    double w = t - static_cast<double>(k_wrap) * period;
    if (w < 0.0) w += period;
    return w + x_min;
}

__global__ void spline2d_many_kernel(const spline2d_handle h,
                                     const double* __restrict__ x,
                                     double* __restrict__ y,
                                     int npts) {
    int ipt = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipt >= npts) return;

    double period1 = h.h1 * static_cast<double>(h.n1 - 1);
    double period2 = h.h2 * static_cast<double>(h.n2 - 1);

    double x1 = x[2 * ipt + 0];
    double x2 = x[2 * ipt + 1];
    double xj1 = h.periodic1 ? wrap_periodic(x1, h.x_min1, period1) : x1;
    double xj2 = h.periodic2 ? wrap_periodic(x2, h.x_min2, period2) : x2;

    double x_norm1 = (xj1 - h.x_min1) / h.h1;
    double x_norm2 = (xj2 - h.x_min2) / h.h2;
    int i1 = static_cast<int>(x_norm1);
    int i2 = static_cast<int>(x_norm2);
    if (i1 < 0) i1 = 0;
    if (i1 > h.n1 - 2) i1 = h.n1 - 2;
    if (i2 < 0) i2 = 0;
    if (i2 > h.n2 - 2) i2 = h.n2 - 2;

    double x_local1 = (x_norm1 - static_cast<double>(i1)) * h.h1;
    double x_local2 = (x_norm2 - static_cast<double>(i2)) * h.h2;

    std::size_t stride_k1 = static_cast<std::size_t>(h.nq);
    std::size_t stride_k2 = stride_k1 * static_cast<std::size_t>(h.order1 + 1);
    std::size_t stride_i1 = stride_k2 * static_cast<std::size_t>(h.order2 + 1);
    std::size_t stride_i2 = stride_i1 * static_cast<std::size_t>(h.n1);
    const double* coeff = h.coeff_dev + static_cast<std::size_t>(i1) * stride_i1 +
                          static_cast<std::size_t>(i2) * stride_i2;

    double* y_out = y + static_cast<std::size_t>(ipt) * static_cast<std::size_t>(h.nq);

    for (int iq = 0; iq < h.nq; ++iq) {
        double v = coeff[static_cast<std::size_t>(h.order2) * stride_k2 +
                         static_cast<std::size_t>(h.order1) * stride_k1 +
                         static_cast<std::size_t>(iq)];
        for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
            v = coeff[static_cast<std::size_t>(h.order2) * stride_k2 +
                      static_cast<std::size_t>(k1) * stride_k1 +
                      static_cast<std::size_t>(iq)] +
                x_local1 * v;
        }
        double yq = v;
        for (int k2 = h.order2 - 1; k2 >= 0; --k2) {
            v = coeff[static_cast<std::size_t>(k2) * stride_k2 +
                      static_cast<std::size_t>(h.order1) * stride_k1 +
                      static_cast<std::size_t>(iq)];
            for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
                v = coeff[static_cast<std::size_t>(k2) * stride_k2 +
                          static_cast<std::size_t>(k1) * stride_k1 +
                          static_cast<std::size_t>(iq)] +
                    x_local1 * v;
            }
            yq = v + x_local2 * yq;
        }
        y_out[iq] = yq;
    }
}

extern "C" void spline2d_many_cuda_c_init(int order1, int order2, int n1, int n2, int nq,
                                         int periodic1, int periodic2, double x_min1,
                                         double x_min2, double h1, double h2, const void* coeff,
                                         void** handle) {
    auto* h = new spline2d_handle();
    h->order1 = order1;
    h->order2 = order2;
    h->n1 = n1;
    h->n2 = n2;
    h->nq = nq;
    h->periodic1 = periodic1;
    h->periodic2 = periodic2;
    h->x_min1 = x_min1;
    h->x_min2 = x_min2;
    h->h1 = h1;
    h->h2 = h2;
    h->x_dev = nullptr;
    h->y_dev = nullptr;
    h->capacity_npts = 0;

    std::size_t n = static_cast<std::size_t>(nq) * static_cast<std::size_t>(order1 + 1) *
                    static_cast<std::size_t>(order2 + 1) * static_cast<std::size_t>(n1) *
                    static_cast<std::size_t>(n2);
    cudaMalloc(&h->coeff_dev, n * sizeof(double));
    cudaMemcpy(h->coeff_dev, coeff, n * sizeof(double), cudaMemcpyHostToDevice);

    *handle = h;
}

extern "C" void spline2d_many_cuda_c_free(void* handle) {
    auto* h = static_cast<spline2d_handle*>(handle);
    if (!h) return;
    if (h->y_dev) cudaFree(h->y_dev);
    if (h->x_dev) cudaFree(h->x_dev);
    cudaFree(h->coeff_dev);
    delete h;
}

extern "C" void spline2d_many_cuda_c_set_x(void* handle, const void* x, int npts) {
    auto* h = static_cast<spline2d_handle*>(handle);
    if (npts > h->capacity_npts) {
        if (h->y_dev) cudaFree(h->y_dev);
        if (h->x_dev) cudaFree(h->x_dev);
        cudaMalloc(&h->x_dev, static_cast<std::size_t>(2 * npts) * sizeof(double));
        cudaMalloc(&h->y_dev,
                   static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->nq) *
                       sizeof(double));
        h->capacity_npts = npts;
    }
    cudaMemcpy(h->x_dev, x, static_cast<std::size_t>(2 * npts) * sizeof(double),
               cudaMemcpyHostToDevice);
}

extern "C" void spline2d_many_cuda_c_eval_device(void* handle, int npts) {
    auto* h = static_cast<spline2d_handle*>(handle);
    int block = 256;
    int grid = (npts + block - 1) / block;
    spline2d_many_kernel<<<grid, block>>>(*h, h->x_dev, h->y_dev, npts);
}

extern "C" void spline2d_many_cuda_c_get_y(void* handle, void* y, int npts) {
    auto* h = static_cast<spline2d_handle*>(handle);
    cudaMemcpy(y, h->y_dev,
               static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->nq) *
                   sizeof(double),
               cudaMemcpyDeviceToHost);
}
