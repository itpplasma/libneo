#include <cuda_runtime.h>
#include <cstddef>

struct spline3d_handle {
    int order1, order2, order3;
    int n1, n2, n3;
    int nq;
    int p1, p2, p3;
    double x_min1, x_min2, x_min3;
    double h1, h2, h3;
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

__global__ void spline3d_many_kernel(const spline3d_handle h, const double* __restrict__ x,
                                     double* __restrict__ y, int npts) {
    int ipt = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipt >= npts) return;

    double period1 = h.h1 * static_cast<double>(h.n1 - 1);
    double period2 = h.h2 * static_cast<double>(h.n2 - 1);
    double period3 = h.h3 * static_cast<double>(h.n3 - 1);

    double x1 = x[3 * ipt + 0];
    double x2 = x[3 * ipt + 1];
    double x3 = x[3 * ipt + 2];
    double xj1 = h.p1 ? wrap_periodic(x1, h.x_min1, period1) : x1;
    double xj2 = h.p2 ? wrap_periodic(x2, h.x_min2, period2) : x2;
    double xj3 = h.p3 ? wrap_periodic(x3, h.x_min3, period3) : x3;

    double x_norm1 = (xj1 - h.x_min1) / h.h1;
    double x_norm2 = (xj2 - h.x_min2) / h.h2;
    double x_norm3 = (xj3 - h.x_min3) / h.h3;
    int i1 = static_cast<int>(x_norm1);
    int i2 = static_cast<int>(x_norm2);
    int i3 = static_cast<int>(x_norm3);
    if (i1 < 0) i1 = 0;
    if (i1 > h.n1 - 2) i1 = h.n1 - 2;
    if (i2 < 0) i2 = 0;
    if (i2 > h.n2 - 2) i2 = h.n2 - 2;
    if (i3 < 0) i3 = 0;
    if (i3 > h.n3 - 2) i3 = h.n3 - 2;

    double x_local1 = (x_norm1 - static_cast<double>(i1)) * h.h1;
    double x_local2 = (x_norm2 - static_cast<double>(i2)) * h.h2;
    double x_local3 = (x_norm3 - static_cast<double>(i3)) * h.h3;

    std::size_t stride_k1 = static_cast<std::size_t>(h.nq);
    std::size_t stride_k2 = stride_k1 * static_cast<std::size_t>(h.order1 + 1);
    std::size_t stride_k3 = stride_k2 * static_cast<std::size_t>(h.order2 + 1);
    std::size_t stride_i1 = stride_k3 * static_cast<std::size_t>(h.order3 + 1);
    std::size_t stride_i2 = stride_i1 * static_cast<std::size_t>(h.n1);
    std::size_t stride_i3 = stride_i2 * static_cast<std::size_t>(h.n2);
    const double* coeff = h.coeff_dev + static_cast<std::size_t>(i1) * stride_i1 +
                          static_cast<std::size_t>(i2) * stride_i2 +
                          static_cast<std::size_t>(i3) * stride_i3;

    double* y_out = y + static_cast<std::size_t>(ipt) * static_cast<std::size_t>(h.nq);

    for (int iq = 0; iq < h.nq; ++iq) {
        double w2 = coeff[static_cast<std::size_t>(h.order3) * stride_k3 +
                          static_cast<std::size_t>(h.order2) * stride_k2 +
                          static_cast<std::size_t>(h.order1) * stride_k1 +
                          static_cast<std::size_t>(iq)];
        for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
            w2 = coeff[static_cast<std::size_t>(h.order3) * stride_k3 +
                       static_cast<std::size_t>(h.order2) * stride_k2 +
                       static_cast<std::size_t>(k1) * stride_k1 +
                       static_cast<std::size_t>(iq)] +
                 x_local1 * w2;
        }
        for (int k2 = h.order2 - 1; k2 >= 0; --k2) {
            double v = coeff[static_cast<std::size_t>(h.order3) * stride_k3 +
                             static_cast<std::size_t>(k2) * stride_k2 +
                             static_cast<std::size_t>(h.order1) * stride_k1 +
                             static_cast<std::size_t>(iq)];
            for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
                v = coeff[static_cast<std::size_t>(h.order3) * stride_k3 +
                          static_cast<std::size_t>(k2) * stride_k2 +
                          static_cast<std::size_t>(k1) * stride_k1 +
                          static_cast<std::size_t>(iq)] +
                    x_local1 * v;
            }
            w2 = v + x_local2 * w2;
        }
        double yq = w2;

        for (int k3 = h.order3 - 1; k3 >= 0; --k3) {
            w2 = coeff[static_cast<std::size_t>(k3) * stride_k3 +
                       static_cast<std::size_t>(h.order2) * stride_k2 +
                       static_cast<std::size_t>(h.order1) * stride_k1 +
                       static_cast<std::size_t>(iq)];
            for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
                w2 = coeff[static_cast<std::size_t>(k3) * stride_k3 +
                           static_cast<std::size_t>(h.order2) * stride_k2 +
                           static_cast<std::size_t>(k1) * stride_k1 +
                           static_cast<std::size_t>(iq)] +
                     x_local1 * w2;
            }
            for (int k2 = h.order2 - 1; k2 >= 0; --k2) {
                double v = coeff[static_cast<std::size_t>(k3) * stride_k3 +
                                 static_cast<std::size_t>(k2) * stride_k2 +
                                 static_cast<std::size_t>(h.order1) * stride_k1 +
                                 static_cast<std::size_t>(iq)];
                for (int k1 = h.order1 - 1; k1 >= 0; --k1) {
                    v = coeff[static_cast<std::size_t>(k3) * stride_k3 +
                              static_cast<std::size_t>(k2) * stride_k2 +
                              static_cast<std::size_t>(k1) * stride_k1 +
                              static_cast<std::size_t>(iq)] +
                        x_local1 * v;
                }
                w2 = v + x_local2 * w2;
            }
            yq = w2 + x_local3 * yq;
        }
        y_out[iq] = yq;
    }
}

extern "C" void spline3d_many_cuda_c_init(int order1, int order2, int order3, int n1, int n2,
                                         int n3, int nq, int p1, int p2, int p3, double x_min1,
                                         double x_min2, double x_min3, double h1, double h2,
                                         double h3, const void* coeff, void** handle) {
    auto* h = new spline3d_handle();
    h->order1 = order1;
    h->order2 = order2;
    h->order3 = order3;
    h->n1 = n1;
    h->n2 = n2;
    h->n3 = n3;
    h->nq = nq;
    h->p1 = p1;
    h->p2 = p2;
    h->p3 = p3;
    h->x_min1 = x_min1;
    h->x_min2 = x_min2;
    h->x_min3 = x_min3;
    h->h1 = h1;
    h->h2 = h2;
    h->h3 = h3;
    h->x_dev = nullptr;
    h->y_dev = nullptr;
    h->capacity_npts = 0;

    std::size_t n = static_cast<std::size_t>(nq) * static_cast<std::size_t>(order1 + 1) *
                    static_cast<std::size_t>(order2 + 1) * static_cast<std::size_t>(order3 + 1) *
                    static_cast<std::size_t>(n1) * static_cast<std::size_t>(n2) *
                    static_cast<std::size_t>(n3);
    cudaMalloc(&h->coeff_dev, n * sizeof(double));
    cudaMemcpy(h->coeff_dev, coeff, n * sizeof(double), cudaMemcpyHostToDevice);
    *handle = h;
}

extern "C" void spline3d_many_cuda_c_free(void* handle) {
    auto* h = static_cast<spline3d_handle*>(handle);
    if (!h) return;
    if (h->y_dev) cudaFree(h->y_dev);
    if (h->x_dev) cudaFree(h->x_dev);
    cudaFree(h->coeff_dev);
    delete h;
}

extern "C" void spline3d_many_cuda_c_set_x(void* handle, const void* x, int npts) {
    auto* h = static_cast<spline3d_handle*>(handle);
    if (npts > h->capacity_npts) {
        if (h->y_dev) cudaFree(h->y_dev);
        if (h->x_dev) cudaFree(h->x_dev);
        cudaMalloc(&h->x_dev, static_cast<std::size_t>(3 * npts) * sizeof(double));
        cudaMalloc(&h->y_dev,
                   static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->nq) *
                       sizeof(double));
        h->capacity_npts = npts;
    }
    cudaMemcpy(h->x_dev, x, static_cast<std::size_t>(3 * npts) * sizeof(double),
               cudaMemcpyHostToDevice);
}

extern "C" void spline3d_many_cuda_c_eval_device(void* handle, int npts) {
    auto* h = static_cast<spline3d_handle*>(handle);
    int block = 256;
    int grid = (npts + block - 1) / block;
    spline3d_many_kernel<<<grid, block>>>(*h, h->x_dev, h->y_dev, npts);
}

extern "C" void spline3d_many_cuda_c_get_y(void* handle, void* y, int npts) {
    auto* h = static_cast<spline3d_handle*>(handle);
    cudaMemcpy(y, h->y_dev,
               static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->nq) *
                   sizeof(double),
               cudaMemcpyDeviceToHost);
}
