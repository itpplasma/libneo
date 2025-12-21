#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>

struct spline1d_handle {
    int order;
    int num_points;
    int num_quantities;
    int periodic;
    double x_min;
    double h_step;
    double* coeff_dev;
    double* x_dev;
    double* y_dev;
    int capacity_npts;
};

static __device__ __forceinline__ double wrap_periodic(double x, double x_min, double period) {
    double t = x - x_min;
    double w = t - floor(t / period) * period;
    return w + x_min;
}

__global__ void spline1d_many_kernel(const spline1d_handle h,
                                     const double* __restrict__ x,
                                     double* __restrict__ y,
                                     int npts)
{
    int ipt = blockIdx.x * blockDim.x + threadIdx.x;
    if (ipt >= npts) return;

    double period = h.h_step * static_cast<double>(h.num_points - 1);
    double xj = h.periodic ? wrap_periodic(x[ipt], h.x_min, period) : x[ipt];

    double x_norm = (xj - h.x_min) / h.h_step;
    int idx = static_cast<int>(x_norm);
    if (idx < 0) idx = 0;
    if (idx > h.num_points - 2) idx = h.num_points - 2;
    double x_local = (x_norm - static_cast<double>(idx)) * h.h_step;

    std::size_t stride_k = static_cast<std::size_t>(h.num_quantities);
    std::size_t stride_i = stride_k * static_cast<std::size_t>(h.order + 1);
    const double* coeff = h.coeff_dev + static_cast<std::size_t>(idx) * stride_i;

    double* y_out = y + static_cast<std::size_t>(ipt) * static_cast<std::size_t>(h.num_quantities);

    for (int iq = 0; iq < h.num_quantities; ++iq) {
        y_out[iq] = coeff[static_cast<std::size_t>(iq) + static_cast<std::size_t>(h.order) * stride_k];
    }
    for (int k = h.order - 1; k >= 0; --k) {
        std::size_t off = static_cast<std::size_t>(k) * stride_k;
        for (int iq = 0; iq < h.num_quantities; ++iq) {
            y_out[iq] = coeff[static_cast<std::size_t>(iq) + off] + x_local * y_out[iq];
        }
    }
}

extern "C" void spline1d_many_cuda_c_init(int order, int num_points, int num_quantities,
                                         int periodic, double x_min, double h_step,
                                         const void* coeff, void** handle)
{
    auto* h = new spline1d_handle();
    h->order = order;
    h->num_points = num_points;
    h->num_quantities = num_quantities;
    h->periodic = periodic;
    h->x_min = x_min;
    h->h_step = h_step;
    h->x_dev = nullptr;
    h->y_dev = nullptr;
    h->capacity_npts = 0;

    std::size_t n = static_cast<std::size_t>(num_quantities) *
                    static_cast<std::size_t>(order + 1) *
                    static_cast<std::size_t>(num_points);
    cudaMalloc(&h->coeff_dev, n * sizeof(double));
    cudaMemcpy(h->coeff_dev, coeff, n * sizeof(double), cudaMemcpyHostToDevice);

    *handle = h;
}

extern "C" void spline1d_many_cuda_c_free(void* handle)
{
    auto* h = static_cast<spline1d_handle*>(handle);
    if (!h) return;
    if (h->y_dev) cudaFree(h->y_dev);
    if (h->x_dev) cudaFree(h->x_dev);
    cudaFree(h->coeff_dev);
    delete h;
}

extern "C" void spline1d_many_cuda_c_eval(void* handle, const void* x, void* y, int npts)
{
    auto* h = static_cast<spline1d_handle*>(handle);
    if (npts > h->capacity_npts) {
        if (h->y_dev) cudaFree(h->y_dev);
        if (h->x_dev) cudaFree(h->x_dev);
        cudaMalloc(&h->x_dev, static_cast<std::size_t>(npts) * sizeof(double));
        cudaMalloc(&h->y_dev, static_cast<std::size_t>(npts) *
                              static_cast<std::size_t>(h->num_quantities) * sizeof(double));
        h->capacity_npts = npts;
    }

    cudaMemcpy(h->x_dev, x, static_cast<std::size_t>(npts) * sizeof(double),
               cudaMemcpyHostToDevice);

    int block = 256;
    int grid = (npts + block - 1) / block;
    spline1d_many_kernel<<<grid, block>>>(*h, h->x_dev, h->y_dev, npts);

    cudaMemcpy(y, h->y_dev,
               static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->num_quantities) *
                   sizeof(double),
               cudaMemcpyDeviceToHost);
}

extern "C" void spline1d_many_cuda_c_set_x(void* handle, const void* x, int npts)
{
    auto* h = static_cast<spline1d_handle*>(handle);
    if (npts > h->capacity_npts) {
        if (h->y_dev) cudaFree(h->y_dev);
        if (h->x_dev) cudaFree(h->x_dev);
        cudaMalloc(&h->x_dev, static_cast<std::size_t>(npts) * sizeof(double));
        cudaMalloc(&h->y_dev, static_cast<std::size_t>(npts) *
                              static_cast<std::size_t>(h->num_quantities) * sizeof(double));
        h->capacity_npts = npts;
    }
    cudaMemcpy(h->x_dev, x, static_cast<std::size_t>(npts) * sizeof(double),
               cudaMemcpyHostToDevice);
}

extern "C" void spline1d_many_cuda_c_eval_device(void* handle, int npts)
{
    auto* h = static_cast<spline1d_handle*>(handle);
    int block = 256;
    int grid = (npts + block - 1) / block;
    spline1d_many_kernel<<<grid, block>>>(*h, h->x_dev, h->y_dev, npts);
}

extern "C" void spline1d_many_cuda_c_get_y(void* handle, void* y, int npts)
{
    auto* h = static_cast<spline1d_handle*>(handle);
    cudaMemcpy(y, h->y_dev,
               static_cast<std::size_t>(npts) * static_cast<std::size_t>(h->num_quantities) *
                   sizeof(double),
               cudaMemcpyDeviceToHost);
}

extern "C" void spline1d_many_cuda_sync()
{
    cudaDeviceSynchronize();
}
