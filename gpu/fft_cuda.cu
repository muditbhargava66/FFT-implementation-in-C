/*
 * CUDA/cuFFT backend.
 *
 * This translation unit is compiled by nvcc as C++, so it deliberately does
 * not include fft_common.h or fft_gpu.h: their C11 `_Complex` typedef is not
 * accepted by the host C++ front end. The public dispatcher declares these
 * functions with complex_t* arguments; C linkage does not encode pointer
 * types, and complex_t has the same two-double layout as cufftDoubleComplex,
 * so the declarations remain ABI-compatible.
 *
 * Forward transforms are unscaled; inverse transforms are scaled by 1/N on
 * the device, matching the CPU kernels.
 */

#include <cuda_runtime.h>
#include <cufft.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct fft_gpu_memory {
  void *device_ptr;
  size_t elements;
};

struct fft_gpu_plan {
  cufftHandle handle;
  int n;       /* transform length for 1-D, rows*cols for 2-D */
  int batch;   /* 1 for 2-D */
  int rows;    /* 0 for a 1-D plan */
  int cols;
  int direction; /* -1 forward, +1 inverse */
  double inverse_scale;
};

static int g_cuda_initialized = 0;
static int g_device_id = 0;

int fft_gpu_dft_1d_batch_cuda(double *in, double *out, int n, int batch,
                              int direction);

static int elements_fit(size_t elements) {
  return elements > 0 && elements <= SIZE_MAX / sizeof(cufftDoubleComplex) &&
         elements <= (size_t)INT32_MAX;
}

int fft_gpu_init_cuda(void) {
  if (g_cuda_initialized)
    return 0;

  int device_count = 0;
  if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count <= 0)
    return -1;

  int best_device = 0;
  int best_multiprocessors = -1;
  for (int i = 0; i < device_count; i++) {
    cudaDeviceProp props;
    if (cudaGetDeviceProperties(&props, i) != cudaSuccess)
      continue;
    if (props.multiProcessorCount > best_multiprocessors) {
      best_multiprocessors = props.multiProcessorCount;
      best_device = i;
    }
  }

  if (cudaSetDevice(best_device) != cudaSuccess)
    return -1;

  g_device_id = best_device;
  g_cuda_initialized = 1;
  return 0;
}

void fft_gpu_cleanup_cuda(void) {
  if (!g_cuda_initialized)
    return;
  /* Do not call cudaDeviceReset(): it destroys the primary context for the
   * whole process, including work owned by other CUDA users in the host
   * application. Finish outstanding work and drop the cached state instead. */
  cudaDeviceSynchronize();
  g_cuda_initialized = 0;
}

int fft_gpu_available_cuda(void) {
  int device_count = 0;
  if (cudaGetDeviceCount(&device_count) != cudaSuccess)
    return 0;
  return device_count > 0;
}

struct fft_gpu_memory *fft_gpu_alloc_cuda(size_t elements) {
  if (!elements_fit(elements))
    return NULL;

  struct fft_gpu_memory *mem =
      (struct fft_gpu_memory *)malloc(sizeof(struct fft_gpu_memory));
  if (!mem)
    return NULL;

  mem->elements = elements;
  if (cudaMalloc(&mem->device_ptr, elements * sizeof(cufftDoubleComplex)) !=
      cudaSuccess) {
    free(mem);
    return NULL;
  }
  return mem;
}

void fft_gpu_free_cuda(struct fft_gpu_memory *mem) {
  if (!mem)
    return;
  if (mem->device_ptr)
    cudaFree(mem->device_ptr);
  free(mem);
}

void fft_gpu_copy_h2d_cuda(struct fft_gpu_memory *dst, const double *src,
                           size_t elements) {
  if (!dst || !src || elements == 0 || elements > dst->elements)
    return;
  cudaMemcpy(dst->device_ptr, src, elements * sizeof(cufftDoubleComplex),
             cudaMemcpyHostToDevice);
}

void fft_gpu_copy_d2h_cuda(double *dst, const struct fft_gpu_memory *src,
                           size_t elements) {
  if (!dst || !src || elements == 0 || elements > src->elements)
    return;
  cudaMemcpy(dst, src->device_ptr, elements * sizeof(cufftDoubleComplex),
             cudaMemcpyDeviceToHost);
}

__global__ void scale_kernel(cufftDoubleComplex *data, int n, double scale) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    data[idx].x *= scale;
    data[idx].y *= scale;
  }
}

struct fft_gpu_plan *fft_gpu_plan_1d_cuda(int n, int batch, int direction) {
  if (n <= 0 || batch <= 0 || !elements_fit((size_t)n * (size_t)batch))
    return NULL;

  struct fft_gpu_plan *plan =
      (struct fft_gpu_plan *)calloc(1, sizeof(struct fft_gpu_plan));
  if (!plan)
    return NULL;

  plan->n = n;
  plan->batch = batch;
  plan->rows = 0;
  plan->cols = 0;
  plan->direction = direction;
  plan->inverse_scale = 1.0 / (double)n;

  int dims[1] = {n};
  if (cufftPlanMany(&plan->handle, 1, dims, NULL, 1, n, NULL, 1, n,
                    CUFFT_Z2Z, batch) != CUFFT_SUCCESS) {
    free(plan);
    return NULL;
  }
  return plan;
}

struct fft_gpu_plan *fft_gpu_plan_2d_cuda(int rows, int cols, int direction) {
  if (rows <= 0 || cols <= 0 ||
      !elements_fit((size_t)rows * (size_t)cols))
    return NULL;

  struct fft_gpu_plan *plan =
      (struct fft_gpu_plan *)calloc(1, sizeof(struct fft_gpu_plan));
  if (!plan)
    return NULL;

  plan->n = rows * cols;
  plan->batch = 1;
  plan->rows = rows;
  plan->cols = cols;
  plan->direction = direction;
  plan->inverse_scale = 1.0 / ((double)rows * (double)cols);

  /* For the row-major buffer used here, data[r*cols + c], cuFFT's 2-D plan
   * takes the row count first: plan2d(rows, cols) matches a direct DFT over r
   * and c. Swapping the arguments transposes the result. Verified against the
   * rectangular 8x16 and 5x7 cases in test_gpu_cuda. */
  if (cufftPlan2d(&plan->handle, rows, cols, CUFFT_Z2Z) != CUFFT_SUCCESS) {
    free(plan);
    return NULL;
  }
  return plan;
}

int fft_gpu_execute_cuda(struct fft_gpu_plan *plan,
                         const struct fft_gpu_memory *in,
                         const struct fft_gpu_memory *out) {
  if (!plan || !in || !out || !in->device_ptr || !out->device_ptr)
    return -1;

  size_t total = plan->rows > 0 ? (size_t)plan->n
                                : (size_t)plan->n * (size_t)plan->batch;
  if (in->elements < total || out->elements < total || total > INT32_MAX)
    return -1;

  cufftDoubleComplex *in_ptr = (cufftDoubleComplex *)in->device_ptr;
  cufftDoubleComplex *out_ptr = (cufftDoubleComplex *)out->device_ptr;
  int direction = plan->direction < 0 ? CUFFT_FORWARD : CUFFT_INVERSE;

  if (cufftExecZ2Z(plan->handle, in_ptr, out_ptr, direction) != CUFFT_SUCCESS)
    return -1;

  if (direction == CUFFT_INVERSE) {
    int threads = 256;
    size_t blocks = (total + (size_t)threads - 1) / (size_t)threads;
    scale_kernel<<<(unsigned int)blocks, threads>>>(out_ptr, (int)total,
                                                    plan->inverse_scale);
    if (cudaGetLastError() != cudaSuccess)
      return -1;
  }

  return cudaDeviceSynchronize() == cudaSuccess ? 0 : -1;
}

void fft_gpu_destroy_plan_cuda(struct fft_gpu_plan *plan) {
  if (!plan)
    return;
  if (plan->handle)
    cufftDestroy(plan->handle);
  free(plan);
}

const char *fft_gpu_get_device_name_cuda(void) {
  static char name[256];
  cudaDeviceProp props;
  if (cudaGetDeviceProperties(&props, g_device_id) != cudaSuccess)
    return "Unknown CUDA Device";
  snprintf(name, sizeof(name), "%s", props.name);
  return name;
}

void fft_gpu_get_memory_info_cuda(size_t *total, size_t *available) {
  if (!total || !available)
    return;
  *total = 0;
  *available = 0;
  cudaMemGetInfo(available, total);
}

int fft_gpu_dft_1d_cuda(double *in, double *out, int n, int direction) {
  return fft_gpu_dft_1d_batch_cuda(in, out, n, 1, direction);
}

int fft_gpu_dft_1d_batch_cuda(double *in, double *out, int n, int batch,
                              int direction) {
  if (!in || !out || n <= 0 || batch <= 0)
    return -1;
  size_t total = (size_t)n * (size_t)batch;
  if (!elements_fit(total))
    return -1;
  if (fft_gpu_init_cuda() != 0)
    return -1;

  struct fft_gpu_memory *gpu_in = fft_gpu_alloc_cuda(total);
  struct fft_gpu_memory *gpu_out = fft_gpu_alloc_cuda(total);
  struct fft_gpu_plan *plan = fft_gpu_plan_1d_cuda(n, batch, direction);
  if (!gpu_in || !gpu_out || !plan) {
    fft_gpu_free_cuda(gpu_in);
    fft_gpu_free_cuda(gpu_out);
    fft_gpu_destroy_plan_cuda(plan);
    return -1;
  }

  int result = -1;
  fft_gpu_copy_h2d_cuda(gpu_in, in, total);
  if (fft_gpu_execute_cuda(plan, gpu_in, gpu_out) == 0) {
    fft_gpu_copy_d2h_cuda(out, gpu_out, total);
    result = 0;
  }

  fft_gpu_destroy_plan_cuda(plan);
  fft_gpu_free_cuda(gpu_in);
  fft_gpu_free_cuda(gpu_out);
  return result;
}

int fft_gpu_dft_2d_cuda(double *in, double *out, int rows, int cols,
                        int direction) {
  if (!in || !out || rows <= 0 || cols <= 0)
    return -1;
  size_t total = (size_t)rows * (size_t)cols;
  if (!elements_fit(total))
    return -1;
  if (fft_gpu_init_cuda() != 0)
    return -1;

  struct fft_gpu_memory *gpu_in = fft_gpu_alloc_cuda(total);
  struct fft_gpu_memory *gpu_out = fft_gpu_alloc_cuda(total);
  struct fft_gpu_plan *plan = fft_gpu_plan_2d_cuda(rows, cols, direction);
  if (!gpu_in || !gpu_out || !plan) {
    fft_gpu_free_cuda(gpu_in);
    fft_gpu_free_cuda(gpu_out);
    fft_gpu_destroy_plan_cuda(plan);
    return -1;
  }

  int result = -1;
  fft_gpu_copy_h2d_cuda(gpu_in, in, total);
  if (fft_gpu_execute_cuda(plan, gpu_in, gpu_out) == 0) {
    fft_gpu_copy_d2h_cuda(out, gpu_out, total);
    result = 0;
  }

  fft_gpu_destroy_plan_cuda(plan);
  fft_gpu_free_cuda(gpu_in);
  fft_gpu_free_cuda(gpu_out);
  return result;
}

#ifdef __cplusplus
}
#endif
