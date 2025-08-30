#ifdef __CUDACC__

#include "../../include/fft_gpu.h"
#include <cuda_runtime.h>
#include <cufft.h>
#include <stdio.h>

/**
 * @file fft_cuda.cu
 * @brief CUDA implementation of GPU-accelerated FFT
 * 
 * New in v2.0.0: Massive performance gains using NVIDIA GPUs
 */

// GPU memory structure
struct fft_gpu_memory {
    void* device_ptr;
    size_t size;
};

// GPU plan structure
struct fft_gpu_plan {
    cufftHandle cufft_plan;
    int n;
    int batch;
    cufftType type;
};

// Global CUDA state
static int g_cuda_initialized = 0;
static int g_device_id = 0;

// Error checking macro
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(error)); \
        return -1; \
    } \
} while(0)

#define CUFFT_CHECK(call) do { \
    cufftResult error = call; \
    if (error != CUFFT_SUCCESS) { \
        fprintf(stderr, "cuFFT error at %s:%d: %d\n", \
                __FILE__, __LINE__, error); \
        return -1; \
    } \
} while(0)

// Initialize CUDA
int fft_gpu_init_cuda(void) {
    if (g_cuda_initialized) return 0;
    
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    
    if (device_count == 0) {
        fprintf(stderr, "No CUDA devices found\n");
        return -1;
    }
    
    // Select best device
    int best_device = 0;
    int max_multiprocessors = 0;
    
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp props;
        CUDA_CHECK(cudaGetDeviceProperties(&props, i));
        
        if (props.multiProcessorCount > max_multiprocessors) {
            max_multiprocessors = props.multiProcessorCount;
            best_device = i;
        }
    }
    
    CUDA_CHECK(cudaSetDevice(best_device));
    g_device_id = best_device;
    g_cuda_initialized = 1;
    
    return 0;
}

// Cleanup CUDA
void fft_gpu_cleanup_cuda(void) {
    if (!g_cuda_initialized) return;
    
    cudaDeviceReset();
    g_cuda_initialized = 0;
}

// Check availability
int fft_gpu_available_cuda(void) {
    int device_count;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess) {
        return 0;
    }
    return device_count > 0;
}

// Memory allocation
fft_gpu_memory_t fft_gpu_alloc_cuda(size_t n) {
    fft_gpu_memory_t mem = malloc(sizeof(struct fft_gpu_memory));
    if (!mem) return NULL;
    
    mem->size = n * sizeof(complex_t);
    
    if (cudaMalloc(&mem->device_ptr, mem->size) != cudaSuccess) {
        free(mem);
        return NULL;
    }
    
    return mem;
}

// Memory free
void fft_gpu_free_cuda(fft_gpu_memory_t mem) {
    if (!mem) return;
    
    cudaFree(mem->device_ptr);
    free(mem);
}

// Memory copy host to device
void fft_gpu_copy_h2d_cuda(fft_gpu_memory_t dst, const complex_t* src, size_t n) {
    cudaMemcpy(dst->device_ptr, src, n * sizeof(complex_t), 
               cudaMemcpyHostToDevice);
}

// Memory copy device to host
void fft_gpu_copy_d2h_cuda(complex_t* dst, fft_gpu_memory_t src, size_t n) {
    cudaMemcpy(dst, src->device_ptr, n * sizeof(complex_t),
               cudaMemcpyDeviceToHost);
}

// Create FFT plan
fft_gpu_plan_t fft_gpu_plan_1d_cuda(int n, int batch, fft_direction dir) {
    fft_gpu_plan_t plan = malloc(sizeof(struct fft_gpu_plan));
    if (!plan) return NULL;
    
    plan->n = n;
    plan->batch = batch;
    plan->type = (dir == FFT_FORWARD) ? CUFFT_Z2Z : CUFFT_Z2Z;
    
    if (batch == 1) {
        if (cufftPlan1d(&plan->cufft_plan, n, plan->type, 1) != CUFFT_SUCCESS) {
            free(plan);
            return NULL;
        }
    } else {
        int dims[] = {n};
        if (cufftPlanMany(&plan->cufft_plan, 1, dims,
                         NULL, 1, n,
                         NULL, 1, n,
                         plan->type, batch) != CUFFT_SUCCESS) {
            free(plan);
            return NULL;
        }
    }
    
    return plan;
}

// Execute FFT
void fft_gpu_execute_cuda(fft_gpu_plan_t plan, fft_gpu_memory_t in, 
                         fft_gpu_memory_t out, fft_direction dir) {
    cufftDoubleComplex* in_ptr = (cufftDoubleComplex*)in->device_ptr;
    cufftDoubleComplex* out_ptr = (cufftDoubleComplex*)out->device_ptr;
    
    cufftExecZ2Z(plan->cufft_plan, in_ptr, out_ptr,
                 (dir == FFT_FORWARD) ? CUFFT_FORWARD : CUFFT_INVERSE);
    
    // Scale for inverse transform
    if (dir == FFT_INVERSE) {
        // Launch scaling kernel
        int threads = 256;
        int blocks = (plan->n * plan->batch + threads - 1) / threads;
        
        // scale_kernel<<<blocks, threads>>>(out_ptr, plan->n * plan->batch, 1.0 / plan->n);
        // TODO: Implement scaling kernel
    }
    
    cudaDeviceSynchronize();
}

// Destroy plan
void fft_gpu_destroy_plan_cuda(fft_gpu_plan_t plan) {
    if (!plan) return;
    
    cufftDestroy(plan->cufft_plan);
    free(plan);
}

// Get device name
const char* fft_gpu_get_device_name_cuda(void) {
    static char name[256];
    cudaDeviceProp props;
    
    if (cudaGetDeviceProperties(&props, g_device_id) == cudaSuccess) {
        snprintf(name, sizeof(name), "%s", props.name);
        return name;
    }
    
    return "Unknown CUDA Device";
}

// Get memory info
void fft_gpu_get_memory_info_cuda(size_t* total, size_t* available) {
    cudaMemGetInfo(available, total);
}

// Convenience function for 1D FFT
int fft_gpu_dft_1d_cuda(complex_t* in, complex_t* out, int n, fft_direction dir) {
    if (!g_cuda_initialized) {
        if (fft_gpu_init_cuda() != 0) return -1;
    }
    
    // Allocate GPU memory
    fft_gpu_memory_t gpu_in = fft_gpu_alloc_cuda(n);
    fft_gpu_memory_t gpu_out = fft_gpu_alloc_cuda(n);
    
    if (!gpu_in || !gpu_out) {
        fft_gpu_free_cuda(gpu_in);
        fft_gpu_free_cuda(gpu_out);
        return -1;
    }
    
    // Create plan
    fft_gpu_plan_t plan = fft_gpu_plan_1d_cuda(n, 1, dir);
    if (!plan) {
        fft_gpu_free_cuda(gpu_in);
        fft_gpu_free_cuda(gpu_out);
        return -1;
    }
    
    // Copy input to GPU
    fft_gpu_copy_h2d_cuda(gpu_in, in, n);
    
    // Execute FFT
    fft_gpu_execute_cuda(plan, gpu_in, gpu_out, dir);
    
    // Copy result back
    fft_gpu_copy_d2h_cuda(out, gpu_out, n);
    
    // Cleanup
    fft_gpu_destroy_plan_cuda(plan);
    fft_gpu_free_cuda(gpu_in);
    fft_gpu_free_cuda(gpu_out);
    
    return 0;
}

// Custom kernels for optimized operations

__global__ void scale_kernel(cufftDoubleComplex* data, int n, double scale) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        data[idx].x *= scale;
        data[idx].y *= scale;
    }
}

__global__ void butterfly_kernel(cufftDoubleComplex* data, int n, int stage) {
    // TODO: Implement custom butterfly kernel for small FFTs
}

#endif // __CUDACC__
