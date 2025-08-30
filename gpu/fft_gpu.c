#include "../include/fft_gpu.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @file fft_gpu.c
 * @brief GPU abstraction layer for FFT
 * 
 * Provides a unified interface for different GPU backends
 */

// Current backend
static fft_gpu_backend_t g_current_backend = FFT_GPU_NONE;

// Backend-specific functions
#ifdef __APPLE__
extern int fft_gpu_init_metal(void);
extern void fft_gpu_cleanup_metal(void);
extern int fft_gpu_available_metal(void);
extern fft_gpu_memory_t fft_gpu_alloc_metal(size_t n);
extern void fft_gpu_free_metal(fft_gpu_memory_t mem);
extern void fft_gpu_copy_h2d_metal(fft_gpu_memory_t dst, const complex_t* src, size_t n);
extern void fft_gpu_copy_d2h_metal(complex_t* dst, fft_gpu_memory_t src, size_t n);
extern fft_gpu_plan_t fft_gpu_plan_1d_metal(int n, int batch, fft_direction dir);
extern void fft_gpu_execute_metal(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out, fft_direction dir);
extern void fft_gpu_destroy_plan_metal(fft_gpu_plan_t plan);
extern const char* fft_gpu_get_device_name_metal(void);
extern void fft_gpu_get_memory_info_metal(size_t* total, size_t* available);
extern int fft_gpu_dft_1d_metal(complex_t* in, complex_t* out, int n, fft_direction dir);
#endif

#ifdef __CUDACC__
extern int fft_gpu_init_cuda(void);
extern void fft_gpu_cleanup_cuda(void);
extern int fft_gpu_available_cuda(void);
extern fft_gpu_memory_t fft_gpu_alloc_cuda(size_t n);
extern void fft_gpu_free_cuda(fft_gpu_memory_t mem);
extern void fft_gpu_copy_h2d_cuda(fft_gpu_memory_t dst, const complex_t* src, size_t n);
extern void fft_gpu_copy_d2h_cuda(complex_t* dst, fft_gpu_memory_t src, size_t n);
extern fft_gpu_plan_t fft_gpu_plan_1d_cuda(int n, int batch, fft_direction dir);
extern void fft_gpu_execute_cuda(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out, fft_direction dir);
extern void fft_gpu_destroy_plan_cuda(fft_gpu_plan_t plan);
extern const char* fft_gpu_get_device_name_cuda(void);
extern void fft_gpu_get_memory_info_cuda(size_t* total, size_t* available);
extern int fft_gpu_dft_1d_cuda(complex_t* in, complex_t* out, int n, fft_direction dir);
#endif

// Initialize GPU backend
int fft_gpu_init(fft_gpu_backend_t backend) {
    if (backend == FFT_GPU_AUTO) {
        // Try CUDA first, then Metal
        #ifdef __CUDACC__
        if (fft_gpu_available_cuda()) {
            backend = FFT_GPU_CUDA;
        }
        #endif
        
        #ifdef __APPLE__
        if (backend == FFT_GPU_AUTO && fft_gpu_available_metal()) {
            backend = FFT_GPU_METAL;
        }
        #endif
        
        if (backend == FFT_GPU_AUTO) {
            return -1;  // No GPU available
        }
    }
    
    int result = -1;
    
    switch (backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            result = fft_gpu_init_cuda();
            #else
            fprintf(stderr, "CUDA support not compiled in\n");
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            result = fft_gpu_init_metal();
            #else
            fprintf(stderr, "Metal support only available on macOS\n");
            #endif
            break;
            
        default:
            break;
    }
    
    if (result == 0) {
        g_current_backend = backend;
    }
    
    return result;
}

// Cleanup GPU resources
void fft_gpu_cleanup(void) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_cleanup_cuda();
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_cleanup_metal();
            #endif
            break;
            
        default:
            break;
    }
    
    g_current_backend = FFT_GPU_NONE;
}

// Check if GPU is available
int fft_gpu_available(void) {
    #ifdef __CUDACC__
    if (fft_gpu_available_cuda()) return 1;
    #endif
    
    #ifdef __APPLE__
    if (fft_gpu_available_metal()) return 1;
    #endif
    
    return 0;
}

// Get current GPU backend
fft_gpu_backend_t fft_gpu_get_backend(void) {
    return g_current_backend;
}

// GPU memory allocation
fft_gpu_memory_t fft_gpu_alloc(size_t size) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            return fft_gpu_alloc_cuda(size);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            return fft_gpu_alloc_metal(size);
            #endif
            break;
            
        default:
            break;
    }
    
    return NULL;
}

// GPU memory free
void fft_gpu_free(fft_gpu_memory_t mem) {
    if (!mem) return;
    
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_free_cuda(mem);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_free_metal(mem);
            #endif
            break;
            
        default:
            break;
    }
}

// Copy host to device
void fft_gpu_copy_h2d(fft_gpu_memory_t dst, const complex_t* src, size_t size) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_copy_h2d_cuda(dst, src, size);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_copy_h2d_metal(dst, src, size);
            #endif
            break;
            
        default:
            break;
    }
}

// Copy device to host
void fft_gpu_copy_d2h(complex_t* dst, fft_gpu_memory_t src, size_t size) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_copy_d2h_cuda(dst, src, size);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_copy_d2h_metal(dst, src, size);
            #endif
            break;
            
        default:
            break;
    }
}

// Create GPU FFT plan
fft_gpu_plan_t fft_gpu_plan_1d(int n, int batch, fft_direction direction) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            return fft_gpu_plan_1d_cuda(n, batch, direction);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            return fft_gpu_plan_1d_metal(n, batch, direction);
            #endif
            break;
            
        default:
            break;
    }
    
    return NULL;
}

// Execute GPU FFT plan
void fft_gpu_execute(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out) {
    if (!plan) return;
    
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_execute_cuda(plan, in, out, FFT_FORWARD);  // Direction from plan
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_execute_metal(plan, in, out, FFT_FORWARD);  // Direction from plan
            #endif
            break;
            
        default:
            break;
    }
}

// Destroy GPU plan
void fft_gpu_destroy_plan(fft_gpu_plan_t plan) {
    if (!plan) return;
    
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_destroy_plan_cuda(plan);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_destroy_plan_metal(plan);
            #endif
            break;
            
        default:
            break;
    }
}

// Get device name
const char* fft_gpu_get_device_name(void) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            return fft_gpu_get_device_name_cuda();
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            return fft_gpu_get_device_name_metal();
            #endif
            break;
            
        default:
            break;
    }
    
    return "No GPU";
}

// Get memory info
void fft_gpu_get_memory_info(size_t* total, size_t* available) {
    if (!total || !available) return;
    
    *total = 0;
    *available = 0;
    
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            fft_gpu_get_memory_info_cuda(total, available);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            fft_gpu_get_memory_info_metal(total, available);
            #endif
            break;
            
        default:
            break;
    }
}

// Convenience function for 1D FFT
int fft_gpu_dft_1d(complex_t* in, complex_t* out, int n, fft_direction direction) {
    switch (g_current_backend) {
        case FFT_GPU_CUDA:
            #ifdef __CUDACC__
            return fft_gpu_dft_1d_cuda(in, out, n, direction);
            #endif
            break;
            
        case FFT_GPU_METAL:
            #ifdef __APPLE__
            return fft_gpu_dft_1d_metal(in, out, n, direction);
            #endif
            break;
            
        default:
            break;
    }
    
    return -1;
}

// Set GPU device
int fft_gpu_set_device(int device) {
    // TODO: Implement multi-GPU support
    (void)device;
    return 0;
}

// Batched FFT
int fft_gpu_dft_1d_batch(complex_t* in, complex_t* out, int n, int batch,
                        fft_direction direction) {
    // Simple implementation using single FFTs
    for (int i = 0; i < batch; i++) {
        int result = fft_gpu_dft_1d(in + i * n, out + i * n, n, direction);
        if (result != 0) return result;
    }
    return 0;
}

// 2D FFT placeholder
fft_gpu_plan_t fft_gpu_plan_2d(int rows, int cols, fft_direction direction) {
    // TODO: Implement 2D GPU FFT
    (void)rows;
    (void)cols;
    (void)direction;
    return NULL;
}

int fft_gpu_dft_2d(complex_t* in, complex_t* out, int rows, int cols,
                  fft_direction direction) {
    // TODO: Implement 2D GPU FFT
    (void)in;
    (void)out;
    (void)rows;
    (void)cols;
    (void)direction;
    return -1;
}
