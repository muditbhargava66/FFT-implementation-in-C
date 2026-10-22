#include "fft_gpu.h"
#include <stddef.h>

/**
 * @file fft_gpu_stubs.c
 * @brief Stub implementations for GPU functions when GPU support is disabled
 */

#ifndef FFT_GPU_ENABLED

int fft_gpu_init(fft_gpu_backend_t backend) {
    (void)backend;
    return -1; // Fail
}

void fft_gpu_cleanup(void) {
    // No-op
}

int fft_gpu_available(void) {
    return 0; // Not available
}

fft_gpu_backend_t fft_gpu_get_backend(void) {
    return FFT_GPU_NONE;
}

fft_gpu_memory_t fft_gpu_alloc(size_t size) {
    (void)size;
    return NULL;
}

void fft_gpu_free(fft_gpu_memory_t mem) {
    (void)mem;
}

void fft_gpu_copy_h2d(fft_gpu_memory_t dst, const complex_t* src, size_t size) {
    (void)dst; (void)src; (void)size;
}

void fft_gpu_copy_d2h(complex_t* dst, fft_gpu_memory_t src, size_t size) {
    (void)dst; (void)src; (void)size;
}

fft_gpu_plan_t fft_gpu_plan_1d(int n, int batch, fft_direction direction) {
    (void)n; (void)batch; (void)direction;
    return NULL;
}

void fft_gpu_execute(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out) {
    (void)plan; (void)in; (void)out;
}

void fft_gpu_destroy_plan(fft_gpu_plan_t plan) {
    (void)plan;
}

int fft_gpu_dft_1d(complex_t* in, complex_t* out, int n, fft_direction direction) {
    (void)in; (void)out; (void)n; (void)direction;
    return -1;
}

int fft_gpu_dft_1d_batch(complex_t* in, complex_t* out, int n, int batch, fft_direction direction) {
    (void)in; (void)out; (void)n; (void)batch; (void)direction;
    return -1;
}

fft_gpu_plan_t fft_gpu_plan_2d(int rows, int cols, fft_direction direction) {
    (void)rows; (void)cols; (void)direction;
    return NULL;
}

int fft_gpu_dft_2d(complex_t* in, complex_t* out, int rows, int cols, fft_direction direction) {
    (void)in; (void)out; (void)rows; (void)cols; (void)direction;
    return -1;
}

const char* fft_gpu_get_device_name(void) {
    return "None (GPU Disabled)";
}

void fft_gpu_get_memory_info(size_t* total, size_t* available) {
    if (total) *total = 0;
    if (available) *available = 0;
}

int fft_gpu_set_device(int device) {
    (void)device;
    return -1;
}

#endif // !FFT_GPU_ENABLED
