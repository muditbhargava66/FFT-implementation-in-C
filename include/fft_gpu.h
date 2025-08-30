#ifndef FFT_GPU_H
#define FFT_GPU_H

#include "fft_common.h"

/**
 * @file fft_gpu.h
 * @brief GPU acceleration for FFT (CUDA and Metal Performance Shaders)
 * 
 * New in v2.0.0: GPU support for massive performance gains
 */

// GPU backend selection
typedef enum {
    FFT_GPU_NONE = 0,
    FFT_GPU_CUDA = 1,
    FFT_GPU_METAL = 2,
    FFT_GPU_OPENCL = 3,
    FFT_GPU_AUTO = -1  // Automatically select best available
} fft_gpu_backend_t;

// GPU memory handle (opaque)
typedef struct fft_gpu_memory* fft_gpu_memory_t;

// GPU plan (opaque)
typedef struct fft_gpu_plan* fft_gpu_plan_t;

// GPU initialization and management

/**
 * @brief Initialize GPU backend
 * @param backend GPU backend to use (or FFT_GPU_AUTO)
 * @return 0 on success, -1 on failure
 */
int fft_gpu_init(fft_gpu_backend_t backend);

/**
 * @brief Cleanup GPU resources
 */
void fft_gpu_cleanup(void);

/**
 * @brief Check if GPU is available
 * @return 1 if available, 0 otherwise
 */
int fft_gpu_available(void);

/**
 * @brief Get current GPU backend
 * @return Current backend type
 */
fft_gpu_backend_t fft_gpu_get_backend(void);

// GPU memory management

/**
 * @brief Allocate memory on GPU
 * @param size Number of complex elements
 * @return GPU memory handle
 */
fft_gpu_memory_t fft_gpu_alloc(size_t size);

/**
 * @brief Free GPU memory
 * @param mem GPU memory handle
 */
void fft_gpu_free(fft_gpu_memory_t mem);

/**
 * @brief Copy data from host to GPU
 * @param dst GPU memory destination
 * @param src Host memory source
 * @param size Number of complex elements
 */
void fft_gpu_copy_h2d(fft_gpu_memory_t dst, const complex_t* src, size_t size);

/**
 * @brief Copy data from GPU to host
 * @param dst Host memory destination
 * @param src GPU memory source
 * @param size Number of complex elements
 */
void fft_gpu_copy_d2h(complex_t* dst, fft_gpu_memory_t src, size_t size);

// GPU FFT planning and execution

/**
 * @brief Create GPU FFT plan
 * @param n Transform size
 * @param batch Number of transforms
 * @param direction FFT direction
 * @return GPU plan handle
 */
fft_gpu_plan_t fft_gpu_plan_1d(int n, int batch, fft_direction direction);

/**
 * @brief Execute GPU FFT plan
 * @param plan GPU plan
 * @param in Input GPU memory
 * @param out Output GPU memory (can be same as in)
 */
void fft_gpu_execute(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out);

/**
 * @brief Destroy GPU plan
 * @param plan GPU plan to destroy
 */
void fft_gpu_destroy_plan(fft_gpu_plan_t plan);

// High-level convenience functions

/**
 * @brief Perform FFT on GPU (handles all memory transfers)
 * @param in Host input array
 * @param out Host output array
 * @param n Transform size
 * @param direction FFT direction
 * @return 0 on success, -1 on failure
 */
int fft_gpu_dft_1d(complex_t* in, complex_t* out, int n, fft_direction direction);

/**
 * @brief Perform batched FFT on GPU
 * @param in Host input array (size n * batch)
 * @param out Host output array
 * @param n Transform size
 * @param batch Number of transforms
 * @param direction FFT direction
 * @return 0 on success, -1 on failure
 */
int fft_gpu_dft_1d_batch(complex_t* in, complex_t* out, int n, int batch, 
                         fft_direction direction);

// 2D FFT on GPU

/**
 * @brief Create 2D GPU FFT plan
 * @param rows Number of rows
 * @param cols Number of columns
 * @param direction FFT direction
 * @return GPU plan handle
 */
fft_gpu_plan_t fft_gpu_plan_2d(int rows, int cols, fft_direction direction);

/**
 * @brief Perform 2D FFT on GPU
 * @param in Host input array (row-major)
 * @param out Host output array
 * @param rows Number of rows
 * @param cols Number of columns
 * @param direction FFT direction
 * @return 0 on success, -1 on failure
 */
int fft_gpu_dft_2d(complex_t* in, complex_t* out, int rows, int cols,
                   fft_direction direction);

// GPU information

/**
 * @brief Get GPU device name
 * @return Device name string
 */
const char* fft_gpu_get_device_name(void);

/**
 * @brief Get GPU memory info
 * @param total Total memory in bytes
 * @param available Available memory in bytes
 */
void fft_gpu_get_memory_info(size_t* total, size_t* available);

/**
 * @brief Set GPU device (for multi-GPU systems)
 * @param device Device index
 * @return 0 on success, -1 on failure
 */
int fft_gpu_set_device(int device);

// Platform-specific features

#ifdef __APPLE__
// Metal Performance Shaders specific
typedef struct {
    int prefer_shared_memory;
    int allow_reduced_precision;
} fft_mps_options_t;

void fft_gpu_set_mps_options(const fft_mps_options_t* options);
#endif

#ifdef __CUDACC__
// CUDA specific
typedef struct {
    int block_size;
    int shared_memory_size;
    cudaStream_t stream;
} fft_cuda_options_t;

void fft_gpu_set_cuda_options(const fft_cuda_options_t* options);
#endif

#endif // FFT_GPU_H
