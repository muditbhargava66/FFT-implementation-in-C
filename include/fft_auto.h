#ifndef FFT_AUTO_H
#define FFT_AUTO_H

#include "fft_common.h"

/**
 * @file fft_auto.h
 * @brief Automatic FFT algorithm selection and planning
 * 
 * New in v2.0.0: Simplified API with automatic optimization
 */

// FFT Planner - automatically selects best algorithm
typedef struct fft_plan* fft_plan_t;

// Execution hints for optimization
typedef enum {
    FFT_ESTIMATE    = 0,     // Quick planning
    FFT_MEASURE     = 1,     // Benchmark algorithms  
    FFT_PATIENT     = 2,     // Extensive search
    FFT_EXHAUSTIVE  = 3,     // Try everything
    FFT_WISDOM_ONLY = 4,     // Use saved wisdom
    FFT_REAL_INPUT  = 1<<5,  // Real-valued input
    FFT_REAL_OUTPUT = 1<<6,  // Real-valued output
    FFT_UNALIGNED   = 1<<7,  // Data may be unaligned
    FFT_CONSERVE_MEMORY = 1<<8,  // Minimize memory usage
    FFT_PREFER_GPU  = 1<<9,  // Use GPU if available
    FFT_THREADED    = 1<<10  // Use multiple threads
} fft_flags_t;

// Main API Functions (v2.0.0)

/**
 * @brief Create an FFT plan with automatic algorithm selection
 * 
 * @param n Size of transform
 * @param in Input array (can be same as out for in-place)
 * @param out Output array
 * @param sign Direction (-1 for forward, +1 for inverse)
 * @param flags Optimization hints
 * @return Optimized FFT plan
 */
fft_plan_t fft_plan_dft_1d(int n, complex_t* in, complex_t* out, 
                           int sign, unsigned flags);

/**
 * @brief Execute an FFT plan
 * 
 * @param plan The FFT plan to execute
 */
void fft_execute(fft_plan_t plan);

/**
 * @brief Execute plan with different arrays
 * 
 * @param plan The FFT plan
 * @param in New input array
 * @param out New output array
 */
void fft_execute_dft(fft_plan_t plan, complex_t* in, complex_t* out);

/**
 * @brief Destroy an FFT plan
 * 
 * @param plan The plan to destroy
 */
void fft_destroy_plan(fft_plan_t plan);

// Simplified one-shot interface

/**
 * @brief Perform FFT with automatic optimization (v2.0.0)
 * 
 * Automatically selects the best algorithm based on:
 * - Transform size
 * - Hardware capabilities
 * - Available memory
 * 
 * @param in Input array
 * @param out Output array (can be same as in)
 * @param n Transform size
 * @param sign Direction (-1 forward, +1 inverse)
 * @return 0 on success, -1 on error
 */
int fft_auto(complex_t* in, complex_t* out, int n, int sign);

// Real-valued FFT interface

/**
 * @brief Plan for real-to-complex FFT
 * 
 * @param n Size of real input
 * @param in Real input array
 * @param out Complex output array (size n/2+1)
 * @param flags Optimization hints
 */
fft_plan_t fft_plan_r2c_1d(int n, double* in, complex_t* out, unsigned flags);

/**
 * @brief Plan for complex-to-real FFT
 * 
 * @param n Size of real output
 * @param in Complex input array (size n/2+1)
 * @param out Real output array
 * @param flags Optimization hints
 */
fft_plan_t fft_plan_c2r_1d(int n, complex_t* in, double* out, unsigned flags);

// 2D FFT interface

/**
 * @brief Plan for 2D complex FFT
 * 
 * @param rows Number of rows
 * @param cols Number of columns
 * @param in Input array (row-major order)
 * @param out Output array
 * @param sign Direction
 * @param flags Optimization hints
 */
fft_plan_t fft_plan_dft_2d(int rows, int cols, complex_t* in, 
                           complex_t* out, int sign, unsigned flags);

// Wisdom (saved optimization data)

/**
 * @brief Export wisdom to string
 * @return Wisdom string (must be freed by caller)
 */
char* fft_export_wisdom_to_string(void);

/**
 * @brief Import wisdom from string
 * @param wisdom Wisdom string
 * @return 1 on success, 0 on failure
 */
int fft_import_wisdom_from_string(const char* wisdom);

// Hardware detection

/**
 * @brief Query available hardware acceleration
 * @return Bitmask of available accelerators
 */
typedef enum {
    FFT_HW_CPU_SSE    = 1<<0,
    FFT_HW_CPU_AVX    = 1<<1,
    FFT_HW_CPU_AVX2   = 1<<2,
    FFT_HW_CPU_AVX512 = 1<<3,
    FFT_HW_CPU_NEON   = 1<<4,
    FFT_HW_GPU_CUDA   = 1<<5,
    FFT_HW_GPU_MPS    = 1<<6,
    FFT_HW_GPU_OPENCL = 1<<7
} fft_hardware_t;

unsigned fft_get_hardware_capabilities(void);

// Thread control

/**
 * @brief Set number of threads for parallel execution
 * @param nthreads Number of threads (0 for automatic)
 */
void fft_plan_with_nthreads(int nthreads);

// Memory allocation with alignment

/**
 * @brief Allocate aligned memory for optimal performance
 * @param n Number of complex elements
 * @return Aligned complex array
 */
complex_t* fft_alloc_complex(size_t n);

/**
 * @brief Allocate aligned memory for real data
 * @param n Number of real elements
 * @return Aligned real array
 */
double* fft_alloc_real(size_t n);

/**
 * @brief Free memory allocated by fft_alloc
 * @param p Pointer to free
 */
void fft_free(void* p);

// Version information

/**
 * @brief Get library version string
 * @return Version string
 */
const char* fft_version(void);

#endif // FFT_AUTO_H
