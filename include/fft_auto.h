#ifndef FFT_AUTO_H
#define FFT_AUTO_H

#include "fft_common.h"

/**
 * @file fft_auto.h
 * @brief Automatic FFT algorithm selection and planning
 * 
 * v3.0.0: Simplified API with automatic optimization, r2c/c2r FFTs, 2D FFT, wisdom cache
 */

// FFT Planner - automatically selects best algorithm
typedef struct fft_plan* fft_plan_t;

// Execution hints for optimization.
//
// These are OR-combined, so each one needs its own bit. FFT_MEASURE and
// FFT_WISDOM_ONLY used to be the consecutive values 1 and 4 alongside two
// others, which meant a flag whose value was 0b11 also read as FFT_MEASURE
// (0b01) in every `flags & FFT_MEASURE` test. They are numbered as bits now.
//
// The gaps below are deliberate. FFT_PATIENT, FFT_EXHAUSTIVE, FFT_UNALIGNED and
// FFT_THREADED were removed because nothing implemented them, but the retained
// flags keep their existing bit values so exported wisdom and compiled callers
// do not silently change meaning. Do not renumber to close the gaps.
typedef enum {
    FFT_ESTIMATE    = 0,     // Quick planning: pick by the size heuristic
    FFT_MEASURE     = 1<<0,  // Time the candidate algorithms and keep the fastest
    FFT_WISDOM_ONLY = 1<<3,  // Use a cached choice, or fail; never measure
    FFT_REAL_INPUT  = 1<<5,  // Real-valued input
    FFT_REAL_OUTPUT = 1<<6,  // Real-valued output
    FFT_CONSERVE_MEMORY = 1<<8,  // Minimize memory usage; overrides FFT_MEASURE
    FFT_PREFER_GPU  = 1<<9   // Use GPU if available
} fft_flags_t;

// Main API Functions (v3.0.0)

/**
 * @brief Create an FFT plan with automatic algorithm selection
 * 
 * @param n Size of transform
 * @param in Input array (can be same as out for in-place)
 * @param out Output array
 * @param sign Direction (-1 for forward, +1 for inverse)
 * @param flags Optimization hints
 * @return Optimized FFT plan, or NULL if n <= 0, a buffer is NULL, allocation
 *         fails, or FFT_WISDOM_ONLY was set and no cached choice exists for
 *         (n, flags). Callers must check for NULL before executing.
 */
fft_plan_t fft_plan_dft_1d(int n, complex_t* in, complex_t* out, 
                           int sign, unsigned flags);

/**
 * @brief Execute an FFT plan
 * 
 * @param plan The FFT plan to execute
 * @return FFT_OK, or the error from the selected algorithm. A null plan gives
 *         FFT_ERR_NULL and leaves every buffer untouched.
 */
fft_error_t fft_execute(fft_plan_t plan);

/**
 * @brief Execute plan with different arrays
 * 
 * @param plan The FFT plan
 * @param in New input array
 * @param out New output array
 * @return FFT_OK, or FFT_ERR_NULL if plan/in/out is null
 */
fft_error_t fft_execute_dft(fft_plan_t plan, complex_t* in, complex_t* out);

/**
 * @brief Destroy an FFT plan
 * 
 * @param plan The plan to destroy
 */
void fft_destroy_plan(fft_plan_t plan);

/**
 * @brief Name of the algorithm a plan actually selected
 *
 * Lets callers (and demos, tests, and logs) report the real choice made by
 * select_algorithm() instead of restating its rules and drifting out of date.
 *
 * @param plan The FFT plan
 * @return Pointer to a static string, e.g. "Radix-4"; "Unknown" for a null plan
 */
const char* fft_plan_algorithm_name(fft_plan_t plan);

// Simplified one-shot interface

/**
 * @brief Perform FFT with automatic optimization (v3.0.0)
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
 * @return FFT_OK, or the error code from planning/execution
 */
fft_error_t fft_auto(complex_t* in, complex_t* out, int n, int sign);

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
//
// A process-local cache that maps (transform size, planning flags) to the
// algorithm chosen for it. Plans consult it first, so repeated planning for the
// same size is free. Building a plan with FFT_MEASURE times each candidate
// algorithm and caches the winner; without FFT_MEASURE the cached entry is the
// size heuristic's answer. The cache crosses process boundaries only through
// the two functions below.

/**
 * @brief Export wisdom to string
 * @return Wisdom string (must be freed by caller), or NULL if it cannot be
 *         allocated
 */
char* fft_export_wisdom_to_string(void);

/**
 * @brief Import wisdom from string, replacing the current cache
 *
 * Malformed lines and entries with a non-positive size, an unknown algorithm
 * id, or an algorithm that cannot transform the recorded size are skipped
 * rather than rejected. Entries naming a GPU algorithm are kept even on a
 * machine without one, so a wisdom file remains portable; they are ignored at
 * planning time when the corresponding hardware is absent.
 *
 * @param wisdom Wisdom string in the format produced by
 *               fft_export_wisdom_to_string()
 * @return Number of entries imported (0 for a null argument or no valid lines)
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
