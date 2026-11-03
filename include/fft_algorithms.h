#ifndef FFT_ALGORITHMS_H
#define FFT_ALGORITHMS_H

#include "fft_common.h"

/* 
 * Function declarations for all FFT algorithms
 * This allows algorithms to use each other without circular dependencies
 *
 * Every transform returns an fft_error_t so that a bad length or a failed
 * allocation is reported to the caller instead of terminating the process.
 * Callers that do not care may discard the value. A transform that fails after
 * it has started may leave its buffer partly modified, so do not use the
 * result after an error.
 */

/* Core FFT algorithms */
fft_error_t radix2_dit_fft(complex_t* x, int n, fft_direction dir);
fft_error_t radix2_dif_fft(complex_t* x, int n, fft_direction dir);
fft_error_t radix4_fft(complex_t* x, int n, fft_direction dir);
fft_error_t radix4_omp_fft(complex_t* x, int n, fft_direction dir);
fft_error_t radix4_avx2_fft(complex_t* x, int n, fft_direction dir);
fft_error_t four_step_fft(complex_t* x, int n, fft_direction dir);
fft_error_t split_radix_fft(complex_t* x, int n, fft_direction dir);
fft_error_t bluestein_fft(complex_t* x, int n, fft_direction dir);
fft_error_t mixed_radix_fft(complex_t* x, int n, fft_direction dir);
fft_error_t fft_recursive(complex_t* x, int n, fft_direction dir);
fft_error_t naive_dft(complex_t* x, int n, fft_direction dir);
fft_error_t optimized_dft(complex_t* x, int n, fft_direction dir);

/* Threaded radix-4 variant. Available only when the library was built with
 * FFT_ENABLE_OPENMP and the process has more than one thread available. */
int radix4_omp_available(void);

/* AVX2 radix-4 variant. The build always contains the dispatcher; it reports
 * unavailable on non-x86 targets or when the CPU lacks AVX2. */
int radix4_avx2_available(void);

/* Cache-blocked four-step transform. The kernel accepts any power-of-two
 * length >= 4, but the planner only considers it from 1M points, and only
 * through FFT_MEASURE on the measured hosts. */
int four_step_available(void);

/* Convenience wrappers */
fft_error_t fft_radix2_dit(complex_t* x, int n);
fft_error_t ifft_radix2_dit(complex_t* x, int n);
fft_error_t fft_radix2_dif(complex_t* x, int n);
fft_error_t ifft_radix2_dif(complex_t* x, int n);
fft_error_t fft_radix4(complex_t* x, int n);
fft_error_t ifft_radix4(complex_t* x, int n);
fft_error_t fft_split_radix(complex_t* x, int n);
fft_error_t ifft_split_radix(complex_t* x, int n);
fft_error_t fft_bluestein(complex_t* x, int n);
fft_error_t ifft_bluestein(complex_t* x, int n);
fft_error_t fft_mixed_radix(complex_t* x, int n);
fft_error_t ifft_mixed_radix(complex_t* x, int n);
fft_error_t fft_recursive_forward(complex_t* x, int n);
fft_error_t ifft_recursive(complex_t* x, int n);
fft_error_t dft_forward(complex_t* x, int n);
fft_error_t dft_inverse(complex_t* x, int n);

/* Strided recursion entry point used by mixed_radix_fft; exposed for callers
 * that already hold a strided view of the input. */
void mixed_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir);

#endif /* FFT_ALGORITHMS_H */
