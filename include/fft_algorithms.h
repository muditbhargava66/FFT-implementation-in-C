#ifndef FFT_ALGORITHMS_H
#define FFT_ALGORITHMS_H

#include "fft_common.h"

/* 
 * Function declarations for all FFT algorithms
 * This allows algorithms to use each other without circular dependencies
 */

/* Core FFT algorithms */
void radix2_dit_fft(complex_t* x, int n, fft_direction dir);
void radix2_dif_fft(complex_t* x, int n, fft_direction dir);
void radix4_fft(complex_t* x, int n, fft_direction dir);
void split_radix_fft(complex_t* x, int n, fft_direction dir);
void bluestein_fft(complex_t* x, int n, fft_direction dir);
void mixed_radix_fft(complex_t* x, int n, fft_direction dir);
void fft_recursive(complex_t* x, int n, fft_direction dir);
void naive_dft(complex_t* x, int n, fft_direction dir);
void optimized_dft(complex_t* x, int n, fft_direction dir);

/* Convenience wrappers */
void fft_radix2_dit(complex_t* x, int n);
void ifft_radix2_dit(complex_t* x, int n);
void fft_radix2_dif(complex_t* x, int n);
void ifft_radix2_dif(complex_t* x, int n);
void fft_radix4(complex_t* x, int n);
void ifft_radix4(complex_t* x, int n);
void fft_split_radix(complex_t* x, int n);
void ifft_split_radix(complex_t* x, int n);
void fft_bluestein(complex_t* x, int n);
void ifft_bluestein(complex_t* x, int n);
void fft_mixed_radix(complex_t* x, int n);
void ifft_mixed_radix(complex_t* x, int n);
void fft_recursive_forward(complex_t* x, int n);
void ifft_recursive(complex_t* x, int n);
void dft_forward(complex_t* x, int n);
void dft_inverse(complex_t* x, int n);

/* Recursive mixed-radix helper (exposed for algorithms that need it) */
void mixed_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir);

#endif /* FFT_ALGORITHMS_H */
