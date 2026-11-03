#ifndef FFT_RADIX4_INTERNAL_H
#define FFT_RADIX4_INTERNAL_H

#include "../../include/fft_common.h"

/* Base-4 digit reversal for the radix-4 family. n must be a power of two.
 * Returns 0 on success and -1 if the permutation buffer cannot be allocated. */
int fft_radix4_reorder(complex_t *x, int n);

/* Same permutation using caller-provided scratch of n complex elements. */
int fft_radix4_reorder_into(complex_t *x, int n, complex_t *temp);

/* Radix-4 transform using caller-provided scratch of n complex elements.
 * Used by the four-step kernel to avoid one allocation per row and column. */
fft_error_t radix4_fft_scratch(complex_t *x, int n, fft_direction dir,
                               complex_t *temp);

#endif /* FFT_RADIX4_INTERNAL_H */
