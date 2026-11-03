#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"
#include "radix4_internal.h"

/*
 * Four-step FFT for large power-of-two lengths.
 *
 * N = N1*N2; the input is an N1-by-N2 row-major matrix. Column FFTs, a twiddle
 * pass and row FFTs operate on cache-sized blocks, then one transpose returns
 * natural order. The twiddle pass W_N^(k1*c) is what makes the decomposition
 * exact; without it the two sets of small DFTs do not compose into DFT_N.
 *
 * Radix-4 sub-transforms share one scratch buffer, so blocking does not add an
 * allocation per row or column.
 */

#define FOUR_STEP_TILE 16

int four_step_available(void) { return 1; }

fft_error_t four_step_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;
    if (n < 4) return radix2_dit_fft(x, n, dir);

    const int log2n = log2_int(n);
    const int k1 = (log2n + 1) / 2;
    const int k2 = log2n - k1;
    const int n1 = 1 << k1;
    const int n2 = 1 << k2;

    const int max_sub = (n1 > n2) ? n1 : n2;
    complex_t* temp = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    complex_t* col = (complex_t*)malloc((size_t)n1 * sizeof(complex_t));
    complex_t* sub = (complex_t*)malloc((size_t)max_sub * sizeof(complex_t));
    if (!temp || !col || !sub) {
        free(temp);
        free(col);
        free(sub);
        return FFT_ERR_ALLOC;
    }

    /* Gather each column into a contiguous buffer so radix-4 can transform
     * it. The result is stored as A[k1][c]; the row pass reads it directly. */
    for (int c = 0; c < n2; c++) {
        for (int r = 0; r < n1; r++)
            col[r] = x[r * n2 + c];

        fft_error_t err = radix4_fft_scratch(col, n1, dir, sub);
        if (err != FFT_OK) {
            free(temp);
            free(col);
            free(sub);
            return err;
        }

        for (int r = 0; r < n1; r++)
            temp[r * n2 + c] = col[r];
    }

    /* Twiddle A[k1][c] by W_N^(k1*c). */
    for (int r = 0; r < n1; r++) {
        const complex_t step = twiddle_factor(r, n, dir);
        complex_t w = 1.0;
        for (int c = 0; c < n2; c++) {
            temp[r * n2 + c] *= w;
            w *= step;
        }
    }

    /* Row transforms. */
    for (int r = 0; r < n1; r++) {
        fft_error_t err = radix4_fft_scratch(&temp[r * n2], n2, dir, sub);
        if (err != FFT_OK) {
            free(temp);
            free(col);
            free(sub);
            return err;
        }
    }

    /* Transpose back to [k2][k1]. Tiling keeps both the source and
     * destination lines resident for a 16x16 block. */
    for (int i0 = 0; i0 < n1; i0 += FOUR_STEP_TILE) {
        const int i_end = (i0 + FOUR_STEP_TILE < n1) ? i0 + FOUR_STEP_TILE : n1;
        for (int j0 = 0; j0 < n2; j0 += FOUR_STEP_TILE) {
            const int j_end =
                (j0 + FOUR_STEP_TILE < n2) ? j0 + FOUR_STEP_TILE : n2;
            for (int i = i0; i < i_end; i++)
                for (int j = j0; j < j_end; j++)
                    x[j * n1 + i] = temp[i * n2 + j];
        }
    }

    free(temp);
    free(col);
    free(sub);
    return FFT_OK;
}
