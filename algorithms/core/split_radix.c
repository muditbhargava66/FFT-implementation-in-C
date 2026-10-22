#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Split-Radix FFT Implementation
 *
 * A split-radix step expresses an n-point DFT as one n/2-point DFT plus two
 * n/4-point DFTs, which is the decomposition with the lowest known arithmetic
 * cost for power-of-two lengths.
 *
 * Starting from the radix-2 split X[k] = E[k] + w_n^k * O[k], where E and O are
 * the DFTs of the even and odd samples, the even branch is split once more:
 *
 *   O[k] = DFT_{n/2}( x[2m+1] )
 *   A[k] = DFT_{n/4}( x[4l] )
 *   B[k] = DFT_{n/4}( x[4l+2] )
 *   E[k] = A[k mod n/4] + w_{n/2}^k * B[k mod n/4]
 *   X[k]     = E[k] + w_n^k * O[k]
 *   X[k+n/2] = E[k] - w_n^k * O[k]
 *
 * for k = 0..n/2-1, with the n=2 butterfly [[1,1],[1,-1]] as the base case (it
 * is its own inverse, so the recursion carries `dir` only through the twiddles).
 *
 * The kernel is unnormalised for both directions; split_radix_fft applies the
 * single 1/n scale at the end.
 *
 * Time Complexity: O(n log n), with the split-radix operation counts.
 */

/*
 * sp is a bump cursor into a caller-owned pool; children allocate above this
 * node's buffers and roll the cursor back when they return, so peak usage is
 * 2.5n + 2.5(n/2) + ... < 5n. Without the pool this kernel issued ~n mallocs
 * per transform, and the allocator - not the arithmetic - was the whole cost.
 */
static int sr_recursion(complex_t* out, const complex_t* in, int n,
                        fft_direction dir, complex_t** sp) {
    if (n == 1) {
        out[0] = in[0];
        return 0;
    }
    if (n == 2) {
        out[0] = in[0] + in[1];
        out[1] = in[0] - in[1];
        return 0;
    }

    int h = n / 2, q = n / 4;

    /* Gather inputs and their transform outputs in separate regions: sr_recursion
     * writes its out argument while later stages still read in, so they must not
     * overlap. */
    complex_t* base = *sp;
    complex_t* in_o = base;                    /* h */
    complex_t* out_o = in_o + h;               /* h */
    complex_t* in_a = out_o + h;               /* q */
    complex_t* out_a = in_a + q;               /* q */
    complex_t* in_b = out_a + q;               /* q */
    complex_t* out_b = in_b + q;               /* q */
    complex_t* e = out_b + q;                  /* h */
    if ((e + h) > *sp + (5 * n)) return -1;    /* pool exhausted */

    for (int m = 0; m < h; m++) in_o[m] = in[2 * m + 1];
    for (int l = 0; l < q; l++) in_a[l] = in[4 * l];
    for (int l = 0; l < q; l++) in_b[l] = in[4 * l + 2];

    complex_t* mark = e + h;
    int rc = sr_recursion(out_o, in_o, h, dir, &mark);
    if (rc == 0) rc = sr_recursion(out_a, in_a, q, dir, &mark);
    if (rc == 0) rc = sr_recursion(out_b, in_b, q, dir, &mark);
    if (rc != 0) return rc;

    /* E[k] = A[k mod q] + w_h^k * B[k mod q], k = 0..h-1 */
    complex_t wh = twiddle_factor(1, h, dir);
    complex_t wk = 1.0;
    for (int k = 0; k < h; k++) {
        int kq = k & (q - 1);          /* k mod q; q is a power of two */
        e[k] = out_a[kq] + out_b[kq] * wk;
        wk *= wh;
    }

    /* X[k] = E[k] + w_n^k * O[k] and X[k+h] = E[k] - w_n^k * O[k] */
    complex_t wn = twiddle_factor(1, n, dir);
    wk = 1.0;
    for (int k = 0; k < h; k++) {
        complex_t t = out_o[k] * wk;
        out[k] = e[k] + t;
        out[k + h] = e[k] - t;
        wk *= wn;
    }

    return 0;
}

// Main split-radix FFT function
fft_error_t split_radix_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    FFT_REQUIRE_POW2(n);
    if (n == 1) return FFT_OK;

    /* The kernel reads in[] while writing out[], so it cannot run in place. */
    complex_t* pool = (complex_t*)malloc(((size_t)6 * n) * sizeof(complex_t));
    if (!pool) return FFT_ERR_ALLOC;
    memcpy(pool, x, (size_t)n * sizeof(complex_t));

    complex_t* cursor = pool + n;        /* pool[0..n-1] holds the input copy */
    int exhausted = sr_recursion(x, pool, n, dir, &cursor);
    free(pool);
    if (exhausted) return FFT_ERR_ALLOC;

    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
    return FFT_OK;
}

// Wrapper functions
fft_error_t fft_split_radix(complex_t* x, int n) {
    return split_radix_fft(x, n, FFT_FORWARD);
}

fft_error_t ifft_split_radix(complex_t* x, int n) {
    return split_radix_fft(x, n, FFT_INVERSE);
}

// Operation counting for comparison
typedef struct {
    int adds;
    int muls;
} op_count_t;

op_count_t count_split_radix_ops(int n) {
    op_count_t count = {0, 0};
    if (n <= 1) return count;
    
    // Based on theoretical analysis
    int log2n = log2_int(n);
    count.muls = (n * log2n - 3 * n + 4) / 3;
    count.adds = n * log2n;
    
    return count;
}

op_count_t count_radix2_ops(int n) {
    op_count_t count = {0, 0};
    int log2n = log2_int(n);
    count.muls = n * log2n / 2;
    count.adds = n * log2n;
    return count;
}

#ifndef LIB_BUILD

// Performance demonstration
int main() {
    printf("Split-Radix FFT Implementation\n");
    printf("==============================\n\n");
    
    // Operation count comparison
    printf("Theoretical Operation Counts:\n");
    printf("Size\tRadix-2 Muls\tSplit-Radix Muls\tReduction\n");
    printf("----\t------------\t----------------\t---------\n");
    
    for (int n = 8; n <= 1024; n *= 2) {
        op_count_t r2_ops = count_radix2_ops(n);
        op_count_t sr_ops = count_split_radix_ops(n);
        double reduction = (1.0 - (double)sr_ops.muls / r2_ops.muls) * 100;
        
        printf("%d\t%d\t\t%d\t\t\t%.1f%%\n", 
               n, r2_ops.muls, sr_ops.muls, reduction);
    }
    
    // Performance comparison
    printf("\n\nPerformance Comparison:\n");
    printf("Size\tRadix-2\t\tRadix-4\t\tSplit-Radix\n");
    printf("----\t-------\t\t-------\t\t-----------\n");
    
    int sizes[] = {128, 512, 2048, 8192};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        complex_t* x_r2 = allocate_complex_array(n);
        complex_t* x_r4 = allocate_complex_array(n);
        complex_t* x_sr = allocate_complex_array(n);
        
        // Generate identical random data
        for (int i = 0; i < n; i++) {
            double real = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double imag = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            x_r2[i] = x_r4[i] = x_sr[i] = real + I * imag;
        }
        
        fft_timer_t timer;
        double time_r2, time_r4, time_sr;
        
        // Time radix-2
        timer_start(&timer);
        radix2_dit_fft(x_r2, n, FFT_FORWARD);
        timer_stop(&timer);
        time_r2 = timer.elapsed_ms;
        
        // Time radix-4
        timer_start(&timer);
        fft_radix4(x_r4, n);
        timer_stop(&timer);
        time_r4 = timer.elapsed_ms;
        
        // Time split-radix
        timer_start(&timer);
        fft_split_radix(x_sr, n);
        timer_stop(&timer);
        time_sr = timer.elapsed_ms;
        
        printf("%d\t%.3f ms\t%.3f ms\t%.3f ms\n",
               n, time_r2, time_r4, time_sr);
        
        // Verify correctness
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = cabs(x_r2[i] - x_sr[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        if (max_diff > 1e-10) {
            printf("  WARNING: Results differ by %.2e\n", max_diff);
        }
        
        free_complex_array(x_r2);
        free_complex_array(x_r4);
        free_complex_array(x_sr);
    }
    
    // Accuracy test
    printf("\n\nAccuracy Test:\n");
    int n = 32;
    complex_t* signal = allocate_complex_array(n);
    complex_t* original = allocate_complex_array(n);
    
    // Generate test signal
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 3 * i / n) + 0.5 * cos(2 * PI * 7 * i / n);
        original[i] = signal[i];
    }
    
    // Forward and inverse transform
    fft_split_radix(signal, n);
    ifft_split_radix(signal, n);
    
    // Check reconstruction
    double error = 0.0;
    for (int i = 0; i < n; i++) {
        error += cabs(signal[i] - original[i]);
    }
    
    printf("Reconstruction error: %.2e %s\n", error / n,
           error / n < 1e-14 ? "(PASS)" : "(FAIL)");
    
    free_complex_array(signal);
    free_complex_array(original);
    
    return 0;
}

#endif /* LIB_BUILD */
