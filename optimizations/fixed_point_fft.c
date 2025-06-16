#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/**
 * Fixed-Point FFT Implementation
 * 
 * Designed for embedded systems and DSPs without floating-point units.
 * Uses Q15 format (1 sign bit, 15 fractional bits) for representing
 * numbers in the range [-1, 1).
 * 
 * Benefits:
 * - No floating-point operations
 * - Predictable precision
 * - Faster on processors without FPU
 * - Lower power consumption
 * 
 * Limitations:
 * - Limited dynamic range
 * - Potential for overflow
 * - Reduced precision compared to floating-point
 */

// Fixed-point configuration
#define FIXED_POINT_BITS 15
#define FIXED_POINT_SCALE (1 << FIXED_POINT_BITS)
#define FIXED_POINT_HALF (1 << (FIXED_POINT_BITS - 1))

// Fixed-point types
typedef int16_t q15_t;
typedef int32_t q31_t;

// Complex number in Q15 format
typedef struct {
    q15_t real;
    q15_t imag;
} complex_q15_t;

// Convert floating-point to Q15
static inline q15_t float_to_q15(float x) {
    if (x >= 1.0f) return 0x7FFF;
    if (x <= -1.0f) return 0x8000;
    return (q15_t)(x * FIXED_POINT_SCALE);
}

// Convert Q15 to floating-point
static inline float q15_to_float(q15_t x) {
    return (float)x / FIXED_POINT_SCALE;
}

// Q15 multiplication with rounding
static inline q15_t q15_mul(q15_t a, q15_t b) {
    q31_t result = ((q31_t)a * b + FIXED_POINT_HALF) >> FIXED_POINT_BITS;
    
    // Saturate on overflow
    if (result > 0x7FFF) return 0x7FFF;
    if (result < -0x8000) return -0x8000;
    
    return (q15_t)result;
}

// Complex multiplication in Q15
static inline complex_q15_t complex_mul_q15(complex_q15_t a, complex_q15_t b) {
    complex_q15_t result;
    
    // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
    q31_t real = (q31_t)a.real * b.real - (q31_t)a.imag * b.imag;
    q31_t imag = (q31_t)a.real * b.imag + (q31_t)a.imag * b.real;
    
    // Round and saturate
    real = (real + FIXED_POINT_HALF) >> FIXED_POINT_BITS;
    imag = (imag + FIXED_POINT_HALF) >> FIXED_POINT_BITS;
    
    if (real > 0x7FFF) real = 0x7FFF;
    if (real < -0x8000) real = -0x8000;
    if (imag > 0x7FFF) imag = 0x7FFF;
    if (imag < -0x8000) imag = -0x8000;
    
    result.real = (q15_t)real;
    result.imag = (q15_t)imag;
    
    return result;
}

// Pre-computed twiddle factors in Q15 format
typedef struct {
    complex_q15_t* factors;
    int size;
} twiddle_table_q15_t;

// Generate twiddle factor table
twiddle_table_q15_t* generate_twiddle_table_q15(int n) {
    twiddle_table_q15_t* table = (twiddle_table_q15_t*)malloc(sizeof(twiddle_table_q15_t));
    table->size = n / 2;
    table->factors = (complex_q15_t*)malloc(table->size * sizeof(complex_q15_t));
    
    for (int k = 0; k < table->size; k++) {
        float angle = -2.0f * 3.14159265358979f * k / n;
        table->factors[k].real = float_to_q15(cosf(angle));
        table->factors[k].imag = float_to_q15(sinf(angle));
    }
    
    return table;
}

void free_twiddle_table_q15(twiddle_table_q15_t* table) {
    free(table->factors);
    free(table);
}

// Bit reversal for fixed-point FFT
void bit_reversal_q15(complex_q15_t* x, int n) {
    int bits = 0;
    int temp = n;
    while (temp >>= 1) bits++;
    
    for (int i = 0; i < n; i++) {
        int reversed = 0;
        for (int j = 0; j < bits; j++) {
            reversed = (reversed << 1) | ((i >> j) & 1);
        }
        
        if (i < reversed) {
            complex_q15_t temp = x[i];
            x[i] = x[reversed];
            x[reversed] = temp;
        }
    }
}

// Fixed-point radix-2 DIT FFT
void fft_radix2_q15(complex_q15_t* x, int n, twiddle_table_q15_t* twiddles) {
    // Check power of 2
    if ((n & (n - 1)) != 0) {
        fprintf(stderr, "FFT size must be power of 2\n");
        return;
    }
    
    // Bit reversal
    bit_reversal_q15(x, n);
    
    // FFT stages
    int stages = 0;
    int temp = n;
    while (temp >>= 1) stages++;
    
    for (int stage = 1; stage <= stages; stage++) {
        int m = 1 << stage;
        int half_m = m >> 1;
        int stride = n >> stage;
        
        for (int k = 0; k < n; k += m) {
            int twiddle_idx = 0;
            
            for (int j = 0; j < half_m; j++) {
                int idx1 = k + j;
                int idx2 = idx1 + half_m;
                
                // Get twiddle factor
                complex_q15_t w = twiddles->factors[twiddle_idx];
                
                // Butterfly computation
                complex_q15_t t = complex_mul_q15(x[idx2], w);
                
                // Prevent overflow in addition/subtraction
                q31_t real_sum = (q31_t)x[idx1].real + t.real;
                q31_t imag_sum = (q31_t)x[idx1].imag + t.imag;
                q31_t real_diff = (q31_t)x[idx1].real - t.real;
                q31_t imag_diff = (q31_t)x[idx1].imag - t.imag;
                
                // Scale down by 2 to prevent overflow (common in fixed-point FFT)
                x[idx1].real = (q15_t)(real_sum >> 1);
                x[idx1].imag = (q15_t)(imag_sum >> 1);
                x[idx2].real = (q15_t)(real_diff >> 1);
                x[idx2].imag = (q15_t)(imag_diff >> 1);
                
                twiddle_idx += stride;
            }
        }
    }
}

// Inverse FFT using conjugate property
void ifft_radix2_q15(complex_q15_t* x, int n, twiddle_table_q15_t* twiddles) {
    // Conjugate input
    for (int i = 0; i < n; i++) {
        x[i].imag = -x[i].imag;
    }
    
    // Forward FFT
    fft_radix2_q15(x, n, twiddles);
    
    // Conjugate output and scale
    int log2n = 0;
    int temp = n;
    while (temp >>= 1) log2n++;
    
    for (int i = 0; i < n; i++) {
        x[i].imag = -x[i].imag;
        // Scale by 1/n (implemented as right shift by log2(n))
        x[i].real >>= log2n;
        x[i].imag >>= log2n;
    }
}

// Block floating-point FFT for better dynamic range
typedef struct {
    complex_q15_t* data;
    int block_exponent;
} block_float_fft_t;

// Normalize block and return exponent
int normalize_block(complex_q15_t* x, int n) {
    // Find maximum absolute value
    int32_t max_val = 0;
    for (int i = 0; i < n; i++) {
        int32_t abs_real = abs(x[i].real);
        int32_t abs_imag = abs(x[i].imag);
        if (abs_real > max_val) max_val = abs_real;
        if (abs_imag > max_val) max_val = abs_imag;
    }
    
    // Calculate shift amount
    int shift = 0;
    while (max_val < 0x4000 && shift < 14) {
        max_val <<= 1;
        shift++;
    }
    
    // Apply normalization
    if (shift > 0) {
        for (int i = 0; i < n; i++) {
            x[i].real <<= shift;
            x[i].imag <<= shift;
        }
    }
    
    return shift;
}

// Performance testing
void test_fixed_point_fft() {
    printf("\nFixed-Point FFT Test:\n");
    printf("====================\n");
    
    int n = 64;
    complex_q15_t* signal = (complex_q15_t*)malloc(n * sizeof(complex_q15_t));
    complex_q15_t* reference = (complex_q15_t*)malloc(n * sizeof(complex_q15_t));
    
    // Generate test signal
    for (int i = 0; i < n; i++) {
        float value = 0.5f * sinf(2 * 3.14159f * 5 * i / n);
        signal[i].real = float_to_q15(value);
        signal[i].imag = 0;
        reference[i] = signal[i];
    }
    
    // Create twiddle table
    twiddle_table_q15_t* twiddles = generate_twiddle_table_q15(n);
    
    // Forward FFT
    fft_radix2_q15(signal, n, twiddles);
    
    // Find peak
    int peak_bin = 0;
    int32_t peak_mag = 0;
    
    for (int i = 0; i < n/2; i++) {
        int32_t mag = (int32_t)signal[i].real * signal[i].real + 
                      (int32_t)signal[i].imag * signal[i].imag;
        if (mag > peak_mag) {
            peak_mag = mag;
            peak_bin = i;
        }
    }
    
    printf("Peak at bin %d (expected: 5)\n", peak_bin);
    
    // Inverse FFT
    ifft_radix2_q15(signal, n, twiddles);
    
    // Calculate error
    int32_t error = 0;
    for (int i = 0; i < n; i++) {
        error += abs(signal[i].real - reference[i].real);
        error += abs(signal[i].imag - reference[i].imag);
    }
    
    printf("Reconstruction error: %d (Q15 units)\n", error);
    printf("Average error per sample: %.6f\n", (float)error / (n * 2 * FIXED_POINT_SCALE));
    
    free(signal);
    free(reference);
    free_twiddle_table_q15(twiddles);
}

// Overflow detection and prevention
void demonstrate_overflow_handling() {
    printf("\n\nOverflow Handling:\n");
    printf("==================\n");
    
    // Example of potential overflow
    q15_t a = 0x6000;  // 0.75
    q15_t b = 0x6000;  // 0.75
    
    // Without saturation: 0.75 * 0.75 = 0.5625
    q31_t result = (q31_t)a * b;
    printf("Raw multiplication: 0x%08X\n", result);
    
    // With proper scaling
    q15_t scaled_result = (q15_t)((result + FIXED_POINT_HALF) >> FIXED_POINT_BITS);
    printf("Scaled result: 0x%04X (%.4f)\n", scaled_result, q15_to_float(scaled_result));
    
    // Demonstrate saturation
    a = 0x7FFF;  // Maximum positive
    b = 0x4000;  // 0.5
    q15_t saturated = q15_mul(a, b);
    printf("Saturated multiplication: %.4f * %.4f = %.4f\n",
           q15_to_float(a), q15_to_float(b), q15_to_float(saturated));
}

// Main demonstration
int main() {
    printf("Fixed-Point FFT Implementation\n");
    printf("==============================\n");
    
    // Display Q15 format info
    printf("\nQ15 Format:\n");
    printf("-----------\n");
    printf("Range: [-1.0, 0.999969] (%.6f)\n", q15_to_float(0x7FFF));
    printf("Resolution: %.6f\n", 1.0f / FIXED_POINT_SCALE);
    printf("Bits: 1 sign + 15 fractional\n");
    
    // Test conversions
    printf("\nConversion Examples:\n");
    printf("-------------------\n");
    float test_values[] = {0.0f, 0.5f, -0.5f, 0.999f, -1.0f};
    
    for (int i = 0; i < 5; i++) {
        q15_t fixed = float_to_q15(test_values[i]);
        float back = q15_to_float(fixed);
        printf("%.3f -> 0x%04X -> %.6f\n", test_values[i], fixed & 0xFFFF, back);
    }
    
    // Run FFT test
    test_fixed_point_fft();
    
    // Demonstrate overflow handling
    demonstrate_overflow_handling();
    
    // Performance comparison
    printf("\n\nPerformance Considerations:\n");
    printf("===========================\n");
    printf("Operation       | Floating-Point | Fixed-Point\n");
    printf("----------------|----------------|-------------\n");
    printf("Addition        | 1 cycle        | 1 cycle\n");
    printf("Multiplication  | 3-5 cycles     | 1-2 cycles\n");
    printf("Memory (complex)| 8 bytes        | 4 bytes\n");
    printf("Dynamic range   | ~10^38         | [-1, 1]\n");
    
    // Best practices
    printf("\n\nBest Practices:\n");
    printf("===============\n");
    printf("1. Pre-compute and store twiddle factors\n");
    printf("2. Use block floating-point for better range\n");
    printf("3. Scale at each stage to prevent overflow\n");
    printf("4. Use saturation arithmetic\n");
    printf("5. Test with worst-case inputs\n");
    printf("6. Consider Q31 for intermediate results\n");
    printf("7. Profile on target hardware\n");
    
    return 0;
}
