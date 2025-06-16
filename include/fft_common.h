#ifndef FFT_COMMON_H
#define FFT_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

// Constants
#define PI 3.14159265358979323846
#define TWO_PI (2.0 * PI)

// Complex number type
typedef double complex complex_t;

// FFT direction
typedef enum {
    FFT_FORWARD = -1,
    FFT_INVERSE = 1
} fft_direction;

// Common utility functions
static inline int is_power_of_two(int n) {
    return n > 0 && (n & (n - 1)) == 0;
}

static inline int next_power_of_two(int n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

static inline int log2_int(int n) {
    int log = 0;
    while (n >>= 1) log++;
    return log;
}

// Bit reversal function
static inline unsigned int bit_reverse(unsigned int x, int log2n) {
    unsigned int r = 0;
    for (int i = 0; i < log2n; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

// Memory allocation helpers
static inline complex_t* allocate_complex_array(int n) {
    return (complex_t*)calloc(n, sizeof(complex_t));
}

static inline void free_complex_array(complex_t* arr) {
    if (arr) free(arr);
}

// Twiddle factor computation
static inline complex_t twiddle_factor(int k, int n, fft_direction dir) {
    double angle = dir * TWO_PI * k / n;
    return cexp(I * angle);
}

// Performance timing
typedef struct {
    clock_t start;
    clock_t end;
    double elapsed_ms;
} timer_t;

static inline void timer_start(timer_t* timer) {
    timer->start = clock();
}

static inline void timer_stop(timer_t* timer) {
    timer->end = clock();
    timer->elapsed_ms = ((double)(timer->end - timer->start) / CLOCKS_PER_SEC) * 1000.0;
}

// Error checking
#define CHECK_NULL(ptr, msg) \
    if (!(ptr)) { \
        fprintf(stderr, "Error: %s\n", msg); \
        exit(EXIT_FAILURE); \
    }

#define CHECK_POWER_OF_TWO(n) \
    if (!is_power_of_two(n)) { \
        fprintf(stderr, "Error: Size %d is not a power of two\n", n); \
        exit(EXIT_FAILURE); \
    }

// Complex number utilities
static inline void print_complex(complex_t c) {
    double real = creal(c);
    double imag = cimag(c);
    if (fabs(real) < 1e-10) real = 0;
    if (fabs(imag) < 1e-10) imag = 0;
    printf("(%.3f, %.3fi)", real, imag);
}

static inline void print_complex_array(const char* label, complex_t* arr, int n) {
    printf("%s: ", label);
    for (int i = 0; i < n; i++) {
        print_complex(arr[i]);
        printf(" ");
    }
    printf("\n");
}

// Signal generation utilities
static inline void generate_sine_wave(complex_t* signal, int n, double freq, double fs) {
    for (int i = 0; i < n; i++) {
        signal[i] = sin(TWO_PI * freq * i / fs) + 0 * I;
    }
}

static inline void generate_square_wave(complex_t* signal, int n, double freq, double fs) {
    int period = (int)(fs / freq);
    for (int i = 0; i < n; i++) {
        signal[i] = (i % period < period / 2) ? 1.0 : -1.0;
    }
}

static inline void generate_impulse(complex_t* signal, int n) {
    memset(signal, 0, n * sizeof(complex_t));
    signal[0] = 1.0;
}

// Magnitude and phase utilities
static inline double* compute_magnitude(complex_t* fft_result, int n) {
    double* mag = (double*)malloc(n * sizeof(double));
    CHECK_NULL(mag, "Failed to allocate magnitude array");
    
    for (int i = 0; i < n; i++) {
        mag[i] = cabs(fft_result[i]);
    }
    return mag;
}

static inline double* compute_phase(complex_t* fft_result, int n) {
    double* phase = (double*)malloc(n * sizeof(double));
    CHECK_NULL(phase, "Failed to allocate phase array");
    
    for (int i = 0; i < n; i++) {
        phase[i] = carg(fft_result[i]);
    }
    return phase;
}

static inline double* compute_power_spectrum(complex_t* fft_result, int n) {
    double* power = (double*)malloc(n * sizeof(double));
    CHECK_NULL(power, "Failed to allocate power spectrum array");
    
    for (int i = 0; i < n; i++) {
        double mag = cabs(fft_result[i]);
        power[i] = mag * mag / n;
    }
    return power;
}

#endif // FFT_COMMON_H
