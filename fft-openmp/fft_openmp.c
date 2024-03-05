#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

#define N 8 // Size of FFT
#define PI 3.1415926535897932384626433832795 // Value of Pi

typedef double complex cplx; // Complex number type

// Perform in-place FFT on input array
void fft(cplx x[N], int inverse) {
    int i, j, k, m;
    cplx u, w, t;

    // Bit reversal permutation
    #pragma omp parallel for private(j, t)
    for (i = 0; i < N; i++) {
        j = 0;
        for (k = 1; k < N; k <<= 1) {
            j = (j << 1) | (i & k ? 1 : 0);
        }
        if (i < j) {
            t = x[i];
            x[i] = x[j];
            x[j] = t;
        }
    }

    // Danielson-Lanczos Algorithm
    for (m = 2; m <= N; m <<= 1) {
        double p = inverse ? PI / (m >> 1) : -PI / (m >> 1);
        w = cexp(I * p);
        #pragma omp parallel for private(j, u, t)
        for (k = 0; k < N; k += m) {
            u = 1;
            for (j = 0; j < m >> 1; j++) {
                t = x[k + j + (m >> 1)] * u;
                x[k + j + (m >> 1)] = x[k + j] - t;
                x[k + j] += t;
                u *= w;
            }
        }
    }

    // Scale the output if inverse FFT
    if (inverse) {
        #pragma omp parallel for
        for (i = 0; i < N; i++) {
            x[i] /= N;
        }
    }
}

// Print the FFT output
void showOutput(cplx x[N], const char* title) {
    printf("%s: ", title);
    for (int i = 0; i < N; i++) {
        double real = creal(x[i]);
        double imag = cimag(x[i]);
        if (fabs(real) < 1e-10) real = 0;
        if (fabs(imag) < 1e-10) imag = 0;
        printf("(%.3lf, %.3lfi) ", real, imag);
    }
    printf("\n");
}

int main() {
    clock_t start, end;
    double timeTaken;

    // Input sequence
    cplx input[N] = {1, 1, 1, 1, -1, -1, -1, -1};
    printf("Input: ");
    for (int i = 0; i < N; i++)
        printf("%2.0lf ", creal(input[i]));
    printf("\n");

    // Perform FFT
    start = clock();
    fft(input, 0);
    end = clock();
    timeTaken = ((double)(end - start)) / CLOCKS_PER_SEC;
    showOutput(input, "FFT");

    // Perform inverse FFT
    start = clock();
    fft(input, 1);
    end = clock();
    timeTaken += ((double)(end - start)) / CLOCKS_PER_SEC;
    showOutput(input, "IFFT");

    printf("Time: %.6lf seconds\n", timeTaken);

    return 0;
}