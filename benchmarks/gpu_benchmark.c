#include "../include/fft_algorithms.h"
#include "../include/fft_gpu.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * CPU-vs-GPU benchmark. The GPU backend is not built by default (see
 * CMakeLists.txt), so this reports "GPU not available" and still exercises the
 * CPU path rather than failing.
 */

/* clock() measures processor time, which does not include a blocking wait for
 * the GPU. Use the monotonic wall clock off Windows; MinGW's clock_gettime is
 * not dependable, and the planner already treats Windows clock() as coarse. */
static double get_time(void) {
#if defined(_WIN32)
  return (double)clock() / CLOCKS_PER_SEC;
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

static void run_benchmark(int n) {
    printf("Benchmarking size %d...\n", n);

    complex_t *in = malloc((size_t)n * sizeof(complex_t));
    complex_t *out_cpu = malloc((size_t)n * sizeof(complex_t));
    complex_t *out_gpu = malloc((size_t)n * sizeof(complex_t));

    if (!in || !out_cpu || !out_gpu) {
        printf("Memory allocation failed\n");
        free(in);
        free(out_cpu);
        free(out_gpu);
        return;
    }

    for (int i = 0; i < n; i++) {
        in[i] = (double)rand() / RAND_MAX + I * ((double)rand() / RAND_MAX);
    }

    double start = get_time();
    memcpy(out_cpu, in, (size_t)n * sizeof(complex_t));
    radix2_dit_fft(out_cpu, n, FFT_FORWARD);
    double cpu_time = get_time() - start;
    printf("CPU Time: %.6f s\n", cpu_time);

    if (fft_gpu_available()) {
        /* Warmup */
        fft_gpu_dft_1d(in, out_gpu, n, FFT_FORWARD);

        start = get_time();
        int rc = fft_gpu_dft_1d(in, out_gpu, n, FFT_FORWARD);
        double gpu_time = get_time() - start;

        if (rc != 0) {
            printf("GPU transform failed (rc=%d)\n", rc);
        } else {
            printf("GPU Time: %.6f s\n", gpu_time);
            if (gpu_time > 0) {
                printf("Speedup: %.2fx\n", cpu_time / gpu_time);
            }

            /* Timing is meaningless if the result disagrees with the CPU.
             * MPSGraph computes in float32 while the CPU kernels use double,
             * so the error floor grows with n: measured 4.5e-6 at 1024 and
             * 5.3e-4 at 1048576 on an M1 Pro. 1e-3 covers the range. */
            double max_err = 0.0;
            for (int i = 0; i < n; i++) {
                double e = cabs(out_gpu[i] - out_cpu[i]);
                double s = cabs(out_cpu[i]);
                e /= (s > 1.0 ? s : 1.0);
                if (e > max_err) max_err = e;
            }
            printf("Max relative error vs CPU: %.3e %s\n", max_err,
                   max_err < 1e-3 ? "(PASS)" : "(FAIL)");
        }
    } else {
        printf("GPU not available\n");
    }

    free(in);
    free(out_cpu);
    free(out_gpu);
    printf("----------------------------------------\n");
}

int main(void) {
    printf("Initializing GPU...\n");
    fft_gpu_init(FFT_GPU_AUTO);

    if (fft_gpu_available()) {
        printf("Using GPU: %s\n", fft_gpu_get_device_name());
    } else {
        printf("No GPU detected, running CPU only benchmarks\n");
    }
    printf("----------------------------------------\n");

    srand(42);
    run_benchmark(1024);
    run_benchmark(4096);
    run_benchmark(65536);
    run_benchmark(1048576);

    fft_gpu_cleanup();
    return 0;
}
