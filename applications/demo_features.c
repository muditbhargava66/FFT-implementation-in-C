#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fft_auto.h"
#include "fft_gpu.h"

/**
 * @file demo_features.c
 * @brief Demonstration of FFT features
 * 
 * Shows:
 * - Automatic algorithm selection
 * - GPU acceleration
 * - New simplified API
 * - Performance comparison
 * - Optimization features (twiddle factors, bit reversal, compiler hints)
 */

// Generate test signal
void generate_test_signal(complex_t* signal, int n) {
    double fs = 1000.0;  // Sample rate
    
    for (int i = 0; i < n; i++) {
        double t = i / fs;
        // Multi-frequency signal
        signal[i] = sin(2 * PI * 50 * t) +       // 50 Hz
                   0.5 * sin(2 * PI * 120 * t) +  // 120 Hz  
                   0.25 * sin(2 * PI * 250 * t);  // 250 Hz
        
        // Add some noise
        signal[i] += 0.1 * ((double)rand() / RAND_MAX - 0.5);
    }
}



// Demo automatic algorithm selection
void demo_automatic_selection(void) {
    printf("\n=== Automatic Algorithm Selection Demo ===\n");
    
    int sizes[] = {64, 256, 1024, 4096, 16384, 97, 360, 1000};
    const char* size_types[] = {"small pow2", "medium pow2", "large pow2", 
                                "larger pow2", "very large pow2", "prime", 
                                "composite", "non-pow2"};
    
    for (int i = 0; i < 8; i++) {
        int n = sizes[i];
        printf("\nSize %d (%s):\n", n, size_types[i]);
        
        complex_t* data = fft_alloc_complex(n);
        fft_plan_t plan = fft_plan_dft_1d(n, data, data, -1, FFT_ESTIMATE);
        if (!plan) {
            printf("  Planning failed for n=%d\n", n);
            fft_free(data);
            continue;
        }

        /* Ask the plan what it chose. This used to be a hand-written copy of
         * the selection rules, which had already drifted from the library. */
        printf("  Selected algorithm: %s\n", fft_plan_algorithm_name(plan));
        
        generate_test_signal(data, n);
        
        /* clock() advances about once a millisecond on Windows, so one small
         * transform measures as 0.000 ms. The batch is grown until the clock
         * records at least ~50 ticks, which keeps the quantisation error under
         * the printed precision. */
        int reps = 32;
        double time_ms = 0.0;
        for (;;) {
            clock_t start = clock();
            for (int r = 0; r < reps; r++)
                fft_execute(plan);
            double ticks = (double)(clock() - start);
            time_ms = ticks * 1000.0 / ((double)CLOCKS_PER_SEC * (double)reps);
            if (ticks >= 50.0 || reps >= 2000000) break;
            int next = (ticks > 0.0) ? (int)(reps * (50.0 / ticks)) : reps * 8;
            reps = (next > reps) ? next : reps + 1;
        }
        printf("  Execution time: %.4f ms (%d runs)\n", time_ms, reps);
        
        fft_destroy_plan(plan);
        fft_free(data);
    }
}

// Demo GPU acceleration
void demo_gpu_acceleration(void) {
    printf("\n\n=== GPU Acceleration Demo ===\n");
    
    // Check GPU availability
    if (!fft_gpu_available()) {
        printf("No GPU available. Skipping GPU demo.\n");
        return;
    }
    
    // Initialize GPU
    if (fft_gpu_init(FFT_GPU_AUTO) != 0) {
        printf("Failed to initialize GPU.\n");
        return;
    }
    
    printf("GPU Device: %s\n", fft_gpu_get_device_name());
    
    size_t total_mem, avail_mem;
    fft_gpu_get_memory_info(&total_mem, &avail_mem);
    printf("GPU Memory: %.1f GB total, %.1f GB available\n", 
           total_mem / 1e9, avail_mem / 1e9);
    
    // Benchmark CPU vs GPU
    printf("\nCPU vs GPU Performance:\n");
    printf("Size    | CPU Time  | GPU Time  | Speedup\n");
    printf("--------|-----------|-----------|--------\n");
    
    int sizes[] = {1024, 4096, 16384, 65536, 262144};
    
    for (int i = 0; i < 5; i++) {
        int n = sizes[i];
        complex_t* signal = fft_alloc_complex(n);
        generate_test_signal(signal, n);
        
        // CPU timing
        fft_plan_t cpu_plan = fft_plan_dft_1d(n, signal, signal, -1, 0);
        clock_t start = clock();
        for (int j = 0; j < 10; j++) {
            fft_execute(cpu_plan);
        }
        clock_t end = clock();
        double cpu_time = ((double)(end - start) / CLOCKS_PER_SEC) * 100.0; // ms per FFT
        
        // GPU timing
        fft_plan_t gpu_plan = fft_plan_dft_1d(n, signal, signal, -1, FFT_PREFER_GPU);
        start = clock();
        for (int j = 0; j < 10; j++) {
            fft_execute(gpu_plan);
        }
        end = clock();
        double gpu_time = ((double)(end - start) / CLOCKS_PER_SEC) * 100.0;
        
        printf("%-7d | %7.2f ms | %7.2f ms | %.1fx\n", 
               n, cpu_time, gpu_time, cpu_time / gpu_time);
        
        fft_destroy_plan(cpu_plan);
        fft_destroy_plan(gpu_plan);
        fft_free(signal);
    }
    
    fft_gpu_cleanup();
}

// Demo new simplified API
void demo_simplified_api(void) {
    printf("\n\n=== Simplified API Demo ===\n");
    
    int n = 1024;
    complex_t* signal = fft_alloc_complex(n);
    
    printf("\nLegacy API:\n");
    printf("  radix2_dit_fft(signal, n, FFT_FORWARD);\n");
    printf("  // Need to know which algorithm to use\n");
    printf("  // No automatic optimization\n");
    
    printf("\nModern API:\n");
    printf("  fft_auto(signal, signal, n, -1);\n");
    printf("  // Automatic algorithm selection\n");
    printf("  // Best CPU kernel for this machine\n");
    printf("  // GPU only with FFT_PREFER_GPU when a GPU build is available\n");
    
    // Demonstrate usage
    generate_test_signal(signal, n);
    
    printf("\nExecuting FFT with automatic optimization...\n");
    /* A single 1024-point transform is far below clock()'s ~1 ms resolution on
     * Windows, so reporting it here printed 0.000 ms. The measured-batch
     * numbers are in the selection demo above and in benchmark_all. */
    fft_auto(signal, signal, n, -1);
    printf("Completed.\n");
    
    // Show some results
    printf("\nSpectrum peaks:\n");
    double* mag = compute_magnitude(signal, n);
    for (int i = 0; i < n/2; i++) {
        if (mag[i] > 100) {  // Significant peaks
            double freq = i * 1000.0 / n;
            printf("  %.1f Hz: magnitude %.1f\n", freq, mag[i]);
        }
    }
    
    free(mag);
    fft_free(signal);
}

// Demo hardware detection
void demo_hardware_detection(void) {
    printf("\n\n=== Hardware Detection ===\n");
    
    unsigned caps = fft_get_hardware_capabilities();
    
    printf("CPU Features:\n");
    if (caps & FFT_HW_CPU_SSE)    printf("  ✓ SSE\n");
    if (caps & FFT_HW_CPU_AVX)    printf("  ✓ AVX\n");
    if (caps & FFT_HW_CPU_AVX2)   printf("  ✓ AVX2\n");
    if (caps & FFT_HW_CPU_AVX512) printf("  ✓ AVX-512\n");
    if (caps & FFT_HW_CPU_NEON)   printf("  ✓ ARM NEON\n");
    
    printf("\nGPU Features:\n");
    if (caps & FFT_HW_GPU_CUDA)   printf("  ✓ NVIDIA CUDA\n");
    if (caps & FFT_HW_GPU_MPS)    printf("  ✓ Apple Metal\n");
    if (caps & FFT_HW_GPU_OPENCL) printf("  ✓ OpenCL\n");
    
    if (!(caps & (FFT_HW_GPU_CUDA | FFT_HW_GPU_MPS | FFT_HW_GPU_OPENCL))) {
        printf("  ✗ No GPU acceleration available\n");
    }
}

// Main demo program
int main() {
    printf("FFT Feature Demonstration\n");
    printf("================================\n");
    
    // Show version
    printf("Library Version: %s\n", fft_version());
    
    // Seed random number generator
    srand(time(NULL));
    
    // Run demonstrations
    demo_hardware_detection();
    demo_automatic_selection();
    demo_simplified_api();
    demo_gpu_acceleration();
    
    // Summary
    printf("\n\n=== Summary of Features ===\n");
    printf("1. Automatic algorithm selection based on size and hardware\n");
    printf("2. GPU backends in source, excluded from the default build (CPU-only)\n");
    printf("3. Simplified API with fft_auto() function\n");
    printf("4. Cross-platform compatibility\n");
    printf("5. Intelligent planning with FFT_MEASURE flag\n");
    printf("6. Hardware detection and optimization\n");
    printf("7. Real-to-complex and complex-to-real FFT (r2c/c2r)\n");
    printf("8. 2D FFT planning with row-column decomposition\n");
    printf("9. Wisdom cache for persistent algorithm selection\n");
    printf("10. Single-threaded planner; see the README's thread-safety note\n");
    
    printf("\nFor more information, see the documentation at:\n");
    printf("https://fft-implementation-in-c.readthedocs.io/\n");
    
    return 0;
}
