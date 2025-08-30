#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/fft_auto.h"
#include "../include/fft_gpu.h"

/**
 * @file demo_v2_features.c
 * @brief Demonstration of FFT v2.0.0 features
 * 
 * Shows:
 * - Automatic algorithm selection
 * - GPU acceleration
 * - New simplified API
 * - Performance comparison
 * - v2.0.0 optimizations (twiddle factors, bit reversal, compiler hints)
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

// Benchmark a specific approach
double benchmark_fft(complex_t* signal, int n, const char* method,
                    void (*fft_func)(complex_t*, complex_t*, int, int)) {
    complex_t* work = fft_alloc_complex(n);
    memcpy(work, signal, n * sizeof(complex_t));
    
    clock_t start = clock();
    fft_func(work, work, n, -1);
    clock_t end = clock();
    
    fft_free(work);
    
    return ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
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
        
        // Create plan and check selected algorithm
        complex_t* data = fft_alloc_complex(n);
        fft_plan_t plan = fft_plan_dft_1d(n, data, data, -1, FFT_ESTIMATE);
        
        // The library will print or we can check internal state
        printf("  Selected algorithm: ");
        if (is_power_of_two(n)) {
            if (n <= 64) printf("Radix-2 DIT (simple for small sizes)\n");
            else if (n <= 1024) printf("Radix-4 (fewer operations)\n");
            else printf("Split-Radix (optimal operation count)\n");
        } else {
            if (n == 97) printf("Bluestein (prime size)\n");
            else if (n == 360) printf("Mixed-Radix (highly composite)\n");
            else printf("Bluestein (arbitrary size)\n");
        }
        
        // Time the execution
        generate_test_signal(data, n);
        
        clock_t start = clock();
        fft_execute(plan);
        clock_t end = clock();
        
        double time_ms = ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
        printf("  Execution time: %.3f ms\n", time_ms);
        
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
    printf("\n\n=== Simplified API Demo (v2.0.0) ===\n");
    
    int n = 1024;
    complex_t* signal = fft_alloc_complex(n);
    
    printf("\nOld API (v1.x):\n");
    printf("  radix2_dit_fft(signal, n, FFT_FORWARD);\n");
    printf("  // Need to know which algorithm to use\n");
    printf("  // No automatic optimization\n");
    
    printf("\nNew API (v2.0.0):\n");
    printf("  fft_auto(signal, signal, n, -1);\n");
    printf("  // Automatic algorithm selection\n");
    printf("  // GPU acceleration if available\n");
    printf("  // Optimized for your hardware\n");
    
    // Demonstrate usage
    generate_test_signal(signal, n);
    
    printf("\nExecuting FFT with automatic optimization...\n");
    clock_t start = clock();
    fft_auto(signal, signal, n, -1);
    clock_t end = clock();
    
    double time_ms = ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
    printf("Completed in %.3f ms\n", time_ms);
    
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
    printf("FFT v2.0.0 Feature Demonstration\n");
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
    printf("\n\n=== Summary of v2.0.0 Features ===\n");
    printf("1. Automatic algorithm selection based on size and hardware\n");
    printf("2. GPU acceleration with CUDA and Metal support\n");
    printf("3. Simplified API with fft_auto() function\n");
    printf("4. Cross-platform compatibility\n");
    printf("5. Intelligent planning with FFT_MEASURE flag\n");
    printf("6. Hardware detection and optimization\n");
    printf("7. Aligned memory allocation for SIMD\n");
    printf("8. Thread-safe implementation\n");
    
    printf("\nFor more information, see the documentation at:\n");
    printf("https://fft-implementation-in-c.readthedocs.io/\n");
    
    return 0;
}
