#include "../../include/fft_auto.h"
#include "../../include/fft_algorithms.h"
#include "../../include/fft_gpu.h"
#include <stdlib.h>
#include <string.h>

/**
 * @file fft_auto.c
 * @brief Automatic FFT algorithm selection and planning (v2.0.0)
 * 
 * This module provides intelligent algorithm selection based on:
 * - Transform size and factorization
 * - Hardware capabilities
 * - Memory constraints
 * - Performance measurements
 */

// Internal plan structure
struct fft_plan {
    int n;                      // Transform size
    complex_t* in;             // Input buffer
    complex_t* out;            // Output buffer  
    fft_direction dir;         // Transform direction
    unsigned flags;            // Planning flags
    
    // Algorithm selection
    enum {
        ALGO_RADIX2_DIT,
        ALGO_RADIX2_DIF,
        ALGO_RADIX4,
        ALGO_SPLIT_RADIX,
        ALGO_BLUESTEIN,
        ALGO_MIXED_RADIX,
        ALGO_GPU_CUDA,
        ALGO_GPU_MPS
    } algorithm;
    
    // Precomputed data
    complex_t* twiddles;       // Twiddle factors
    int* bit_reverse_table;    // Bit reversal indices
    void* algorithm_data;      // Algorithm-specific data
    
    // GPU resources
    fft_gpu_plan_t gpu_plan;
    fft_gpu_memory_t gpu_in;
    fft_gpu_memory_t gpu_out;
};

// Global state
static unsigned g_hardware_caps = 0;
static int g_num_threads = 0;
static int g_hardware_detected = 0;

// Detect hardware capabilities
static void detect_hardware(void) {
    if (g_hardware_detected) return;
    
    g_hardware_caps = 0;
    
    #ifdef __x86_64__
        // Check CPU features
        unsigned int eax, ebx, ecx, edx;
        
        // Check SSE
        __asm__ __volatile__ ("cpuid" : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) : "a"(1));
        if (edx & (1 << 25)) g_hardware_caps |= FFT_HW_CPU_SSE;
        
        // Check AVX
        if (ecx & (1 << 28)) g_hardware_caps |= FFT_HW_CPU_AVX;
        
        // Check AVX2
        __asm__ __volatile__ ("cpuid" : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx) : "a"(7), "c"(0));
        if (ebx & (1 << 5)) g_hardware_caps |= FFT_HW_CPU_AVX2;
    #endif
    
    #ifdef __ARM_NEON
        g_hardware_caps |= FFT_HW_CPU_NEON;
    #endif
    
    // Check GPU availability
    #if defined(__APPLE__) || defined(__CUDACC__)
    if (fft_gpu_available()) {
        fft_gpu_backend_t backend = fft_gpu_get_backend();
        if (backend == FFT_GPU_CUDA) {
            g_hardware_caps |= FFT_HW_GPU_CUDA;
        } else if (backend == FFT_GPU_METAL) {
            g_hardware_caps |= FFT_HW_GPU_MPS;
        }
    }
    #endif
    
    g_hardware_detected = 1;
}

// Factor analysis for algorithm selection
static int count_factors(int n, int factor) {
    int count = 0;
    while (n % factor == 0) {
        count++;
        n /= factor;
    }
    return count;
}

static int is_prime(int n) {
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    
    for (int i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

static int is_highly_composite(int n) {
    int factors = 0;
    int temp = n;
    
    factors += count_factors(temp, 2);
    temp /= (1 << count_factors(temp, 2));
    
    factors += count_factors(temp, 3);
    temp /= pow(3, count_factors(temp, 3));
    
    factors += count_factors(temp, 5);
    temp /= pow(5, count_factors(temp, 5));
    
    factors += count_factors(temp, 7);
    temp /= pow(7, count_factors(temp, 7));
    
    return (temp == 1 && factors >= 3);
}

// Select best algorithm for given size
static int select_algorithm(int n, unsigned flags) {
    // GPU preference
    if (flags & FFT_PREFER_GPU) {
        if (g_hardware_caps & FFT_HW_GPU_CUDA) {
            return ALGO_GPU_CUDA;
        } else if (g_hardware_caps & FFT_HW_GPU_MPS) {
            return ALGO_GPU_MPS;
        }
    }
    
    // CPU algorithms
    if (is_power_of_two(n)) {
        if (n <= 64) {
            return ALGO_RADIX2_DIT;  // Simple for small sizes
        } else if (n <= 1024) {
            if (is_power_of_two(n) && (n & 3) == 0) {
                return ALGO_RADIX4;  // Power of 4
            }
            return ALGO_RADIX2_DIT;
        } else {
            // Large sizes - use split-radix for optimal operation count
            if (flags & FFT_CONSERVE_MEMORY) {
                return ALGO_RADIX2_DIT;  // In-place
            }
            return ALGO_SPLIT_RADIX;
        }
    } else {
        // Non-power-of-2
        if (is_prime(n)) {
            return ALGO_BLUESTEIN;
        } else if (is_highly_composite(n)) {
            return ALGO_MIXED_RADIX;
        } else {
            return ALGO_BLUESTEIN;
        }
    }
}

// Create optimized plan
fft_plan_t fft_plan_dft_1d(int n, complex_t* in, complex_t* out,
                           int sign, unsigned flags) {
    if (n <= 0 || !in || !out) return NULL;
    
    detect_hardware();
    
    fft_plan_t plan = calloc(1, sizeof(struct fft_plan));
    if (!plan) return NULL;
    
    plan->n = n;
    plan->in = in;
    plan->out = out;
    plan->dir = (sign < 0) ? FFT_FORWARD : FFT_INVERSE;
    plan->flags = flags;
    
    // Select algorithm
    plan->algorithm = select_algorithm(n, flags);
    
    // Precompute data based on algorithm
    switch (plan->algorithm) {
        case ALGO_RADIX2_DIT:
        case ALGO_RADIX2_DIF:
        case ALGO_RADIX4:
        case ALGO_SPLIT_RADIX:
            // Precompute twiddle factors
            plan->twiddles = allocate_complex_array(n);
            for (int k = 0; k < n; k++) {
                plan->twiddles[k] = cexp(-2.0 * PI * I * k / n);
            }
            
            // Precompute bit reversal
            if (is_power_of_two(n)) {
                plan->bit_reverse_table = malloc(n * sizeof(int));
                int log2n = log2_int(n);
                for (int i = 0; i < n; i++) {
                    plan->bit_reverse_table[i] = bit_reverse(i, log2n);
                }
            }
            break;
            
        case ALGO_BLUESTEIN:
        case ALGO_MIXED_RADIX:
            // These algorithms handle their own precomputation
            break;
            
        case ALGO_GPU_CUDA:
        case ALGO_GPU_MPS:
            #if defined(__APPLE__) || defined(__CUDACC__)
            // Initialize GPU resources
            fft_gpu_init(FFT_GPU_AUTO);
            plan->gpu_plan = fft_gpu_plan_1d(n, 1, plan->dir);
            plan->gpu_in = fft_gpu_alloc(n);
            plan->gpu_out = (in == out) ? plan->gpu_in : fft_gpu_alloc(n);
            #endif
            break;
    }
    
    // Measure performance if requested
    if (flags & FFT_MEASURE) {
        // TODO: Benchmark different algorithms and select best
    }
    
    return plan;
}

// Execute plan
void fft_execute(fft_plan_t plan) {
    if (!plan) return;
    
    // Copy input if not in-place
    if (plan->in != plan->out) {
        memcpy(plan->out, plan->in, plan->n * sizeof(complex_t));
    }
    
    // Execute based on algorithm
    switch (plan->algorithm) {
        case ALGO_RADIX2_DIT:
            radix2_dit_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_RADIX2_DIF:
            radix2_dif_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_RADIX4:
            radix4_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_SPLIT_RADIX:
            split_radix_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_BLUESTEIN:
            bluestein_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_MIXED_RADIX:
            mixed_radix_fft(plan->out, plan->n, plan->dir);
            break;
            
        case ALGO_GPU_CUDA:
        case ALGO_GPU_MPS:
            #if defined(__APPLE__) || defined(__CUDACC__)
            fft_gpu_copy_h2d(plan->gpu_in, plan->in, plan->n);
            fft_gpu_execute(plan->gpu_plan, plan->gpu_in, plan->gpu_out);
            fft_gpu_copy_d2h(plan->out, plan->gpu_out, plan->n);
            #endif
            break;
    }
}

// Execute with different arrays
void fft_execute_dft(fft_plan_t plan, complex_t* in, complex_t* out) {
    if (!plan || !in || !out) return;
    
    // Temporarily swap arrays
    complex_t* old_in = plan->in;
    complex_t* old_out = plan->out;
    
    plan->in = in;
    plan->out = out;
    
    fft_execute(plan);
    
    // Restore original arrays
    plan->in = old_in;
    plan->out = old_out;
}

// Destroy plan
void fft_destroy_plan(fft_plan_t plan) {
    if (!plan) return;
    
    free_complex_array(plan->twiddles);
    free(plan->bit_reverse_table);
    
    if (plan->algorithm == ALGO_GPU_CUDA || plan->algorithm == ALGO_GPU_MPS) {
        #if defined(__APPLE__) || defined(__CUDACC__)
        fft_gpu_destroy_plan(plan->gpu_plan);
        fft_gpu_free(plan->gpu_in);
        if (plan->gpu_out != plan->gpu_in) {
            fft_gpu_free(plan->gpu_out);
        }
        #endif
    }
    
    free(plan);
}

// Simple one-shot interface
int fft_auto(complex_t* in, complex_t* out, int n, int sign) {
    fft_plan_t plan = fft_plan_dft_1d(n, in, out, sign, FFT_ESTIMATE);
    if (!plan) return -1;
    
    fft_execute(plan);
    fft_destroy_plan(plan);
    
    return 0;
}

// Hardware capabilities
unsigned fft_get_hardware_capabilities(void) {
    detect_hardware();
    return g_hardware_caps;
}

// Thread control
void fft_plan_with_nthreads(int nthreads) {
    g_num_threads = nthreads;
    #ifdef _OPENMP
    if (nthreads > 0) {
        omp_set_num_threads(nthreads);
    }
    #endif
}

// Aligned memory allocation
complex_t* fft_alloc_complex(size_t n) {
    void* ptr;
    #ifdef _WIN32
        ptr = _aligned_malloc(n * sizeof(complex_t), 64);
    #else
        if (posix_memalign(&ptr, 64, n * sizeof(complex_t)) != 0) {
            return NULL;
        }
    #endif
    return (complex_t*)ptr;
}

double* fft_alloc_real(size_t n) {
    void* ptr;
    #ifdef _WIN32
        ptr = _aligned_malloc(n * sizeof(double), 64);
    #else
        if (posix_memalign(&ptr, 64, n * sizeof(double)) != 0) {
            return NULL;
        }
    #endif
    return (double*)ptr;
}

void fft_free(void* p) {
    if (!p) return;
    #ifdef _WIN32
        _aligned_free(p);
    #else
        free(p);
    #endif
}

// Version string
const char* fft_version(void) {
    return "FFT Library v2.0.0 - Automatic Algorithm Selection";
}

// Real FFT planning
fft_plan_t fft_plan_r2c_1d(int n, double* in, complex_t* out, unsigned flags) {
    // For now, convert to complex and use regular FFT
    // TODO: Implement optimized real FFT
    complex_t* temp = fft_alloc_complex(n);
    for (int i = 0; i < n; i++) {
        temp[i] = in[i];
    }
    
    fft_plan_t plan = fft_plan_dft_1d(n, temp, out, -1, flags | FFT_REAL_INPUT);
    fft_free(temp);
    
    return plan;
}

fft_plan_t fft_plan_c2r_1d(int n, complex_t* in, double* out, unsigned flags) {
    // TODO: Implement optimized inverse real FFT
    return NULL;
}

// 2D FFT planning
fft_plan_t fft_plan_dft_2d(int rows, int cols, complex_t* in,
                           complex_t* out, int sign, unsigned flags) {
    // TODO: Implement 2D FFT planning
    return NULL;
}

// Wisdom export/import
char* fft_export_wisdom_to_string(void) {
    // TODO: Implement wisdom export
    return strdup("# FFT Wisdom v2.0.0\n");
}

int fft_import_wisdom_from_string(const char* wisdom) {
    // TODO: Implement wisdom import
    return wisdom != NULL;
}

#ifndef LIB_BUILD
// Demo program
int main() {
    printf("FFT Auto - Automatic Algorithm Selection v2.0.0\n");
    printf("===============================================\n\n");
    
    // Test automatic selection for various sizes
    int test_sizes[] = {16, 64, 256, 1024, 97, 360, 1000};
    const char* size_types[] = {"tiny", "small", "medium", "large", 
                                "prime", "composite", "arbitrary"};
    
    printf("Testing automatic algorithm selection:\n");
    printf("Size\tType\t\tAlgorithm\tTime (ms)\n");
    printf("----\t----\t\t---------\t---------\n");
    
    for (int i = 0; i < 7; i++) {
        int n = test_sizes[i];
        complex_t* data = fft_alloc_complex(n);
        
        // Generate test data
        for (int j = 0; j < n; j++) {
            data[j] = sin(2 * PI * j / n) + I * cos(2 * PI * j / n);
        }
        
        // Time the automatic FFT
        clock_t start = clock();
        fft_auto(data, data, n, -1);
        clock_t end = clock();
        
        double time_ms = ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0;
        
        // Determine which algorithm was selected
        const char* algo_name = "Unknown";
        if (is_power_of_two(n)) {
            if (n <= 64) algo_name = "Radix-2 DIT";
            else if (n <= 1024) algo_name = (n % 4 == 0) ? "Radix-4" : "Radix-2 DIT";
            else algo_name = "Split-Radix";
        } else {
            if (is_prime(n)) algo_name = "Bluestein";
            else if (is_highly_composite(n)) algo_name = "Mixed-Radix";
            else algo_name = "Bluestein";
        }
        
        printf("%d\t%s\t\t%s\t%.3f\n", n, size_types[i], algo_name, time_ms);
        
        fft_free(data);
    }
    
    // Test planning API
    printf("\n\nTesting planning API for repeated transforms:\n");
    int n = 1024;
    complex_t* in = fft_alloc_complex(n);
    complex_t* out = fft_alloc_complex(n);
    
    // Create plan
    fft_plan_t plan = fft_plan_dft_1d(n, in, out, -1, FFT_MEASURE);
    
    // Execute multiple times
    clock_t start = clock();
    for (int i = 0; i < 1000; i++) {
        fft_execute(plan);
    }
    clock_t end = clock();
    
    double avg_time = ((double)(end - start) / CLOCKS_PER_SEC) * 1000.0 / 1000.0;
    printf("Average time per FFT (1000 iterations): %.3f ms\n", avg_time);
    
    fft_destroy_plan(plan);
    fft_free(in);
    fft_free(out);
    
    // Show hardware capabilities
    printf("\n\nDetected Hardware Capabilities:\n");
    unsigned caps = fft_get_hardware_capabilities();
    
    printf("CPU: ");
    if (caps & FFT_HW_CPU_SSE) printf("SSE ");
    if (caps & FFT_HW_CPU_AVX) printf("AVX ");
    if (caps & FFT_HW_CPU_AVX2) printf("AVX2 ");
    if (caps & FFT_HW_CPU_NEON) printf("NEON ");
    printf("\n");
    
    printf("GPU: ");
    if (caps & FFT_HW_GPU_CUDA) printf("CUDA ");
    if (caps & FFT_HW_GPU_MPS) printf("Metal ");
    if (!(caps & (FFT_HW_GPU_CUDA | FFT_HW_GPU_MPS))) printf("None");
    printf("\n");
    
    return 0;
}
#endif
