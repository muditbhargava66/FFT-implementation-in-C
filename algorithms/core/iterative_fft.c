#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Iterative FFT Implementation
 * 
 * Memory-efficient iterative implementation of the FFT algorithm.
 * Shows explicit loop structure and in-place computation.
 * 
 * Benefits over recursive:
 * - O(1) space complexity (in-place)
 * - Better cache locality
 * - No function call overhead
 * - Predictable memory access patterns
 * 
 * Time Complexity: O(n log n)
 * Space Complexity: O(1)
 */

// Iterative bit-reversal with optimization
void optimized_bit_reversal(complex_t* x, int n) {
    int log2n = log2_int(n);
    
    // Pre-compute bit-reversal table for small sizes
    if (n <= 256) {
        int* rev_table = (int*)malloc(n * sizeof(int));
        
        for (int i = 0; i < n; i++) {
            rev_table[i] = bit_reverse(i, log2n);
        }
        
        // Apply permutation
        for (int i = 0; i < n; i++) {
            int j = rev_table[i];
            if (i < j) {
                complex_t temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
        
        free(rev_table);
    } else {
        // For large sizes, compute on the fly
        for (int i = 0; i < n; i++) {
            int j = bit_reverse(i, log2n);
            if (i < j) {
                complex_t temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }
}

// Iterative FFT with detailed loop structure
void iterative_fft_detailed(complex_t* x, int n, fft_direction dir) {
    CHECK_POWER_OF_TWO(n);
    
    // Step 1: Bit-reversal permutation
    optimized_bit_reversal(x, n);
    
    // Step 2: Iterative computation
    int log2n = log2_int(n);
    
    // For each stage of the FFT
    for (int stage = 1; stage <= log2n; stage++) {
        int m = 1 << stage;        // Current FFT size (2, 4, 8, ...)
        int half_m = m >> 1;        // Half of current size
        
        // Twiddle factor for this stage
        complex_t w_m = twiddle_factor(1, m, dir);
        
        // Process each FFT of size m
        for (int k = 0; k < n; k += m) {
            complex_t w = 1.0;
            
            // Butterfly operations within this FFT
            for (int j = 0; j < half_m; j++) {
                int idx1 = k + j;
                int idx2 = idx1 + half_m;
                
                complex_t t = x[idx2] * w;
                x[idx2] = x[idx1] - t;
                x[idx1] = x[idx1] + t;
                
                w *= w_m;
            }
        }
    }
    
    // Step 3: Scale for inverse transform
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Memory access pattern visualization
void visualize_memory_access(int n) {
    printf("\nMemory Access Pattern for n=%d:\n", n);
    printf("================================\n");
    
    int log2n = log2_int(n);
    
    // Show bit-reversal access pattern
    printf("\nBit-reversal stage:\n");
    for (int i = 0; i < n && i < 16; i++) {
        int j = bit_reverse(i, log2n);
        printf("  x[%2d] <-> x[%2d]\n", i, j);
    }
    if (n > 16) printf("  ...\n");
    
    // Show butterfly access patterns for each stage
    for (int stage = 1; stage <= log2n && stage <= 3; stage++) {
        int m = 1 << stage;
        int half_m = m >> 1;
        
        printf("\nStage %d (m=%d):\n", stage, m);
        int count = 0;
        
        for (int k = 0; k < n && count < 8; k += m) {
            for (int j = 0; j < half_m && count < 8; j++) {
                int idx1 = k + j;
                int idx2 = idx1 + half_m;
                printf("  Butterfly: x[%2d] <-> x[%2d]\n", idx1, idx2);
                count++;
            }
        }
        if (count == 8 && n > 16) printf("  ...\n");
    }
}

// Cache analysis
typedef struct {
    int cache_misses;
    int cache_hits;
    int* last_access;
    int time;
    int cache_line_size;
} cache_sim_t;

void simulate_cache_behavior(int n, int cache_size) {
    cache_sim_t sim;
    sim.cache_misses = 0;
    sim.cache_hits = 0;
    sim.time = 0;
    sim.cache_line_size = 64 / sizeof(complex_t);  // 64-byte cache lines
    sim.last_access = (int*)calloc(n, sizeof(int));
    
    printf("\nCache Simulation (cache size = %d elements):\n", cache_size);
    printf("============================================\n");
    
    // Simulate bit-reversal
    for (int i = 0; i < n; i++) {
        int j = bit_reverse(i, log2_int(n));
        
        // Check cache for i and j
        if (sim.time - sim.last_access[i] > cache_size) sim.cache_misses++;
        else sim.cache_hits++;
        
        if (sim.time - sim.last_access[j] > cache_size) sim.cache_misses++;
        else sim.cache_hits++;
        
        sim.last_access[i] = sim.time++;
        sim.last_access[j] = sim.time++;
    }
    
    printf("Bit-reversal: %d hits, %d misses (%.1f%% hit rate)\n",
           sim.cache_hits, sim.cache_misses,
           100.0 * sim.cache_hits / (sim.cache_hits + sim.cache_misses));
    
    free(sim.last_access);
}

// Performance analysis function
void analyze_iterative_performance() {
    printf("\nIterative FFT Performance Analysis:\n");
    printf("===================================\n");
    
    // Operation count analysis
    printf("\nOperation counts per stage:\n");
    printf("Stage\tSize\tButterflies\tComplex Muls\n");
    printf("-----\t----\t-----------\t------------\n");
    
    int n = 64;
    int log2n = log2_int(n);
    int total_muls = 0;
    
    for (int stage = 1; stage <= log2n; stage++) {
        int m = 1 << stage;
        int butterflies = n / 2;
        int muls_per_stage = butterflies;
        total_muls += muls_per_stage;
        
        printf("%d\t%d\t%d\t\t%d\n", stage, m, butterflies, muls_per_stage);
    }
    
    printf("\nTotal: %d complex multiplications\n", total_muls);
    printf("Formula: (n/2) * log2(n) = %d\n", (n/2) * log2n);
}

// Main demonstration
int main() {
    printf("Iterative FFT Implementation\n");
    printf("============================\n");
    
    // Demonstrate with small example
    int n = 8;
    complex_t* signal = allocate_complex_array(n);
    
    printf("\nStep-by-step FFT computation for n=%d:\n", n);
    printf("-------------------------------------\n");
    
    // Initialize with simple sequence
    printf("Input: ");
    for (int i = 0; i < n; i++) {
        signal[i] = i;
        printf("%.0f ", creal(signal[i]));
    }
    printf("\n");
    
    // Apply FFT with detailed output
    iterative_fft_detailed(signal, n, FFT_FORWARD);
    
    printf("\nOutput: ");
    for (int i = 0; i < n; i++) {
        print_complex(signal[i]);
        printf(" ");
    }
    printf("\n");
    
    // Visualize memory access patterns
    visualize_memory_access(16);
    
    // Cache behavior analysis
    simulate_cache_behavior(64, 16);
    
    // Performance analysis
    analyze_iterative_performance();
    
    // Compare with recursive implementation
    printf("\n\nIterative vs Recursive Comparison:\n");
    printf("==================================\n");
    printf("Aspect\t\t\tIterative\tRecursive\n");
    printf("------\t\t\t---------\t---------\n");
    printf("Space complexity\tO(1)\t\tO(n log n)\n");
    printf("Cache efficiency\tBetter\t\tWorse\n");
    printf("Function overhead\tNone\t\tO(n log n) calls\n");
    printf("Code complexity\t\tModerate\tSimple\n");
    printf("Parallelization\t\tEasier\t\tHarder\n");
    
    // Large-scale performance test
    printf("\n\nLarge-scale Performance Test:\n");
    printf("=============================\n");
    printf("Size\tTime (ms)\tThroughput (M samples/s)\n");
    printf("----\t---------\t------------------------\n");
    
    for (int size = 1024; size <= 65536; size *= 2) {
        complex_t* x = allocate_complex_array(size);
        
        // Fill with random data
        for (int i = 0; i < size; i++) {
            x[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        
        timer_t timer;
        timer_start(&timer);
        
        // Run multiple iterations for accuracy
        int iterations = 1000000 / size;
        if (iterations < 1) iterations = 1;
        
        for (int iter = 0; iter < iterations; iter++) {
            iterative_fft_detailed(x, size, FFT_FORWARD);
        }
        
        timer_stop(&timer);
        
        double time_per_fft = timer.elapsed_ms / iterations;
        double throughput = size / (time_per_fft * 1000.0);  // M samples/s
        
        printf("%d\t%.3f\t\t%.2f\n", size, time_per_fft, throughput);
        
        free_complex_array(x);
    }
    
    free_complex_array(signal);
    
    return 0;
}
