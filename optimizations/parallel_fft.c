#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <pthread.h>
#include <unistd.h>
#include <sys/sysctl.h>
#include <errno.h>

// macOS doesn't have pthread_barrier_t, so we'll implement our own
#ifdef __APPLE__
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;

typedef void* pthread_barrierattr_t;

int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count) {
    if(count == 0) {
        errno = EINVAL;
        return -1;
    }
    if(pthread_mutex_init(&barrier->mutex, 0) < 0) {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0) {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }
    barrier->tripCount = count;
    barrier->count = 0;
    return 0;
}

int pthread_barrier_destroy(pthread_barrier_t *barrier) {
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *barrier) {
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount) {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    } else {
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}
#endif

/**
 * Multi-threaded FFT Implementation
 * 
 * Parallelizes FFT computation using pthreads for multi-core processors.
 * Implements both data-parallel and task-parallel approaches.
 * 
 * Strategies:
 * 1. Parallel butterflies within stages
 * 2. Parallel independent FFTs
 * 3. Four-step FFT algorithm for cache efficiency
 * 
 * Benefits:
 * - Near-linear speedup with core count
 * - Better cache utilization
 * - Scalable to large FFT sizes
 */

// Thread configuration
typedef struct {
    int num_threads;
    int min_parallel_size;  // Minimum FFT size to parallelize
} thread_config_t;

// Thread work unit
typedef struct {
    complex_t* data;
    int start;
    int end;
    int stage;
    int n;
    fft_direction dir;
    pthread_barrier_t* barrier;
} thread_work_t;

// Global thread pool
typedef struct {
    pthread_t* threads;
    thread_work_t* work_units;
    int num_threads;
    volatile int shutdown;
} thread_pool_t;

// Parallel butterfly computation for one stage
void* parallel_butterfly_stage(void* arg) {
    thread_work_t* work = (thread_work_t*)arg;
    int n = work->n;
    int stage = work->stage;
    
    int m = 1 << stage;
    int half_m = m >> 1;
    
    // Process assigned range of butterflies
    for (int k = work->start; k < work->end; k += m) {
        complex_t w_m = twiddle_factor(1, m, work->dir);
        complex_t w = 1.0;
        
        for (int j = 0; j < half_m; j++) {
            int idx1 = k + j;
            int idx2 = idx1 + half_m;
            
            complex_t t = work->data[idx2] * w;
            work->data[idx2] = work->data[idx1] - t;
            work->data[idx1] = work->data[idx1] + t;
            
            w *= w_m;
        }
    }
    
    return NULL;
}

// Parallel radix-2 FFT
void fft_radix2_parallel(complex_t* x, int n, fft_direction dir, int num_threads) {
    CHECK_POWER_OF_TWO(n);
    
    // Sequential bit reversal (could be parallelized too)
    int log2n = log2_int(n);
    for (int i = 0; i < n; i++) {
        int j = bit_reverse(i, log2n);
        if (i < j) {
            complex_t temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }
    
    // Create threads
    pthread_t* threads = (pthread_t*)malloc((num_threads - 1) * sizeof(pthread_t));
    thread_work_t* work_units = (thread_work_t*)malloc(num_threads * sizeof(thread_work_t));
    
    // Process each stage
    for (int stage = 1; stage <= log2n; stage++) {
        int m = 1 << stage;
        
        // Divide work among threads
        int butterflies_per_thread = n / num_threads;
        butterflies_per_thread = (butterflies_per_thread / m) * m;  // Align to stage size
        
        if (butterflies_per_thread < m) {
            // Not enough work to parallelize this stage
            thread_work_t single_work = {
                .data = x, .start = 0, .end = n,
                .stage = stage, .n = n, .dir = dir
            };
            parallel_butterfly_stage(&single_work);
        } else {
            // Launch parallel threads
            for (int t = 0; t < num_threads; t++) {
                work_units[t].data = x;
                work_units[t].start = t * butterflies_per_thread;
                work_units[t].end = (t == num_threads - 1) ? n : (t + 1) * butterflies_per_thread;
                work_units[t].stage = stage;
                work_units[t].n = n;
                work_units[t].dir = dir;
                
                if (t < num_threads - 1) {
                    pthread_create(&threads[t], NULL, parallel_butterfly_stage, &work_units[t]);
                }
            }
            
            // Main thread does work too
            parallel_butterfly_stage(&work_units[num_threads - 1]);
            
            // Wait for threads
            for (int t = 0; t < num_threads - 1; t++) {
                pthread_join(threads[t], NULL);
            }
        }
    }
    
    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        int elements_per_thread = n / num_threads;
        
        for (int t = 0; t < num_threads - 1; t++) {
            int start = t * elements_per_thread;
            int end = (t + 1) * elements_per_thread;
            
            // Simple inline scaling work
            for (int i = start; i < end; i++) {
                x[i] /= n;
            }
        }
        
        // Main thread handles remainder
        for (int i = (num_threads - 1) * elements_per_thread; i < n; i++) {
            x[i] /= n;
        }
    }
    
    free(threads);
    free(work_units);
}

// Four-step FFT algorithm for cache efficiency
void four_step_fft(complex_t* x, int n, fft_direction dir, int num_threads) {
    // For n = n1 * n2, compute:
    // 1. n2 FFTs of size n1 (columns)
    // 2. Multiply by twiddle factors
    // 3. n1 FFTs of size n2 (rows)
    // 4. Transpose
    
    int n1 = 1;
    while (n1 * n1 < n) n1 <<= 1;
    int n2 = n / n1;
    
    printf("Four-step FFT: n=%d, n1=%d, n2=%d\n", n, n1, n2);
    
    // Step 1: Column FFTs (parallelizable)
    #pragma omp parallel for num_threads(num_threads)
    for (int j = 0; j < n2; j++) {
        complex_t* column = allocate_complex_array(n1);
        
        // Extract column
        for (int i = 0; i < n1; i++) {
            column[i] = x[i * n2 + j];
        }
        
        // FFT of column
        radix2_dit_fft(column, n1, dir);
        
        // Put back
        for (int i = 0; i < n1; i++) {
            x[i * n2 + j] = column[i];
        }
        
        free_complex_array(column);
    }
    
    // Step 2: Twiddle factor multiplication
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            if (i > 0 || j > 0) {
                x[i * n2 + j] *= twiddle_factor(i * j, n, dir);
            }
        }
    }
    
    // Step 3: Row FFTs
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < n1; i++) {
        radix2_dit_fft(&x[i * n2], n2, dir);
    }
    
    // Step 4: Transpose (could be optimized)
    complex_t* temp = allocate_complex_array(n);
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            temp[j * n1 + i] = x[i * n2 + j];
        }
    }
    memcpy(x, temp, n * sizeof(complex_t));
    free_complex_array(temp);
}

// Benchmark parallel vs sequential
void benchmark_parallel_fft() {
    printf("\nParallel FFT Benchmarks:\n");
    printf("========================\n");
    
    int sizes[] = {1024, 8192, 65536, 262144};
    int thread_counts[] = {1, 2, 4, 8};
    
    for (int s = 0; s < 4; s++) {
        int n = sizes[s];
        printf("\nSize: %d\n", n);
        printf("Threads\tTime (ms)\tSpeedup\n");
        printf("-------\t---------\t-------\n");
        
        complex_t* data = allocate_complex_array(n);
        complex_t* reference = allocate_complex_array(n);
        
        // Generate test data
        for (int i = 0; i < n; i++) {
            data[i] = ((double)rand() / RAND_MAX) + I * ((double)rand() / RAND_MAX);
            reference[i] = data[i];
        }
        
        fft_timer_t timer;
        double sequential_time = 0;
        
        for (int t = 0; t < 4; t++) {
            int num_threads = thread_counts[t];
            
            // Reset data
            memcpy(data, reference, n * sizeof(complex_t));
            
            timer_start(&timer);
            if (num_threads == 1) {
                radix2_dit_fft(data, n, FFT_FORWARD);
            } else {
                fft_radix2_parallel(data, n, FFT_FORWARD, num_threads);
            }
            timer_stop(&timer);
            
            if (num_threads == 1) {
                sequential_time = timer.elapsed_ms;
            }
            
            double speedup = sequential_time / timer.elapsed_ms;
            printf("%d\t%.3f\t\t%.2fx\n", num_threads, timer.elapsed_ms, speedup);
        }
        
        free_complex_array(data);
        free_complex_array(reference);
    }
}

// Thread pool for batch FFT processing
void* fft_worker_thread(void* arg) {
    thread_pool_t* pool = (thread_pool_t*)arg;
    
    while (!pool->shutdown) {
        // In a real implementation, would wait on condition variable
        // for work to be available
        usleep(1000);
    }
    
    return NULL;
}

// Initialize thread pool
thread_pool_t* create_thread_pool(int num_threads) {
    thread_pool_t* pool = (thread_pool_t*)malloc(sizeof(thread_pool_t));
    pool->num_threads = num_threads;
    pool->shutdown = 0;
    pool->threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&pool->threads[i], NULL, fft_worker_thread, pool);
    }
    
    return pool;
}

void destroy_thread_pool(thread_pool_t* pool) {
    pool->shutdown = 1;
    
    for (int i = 0; i < pool->num_threads; i++) {
        pthread_join(pool->threads[i], NULL);
    }
    
    free(pool->threads);
    free(pool);
}

// Main demonstration
int main() {
    printf("Multi-threaded FFT Implementation\n");
    printf("=================================\n");
    
    // Get number of CPU cores
#ifdef __APPLE__
    int num_cores;
    size_t size = sizeof(num_cores);
    sysctlbyname("hw.ncpu", &num_cores, &size, NULL, 0);
#else
    int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
#endif
    printf("\nSystem has %d CPU cores\n", num_cores);
    
    // Simple correctness test
    printf("\nCorrectness Test:\n");
    printf("-----------------\n");
    
    int n = 1024;
    complex_t* test_seq = allocate_complex_array(n);
    complex_t* test_par = allocate_complex_array(n);
    
    // Generate test signal
    for (int i = 0; i < n; i++) {
        double value = sin(2 * PI * 10 * i / n);
        test_seq[i] = test_par[i] = value;
    }
    
    // Sequential FFT
    radix2_dit_fft(test_seq, n, FFT_FORWARD);
    
    // Parallel FFT
    fft_radix2_parallel(test_par, n, FFT_FORWARD, 4);
    
    // Compare results
    double max_error = 0;
    for (int i = 0; i < n; i++) {
        double error = cabs(test_seq[i] - test_par[i]);
        if (error > max_error) max_error = error;
    }
    
    printf("Maximum error: %.2e %s\n", max_error,
           max_error < 1e-10 ? "(PASS)" : "(FAIL)");
    
    free_complex_array(test_seq);
    free_complex_array(test_par);
    
    // Run benchmarks
    benchmark_parallel_fft();
    
    // Parallel strategies
    printf("\n\nParallelization Strategies:\n");
    printf("===========================\n");
    printf("1. Data Parallel:\n");
    printf("   - Divide butterflies among threads\n");
    printf("   - Best for large FFTs\n");
    printf("   - Limited by memory bandwidth\n\n");
    
    printf("2. Task Parallel:\n");
    printf("   - Assign complete FFTs to threads\n");
    printf("   - Best for batch processing\n");
    printf("   - Good cache locality\n\n");
    
    printf("3. Hybrid:\n");
    printf("   - Four-step algorithm\n");
    printf("   - Recursive decomposition\n");
    printf("   - Balances computation and communication\n");
    
    // Best practices
    printf("\n\nBest Practices:\n");
    printf("===============\n");
    printf("1. Use power-of-2 thread counts\n");
    printf("2. Consider NUMA effects on large systems\n");
    printf("3. Minimize false sharing (pad data structures)\n");
    printf("4. Use thread pools for batch processing\n");
    printf("5. Profile to find optimal parallelization threshold\n");
    printf("6. Consider OpenMP for simpler code\n");
    printf("7. Be aware of Amdahl's Law limitations\n");
    
    return 0;
}
