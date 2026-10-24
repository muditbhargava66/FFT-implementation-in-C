#include "../../include/fft_auto.h"
#include "../../include/fft_algorithms.h"
#include "../../include/fft_gpu.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file fft_auto.c
 * @brief Planner, algorithm selection, wisdom cache, real and 2-D paths.
 */

/* Algorithm identifiers used both in plans and in the wisdom cache */
enum fft_algo {
  ALGO_RADIX2_DIT,
  ALGO_RADIX2_DIF,
  ALGO_RADIX4,
  ALGO_SPLIT_RADIX,
  ALGO_BLUESTEIN,
  ALGO_MIXED_RADIX,
  ALGO_GPU_CUDA,
  ALGO_GPU_MPS,
  ALGO_COUNT
};

/* Defined below with the execution paths; the measurement in choose_algorithm()
 * needs to run transforms while planning. */
static fft_error_t execute_1d(complex_t *buf, int len, fft_direction dir,
                              enum fft_algo algo);

/* Human-readable name for an algorithm */
static const char *algorithm_name(enum fft_algo a) {
  static const char *names[] = {"Radix-2 DIT", "Radix-2 DIF", "Radix-4",
                                "Split-Radix", "Bluestein",   "Mixed-Radix",
                                "GPU (CUDA)",  "GPU (Metal)"};
  return (a >= 0 && a < ALGO_COUNT) ? names[a] : "Unknown";
}

/* Internal plan structure */
struct fft_plan {
  int n;             /* Transform size (1-D), or rows*cols for 2-D */
  complex_t *in;     /* Input buffer */
  complex_t *out;    /* Output buffer */
  fft_direction dir; /* Transform direction */
  unsigned flags;    /* Planning flags */
  enum fft_algo algorithm;
  enum fft_algo col_algorithm; /* 2-D column pass only */

  /* Optional real-data pointers (r2c / c2r) */
  double *real_in;   /* non-NULL for r2c plans */
  double *real_out;  /* non-NULL for c2r plans */

  /* 2-D shape (zero for 1-D plans) */
  int rows;
  int cols;

  /* Scratch buffer for r2c / c2r plans */
  void *algorithm_data;

  /* GPU resources */
  fft_gpu_plan_t gpu_plan;
  fft_gpu_memory_t gpu_in;
  fft_gpu_memory_t gpu_out;
};

// Global state
static unsigned g_hardware_caps = 0;
static int g_hardware_detected = 0;

// Detect hardware capabilities
static void detect_hardware(void) {
  if (g_hardware_detected)
    return;

  g_hardware_caps = 0;

#ifdef __x86_64__
  // Check CPU features
  unsigned int eax, ebx, ecx, edx;

  // Highest leaf this CPU answers. Leaf 7 does not exist on older parts, and
  // reading it returns undefined data rather than zeroes.
  __asm__ __volatile__("cpuid"
                       : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                       : "a"(0));
  unsigned int max_leaf = eax;

  // Leaf 1: EDX bit 25 = SSE, ECX bit 27 = OSXSAVE, ECX bit 28 = AVX
  __asm__ __volatile__("cpuid"
                       : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                       : "a"(1));

  /* CPUID says the silicon can run the instruction; only the operating system
   * can say whether it preserves the vector register state across context
   * switches. Without XCR0 bits 1 (SSE) and 2 (AVX) set, executing SSE/AVX
   * raises #UD, so advertising support here would hand callers a crash. */
  unsigned int xcr0 = 0;
  int osxsave = (ecx & (1u << 27)) != 0;
  if (osxsave) {
    unsigned int hi;
    __asm__ __volatile__("xgetbv" : "=a"(xcr0), "=d"(hi) : "c"(0));
    (void)hi;
  }

  if ((edx & (1u << 25)) && (!osxsave || (xcr0 & 0x2)))
    g_hardware_caps |= FFT_HW_CPU_SSE;

  int avx_ok = (ecx & (1u << 28)) && osxsave && (xcr0 & 0x6) == 0x6;
  if (avx_ok)
    g_hardware_caps |= FFT_HW_CPU_AVX;

  // Check AVX2
  if (avx_ok && max_leaf >= 7) {
    __asm__ __volatile__("cpuid"
                         : "=a"(eax), "=b"(ebx), "=c"(ecx), "=d"(edx)
                         : "a"(7), "c"(0));
    if (ebx & (1 << 5))
      g_hardware_caps |= FFT_HW_CPU_AVX2;
  }
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

/* Factor analysis for algorithm selection */
static int is_highly_composite(int n) {
  /* A number is deemed "highly composite" for our purposes if it
     factors completely into small primes (2, 3, 5, 7) and has >= 3
     total prime factors. */
  static const int small_primes[] = {2, 3, 5, 7};
  int total_factors = 0;
  for (int p = 0; p < 4; p++) {
    while (n % small_primes[p] == 0) {
      total_factors++;
      n /= small_primes[p];
    }
  }
  return (n == 1 && total_factors >= 3);
}

/*
 * Whether `algo` can execute a transform of length n. This mirrors the size
 * checks inside the kernels so that a wisdom entry can be validated before it
 * becomes a plan. Size validity is machine-independent; hardware availability
 * is a separate question, handled below, so that an exported wisdom file stays
 * portable between machines.
 */
static int algorithm_handles_size(enum fft_algo algo, int n) {
  switch (algo) {
  case ALGO_RADIX2_DIT:
  case ALGO_RADIX2_DIF:
  case ALGO_RADIX4:
  case ALGO_SPLIT_RADIX:
    return is_power_of_two(n);
  case ALGO_BLUESTEIN:
    /* Same bound as bluestein_fft(): 2n-1 must not overflow. */
    return n > 0 && n <= (1 << 29);
  case ALGO_MIXED_RADIX:
    return n > 0;
  case ALGO_GPU_CUDA:
  case ALGO_GPU_MPS:
    return n > 0; /* cuFFT and MPS accept non-power-of-two lengths */
  default:
    return 0; /* ALGO_COUNT, or an id this build does not know */
  }
}

/* Whether this build and machine can dispatch `algo` at all. */
static int algorithm_available(enum fft_algo algo) {
  switch (algo) {
  case ALGO_GPU_CUDA:
    return (g_hardware_caps & FFT_HW_GPU_CUDA) != 0;
  case ALGO_GPU_MPS:
    return (g_hardware_caps & FFT_HW_GPU_MPS) != 0;
  case ALGO_COUNT:
    return 0;
  default:
    return 1;
  }
}

/* Select best algorithm for given size */
static enum fft_algo select_algorithm(int n, unsigned flags) {
  /* GPU preference */
  if (flags & FFT_PREFER_GPU) {
    if (g_hardware_caps & FFT_HW_GPU_CUDA)
      return ALGO_GPU_CUDA;
    if (g_hardware_caps & FFT_HW_GPU_MPS)
      return ALGO_GPU_MPS;
  }

  /* The size rules come from the ratios in docs/performance.md: radix-4 loses
   * to radix-2 DIT at 16 and 64, ties at 256, and wins from 512 up, so 256 is
   * the crossover. Bluestein won every non-power-of-two length measured, but
   * needs about 9n complex numbers of scratch against mixed-radix's n. */
  if (is_power_of_two(n)) {
    if (flags & FFT_CONSERVE_MEMORY)
      return ALGO_RADIX2_DIT;
    return (n <= 256) ? ALGO_RADIX2_DIT : ALGO_RADIX4;
  }

  if ((flags & FFT_CONSERVE_MEMORY) && is_highly_composite(n))
    return ALGO_MIXED_RADIX;
  return ALGO_BLUESTEIN;
}

/* ---- Wisdom cache: (n, flags) -> algorithm, filled by FFT_MEASURE ---- */
#define WISDOM_MAX 256

typedef struct {
  int n;
  unsigned flags;
  enum fft_algo algo;
} wisdom_entry_t;

static wisdom_entry_t g_wisdom[WISDOM_MAX];
static int g_wisdom_count = 0;

static int wisdom_lookup(int n, unsigned flags, enum fft_algo *out) {
  for (int i = 0; i < g_wisdom_count; i++) {
    if (g_wisdom[i].n == n && g_wisdom[i].flags == flags) {
      *out = g_wisdom[i].algo;
      return 1;
    }
  }
  return 0;
}

static void wisdom_store(int n, unsigned flags, enum fft_algo algo) {
  for (int i = 0; i < g_wisdom_count; i++) {
    if (g_wisdom[i].n == n && g_wisdom[i].flags == flags) {
      g_wisdom[i].algo = algo;
      return;
    }
  }
  if (g_wisdom_count < WISDOM_MAX) {
    g_wisdom[g_wisdom_count].n = n;
    g_wisdom[g_wisdom_count].flags = flags;
    g_wisdom[g_wisdom_count].algo = algo;
    g_wisdom_count++;
  }
  /* A full cache stops recording; the heuristic still answers every plan. */
}

/*
 * Return the per-transform time of one batch in milliseconds, growing the batch
 * until clock() can resolve it. On Windows the clock advances about once a
 * millisecond, so a single small transform would measure as 0.000 ms and leave
 * the ranking to list order. *reps persists so later rounds reuse the size.
 */
static double time_batch(enum fft_algo algo, complex_t *scratch,
                         const complex_t *seed, int n, int *reps) {
  /* 20 ms in clock ticks. A literal tick count is not portable: CLOCKS_PER_SEC
   * is 1000 on Windows and 1000000 on POSIX. */
  const double target_ticks = 0.020 * (double)CLOCKS_PER_SEC;
  const int max_reps = 2000000;

  for (;;) {
    clock_t start = clock();
    for (int r = 0; r < *reps; r++) {
      /* A forward transform scales by n, so replaying one buffer climbs to
       * Inf/NaN and times the wrong work. Restore the seed each rep. */
      memcpy(scratch, seed, (size_t)n * sizeof(complex_t));
      execute_1d(scratch, n, FFT_FORWARD, algo);
    }
    double ticks = (double)(clock() - start);
    double ms = ticks * 1000.0 / ((double)CLOCKS_PER_SEC * (double)*reps);

    if (ticks >= target_ticks || *reps >= max_reps)
      return ms;

    /* Scale towards the target, or by 8x while the clock reports no elapsed
     * time at all. */
    double want = (ticks > 0.0) ? (double)(*reps) * (target_ticks / ticks) * 1.25
                                : (double)(*reps) * 8.0;
    int next = (want >= (double)max_reps) ? max_reps : (int)want;
    if (next <= *reps)
      next = *reps + 1;
    *reps = next;
  }
}

/*
 * Compare a challenger and the incumbent in alternating batches. The power-of-
 * two kernels are within a few percent of each other, and clock drift over a
 * one-pass sweep is larger than that. Alternating puts both at the same
 * moments. Returns 1 if the challenger wins all but one round and is faster by
 * more than `margin` overall.
 */
static int confirm_choice(enum fft_algo challenger, enum fft_algo incumbent,
                          complex_t *scratch, const complex_t *seed, int n) {
  const int rounds = 5;
  const double margin = 0.95;

  int rep_a = 1, rep_b = 1;
  double sum_a = 0.0, sum_b = 0.0;
  int wins = 0;

  for (int r = 0; r < rounds; r++) {
    double a = time_batch(challenger, scratch, seed, n, &rep_a);
    double b = time_batch(incumbent, scratch, seed, n, &rep_b);
    if (a < b)
      wins++;
    sum_a += a;
    sum_b += b;
  }

  return wins >= rounds - 1 && sum_a < sum_b * margin;
}

/*
 * Time the CPU algorithms that accept length n. Each is sampled three times
 * with the minimum kept, and the samples are interleaved so clock drift moves
 * all candidates together. `default_algo` is the heuristic's answer; a
 * challenger must also survive confirm_choice() before it is used.
 */
static enum fft_algo measure_best_algorithm(int n, enum fft_algo default_algo) {
  static const enum fft_algo cpu_algorithms[] = {
      ALGO_RADIX2_DIT, ALGO_RADIX2_DIF, ALGO_RADIX4,
      ALGO_SPLIT_RADIX, ALGO_BLUESTEIN, ALGO_MIXED_RADIX};
  const int num_algorithms =
      (int)(sizeof(cpu_algorithms) / sizeof(cpu_algorithms[0]));
  const int trials = 3;

  complex_t *scratch = allocate_complex_array(n);
  complex_t *seed = allocate_complex_array(n);
  if (!scratch || !seed) {
    free_complex_array(scratch);
    free_complex_array(seed);
    return ALGO_COUNT;
  }

  /* Deterministic data, so planning twice reaches the same conclusion. Held in
   * its own buffer because every timed repetition overwrites scratch. */
  for (int i = 0; i < n; i++)
    seed[i] = (double)((i * 31) % 97) - I * (double)((i * 17) % 89);

  int reps[sizeof(cpu_algorithms) / sizeof(cpu_algorithms[0])];
  double best_ms[sizeof(cpu_algorithms) / sizeof(cpu_algorithms[0])];

  for (int c = 0; c < num_algorithms; c++) {
    enum fft_algo algo = cpu_algorithms[c];

    /* Also the validity check: radix-4 and the radix-2 family reject lengths
     * they cannot factor, so non-power-of-two sizes drop out here.
     *
     * mixed_radix_fft() forwards power-of-two lengths straight to
     * radix2_dit_fft(), which would make the sweep measure the same code twice
     * and let noise name a winner between them. */
    memcpy(scratch, seed, (size_t)n * sizeof(complex_t));
    if (execute_1d(scratch, n, FFT_FORWARD, algo) != FFT_OK ||
        (algo == ALGO_MIXED_RADIX && is_power_of_two(n))) {
      best_ms[c] = -1.0;
      continue;
    }
    reps[c] = 1;
    best_ms[c] = time_batch(algo, scratch, seed, n, &reps[c]);
  }

  for (int t = 1; t < trials; t++) {
    for (int c = 0; c < num_algorithms; c++) {
      if (best_ms[c] < 0.0)
        continue;
      double ms = time_batch(cpu_algorithms[c], scratch, seed, n, &reps[c]);
      if (ms < best_ms[c])
        best_ms[c] = ms;
    }
  }

  enum fft_algo winner = ALGO_COUNT;
  double winner_ms = -1.0;
  double default_ms = -1.0;

  for (int c = 0; c < num_algorithms; c++) {
    if (best_ms[c] < 0.0)
      continue;
    if (cpu_algorithms[c] == default_algo)
      default_ms = best_ms[c];
    if (winner_ms < 0.0 || best_ms[c] < winner_ms) {
      winner_ms = best_ms[c];
      winner = cpu_algorithms[c];
    }
  }

  enum fft_algo chosen = winner;
  if (chosen != default_algo) {
    /* Trust the challenge only when the sweep saw a >15% win *and* a paired
     * re-run repeats it. Repeated runs of one kernel move by about 10%, so
     * smaller gaps are noise, and caching noise would make the choice depend on
     * which process happened to be running. */
    if (default_ms <= 0.0 || winner_ms >= default_ms * 0.85 ||
        !confirm_choice(winner, default_algo, scratch, seed, n))
      chosen = default_algo;
  }

  free_complex_array(scratch);
  free_complex_array(seed);
  return chosen;
}

/*
 * Cached choice if there is one, otherwise the size heuristic, optionally
 * upgraded by FFT_MEASURE. Returns ALGO_COUNT when FFT_WISDOM_ONLY was
 * requested and no usable entry exists.
 */
static enum fft_algo choose_algorithm(int n, unsigned flags) {
  /* Mask FFT_WISDOM_ONLY out of the key so a decision cached under
   * FFT_MEASURE can be found by FFT_MEASURE | FFT_WISDOM_ONLY. */
  const unsigned key = flags & ~(unsigned)FFT_WISDOM_ONLY;

  enum fft_algo cached;
  /* Imported entries are untrusted: check that the algorithm can handle this
   * size and is available here. A bad entry is a cache miss. */
  if (wisdom_lookup(n, key, &cached) && algorithm_handles_size(cached, n) &&
      algorithm_available(cached))
    return cached;

  if (flags & FFT_WISDOM_ONLY)
    return ALGO_COUNT;

  enum fft_algo algo = select_algorithm(n, flags);

  /* Memory is a constraint, not a tie-breaker for measurement. */
  if ((flags & FFT_MEASURE) && !(flags & FFT_CONSERVE_MEMORY) &&
      algo != ALGO_GPU_CUDA && algo != ALGO_GPU_MPS) {
    enum fft_algo timed = measure_best_algorithm(n, algo);
    if (timed != ALGO_COUNT && timed != ALGO_GPU_CUDA && timed != ALGO_GPU_MPS)
      algo = timed;
  }

  wisdom_store(n, key, algo);
  return algo;
}

// Create optimized plan
fft_plan_t fft_plan_dft_1d(int n, complex_t *in, complex_t *out, int sign,
                           unsigned flags) {
  if (n <= 0 || !in || !out)
    return NULL;

  detect_hardware();

  fft_plan_t plan = calloc(1, sizeof(struct fft_plan));
  if (!plan)
    return NULL;

  plan->n = n;
  plan->in = in;
  plan->out = out;
  plan->dir = (sign < 0) ? FFT_FORWARD : FFT_INVERSE;
  plan->flags = flags;

  // Select algorithm, consulting and refilling the wisdom cache.
  plan->algorithm = choose_algorithm(n, flags);

  /* ALGO_COUNT means FFT_WISDOM_ONLY found no entry. Return NULL rather than a
   * plan that fails on every execute. */
  if (plan->algorithm == ALGO_COUNT) {
    free(plan);
    return NULL;
  }

  // Allocate per-algorithm resources.
  switch (plan->algorithm) {
  case ALGO_RADIX2_DIT:
  case ALGO_RADIX2_DIF:
  case ALGO_RADIX4:
  case ALGO_SPLIT_RADIX:
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

  case ALGO_COUNT:
    break;
  }

  return plan;
}

/* Dispatch a 1-D FFT on buffer buf of length len */
static fft_error_t execute_1d(complex_t *buf, int len, fft_direction dir,
                       enum fft_algo algo) {
  switch (algo) {
  case ALGO_RADIX2_DIT: return radix2_dit_fft(buf, len, dir);
  case ALGO_RADIX2_DIF: return radix2_dif_fft(buf, len, dir);
  case ALGO_RADIX4:     return radix4_fft(buf, len, dir);
  case ALGO_SPLIT_RADIX:return split_radix_fft(buf, len, dir);
  case ALGO_BLUESTEIN:  return bluestein_fft(buf, len, dir);
  case ALGO_MIXED_RADIX:return mixed_radix_fft(buf, len, dir);
  case ALGO_COUNT:      return FFT_ERR_SIZE; /* sentinel — never selected */
  default:              return FFT_ERR_SIZE; /* unknown algorithm id */
  }
}

/* Execute plan */
fft_error_t fft_execute(fft_plan_t plan) {
  if (!plan)
    return FFT_ERR_NULL;

  /* ---------- r2c: promote real input, execute, copy n/2+1 bins ---------- */
  if (plan->real_in) {
    complex_t *cbuf = (complex_t *)plan->algorithm_data;
    if (!cbuf) return FFT_ERR_ALLOC;
    for (int i = 0; i < plan->n; i++)
      cbuf[i] = plan->real_in[i];
    fft_error_t err = execute_1d(cbuf, plan->n, plan->dir, plan->algorithm);
    if (err) return err;
    int half = plan->n / 2 + 1;
    memcpy(plan->out, cbuf, (size_t)half * sizeof(complex_t));
    return FFT_OK;
  }

  /* ---------- c2r: expand Hermitian spectrum, inverse, extract real ---------- */
  if (plan->real_out) {
    complex_t *cbuf = (complex_t *)plan->algorithm_data;
    if (!cbuf) return FFT_ERR_ALLOC;
    /* in[] holds n/2+1 bins per fft_plan_c2r_1d's contract; the remaining bins
     * are the conjugate mirror, not extra caller data. */
    int half = plan->n / 2 + 1;
    memcpy(cbuf, plan->in, (size_t)half * sizeof(complex_t));
    for (int k = half; k < plan->n; k++)
      cbuf[k] = conj(plan->in[plan->n - k]);
    fft_error_t err = execute_1d(cbuf, plan->n, plan->dir, plan->algorithm);
    if (err) return err;
    for (int i = 0; i < plan->n; i++)
      plan->real_out[i] = creal(cbuf[i]);
    return FFT_OK;
  }

  /* ---------- 2-D: row-column decomposition ---------- */
  if (plan->rows > 0 && plan->cols > 0) {
    int rows = plan->rows;
    int cols = plan->cols;
    complex_t *buf = plan->out;

    /* Allocate the transpose buffer before touching either data buffer, so a
     * failed allocation leaves the caller's data alone. */
    complex_t *tmp = allocate_complex_array(rows * cols);
    if (!tmp) return FFT_ERR_ALLOC;

    if (plan->in != plan->out)
      memcpy(buf, plan->in, (size_t)rows * cols * sizeof(complex_t));

    /* Row-wise FFTs */
    for (int r = 0; r < rows; r++) {
      fft_error_t err = execute_1d(&buf[r * cols], cols, plan->dir, plan->algorithm);
      if (err) {
        free_complex_array(tmp);
        return err;
      }
    }

    /* Transpose into a temp buffer */
    for (int r = 0; r < rows; r++)
      for (int c = 0; c < cols; c++)
        tmp[c * rows + r] = buf[r * cols + c];

    /* Column-wise FFTs (now contiguous after transpose). The column length
     * was planned separately and stored, so execution does not consult the
     * wisdom cache again. */
    for (int c = 0; c < cols; c++) {
      fft_error_t err =
          execute_1d(&tmp[c * rows], rows, plan->dir, plan->col_algorithm);
      if (err) {
        free_complex_array(tmp);
        return err;
      }
    }

    /* Transpose back */
    for (int r = 0; r < rows; r++)
      for (int c = 0; c < cols; c++)
        buf[r * cols + c] = tmp[c * rows + r];

    free_complex_array(tmp);
    return FFT_OK;
  }

  /* ---------- Regular 1-D ---------- */
  if (plan->in != plan->out)
    memcpy(plan->out, plan->in, (size_t)plan->n * sizeof(complex_t));

  switch (plan->algorithm) {
  case ALGO_GPU_CUDA:
  case ALGO_GPU_MPS:
#if defined(__APPLE__) || defined(__CUDACC__)
    fft_gpu_copy_h2d(plan->gpu_in, plan->in, plan->n);
    fft_gpu_execute(plan->gpu_plan, plan->gpu_in, plan->gpu_out);
    fft_gpu_copy_d2h(plan->out, plan->gpu_out, plan->n);
#endif
    return FFT_OK;
  default:
    return execute_1d(plan->out, plan->n, plan->dir, plan->algorithm);
  }
}

// Execute with different arrays
fft_error_t fft_execute_dft(fft_plan_t plan, complex_t *in, complex_t *out) {
  if (!plan || !in || !out)
    return FFT_ERR_NULL;

  /* Real-input/output plans bind a double* at planning time and have no
   * complex in/out pair to substitute. Reject the call instead of quietly
   * transforming the buffers the plan was created with. */
  if (plan->real_in || plan->real_out)
    return FFT_ERR_UNSUPPORTED;

  // Temporarily swap arrays
  complex_t *old_in = plan->in;
  complex_t *old_out = plan->out;

  plan->in = in;
  plan->out = out;

  fft_error_t err = fft_execute(plan);

  // Restore original arrays
  plan->in = old_in;
  plan->out = old_out;

  return err;
}

/* Destroy plan */
void fft_destroy_plan(fft_plan_t plan) {
  if (!plan)
    return;

  /* Free internal complex buffer used by r2c / c2r plans */
  if (plan->algorithm_data)
    fft_free(plan->algorithm_data);

  if (plan->algorithm == ALGO_GPU_CUDA || plan->algorithm == ALGO_GPU_MPS) {
#if defined(__APPLE__) || defined(__CUDACC__)
    fft_gpu_destroy_plan(plan->gpu_plan);
    fft_gpu_free(plan->gpu_in);
    if (plan->gpu_out != plan->gpu_in)
      fft_gpu_free(plan->gpu_out);
#endif
  }

  free(plan);
}

// Simple one-shot interface
fft_error_t fft_auto(complex_t *in, complex_t *out, int n, int sign) {
  if (!in || !out) return FFT_ERR_NULL;
  if (n <= 0) return FFT_ERR_SIZE;

  fft_plan_t plan = fft_plan_dft_1d(n, in, out, sign, FFT_ESTIMATE);
  if (!plan) return FFT_ERR_ALLOC;

  fft_error_t err = fft_execute(plan);
  fft_destroy_plan(plan);

  return err;
}

const char *fft_plan_algorithm_name(fft_plan_t plan) {
  return plan ? algorithm_name(plan->algorithm) : "Unknown";
}

// Hardware capabilities
unsigned fft_get_hardware_capabilities(void) {
  detect_hardware();
  return g_hardware_caps;
}

// Aligned memory allocation
complex_t *fft_alloc_complex(size_t n) {
  void *ptr;
#ifdef _WIN32
  ptr = _aligned_malloc(n * sizeof(complex_t), 64);
#else
  if (posix_memalign(&ptr, 64, n * sizeof(complex_t)) != 0) {
    return NULL;
  }
#endif
  return (complex_t *)ptr;
}

double *fft_alloc_real(size_t n) {
  void *ptr;
#ifdef _WIN32
  ptr = _aligned_malloc(n * sizeof(double), 64);
#else
  if (posix_memalign(&ptr, 64, n * sizeof(double)) != 0) {
    return NULL;
  }
#endif
  return (double *)ptr;
}

void fft_free(void *p) {
  if (!p)
    return;
#ifdef _WIN32
  _aligned_free(p);
#else
  free(p);
#endif
}

// Version string
const char *fft_version(void) {
  return "FFT Library v3.1.0 (2026-10-24) - Automatic Algorithm Selection";
}

/*
 * Real-to-complex FFT plan.
 * Stores the real-input pointer and promotes to complex at execute time.
 * Produces n/2+1 unique complex bins (Hermitian symmetry).
 */
fft_plan_t fft_plan_r2c_1d(int n, double *in, complex_t *out, unsigned flags) {
  if (n <= 0 || !in || !out)
    return NULL;

  /* Use an internal complex buffer for the full-size FFT */
  complex_t *cbuf = fft_alloc_complex(n);
  if (!cbuf)
    return NULL;

  fft_plan_t plan =
      fft_plan_dft_1d(n, cbuf, cbuf, -1, flags | FFT_REAL_INPUT);
  if (!plan) {
    fft_free(cbuf);
    return NULL;
  }

  plan->real_in = in;
  plan->out = out;
  plan->algorithm_data = cbuf; /* freed in fft_destroy_plan */
  return plan;
}

/*
 * Complex-to-real FFT plan.
 * Stores the real-output pointer; execute runs inverse FFT then extracts real
 * parts.
 */
fft_plan_t fft_plan_c2r_1d(int n, complex_t *in, double *out, unsigned flags) {
  if (n <= 0 || !in || !out)
    return NULL;

  complex_t *cbuf = fft_alloc_complex(n);
  if (!cbuf)
    return NULL;

  fft_plan_t plan =
      fft_plan_dft_1d(n, in, cbuf, 1, flags | FFT_REAL_OUTPUT);
  if (!plan) {
    fft_free(cbuf);
    return NULL;
  }

  plan->real_out = out;
  plan->algorithm_data = cbuf;
  return plan;
}

/*
 * 2-D FFT plan via row-column decomposition.
 * Stores rows/cols in the plan; execute performs row FFTs then column FFTs
 * (with an internal transpose).
 */
fft_plan_t fft_plan_dft_2d(int rows, int cols, complex_t *in, complex_t *out,
                           int sign, unsigned flags) {
  if (rows <= 0 || cols <= 0 || !in || !out)
    return NULL;
  if ((long long)rows * cols > INT_MAX)
    return NULL;
  if (!is_power_of_two(rows) || !is_power_of_two(cols))
    return NULL;

  detect_hardware();

  fft_plan_t plan = calloc(1, sizeof(struct fft_plan));
  if (!plan)
    return NULL;

  plan->n = rows * cols;
  plan->rows = rows;
  plan->cols = cols;
  plan->in = in;
  plan->out = out;
  plan->dir = (sign < 0) ? FFT_FORWARD : FFT_INVERSE;
  plan->flags = flags;
  plan->algorithm = choose_algorithm(cols, flags); /* per-row algo */
  plan->col_algorithm = choose_algorithm(rows, flags);

  /* Under FFT_WISDOM_ONLY a missing entry makes choose_algorithm() return the
   * ALGO_COUNT sentinel. Reject the plan now rather than failing partway
   * through execution. */
  if (plan->algorithm == ALGO_COUNT || plan->col_algorithm == ALGO_COUNT) {
    free(plan);
    return NULL;
  }

  return plan;
}

/* Export wisdom as a human-readable string.
   Format: one line per entry: "<n> <flags> <algo>\n".
   The header carries the wisdom-format version (v3.0.0), not the library
   version, so files written by 3.0 remain valid in 3.1. */
char *fft_export_wisdom_to_string(void) {
  /* Header + up to WISDOM_MAX entries of ~40 chars each */
  size_t cap = 64 + (size_t)g_wisdom_count * 48;
  char *buf = (char *)malloc(cap);
  if (!buf)
    return NULL;

  int pos = snprintf(buf, cap, "# FFT Wisdom v3.0.0 (%d entries)\n",
                     g_wisdom_count);
  if (pos < 0 || (size_t)pos >= cap) {
    free(buf);
    return NULL;
  }
  /* snprintf returns the length it *would* have written, so accumulating it
   * unguarded could push pos past cap and make cap - pos a huge size_t. */
  for (int i = 0; i < g_wisdom_count; i++) {
    int written = snprintf(buf + pos, cap - (size_t)pos, "%d %u %d\n",
                           g_wisdom[i].n, g_wisdom[i].flags,
                           (int)g_wisdom[i].algo);
    if (written < 0 || (size_t)written >= cap - (size_t)pos) break;
    pos += written;
  }
  return buf;
}

/* Import wisdom from a previously exported string.
   Returns the number of entries imported. */
int fft_import_wisdom_from_string(const char *wisdom) {
  if (!wisdom)
    return 0;

  g_wisdom_count = 0;
  const char *p = wisdom;

  /* Skip header line */
  while (*p && *p != '\n') p++;
  if (*p) p++;

  while (*p && g_wisdom_count < WISDOM_MAX) {
    int n;
    unsigned flags;
    int algo;
    if (sscanf(p, "%d %u %d", &n, &flags, &algo) == 3 && n > 0 &&
        algo >= 0 && algo < ALGO_COUNT &&
        algorithm_handles_size((enum fft_algo)algo, n)) {
      g_wisdom[g_wisdom_count].n = n;
      g_wisdom[g_wisdom_count].flags = flags;
      g_wisdom[g_wisdom_count].algo = (enum fft_algo)algo;
      g_wisdom_count++;
    }
    while (*p && *p != '\n') p++;
    if (*p) p++;
  }
  return g_wisdom_count;
}

#ifndef LIB_BUILD

int main() {
  printf("FFT Auto - Automatic Algorithm Selection v3.1.0\n");
  printf("===============================================\n\n");

  int test_sizes[] = {16, 64, 256, 1024, 97, 360, 1000};
  const char *size_types[] = {"tiny",  "small",     "medium",   "large",
                              "prime", "composite", "arbitrary"};

  printf("Testing automatic algorithm selection:\n");
  printf("Size\tType\t\tAlgorithm\tTime (ms)\n");
  printf("----\t----\t\t---------\t---------\n");

  for (int i = 0; i < 7; i++) {
    int n = test_sizes[i];
    complex_t *data = fft_alloc_complex(n);

    for (int j = 0; j < n; j++)
      data[j] = sin(2 * PI * j / n) + I * cos(2 * PI * j / n);

    clock_t t0 = clock();
    fft_auto(data, data, n, -1);
    clock_t t1 = clock();
    double time_ms = ((double)(t1 - t0) / CLOCKS_PER_SEC) * 1000.0;

    enum fft_algo algo = select_algorithm(n, FFT_ESTIMATE);
    printf("%d\t%s\t\t%s\t%.3f\n", n, size_types[i],
           algorithm_name(algo), time_ms);

    fft_free(data);
  }

  /* Repeated-transform benchmark */
  printf("\n\nTesting planning API for repeated transforms:\n");
  int n = 1024;
  complex_t *in = fft_alloc_complex(n);
  complex_t *out = fft_alloc_complex(n);

  fft_plan_t plan = fft_plan_dft_1d(n, in, out, -1, FFT_MEASURE);
  clock_t t0 = clock();
  for (int i = 0; i < 1000; i++)
    fft_execute(plan);
  clock_t t1 = clock();
  printf("Average time per FFT (1000 iterations): %.3f ms\n",
         ((double)(t1 - t0) / CLOCKS_PER_SEC));

  fft_destroy_plan(plan);
  fft_free(in);
  fft_free(out);

  /* r2c / c2r round-trip test */
  printf("\nReal-to-complex / complex-to-real round-trip (n=%d):\n", n);
  double *rin = fft_alloc_real(n);
  complex_t *cout = fft_alloc_complex(n / 2 + 1);
  double *rout = fft_alloc_real(n);

  for (int i = 0; i < n; i++)
    rin[i] = sin(2 * PI * 3 * i / n);

  fft_plan_t pr2c = fft_plan_r2c_1d(n, rin, cout, FFT_ESTIMATE);

  /* c2r consumes the n/2+1 bins that r2c produced; the mirror is implied. */
  fft_execute(pr2c);
  fft_plan_t pc2r = fft_plan_c2r_1d(n, cout, rout, FFT_ESTIMATE);
  fft_execute(pc2r);

  double max_err = 0;
  for (int i = 0; i < n; i++) {
    double e = fabs(rin[i] - rout[i]);
    if (e > max_err) max_err = e;
  }
  printf("  Max reconstruction error: %.2e\n", max_err);

  fft_destroy_plan(pr2c);
  fft_destroy_plan(pc2r);
  fft_free(rin);
  fft_free(cout);
  fft_free(rout);

  /* Wisdom test */
  printf("\nWisdom export/import test:\n");
  char *wis = fft_export_wisdom_to_string();
  printf("  Exported: %s", wis);
  int imported = fft_import_wisdom_from_string(wis);
  printf("  Imported %d entries\n", imported);
  free(wis);

  /* Hardware capabilities */
  printf("\nDetected Hardware Capabilities:\n");
  unsigned caps = fft_get_hardware_capabilities();
  printf("CPU:");
  if (caps & FFT_HW_CPU_SSE)  printf(" SSE");
  if (caps & FFT_HW_CPU_AVX)  printf(" AVX");
  if (caps & FFT_HW_CPU_AVX2) printf(" AVX2");
  if (caps & FFT_HW_CPU_NEON) printf(" NEON");
  printf("\nGPU:");
  if (caps & FFT_HW_GPU_CUDA) printf(" CUDA");
  if (caps & FFT_HW_GPU_MPS)  printf(" Metal");
  if (!(caps & (FFT_HW_GPU_CUDA | FFT_HW_GPU_MPS))) printf(" None");
  printf("\n");

  return 0;
}
#endif
