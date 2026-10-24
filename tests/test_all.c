#include "../include/fft_algorithms.h"
#include "../include/fft_auto.h"
#include "../include/fft_common.h"
#include <float.h>
#include <stdlib.h>

/**
 * FFT Test Suite
 *
 * Each implementation is compared against a direct DFT computed here, sharing
 * no code with the library, plus algebraic property tests. The reference
 * comparison catches ordering and twiddle-index defects that property tests
 * cannot see.
 */

/* ---------- configuration ---------- */

typedef struct {
  double tolerance;   /* relative tolerance for reference agreement */
  int verbose;
  int stop_on_failure;
} test_config_t;

/* ---------- implementation registry ---------- */

typedef struct {
  const char *name;
  fft_error_t (*fft_func)(complex_t *, int, fft_direction);
  int requires_power_of_2;
  int max_size; /* skip sizes above this; direct DFTs are O(n^2) */
} fft_implementation_t;

static fft_implementation_t implementations[] = {
    {"Radix-2 DIT", radix2_dit_fft, 1, 16384},
    {"Radix-2 DIF", radix2_dif_fft, 1, 16384},
    {"Radix-4", radix4_fft, 1, 16384},
    {"Split-Radix", split_radix_fft, 1, 16384},
    {"Bluestein", bluestein_fft, 0, 4096},
    {"Mixed-Radix", mixed_radix_fft, 0, 4096},
    {"Recursive", fft_recursive, 1, 4096},
    {"Naive DFT", naive_dft, 0, 256}};

static const int num_implementations =
    sizeof(implementations) / sizeof(implementations[0]);

/* ---------- counters ---------- */

static int g_total_passed = 0;
static int g_total_tests = 0;
static int g_failures = 0;

static void report(const char *scope, const char *test, int passed,
                   double error, const char *detail) {
  g_total_tests++;
  if (passed) {
    g_total_passed++;
  } else {
    g_failures++;
    printf("  [FAIL] %-14s %-24s err=%.3e  %s\n", scope, test, error,
           detail ? detail : "");
  }
}

/* ---------- independent reference DFT ---------- */
/*
 * Computed by direct summation with trigonometry obtained straight from the
 * definition of the DFT. Deliberately does not use twiddle_factor(), the
 * library's bit_reverse(), or any library allocation helper, so a defect in
 * those cannot hide by being reproduced on both sides of the comparison.
 */

/* One DFT bin in O(n). Checking selected bins gives an independent reference
 * without paying the O(n^2) cost of a full direct transform at large n. */
static complex_t reference_bin(const complex_t *in, int n, int k,
                               fft_direction dir) {
  complex_t acc = 0.0;
  for (int i = 0; i < n; i++) {
    double angle =
        (double)dir * 2.0 * PI * (double)((long)i * k % (long)n) / (double)n;
    acc += in[i] * (cos(angle) + I * sin(angle));
  }
  return (dir == FFT_INVERSE) ? acc / (double)n : acc;
}

static double rel_err(complex_t a, complex_t b) {
  double s = cabs(b);
  return cabs(a - b) / (s > 1.0 ? s : 1.0);
}

/* Independent of the library's log2_int() on purpose. */
static int ilog2i(int n) {
  int l = 0;
  while (n > 1) { n >>= 1; l++; }
  return l;
}

/*
 * Cooley-Tukey is backward stable with rms error O(eps*sqrt(log n)), so the
 * expected deviation from a direct DFT grows slowly with n. A fixed 1e-10
 * would flag perfectly good double-precision results at large n; scaling the
 * bound keeps the test sensitive to the defect class that matters (ordering
 * and twiddle-index errors produce O(1) deviations, not O(1e-9)).
 */
static double tol_for(const test_config_t *c, int n) {
  return c->tolerance * (1.0 + (double)ilog2i(n > 1 ? n : 2));
}

static void fill_random(complex_t *x, int n, unsigned seed) {
  srand(seed);
  for (int i = 0; i < n; i++) {
    x[i] = ((double)rand() / RAND_MAX - 0.5) +
           I * ((double)rand() / RAND_MAX - 0.5);
  }
}

/* ---------- test: agreement with the reference DFT ---------- */

/*
 * Checks a specific set of output bins against directly-computed values.
 * For n <= 1024 every bin is checked; above that a strided sample plus the
 * edges, which keeps the reference at O(n * bins) instead of O(n^2) while
 * still pinning the bin indices (a permutation error moves a checked bin).
 */
static void test_reference_agreement(fft_implementation_t *impl, int n,
                                     test_config_t *config) {
  complex_t *work = malloc((size_t)n * sizeof(complex_t));
  complex_t *orig = malloc((size_t)n * sizeof(complex_t));
  if (!work || !orig) {
    report(impl->name, "reference alloc failed", 0, 0.0, "out of memory");
    free(work);
    free(orig);
    return;
  }

  fill_random(orig, n, 100000u + (unsigned)n);

  int step = (n <= 1024) ? 1 : n / 64;
  if (step < 1) step = 1;

  struct { int idx; fft_direction dir; } checks[512];
  int nchecks = 0;
  for (int k = 0; k < n && nchecks < 480; k += step) {
    checks[nchecks].idx = k;
    checks[nchecks].dir = FFT_FORWARD;
    nchecks++;
  }
  /* Edges and the Nyquist bin, where ordering defects show up first. */
  int extra[] = {0, 1, n / 2, n / 2 + 1, n - 2, n - 1};
  for (int e = 0; e < 6 && nchecks < 510; e++) {
    int k = extra[e];
    if (k < 0 || k >= n) continue;
    int dup = 0;
    for (int q = 0; q < nchecks; q++)
      if (checks[q].idx == k && checks[q].dir == FFT_FORWARD) dup = 1;
    if (!dup) {
      checks[nchecks].idx = k;
      checks[nchecks].dir = FFT_FORWARD;
      nchecks++;
    }
  }

  /* Forward */
  memcpy(work, orig, (size_t)n * sizeof(complex_t));
  impl->fft_func(work, n, FFT_FORWARD);
  double worst = 0.0;
  int worst_bin = -1;
  for (int q = 0; q < nchecks; q++) {
    int k = checks[q].idx;
    double e = rel_err(work[k], reference_bin(orig, n, k, FFT_FORWARD));
    if (e > worst) {
      worst = e;
      worst_bin = k;
    }
  }
  char detail[128];
  snprintf(detail, sizeof(detail), "forward, worst bin %d of %d checked",
           worst_bin, nchecks);
  double bound = tol_for(config, n);
  int ok = worst <= bound;
  report(impl->name, "reference agreement fwd", ok, worst, ok ? NULL : detail);

  /* Inverse */
  memcpy(work, orig, (size_t)n * sizeof(complex_t));
  impl->fft_func(work, n, FFT_INVERSE);
  worst = 0.0;
  for (int q = 0; q < nchecks; q++) {
    int k = checks[q].idx;
    double e = rel_err(work[k], reference_bin(orig, n, k, FFT_INVERSE));
    if (e > worst) worst = e;
  }
  ok = worst <= bound;
  report(impl->name, "reference agreement inv", ok, worst,
         ok ? NULL : "inverse transform disagrees");

  free(work);
  free(orig);
}

/* ---------- property tests ---------- */

static void test_impulse_response(fft_implementation_t *impl, int n,
                                  test_config_t *config) {
  complex_t *signal = calloc((size_t)n, sizeof(complex_t));
  signal[0] = 1.0;
  impl->fft_func(signal, n, FFT_FORWARD);

  double worst = 0.0;
  for (int i = 0; i < n; i++) {
    double e = fabs(cabs(signal[i]) - 1.0);
    if (e > worst) worst = e;
  }
  int ok = worst <= tol_for(config, n);
  report(impl->name, "impulse response", ok, worst, ok ? NULL : "|X[k]| != 1");
  free(signal);
}

static void test_dc_signal(fft_implementation_t *impl, int n,
                           test_config_t *config) {
  complex_t *signal = malloc((size_t)n * sizeof(complex_t));
  for (int i = 0; i < n; i++) signal[i] = 1.0;
  impl->fft_func(signal, n, FFT_FORWARD);

  /* Normalise by the transform scale (the DC bin is n), so round-off on the
   * near-zero bins is measured relative to what the transform can resolve. */
  double worst = 0.0;
  for (int i = 0; i < n; i++) {
    double want = (i == 0) ? (double)n : 0.0;
    double e = cabs(signal[i] - want) / (double)n;
    if (e > worst) worst = e;
  }
  int ok = worst <= tol_for(config, n);
  report(impl->name, "dc signal", ok, worst, ok ? NULL : "DC bin wrong");
  free(signal);
}

static void test_linearity(fft_implementation_t *impl, int n,
                           test_config_t *config) {
  complex_t *a = malloc((size_t)n * sizeof(complex_t));
  complex_t *b = malloc((size_t)n * sizeof(complex_t));
  complex_t *s = malloc((size_t)n * sizeof(complex_t));
  complex_t *fs = malloc((size_t)n * sizeof(complex_t));
  fill_random(a, n, 5000u + (unsigned)n);
  fill_random(b, n, 9000u + (unsigned)n);
  for (int i = 0; i < n; i++) s[i] = 2 * a[i] + 3 * b[i];

  memcpy(fs, s, (size_t)n * sizeof(complex_t));
  impl->fft_func(fs, n, FFT_FORWARD);
  impl->fft_func(a, n, FFT_FORWARD);
  impl->fft_func(b, n, FFT_FORWARD);

  double worst = 0.0;
  for (int i = 0; i < n; i++) {
    double e = rel_err(fs[i], 2 * a[i] + 3 * b[i]);
    if (e > worst) worst = e;
  }
  int ok = worst <= tol_for(config, n);
  report(impl->name, "linearity", ok, worst, ok ? NULL : "FFT(2a+3b) mismatch");
  free(a); free(b); free(s); free(fs);
}

static void test_parsevals(fft_implementation_t *impl, int n,
                           test_config_t *config) {
  complex_t *signal = malloc((size_t)n * sizeof(complex_t));
  fill_random(signal, n, 6000u + (unsigned)n);

  double time_energy = 0.0;
  for (int i = 0; i < n; i++) time_energy += cabs(signal[i]) * cabs(signal[i]);

  impl->fft_func(signal, n, FFT_FORWARD);
  double freq_energy = 0.0;
  for (int i = 0; i < n; i++) freq_energy += cabs(signal[i]) * cabs(signal[i]);
  freq_energy /= n;

  double rel = fabs(time_energy - freq_energy) / time_energy;
  /* Energy is a sum of magnitudes, so it tolerates no ordering error either. */
  int ok = rel <= tol_for(config, n);
  report(impl->name, "parseval theorem", ok, rel, ok ? NULL : "energy mismatch");
  free(signal);
}

static void test_inverse_roundtrip(fft_implementation_t *impl, int n,
                                   test_config_t *config) {
  complex_t *original = malloc((size_t)n * sizeof(complex_t));
  complex_t *work = malloc((size_t)n * sizeof(complex_t));
  fill_random(original, n, 7000u + (unsigned)n);
  memcpy(work, original, (size_t)n * sizeof(complex_t));

  impl->fft_func(work, n, FFT_FORWARD);
  impl->fft_func(work, n, FFT_INVERSE);

  double worst = 0.0;
  for (int i = 0; i < n; i++) {
    double e = rel_err(work[i], original[i]);
    if (e > worst) worst = e;
  }
  double tol = tol_for(config, n);
  int ok = worst <= tol;
  report(impl->name, "inverse roundtrip", ok, worst,
         ok ? NULL : "x != IFFT(FFT(x))");
  free(original); free(work);
}

static void test_known_pair(fft_implementation_t *impl, int n,
                            test_config_t *config) {
  if (n < 4) return; /* cos at bin 1 needs at least 4 points to be distinct */
  complex_t *signal = malloc((size_t)n * sizeof(complex_t));
  int freq = 1;
  for (int i = 0; i < n; i++) signal[i] = cos(2 * PI * freq * i / n);
  impl->fft_func(signal, n, FFT_FORWARD);

  double expected = n / 2.0;
  double worst = 0.0;
  for (int i = 0; i < n; i++) {
    double want = (i == freq || i == n - freq) ? expected : 0.0;
    double e = cabs(signal[i] - want) / expected;
    if (e > worst) worst = e;
  }
  /* Property-only check; the reference test is what pins the bin indices. */
  int ok = worst <= tol_for(config, n);
  report(impl->name, "known pair cos", ok, worst,
         ok ? NULL : "energy in wrong bins");
  free(signal);
}

static void test_stability(fft_implementation_t *impl, int n,
                           test_config_t *config) {
  if (n > 1024) return; /* repeated round trips are costly for O(n^2) paths */
  complex_t *signal = malloc((size_t)n * sizeof(complex_t));
  complex_t *expect = malloc((size_t)n * sizeof(complex_t));
  fill_random(expect, n, 8000u + (unsigned)n);

  /* Wide dynamic range: 1e-5 .. 1e5. expect[] holds the same scaling so the
   * comparison is like-for-like rather than against unscaled values. */
  for (int i = 0; i < n; i++) {
    double scale = pow(10, -5 + 10.0 * i / n);
    expect[i] *= scale;
    signal[i] = expect[i];
  }

  for (int iter = 0; iter < 10; iter++) {
    impl->fft_func(signal, n, FFT_FORWARD);
    impl->fft_func(signal, n, FFT_INVERSE);
  }

  double worst = 0.0;
  for (int i = 0; i < n; i++) worst = fmax(worst, rel_err(signal[i], expect[i]));
  /* Ten round trips multiply the round-off, hence the looser bound than the
   * single-pass reference test uses. */
  double tol = fmax(1e-6, tol_for(config, n) * 10.0);
  int ok = worst <= tol;
  report(impl->name, "10x roundtrip stability", ok, worst,
         ok ? NULL : "error accumulates across round trips");
  free(signal); free(expect);
}

/* ---------- tests that are not per-implementation ---------- */

static void test_bit_reverse_network(void) {
  printf("\n=== bit_reverse() permutation network ===\n");
  for (int lg = 1; lg <= 16; lg++) {
    int n = 1 << lg, bad = 0;
    unsigned ex_i = 0, ex_got = 0, ex_want = 0;
    for (int i = 0; i < n; i++) {
      unsigned x = (unsigned)i, want = 0;
      for (int b = 0; b < lg; b++) {
        want = (want << 1) | (x & 1u);
        x >>= 1;
      }
      unsigned got = bit_reverse((unsigned)i, lg);
      if (got != want) {
        if (!bad) { ex_i = (unsigned)i; ex_got = got; ex_want = want; }
        bad++;
      }
    }
    char detail[128];
    snprintf(detail, sizeof(detail), "%d/%d wrong, e.g. bit_reverse(%u,%d)=%u want %u",
             bad, n, ex_i, lg, ex_got, ex_want);
    report("bit_reverse", detail, bad == 0, bad ? 1.0 : 0.0,
           bad ? "index permutation broken" : NULL);
  }
}

/*
 * log2_int() must not loop on negative input: an arithmetic right shift of a
 * negative value never becomes zero. A regression here hangs the process, which
 * the harness timeout catches.
 */
static void test_index_helpers(void) {
  printf("\n=== index/size helpers ===\n");
  struct { int in; int want; } lg[] = {
      {-7, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 1}, {3, 1}, {4, 2},
      {1024, 10}, {1 << 20, 20}};
  for (unsigned q = 0; q < sizeof(lg) / sizeof(lg[0]); q++) {
    int got = log2_int(lg[q].in);
    char d[96];
    snprintf(d, sizeof(d), "log2_int(%d)=%d want %d", lg[q].in, got, lg[q].want);
    report("log2_int", d, got == lg[q].want, 0.0, got == lg[q].want ? NULL : d);
  }

  struct { int in; int want; } np[] = {
      {-5, 1}, {0, 1}, {1, 1}, {2, 2}, {3, 4}, {5, 8}, {17, 32}, {1024, 1024},
      {1025, 2048}};
  for (unsigned q = 0; q < sizeof(np) / sizeof(np[0]); q++) {
    int got = next_power_of_two(np[q].in);
    int ok = (got == np[q].want) ||
             (got > np[q].in - 1 && (got & (got - 1)) == 0 && got >= np[q].in);
    char d[96];
    snprintf(d, sizeof(d), "next_power_of_two(%d)=%d want %d", np[q].in, got,
             np[q].want);
    report("next_power_of_two", d, ok, 0.0, ok ? NULL : d);
  }

  /* n <= 0 must not attempt a huge allocation via int -> size_t underflow. */
  int ok_null = (allocate_complex_array(0) == NULL) &&
                (allocate_complex_array(-4) == NULL);
  report("allocator", "n<=0 returns NULL", ok_null, 0.0,
         ok_null ? NULL : "negative size reached calloc");
}

/* A size that is not supported must be reported, never silently ignored. */
static void test_rejects_invalid_size(void) {
  printf("\n=== error-code reporting ===\n");
  complex_t buf[8];
  for (int i = 0; i < 8; i++) buf[i] = 1.0;

  /*
   * Before this API change the transforms reported a bad length by calling
   * exit(), which killed whatever program linked the library. Each one now
   * returns a code and must leave the caller's buffer untouched.
   */
  for (int q = 0; q < num_implementations; q++) {
    fft_implementation_t *impl = &implementations[q];
    char scope[64];
    snprintf(scope, sizeof(scope), "%s errors", impl->name);

    int ok = impl->fft_func(buf, 8, FFT_FORWARD) == FFT_OK;

    /* Every rejection below happens before any write, so the buffer must still
     * hold the values set above. */
    for (int i = 0; i < 8; i++) buf[i] = 1.0;
    ok = ok && impl->fft_func(NULL, 8, FFT_FORWARD) == FFT_ERR_NULL;
    ok = ok && impl->fft_func(buf, 0, FFT_FORWARD) == FFT_ERR_SIZE;
    ok = ok && impl->fft_func(buf, -1, FFT_FORWARD) == FFT_ERR_SIZE;
    if (impl->requires_power_of_2)
      ok = ok && impl->fft_func(buf, 3, FFT_FORWARD) == FFT_ERR_SIZE;

    for (int i = 0; i < 8; i++)
      ok = ok && creal(buf[i]) == 1.0 && cimag(buf[i]) == 0.0;

    report(scope, "returns code, not exit()", ok, 0.0, NULL);
  }
}

/* ---------- test: planning API, wisdom cache, real and 2-D paths ---------- */

static int algorithm_name_is_known(const char *name) {
  static const char *const known[] = {"Radix-2 DIT", "Radix-2 DIF", "Radix-4",
                                      "Split-Radix", "Bluestein", "Mixed-Radix",
                                      "GPU (CUDA)",  "GPU (Metal)"};
  if (!name) return 0;
  for (size_t i = 0; i < sizeof(known) / sizeof(known[0]); i++)
    if (strcmp(name, known[i]) == 0) return 1;
  return 0;
}

static void test_plan_api(test_config_t *config) {
  printf("\n=== planning API and wisdom ===\n");
  const int n = 256;
  const double tol = tol_for(config, n);

  complex_t *src = malloc((size_t)n * sizeof(complex_t));
  complex_t *work = malloc((size_t)n * sizeof(complex_t));
  if (!src || !work) {
    report("plan api", "alloc", 0, 0.0, "out of memory");
    free(src);
    free(work);
    return;
  }
  fill_random(src, n, 4242u);

  /*
   * The flags are OR-combined and tested with `flags & FFT_MEASURE`, so each
   * must be a single distinct bit and no two may alias.
   */
  {
    static const unsigned planning_flags[] = {
        FFT_MEASURE,     FFT_WISDOM_ONLY, FFT_REAL_INPUT,
        FFT_REAL_OUTPUT, FFT_CONSERVE_MEMORY, FFT_PREFER_GPU};
    const int nf = (int)(sizeof(planning_flags) / sizeof(planning_flags[0]));
    int single_bit = 1, disjoint = 1;
    unsigned seen = 0;
    for (int i = 0; i < nf; i++) {
      unsigned f = planning_flags[i];
      if (f == 0 || (f & (f - 1)) != 0) single_bit = 0; /* not one bit */
      if (seen & f) disjoint = 0;                        /* overlaps another */
      seen |= f;
    }
    report("plan api", "planning flags are distinct single bits",
           single_bit && disjoint, 0.0, NULL);
  }

  /* A null plan and an unusable length report instead of crashing. */
  report("plan api", "null plan rejected", fft_execute(NULL) == FFT_ERR_NULL &&
                                             fft_execute_dft(NULL, src, work) ==
                                                 FFT_ERR_NULL,
         0.0, NULL);
  report("plan api", "n<=0 not plannable",
         fft_plan_dft_1d(0, src, work, -1, FFT_ESTIMATE) == NULL &&
             fft_plan_dft_1d(-4, src, work, -1, FFT_ESTIMATE) == NULL,
         0.0, NULL);

  static const struct {
    unsigned flags;
    const char *what;
  } modes[] = {{FFT_ESTIMATE, "heuristic plan"}, {FFT_MEASURE, "measured plan"}};

  char measured_name[64];
  measured_name[0] = '\0';

  for (size_t m = 0; m < sizeof(modes) / sizeof(modes[0]); m++) {
    memcpy(work, src, (size_t)n * sizeof(complex_t));
    fft_plan_t plan = fft_plan_dft_1d(n, work, work, -1, modes[m].flags);
    if (!plan) {
      report("plan api", modes[m].what, 0, 0.0, "planning returned NULL");
      continue;
    }

    int ok = fft_execute(plan) == FFT_OK;
    ok = ok && algorithm_name_is_known(fft_plan_algorithm_name(plan));

    /* The planned transform must agree with the direct DFT just like the
     * hand-selected kernels do; this is the only end-to-end check of the
     * dispatch layer. */
    for (int k = 0; k < n && ok; k += 17)
      ok = rel_err(work[k], reference_bin(src, n, k, FFT_FORWARD)) < tol;

    if (modes[m].flags == FFT_MEASURE)
      snprintf(measured_name, sizeof(measured_name), "%s",
               fft_plan_algorithm_name(plan));

    fft_destroy_plan(plan);
    report("plan api", modes[m].what, ok, 0.0, NULL);
  }

  /* Wisdom: the measured choice survives an export/import round trip. */
  {
    int ok = measured_name[0] != '\0';
    char *wis = fft_export_wisdom_to_string();
    ok = ok && wis != NULL;
    if (ok) ok = strncmp(wis, "# FFT Wisdom", 12) == 0;
    if (ok) ok = fft_import_wisdom_from_string(wis) >= 1;
    free(wis);

    fft_plan_t plan = fft_plan_dft_1d(n, work, work, -1, FFT_MEASURE);
    if (ok)
      ok = plan != NULL &&
           strcmp(fft_plan_algorithm_name(plan), measured_name) == 0;
    fft_destroy_plan(plan);
    report("plan api", "wisdom round trip", ok, 0.0, NULL);
  }

  /*
   * FFT_WISDOM_ONLY must refuse to plan rather than silently fall back to a
   * fresh choice -- the point of the flag is that planning may not measure.
   *
   * The cold case uses FFT_CONSERVE_MEMORY, which no plan anywhere in this
   * suite is built with, so its wisdom key is provably absent. Without that the
   * test would depend on the global cache being empty, which it is not.
   */
  {
    fft_plan_t cold = fft_plan_dft_1d(
        n, work, work, -1, FFT_CONSERVE_MEMORY | FFT_WISDOM_ONLY);
    report("plan api", "wisdom-only on a cold cache returns NULL", cold == NULL,
           0.0, NULL);
    fft_destroy_plan(cold); /* must tolerate NULL */
  }

  /*
   * The warm case is the one that exercises the subtlety: the flag is masked out
   * of the cache key, so a choice made under FFT_MEASURE is found by a caller
   * asking for FFT_MEASURE|FFT_WISDOM_ONLY. Keying on the flags verbatim would
   * make that combination always miss, and the flag would be useless.
   */
  {
    int ok = measured_name[0] != '\0';
    memcpy(work, src, (size_t)n * sizeof(complex_t));
    fft_plan_t warm =
        fft_plan_dft_1d(n, work, work, -1, FFT_MEASURE | FFT_WISDOM_ONLY);
    if (!warm) {
      report("plan api", "wisdom-only on a warm cache plans and runs", 0, 0.0,
             "planning returned NULL");
    } else {
      ok = ok && strcmp(fft_plan_algorithm_name(warm), measured_name) == 0;
      ok = ok && fft_execute(warm) == FFT_OK;
      for (int k = 0; k < n && ok; k += 17)
        ok = rel_err(work[k], reference_bin(src, n, k, FFT_FORWARD)) < tol;
      fft_destroy_plan(warm);
      report("plan api", "wisdom-only on a warm cache plans and runs", ok, 0.0,
             NULL);
    }
  }

  /* A hostile cache file must not install an algorithm id that execute_1d
   * would dispatch on, nor a length no buffer could satisfy. */
  {
    int imported =
        fft_import_wisdom_from_string("# FFT Wisdom v3.0.0 (2 entries)\n"
                                      "-5 0 1\n"
                                      "256 0 9999\n");
    report("plan api", "rejects corrupt wisdom", imported == 0, 0.0, NULL);
  }

  /*
   * Wisdom is an external input format: fft_import_wisdom_from_string() reads
   * text that may pair any algorithm with any n. An entry whose algorithm
   * cannot transform the recorded length must be dropped at the door. Radix-4
   * (id 2) on a prime length is the concrete case: cached unchecked, it yields a
   * non-NULL plan whose every execute returns FFT_ERR_SIZE.
   */
  {
    int imported = fft_import_wisdom_from_string(
        "# FFT Wisdom v3.0.0 (1 entries)\n7 0 2\n");
    complex_t seven[7], original[7];
    for (int i = 0; i < 7; i++)
      seven[i] = original[i] = (double)(i + 1) + I * (double)i;

    fft_plan_t bad = fft_plan_dft_1d(7, seven, seven, -1, FFT_ESTIMATE);
    int ok = imported == 0 && bad != NULL;
    if (bad) {
      ok = ok && strcmp(fft_plan_algorithm_name(bad), "Bluestein") == 0;
      ok = ok && fft_execute(bad) == FFT_OK;
    }
    for (int k = 0; k < 7 && ok; k++)
      ok = rel_err(seven[k], reference_bin(original, 7, k, FFT_FORWARD)) <
           tol_for(config, 7);
    fft_destroy_plan(bad);
    report("plan api", "unexecutable wisdom pairing rejected", ok, 0.0, NULL);
  }

  /*
   * A GPU id is well-formed and portable -- another machine may have a GPU --
   * so import accepts it. It must not be returned here, where detect_hardware()
   * found none: fft_execute()'s GPU case is compiled out in this build and
   * would report FFT_OK without touching the buffer.
   */
  {
    int imported = fft_import_wisdom_from_string(
        "# FFT Wisdom v3.0.0 (1 entries)\n16 0 6\n");
    unsigned caps = fft_get_hardware_capabilities();
    int has_gpu = (caps & (FFT_HW_GPU_CUDA | FFT_HW_GPU_MPS)) != 0;

    complex_t gpu_buf[16];
    for (int i = 0; i < 16; i++) gpu_buf[i] = 100.0 + (double)i;
    fft_plan_t gpu = fft_plan_dft_1d(16, gpu_buf, gpu_buf, -1, FFT_WISDOM_ONLY);

    int ok = imported == 1;
    if (has_gpu)
      ok = ok && gpu != NULL && fft_execute(gpu) == FFT_OK;
    else
      ok = ok && gpu == NULL; /* unavailable entry must be a cache miss */
    fft_destroy_plan(gpu);
    report("plan api", "unavailable GPU wisdom not executed", ok, 0.0, NULL);
  }

  /*
   * FFT_CONSERVE_MEMORY is a constraint. Radix-4 is faster at 512 and up but
   * allocates an O(n) digit-reversal buffer; radix-2 DIT is in place. The flag
   * must pick radix-2 for power-of-two sizes, and FFT_MEASURE must not override
   * it back to the allocating kernel.
   */
  {
    const int cn = 512;
    complex_t *csrc = malloc((size_t)cn * sizeof(complex_t));
    complex_t *cwork = malloc((size_t)cn * sizeof(complex_t));
    if (!csrc || !cwork) {
      report("plan api", "conserve-memory alloc", 0, 0.0, "out of memory");
    } else {
      static const char *const inplace[] = {"Radix-2 DIT", "Radix-2 DIF"};
      const double ctol = tol_for(config, cn);
      fill_random(csrc, cn, 99u);
      int ok = 1;

      /* The GPU test above left an entry in the cache. */
      fft_import_wisdom_from_string("# FFT Wisdom v3.0.0 (0 entries)\n");
      for (int mode = 0; mode < 2; mode++) {
        unsigned f = FFT_CONSERVE_MEMORY;
        if (mode == 1) f |= FFT_MEASURE;
        memcpy(cwork, csrc, (size_t)cn * sizeof(complex_t));
        fft_plan_t cplan = fft_plan_dft_1d(cn, cwork, cwork, -1, f);
        int inplace_choice = 0;
        if (cplan)
          for (size_t i = 0; i < sizeof(inplace) / sizeof(inplace[0]); i++)
            inplace_choice |=
                strcmp(fft_plan_algorithm_name(cplan), inplace[i]) == 0;
        ok = ok && cplan != NULL && inplace_choice;
        ok = ok && fft_execute(cplan) == FFT_OK;
        for (int k = 0; k < cn && ok; k += 37)
          ok = rel_err(cwork[k], reference_bin(csrc, cn, k, FFT_FORWARD)) < ctol;
        fft_destroy_plan(cplan);
      }
      report("plan api", "conserve-memory picks in-place kernel", ok, 0.0, NULL);
    }
    free(csrc);
    free(cwork);
  }

  /* Real-to-complex then complex-to-real on a constant signal: the spectrum is
   * one non-zero bin and the round trip returns the input. */
  {
    const int rn = 64;
    double *rin = malloc((size_t)rn * sizeof(double));
    double *rout = malloc((size_t)rn * sizeof(double));
    complex_t *cout = malloc((size_t)(rn / 2 + 1) * sizeof(complex_t));
    if (!rin || !rout || !cout) {
      report("plan api", "r2c alloc", 0, 0.0, "out of memory");
    } else {
      for (int i = 0; i < rn; i++) rin[i] = 1.0;
      int ok = 1;

      fft_plan_t p = fft_plan_r2c_1d(rn, rin, cout, FFT_ESTIMATE);
      ok = ok && p != NULL && fft_execute(p) == FFT_OK;
      /* A real plan has no complex in/out pair to substitute; execute_dft
       * must refuse rather than transform the arrays it was planned with. */
      ok = ok && fft_execute_dft(p, cout, cout) == FFT_ERR_UNSUPPORTED;
      fft_destroy_plan(p);

      ok = ok && rel_err(cout[0], (double)rn) < 1e-12;
      ok = ok && cabs(cout[3]) < 1e-12;

      /* in[] carries the n/2+1 bins; c2r mirrors the rest. */
      for (int k = 1; k <= rn / 2; k++) cout[k] = 0.0;
      p = fft_plan_c2r_1d(rn, cout, rout, FFT_ESTIMATE);
      ok = ok && p != NULL && fft_execute(p) == FFT_OK;
      ok = ok && fft_execute_dft(p, cout, cout) == FFT_ERR_UNSUPPORTED;
      fft_destroy_plan(p);

      for (int i = 0; i < rn && ok; i++)
        ok = fabs(rout[i] - 1.0) < 1e-12;

      report("plan api", "r2c / c2r round trip", ok, 0.0, NULL);
    }
    free(rin);
    free(rout);
    free(cout);
  }

  /*
   * 2-D on a 4x8 grid. The dimensions are deliberately different: a square
   * plan uses the same algorithm for both passes, so it cannot show whether
   * execution uses the choice stored for the column pass. Wisdom supplies
   * Radix-4 for 4 points and Split-Radix for 8, then is cleared before
   * execute. The plan must still produce the direct 2-D DFT.
   */
  {
    const int rows = 4;
    const int cols = 8;
    const int total = rows * cols;
    complex_t *gin = malloc((size_t)total * sizeof(complex_t));
    complex_t *gout = malloc((size_t)total * sizeof(complex_t));
    if (!gin || !gout) {
      report("plan api", "2d alloc", 0, 0.0, "out of memory");
    } else {
      fill_random(gin, total, 777u);

      fft_import_wisdom_from_string(
          "# FFT Wisdom v3.0.0 (2 entries)\n"
          "4 0 2\n"
          "8 0 3\n");
      fft_plan_t p =
          fft_plan_dft_2d(rows, cols, gin, gout, -1, FFT_WISDOM_ONLY);
      int ok = p != NULL && fft_execute(p) == FFT_OK;

      /* Clearing the cache used to make execute re-select the column pass and
       * return FFT_ERR_SIZE under FFT_WISDOM_ONLY. */
      fft_import_wisdom_from_string("# FFT Wisdom v3.0.0 (0 entries)\n");
      ok = ok && fft_execute(p) == FFT_OK;

      if (ok) {
        const double tol = tol_for(config, total);
        for (int kr = 0; kr < rows && ok; kr++) {
          for (int kc = 0; kc < cols && ok; kc++) {
            complex_t expected = 0.0;
            for (int r = 0; r < rows; r++) {
              for (int c = 0; c < cols; c++) {
                double angle =
                    -TWO_PI * ((double)kr * r / rows + (double)kc * c / cols);
                expected +=
                    gin[r * cols + c] * (cos(angle) + I * sin(angle));
              }
            }
            ok = rel_err(gout[kr * cols + kc], expected) < tol;
          }
        }
      }

      ok = ok && fft_plan_dft_2d(65536, 65536, gout, gout, -1,
                                 FFT_ESTIMATE) == NULL;
      fft_destroy_plan(p);

      report("plan api", "2-D transform", ok, 0.0, NULL);
    }
    free(gin);
    free(gout);
  }

  free(src);
  free(work);
}

/* ---------- driver ---------- */

static int size_supported(fft_implementation_t *impl, int n) {
  if (n > impl->max_size) return 0;
  if (impl->requires_power_of_2 && !is_power_of_two(n)) return 0;
  return 1;
}

int main(int argc, char *argv[]) {
  printf("FFT Implementation Test Suite\n");
  printf("=============================\n");

  test_config_t config;
  config.tolerance = 1e-10;
  config.verbose = 0;
  config.stop_on_failure = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
      config.verbose = 1;
    } else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--stop") == 0) {
      config.stop_on_failure = 1;
    } else if (strcmp(argv[i], "-t") == 0 ||
               strcmp(argv[i], "--tolerance") == 0) {
      if (i + 1 < argc) config.tolerance = atof(argv[++i]);
    } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      printf("Usage: %s [-v] [-s] [-t tolerance]\n", argv[0]);
      return 0;
    }
  }

  printf("\n  tolerance: %.2e (relative)\n", config.tolerance);

  test_bit_reverse_network();
  test_index_helpers();
  test_rejects_invalid_size();
  test_plan_api(&config);

  /* Sizes outer so the O(n^2) reference is computed once per size.
   * Includes 1, primes, 2*4^k (which only a true mixed-radix-4 handles),
   * non-power-of-two composites, and sizes above 1024. */
  static const int sizes[] = {1,    2,    3,    4,    5,    6,    7,    8,
                              12,   15,   16,   20,   24,   30,   31,   32,
                              60,   64,   97,   128,  251,  256,  360,  512,
                              1000, 1021, 1024, 2048, 4096, 16384};
  const int num_sizes = (int)(sizeof(sizes) / sizeof(sizes[0]));

  for (int s = 0; s < num_sizes; s++) {
    int n = sizes[s];
    if (config.verbose) printf("\n--- n=%d ---\n", n);

    for (int q = 0; q < num_implementations; q++) {
      fft_implementation_t *impl = &implementations[q];
      if (!size_supported(impl, n)) continue;

      char scope[64];
      snprintf(scope, sizeof(scope), "%s n=%d", impl->name, n);
      (void)scope;

      test_reference_agreement(impl, n, &config);
      test_impulse_response(impl, n, &config);
      test_dc_signal(impl, n, &config);
      test_linearity(impl, n, &config);
      test_parsevals(impl, n, &config);
      test_inverse_roundtrip(impl, n, &config);
      test_known_pair(impl, n, &config);
      test_stability(impl, n, &config);

      if (config.stop_on_failure && g_failures > 0) {
        printf("\nStopping on first failure.\n");
        return 1;
      }
    }
  }

  printf("\n=== Overall Summary ===\n");
  printf("%d/%d checks passed (%.1f%%)\n", g_total_passed, g_total_tests,
         g_total_tests ? 100.0 * g_total_passed / g_total_tests : 0.0);
  if (g_failures == 0) {
    printf("ALL CHECKS PASSED\n");
    return 0;
  }
  printf("%d CHECKS FAILED\n", g_failures);
  return 1;
}
