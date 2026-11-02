#include "../include/fft_gpu.h"
#include "../include/fft_auto.h"

/* Metal backend checks. The GPU works in float32 while complex_t is double, so
 * results are compared to a direct DFT with a tolerance rather than exactly.
 * Built only when FFT_ENABLE_METAL is on; at runtime a machine without a Metal
 * device (or older than macOS 14) reports a skip instead of failing. */

static int g_checks = 0;
static int g_failures = 0;

static void check(int condition, const char *what) {
  g_checks++;
  if (!condition) {
    g_failures++;
    printf("  [FAIL] %s\n", what);
  }
}

static void check_relative(const char *what, double max_err, double scale,
                           double tol) {
  g_checks++;
  if (scale < 1.0) scale = 1.0;
  if (!(max_err <= tol * scale)) {
    g_failures++;
    printf("  [FAIL] %s (max error %.3e, scale %.3e, tolerance %.3e)\n", what,
           max_err, scale, tol * scale);
  }
}

static unsigned g_seed = 12345u;

static double next_value(void) {
  g_seed = g_seed * 1664525u + 1013904223u;
  return (double)((g_seed >> 8) & 0xffffu) / 32768.0 - 1.0;
}

static void fill_signal(complex_t *x, int n) {
  for (int i = 0; i < n; i++)
    x[i] = next_value() + I * next_value();
}

static double max_magnitude(const complex_t *x, int n) {
  double m = 0.0;
  for (int i = 0; i < n; i++) {
    double a = cabs(x[i]);
    if (a > m) m = a;
  }
  return m;
}

static double max_difference(const complex_t *a, const complex_t *b, int n) {
  double m = 0.0;
  for (int i = 0; i < n; i++) {
    double d = cabs(a[i] - b[i]);
    if (d > m) m = d;
  }
  return m;
}

/* Direct transform used as the reference. inverse includes 1/n, matching the
 * library's convention. */
static void dft_1d_reference(const complex_t *in, complex_t *out, int n,
                             fft_direction direction) {
  for (int k = 0; k < n; k++) {
    complex_t sum = 0.0;
    for (int t = 0; t < n; t++) {
      double angle = TWO_PI * (double)k * (double)t / (double)n;
      sum += in[t] * (direction == FFT_FORWARD ? cexp(-I * angle)
                                               : cexp(I * angle));
    }
    out[k] = (direction == FFT_INVERSE) ? sum / n : sum;
  }
}

static void dft_2d_reference(const complex_t *in, complex_t *out, int rows,
                             int cols, fft_direction direction) {
  double norm = (direction == FFT_INVERSE) ? 1.0 / ((double)rows * cols) : 1.0;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      complex_t sum = 0.0;
      for (int u = 0; u < rows; u++) {
        for (int v = 0; v < cols; v++) {
          double angle = TWO_PI * ((double)r * u / rows + (double)c * v / cols);
          complex_t phase = direction == FFT_FORWARD ? cexp(-I * angle)
                                                     : cexp(I * angle);
          sum += in[u * cols + v] * phase;
        }
      }
      out[r * cols + c] = sum * norm;
    }
  }
}

static void test_initialization(void) {
  printf("Initialization\n");
  check(fft_gpu_available() == 1, "fft_gpu_available reports the Metal device");
  check(fft_gpu_init(FFT_GPU_METAL) == 0, "fft_gpu_init succeeds");
  check(fft_gpu_get_backend() == FFT_GPU_METAL, "backend is FFT_GPU_METAL");

  const char *name = fft_gpu_get_device_name();
  check(name != NULL && name[0] != '\0', "device name is not empty");
  printf("  device: %s\n", name != NULL ? name : "(null)");

  size_t total = 0, available = 0;
  fft_gpu_get_memory_info(&total, &available);
  check(total > 0, "recommended working set is reported");
  check(available > 0 && available <= total, "available memory is bounded");
}

static void test_host_device_copy(void) {
  printf("Host/device copy\n");
  enum { N = 5 };
  complex_t in[N], out[N];
  for (int i = 0; i < N; i++)
    in[i] = (double)i - 2.0 + I * (1.0 - (double)i);

  fft_gpu_memory_t mem = fft_gpu_alloc(N);
  check(mem != NULL, "fft_gpu_alloc returns a handle");
  if (!mem) return;

  fft_gpu_copy_h2d(mem, in, N);
  for (int i = 0; i < N; i++) out[i] = 0.0;
  fft_gpu_copy_d2h(out, mem, N);

  double err = max_difference(in, out, N);
  check_relative("copy round-trip preserves values", err, max_magnitude(in, N),
                 1e-6);
  fft_gpu_free(mem);
}

static void test_1d_forward(void) {
  printf("1-D forward\n");
  static const int sizes[] = {1, 8, 16, 100, 1000};

  for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
    int n = sizes[s];
    complex_t *in = malloc((size_t)n * sizeof(complex_t));
    complex_t *got = malloc((size_t)n * sizeof(complex_t));
    complex_t *want = malloc((size_t)n * sizeof(complex_t));
    if (!in || !got || !want) {
      check(0, "forward test allocation");
      free(in);
      free(got);
      free(want);
      return;
    }

    fill_signal(in, n);
    dft_1d_reference(in, want, n, FFT_FORWARD);
    int rc = fft_gpu_dft_1d(in, got, n, FFT_FORWARD);

    char what[96];
    snprintf(what, sizeof(what), "forward n=%d returns 0", n);
    check(rc == 0, what);
    snprintf(what, sizeof(what), "forward n=%d matches direct DFT", n);
    check_relative(what, max_difference(got, want, n), max_magnitude(want, n),
                   1e-4);

    free(in);
    free(got);
    free(want);
  }
}

static void test_1d_inverse(void) {
  printf("1-D inverse\n");
  static const int sizes[] = {8, 100};

  for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
    int n = sizes[s];
    complex_t *in = malloc((size_t)n * sizeof(complex_t));
    complex_t *got = malloc((size_t)n * sizeof(complex_t));
    complex_t *want = malloc((size_t)n * sizeof(complex_t));
    if (!in || !got || !want) {
      check(0, "inverse test allocation");
      free(in);
      free(got);
      free(want);
      return;
    }

    fill_signal(in, n);
    dft_1d_reference(in, want, n, FFT_INVERSE);
    int rc = fft_gpu_dft_1d(in, got, n, FFT_INVERSE);

    char what[96];
    snprintf(what, sizeof(what), "inverse n=%d returns 0", n);
    check(rc == 0, what);
    snprintf(what, sizeof(what), "inverse n=%d matches direct inverse DFT", n);
    check_relative(what, max_difference(got, want, n), max_magnitude(want, n),
                   1e-4);

    free(in);
    free(got);
    free(want);
  }
}

static void test_1d_roundtrip(void) {
  printf("1-D round-trip\n");
  enum { N = 64 };
  complex_t in[N], spectrum[N], back[N];
  fill_signal(in, N);

  check(fft_gpu_dft_1d(in, spectrum, N, FFT_FORWARD) == 0,
        "forward returns 0");
  check(fft_gpu_dft_1d(spectrum, back, N, FFT_INVERSE) == 0,
        "inverse returns 0");
  check_relative("forward then inverse reconstructs the input",
                 max_difference(in, back, N), max_magnitude(in, N), 1e-4);
}

static void test_batch(void) {
  printf("Batched 1-D\n");
  enum { N = 64, BATCH = 4 };
  complex_t in[N * BATCH], got[N * BATCH], want[N * BATCH];

  fill_signal(in, N * BATCH);
  for (int b = 0; b < BATCH; b++)
    dft_1d_reference(in + b * N, want + b * N, N, FFT_FORWARD);

  int rc = fft_gpu_dft_1d_batch(in, got, N, BATCH, FFT_FORWARD);
  check(rc == 0, "batched forward returns 0");
  check_relative("batched forward matches direct DFT",
                 max_difference(got, want, N * BATCH),
                 max_magnitude(want, N * BATCH), 1e-4);
}

static void test_plan_and_execute(void) {
  printf("Plan execute\n");
  enum { N = 64 };
  complex_t in[N], out[N], want[N];
  fill_signal(in, N);
  dft_1d_reference(in, want, N, FFT_FORWARD);

  fft_gpu_plan_t plan = fft_gpu_plan_1d(N, 1, FFT_FORWARD);
  check(plan != NULL, "fft_gpu_plan_1d returns a plan");
  if (!plan) return;

  fft_gpu_memory_t mem = fft_gpu_alloc(N);
  check(mem != NULL, "plan test allocation");
  if (!mem) {
    fft_gpu_destroy_plan(plan);
    return;
  }

  fft_gpu_copy_h2d(mem, in, N);
  check(fft_gpu_execute(plan, mem, mem) == 0, "in-place execute returns 0");
  fft_gpu_copy_d2h(out, mem, N);
  check_relative("in-place execute matches direct DFT",
                 max_difference(out, want, N), max_magnitude(want, N), 1e-4);

  fft_gpu_free(mem);
  fft_gpu_destroy_plan(plan);
}

static void test_plan_reuse(void) {
  printf("Plan reuse\n");
  enum { N = 32, ROUNDS = 3 };
  fft_gpu_plan_t plan = fft_gpu_plan_1d(N, 1, FFT_FORWARD);
  check(plan != NULL, "reuse plan is created");
  if (!plan) return;

  fft_gpu_memory_t src = fft_gpu_alloc(N);
  fft_gpu_memory_t dst = fft_gpu_alloc(N);
  check(src != NULL && dst != NULL, "reuse buffers allocated");
  if (!src || !dst) {
    fft_gpu_free(src);
    fft_gpu_free(dst);
    fft_gpu_destroy_plan(plan);
    return;
  }

  int ok = 1;
  for (int round = 0; round < ROUNDS; round++) {
    complex_t in[N], out[N], want[N];
    fill_signal(in, N);
    dft_1d_reference(in, want, N, FFT_FORWARD);
    fft_gpu_copy_h2d(src, in, N);
    if (fft_gpu_execute(plan, src, dst) != 0) ok = 0;
    fft_gpu_copy_d2h(out, dst, N);
    double err = max_difference(out, want, N);
    double scale = max_magnitude(want, N);
    if (err > 1e-4 * (scale < 1.0 ? 1.0 : scale)) ok = 0;
  }
  check(ok, "one plan produces correct results across three executions");

  fft_gpu_free(src);
  fft_gpu_free(dst);
  fft_gpu_destroy_plan(plan);
}

static void test_plan_uses_its_direction(void) {
  printf("Plan direction\n");
  enum { N = 32 };
  complex_t in[N], out[N], want[N];
  fill_signal(in, N);
  dft_1d_reference(in, want, N, FFT_INVERSE);

  fft_gpu_plan_t plan = fft_gpu_plan_1d(N, 1, FFT_INVERSE);
  check(plan != NULL, "inverse plan is created");
  if (!plan) return;

  fft_gpu_memory_t src = fft_gpu_alloc(N);
  fft_gpu_memory_t dst = fft_gpu_alloc(N);
  check(src != NULL && dst != NULL, "inverse plan buffers allocated");
  if (!src || !dst) {
    fft_gpu_free(src);
    fft_gpu_free(dst);
    fft_gpu_destroy_plan(plan);
    return;
  }

  fft_gpu_copy_h2d(src, in, N);
  check(fft_gpu_execute(plan, src, dst) == 0, "inverse execute returns 0");
  fft_gpu_copy_d2h(out, dst, N);
  check_relative("execute honours the plan's inverse direction",
                 max_difference(out, want, N), max_magnitude(want, N), 1e-4);

  fft_gpu_free(src);
  fft_gpu_free(dst);
  fft_gpu_destroy_plan(plan);
}

static void test_2d(void) {
  printf("2-D\n");
  enum { ROWS = 8, COLS = 16 };
  complex_t in[ROWS * COLS], got[ROWS * COLS], want[ROWS * COLS];

  fill_signal(in, ROWS * COLS);
  dft_2d_reference(in, want, ROWS, COLS, FFT_FORWARD);

  int rc = fft_gpu_dft_2d(in, got, ROWS, COLS, FFT_FORWARD);
  check(rc == 0, "fft_gpu_dft_2d forward returns 0");
  check_relative("2-D forward matches direct DFT",
                 max_difference(got, want, ROWS * COLS),
                 max_magnitude(want, ROWS * COLS), 1e-4);

  dft_2d_reference(in, want, ROWS, COLS, FFT_INVERSE);
  rc = fft_gpu_dft_2d(in, got, ROWS, COLS, FFT_INVERSE);
  check(rc == 0, "fft_gpu_dft_2d inverse returns 0");
  check_relative("2-D inverse matches direct inverse DFT",
                 max_difference(got, want, ROWS * COLS),
                 max_magnitude(want, ROWS * COLS), 1e-4);

  complex_t back[ROWS * COLS], spectrum[ROWS * COLS];
  check(fft_gpu_dft_2d(in, spectrum, ROWS, COLS, FFT_FORWARD) == 0,
        "2-D forward for round-trip");
  check(fft_gpu_dft_2d(spectrum, back, ROWS, COLS, FFT_INVERSE) == 0,
        "2-D inverse for round-trip");
  check_relative("2-D round-trip reconstructs the input",
                 max_difference(in, back, ROWS * COLS),
                 max_magnitude(in, ROWS * COLS), 1e-4);

  fft_gpu_plan_t plan = fft_gpu_plan_2d(ROWS, COLS, FFT_FORWARD);
  check(plan != NULL, "fft_gpu_plan_2d returns a plan");
  if (!plan) return;

  fft_gpu_memory_t mem = fft_gpu_alloc(ROWS * COLS);
  check(mem != NULL, "2-D plan allocation");
  if (mem) {
    fft_gpu_copy_h2d(mem, in, ROWS * COLS);
    check(fft_gpu_execute(plan, mem, mem) == 0, "2-D in-place execute returns 0");
    fft_gpu_copy_d2h(got, mem, ROWS * COLS);
    dft_2d_reference(in, want, ROWS, COLS, FFT_FORWARD);
    check_relative("2-D in-place plan matches direct DFT",
                   max_difference(got, want, ROWS * COLS),
                   max_magnitude(want, ROWS * COLS), 1e-4);
    fft_gpu_free(mem);
  }
  fft_gpu_destroy_plan(plan);
}

static void test_2d_non_power_of_two(void) {
  printf("2-D non-power-of-two\n");
  enum { ROWS = 5, COLS = 7 };
  complex_t in[ROWS * COLS], got[ROWS * COLS], want[ROWS * COLS];

  fill_signal(in, ROWS * COLS);
  int rc = fft_gpu_dft_2d(in, got, ROWS, COLS, FFT_FORWARD);
  dft_2d_reference(in, want, ROWS, COLS, FFT_FORWARD);
  check(rc == 0, "5x7 forward returns 0");
  check_relative("5x7 forward matches direct DFT",
                 max_difference(got, want, ROWS * COLS),
                 max_magnitude(want, ROWS * COLS), 1e-4);
}

static void test_planner_gpu(void) {
  printf("Planner integration\n");
  enum { N = 64 };
  complex_t in[N], out[N], want[N];
  fill_signal(in, N);
  dft_1d_reference(in, want, N, FFT_FORWARD);

  unsigned caps = fft_get_hardware_capabilities();
  check((caps & FFT_HW_GPU_MPS) != 0, "planner detects the Metal device");

  fft_plan_t plan = fft_plan_dft_1d(N, in, out, -1, FFT_PREFER_GPU);
  check(plan != NULL, "FFT_PREFER_GPU creates a plan");
  if (!plan) return;

  const char *name = fft_plan_algorithm_name(plan);
  check(name != NULL && strcmp(name, "GPU (Metal)") == 0,
        "planner selects GPU (Metal)");
  check(fft_execute(plan) == FFT_OK, "planner GPU execute returns FFT_OK");
  check_relative("planner GPU result matches direct DFT",
                 max_difference(out, want, N), max_magnitude(want, N), 1e-4);
  fft_destroy_plan(plan);
}

static void test_planner_real_gpu(void) {
  printf("Planner real-input GPU\n");
  enum { N = 64, HALF = N / 2 + 1 };
  double rin[N];
  complex_t spectrum[HALF];
  complex_t full[N], full_spectrum[N], want[HALF];

  for (int i = 0; i < N; i++)
    rin[i] = next_value();
  for (int i = 0; i < N; i++)
    full[i] = rin[i] + 0.0 * I;
  dft_1d_reference(full, full_spectrum, N, FFT_FORWARD);
  memcpy(want, full_spectrum, sizeof(want));

  fft_plan_t plan = fft_plan_r2c_1d(N, rin, spectrum, FFT_PREFER_GPU);
  check(plan != NULL, "FFT_PREFER_GPU r2c plan is created");
  if (!plan) return;

  check(fft_execute(plan) == FFT_OK, "r2c execute returns FFT_OK");
  check_relative("r2c result matches direct DFT bins",
                 max_difference(spectrum, want, HALF),
                 max_magnitude(want, HALF), 1e-4);
  fft_destroy_plan(plan);
}

static void test_planner_2d_gpu(void) {
  printf("Planner 2-D GPU\n");
  enum { ROWS = 8, COLS = 8 };
  complex_t in[ROWS * COLS], out[ROWS * COLS], want[ROWS * COLS];

  fill_signal(in, ROWS * COLS);
  dft_2d_reference(in, want, ROWS, COLS, FFT_FORWARD);

  fft_plan_t plan = fft_plan_dft_2d(ROWS, COLS, in, out, -1, FFT_PREFER_GPU);
  check(plan != NULL, "FFT_PREFER_GPU 2-D plan is created");
  if (!plan) return;

  check(fft_execute(plan) == FFT_OK, "2-D planner execute returns FFT_OK");
  check_relative("2-D planner result matches direct DFT",
                 max_difference(out, want, ROWS * COLS),
                 max_magnitude(want, ROWS * COLS), 1e-4);
  fft_destroy_plan(plan);
}

static void test_invalid_arguments(void) {
  printf("Invalid arguments\n");
  complex_t x[4] = {1.0, 2.0, 3.0, 4.0};

  check(fft_gpu_plan_1d(0, 1, FFT_FORWARD) == NULL,
        "zero-length 1-D plan is rejected");
  check(fft_gpu_plan_1d(8, 0, FFT_FORWARD) == NULL,
        "zero-batch 1-D plan is rejected");
  check(fft_gpu_plan_2d(0, 4, FFT_FORWARD) == NULL,
        "zero-row 2-D plan is rejected");
  check(fft_gpu_dft_1d(NULL, x, 4, FFT_FORWARD) == -1,
        "null 1-D input is rejected");
  check(fft_gpu_dft_1d(x, x, 0, FFT_FORWARD) == -1,
        "zero-length 1-D transform is rejected");
  check(fft_gpu_dft_2d(x, x, -1, 4, FFT_FORWARD) == -1,
        "negative 2-D size is rejected");
  check(fft_gpu_alloc(0) == NULL, "zero-length allocation is rejected");
}

int main(void) {
  if (!fft_gpu_available()) {
    printf("Metal backend not available; skipping GPU checks.\n");
    return 77;
  }

  test_initialization();
  test_host_device_copy();
  test_1d_forward();
  test_1d_inverse();
  test_1d_roundtrip();
  test_batch();
  test_plan_and_execute();
  test_plan_reuse();
  test_plan_uses_its_direction();
  test_2d();
  test_2d_non_power_of_two();
  test_planner_gpu();
  test_planner_real_gpu();
  test_planner_2d_gpu();
  test_invalid_arguments();

  fft_gpu_cleanup();
  check(fft_gpu_get_backend() == FFT_GPU_NONE, "cleanup clears the backend");

  printf("%d/%d checks passed\n", g_checks - g_failures, g_checks);
  if (g_failures == 0) {
    printf("ALL CHECKS PASSED\n");
    return 0;
  }
  printf("%d CHECKS FAILED\n", g_failures);
  return 1;
}
