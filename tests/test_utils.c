#include "../include/fft_utils.h"

/* Utility-function tests. File I/O uses fixed names in the working directory
 * because the API takes a path; CTest runs this binary alone in the build
 * directory, and each test removes its files. */

static int g_checks = 0;
static int g_failures = 0;

static void check(int condition, const char *what) {
  g_checks++;
  if (!condition) {
    g_failures++;
    printf("  [FAIL] %s\n", what);
  }
}

static void check_close(double got, double want, double tol, const char *what) {
  double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
  g_checks++;
  if (!(fabs(got - want) <= tol * scale)) {
    g_failures++;
    printf("  [FAIL] %s (got %.17g, want %.17g)\n", what, got, want);
  }
}

static void write_text(const char *path, const char *text) {
  FILE *fp = fopen(path, "w");
  if (!fp) {
    printf("  [FAIL] could not create %s\n", path);
    exit(2);
  }
  fputs(text, fp);
  fclose(fp);
}

static void test_chirp_phase(void) {
  enum { N = 800 };
  complex_t signal[N];
  const double f_start = 100.0, f_end = 500.0, fs = 4000.0;

  generate_chirp_signal(signal, N, f_start, f_end, fs);

  /* The phase is the integral of the linear frequency sweep:
   * phi(t) = 2*pi*(f0*t + 0.5*rate*t^2). Multiplying the instantaneous
   * frequency by t, which reads plausibly, doubles the sweep rate. */
  double duration = (double)N / fs;
  double rate = (f_end - f_start) / duration;
  double max_err = 0.0;
  double max_imag = 0.0;
  for (int i = 0; i < N; i++) {
    double t = (double)i / fs;
    double want = sin(TWO_PI * (f_start * t + 0.5 * rate * t * t));
    double err = fabs(creal(signal[i]) - want);
    if (err > max_err) max_err = err;
    if (fabs(cimag(signal[i])) > max_imag) max_imag = fabs(cimag(signal[i]));
  }
  check(max_err < 1e-9, "chirp follows the integrated phase");
  check(max_imag < 1e-12, "chirp is real-valued");
}

static void test_windows(void) {
  complex_t signal[8];
  for (int i = 0; i < 8; i++) signal[i] = 1.0 + 0.0 * I;

  apply_window_kaiser(signal, 8, 0.0);
  double kaiser_err = 0.0;
  for (int i = 0; i < 8; i++) {
    double err = cabs(signal[i] - 1.0);
    if (err > kaiser_err) kaiser_err = err;
  }
  check(kaiser_err < 1e-12, "Kaiser beta=0 leaves the signal unchanged");

  for (int i = 0; i < 8; i++) signal[i] = 1.0 + 0.0 * I;
  apply_window_tukey(signal, 8, 0.0);
  double tukey_err = 0.0;
  for (int i = 0; i < 8; i++) {
    double err = cabs(signal[i] - 1.0);
    if (err > tukey_err) tukey_err = err;
  }
  check(tukey_err < 1e-12, "Tukey alpha=0 leaves the signal unchanged");

  for (int i = 0; i < 8; i++) signal[i] = 1.0 + 0.0 * I;
  apply_window_tukey(signal, 8, 1.0);
  check(cabs(signal[0]) < 1e-12 && cabs(signal[7]) < 1e-12,
        "Tukey alpha=1 tapers both ends");

  for (int i = 0; i < 8; i++) signal[i] = 1.0 + 0.0 * I;
  apply_window_tukey(signal, 8, 0.5);
  check(cabs(signal[3] - 1.0) < 1e-12,
        "Tukey leaves the middle at unity");

  /* These calls used to dereference the pointer immediately. Reaching the
   * next check is the assertion. */
  apply_window_kaiser(NULL, 8, 1.0);
  apply_window_tukey(NULL, 8, 1.0);
  apply_window_kaiser(signal, 0, 1.0);
  apply_window_tukey(signal, 1, 1.0);
  check(1, "window functions ignore null or short input");
}

static void test_snr(void) {
  complex_t signal[4] = {1.0, 1.0, 1.0, 1.0};
  complex_t noise[4] = {0.1, 0.1, 0.1, 0.1};

  check_close(compute_snr(signal, noise, 4), 20.0, 1e-12,
              "10:1 amplitude ratio is 20 dB");
  check(isnan(compute_snr(NULL, noise, 4)), "SNR rejects a null signal");
  check(isnan(compute_snr(signal, NULL, 4)), "SNR rejects null noise");
  check(isnan(compute_snr(signal, noise, 0)), "SNR rejects an empty array");
}

static void test_zero_pad(void) {
  complex_t source[4] = {1.0 + 1.0 * I, 2.0 + 2.0 * I, 3.0 + 3.0 * I,
                         4.0 + 4.0 * I};

  check(zero_pad(NULL, 4, 8) == NULL, "zero_pad rejects a null signal");
  check(zero_pad(source, 8, 4) == NULL,
        "zero_pad rejects a shorter output size");

  complex_t *padded = zero_pad(source, 4, 8);
  check(padded != NULL, "zero_pad allocates");
  if (padded) {
    double copy_err = 0.0;
    for (int i = 0; i < 4; i++) {
      double err = cabs(padded[i] - source[i]);
      if (err > copy_err) copy_err = err;
    }
    check(copy_err < 1e-12, "zero_pad copies the original samples");
    check(cabs(padded[4]) == 0.0 && cabs(padded[7]) == 0.0,
          "zero_pad zeroes the tail");
    free_complex_array(padded);
  }
}

static void test_spectrum_guards(void) {
  double freqs[2] = {0.0, 0.0};
  double mags[2] = {0.0, 0.0};
  int peaks = 7;

  find_spectral_peaks(NULL, 16, 1000.0, freqs, mags, &peaks, 2);
  check(peaks == 0, "peak finder rejects a null spectrum");

  peaks = 7;
  find_spectral_peaks(NULL, 16, 1000.0, freqs, mags, NULL, 2);
  check(1, "peak finder accepts a null peak counter");

  print_spectrum_plot(NULL, 16, 1000.0, 20);
  save_gnuplot_script(NULL, NULL, NULL, NULL);

  complex_t signal[4] = {1.0, 2.0, 3.0, 4.0};
  frequency_shift(NULL, 4, 100.0, 1000.0);
  complex_t before[4];
  memcpy(before, signal, sizeof(before));
  frequency_shift(signal, 4, 0.0, 1000.0);
  check(cabs(signal[0] - before[0]) < 1e-12 &&
            cabs(signal[3] - before[3]) < 1e-12,
        "zero frequency shift leaves the signal unchanged");
  check(1, "plot, script and shift helpers reject null input");
}

static void test_file_roundtrip(void) {
  const char *path = "test_utils_roundtrip.txt";
  complex_t input[4] = {1.0 + 2.0 * I, -3.0 + 4.0 * I, 0.0, 5.0 - 6.0 * I};

  check(save_complex_array(path, input, 4) == 0, "save_complex_array writes");

  complex_t *loaded = NULL;
  int n = -1;
  int rc = load_complex_array(path, &loaded, &n);
  check(rc == 0, "load_complex_array reads what save wrote");
  check(n == 4, "round trip preserves the sample count");
  if (rc == 0 && loaded && n == 4) {
    for (int i = 0; i < 4; i++) {
      check_close(creal(loaded[i]), creal(input[i]), 1e-6,
                  "round trip real part");
      check_close(cimag(loaded[i]), cimag(input[i]), 1e-6,
                  "round trip imaginary part");
    }
  }
  free_complex_array(loaded);
  remove(path);
}

static void test_file_errors(void) {
  const char *missing = "test_utils_missing.txt";
  const char *malformed = "test_utils_malformed.txt";
  const char *truncated = "test_utils_truncated.txt";
  complex_t sample[1] = {1.0 + 0.0 * I};
  complex_t *loaded = (complex_t *)0x1;
  int n = 7;

  check(save_complex_array(NULL, sample, 1) == -1,
        "save rejects a null filename");
  check(save_complex_array("test_utils_never.txt", NULL, 1) == -1,
        "save rejects a null array");
  check(save_complex_array("test_utils_never.txt", sample, 0) == -1,
        "save rejects an empty array");

  check(load_complex_array(NULL, &loaded, &n) == -1 && loaded == NULL && n == 0,
        "load rejects a null filename and clears the outputs");
  check(load_complex_array(missing, &loaded, &n) == -1 && loaded == NULL &&
            n == 0,
        "load reports a missing file and clears the outputs");
  check(load_complex_array("test_utils_never.txt", NULL, &n) == -1,
        "load rejects a null output pointer");
  check(load_complex_array("test_utils_never.txt", &loaded, NULL) == -1,
        "load rejects a null count pointer");

  write_text(malformed, "# Size: banana\n0 1.0 2.0\n");
  loaded = (complex_t *)0x1;
  n = 7;
  check(load_complex_array(malformed, &loaded, &n) == -1 && loaded == NULL &&
            n == 0,
        "load rejects a malformed size header");
  remove(malformed);

  write_text(truncated, "# Size: 4\n0 1.0 0.0\n1 2.0 0.0\n");
  loaded = (complex_t *)0x1;
  n = 7;
  check(load_complex_array(truncated, &loaded, &n) == -1 && loaded == NULL &&
            n == 0,
        "load rejects a truncated sample list");
  remove(truncated);
}

int main(void) {
  printf("FFT Utility Test Suite\n");
  printf("======================\n");

  test_chirp_phase();
  test_windows();
  test_snr();
  test_zero_pad();
  test_spectrum_guards();
  test_file_roundtrip();
  test_file_errors();

  printf("\n=== Utility Summary ===\n");
  if (g_failures == 0) {
    printf("%d/%d checks passed (100.0%%)\n", g_checks, g_checks);
    printf("ALL UTILITY CHECKS PASSED\n");
    return 0;
  }
  printf("%d/%d checks passed (%.1f%%)\n", g_checks - g_failures, g_checks,
         100.0 * (g_checks - g_failures) / g_checks);
  printf("%d checks failed\n", g_failures);
  return 1;
}
