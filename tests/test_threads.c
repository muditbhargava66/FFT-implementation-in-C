#include "../include/fft_auto.h"
#include "../include/fft_algorithms.h"
#include "../include/fft_common.h"
#include <omp.h>
#include <stdatomic.h>

/* Exercises the process-wide planner state from several threads at once:
 * first-use hardware detection, wisdom lookup/store, import and export. Each
 * thread owns its plans and buffers, so a race can only come from the cache or
 * the hardware probe. */

static atomic_int g_checks;
static atomic_int g_failures;

static void check(int condition) {
    atomic_fetch_add_explicit(&g_checks, 1, memory_order_relaxed);
    if (!condition)
        atomic_fetch_add_explicit(&g_failures, 1, memory_order_relaxed);
}

static complex_t reference_bin(const complex_t *in, int n, int k) {
    complex_t acc = 0.0;
    for (int i = 0; i < n; i++) {
        double angle = -TWO_PI * (double)((long)i * k % (long)n) / (double)n;
        acc += in[i] * (cos(angle) + I * sin(angle));
    }
    return acc;
}

static double rel_err(complex_t a, complex_t b) {
    double scale = cabs(b);
    return cabs(a - b) / (scale > 1.0 ? scale : 1.0);
}

static void run_transform(int n, unsigned seed) {
    complex_t *in = fft_alloc_complex((size_t)n);
    complex_t *out = fft_alloc_complex((size_t)n);
    if (!in || !out) {
        check(0);
        fft_free(in);
        fft_free(out);
        return;
    }

    for (int i = 0; i < n; i++) {
        in[i] = sin(0.17 * (double)(i + seed)) + I * cos(0.11 * (double)(i + seed));
        out[i] = 0.0;
    }

    fft_plan_t plan = fft_plan_dft_1d(n, in, out, -1, FFT_ESTIMATE);
    check(plan != NULL);
    if (plan) {
        check(fft_execute(plan) == FFT_OK);
        for (int k = 0; k < n; k++) {
            complex_t want = reference_bin(in, n, k);
            check(rel_err(out[k], want) < 1e-10 * (1.0 + log2_int(n)));
        }
        fft_destroy_plan(plan);
    }

    fft_free(in);
    fft_free(out);
}

int main(void) {
    const int max_threads = omp_get_max_threads();
    if (max_threads < 2) {
        printf("Only one OpenMP thread is available; skipping concurrency test.\n");
        return 77;
    }

    /* First call from every thread: detect_hardware() must publish the mask
     * exactly once and every caller must see it. */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 64; i++)
        check(fft_get_hardware_capabilities() != 0);

    /* Planning and execution of independent plans. The sizes overlap so the
     * threads contend on the same wisdom keys. */
    static const int sizes[] = {64, 128, 256, 512};
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 64; i++)
        run_transform(sizes[i % 4], (unsigned)i * 17u + 1u);

    /* Export while other threads plan; the snapshot must stay well formed. */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 64; i++) {
        char *wisdom = fft_export_wisdom_to_string();
        check(wisdom != NULL && wisdom[0] == '#');
        free(wisdom);
    }

    /* Import replaces the cache. Plans may see the old or the new generation;
     * both are valid, and planning must not crash or return a bad plan. */
    static const char *const wisdom_text =
        "# FFT Wisdom v3.0.0 (4 entries)\n"
        "64 0 0\n128 0 0\n256 0 2\n512 0 2\n";
#pragma omp parallel for schedule(static)
    for (int i = 0; i < 32; i++) {
        if ((i & 1) == 0)
            check(fft_import_wisdom_from_string(wisdom_text) == 4);
        else
            run_transform(128, (unsigned)i + 3u);
    }

    const int checks = atomic_load_explicit(&g_checks, memory_order_relaxed);
    const int failures = atomic_load_explicit(&g_failures, memory_order_relaxed);
    printf("%d/%d threaded checks passed\n", checks - failures, checks);
    return failures == 0 ? 0 : 1;
}
