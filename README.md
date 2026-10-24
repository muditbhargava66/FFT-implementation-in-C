<div align="center">

# FFT Implementation in C

[![Documentation](https://app.readthedocs.org/projects/fft-implementation-in-c/badge/?version=latest)](https://fft-implementation-in-c.readthedocs.io/en/latest/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue)]()
[![GPU Support](https://img.shields.io/badge/GPU-source--only-lightgrey)]()
[![Release](https://img.shields.io/badge/release-3.1.0-success)]()

</div>

A C11 FFT library with eight algorithm implementations, a planning API, and
algorithm selection based on measured cost. Builds with GCC and Clang on Linux,
macOS and Windows (MinGW-w64). MSVC's `cl.exe` is not supported.

## Contents

- [What's in the library](#whats-in-the-library)
- [Build and test](#build-and-test)
- [Usage](#usage)
- [Performance](#performance)
- [Documentation](#documentation)
- [Project layout](#project-layout)
- [Contributing](#contributing)
- [License](#license)

## What's in the library

**Algorithms:** radix-2 DIT, radix-2 DIF, radix-4, split-radix, Bluestein,
mixed-radix, a recursive radix-2 implementation, and naive/optimized DFT
references.

**Planning API:** one-shot `fft_auto()`, reusable plans for 1-D complex,
real-to-complex, complex-to-real and 2-D transforms, a wisdom cache, and
`fft_plan_algorithm_name()` to see what a plan selected.

**Error handling:** every transform returns `fft_error_t`. On failure, do not
use the output buffer: a 1-D transform may already have copied or partly
transformed it. Out-of-place calls leave the input buffer unchanged, and the
library never calls `exit()`.

What it does not do:

- No SIMD, threaded or fixed-point kernels.
- No GPU execution. The CUDA and Metal sources are in `gpu/` but outside the
  default build, and the planner resolves `fft_gpu_*` to CPU stubs.
- Not thread-safe. Algorithm selection and the wisdom cache share global state,
  so do not plan from multiple threads.
- 2-D plans accept power-of-two row and column counts only.

## Build and test

Requires CMake 3.20+ and a C11 compiler with `_Complex` support.

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The test binary is `build/bin/test_all`. It compares every implementation
against a direct DFT computed independently inside the test file, across 30
sizes from 1 to 16384. The benchmark suite is `build/bin/benchmark_all`.

To build with a specific compiler:

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=clang
```

`-march=native` is off by default; enable it with `-DFFT_NATIVE_ARCH=ON` if the
binary only needs to run on the build machine.

### Verified builds

| Platform | Compiler | Result |
| --- | --- | --- |
| Windows 11 x86_64 | MinGW-w64 GCC 13.2.0 | 3.0.0 release: 1361/1361 checks, 0 warnings; benchmarks 52 passed, 0 failed, 2 skipped |
| Windows (CI, windows-latest) | MinGW-w64 GCC 16.2.0 | 3.1.0 CI: ctest 1/1, 0 warnings; benchmarks 52/0/2; demos 8/8 |
| Windows (CI, windows-latest) | MinGW-w64 Clang 23.1.3 | 3.1.0 CI: ctest 1/1, 0 warnings; benchmarks 52/0/2; demos 8/8 |
| Ubuntu 24.04 (x86_64) | GCC 13.3.0 | 3.1.0: 1361/1361 checks, 0 warnings with `-Werror`; benchmarks 52/0/2 |
| Ubuntu 24.04 (x86_64) | Clang 18.1.3 | 3.1.0: 1361/1361 checks, 0 warnings with `-Werror` |
| macOS 27, arm64 | Apple Clang 21.0.0 | 3.1.0: 1361/1361 checks, 0 warnings; benchmarks 52/0/2 |
| macOS 27, arm64 | Homebrew GCC 16.2.0 | 3.1.0: 1361/1361 checks, 0 warnings |
| Windows / Linux / macOS | MSVC `cl.exe` | Not supported (no C11 `_Complex`) |

The Windows 11 row is the 3.0.0 release measurement. The other rows were run on
the 3.1.0 source: Ubuntu 24.04 in a container, macOS on an Apple M1 Pro, and
the two Windows rows by the 3.1.0 CI workflow on `windows-latest`. The eight
demo binaries used by CI exit 0 on all three platforms; `realtime_analyzer` is
exercised on macOS and Linux. Performance figures below are from the Windows
3.0.0 machine; they are not a cross-platform comparison.

## Usage

### One-shot transform

```c
#include <fft_auto.h>

int main(void) {
    int n = 1024;
    complex_t *signal = fft_alloc_complex(n);

    for (int i = 0; i < n; i++)
        signal[i] = sin(2 * PI * 50 * i / 1000.0);

    fft_error_t err = fft_auto(signal, signal, n, FFT_FORWARD);
    if (err != FFT_OK) {
        fprintf(stderr, "fft_auto failed: %d\n", (int)err);
        fft_free(signal);
        return 1;
    }

    fft_free(signal);
    return 0;
}
```

### Reusable plan

```c
fft_plan_t plan = fft_plan_dft_1d(n, signal, signal, -1, FFT_MEASURE);
if (!plan) return 1;

for (int i = 0; i < 1000; i++) {
    generate_signal(signal, n);
    fft_execute(plan);
    process_spectrum(signal, n);
}

printf("algorithm: %s\n", fft_plan_algorithm_name(plan));
fft_destroy_plan(plan);
```

`FFT_ESTIMATE` picks by size alone. `FFT_MEASURE` times the candidates that
accept the length and caches the winner. `FFT_WISDOM_ONLY` reuses a cached
choice or returns a NULL plan, so it never measures. `FFT_CONSERVE_MEMORY`
prefers in-place kernels and takes precedence over `FFT_MEASURE`.

### Real-valued transforms

```c
double *in = fft_alloc_real(1024);
complex_t *spectrum = fft_alloc_complex(513);   /* n/2 + 1 bins */

fft_plan_t plan = fft_plan_r2c_1d(1024, in, spectrum, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

`fft_plan_c2r_1d()` is the inverse. It expects `n/2+1` complex bins and writes
`n` real samples.

### 2-D transform

```c
int rows = 256, cols = 256;
complex_t *image = fft_alloc_complex(rows * cols);

fft_plan_t plan = fft_plan_dft_2d(rows, cols, image, image, -1, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

Rows and columns must both be powers of two. The transform is implemented as a
row pass, a transpose, a column pass, and a transpose back.

### Wisdom cache

```c
char *wisdom = fft_export_wisdom_to_string();
if (wisdom) {
    save_to_file("fft_wisdom.txt", wisdom);
    free(wisdom);
}

/* Later, in another process: */
fft_import_wisdom_from_string(saved);
fft_plan_t plan = fft_plan_dft_1d(n, in, out, -1,
                                  FFT_MEASURE | FFT_WISDOM_ONLY);
```

Imported entries are validated against the transform size and the current
machine before use. A bad entry is treated as a cache miss, not as a plan that
fails later.

## Performance

Costs below are relative to radix-2 DIT at the same `n`, so below 1.00 is
faster. They were measured on an Intel i7-9750H (Windows 11, MinGW-w64 GCC
13.2.0) with the default `-O3` build. Ratios are used because absolute times on
that laptop move by up to 2.9x between runs of the same binary; the ranges are
in [docs/performance.md](docs/performance.md).

| n | Radix-2 DIF | Radix-4 | Split-Radix | Mixed-Radix | Bluestein |
| --- | --- | --- | --- | --- | --- |
| 16 | 1.06 | 1.41 | 1.64 | 0.99 | 11.11 |
| 64 | 1.06 | 1.21 | 2.44 | 1.02 | 11.30 |
| 256 | 1.00 | 1.10 | 2.38 | 0.99 | 11.90 |
| 512 | 0.95 | **0.83** | 2.00 | 1.01 | 10.34 |
| 1024 | 0.97 | 0.97 | 2.30 | 1.04 | 11.41 |
| 4096 | 1.02 | **0.89** | 1.99 | 1.03 | 9.54 |
| 16384 | 1.00 | **0.93** | 2.86 | 1.01 | 9.68 |

Mixed-radix forwards power-of-two lengths to radix-2 DIT, so its column is the
control: it must read 1.00, and it measures 0.99-1.04. That is the resolution of
the method. Radix-4 wins from 512 up, which is where `select_algorithm()` puts
the crossover. Split-radix has the lowest multiplication count but is 1.6-2.9x
slower in wall time because of its memory traffic.

For non-power-of-two lengths Bluestein beat mixed-radix at every size tested
(97, 251, 360, 1000, 1021, 4095, 6561), by 2.1x at 360 and 241x at the prime
1021. Mixed-radix stays available for callers that pass `FFT_CONSERVE_MEMORY`,
where Bluestein's roughly 9n complex-numbers of scratch do not fit.

## Documentation

The full docs are at <https://fft-implementation-in-c.readthedocs.io/>:

- [Getting started](https://fft-implementation-in-c.readthedocs.io/en/latest/getting-started/)
- [API reference](https://fft-implementation-in-c.readthedocs.io/en/latest/api-reference/)
- [Algorithm guide](https://fft-implementation-in-c.readthedocs.io/en/latest/algorithms/)
- [Performance guide](https://fft-implementation-in-c.readthedocs.io/en/latest/performance/)

## Project layout

```
include/        public headers
algorithms/     core kernels, DFT references, planner (core/, dft/, auto/)
utils/          windows, helpers, GPU stubs
applications/   runnable demos
benchmarks/     benchmark_all and gpu_benchmark
tests/          test_all.c
gpu/            CUDA and Metal sources, not in the default build
docs/           MkDocs site
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). For security reports, use the private
reporting described in [SECURITY.md](SECURITY.md) rather than a public issue.

## Migration from v2

The one-shot and planning calls from v2 still compile:

```c
fft_auto(signal, signal, n, -1);
fft_plan_t plan = fft_plan_dft_1d(n, signal, signal, -1, FFT_MEASURE);
```

Four v2 flags were removed because nothing implemented them: `FFT_PATIENT`,
`FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and `FFT_THREADED`. If your code passed one
of them, delete it.

3.1 removes `fft_plan_with_nthreads()`, which only stored a thread count that no
kernel read, and `FFT_HW_CPU_AVX512`, which the detection code never reported.
`fft_execute_dft()` now returns `FFT_ERR_UNSUPPORTED` for real-valued plans
instead of transforming the buffers stored in the plan. `fft_plan_c2r_1d()`
now records its wisdom under `FFT_REAL_OUTPUT`; cached c2r choices from older
versions no longer satisfy `FFT_WISDOM_ONLY`, so re-plan with `FFT_MEASURE` if
you depend on them.

## License

Apache License 2.0. See [LICENSE](LICENSE).

---

Copyright 2024-2026 Mudit Bhargava.
