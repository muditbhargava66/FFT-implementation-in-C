<div align="center">

# FFT Implementation in C

[![Documentation](https://app.readthedocs.org/projects/fft-implementation-in-c/badge/?version=stable)](https://fft-implementation-in-c.readthedocs.io/en/stable/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue)]()
[![GPU](https://img.shields.io/badge/GPU-Metal%20opt--in-lightgrey)]()
[![Release](https://img.shields.io/badge/release-3.2.0-success)]()

</div>

A C11 FFT library with seven FFT kernels, two DFT references, a planning API,
algorithm selection based on measured cost, and an optional Metal backend for
Apple Silicon. Builds with GCC and Clang on Linux, macOS and Windows
(MinGW-w64).

> [!NOTE]
> MSVC's `cl.exe` is not supported.

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

**Error handling:** every transform returns `fft_error_t`. A failed transform
may leave its buffer partly modified, so do not use the output. The library
never calls `exit()`. The failure contract and error codes are in the
[API reference](https://fft-implementation-in-c.readthedocs.io/en/stable/api-reference/#error-handling).

What it does not do:

- No SIMD, threaded or fixed-point kernels.
- No CUDA execution. The CUDA sources compile but are not wired to the planner.
- GPU execution outside macOS. The Metal backend is opt-in
  (`-DFFT_ENABLE_METAL=ON`, macOS 14+) and computes in float32; default builds
  resolve `fft_gpu_*` to CPU stubs.
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
sizes from 1 to 16384. `build/bin/test_utils` covers the utility functions and
file I/O. The benchmark suite is `build/bin/benchmark_all`.

To build with a specific compiler:

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=clang
```

`-march=native` is off by default; enable it with `-DFFT_NATIVE_ARCH=ON` if the
binary only needs to run on the build machine.

### Metal backend

On macOS 14 or newer, with a Metal-capable device:

```bash
cmake -S . -B build-metal -DFFT_ENABLE_METAL=ON
cmake --build build-metal
ctest --test-dir build-metal --output-on-failure
```

The backend is Objective-C with automatic reference counting, so it needs a
Clang toolchain; configure stops with an explanation if the C compiler is GCC.

This compiles `gpu/fft_gpu.c` and the MPSGraph implementation in
`gpu/fft_metal.m`, links MetalPerformanceShadersGraph, and adds
`build/bin/test_gpu_metal`. That test checks 1-D, batched, in-place, 2-D and
planner-driven transforms against a direct DFT. The GPU computes in float32, so
results are compared with a tolerance; the host-side `complex_t` stays double.
The GPU API's 2-D transform accepts any positive size; the planner's 2-D plan
still requires powers of two.

The backend also powers `FFT_PREFER_GPU` in the planner and
`build/bin/gpu_benchmark`. The CUDA sources are still not part of any working
dispatch path.

### Verified builds

| Platform | Compiler | Result |
| --- | --- | --- |
| macOS 27, arm64 | Apple Clang 21.0.0 | 3.2.0: 1361/1361 + 40/40 checks, 0 warnings with `-Werror`; benchmarks 52/0/2; demos 9/9 |
| macOS 27, arm64, `FFT_ENABLE_METAL=ON` | Apple Clang 21.0.0 | 3.2.0: 1361/1361 + 40/40 + 70/70 GPU checks, 0 warnings with `-Werror`; benchmarks 52/0/2; demos 9/9 |
| macOS (CI, macos-latest) | Apple Clang | 3.2.0 CI: 2/2 tests, 0 warnings; benchmarks 52/0/2; demos 9/9; Metal test passed |
| macOS 27, arm64 | Homebrew GCC 16.2.0 | 3.2.0: 1361/1361 + 40/40 checks, 0 warnings with `-Werror`; benchmarks 52/0/2; demos 9/9 |
| Ubuntu 24.04, arm64 (container) | GCC 13.3.0 | 3.2.0: 1361/1361 + 40/40 checks, 0 warnings with `-Werror`; benchmarks 52/0/2; demos 9/9 |
| Ubuntu 24.04, arm64 (container) | Clang 18.1.3 | 3.2.0: 1361/1361 + 40/40 checks, 0 warnings with `-Werror`; benchmarks 52/0/2; demos 9/9 |
| Ubuntu 24.04, arm64 (container) | Clang 18.1.3, ASan + UBSan | 3.2.0: 1361/1361 + 40/40 checks, no sanitizer reports; demos 9/9 |
| Windows x86_64 (cross-build) | MinGW-w64 GCC 16.2.0 | 3.2.0: full tree, including `test_utils.exe` and `realtime_analyzer.exe`, builds with `-Werror`; binaries not executed here |
| Windows (CI, windows-latest) | MinGW-w64 GCC 16.2.0 | 3.2.0 CI: 2/2 tests, 0 warnings; benchmarks 52/0/2; demos 9/9 |
| Windows (CI, windows-latest) | MinGW-w64 Clang 23.1.4 | 3.2.0 CI: 2/2 tests, 0 warnings; benchmarks 52/0/2; demos 9/9 |
| Windows 11 x86_64 | MinGW-w64 GCC 13.2.0 | 3.0.0 timing machine: 1361/1361 checks; benchmarks 52/0/2 |
| Windows / Linux / macOS | MSVC `cl.exe` | Not supported (no C11 `_Complex`) |

The 3.2.0 rows were run locally on the release branch: Ubuntu 24.04 in a
container on arm64, and macOS on an Apple M1 Pro. The Metal row uses the same
M1 Pro with `-DFFT_ENABLE_METAL=ON` and `-Werror` for both C and Objective-C.
The Windows x86_64 row is a cross-build with the same compiler version as CI;
it proves the sources compile and link, not that the binaries run. The CI rows
are the 3.2.0 results; all seven jobs passed, including the Linux ASan/UBSan
lane, the Metal build and test on macOS, and `realtime_analyzer` on Windows.
Performance figures below are from the Windows 3.0.0 machine; they are not a
cross-platform comparison.

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

The [performance guide](docs/performance.md) has the measured ratios, the
harness and the noise floor. In short, radix-4 beats radix-2 DIT from n=512 on
the measured machine, radix-2 DIF matches DIT, and split-radix is 1.6-2.9x
slower in wall time despite the lowest multiplication count. For
non-power-of-two lengths Bluestein beat mixed-radix at every size tested, by
2.1x at 360 and 241x at the prime 1021. The numbers come from a Windows 11
machine running the 3.0.0 kernels; later releases do not change the kernels, so
they still describe this code.

The Metal backend has a different profile. Each one-shot MPSGraph transform
includes buffer setup and a blocking read-back, so on an M1 Pro
`gpu_benchmark` measured GPU times of 6.9-11.0 ms across 1024 to 65536 points
against 20 us to 2.1 ms on one CPU core: slower below 65536. At 1048576 points
three runs gave 13.1-16.0 ms on the GPU against 48.0-60.1 ms on the CPU, or
3.3-3.8x. Those are single-machine, one-shot measurements, not a batched
comparison. The GPU path also computes in float32, with relative error growing
from 4.5e-6 at 1024 to 5.3e-4 at 1048576.

## Documentation

The full docs are at <https://fft-implementation-in-c.readthedocs.io/>:

- [Getting started](https://fft-implementation-in-c.readthedocs.io/en/stable/getting-started/)
- [API reference](https://fft-implementation-in-c.readthedocs.io/en/stable/api-reference/)
- [Algorithm guide](https://fft-implementation-in-c.readthedocs.io/en/stable/algorithms/)
- [Performance guide](https://fft-implementation-in-c.readthedocs.io/en/stable/performance/)

## Project layout

```
include/        public headers
algorithms/     core kernels, DFT references, planner (core/, dft/, auto/)
utils/          windows, helpers, GPU stubs
applications/   runnable demos
benchmarks/     benchmark_all and gpu_benchmark
tests/          test_all, test_utils, test_gpu_metal
gpu/            CUDA sources and the opt-in MPSGraph Metal backend
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

3.2 adds `include/fft_utils.h`, which declares the utilities that
`utils/fft_utils.c` already compiled: signal generation, the Kaiser and Tukey
windows, file I/O, peak finding, SNR, zero padding and frequency shifting.
`load_complex_array()` and `save_complex_array()` now return -1 without printing
when a file cannot be opened, `compute_snr()` returns `NAN` for invalid
arguments, and the chirp generator integrates the sweep phase correctly. No
kernel, planner or wisdom-format change; 3.1 wisdom files load unchanged. The
public `CHECK_POWER_OF_TWO` macro is removed because it called `exit()`; use
`is_power_of_two()` or `FFT_REQUIRE_POW2` instead.

## License

Apache License 2.0. See [LICENSE](LICENSE).

---

<div align="center">

Copyright 2024-2026 Mudit Bhargava.

</div>
