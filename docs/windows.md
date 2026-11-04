# Windows Support

How to build and test this library on Windows, and what is actually verified
there. Current release: **3.5.0 (4 November 2026)**.

## What works

| Toolchain | Status |
| --------- | ------ |
| MinGW-w64 GCC (MSYS2 or Strawberry) | **Verified locally**: 3.4.0 with GCC 13.2.0, 1406/1406 + 54/54 checks, and 1415/1415 + 17712 with `-DFFT_ENABLE_OPENMP=ON`; the CI lane uses GCC 16.2.0 and passes CTest 2/2 with benchmarks 66/0/54 |
| MinGW-w64 clang (MSYS2 UCRT64) | **Verified in CI for 3.4.0** with Clang 23.1.4, CTest 2/2, benchmarks 66/0/54 and no compiler warnings |
| clang-cl | Expected to work; not verified |
| MSVC `cl.exe` | **Cannot build this library** |
| CUDA backend | **Not available natively**; verified under WSL on a GTX 1050 (CUDA 11.5) |

The 3.4.0 Windows run executes both test binaries, the benchmarks and all nine
demo programs (eight applications plus `gpu_benchmark`) with zero warnings.
Release 3.5.0 changes documentation and reporting only, so those results still
describe the code. Windows application control on the verification machine
blocked three of the optimized demo binaries by hash; the full demo sweep was
run from a Release build with `-g` added, which changes the binaries but not
the code paths. The tests and benchmarks ran from the plain Release build. The
3.3.0 CI rows are retained for history. CUDA cannot be built natively: nvcc on
Windows requires MSVC as its host compiler, and MSVC cannot build the library's
C11 `_Complex` kernels. Use WSL or Linux instead.

### Why `cl.exe` is excluded

Every kernel is written against C11 complex arithmetic:

```c
x[t] = x[t] + temp;      /* + on _Complex  */
w *= w_m;                /* *= on _Complex  */
x[i] /= n;               /* /= on _Complex  */
```

That needs the `_Complex` keyword plus the operators the standard defines for
it. Microsoft's C compiler has no `_Complex` keyword, and its own documentation
for `complex.h` states plainly that "the complex math library functions are NOT
overloaded with the +, -, *, and / operators" because complex numbers are not a
native type there — it offers `_Dcomplex` structures and function forms such as
`_Cmulcc` instead.

Making MSVC work would mean rewriting all the kernels against those helper
calls. Until that is done, `CMakeLists.txt` fails at configure time with this
explanation, and `include/fft_common.h` repeats it as an `#error`, rather than
producing hundreds of diagnostics halfway through the build.

## Build

```cmd
cmake -S . -B build
cmake --build build
.\build\bin\test_all.exe
ctest --test-dir build
```

CMake picks the compiler the same way it does everywhere else. To select MinGW
explicitly, run from an MSYS2/MinGW shell, or pass `-DCMAKE_C_COMPILER=path/to/gcc.exe`.

## Timing on Windows

`clock()` advances about once a millisecond here, so a single small transform
measures as 0.000 ms and any number derived from one looks like a timing of zero
rather than a fast function. Both the benchmark helpers and the planner's
`FFT_MEASURE` sweep work around this by growing the batch until the clock can
resolve it and reporting per-transform time. If you write your own timing on
Windows, do the same, or use `QueryPerformanceCounter`.

Consequence for the numbers in this repository: measured differences below about
10% between two of these kernels are not distinguishable from repeated runs of
one kernel. [Performance Guide](performance.md) records the measurements and
this limit together.

## Tests and demos

```cmd
.\build\bin\test_all.exe        :: correctness against an independent direct DFT
.\build\bin\test_utils.exe     :: file I/O and utility helpers
.\build\bin\benchmark_all.exe  :: timings for every implementation
.\build\bin\demo_features.exe  :: planning API and algorithm selection
.\build\bin\audio_spectrum.exe :: application example
```

## GPU support

Native Windows builds are CPU-only. Every `fft_gpu_*` symbol resolves to
`utils/fft_gpu_stubs.c`, so `fft_gpu_available()` returns 0 and
`FFT_PREFER_GPU` falls back to a CPU algorithm. `FFT_ENABLE_METAL` prints a
warning on Windows and leaves the stubs in place.

### CUDA through WSL

The CUDA backend builds and runs in WSL against the Windows display driver.
CUDA 11.5.119 was verified on a GTX 1050:

```bash
cmake -S . -B build-cuda -G Ninja -DFFT_ENABLE_CUDA=ON
cmake --build build-cuda
ctest --test-dir build-cuda --output-on-failure
```

The run passed 1406/1406 `test_all` checks, 54/54 utility checks and 74/74
`test_gpu_cuda` checks, with CTest reporting 3/3. `gpu_benchmark` measured
0.22x CPU speed at 1024, 0.67x at 4096, 3.56x at 65536 and 13.41x at 1048576;
small transforms are dominated by transfer and Pascal's FP64 rate. The CUDA
device path is double precision.

CUDA 13 dropped Pascal (sm_61), so this GPU needs CUDA 12.x or older. The
installed Windows toolkit (13.1) cannot target it; run nvcc inside WSL instead.
`FFT_ENABLE_GPU` remains as a deprecated alias for `FFT_ENABLE_CUDA`, and
enabling either one on native Windows stops at configure time with the
explanation above.

## Troubleshooting

**`cmake` not found** — install CMake and add it to `PATH`.

**`gcc` not found** — install MinGW-w64 and add its `bin` directory to `PATH`, or
configure with an explicit compiler as above.

**A freshly built `.exe` reports "Permission denied" when run** — this is
Windows application control (Smart App Control / Device Guard) refusing an
unknown binary, not a build failure. Waiting, or running the executable from a
signed/allowed location, resolves it.

**Linker reports multiple definitions of `main`** — each algorithm file keeps its
demo `main()` behind `#ifndef LIB_BUILD`, which the library target defines. Build
through CMake rather than compiling files together by hand.
