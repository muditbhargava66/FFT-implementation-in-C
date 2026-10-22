# Windows Support

How to build and test this library on Windows, and what is actually verified
there. Current release: **3.0.0 (23 October 2026)**.

## What works

| Toolchain | Status |
| --------- | ------ |
| MinGW-w64 GCC (MSYS2 or Strawberry) | **Verified** — 1361/1361 checks pass, no compiler warnings, on Windows 11 x86_64 with GCC 13.2.0 |
| clang / clang-cl | Expected to work; not verified here |
| MSVC `cl.exe` | **Cannot build this library** |

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

Making MSVC work would mean rewriting all eight algorithms against those helper
calls. Until that is done, `CMakeLists.txt` fails at configure time with this
explanation, and `include/fft_common.h` repeats it as an `#error`, rather than
producing hundreds of unhelpful diagnostics halfway through the build. An
earlier revision carried a `#ifdef _MSC_VER` branch that typedef'd `_Dcomplex`
and included `<complex>` — a C++ header — so it could not compile either.

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
.\build\bin\benchmark_all.exe   :: timings for every implementation
.\build\bin\demo_features.exe   :: planning API, algorithm selection, wisdom
.\build\bin\audio_spectrum.exe  :: application example
```

## GPU support

Not available on a Windows build. The CUDA sources under `gpu/` are outside the
default build, and the dispatch layer that would call them is not wired up —
every `fft_gpu_*` symbol resolves to `utils/fft_gpu_stubs.c`, so
`fft_gpu_available()` returns 0 and `FFT_PREFER_GPU` falls back to a CPU
algorithm. Earlier revisions of this page listed a `-DENABLE_CUDA=ON` option
that has never existed in `CMakeLists.txt`; the real option is
`FFT_ENABLE_GPU`, and turning it on compiles the kernels without connecting
them.

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
