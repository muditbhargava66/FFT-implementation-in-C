# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.0] - 2026-10-24

3.1.0 is a correctness and portability release. It fixes the planner and 2-D
paths found during a full macOS/Linux audit, aligns the documentation with the
code, and adds CI coverage for macOS and both Windows compilers. The CI run on
the release branch passes on Linux GCC and Clang, macOS Apple Clang, Windows
MinGW-w64 GCC 16.2.0, Windows MinGW-w64 clang 23.1.3, and the strict docs build.

### Added
- macOS CI lane (Apple Clang) running the test suite, the benchmark pass/fail
  check, and the non-interactive application demos.
- Windows CI now covers both MinGW-w64 GCC and clang through the MSYS2 UCRT64
  toolchain, and runs the application demos as well as the tests.
- Linux CI runs the application demos in addition to `test_all` and
  `benchmark_all`.
- `FFT_ERR_UNSUPPORTED` for an operation a plan cannot perform. It is used by
  `fft_execute_dft()` when handed a real-valued plan.

### Fixed
- 2-D execution allocated its transpose buffer after the row pass. An
  allocation failure returned `FFT_ERR_ALLOC` with the output already partly
  transformed; the buffer is now allocated first.
- `load_complex_array()` in `utils/fft_utils.c` did not check the allocation
  result and ignored read errors. It now fails cleanly, including when the
  stream cannot be rewound.
- The filtering demo's `design_fir_filter()` divides by `taps - 1` when it
  builds the Hamming window and did not check its allocations. It now rejects
  fewer than two taps and fails cleanly. The realtime analyzer also checks its
  allocations and zero-initializes the spectrum buffers used for peak
  detection.
- 2-D plans selected the column-pass algorithm at execute time by consulting
  the global wisdom cache again. Replacing the cache between planning and
  execution could therefore change the plan, and under `FFT_WISDOM_ONLY` could
  make an existing plan fail mid-transform. Both passes are now chosen and
  validated at planning time. The planner also rejects grids whose element
  count overflows `int`.
- `fft_execute_dft()` accepted real-to-complex and complex-to-real plans and
  silently transformed the buffers the plan was created with. It now returns
  `FFT_ERR_UNSUPPORTED` for real-valued plans; the documented contract covers
  complex 1-D and 2-D plans.
- The `FFT_MEASURE` timing sweep calibrated its batch to 60 raw `clock()` ticks.
  That is 60 ms on Windows but 60 microseconds on POSIX, where `CLOCKS_PER_SEC`
  is 1000000, so the measurement was noise on Linux and macOS. The target is
  now 20 ms expressed as `0.020 * CLOCKS_PER_SEC`.
- `fft_plan_dft_2d()` wrote a diagnostic to stderr from inside the library. It
  returns `NULL` for a non-power-of-two grid without printing.
- The 2-D regression test used an 8x8 plan and only checked that execution
  succeeded. It now plans a 4x8 transform under `FFT_WISDOM_ONLY`, clears the
  cache, and compares every output bin with a direct 2-D DFT, which catches
  execution re-reading wisdom for the column pass.

### Removed
- `fft_plan_with_nthreads()` and the unused `g_num_threads` global. The
  function stored a count that no kernel read, and removing it also removed the
  only warning in a debug build.
- `FFT_HW_CPU_AVX512` from the hardware capability flags. The detection code
  never set it.

### Changed
- Version bumped to 3.1.0 in CMake, `fft_version()`, the documentation, the
  issue template and the release badge.
- Documentation corrected against the code: removed the deleted `CHECK_NULL`
  example, documented the 2-D power-of-two restriction and the
  `fft_execute_dft()` limitation, replaced code samples that called functions
  which do not exist, described the demos as they are actually implemented
  (frequency-sampling FIR, no per-hop plan, four filter types), removed the
  stale standalone `gcc` command, and fixed the Bluestein/mixed-radix ratios in
  the README.
- The documented error contract now says that a failed 1-D transform may leave
  the output buffer partly written. The old wording promised every caller
  buffer was untouched, which `fft_execute()` did not provide.
- The CMake minimum is 3.20. The project uses `CONFIGURE_DEPENDS`,
  `add_compile_definitions`, `CMAKE_CUDA_ARCHITECTURES`, and the documented
  `ctest --test-dir` command, none of which work on 3.10.
- `fft_plan_c2r_1d()` now sets `FFT_REAL_OUTPUT` in the plan flags, matching
  the documented behavior and `fft_plan_r2c_1d()`.
- Removed the unused `mkdocstrings` plugin and its dependencies. The project has
  no Python modules to document, the plugin required Python 3.11, and it made
  the docs build fail on Python 3.10.
- Removed the unimplemented `FILTER_CUSTOM` mode from the filtering demo and
  added the missing raised-cosine transitions for band-pass and band-stop
  filters.
- Removed the unused `cmake/FFTConfig.cmake.in`; the project does not install or
  export a CMake package, so the file could never be used.
- Added Linux and macOS results and the actual compiler versions to the
  platform table. These are builds that were run, not expectations.

## [3.0.0] - 2026-10-23

Release 3.0.0 adds the planning API, real and 2-D transforms, a wisdom cache,
and a reference-comparison test suite. The release was verified on Windows 11
with MinGW-w64 GCC 13.2.0: 1361/1361 checks, no compiler warnings, and 52 passed
/ 0 failed / 2 skipped benchmark rows. The same suite and samples were re-run
after release on Ubuntu 24.04 (GCC 13.3.0, Clang 18.1.3) and macOS 27 arm64
(Apple Clang 21.0.0, Homebrew GCC 16.2.0) with the same pass counts.

### Added
- **Planning API:** `fft_plan_dft_1d`, reusable plans, `fft_execute`,
  `fft_execute_dft`, `fft_destroy_plan`, and `fft_plan_algorithm_name()`.
- **Real-valued transforms:** `fft_plan_r2c_1d` and `fft_plan_c2r_1d`, with
  Hermitian-symmetric output and `n/2+1` input bins for the inverse.
- **2-D transforms:** `fft_plan_dft_2d` using a row pass, transpose, column
  pass, and transpose back. Both dimensions must be powers of two.
- **Wisdom cache:** export/import as text, with validation of imported entries
  against the transform size and the current machine. `FFT_WISDOM_ONLY` uses a
  cached choice or refuses to plan; it never measures.
- **Measured selection:** with `FFT_MEASURE`, candidate kernels are timed on a
  private scratch buffer, interleaved, and a challenger needs a repeatable win
  before it replaces the size heuristic. `FFT_CONSERVE_MEMORY` overrides
  measurement and selects in-place or low-scratch kernels.
- **Error codes:** all transforms return `fft_error_t` and leave the buffer
  untouched on failure. The library no longer terminates the host process.
- **Applications:** 2-D convolution, coherence estimation, and the existing
  audio, filtering, pitch, power-spectrum and realtime demos.
- **CI and maintenance:** Linux GCC/Clang, Windows MinGW-w64 and a strict MkDocs
  job; Dependabot; issue forms and a pull-request template; `SECURITY.md`.
- **Tests:** an independent direct-DFT reference over 30 sizes, error-code and
  untouched-buffer checks, planning and wisdom tests, r2c/c2r round trips and a
  2-D transform check. 1361 checks in total.

### Fixed
- **Correctness:** `bit_reverse()` produced a partial permutation for
  `log2n <= 4`, which corrupted radix-2, Bluestein and mixed-radix results at
  N = 4, 8 and 16. `mixed_radix_fft_recursive()` ignored its stride argument.
  `optimized_dft_real()` left the middle bin uninitialised for odd lengths.
  `goertzel_single_bin()` returned the wrong phase.
  `frequency_to_note()` indexed the note table out of bounds below C0.
- **Memory safety:** fixed an over-read of the filter in the convolution demo,
  an over-read of the caller's spectrum in the c2r path, a buffer overflow in
  `fft_export_wisdom_to_string()`, unchecked allocations in the convolution,
  correlation and recursive paths, and a null dereference in
  `split_radix_fft()` and `mixed_radix_fft()`.
- **Edge cases:** `log2_int()` looped forever on negative input,
  `next_power_of_two(0)` returned 0, `allocate_complex_array(n <= 0)` requested
  a huge allocation, window functions produced NaN for `n <= 1`, and Bluestein
  overflowed its signed scratch size above `2^29`.
- **Planner:** a failed algorithm choice could return a non-NULL plan that
  failed every execute. `FFT_MEASURE`, `FFT_PATIENT`, `FFT_EXHAUSTIVE` and
  `FFT_WISDOM_ONLY` shared consecutive values in an OR-combined bitmask, so a
  test for one matched another. `FFT_CONSERVE_MEMORY` was ignored at
  power-of-two sizes and overridden by `FFT_MEASURE`. Imported wisdom entries
  are now validated before use.
- **Hardware detection:** AVX and AVX2 were reported from CPUID alone, without
  checking that the OS saves YMM state. Detection now checks the `XCR0` bits
  through `XGETBV` and requires CPUID leaf 7 for AVX2.
- **Build:** Linux consumers failed to link without `libm`; the library target
  now propagates it. Clang's `-Wunused-function` rejected a demo helper; the
  helper is compiled only into its demo. The feature-test macros for
  `posix_memalign` and `sysconf` are set in `fft_common.h` so strict C11/C99
  builds see the declarations.
- **Test harness:** `test_all` returned 0 even when checks failed. The benchmark
  recomputed the O(n^2) reference once per implementation, reported `0.000 ms`
  for fast transforms, and timed `Inf`/`NaN` data because it re-transformed one
  buffer without restoring it. The inverse timing was also seeded from a partial
  spectrum. All three are fixed; benchmark rows now use batched timings and a
  restored seed.
- **Applications:** `realtime_analyzer` crashed intermittently when built with
  `-march=native` on MinGW-w64 GCC 13.2.0 because of a misaligned 32-byte AVX
  store. `-march=native` is now opt-in. `demo_features` prints the algorithm the
  planner actually selected instead of a hand-written copy of the rules.
  Window functions were defined three times; the duplicate application copies
  are gone. `utils/fft_utils.c` no longer defines an unguarded `main()`.

### Changed
- **License:** Apache 2.0 instead of MIT.
- **Build:** CMake only. The legacy `Makefile`, batch scripts, `fft/`, `bin/`
  and `lib/` directories are gone. `-march=native` is behind
  `-DFFT_NATIVE_ARCH=ON`, off by default. MSVC is refused at configure time with
  an explanation instead of failing partway through.
- **Dead code removed:** the `optimizations/` directory (SSE2/AVX, pthreads,
  OpenMP and fixed-point kernels that no header declared and nothing called),
  the duplicate `iterative_fft.c`, and the unused twiddle and bit-reversal
  tables in the plan structure.
- **Kernels:** split-radix is now a real split-radix decomposition, mixed-radix
  is a real Cooley-Tukey implementation for composite sizes, and radix-4
  accumulates twiddles per stage instead of calling `cexp()` in the butterfly.
- **Performance reporting:** absolute milliseconds were not reproducible on the
  reference machine, so the README and performance guide publish ratios against
  radix-2 DIT with the observed range, plus the mixed-radix control column.
- **Documentation:** corrected claims about SIMD, threading, fixed-point, GPU
  execution and MSVC; removed the unreproducible benchmark and GPU speedup
  tables; documented the real build targets and library name.
- **Flags removed:** `FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and
  `FFT_THREADED` were accepted and ignored. Their bit positions are left unused
  so exported wisdom keeps its meaning.

## [2.0.0] - 2025-08-30

*Historical entry, kept as released. The GPU acceleration, SIMD and OpenMP
items below were later found unreachable from the planner or unbuildable and
were removed in 3.0.0; see that entry for the details.*

### Added
- **Automatic Algorithm Selection**: New `fft_auto()` API that intelligently selects the best algorithm based on transform size and hardware capabilities
- **GPU Acceleration**: 
  - CUDA support for NVIDIA GPUs with up to 33x speedup
  - Metal Performance Shaders support for Apple Silicon (M1/M2/M3)
- **Breaking: Direction Constants**: 
  - Now use -1 for forward FFT, +1 for inverse (was enum)
- **Improved Algorithm Selection**: More intelligent selection based on factorization
- **Better Memory Management**: All algorithms now use aligned memory for better performance
- **Enhanced Documentation**: Complete rewrite with ReadTheDocs support

### Fixed
- **GitHub Issue #1**: Fixed timer_t type conflict on Linux systems by renaming custom timer type to fft_timer_t
- **Algorithm Reliability**: 
  - Fixed Radix-4 FFT implementation with proper radix-4 butterflies for optimal performance
  - Fixed Split-Radix FFT algorithm with simplified reliable implementation
  - Fixed Mixed-Radix FFT scaling issues causing reconstruction errors
  - Fixed Bluestein FFT with hybrid approach: uses fft_recursive for small convolution sizes (M≤32) ensuring correctness, and radix2_dit_fft for large sizes (M>32) for optimal performance
- **Build System**: 
  - Resolved all compiler warnings across the codebase
  - Fixed platform-specific compilation issues (pthread barriers, sysconf)
  - Made SIMD code platform-aware (x86 vs ARM64)
  - Fixed duplicate main function conflicts with LIB_BUILD guards
- **Memory Management**: 
  - Fixed memory corruption issues in Split-Radix implementation
  - Proper memory alignment for SIMD operations

### Deprecated
- Direct algorithm function calls (e.g., `radix2_dit_fft()`) - use `fft_auto()` or planning API instead
- Non-aligned memory allocation - use `fft_alloc_complex()` instead of `malloc()`

## [1.0.0] - 2025-05-16

### Added
- Initial release with core FFT algorithms:
  - Radix-2 DIT (Decimation-in-Time)
  - Radix-2 DIF (Decimation-in-Frequency)
  - Radix-4 FFT
  - Split-Radix FFT
  - Bluestein FFT (arbitrary sizes)
  - Mixed-Radix FFT
  - Recursive FFT
  - Iterative FFT
- Applications:
  - Audio spectrum analyzer
  - Digital signal filtering
  - Fast convolution
  - Power spectrum analysis
  - Image FFT processing
- Optimizations:
  - SIMD implementation (SSE/AVX)
  - Parallel FFT (OpenMP)
  - Fixed-point FFT
- Benchmark suite
- Unit tests
- Example programs
- Basic documentation

[Unreleased]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.1.0...HEAD
[3.1.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v2.0.0...v3.0.0
[2.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/releases/tag/v1.0.0
