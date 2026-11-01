# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.2.0] - 2026-11-01

3.2.0 exposes the utility module that the library already compiled, fixes the
bugs in it, adds an opt-in Metal backend, and tightens CI. The CPU kernels and
the wisdom format are unchanged, so the 3.0.0 CPU timing results still describe
them; the GPU has separate one-shot numbers.

### Added
- `include/fft_utils.h`, declaring signal generation, the Kaiser and Tukey
  windows, file I/O, peak finding, SNR, plotting, zero padding and frequency
  shifting. These functions were compiled into the library without a
  declaration, so outside code could not call them.
- `tests/test_utils.c` with 40 checks, registered as a second CTest test.
- Linux ASan/UBSan CI job with leak detection, running both test binaries and
  every application demo.
- `-DFFT_ENABLE_METAL=ON` builds a working Metal backend from
  `gpu/fft_metal.m` through MPSGraph (macOS 14+). It covers 1-D, batched,
  in-place, 2-D and one-shot transforms, and `FFT_PREFER_GPU` reaches it for
  complex, r2c, c2r and 2-D plans. `tests/test_gpu_metal.c` adds 70 checks
  against a direct DFT. The device path is float32, so host doubles are
  narrowed on upload and widened on download.

### Fixed
- The Metal file did not compute an FFT. It created `MPSImageConversion`, a
  color-space converter, and its convenience function called the CPU
  `radix2_dit_fft`; 2-D planning returned `NULL`. It is replaced with MPSGraph
  transforms, and `gpu/fft_gpu.c` is now compiled into the library when the
  backend is enabled.
- Hardware detection called `fft_gpu_get_backend()` without initializing the
  backend, so the GPU capability bits were never set and `FFT_PREFER_GPU` could
  not select a GPU. Detection now initializes first.
- A GPU plan that failed to allocate its device buffers used to execute as a
  silent no-op and report success. Plan creation now fails instead.
- The 2-D and r2c/c2r planner paths called `execute_1d()` with a GPU algorithm
  id, which returned `FFT_ERR_SIZE`. r2c/c2r now run through `fft_gpu_dft_1d()`,
  and the 2-D path issues one `fft_gpu_dft_2d()` call instead of one launch per
  row and column.
- `fft_gpu_execute()` ignored the plan's direction for the Metal backend, and
  the 1-D batch and 2-D entry points did not dispatch to Metal at all. It now
  returns 0 or -1 so the planner can report a failed graph run instead of
  reporting success.
- `gpu_benchmark` timed a blocking GPU call with `clock()` (processor time) and
  rejected float32 output at a 1e-6 threshold that the backend cannot meet. It
  now uses the monotonic wall clock off Windows and a 1e-3 threshold.
- `generate_chirp_signal()` multiplied the instantaneous frequency by time
  instead of integrating it, which swept at twice the requested rate.
- `apply_window_kaiser()` divided by zero for a one-sample window;
  `apply_window_tukey()` did the same and mishandled `alpha` outside 0..1.
- `find_spectral_peaks()` dereferenced a failed `compute_magnitude()` result and
  divided by a zero interpolation denominator. It now clears `*num_peaks`,
  rejects invalid input, and clamps the interpolated offset to one bin.
- `compute_snr()` divided by zero for an empty array; it returns `NAN` for
  invalid arguments and keeps `INFINITY` for zero noise.
- `save_complex_array()` and `load_complex_array()` dereferenced null pointers
  and printed to stderr from inside the library. They now reject null or empty
  input, return -1 on failure, and leave `*data` as `NULL` and `*n` as 0.
  `load_complex_array()` also rejects a malformed `# Size:` header.
- `zero_pad()`, the first-pass plotting helpers, the gnuplot script writer and
  `frequency_shift()` now return without touching memory when handed null or
  invalid arguments.
- The recursive FFT's header comment claimed O(n log n) space; peak usage is
  O(n).
- The application demos check every allocation before use. Convolution
  evaluates output lengths in 64 bits and validates 2-D FFT dimensions before
  allocating; the DSP demos reject FFT lengths the kernels cannot transform.
- `realtime_analyzer` printed a hard-coded 0.1 ms FFT time; it now measures the
  transform, and its real-time factor uses the hop duration. `demo_features`
  no longer prints a CPU-versus-GPU speedup from `clock()` timings that omit
  the blocking GPU wait; it compares the two outputs instead.
- `power_spectrum` avoids `log10(0)` and rejects lengths where `2 * n`
  overflows; `fft_filtering` resets the working buffer before its timing pass.

### Changed
- Windows CI runs `realtime_analyzer`, so all nine demos are exercised there.
- The macOS CI job also builds the Metal backend with `-Werror` for C and
  Objective-C and runs `test_gpu_metal`; the test reports a skip when the
  runner has no Metal device.
- Linux and macOS CI configure with `-Werror`; the warning grep stays as a
  second check.
- Removed the `fft_gpu_set_mps_options()` and `fft_gpu_set_cuda_options()`
  declarations. Neither had a definition, and the Metal backend has no runtime
  precision or storage-mode switch to expose.
- `docs/requirements.txt` raises the floors to mkdocs 1.6.1 and
  mkdocs-material 9.8.0. Dependabot's commit prefix no longer doubles the
  `(deps)` scope.
- Removed the public `CHECK_POWER_OF_TWO` macro. It printed a message and
  called `exit()`, no library code used it, and `image_fft` was its only caller.
  Use `FFT_REQUIRE_POW2` or `is_power_of_two()` instead.

### Documentation
- Corrected the algorithm count, the number of convenience wrappers, the DIF
  output order, the Bluestein length limit, the planner's highly-composite
  rule, the forward/inverse scaling note and the recursive FFT's memory claim.
- Removed duplicated performance tables and repeated error-contract text; the
  README links to the stable documentation version.
- The README, API reference, performance guide, optimizations page, Windows
  notes and security policy now describe the opt-in Metal backend instead of
  saying no GPU path exists.
- Rewrote the core-algorithm README, which still used repeated template blocks
  and star ratings that contradicted the measured results.
- `docs/applications.md` now lists only the helpers each demo actually calls
  and says which routines are present but unused.

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
  fewer than two taps and fails cleanly; the demo also checks its test-signal
  allocations and frees them on the FIR failure path. The realtime analyzer
  checks its allocations and zero-initializes the spectrum buffers used for
  peak detection.
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
  the output buffer partly written, and the core-transform header no longer
  promises an unchanged buffer for every failure. The old wording was broader
  than the implementation.
- The CMake minimum is 3.20. The project uses `CONFIGURE_DEPENDS`,
  `add_compile_definitions`, `CMAKE_CUDA_ARCHITECTURES`, and the documented
  `ctest --test-dir` command, none of which work on 3.10.
- `fft_plan_c2r_1d()` now sets `FFT_REAL_OUTPUT` in the plan flags, matching
  the documented behavior and `fft_plan_r2c_1d()`. This changes the c2r wisdom
  key from the complex-plan key to the real-output key. C2r entries written by
  older versions, or entries for complex plans, no longer satisfy
  `FFT_WISDOM_ONLY`; re-plan with `FFT_MEASURE` if you depend on a cached c2r
  choice.
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

[Unreleased]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.2.0...HEAD
[3.2.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.1.0...v3.2.0
[3.1.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.0.0...v3.1.0
[3.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v2.0.0...v3.0.0
[2.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/releases/tag/v1.0.0
