# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]


## [3.0.0] - 2026-10-23

Release 3.0.0 consolidates the v3 planning API with the correctness,
performance and build work completed for it. Verified on 2026-10-23 with
MinGW-w64 GCC 13.2.0: **1361/1361 checks**, no compiler warnings, and
**52 passed / 0 failed / 2 skipped** benchmark rows in 17.1 s.

### Fixed
- **Linux links failed without `libm`.** The CMake library target never linked
  `m`, so `test_all`, the benchmarks and the examples failed with undefined
  references to `cexp`, `cabs` and `sincos` on both Linux lanes. MinGW and
  macOS resolve those from the C runtime, which is why the local Windows
  verification missed it. `fft_lib_c` now propagates `libm` to its consumers on
  Unix.
- **Clang's warning-free gate rejected an unused demonstration helper.**
  `radix4_butterfly_demo()` is called only by the guarded demonstration
  `main()`; it is now compiled only into that demo, so the library build is
  clean under `-Wunused-function`.
- **GitHub issue #3**: build conflicts and documentation inconsistencies
  resolved, including the stale radix-4 description that claimed it used
  radix-2 butterflies for every stage.
- **Multiple-definition linker errors**: every demonstration `main()` and its
  file-local helpers are guarded by `LIB_BUILD`, which the library target
  defines.
- **Dead static helpers** (`dft_2`, `dft_5`, `count_factors`) removed from
  `mixed_radix.c` and `fft_auto.c`.
- **`factorize` and `free_factorization`** moved behind `#ifndef LIB_BUILD` in
  `mixed_radix.c`; only the demo uses them.
- **`ALGO_COUNT` cases** added to the dispatch switches to silence `-Wswitch`
  and make an unknown algorithm fail visibly.
- **`bit_reverse()` produced a wrong permutation for `log2n <= 4`.** The final
  byte-swap of the 16-bit reversal network was guarded by `if (log2n > 4)`, so
  `radix2_dit_fft` and `radix2_dif_fft` returned incorrect spectra for N = 4, 8
  and 16 (measured relative error 0.9-2.7 against a direct DFT), and Bluestein
  and mixed-radix inherited the fault at those sizes. Not previously reported.
- **`test_all` could not report failure.** `main()` returned 0 unconditionally,
  so CTest passed regardless of outcome. It now returns 1 when any check fails.
- **The test suite had no reference comparison.** No test compared an FFT
  against an independently computed DFT, so bin-ordering and twiddle-index
  defects were undetectable. Added a direct-summation reference over a checked
  bin set (every bin below N=1024, a strided sample plus edges above), with
  tolerance scaled as O(sqrt(log n)) to match backward-stable round-off.
  Coverage now includes N=1, primes (3,5,7,31,97,251,1021), non-power-of-two
  composites, 2*4^k sizes for radix-4, and N up to 16384: 1361 checks.
- **`log2_int()` looped forever on negative input** (`while (n >>= 1)` never
  terminates for a negative value under arithmetic shift). Now returns 0 for
  `n <= 1`, with regression tests.
- **`next_power_of_two(0)` returned 0**, which callers then allocated as a
  zero-length FFT. Now returns 1 and computes in unsigned to avoid signed
  overflow.
- **`allocate_complex_array()` with `n <= 0`** passed a negative count through
  `int` to `size_t`, producing a huge allocation request. Returns `NULL` now.
- **`circular_convolution()` over-read its filter argument by 64 bytes** in the
  convolution demo: `h` holds `nh = 4` elements but was copied as if it held
  `nx = 8`. The caller now zero-extends `h`.
- **`fft_execute()`'s c2r path over-read the caller's spectrum.** It copied `n`
  complex values from an input documented (in `fft_auto.h` and
  `docs/api-reference.md`) as holding `n/2+1`; the mirror bins are now
  generated from the conjugate symmetry instead of being read out of bounds.
- **`frequency_to_note()` indexed `notes[]` out of bounds** for frequencies
  below C0 (~16.35 Hz), because C's `%` keeps the sign of the dividend. The
  pitch class is now normalised and the octave floored.
- **`fft_export_wisdom_to_string()` could overflow its buffer** by accumulating
  `snprintf`'s would-be length rather than the bytes written. Truncation now
  stops the loop instead of driving `cap - pos` negative.
- **`mixed_radix_fft_recursive()` ignored its `stride` argument** despite being
  declared in `fft_algorithms.h`; strided calls transformed the wrong elements
  and overwrote `x[0..n-1]`.
- **`optimized_dft_real()` left bin `(n+1)/2` uninitialised for odd `n`.** The
  conjugate-mirror guard tested `k < n/2`, which skips the write for every odd
  length.
- **`goertzel_single_bin()` returned the wrong phase.** It computed
  `s1 - s2*e^{jw}` instead of `e^{jw}*s1 - s2`.
- **Window functions turned the signal into NaN for `n <= 1`** by dividing by
  `n-1`; a shorter-than-two window is now left untouched.
- **`apply_window_hann/hamming/blackman` were defined three times** (library
  plus two applications), leaving duplicate external symbols resolved only by
  archive-membership luck. The application copies are gone.
- **`utils/fft_utils.c` defined an unguarded `main()`**, so any program calling
  one of its helpers failed to link with "multiple definition of main". This
  was the one file missed by the v3.0.0 `LIB_BUILD` sweep.
- **`fft_radix2_sse2()` was mathematically wrong for all n >= 8**: it broadcast
  a single twiddle across four lanes that each need a different power, and
  advanced it once per group of four. Each lane now carries its own twiddle and
  the vector is advanced by w^4 per group. Verified against double-precision
  radix-2 to within 1.5e-08 (N=2) .. 9.1e-05 (N=4096), consistent with
  single-precision round-off; previously it produced absolute errors of 3-22.
- **Feature-test macros are now set in `fft_common.h`**, so `posix_memalign`
  and `sysconf` are declared under strict `-std=c11`/`-std=c99` regardless of
  compiler flags (GCC 14+ makes the resulting implicit declarations hard
  errors). This is the remaining cause of the Debian build failures in issue #4.
- **`M_PI` is defined on every platform**, not only under `_WIN32`.
- **Unchecked allocations in `fft_convolution`, `autocorrelation_fft`,
  `cross_correlation_fft` and `circular_convolution`** were dereferenced
  immediately; they now fail cleanly.
- **`fft_recursive()` dereferenced a failed allocation.** Its scratch buffer and
  the even/odd buffers inside the recursion were used without a null check.
- **`optimized_dft()` and `create_twiddle_cache()`** dereferenced three
  allocations and an `n*n` twiddle table without checking, and
  `optimized_dft()` accepted `n <= 0`.
- **`bluestein_fft()` computed `2 * n - 1` as a signed `int`.** For
  `n > 2^30` that overflows, and the negative result fed `next_power_of_two()`
  into a path where the scratch arrays were smaller than the loops indexing
  them. Lengths above `2^29` now return `FFT_ERR_RANGE`.
- **`split_radix_fft()` and `mixed_radix_fft()` dereferenced a null `x`.**
  Neither checked the buffer before `memcpy`.
- **Planning flags aliased.** `FFT_MEASURE`, `FFT_PATIENT`, `FFT_EXHAUSTIVE`
  and `FFT_WISDOM_ONLY` were the consecutive values 1, 2, 3, 4 in an
  OR-combined bitmask, so `flags & FFT_MEASURE` was also true for
  `FFT_EXHAUSTIVE`. They are distinct bits now.
- **`detect_hardware()` advertised AVX and AVX2 on CPUID alone.** Executing them
  raises `#UD` unless the operating system preserves YMM state, so SSE/AVX are
  now gated on the `XCR0` bits checked through `XGETBV`, AVX2 additionally on
  CPUID leaf 7 existing.
- **The demo in `applications/demo_features.c` printed a hand-written copy of
  the algorithm-selection rules** instead of asking the plan, and had already
  drifted: it claimed 4096 and 16384 used Split-Radix and 360 used Mixed-Radix,
  while the library chose other algorithms. It now reports
  `fft_plan_algorithm_name()`.
- **Dead plan data allocated and never read.** `fft_plan_dft_1d()` built an
  `n`-entry twiddle table (with the forward sign hardcoded, even for inverse
  plans) and an `n`-entry bit-reversal table that no kernel ever consulted,
  costing two allocations and `n` `cexp()` calls per plan.
- **`applications/realtime_analyzer` crashed on about half of all launches.**
  Built with the project's then-default `-march=native`, MinGW-w64 GCC 13.2.0
  materialised the 32-byte `analyzer_config_t` passed by value to
  `create_analyzer()` with `vmovdqa %ymm4, 48(%rsp)` — a 32-byte *aligned*
  store — into a frame it had only aligned to 16 bytes. The x64 Windows ABI
  guarantees 16, and `-mpreferred-stack-boundary` is capped at 4 on this target,
  so 32 cannot be assumed; whether the store faults depends on the mod-32 phase
  of `%rsp` when the CRT calls `main()`, which varies per process launch. An
  exception filter confirmed it: 18 of 18 captured faults were that instruction
  with `(%rsp+48) % 32 == 16`. Measured 14/25 launches SIGSEGV with
  `-march=native`, 0/40 without; `-mstackrealign` still crashed 10/25 because it
  only realigns to 16. Fixed by making `-march=native` opt-in (see Changed).
- **`benchmarks/benchmark_all` took minutes to run and reported `0.000 ms`.**
  Two separate defects. It recomputed the O(n^2) reference DFT inside
  `run_single_benchmark()`, i.e. once per implementation — 9 times per size,
  which is 9 x 14.3 s = 129 s at N=16384 before the O(n^2) `naive_dft` and
  `optimized_dft` rows. The reference is now computed once per size, and above
  128 bins it is evaluated on an evenly spaced sample, which is what makes the
  whole suite finish in 17 s instead of stalling. Separately, a fixed iteration
  count against `clock()` — which advances about once a millisecond here — made
  every fast transform print `0.000 ms` and made the scaling analysis take
  `log(0)`. Batches now grow until they span at least 20 ms, as the planner
  already did.
- **`benchmark_all`'s inverse-FFT timing fed a partially computed spectrum.**
  The round trip is now seeded with the implementation's own forward output, so
  reconstruction can be checked over all `n` samples at every size rather than
  only where the full reference happened to exist.
- **`measure_best_algorithm()` timed algorithms on data that overflowed to Inf.**
  `time_batch()` re-ran `execute_1d()` on one scratch buffer without restoring it,
  and a forward transform scales its input by `n`, so the buffer reached `Inf`
  after a few dozen reps and then `NaN`. That measured about 2x slow at n=16384
  and, more damagingly, inverted the ranking: it put radix-4 *behind* radix-2 DIT
  where a correct harness consistently puts it 6-17% ahead from 512 up. Each
  repetition now restores the buffer from a seed. The same defect was in every
  timing probe in this repository, which is where the table in
  `docs/performance.md` came from.
- **The published performance tables were absolute milliseconds, and absolute
  milliseconds are not reproducible here.** Three runs of one probe, minutes
  apart, measured radix-2 DIT at n=16384 at 0.501, 1.372 and 1.433 ms — a 2.9x
  spread for byte-identical code, from the machine's power and thermal state.
  `README.md` and `docs/performance.md` now publish cost *relative to radix-2 DIT
  at the same n*, measured within a run, with the observed min-max range per cell.
  The control column is mixed-radix, which forwards power-of-two lengths straight
  to radix-2 DIT and so must read 1.00; it measures 0.99-1.04, which bounds the
  resolution of every other figure.
- **The planners could return a non-NULL plan that could never execute.**
  `choose_algorithm()` returns `ALGO_COUNT` when it has no algorithm to offer,
  and both `fft_plan_dft_1d` and `fft_plan_dft_2d` fell straight through their
  `switch` on that value and returned the half-built plan. `execute_1d()` then
  failed every call with `FFT_ERR_SIZE`, labelled in the source as
  "sentinel — never selected". A caller checking only for `NULL`, as the API
  invites, got a plan that reported a size error on a valid size. Both planners
  now free the partial plan and return `NULL`. `fft_plan_dft_2d` checks the
  column length as well as the row length, because the column pass chooses its
  own algorithm at execute time and would otherwise fail partway through a
  transform rather than at planning.
- **`benchmark_all` dereferenced unchecked allocations and divided by a count
  that could be zero.** Its per-size loop passed `allocate_complex_array()`
  results straight to `memcpy` and the signal generator, and the summary divided
  the pass/fail percentages by `result_idx`, which is 0 if that loop breaks
  before recording anything.
- **Wisdom entries were believed without validation.** `choose_algorithm()`
  returned whatever the `(n, flags)` cache held, and that cache can be loaded
  with `fft_import_wisdom_from_string()`. An entry such as `7 0 2` (radix-4 at a
  prime length) produced a non-NULL plan that failed every execute with
  `FFT_ERR_SIZE`; an imported GPU id on a build without a GPU produced a plan
  whose execute returned `FFT_OK` *without transforming anything*, because the
  GPU dispatch is compiled out. Cache hits and import now check that the
  algorithm can transform the recorded size and is available on this machine;
  an unusable entry is skipped at import and treated as a miss at planning time,
  so `FFT_WISDOM_ONLY` fails closed instead of returning a dead plan.
- **`FFT_CONSERVE_MEMORY` was ignored at power-of-two lengths and overridden by
  `FFT_MEASURE`.** The power-of-two branch of `select_algorithm()` returned
  before consulting the flag, choosing radix-4 at n >= 512 — which digit-reverses
  through an O(n) scratch array — while in-place radix-2 DIT was available. The
  flag now selects radix-2 DIT for every power-of-two length and suppresses
  measurement, which would otherwise hand the choice back to the allocating
  kernel: at n=4096, `FFT_CONSERVE_MEMORY|FFT_MEASURE` chose radix-4 too.

### Performance
- **`split_radix_fft()` is now a real split-radix FFT.** The previous version
  was radix-2 DIT under a split-radix name and docstring; it now performs the
  one `n/2` plus two `n/4` decomposition that gives split-radix its
  multiplication count, using a bump-allocated scratch pool rather than a
  `malloc` per node.
- **`mixed_radix_fft()` is now a real Cooley-Tukey mixed-radix FFT.** Non
  power-of-two lengths previously fell through to a direct O(n^2) DFT. Measured
  at N=4095: **1513 ms -> 5.90 ms**. The 1513 ms is of the superseded code and
  predates the harness correction below, so the ratio is the meaningful part.
- **Radix-4 no longer calls `cexp()` inside the innermost loop.** It evaluated
  three transcendental functions per butterfly; twiddles are now accumulated
  multiplicatively per stage, as `radix2_dit` already did. That moved it from
  3.7-4.0x slower than radix-2 DIT to between 0.83x and 0.97x of it at every
  power-of-two length from 512 up, with the cost of its O(n) digit-reversal
  buffer still dominating below that.
- **Measured, then reported honestly.** Costs are now published as ratios
  against radix-2 DIT at the same length, measured inside a single run, because
  absolute milliseconds on this machine vary by up to 2.9x between runs of one
  binary. Radix-2 DIF is indistinguishable from DIT (0.95-1.07 at every size).
  Radix-4 is 1.41x and 1.21x slower at 16 and 64 and then 0.83x, 0.89x and 0.93x
  faster at 512, 4096 and 16384. Split-radix is reproducibly 1.6-2.9x slower
  despite the lower flop count, because its recursion moves about 2.5n elements
  per level. Bluestein costs about 10x where a radix kernel applies. See
  `docs/performance.md` for the full table with per-cell ranges.
- **`select_algorithm()` retuned against those measurements**: the radix-2 DIT /
  radix-4 crossover moved from 64 to 256, which is where the corrected ratios
  put it. Bluestein stays the choice for every non-power-of-two length, measured
  2.1x to 241x faster than mixed-radix, which degrades to O(n^2) at a prime. The
  previous rules cited a power-of-four test that no longer decides anything.
  Planning with `FFT_MEASURE` now reaches the same answer as the heuristic at
  every size checked (64 through 16384, and 360/1000/1021/4095), which is the
  evidence that the retuned threshold is right rather than merely re-fitted.

### Added
- **Repository automation and security metadata**:
  `.github/dependabot.yml` (GitHub Actions weekly, docs pip dependencies
  monthly), `.github/pull_request_template.md`, bug and feature issue forms,
  and `SECURITY.md` with private reporting guidance and an explicit scope. CI
  now builds the MkDocs site in strict mode and runs C jobs on `release/**`
  pushes as well as `main`.
- **Real-to-Complex FFT** (`fft_plan_r2c_1d`): real input is promoted at execute
  time and produces `n/2+1` Hermitian-symmetric complex bins.
- **Complex-to-Real FFT** (`fft_plan_c2r_1d`): the inverse counterpart, which
  expands the Hermitian half-spectrum and returns `n` real samples.
- **2-D FFT planning** (`fft_plan_dft_2d`): row-column decomposition with an
  internal transpose for power-of-two grids.
- **2-D FFT convolution** (`fft_convolution_2d`): zero-padded row-column FFT,
  pointwise multiplication, valid-region extraction.
- **Cross-PSD coherence** (`compute_coherence`): magnitude-squared coherence
  via Welch's method with 50% overlap and a Hann window.
- **Modern CMake build** and **CONTRIBUTING.md**; the legacy build scripts are
  gone.
- **`FFT_WISDOM_ONLY` is implemented.** It was declared, documented and ignored.
  A plan built with it now uses a cached algorithm choice or refuses, so a
  latency-bounded caller can get `NULL` instead of paying for a timing sweep it
  did not ask for. The flag is masked out of the cache key before the lookup,
  which is what makes it usable: a decision recorded under `FFT_MEASURE` is found
  by a later plan passing `FFT_MEASURE | FFT_WISDOM_ONLY`. Keying on the flags
  verbatim would make that combination always miss and the flag pointless.
- **Continuous integration** (`.github/workflows/ci.yml`): Linux under gcc and
  clang, and Windows under MinGW-w64 ucrt64 via `msys2/setup-msys2`. MSVC is not
  a lane because `CMakeLists.txt` refuses it at configure time. Each lane
  configures, builds, fails on any compiler warning, runs `ctest`, and runs
  `benchmark_all` — checking its `Failed: 0` line, because that binary exits 0
  even when rows report FAIL. The Windows lane needs the MSYS2 shell for every
  step, which is why this is two jobs rather than one matrix: `shell:` cannot be
  varied per matrix entry reliably. The gcc lanes are validated locally; the
  clang lane and both runners are not, since neither clang nor a Linux host is
  available here.
- **Error codes instead of `exit()`.** Every transform now returns `fft_error_t`
  (`FFT_OK`, `FFT_ERR_NULL`, `FFT_ERR_SIZE`, `FFT_ERR_ALLOC`, `FFT_ERR_RANGE`)
  and leaves the caller's buffer untouched when it fails. Previously a
  non-power-of-two length reached a library function via `CHECK_POWER_OF_TWO`
  and terminated the host process. The API change covers the nine transforms in
  `fft_algorithms.h`, the sixteen convenience wrappers, `fft_execute()`,
  `fft_execute_dft()` and `fft_auto()`. Discarding the value remains valid C, so
  existing callers compile unchanged.
- **`fft_plan_algorithm_name(plan)`** reports which algorithm a plan selected.
  `select_algorithm()` is static, and without an accessor the only way to show
  the choice was to re-implement the rules — which `demo_features.c` had done and
  gotten wrong.
- **The wisdom cache is now real.** `fft_export_wisdom_to_string()` and
  `fft_import_wisdom_from_string()` existed since v3.0.0 but nothing populated
  `g_wisdom` and nothing read it, so planning never consulted it and
  `FFT_MEASURE` was an empty `if` block with a `TODO`. Planning now goes through
  `choose_algorithm()`: cached answer if present, otherwise the size heuristic,
  upgraded to a measured choice under `FFT_MEASURE`.
- **`FFT_MEASURE` performs a timed sweep.** Each CPU algorithm that accepts the
  length is batch-timed against a private scratch buffer (never the caller's,
  since an FFT destroys its input), sampled min-of-three with the candidates
  interleaved, and a challenger must beat the heuristic by 15% *and* repeat that
  under a paired A/B confirmation before it is selected. Interleaving and the
  margin exist because the kernels are closer together than this machine's
  run-to-run spread; without them the choice flipped between runs.
- **Batch-calibrated timing throughout.** `clock()` ticks about once a
  millisecond on Windows, so a single small transform measured as exactly
  `0.000 ms` and the planner ranked candidates by their position in a list.
  Batches now grow until the clock can resolve them.
- **Regression tests for the new surface** (`tests/test_all.c`): every
  implementation's error codes and untouched-on-failure buffer, planning-flag
  non-aliasing, heuristic and measured plans compared against the independent
  reference DFT, a wisdom export/import round trip, rejection of a corrupt
  wisdom string, the r2c/c2r round trip, and a 2-D DC-bin check. 1361 checks,
  up from 1340.

### Changed
- **`.gitignore`** now covers MkDocs output (`site/`, `docs/site/`), in-source
  CMake artifacts, Python caches and coverage files, local `.env` files, OS
  metadata, and crash dumps.
- **License**: moved from MIT to Apache License 2.0 for patent protection.
- **`fft_execute` refactor**: an `execute_1d()` helper now backs the 1-D, r2c,
  c2r and 2-D paths, with the 2-D path transposing between row and column
  passes.
- **`fft_destroy_plan`** now frees the internal complex buffer used by r2c/c2r
  plans.
- **Window functions** use standard `creal()`/`cimag()` and `complex_t`
  assignment instead of the non-portable `.real`/`.imag` members.
- **Project layout** streamlined around `include/`, `algorithms/{core,dft,auto}`,
  `utils/`, `applications/`, `benchmarks/`, `tests/`, `gpu/` and `docs/`.
- **Documentation corrected where it contradicted the code.** `radix4.c` no
  longer claims O(1) space (it allocates an O(n) digit-reversal buffer) or that
  the implementation "uses reliable radix-2 butterflies for all stages" — a
  stale comment and the likely origin of issue #3. `split_radix.c`'s header
  comment now describes the decomposition it actually performs.
- **Unsupported claims removed from the documentation**, having been checked
  against the build:
  - *MSVC supported*: it is not, and never was — see below.
  - *GPU acceleration, with a table of CUDA speedups*: the GPU sources are
    outside the default build behind `if(FALSE AND CMAKE_CUDA_COMPILER)`, and no
    dispatch layer reaches them, so no GPU number in this repository was
    measured. The speedup table is gone; the badge and feature list say
    "source only".
  - *"100% benchmark pass rate (543/543)"* and a `performance.md` results table
    whose figures (radix-2 DIT at 0.020 ms for N=1024, split-radix at 0.079 ms
    for N=4096) disagree with the current code by 4x and 6.5x, including rows
    reading `0.000 ms`. Replaced with measured numbers, the method, the machine,
    and the noise floor.
  - *"OpenMP multi-core support"* and *"AVX-512 vectorization"*: neither is
    built or reachable; `FFT_THREADED` is accepted and ignored.
  - `-DENABLE_CUDA=ON` in `docs/windows.md`: no such option existed.
- **MSVC is now refused with an explanation rather than half-declared.**
  `fft_common.h` had an `#ifdef _MSC_VER` branch that included `<complex>` — a
  C++ header — and typedef'd `_Dcomplex`, a structure with no arithmetic
  operators in C. It could not compile either way. `CMakeLists.txt` now stops at
  configure time with the reason, and the header repeats it as an `#error`.
  Clang and clang-cl are explicitly not caught by this.
- **`docs/optimizations.md`'s manual build command now lists every source
  directory and `-DLIB_BUILD`.** The old one omitted `algorithms/dft`,
  `algorithms/auto` and `optimizations`, and would have failed to link.
- `benchmarks/gpu_benchmark.c` is built as a CMake target so it cannot silently
  drift from the GPU API again; it had been calling `fft_dft_1d()` and
  `complex_t.re/.im`, neither of which exists.
- Build documentation no longer prescribes `-ffast-math` (it reassociates
  butterfly sums and assumes no NaN/Inf, which is actively harmful for an FFT),
  references the real library name `libfft_lib_c.a` instead of `-lfft`, and
  drops instructions for the removed `Makefile`.
- **`-march=native` is no longer a default compile option.** It is behind
  `option(FFT_NATIVE_ARCH ... OFF)`, which warns when enabled. Two reasons: it
  makes binaries non-portable to other CPUs, and it triggers the MinGW-w64 GCC
  misaligned-AVX-store bug that made `realtime_analyzer` crash on about half of
  all launches (see Fixed). Every timing in `README.md` and
  `docs/performance.md` was re-measured without it, so the published numbers now
  describe what the default build actually produces.
- **The source glob uses `CONFIGURE_DEPENDS`.** Without it CMake caches the file
  list, so a stale build directory keeps trying to compile sources that were
  deleted — which is exactly what removing `optimizations/` would have hit.
- **Documentation no longer describes code that does not exist.** `README.md`,
  `docs/index.md`, `docs/optimizations.md`, `docs/performance.md`,
  `docs/algorithms.md` and `CONTRIBUTING.md` all advertised the SIMD, threading
  and fixed-point kernels as features. The AVX2 and OpenMP listings in
  `docs/performance.md` are kept but relabelled as unimplemented sketches, since
  their shape is useful to whoever picks the work up and leaving them looking
  like shipped code is how the original claim started. The `bit_reverse()`
  listings in `docs/optimizations.md` and `docs/performance.md` reproduced the
  buggy `if (log2n > 4)` version under a function name that never existed; both
  now show the real code and say what the guard cost. The macOS OpenMP
  troubleshooting section in `docs/getting-started.md` is gone, along with the
  claim that "some warnings are expected" — the build is warning-clean and CI
  fails on the first one.

### Removed
- Legacy build scripts (`Makefile`, `.bat`, `.ps1`) and the redundant `fft/`,
  `bin/` and `lib/` directories.
- **`CHECK_NULL`** from `fft_common.h`. A macro in a public header that calls
  `exit()` is not something a library should export; its last users were the
  three `compute_*` helpers, which now return `NULL`. `CHECK_POWER_OF_TWO`
  remains for standalone programs (`applications/image_fft.c`) and is documented
  as such.
- The dead `twiddles` and `bit_reverse_table` fields from `struct fft_plan`.
- `is_power_of_four()` and `is_prime()` from `fft_auto.c`; nothing selects on
  them any more. (`is_power_of_four` was listed as a v3.0.0 feature.)
- The stale `extern` re-declaration of all nine transforms in
  `benchmarks/benchmark_all.c`, which duplicated `fft_algorithms.h` and drifted.
- The MSVC `_Dcomplex` compatibility branch in `fft_common.h`.
- **The whole `optimizations/` directory**: `simd_fft.c`, `parallel_fft.c` and
  `fixed_point_fft.c`. All three were compiled into `libfft_lib_c.a` by the
  source glob, but no header declared any of their functions and nothing called
  them, so no plan could ever reach them. Verified before deleting: zero
  declarations in `include/*.h` and zero references outside their own files for
  `fft_radix2_sse2`, `fft_radix2_avx`, `detect_cpu_features`,
  `fft_radix2_parallel`, `four_step_fft`, `create_thread_pool`,
  `fft_radix2_q15`, `ifft_radix2_q15` and `normalize_block`. On Windows
  `four_step_fft` was not a four-step transform at all but
  `return radix2_dit_fft(...)`. Restoring any of this means writing kernels that
  `select_algorithm()` can choose and `tests/test_all.c` can verify, not
  reinstating the files.
- **`algorithms/core/iterative_fft.c`.** `iterative_fft_detailed()` was a second
  radix-2 DIT, algorithmically identical to `radix2_dit_fft()`, declared in no
  header and called only by that file's own demonstration `main()`. The teaching
  material it existed for — memory-access diagrams and cache simulations — is
  prose in `docs/performance.md` and does not need to ship as compiled objects.
- **`FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and `FFT_THREADED`.** All
  four were accepted by every planner and silently ignored. `FFT_THREADED` could
  not have been honoured anyway once `parallel_fft.c` went. Their bit positions
  are deliberately left unused rather than renumbered, so exported wisdom and
  already-compiled callers keep their meaning. `FFT_WISDOM_ONLY` survives because
  it is now implemented — see Added.

## [2.0.0] - 2025-08-30

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
- Comprehensive benchmark suite
- Unit tests
- Example programs
- Basic documentation

[Unreleased]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v3.0.0...HEAD
[3.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v2.0.0...v3.0.0
[2.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/muditbhargava66/FFT-implementation-in-C/releases/tag/v1.0.0
