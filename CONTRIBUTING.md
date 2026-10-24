# Contributing

Bug reports, fixes, tests and documentation improvements are welcome. This page
covers the build, the checks a change has to pass, and the house style.

## Getting started

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

If you are working from a fork, add the original as a second remote:

```bash
git remote add upstream https://github.com/muditbhargava66/FFT-implementation-in-C.git
```

Create a branch for your change:

```bash
git checkout -b fix/bit-reversal-size
```

## Checks before opening a pull request

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/bin/benchmark_all
```

The requirements:

- The build must be warning-free under `-Wall -Wextra`. CI fails on the first
  `warning:` line.
- `test_all` must pass every check. It compares each implementation against a
  direct DFT computed independently inside the test file.
- `benchmark_all` must report `Failed: 0`. The binary exits 0 even when rows
  fail, so read the summary.
- Documentation changes must pass the strict MkDocs build:

  ```bash
  python -m pip install -r docs/requirements.txt
  python -m mkdocs build --strict -f docs/mkdocs.yml -d /tmp/docs-site
  ```

  The docs build needs Python 3.9 or newer.

## Adding an algorithm

1. Put core FFT kernels in `algorithms/core/`, DFT references in
   `algorithms/dft/`, and planner code in `algorithms/auto/`.
2. Include `../../include/fft_common.h` and `../../include/fft_algorithms.h`.
3. Return `fft_error_t` and leave the caller's buffer untouched on failure.
4. Declare the function in `include/fft_algorithms.h`. A kernel that no header
   declares cannot be used and should not be committed.
5. If the planner should be able to choose it, add it to `select_algorithm()`
   and to the timing sweep in `algorithms/auto/fft_auto.c`.
6. Add a reference comparison and edge-case checks to `tests/test_all.c`.
7. Add it to the benchmark registry in `benchmarks/benchmark_all.c`.

Every `.c` file under `algorithms/`, `utils/` and `applications/` is picked up
by the CMake glob. That means an unreachable kernel still compiles; the tests
and the planner integration are what make it part of the library.

### Demonstration code

Most algorithm files end with a `main()` guarded by `#ifndef LIB_BUILD`. The
library target defines `LIB_BUILD`, so the demo is excluded from
`libfft_lib_c.a`. Keep new demo code inside that guard, and do not define a
second unguarded `main()` anywhere under `utils/` or `algorithms/`.

## C style

- C11, 4-space indentation, no tabs.
- Functions and variables in `snake_case`; types in `snake_case_t`; macros in
  `UPPER_CASE`.
- Declare at the top of a block only if it helps readability; C11 permits
  declarations where they are used.
- Check allocation results. Return an error code rather than calling `exit()` or
  `abort()`.
- Keep comments about constraints and non-obvious arithmetic. Avoid comments
  that restate the next line or narrate what an earlier version did; that
  belongs in the commit message.

## Tests

`tests/test_all.c` is the single test binary. Add checks with `report()`:

```c
report("scope", "what is checked", condition, error_value, "detail");
```

A new implementation should be checked against `reference_bin()`, which
computes a DFT bin by direct summation and shares no code with the library.
Property tests are useful, but they cannot catch a consistently wrong bin
ordering; the reference comparison can.

## Documentation

- Public functions belong in `docs/api-reference.md`, with the signature copied
  from the header.
- Performance claims need a reproducible method and the machine they were
  measured on.
- Code samples are compiled before merging. If an example uses a helper, make
  sure that helper exists in `include/` or is defined in the same snippet.

## Commits and pull requests

Use a short imperative subject, optionally with a conventional prefix
(`fix:`, `docs:`, `test:`):

```
fix: reject execute_dft on real-valued plans
```

Explain the reason in the body if it is not obvious from the diff, and reference
an issue when one exists. In the pull request, list the commands you ran and
their results. The template at `.github/pull_request_template.md` covers the
build, the test suite, the warning gate, the benchmark smoke test and the strict
docs build.

## Security

Do not open a public issue for a suspected vulnerability. Follow
[SECURITY.md](SECURITY.md) and use GitHub's private vulnerability reporting from
the repository's Security tab. Include the affected version, platform, compiler,
build flags and a minimal reproducer.
