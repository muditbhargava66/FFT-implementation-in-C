## Summary

<!-- What changed, and why. Link the issue this closes: "Fixes #123". -->

## Type of change

- [ ] Bug fix (correctness, memory safety, build)
- [ ] New feature or algorithm
- [ ] Performance change
- [ ] Tests
- [ ] Documentation
- [ ] CI / tooling
- [ ] Release metadata

## Verification

Paste the commands you ran and their results. Do not tick a box you did not
run; a reviewer will ask for the output.

- [ ] Configured and built cleanly:

  ```bash
  cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
  cmake --build build
  ```

- [ ] `ctest --test-dir build --output-on-failure` — expected **1406/1406**
      checks in `test_all` on an AVX2-capable x86 build and **54/54** in
      `test_utils`, 2/2 CTest tests. Hosts without AVX2 skip that kernel's
      entries and report fewer checks. An OpenMP build expects **1415/1415**
      plus **17712** threaded checks and 3/3 CTest tests; a GPU backend adds one
      more test. GPU tests report a skip without a device.
- [ ] No new compiler warnings: Linux and macOS CI builds with `-Werror`, and
      every C job also fails on a `warning:` line in the build log.
- [ ] `./build/bin/benchmark_all` reports `Failed: 0` (include the row count and
      elapsed time for performance changes).
- [ ] Example programs affected by this change were run; paste the relevant
      output.
- [ ] Docs build, when documentation changed:

  ```bash
  python -m mkdocs build --strict -f docs/mkdocs.yml -d "$RUNNER_TEMP/docs-site"
  ```

## If a kernel or planner changed

- [ ] Compared against the independent direct DFT in `tests/test_all.c`, not
      against another library kernel.
- [ ] Wisdom/cache paths were checked with an invalid entry as well as a valid
      one.
- [ ] Performance figures are ratios from a seeded, batched run. Single-shot
      `clock()` readings are quantisation on Windows, not measurements.
- [ ] Numerical claims state the harness, compiler, flags, and whether
      `-march=native` was enabled (it is off by default).

## Checklist

- [ ] The change is focused; unrelated refactors and formatting are separate.
- [ ] No build directories, binaries, logs, generated docs, or local paths are
      committed.
- [ ] No secrets, credentials, API keys, or tokens are included.
- [ ] Public headers and `docs/api-reference.md` are updated if the API changed.
- [ ] `CHANGELOG.md` has an entry under `[Unreleased]` (or the release section
      for a release PR).
- [ ] New dependencies are justified and reflected in the relevant manifest.
- [ ] Security-sensitive changes follow [SECURITY.md](SECURITY.md).
