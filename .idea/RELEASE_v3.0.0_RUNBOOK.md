# v3.0.0 — PR, merge and release runbook

Release working notes. `.idea/` is listed in `.gitignore`; the release documents
in this folder are force-added to git deliberately, while `.idea/.agents` and
unrelated IDE files stay local. The PR body and release notes are read from
these files with `--body-file` / `--notes-file`.

Current state:

- Branch: `release/v3.0.0`, pushed and tracking `origin/release/v3.0.0`
- Tip: the audit/docs commit that adds this section (see `git log -1`)
- Pre-audit tip: `56b9351` — `docs: record the Linux link and Clang gate fixes`
- Branch commits: `56b9351`, `5a3948a`, `20edd3d`, `87e437f`, plus the audit/docs
  commit
- CI: workflow `CI` run `39836257536` at `56b9351` — all four jobs green
  (Linux GCC, Linux Clang, Windows MinGW-w64, docs strict). Runs `39834532847`
  and `39835342504` failed first: Linux lacked `-lm`, then Clang rejected an
  unused demo helper. Both are fixed on the branch. The audit/docs commit is a
  documentation-only follow-up; check its run with
  `gh run list --branch release/v3.0.0 --limit 1`.
- Repository description and topics are already updated on GitHub
- Homepage remains https://fft-implementation-in-c.readthedocs.io/
- Working tree clean; the two docs files that an editor had reverted to
  pre-release text were restored from `HEAD`

Files used below:

- `.idea/PR_BODY_v3.0.0.md` — pull request description
- `.idea/RELEASE_NOTES_v3.0.0.md` — GitHub Release notes
- `.idea/github_issue_replies.md` — ready-to-post replies for #3, #4, #5 and PR #7

---

## Audit: PR #7 and open issues (2026-10-23)

| Item | Status in v3.0.0 | After merge |
| --- | --- | --- |
| Issue #3 (radix-4 is not radix-4) | Fixed: real radix-4 kernel for every power of two | Close |
| Issue #4 (Debian compile errors) | Fixed: feature macros before headers, no sysctl; GCC 15.3 strict C99 passes 1361/1361 | Close |
| Issue #5 (Meson support) | Not in this release; CMake only | Keep open |
| PR #7 (true radix-4) | Superseded; do not merge | Close without merge |

**PR #7 verdict.** The submitted butterfly is a real radix-4 DIT and its n=16
"mismatch" was caused by the base branch's buggy `bit_reverse()` at n <= 16,
not by the PR. It is still not mergeable here:

- it accepts powers of four only, so n = 8, 32, 128, 512, 2048 and every other
  `2 * 4^k` length loses the support v3.0.0 provides through a leading radix-2
  stage;
- it reverts `radix4_fft()` to `void` and reports on stderr, while the v3.0.0
  API returns `fft_error_t` and must not terminate or print from the library;
- it allocates a plan and an n-entry twiddle table inside every call, whereas
  v3.0.0 has the public `fft_plan_*` API and accumulates stage twiddles;
- it ships no independent-DFT test, and its header change removes a correct
  `twiddle_factor()` fast path used by radix-2.

No code from PR #7 is incorporated. Close it with the draft in
`github_issue_replies.md`; do not merge it after the release PR.

**Issue #3 evidence.** `algorithms/core/radix4.c` has base-4 digit reversal
(lines 91-127), a leading radix-2 stage for odd `log2(n)` (138-152), and
`m *= 4` radix-4 stages with the four-point butterfly (155-202). The suite runs
radix-4 at n = 4, 8, 16, 32, 64, 256, 512, 1024, 4096 and 16384 and compares it
with an independent direct DFT (`tests/test_all.c` size list at 819-822). The
planner uses radix-4 above 256 points (`algorithms/auto/fft_auto.c:257`).

**Issue #4 evidence.** `include/fft_common.h` defines `_GNU_SOURCE` before any
system header, so `strdup`/`posix_memalign` are visible under strict C99; there
is no `strdup` call and no `<sys/sysctl.h>` or `sysctl` usage left. A container
build with GCC 15.3 and `-Werror=implicit-function-declaration` passed
1361/1361. CI's Linux GCC and Clang lanes also pass.

**Issue #5.** Meson is not part of v3.0.0 and is called out as out of scope in
the release notes. Keep it open as a future contribution; the draft reply
corrects the stale "optional pthreads" wording.

**Doc cleanup before closing #3.** Three stale pages were corrected with this
audit: `algorithms/core/README.md` (power-of-four requirement, O(1) memory, and
a removed iterative-FFT entry), `docs/algorithms.md` (comparison table said
Power of 4), and `docs/optimizations.md` (old `void radix4_fft` snippet).

---

## 1. Pre-flight

```bash
gh auth status
git fetch origin
git status -sb
git log --oneline --decorate origin/main..release/v3.0.0
```

Expected: authenticated as `muditbhargava66`, branch clean and ahead of
`origin/main`, and the only commit listed is `87e437f`.

Confirm the local build still passes (optional before opening the PR):

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/benchmarks/benchmark_all
python -m mkdocs build --strict -f docs/mkdocs.yml -d build/docs-site-final
```

The recorded release numbers are in `README.md` and `CHANGELOG.md`: 1361/1361
checks, 52 passed / 0 failed / 2 skipped benchmark rows in 17.1 s, 9 examples,
realtime 5/5, docs strict 1.36 s.

## 2. Open the pull request

CLI:

```bash
gh pr create \
  --base main \
  --head release/v3.0.0 \
  --title "Release v3.0.0 (2026-10-23)" \
  --body-file .idea/PR_BODY_v3.0.0.md
```

Browser alternative: open
https://github.com/muditbhargava66/FFT-implementation-in-C/pull/new/release/v3.0.0
and paste the contents of `.idea/PR_BODY_v3.0.0.md`.

Then capture the PR number and open it:

```bash
gh pr view --web
```

## 3. Wait for CI and review

```bash
gh pr checks --watch
gh run list --branch release/v3.0.0 --limit 5
```

Required jobs:

| Job | What it proves |
| --- | --- |
| `linux-gcc` | configure, warning-free build, `ctest`, benchmark `Failed: 0` |
| `linux-clang` | same under Clang |
| `windows-mingw` | MSYS2 UCRT64 build and tests |
| `docs` | `mkdocs build --strict` with the pinned docs requirements |

If a job fails:

```bash
gh run view <run-id> --log-failed
```

Fix on the branch, commit, and push; the open PR updates automatically. Do not
merge with a red required check.

Review the rendered PR description once. Confirm it still says "Fixes #3",
"Fixes #4", and explicitly leaves #5 (Meson) open.

## 4. Merge

Use a merge commit, matching the repository's existing history (the v2.0.0
release was merged the same way). Squashing loses the individual fix commits
that the changelog references.

```bash
gh pr merge --merge --delete-branch
```

If branch protection requires an approval or a different strategy, use
**Merge pull request** in the web UI instead; do not force-merge.

## 5. Post-merge: sync main and tag

```bash
git switch main
git pull --ff-only origin main
git log --oneline -5
git tag -a v3.0.0 -m "FFT Implementation in C v3.0.0 (2026-10-23)"
git push origin v3.0.0
```

Verify the tag points at the merge (or at the release commit reachable from
`main`):

```bash
git show --no-patch --oneline v3.0.0
git ls-remote --tags origin v3.0.0
```

If the branch was not deleted by `gh pr merge`:

```bash
git push origin --delete release/v3.0.0
git branch -d release/v3.0.0
git fetch --prune
```

## 6. Create the GitHub Release

```bash
gh release create v3.0.0 \
  --title "v3.0.0 — planner correctness, real split-radix and mixed-radix" \
  --notes-file .idea/RELEASE_NOTES_v3.0.0.md \
  --latest \
  --verify-tag
```

Verify:

```bash
gh release view v3.0.0
gh release list --limit 3
```

`--latest` moves the "Latest" badge from v2.0.0 to v3.0.0.

## 7. Read the Docs

Read the Docs builds `latest` from `main`, so the site updates only after the
merge, not from the release branch. A tag build is a separate version and is not
automatic unless the project is configured to build tags.

1. Open https://readthedocs.org/projects/fft-implementation-in-c/builds/ and
   confirm a build for `latest` started from the merge commit. If none started,
   use **Build version** for `latest`.
2. Open https://readthedocs.org/projects/fft-implementation-in-c/versions/ and
   activate `v3.0.0` if it appears. Set it as the default/stable version if that
   matches the project's policy.
3. Confirm the live pages:
   - https://fft-implementation-in-c.readthedocs.io/ shows Release 3.0.0 and
     1361/1361 checks
   - navigation has Home, Getting Started, Algorithms, Applications,
     API Reference, Performance/Benchmarks, Optimization Guide
   - no 404s, no raw Markdown in navigation, no missing assets
   - MathJax renders the DFT equation on the home page
4. If the build fails, read the build log. `mkdocs build --strict` passes
   locally, so a failure is usually an environment or requirements issue in
   `.readthedocs.yaml` (Ubuntu 24.04, Python 3.12).

The version selector on the live site will keep showing 2.0.0 until RTD finishes
the `latest` rebuild; that is expected while the PR is open.

## 8. Close the fixed issues

Only after the release is live. Paste the full drafts from
`.idea/github_issue_replies.md`; the short commands below are the minimum:

```bash
gh issue close 3 --comment "Fixed in v3.0.0 (merge of release/v3.0.0). radix-4 now uses real base-4 stages for every power-of-two length, and the stale power-of-four/O(1) docs are corrected."
gh issue close 4 --comment "Fixed in v3.0.0 (merge of release/v3.0.0). Feature-test macros are set before system headers, strdup is gone, and no sysctl header remains. GCC 15.3 strict C99 passes 1361/1361."
gh pr close 7 --comment "Superseded by the v3.0.0 radix-4 implementation. The new kernel supports every power of two and keeps the fft_error_t API; the PR accepts powers of four only and reverts the API to void. The n=16 mismatch came from the base branch's bit_reverse bug, which v3.0.0 fixes."
```

Leave #5 open (Meson build support is not part of this release); optionally add
a comment saying it is still welcome and correct the old "optional pthreads"
wording using the draft.

## 9. Post-release verification checklist

- [ ] `git log --oneline -1 origin/main` shows the v3.0.0 merge
- [ ] `git ls-remote --tags origin v3.0.0` returns the tag
- [ ] GitHub Release `v3.0.0` is published and marked Latest
- [ ] CI on `main` is green after the merge
- [ ] Read the Docs `latest` build is green and shows 3.0.0
- [ ] Repository description and 20 topics are unchanged from this run
- [ ] Issues #3 and #4 closed, #5 still open
- [ ] PR #7 closed without merge, with the superseded explanation
- [ ] `release/v3.0.0` deleted locally and remotely
- [ ] `README.md` release badge and `CHANGELOG.md` version match the tag
- [ ] `algorithms/core/README.md`, `docs/algorithms.md` and
      `docs/optimizations.md` no longer claim power-of-four-only radix-4

## 10. Rollback notes

If a problem is found before the release is announced:

- Revert the merge on `main` with a new PR:
  `git revert -m 1 <merge-sha>` on a fresh branch, then open a PR. Do not
  rewrite `main`.
- If the tag was pushed but the release is wrong, delete the GitHub Release
  first, then delete and re-push the tag only with explicit intent:
  `gh release delete v3.0.0`, `git push origin :refs/tags/v3.0.0`,
  `git tag -d v3.0.0`. This is destructive for anyone who already fetched the
  tag; prefer a follow-up `v3.0.1` once users may have seen it.
- Read the Docs can rebuild an older commit from the Versions page if the
  `latest` build is bad; fix forward on `main`.
