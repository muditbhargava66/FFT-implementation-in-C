# FFT Algorithms

## Comparison

| Algorithm | Complexity | Length | Notes |
| --- | --- | --- | --- |
| Radix-2 DIT | O(N log N) | power of 2 | in-place, no scratch |
| Radix-2 DIF | O(N log N) | power of 2 | same cost as DIT here |
| Radix-4 | O(N log N) | power of 2 | faster from N=256 up |
| Radix-4 AVX2 | O(N log N) | power of 2 | two complex values per AVX2 register; selected from N=4096 |
| Radix-4 OpenMP | O(N log N) | power of 2 | opt-in; parallel butterflies; selected from N=524288 |
| Four-Step | O(N log N) | power of 2 | cache-blocked; `FFT_MEASURE` candidate |
| Split-Radix | O(N log N) | power of 2 | fewest multiplications, still slower in wall time |
| Bluestein | O(N log N) | any | about 9N complex scratch |
| Mixed-Radix | O(N log N) | any | recursive Cooley-Tukey for composite lengths; direct DFT at a prime |
| Recursive | O(N log N) | power of 2 | clear structure, not the fastest |
| Naive DFT | O(N²) | any | reference |
| Optimized DFT | O(N²) | any | reference with twiddle reuse |

## Radix-2 DIT

The standard Cooley-Tukey decomposition into even and odd samples:

$$X[k] = \sum_{n=0}^{N/2-1} x[2n]\,W_N^{2nk}
      + W_N^k \sum_{n=0}^{N/2-1} x[2n+1]\,W_N^{2nk}.$$

The kernel permutes the input into bit-reversed order and then runs
log₂N butterfly stages in place. It allocates no scratch and is the reference
column in the performance tables.

```c
radix2_dit_fft(signal, n, FFT_FORWARD);
```

## Radix-2 DIF

Decimation in frequency. The butterflies run in the opposite order and the
bit-reversal permutation is applied at the end. Its measured cost is the same as
DIT at every size tested, so choose between them for the data layout, not for
speed.

```c
radix2_dif_fft(signal, n, FFT_FORWARD);
```

## Radix-4

Four-way decomposition. It accepts any power of two and uses one leading radix-2
stage when `log2(n)` is odd.

$$X[k] = \sum_{r=0}^{3} W_N^{rk} \sum_{n=0}^{N/4-1} x[4n+r]\,W_{N/4}^{nk}.$$

Twiddles advance multiplicatively per stage. The digit-reversal permutation
needs an O(n) buffer, which is why radix-4 loses to radix-2 below 256 points
despite doing fewer multiplications per stage.

```c
radix4_fft(signal, n, FFT_FORWARD);
```

## Radix-4 AVX2

The same decomposition with the stage arithmetic vectorized. Two complex values
fit in one 256-bit register, and the twiddle and butterfly arithmetic use
`__m256d` throughout. Dispatch is runtime-only: `radix4_avx2_available()` calls
`__builtin_cpu_supports("avx2")`, and the planner selects the kernel from 4096
points up when it returns true. The scalar kernel remains the fallback on
non-x86 targets and on builds without GCC/Clang vector support; those builds get
a stub that reports unavailable.

```c
radix4_avx2_fft(signal, n, FFT_FORWARD);
```

## Radix-4 OpenMP

An opt-in threaded version of the radix-4 stage, built with
`-DFFT_ENABLE_OPENMP=ON`. The butterfly loop is parallelized over blocks and
indices, and one twiddle table per stage is shared by the threads. Below 524288
points, or with a single thread, the kernel calls scalar radix-4 instead. When
the build has no OpenMP, the entry point returns `FFT_ERR_UNSUPPORTED` and
`radix4_omp_available()` returns 0. The planner selects it from 524288 points
when it is available.

```c
radix4_omp_fft(signal, n, FFT_FORWARD);
```

## Four-Step

For large power-of-two lengths, `N = N1*N2` with the two factors near
`sqrt(N)`. The input is treated as an N1-by-N2 row-major matrix: columns are
gathered into cache-sized scratch and transformed, the result is multiplied by
`W_N^(k1*c)`, rows are transformed, and a 16x16 tiled transpose returns natural
order. Sub-transforms use the radix-4 kernel and share one scratch buffer.

```c
four_step_fft(signal, n, FFT_FORWARD);
```

## Split-radix

Combines radix-2 and radix-4 decompositions: one `n/2` subtransform and two
`n/4` subtransforms per level. It has the lowest multiplication count of the
power-of-two kernels. The recursion reads its input with a stride instead of
gathering even and odd samples into scratch at every level, and the two combine
passes share one loop, so peak scratch is about 3n (the pool is sized 4n).
Wall-time ratios are in the [Performance Guide](performance.md).

```c
split_radix_fft(signal, n, FFT_FORWARD);
```

## Bluestein

Reformulates the DFT as a convolution so any length can use power-of-two
transforms:

$$X[k] = W_N^{k^2/2} \sum_{n=0}^{N-1}
        \left(x[n]\,W_N^{n^2/2}\right) W_N^{-(k-n)^2/2}.$$

The convolution is padded to `next_power_of_two(2n-1)`, which is roughly 4n,
and three transforms of that length are run. That is why it costs about 10x a
radix kernel at a power-of-two length, and why it still wins for lengths where
mixed-radix has no factor to use.

```c
bluestein_fft(signal, n, FFT_FORWARD);   /* any n up to 2^29 */
```

## Mixed-radix

Recursive Cooley-Tukey using the smallest prime factor at each step. Prime
leaves fall back to a direct DFT. Under `FFT_CONSERVE_MEMORY` the planner
chooses it for lengths whose prime factors are all 2, 3, 5 or 7 and number at
least three; smaller composites still go to Bluestein.

```c
mixed_radix_fft(signal, n, FFT_FORWARD);
```

Measured against Bluestein it is slower at every non-power-of-two length tested
(2.1x at N=360, 241x at N=1021). It remains the right choice when memory is the
constraint; pass `FFT_CONSERVE_MEMORY` to the planner.

## Choosing a kernel

The planner uses these rules by default (`FFT_ESTIMATE`):

- Power-of-two lengths up to 256: radix-2 DIT.
- Power-of-two lengths from 257 to 4095: radix-4.
- From 4096 to 524287: AVX2 radix-4 when the CPU supports AVX2, else radix-4.
- From 524288 up: radix-4 OpenMP when built with OpenMP and more than one
  thread is available, else AVX2 radix-4 when available, else radix-4.
- `FFT_CONSERVE_MEMORY`: radix-2 DIT for every power-of-two length.
- Non-power-of-two lengths: Bluestein, unless `FFT_CONSERVE_MEMORY` is set and
  the length is highly composite (at least three factors drawn from 2, 3, 5 and
  7), in which case mixed-radix.

The four-step kernel is not in the heuristic. It is a candidate for
`FFT_MEASURE`, which picks it only when it actually wins; on the verification
machines it was 17-30% slower than AVX2/OpenMP at 1M-4M points (see the
[Performance Guide](performance.md)).

The fixed thresholds sit inside the measured crossovers on the verification
machine; a different CPU can move them. `FFT_MEASURE` times the candidates and
keeps the heuristic unless a challenger wins the sweep by more than 15% and
repeats the win in a paired re-run. The choice is then cached.

`FFT_CONSERVE_MEMORY` takes precedence over measurement and selects in-place or
low-scratch kernels.

## Memory and stability

- Radix-2 DIT and DIF are in place.
- Radix-4 allocates an O(n) digit-reversal buffer; the AVX2 and OpenMP versions
  use the same reorder.
- Four-step allocates n elements for the reordered output plus about
  `2*sqrt(n)` for the column and sub-transform scratch.
- Split-radix allocates a 4n-element pool and peaks near 3n live.
- Bluestein allocates roughly 9n complex numbers.
- Mixed-radix allocates O(n) scratch for its recursion.

All kernels use double-precision complex arithmetic and are compared against an
independent direct DFT in `tests/test_all.c`. Bluestein accumulates somewhat
more round-off than a single radix kernel because it runs three transforms of
about 4n, which the test tolerance accounts for.

There is no fixed-point path. SIMD and threading exist as the AVX2 and OpenMP
kernels above; both have their own availability checks and test entries.
