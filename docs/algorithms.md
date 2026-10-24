# FFT Algorithms

## Comparison

| Algorithm | Complexity | Length | Notes |
| --- | --- | --- | --- |
| Radix-2 DIT | O(N log N) | power of 2 | in-place, no scratch |
| Radix-2 DIF | O(N log N) | power of 2 | same cost as DIT here |
| Radix-4 | O(N log N) | power of 2 | faster from N=512 up |
| Split-Radix | O(N log N) | power of 2 | fewest multiplications, slowest of the three in wall time |
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
needs an O(n) buffer, which is why radix-4 loses to radix-2 below 512 despite
doing fewer multiplications per stage.

```c
radix4_fft(signal, n, FFT_FORWARD);
```

## Split-radix

Combines radix-2 and radix-4 decompositions: one `n/2` subtransform and two
`n/4` subtransforms per level. It has the lowest multiplication count of the
power-of-two kernels, but this implementation moves about 2.5n elements per
level and measures 1.6-2.9x slower than radix-2 DIT. It allocates a scratch
pool; it does not run in place.

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
leaves fall back to a direct DFT. The planner chooses it for lengths that
factor into 2, 3, 5 and 7 when `FFT_CONSERVE_MEMORY` is set.

```c
mixed_radix_fft(signal, n, FFT_FORWARD);
```

Measured against Bluestein it is slower at every non-power-of-two length tested
(2.1x at N=360, 241x at N=1021). It remains the right choice when memory is the
constraint; pass `FFT_CONSERVE_MEMORY` to the planner.

## Choosing a kernel

The planner uses these rules by default (`FFT_ESTIMATE`):

- Power-of-two lengths below 512: radix-2 DIT.
- Power-of-two lengths from 512 up: radix-4.
- Non-power-of-two lengths that factor into small primes: Bluestein, unless
  `FFT_CONSERVE_MEMORY` selects mixed-radix.
- Other non-power-of-two lengths: Bluestein.

`FFT_MEASURE` times the candidates instead of using the rules, then caches the
result. It is worth using when a program plans once and executes many times.

`FFT_CONSERVE_MEMORY` takes precedence over measurement and selects in-place or
low-scratch kernels.

## Memory and stability

- Radix-2 DIT and DIF are in place.
- Radix-4 allocates an O(n) digit-reversal buffer.
- Split-radix allocates a scratch pool.
- Bluestein allocates roughly 9n complex numbers.
- Mixed-radix allocates O(n) scratch for its recursion.

All kernels use double-precision complex arithmetic and are compared against an
independent direct DFT in `tests/test_all.c`. Bluestein accumulates somewhat
more round-off than a single radix kernel because it runs three transforms of
about 4n, which the test tolerance accounts for.

There is no fixed-point path. There is no SIMD or threaded kernel; the
`docs/optimizations.md` page lists what would be required to add one.
