# FFT Algorithms Overview

This guide provides an overview of all FFT algorithms implemented in this project, their characteristics, and when to use each one.

## Algorithm Comparison

| Algorithm | Complexity | Size Requirement | Best Use Case |
|-----------|------------|------------------|---------------|
| Radix-2 DIT | O(N log N) | Power of 2 | General purpose, industry standard |
| Radix-2 DIF | O(N log N) | Power of 2 | Hardware implementations |
| Radix-4 | O(N log N) | Power of 2 | Reduced operations vs Radix-2 |
| Split-Radix | O(N log N) | Power of 2 | Minimum operation count |
| Bluestein | O(N log N) | Any size | Non-power-of-2 sizes |
| Mixed-Radix | O(N log N) | Composite | Efficient for specific sizes |
| Recursive | O(N log N) | Power of 2 | Educational, clear structure |
| Naive DFT | O(N²) | Any size | Reference, small sizes |

## Core Algorithms

### Radix-2 Decimation-in-Time (DIT)

The most common FFT algorithm, implementing the Cooley-Tukey approach.

**Key Features:**
- In-place computation
- Bit-reversal permutation
- Butterfly operations
- Well-suited for vectorization

**Mathematical Foundation:**
The DFT is decomposed into even and odd indexed samples:

$$X[k] = \sum_{n=0}^{N/2-1} x[2n] W_N^{2nk} + W_N^k \sum_{n=0}^{N/2-1} x[2n+1] W_N^{2nk}$$

**Usage:**
```c
radix2_dit_fft(signal, n, FFT_FORWARD);
```

### Radix-2 Decimation-in-Frequency (DIF)

Alternative formulation that decimates in the frequency domain.

**Key Features:**
- Reversed butterfly structure compared to DIT
- No bit-reversal needed at output for some applications
- Natural for hardware pipelines

**Usage:**
```c
radix2_dif_fft(signal, n, FFT_FORWARD);
```

### Radix-4 FFT

Processes 4 samples at a time, reducing the number of stages.

**Key Features:**
- 25% fewer complex multiplications than Radix-2
- Works for any power of 2 (one leading radix-2 stage for odd powers)
- True radix-4 butterflies; twiddles advance multiplicatively per stage
- Costs an O(n) digit-reversal buffer, which is why `FFT_CONSERVE_MEMORY`
  selects in-place radix-2 DIT instead

**Mathematical Basis:**
Decomposes DFT into 4 interleaved subsequences:

$$X[k] = \sum_{r=0}^{3} W_N^{rk} \sum_{n=0}^{N/4-1} x[4n+r] W_{N/4}^{nk}$$

**Usage:**
```c
radix4_fft(signal, n, FFT_FORWARD);
```

### Split-Radix FFT

Combines Radix-2 and Radix-4 decompositions optimally.

**Key Features:**
- Lowest operation count among power-of-2 algorithms
- About 33% fewer multiplications than Radix-2
- More complex implementation

**Usage:**
```c
split_radix_fft(signal, n, FFT_FORWARD);
```

### Bluestein's Algorithm (Chirp Z-Transform)

Computes FFT for arbitrary sizes using convolution.

**Key Features:**
- Works for any size (not just powers of 2)
- Uses convolution theorem
- Higher memory requirements

**Mathematical Approach:**
Reformulates DFT as convolution:

$$X[k] = W_N^{k^2/2} \sum_{n=0}^{N-1} (x[n] W_N^{n^2/2}) W_N^{-(k-n)^2/2}$$

**Usage:**
```c
bluestein_fft(signal, n, FFT_FORWARD);  // n can be any size
```

### Mixed-Radix FFT

Handles composite sizes by factoring into smaller DFTs.

**Key Features:**
- Efficient for highly composite numbers
- Good for sizes like 12, 24, 36, etc.
- Recursive decomposition

**Usage:**
```c
mixed_radix_fft(signal, n, FFT_FORWARD);  // n = product of small primes
```

## Algorithm Selection Guide

The recommendations below follow the measurements in the
[Performance Guide](performance.md), not the operation counts. Where the two
disagree, the measurements are the ones that describe what you will observe.

### For Power-of-2 Sizes

1. **General Purpose**: Use Radix-2 DIT
   - Well-tested and understood
   - In-place, no scratch allocation
   - The reference column every other kernel is measured against; nothing beats
     it by more than ~17% and only radix-4 does so above 256

2. **Fewer Multiplications**: Split-Radix, if you need the count rather than the
   time
   - About a third fewer complex multiplications than radix-2
   - **1.6-2.9x slower in wall time**, because the recursion moves roughly 2.5n
     elements per level and memory traffic, not flops, is the constraint
   - It also cannot run in place, so it allocates 6n complex numbers

3. **Hardware/FPGA**: Consider Radix-2 DIF or Radix-4
   - Better pipeline structure
   - Reduced memory access patterns

### For Arbitrary Sizes

1. **Prime Sizes**: Use Bluestein
   - The only O(N log N) option here: mixed-radix has no factor to split on and
     falls back to a direct O(N²) DFT, measured 241x slower at N=1021
   - Costs roughly 9N complex numbers of scratch

2. **Composite Sizes**: Still Bluestein, unless memory is the binding constraint
   - Measured faster than mixed-radix at every composite length tried: 2.1x at
     N=360, 3.2x at N=1000, 6.0x at N=4095
   - Mixed-radix is the right call when the 2N padding Bluestein needs does not
     fit; pass `FFT_CONSERVE_MEMORY` and the planner selects it for lengths that
     factor into 2, 3, 5 and 7
   - The same flag keeps power-of-two lengths on in-place radix-2 DIT instead of
     radix-4, whose digit-reversal permutation allocates an O(n) scratch buffer

3. **Small Sizes**: Radix-2 DIT still wins
   - Naive DFT is not competitive at any size: measured 32x slower than radix-2
     DIT at N=8, rising to 498x at N=128, because it evaluates a complex
     exponential per term
   - `optimized_dft` precomputes or accumulates its twiddles and does pull ahead
     of naive DFT from N=64 (15x radix-2, where naive is 261x), but it is still
     15-25x behind radix-2 throughout
   - Use the DFTs as references for checking the FFTs, which is what the test
     suite does

## Performance Characteristics

### Memory Access Patterns

- **Radix-2**: Good locality in later stages
- **Radix-4**: Fewer passes through data
- **Split-Radix**: Complex access pattern
- **Bluestein**: Requires additional workspace

### Parallelization Potential

How well each algorithm *would* parallelise, if a threaded kernel existed. None
does today — see `docs/optimizations.md`.

- **High**: Radix-2 (independent butterflies)
- **Medium**: Radix-4, Split-Radix
- **Low**: Recursive implementations

### Numerical Stability

All algorithms are numerically stable, but:
- Radix-2 has well-understood error propagation
- Bluestein may accumulate more rounding errors, because it runs three transforms
  of about 4n instead of one of n
- A fixed-point implementation would need careful scaling, but there is no
  fixed-point path in this library

## Implementation Details

### Common Optimizations

1. **Twiddle Factor Optimization**
   ```c
   // Optimized twiddle factor computation with special cases
   static inline complex_t twiddle_factor(int k, int n, fft_direction dir) {
       // Fast paths for common angles
       if (k == 0) return 1.0;
       if (k * 4 == n) return (dir == FFT_FORWARD) ? -I : I;  // ±j
       if (k * 2 == n) return -1.0;  // -1
       
       double angle = dir * TWO_PI * k / n;
       return cexp(I * angle);
   }
   ```

2. **Bit Reversal Optimization**
   - A four-stage swap network for `log2n <= 8`, a shift-and-accumulate loop
     above that. No lookup table and no SIMD — it is plain integer bit
     manipulation

3. **Compiler Optimization Hints**
   ```c
   // Branch prediction hints
   #define LIKELY(x)   __builtin_expect(!!(x), 1)
   #define UNLIKELY(x) __builtin_expect(!!(x), 0)
   
   // Force inlining for critical functions
   #define FORCE_INLINE __attribute__((always_inline)) inline
   ```

4. **Memory Management**
   - Proper error checking and validation
   - Cache-optimized memory access patterns
   - Aligned allocation (64-byte from `fft_alloc_complex`, 32-byte from
     `allocate_complex_array` on Windows), which helps cache-line and page
     behaviour; there is no SIMD path that would need it

## Next Steps

- Explore individual algorithm pages for detailed explanations
- Run benchmarks to compare performance on your hardware
- Check the [API Reference](api-reference.md) for function details
- See [Applications](applications.md) for practical usage examples
