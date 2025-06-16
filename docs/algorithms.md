# FFT Algorithms Overview

This guide provides an overview of all FFT algorithms implemented in this project, their characteristics, and when to use each one.

## Algorithm Comparison

| Algorithm | Complexity | Size Requirement | Best Use Case |
|-----------|------------|------------------|---------------|
| Radix-2 DIT | O(N log N) | Power of 2 | General purpose, industry standard |
| Radix-2 DIF | O(N log N) | Power of 2 | Hardware implementations |
| Radix-4 | O(N log N) | Power of 4 | Reduced operations vs Radix-2 |
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
- Requires size to be power of 4
- More complex butterfly structure

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

### For Power-of-2 Sizes

1. **General Purpose**: Use Radix-2 DIT
   - Well-tested and understood
   - Good cache behavior
   - Easy to optimize

2. **Minimum Operations**: Use Split-Radix
   - Best theoretical complexity
   - Worth it for large transforms

3. **Hardware/FPGA**: Consider Radix-2 DIF or Radix-4
   - Better pipeline structure
   - Reduced memory access patterns

### For Arbitrary Sizes

1. **Prime Sizes**: Use Bluestein
   - Only option for prime sizes
   - Reasonable performance

2. **Composite Sizes**: Use Mixed-Radix
   - Better than Bluestein for highly composite numbers
   - Exploits factorization

3. **Small Sizes**: Consider Naive DFT
   - Simple and direct
   - Competitive for N < 32

## Performance Characteristics

### Memory Access Patterns

- **Radix-2**: Good locality in later stages
- **Radix-4**: Fewer passes through data
- **Split-Radix**: Complex access pattern
- **Bluestein**: Requires additional workspace

### Parallelization Potential

- **High**: Radix-2 (independent butterflies)
- **Medium**: Radix-4, Split-Radix
- **Low**: Recursive implementations

### Numerical Stability

All algorithms are numerically stable, but:
- Radix-2 has well-understood error propagation
- Bluestein may accumulate more rounding errors
- Fixed-point implementations need careful scaling

## Implementation Details

### Common Optimizations

1. **Twiddle Factor Precomputation**
   ```c
   // Precompute once
   complex_t* twiddles = precompute_twiddles(n);
   
   // Reuse for multiple FFTs
   radix2_dit_fft_with_twiddles(signal, n, twiddles, FFT_FORWARD);
   ```

2. **Real-valued FFT**
   - Exploit conjugate symmetry
   - Half the computation and storage

3. **Cache Optimization**
   - Block algorithms for large sizes
   - Minimize cache misses

## Next Steps

- Explore individual algorithm pages for detailed explanations
- Run benchmarks to compare performance on your hardware
- Check the [API Reference](api-reference.md) for function details
- See [Applications](applications.md) for practical usage examples
