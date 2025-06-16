# Core FFT Algorithms

This directory contains implementations of fundamental FFT algorithms.

## Algorithms Included

### 1. Radix-2 DIT (radix2_dit.c)
- **Description**: Decimation-in-Time radix-2 FFT
- **Complexity**: O(n log n)
- **Requirements**: Input size must be power of 2
- **Best for**: General purpose, well-balanced performance

### 2. Radix-2 DIF (radix2_dif.c)
- **Description**: Decimation-in-Frequency radix-2 FFT
- **Complexity**: O(n log n)
- **Requirements**: Input size must be power of 2
- **Best for**: Systems where output can be in bit-reversed order

### 3. Radix-4 (radix4.c)
- **Description**: Processes 4 points at a time
- **Complexity**: O(n log n) with 25% fewer multiplications than radix-2
- **Requirements**: Input size must be power of 4
- **Best for**: When multiplication is expensive

### 4. Split-Radix (split_radix.c)
- **Description**: Hybrid of radix-2 and radix-4
- **Complexity**: O(n log n) with lowest operation count
- **Requirements**: Input size must be power of 2
- **Best for**: Maximum computational efficiency

### 5. Bluestein's Algorithm (bluestein.c)
- **Description**: FFT for arbitrary input sizes using chirp z-transform
- **Complexity**: O(n log n)
- **Requirements**: None (works for any size)
- **Best for**: Non-power-of-2 sizes

### 6. Mixed-Radix (mixed_radix.c)
- **Description**: Handles composite sizes by factorization
- **Complexity**: O(n log n) for highly composite n
- **Requirements**: None
- **Best for**: Sizes with small prime factors

### 7. Recursive FFT (recursive_fft.c)
- **Description**: Educational implementation showing divide-and-conquer
- **Complexity**: O(n log n)
- **Requirements**: Input size must be power of 2
- **Best for**: Understanding FFT algorithm structure

### 8. Iterative FFT (iterative_fft.c)
- **Description**: Memory-efficient iterative implementation
- **Complexity**: O(n log n)
- **Requirements**: Input size must be power of 2
- **Best for**: Embedded systems, cache efficiency

## Quick Comparison

| Algorithm | Multiplications | Memory | Flexibility | Speed |
|-----------|----------------|---------|-------------|--------|
| Radix-2 DIT | n/2 log₂(n) | O(1) | Power of 2 | ★★★★☆ |
| Radix-2 DIF | n/2 log₂(n) | O(1) | Power of 2 | ★★★★☆ |
| Radix-4 | 3n/8 log₂(n) | O(1) | Power of 4 | ★★★★☆ |
| Split-Radix | ~n/3 log₂(n) | O(n) | Power of 2 | ★★★★★ |
| Bluestein | O(n log n) | O(n) | Any size | ★★★☆☆ |
| Mixed-Radix | Varies | O(n) | Any size | ★★★☆☆ |
| Recursive | n/2 log₂(n) | O(n log n) | Power of 2 | ★★☆☆☆ |
| Iterative | n/2 log₂(n) | O(1) | Power of 2 | ★★★★☆ |

## Usage Example

```c
#include "fft_common.h"

int main() {
    int n = 1024;
    complex_t* signal = allocate_complex_array(n);
    
    // Generate signal
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 10 * i / n);
    }
    
    // Apply FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Process spectrum...
    
    // Inverse FFT
    radix2_dit_fft(signal, n, FFT_INVERSE);
    
    free_complex_array(signal);
    return 0;
}
```

## Building

Each algorithm can be built independently:

```bash
make radix2_dit
make radix4
# etc...
```

Or build all core algorithms:

```bash
make core
```

## Testing

Run the built-in tests in each file:

```bash
./radix2_dit
./split_radix
```

Or run the comprehensive test suite:

```bash
make test
./bin/test_all
```
