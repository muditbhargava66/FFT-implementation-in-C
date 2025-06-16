# Contributing to FFT Study Repository

Thank you for your interest in contributing to this educational FFT repository! This document provides guidelines for contributing code, documentation, and improvements.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How to Contribute](#how-to-contribute)
- [Development Process](#development-process)
- [Code Standards](#code-standards)
- [Testing Requirements](#testing-requirements)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

This project adheres to a code of conduct that all contributors are expected to follow:

- Be respectful and inclusive
- Welcome newcomers and help them get started
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/yourusername/FFT-implementation-in-C.git
   cd FFT-implementation-in-C
   ```
3. Add the upstream repository:
   ```bash
   git remote add upstream https://github.com/originalowner/FFT-implementation-in-C.git
   ```
4. Create a new branch for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## How to Contribute

### Types of Contributions

We welcome the following types of contributions:

1. **New FFT Algorithms**
   - Implement additional FFT variants (e.g., Prime Factor Algorithm, Rader's Algorithm)
   - Add optimized versions for specific use cases

2. **Performance Optimizations**
   - SIMD implementations for different architectures
   - Cache-aware algorithms
   - GPU implementations (CUDA/OpenCL)

3. **Applications**
   - Signal processing examples
   - Real-world use cases
   - Educational demonstrations

4. **Documentation**
   - Algorithm explanations
   - Tutorial improvements
   - Code comments and examples

5. **Bug Fixes**
   - Fix implementation errors
   - Improve numerical stability
   - Address edge cases

6. **Testing**
   - Add test cases
   - Improve test coverage
   - Performance regression tests

### Adding a New Algorithm

When adding a new FFT algorithm:

1. Create a new file in the appropriate directory:
   - `algorithms/core/` for FFT variants
   - `algorithms/dft/` for DFT implementations
   - `applications/` for applications
   - `optimizations/` for optimized versions

2. Include the common header:
   ```c
   #include "../../include/fft_common.h"
   ```

3. Implement the algorithm with clear comments explaining the approach

4. Add a comprehensive main() function demonstrating:
   - Basic usage
   - Performance characteristics
   - Accuracy verification
   - Comparison with other methods

5. Update the benchmark suite to include your implementation

6. Add appropriate test cases

## Development Process

### Code Structure

```
algorithm_name.c
├── File header comment (algorithm description)
├── Include statements
├── Type definitions (if needed)
├── Helper functions
├── Main algorithm implementation
├── Wrapper functions
├── Demonstration/test functions
└── main() with examples
```

### Example Template

```c
/**
 * Algorithm Name Implementation
 * 
 * Brief description of the algorithm and its characteristics.
 * 
 * Time Complexity: O(?)
 * Space Complexity: O(?)
 * 
 * References:
 * [1] Author, "Paper Title", Journal, Year
 */

#include "../../include/fft_common.h"

// Implementation-specific definitions
typedef struct {
    // ...
} algorithm_data_t;

// Main algorithm function
void algorithm_name_fft(complex_t* x, int n, fft_direction dir) {
    // Input validation
    if (!validate_input(x, n)) {
        return;
    }
    
    // Algorithm implementation
    // ...
    
    // Handle inverse transform
    if (dir == FFT_INVERSE) {
        // Scale output
    }
}

// Convenience wrappers
void fft_algorithm_name(complex_t* x, int n) {
    algorithm_name_fft(x, n, FFT_FORWARD);
}

void ifft_algorithm_name(complex_t* x, int n) {
    algorithm_name_fft(x, n, FFT_INVERSE);
}

// Demonstration and testing
int main() {
    printf("Algorithm Name Implementation\n");
    printf("=============================\n");
    
    // 1. Basic functionality demo
    // 2. Performance measurement
    // 3. Accuracy verification
    // 4. Comparison with reference
    // 5. Special cases
    
    return 0;
}
```

## Code Standards

### C Style Guide

1. **Indentation**: 4 spaces (no tabs)

2. **Naming Conventions**:
   - Functions: `snake_case`
   - Types: `snake_case_t`
   - Constants: `UPPER_CASE`
   - Macros: `UPPER_CASE`

3. **Comments**:
   - Use block comments for function descriptions
   - Use line comments for complex logic
   - Document all parameters and return values

4. **Functions**:
   - Keep functions focused and concise
   - Validate inputs
   - Handle errors gracefully
   - Free allocated memory

5. **Memory Management**:
   - Always check allocation success
   - Free memory in reverse order of allocation
   - Use provided allocation helpers

### Example Code Style

```c
/**
 * Compute the FFT of a complex array.
 * 
 * @param x     Input/output array
 * @param n     Size of array (must be power of 2)
 * @param dir   Transform direction (FFT_FORWARD or FFT_INVERSE)
 * @return      0 on success, -1 on error
 */
int compute_fft(complex_t* x, int n, fft_direction dir) {
    // Validate inputs
    if (!x || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return -1;
    }
    
    if (!is_power_of_two(n)) {
        fprintf(stderr, "Error: Size must be power of 2\n");
        return -1;
    }
    
    // Algorithm implementation
    for (int stage = 1; stage <= log2_int(n); stage++) {
        // Process butterflies
        process_stage(x, n, stage, dir);
    }
    
    return 0;
}
```

## Testing Requirements

All contributions must include appropriate tests:

1. **Unit Tests**:
   - Test basic functionality
   - Test edge cases (n=1, n=2, large n)
   - Test error handling

2. **Accuracy Tests**:
   - Compare against reference implementation
   - Verify mathematical properties
   - Test special inputs (impulse, DC, etc.)

3. **Performance Tests**:
   - Measure execution time
   - Compare with existing implementations
   - Test scaling behavior

4. **Integration**:
   - Update `tests/test_all.c` to include new tests
   - Ensure all existing tests still pass

## Documentation

### Code Documentation

Every source file must include:

1. **File Header**:
   ```c
   /**
    * Brief description
    * 
    * Detailed explanation of the algorithm/application
    * 
    * Features:
    * - Feature 1
    * - Feature 2
    * 
    * Limitations:
    * - Limitation 1
    * 
    * References:
    * [1] Citation
    */
   ```

2. **Function Documentation**:
   ```c
   /**
    * Brief function description.
    * 
    * Detailed explanation if needed.
    * 
    * @param param1  Description
    * @param param2  Description
    * @return        Description
    */
   ```

3. **Inline Comments**:
   - Explain complex algorithms
   - Clarify non-obvious code
   - Reference equations or papers

### README Files

When adding a new category or significant feature, include a README.md:

```markdown
# Feature Name

Brief description of the feature.

## Overview

Detailed explanation of what this provides.

## Usage

```c
// Example code
```

## Theory

Mathematical or algorithmic background.

## Performance

Performance characteristics and benchmarks.

## References

1. [Paper/Book citations]
```

## Submitting Changes

### Pull Request Process

1. **Ensure your code**:
   - Follows the style guide
   - Includes tests
   - Has proper documentation
   - Builds without warnings
   - Passes all tests

2. **Update documentation**:
   - Add your algorithm to the main README
   - Update relevant section READMEs
   - Include usage examples

3. **Commit messages**:
   ```
   Short summary (50 chars or less)
   
   More detailed explanation if needed. Explain what changed
   and why. Reference any issues being fixed.
   
   Fixes #123
   ```

4. **Create pull request**:
   - Use a descriptive title
   - Fill out the PR template
   - Reference related issues
   - Include benchmark results

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] New algorithm implementation
- [ ] Bug fix
- [ ] Performance optimization
- [ ] Documentation update
- [ ] New application/example

## Testing
- [ ] All tests pass
- [ ] Added new tests
- [ ] Tested on multiple platforms

## Performance Impact
Include benchmark results comparing before/after

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Comments added for complex code
- [ ] Documentation updated
- [ ] No new warnings
```

## Questions?

If you have questions about contributing:

1. Check existing issues and discussions
2. Read the documentation thoroughly
3. Open a new issue with the question tag
4. Join our community discussions

Thank you for contributing to make this the best FFT educational resource!
