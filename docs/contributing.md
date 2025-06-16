# Contributing to FFT Implementation in C

Thank you for your interest in contributing to this project! This document provides guidelines and instructions for contributing.

## Code of Conduct

By participating in this project, you agree to maintain a respectful and inclusive environment for all contributors.

## How to Contribute

### 1. Reporting Issues

- Use the GitHub issue tracker to report bugs
- Check existing issues before creating a new one
- Include minimal reproducible examples
- Provide system information (OS, compiler, etc.)

#### Bug Report Template

```markdown
**Description**
Clear description of the bug

**To Reproduce**
1. Build with `make ...`
2. Run `./bin/...`
3. See error

**Expected behavior**
What should happen

**System Info**
- OS: [e.g., Ubuntu 20.04]
- Compiler: [e.g., GCC 9.3.0]
- Architecture: [e.g., x86_64]

**Additional context**
Any other relevant information
```

### 2. Suggesting Enhancements

- Open an issue with the "enhancement" label
- Clearly describe the feature and its benefits
- Provide use cases and examples

### 3. Code Contributions

#### Development Setup

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/FFT-implementation-in-C.git
   cd FFT-implementation-in-C
   ```

3. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

4. Make your changes following the coding standards

5. Build and test:
   ```bash
   make clean all
   make tests
   ./bin/test_all
   ```

#### Coding Standards

**C Style Guide:**

```c
/**
 * @brief Function description
 * 
 * @param x Parameter description
 * @return Return value description
 */
return_type function_name(param_type param) {
    // 4-space indentation
    if (condition) {
        // Braces on same line for statements
        do_something();
    } else {
        do_something_else();
    }
    
    // Meaningful variable names
    int sample_count = get_sample_count();
    
    // Constants in UPPER_CASE
    const int MAX_FFT_SIZE = 65536;
    
    return result;
}
```

**File Organization:**
- Header files in `include/`
- Implementation files in appropriate subdirectories
- Keep related functions together
- One algorithm per file

**Documentation:**
- Use Doxygen-style comments for all public functions
- Include mathematical references where applicable
- Provide usage examples in comments

#### Testing Requirements

All contributions must include appropriate tests:

1. **Unit Tests**: Test individual functions
2. **Integration Tests**: Test algorithm combinations
3. **Performance Tests**: Benchmark if claiming improvements

Example test:

```c
void test_fft_properties() {
    // Parseval's theorem test
    complex_t* signal = generate_random_signal(1024);
    complex_t* spectrum = allocate_complex_array(1024);
    memcpy(spectrum, signal, 1024 * sizeof(complex_t));
    
    double time_energy = compute_energy(signal, 1024);
    radix2_dit_fft(spectrum, 1024, FFT_FORWARD);
    double freq_energy = compute_energy(spectrum, 1024) / 1024;
    
    assert(fabs(time_energy - freq_energy) < 1e-10);
}
```

### 4. Pull Request Process

1. Ensure your code follows the style guide
2. Add tests for new functionality
3. Update documentation as needed
4. Ensure all tests pass
5. Update the README.md if needed
6. Submit PR with clear description

#### PR Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] All tests pass
- [ ] Added new tests
- [ ] Tested on multiple platforms

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No new warnings
```

## Areas for Contribution

### Priority Areas

1. **New Algorithms**
   - Prime Factor Algorithm
   - Winograd FFT
   - Number Theoretic Transform

2. **Optimizations**
   - GPU implementations (CUDA/OpenCL)
   - ARM NEON optimizations
   - AVX-512 support

3. **Applications**
   - Real-time audio effects
   - Radar signal processing
   - Medical imaging (MRI reconstruction)

4. **Tools and Utilities**
   - FFT plan generator
   - Automatic algorithm selection
   - Visualization tools

5. **Documentation**
   - Tutorial notebooks
   - Video explanations
   - More examples

### Good First Issues

Look for issues labeled "good first issue":
- Fix compiler warnings
- Improve error messages
- Add input validation
- Write missing tests
- Improve documentation

## Development Guidelines

### 1. Algorithm Implementation

When implementing a new FFT algorithm:

```c
// Template for new algorithm
#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"

/**
 * @file algorithm_name.c
 * @brief One-line description
 * 
 * @details
 * Detailed explanation including:
 * - Mathematical background
 * - Algorithm