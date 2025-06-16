# FFT Implementation Makefile
# Comprehensive build system for all FFT algorithms

CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native -ffast-math
LDFLAGS = -lm -lpthread
OPENMP_FLAGS = -fopenmp

# Directories
SRC_DIR = .
INC_DIR = include
ALGO_DIR = algorithms
APP_DIR = applications
OPT_DIR = optimizations
BENCH_DIR = benchmarks
TEST_DIR = tests
BIN_DIR = bin
OBJ_DIR = obj
LIB_DIR = lib

# Create directories
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR))
$(shell mkdir -p $(OBJ_DIR)/algorithms/core $(OBJ_DIR)/algorithms/dft)
$(shell mkdir -p $(OBJ_DIR)/applications $(OBJ_DIR)/optimizations)

# Include paths
INCLUDES = -I$(INC_DIR)

# Source files
CORE_SRCS = $(wildcard $(ALGO_DIR)/core/*.c)
DFT_SRCS = $(wildcard $(ALGO_DIR)/dft/*.c)
APP_SRCS = $(wildcard $(APP_DIR)/*.c)
OPT_SRCS = $(wildcard $(OPT_DIR)/*.c)
BENCH_SRCS = $(wildcard $(BENCH_DIR)/*.c)
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)

# Object files for library
RADIX2_DIT_OBJ = $(OBJ_DIR)/algorithms/core/radix2_dit_lib.o
RADIX2_DIF_OBJ = $(OBJ_DIR)/algorithms/core/radix2_dif_lib.o
RADIX4_OBJ = $(OBJ_DIR)/algorithms/core/radix4_lib.o
SPLIT_RADIX_OBJ = $(OBJ_DIR)/algorithms/core/split_radix_lib.o
BLUESTEIN_OBJ = $(OBJ_DIR)/algorithms/core/bluestein_lib.o
MIXED_RADIX_OBJ = $(OBJ_DIR)/algorithms/core/mixed_radix_lib.o
RECURSIVE_OBJ = $(OBJ_DIR)/algorithms/core/recursive_fft_lib.o
ITERATIVE_OBJ = $(OBJ_DIR)/algorithms/core/iterative_fft_lib.o
NAIVE_DFT_OBJ = $(OBJ_DIR)/algorithms/dft/naive_dft_lib.o
OPT_DFT_OBJ = $(OBJ_DIR)/algorithms/dft/optimized_dft_lib.o

# All library objects
LIB_OBJS = $(RADIX2_DIT_OBJ) $(RADIX2_DIF_OBJ) $(RADIX4_OBJ) $(SPLIT_RADIX_OBJ) \
           $(BLUESTEIN_OBJ) $(MIXED_RADIX_OBJ) $(RECURSIVE_OBJ) $(ITERATIVE_OBJ) \
           $(NAIVE_DFT_OBJ) $(OPT_DFT_OBJ)

# Static library
STATIC_LIB = $(LIB_DIR)/libfft.a

# Executables
CORE_BINS = $(CORE_SRCS:$(ALGO_DIR)/core/%.c=$(BIN_DIR)/%)
DFT_BINS = $(DFT_SRCS:$(ALGO_DIR)/dft/%.c=$(BIN_DIR)/%)
APP_BINS = $(APP_SRCS:$(APP_DIR)/%.c=$(BIN_DIR)/%)
OPT_BINS = $(OPT_SRCS:$(OPT_DIR)/%.c=$(BIN_DIR)/%)
BENCH_BIN = $(BIN_DIR)/benchmark_all
TEST_BIN = $(BIN_DIR)/test_all

# Default target
all: $(STATIC_LIB) algorithms applications optimizations benchmarks tests

# Build static library
$(STATIC_LIB): $(LIB_OBJS)
	ar rcs $@ $^

# Algorithm library objects (compile without main)
$(OBJ_DIR)/algorithms/core/%_lib.o: $(ALGO_DIR)/core/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

$(OBJ_DIR)/algorithms/dft/%_lib.o: $(ALGO_DIR)/dft/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

# Algorithm targets
algorithms: $(STATIC_LIB) core dft

core: $(CORE_BINS)
	@echo "Built core FFT algorithms"

dft: $(DFT_BINS)
	@echo "Built DFT implementations"

# Application targets
applications: $(STATIC_LIB) $(APP_BINS)
	@echo "Built FFT applications"

# Optimization targets
optimizations: $(STATIC_LIB) $(OPT_BINS)
	@echo "Built optimized implementations"

# Benchmark target
benchmarks: $(STATIC_LIB) $(BENCH_BIN)
	@echo "Built benchmark suite"

# Test target
tests: $(STATIC_LIB) $(TEST_BIN)
	@echo "Built test suite"

# Build executables for algorithms (with library)
$(BIN_DIR)/%: $(ALGO_DIR)/core/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(ALGO_DIR)/dft/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(APP_DIR)/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(OPT_DIR)/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Special rules for parallel implementations
$(BIN_DIR)/parallel_fft: $(OPT_DIR)/parallel_fft.c $(STATIC_LIB)
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Benchmark compilation
$(BENCH_BIN): $(BENCH_DIR)/benchmark_all.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Test compilation
$(TEST_BIN): $(TEST_DIR)/test_all.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Run targets
run-benchmarks: $(BENCH_BIN)
	./$(BENCH_BIN)

run-tests: $(TEST_BIN)
	./$(TEST_BIN)

# Demo target - runs a selection of implementations
demo: $(BIN_DIR)/radix2_dit $(BIN_DIR)/audio_spectrum $(BIN_DIR)/convolution
	@echo "=== Radix-2 DIT FFT Demo ==="
	./$(BIN_DIR)/radix2_dit
	@echo "\n=== Audio Spectrum Analysis Demo ==="
	./$(BIN_DIR)/audio_spectrum
	@echo "\n=== Convolution Demo ==="
	./$(BIN_DIR)/convolution

# Performance profiling
profile: CFLAGS += -pg
profile: clean all
	@echo "Built with profiling enabled"

# Debug build
debug: CFLAGS = -Wall -Wextra -g -O0 -DDEBUG
debug: clean all
	@echo "Built with debug symbols"

# Static analysis
analyze:
	@echo "Running static analysis..."
	cppcheck --enable=all --suppress=missingIncludeSystem -I$(INC_DIR) $(SRC_DIR)

# Clean
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)
	rm -f gmon.out

# Install (copy binaries to system location)
install: all
	@echo "Installing FFT tools..."
	mkdir -p /usr/local/bin/fft-tools
	cp $(BIN_DIR)/* /usr/local/bin/fft-tools/
	@echo "Installation complete. Add /usr/local/bin/fft-tools to your PATH"

# Documentation
docs:
	@echo "Generating documentation..."
	doxygen Doxyfile

# Help
help:
	@echo "FFT Implementation Build System"
	@echo "=============================="
	@echo "Targets:"
	@echo "  all          - Build everything"
	@echo "  algorithms   - Build core FFT algorithms"
	@echo "  applications - Build FFT applications"
	@echo "  optimizations- Build optimized implementations"
	@echo "  benchmarks   - Build benchmark suite"
	@echo "  tests        - Build test suite"
	@echo "  demo         - Run demonstrations"
	@echo "  run-benchmarks - Run all benchmarks"
	@echo "  run-tests    - Run all tests"
	@echo "  profile      - Build with profiling"
	@echo "  debug        - Build with debug symbols"
	@echo "  analyze      - Run static analysis"
	@echo "  docs         - Generate documentation"
	@echo "  clean        - Remove build artifacts"
	@echo "  install      - Install tools system-wide"

.PHONY: all algorithms core dft applications optimizations benchmarks tests \
        demo run-benchmarks run-tests profile debug analyze docs clean install help
