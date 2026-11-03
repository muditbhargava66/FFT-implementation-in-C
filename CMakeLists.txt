cmake_minimum_required(VERSION 3.20)
project(FFT_Implementation_C VERSION 3.4.0 LANGUAGES C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# The kernels use C11 _Complex arithmetic, which the MSVC C compiler does not
# provide. Fail at configure time with the reason instead of partway through
# the build, and do it before any CUDA probing: check_language(CUDA) invokes
# nvcc, which needs cl.exe as its host compiler on Windows and otherwise fails
# with a less useful message. Gate on the compiler id because clang-cl sets
# MSVC and does support _Complex.
if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    message(FATAL_ERROR
        "MSVC's C compiler cannot build this library: it requires C11 _Complex "
        "arithmetic, which cl.exe does not provide. Configure with MinGW-w64 "
        "(for example `cmake -G \"MinGW Makefiles\" ..`) or with clang.")
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(include)

# -march=native is opt-in: it produces binaries that only run on the build CPU,
# and MinGW-w64 GCC 13.2.0 miscompiles 32-byte AVX stack stores at some
# alignment phases, which crashed realtime_analyzer intermittently.
option(FFT_NATIVE_ARCH "Tune for the build host CPU (-march=native); non-portable" OFF)
option(FFT_ENABLE_OPENMP "Build the threaded radix-4 kernel (OpenMP)" OFF)

# Compiler options. Language-scoped so nvcc does not receive GCC-only flags
# such as -march=native; CUDA gets its warnings through -Xcompiler.
add_compile_options(
    $<$<COMPILE_LANGUAGE:C>:-Wall>
    $<$<COMPILE_LANGUAGE:C>:-Wextra>
    $<$<COMPILE_LANGUAGE:C>:-O3>
    $<$<COMPILE_LANGUAGE:OBJC>:-Wall>
    $<$<COMPILE_LANGUAGE:OBJC>:-Wextra>
    $<$<COMPILE_LANGUAGE:OBJC>:-O3>
    $<$<COMPILE_LANGUAGE:CUDA>:-O3>
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra>
)
if(FFT_NATIVE_ARCH)
    message(WARNING
        "FFT_NATIVE_ARCH is ON: binaries will only run on CPUs compatible "
        "with the build host, and MinGW-w64 GCC may emit misaligned 32-byte "
        "AVX stack stores that crash intermittently.")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:C>:-march=native>
        $<$<COMPILE_LANGUAGE:OBJC>:-march=native>
    )
endif()
# Expose POSIX functions (posix_memalign, sysconf) in strict C11 mode
add_compile_definitions(_GNU_SOURCE _POSIX_C_SOURCE=200809L)

# CONFIGURE_DEPENDS makes the glob re-run when files are added or removed.
file(GLOB_RECURSE FFT_C_SOURCES CONFIGURE_DEPENDS
    "algorithms/core/*.c"
    "algorithms/dft/*.c"
    "algorithms/auto/*.c"
    "utils/*.c"
)

# ---- GPU backends -----------------------------------------------------------
# Both are opt-in. Exactly one implementation of the fft_gpu_* symbols is
# linked: utils/fft_gpu_stubs.c when neither backend is enabled, otherwise
# gpu/fft_gpu.c plus the selected backend.
option(FFT_ENABLE_CUDA
    "Build the CUDA/cuFFT GPU backend (Linux; use WSL on Windows)" OFF)
option(FFT_ENABLE_GPU
    "Deprecated alias for FFT_ENABLE_CUDA" OFF)
if(FFT_ENABLE_GPU AND NOT FFT_ENABLE_CUDA)
    message(WARNING "FFT_ENABLE_GPU is deprecated; enabling FFT_ENABLE_CUDA.")
    set(FFT_ENABLE_CUDA ON)
endif()

set(FFT_CUDA_ENABLED OFF)
if(FFT_ENABLE_CUDA)
    if(WIN32)
        message(FATAL_ERROR
            "The CUDA backend is not supported in native Windows builds. This "
            "library requires C11 _Complex, which MSVC does not provide, and "
            "nvcc on Windows requires MSVC as its host compiler. Use WSL or "
            "Linux instead. Note that CUDA 13 also dropped Pascal (sm_61), so "
            "a GTX 1050 needs CUDA 12.x in WSL.")
    endif()
    include(CheckLanguage)
    check_language(CUDA)
    if(NOT CMAKE_CUDA_COMPILER)
        message(FATAL_ERROR
            "FFT_ENABLE_CUDA=ON but no CUDA compiler was found. Install the "
            "CUDA toolkit (nvcc) or set CMAKE_CUDA_COMPILER to its path.")
    endif()
    # Choose a default architecture list before enable_language(), otherwise
    # CMake fills CMAKE_CUDA_ARCHITECTURES with its own default and this block
    # never runs. Parse nvcc's banner because the compiler version is not
    # available until the language is enabled.
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        execute_process(
            COMMAND "${CMAKE_CUDA_COMPILER}" --version
            OUTPUT_VARIABLE _fft_nvcc_version
            ERROR_QUIET)
        if(_fft_nvcc_version MATCHES "release ([0-9]+)")
            if(CMAKE_MATCH_1 GREATER_EQUAL 13)
                # CUDA 13 dropped Maxwell, Pascal and Volta.
                set(CMAKE_CUDA_ARCHITECTURES "75;80;90" CACHE STRING "CUDA architectures")
            else()
                # Includes Pascal (sm_61) for the GTX 1050 used in local testing.
                set(CMAKE_CUDA_ARCHITECTURES "61;70;75;80" CACHE STRING "CUDA architectures")
            endif()
        else()
            set(CMAKE_CUDA_ARCHITECTURES "61;70;75;80" CACHE STRING "CUDA architectures")
        endif()
    endif()
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS
        "CUDA backend enabled (architectures: ${CMAKE_CUDA_ARCHITECTURES})")
    set(FFT_CUDA_ENABLED ON)
endif()

# Metal is opt-in because it requires macOS 14 and pulls in
# MetalPerformanceShadersGraph.
option(FFT_ENABLE_METAL
    "Build the Metal/MPSGraph GPU backend (Apple Silicon, macOS 14+)" OFF)

if(FFT_ENABLE_METAL AND NOT APPLE)
    message(WARNING
        "FFT_ENABLE_METAL is ON but this is not macOS; the Metal backend will "
        "not be built and fft_gpu_* stays on the CPU stubs.")
endif()

if(APPLE AND FFT_ENABLE_METAL AND NOT CMAKE_C_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR
        "FFT_ENABLE_METAL requires Clang: gpu/fft_metal.m is Objective-C with "
        "automatic reference counting, which GCC does not implement.")
endif()

set(FFT_METAL_ENABLED OFF)
if(APPLE AND FFT_ENABLE_METAL)
    enable_language(OBJC)
    set_source_files_properties(gpu/fft_metal.m PROPERTIES
        COMPILE_OPTIONS "-fobjc-arc")
    message(STATUS "Metal backend enabled (MetalPerformanceShadersGraph)")
    set(FFT_METAL_ENABLED ON)
endif()

# Select exactly one implementation of the fft_gpu_* symbols.
if(FFT_CUDA_ENABLED OR FFT_METAL_ENABLED)
    list(FILTER FFT_C_SOURCES EXCLUDE REGEX "/utils/fft_gpu_stubs\\.c$")
    list(APPEND FFT_C_SOURCES "gpu/fft_gpu.c")
endif()
if(FFT_CUDA_ENABLED)
    list(APPEND FFT_C_SOURCES "gpu/fft_cuda.cu")
endif()
if(FFT_METAL_ENABLED)
    list(APPEND FFT_C_SOURCES "gpu/fft_metal.m")
endif()

add_library(fft_lib_c STATIC ${FFT_C_SOURCES})
target_compile_definitions(fft_lib_c PRIVATE LIB_BUILD)

if(FFT_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS C)
    target_link_libraries(fft_lib_c PUBLIC OpenMP::OpenMP_C)
    message(STATUS "OpenMP radix-4 kernel enabled")
endif()

if(FFT_CUDA_ENABLED)
    target_compile_definitions(fft_lib_c PUBLIC FFT_GPU_CUDA_ENABLED)
    target_link_libraries(fft_lib_c PUBLIC CUDA::cudart CUDA::cufft)
endif()

if(FFT_METAL_ENABLED)
    target_compile_definitions(fft_lib_c PUBLIC FFT_GPU_METAL_ENABLED)
    target_link_libraries(fft_lib_c PUBLIC
        "-framework Foundation"
        "-framework Metal"
        "-framework MetalPerformanceShaders"
        "-framework MetalPerformanceShadersGraph")
endif()

# libm is separate on Linux. PUBLIC so tests, benchmarks and applications
# inherit it; MinGW and macOS resolve these symbols from the C runtime.
if(UNIX AND NOT APPLE)
    target_link_libraries(fft_lib_c PUBLIC m)
endif()

add_library(fft_lib ALIAS fft_lib_c)

# Tests
enable_testing()

# Helper to add test executables
function(add_fft_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE fft_lib)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# Add test_all
add_executable(test_all tests/test_all.c)
target_link_libraries(test_all PRIVATE fft_lib)
add_test(NAME test_all COMMAND test_all)

# Utility functions: file I/O and the helpers that were previously compiled
# without a public declaration.
add_executable(test_utils tests/test_utils.c)
target_link_libraries(test_utils PRIVATE fft_lib)
add_test(NAME test_utils COMMAND test_utils)

if(FFT_ENABLE_OPENMP)
    add_executable(test_threads tests/test_threads.c)
    target_link_libraries(test_threads PRIVATE fft_lib)
    add_test(NAME test_threads COMMAND test_threads)
    set_tests_properties(test_threads PROPERTIES SKIP_RETURN_CODE 77)
endif()

if(FFT_METAL_ENABLED)
    add_executable(test_gpu_metal tests/test_gpu_metal.c)
    target_link_libraries(test_gpu_metal PRIVATE fft_lib)
    add_test(NAME test_gpu_metal COMMAND test_gpu_metal)
    set_tests_properties(test_gpu_metal PROPERTIES SKIP_RETURN_CODE 77)
endif()

if(FFT_CUDA_ENABLED)
    add_executable(test_gpu_cuda tests/test_gpu_cuda.c)
    target_link_libraries(test_gpu_cuda PRIVATE fft_lib)
    add_test(NAME test_gpu_cuda COMMAND test_gpu_cuda)
    set_tests_properties(test_gpu_cuda PROPERTIES SKIP_RETURN_CODE 77)
endif()

# Add benchmarks
add_executable(benchmark_all benchmarks/benchmark_all.c)
target_link_libraries(benchmark_all PRIVATE fft_lib)

# Built as a target so it cannot silently drift from the GPU API again.
add_executable(gpu_benchmark benchmarks/gpu_benchmark.c)
target_link_libraries(gpu_benchmark PRIVATE fft_lib)

# Add applications
file(GLOB APP_SOURCES "applications/*.c")
foreach(APP_SRC ${APP_SOURCES})
    get_filename_component(APP_NAME ${APP_SRC} NAME_WE)
    add_executable(${APP_NAME} ${APP_SRC})
    target_link_libraries(${APP_NAME} PRIVATE fft_lib)
endforeach()

# CLI tools (individual algorithm demos) are disabled in the CMake build.
# Each algorithm .c file defines both the algorithm functions and a main()
# guarded by #ifndef LIB_BUILD. Building standalone CLIs would cause duplicate
# symbols against the library, and the demos depend on each other. Use test_all
# or benchmark_all instead.
