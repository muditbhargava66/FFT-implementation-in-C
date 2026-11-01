#ifdef __APPLE__

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#import <MetalPerformanceShadersGraph/MetalPerformanceShadersGraph.h>
#include "../include/fft_gpu.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/**
 * @file fft_metal.m
 * @brief Metal backend for the GPU interface, built on MPSGraph
 *
 * MPSGraph exposes a complex-valued fast Fourier transform from macOS 14.
 * Apple's metal-cpp wrappers cover Metal, Foundation, QuartzCore and MetalFX
 * but contain no FFT, so the transform comes from MetalPerformanceShadersGraph.
 * The transform runs in float32, so host double-precision complex_t values are
 * narrowed on upload and widened again on download. Forward transforms are
 * unscaled and inverse transforms carry 1/N, matching the CPU kernels.
 *
 * The real and imaginary parts are stored in separate MTLBuffers. Keeping them
 * planar means the buffers can be fed to MPSGraph without a repack, and the
 * result can be read straight back into them.
 */

#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 140000

@interface FFTMetalPlan : NSObject
@property(nonatomic, strong) MPSGraph *graph;
@property(nonatomic, strong) MPSGraphTensor *realIn;
@property(nonatomic, strong) MPSGraphTensor *imagIn;
@property(nonatomic, strong) MPSGraphTensor *realOut;
@property(nonatomic, strong) MPSGraphTensor *imagOut;
@property(nonatomic, assign) int rows;
@property(nonatomic, assign) int cols;
@property(nonatomic, assign) size_t elements;
@end

@implementation FFTMetalPlan
@end

struct fft_gpu_memory {
    void *realBuffer; /* +1 MTLBuffer, one float per element */
    void *imagBuffer; /* +1 MTLBuffer, one float per element */
    size_t elements;
};

struct fft_gpu_plan {
    void *state; /* +1 FFTMetalPlan */
};

static id<MTLDevice> g_device = nil;
static id<MTLCommandQueue> g_command_queue = nil;
static int g_metal_initialized = 0;
static char g_device_name[256] = "No Metal device";

static id<MTLBuffer> memory_real(fft_gpu_memory_t mem) {
    return (__bridge id<MTLBuffer>)mem->realBuffer;
}

static id<MTLBuffer> memory_imag(fft_gpu_memory_t mem) {
    return (__bridge id<MTLBuffer>)mem->imagBuffer;
}

static int metal_available(void) {
    if (@available(macOS 14.0, *)) {
        return 1;
    }
    return 0;
}

int fft_gpu_init_metal(void) {
    if (g_metal_initialized) return 0;
    if (!metal_available()) return -1;

    @autoreleasepool {
        g_device = MTLCreateSystemDefaultDevice();
        if (!g_device) return -1;

        g_command_queue = [g_device newCommandQueue];
        if (!g_command_queue) {
            g_device = nil;
            return -1;
        }

        const char *name = [[g_device name] UTF8String];
        strncpy(g_device_name, name ? name : "Metal device",
                sizeof(g_device_name) - 1);
        g_device_name[sizeof(g_device_name) - 1] = '\0';
        g_metal_initialized = 1;
    }
    return 0;
}

void fft_gpu_cleanup_metal(void) {
    if (!g_metal_initialized) return;

    @autoreleasepool {
        g_command_queue = nil;
        g_device = nil;
        g_metal_initialized = 0;
        strncpy(g_device_name, "No Metal device", sizeof(g_device_name) - 1);
    }
}

int fft_gpu_available_metal(void) {
    if (!metal_available()) return 0;

    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return device != nil ? 1 : 0;
    }
}

fft_gpu_memory_t fft_gpu_alloc_metal(size_t n) {
    if (n == 0 || n > SIZE_MAX / sizeof(float)) return NULL;
    if (fft_gpu_init_metal() != 0) return NULL;

    fft_gpu_memory_t mem = calloc(1, sizeof(*mem));
    if (!mem) return NULL;

    @autoreleasepool {
        id<MTLBuffer> real = [g_device newBufferWithLength:n * sizeof(float)
                                                  options:MTLResourceStorageModeShared];
        id<MTLBuffer> imag = [g_device newBufferWithLength:n * sizeof(float)
                                                  options:MTLResourceStorageModeShared];
        if (!real || !imag) {
            free(mem);
            return NULL;
        }
        mem->realBuffer = (void *)CFBridgingRetain(real);
        mem->imagBuffer = (void *)CFBridgingRetain(imag);
        mem->elements = n;
    }
    return mem;
}

void fft_gpu_free_metal(fft_gpu_memory_t mem) {
    if (!mem) return;

    if (mem->realBuffer) CFBridgingRelease(mem->realBuffer);
    if (mem->imagBuffer) CFBridgingRelease(mem->imagBuffer);
    free(mem);
}

void fft_gpu_copy_h2d_metal(fft_gpu_memory_t dst, const complex_t *src,
                            size_t n) {
    if (!dst || !src) return;

    size_t count = n < dst->elements ? n : dst->elements;
    float *re = (float *)[memory_real(dst) contents];
    float *im = (float *)[memory_imag(dst) contents];
    for (size_t i = 0; i < count; i++) {
        re[i] = (float)creal(src[i]);
        im[i] = (float)cimag(src[i]);
    }
}

void fft_gpu_copy_d2h_metal(complex_t *dst, fft_gpu_memory_t src, size_t n) {
    if (!dst || !src) return;

    size_t count = n < src->elements ? n : src->elements;
    const float *re = (const float *)[memory_real(src) contents];
    const float *im = (const float *)[memory_imag(src) contents];
    for (size_t i = 0; i < count; i++) {
        dst[i] = (double)re[i] + I * (double)im[i];
    }
}

/* Build the graph for a rank-2 shape. A 1-D batch is [batch, n] transformed
 * along the last axis; a 2-D transform is [rows, cols] along both axes. */
static FFTMetalPlan *create_plan_state(int rows, int cols, BOOL is2D,
                                       fft_direction direction) {
    if (@available(macOS 14.0, *)) {
        @autoreleasepool {
            MPSGraph *graph = [[MPSGraph alloc] init];
            if (!graph) return nil;

            MPSShape *shape = @[ @(rows), @(cols) ];
            MPSGraphTensor *realIn =
                [graph placeholderWithShape:shape
                                   dataType:MPSDataTypeFloat32
                                       name:@"real"];
            MPSGraphTensor *imagIn =
                [graph placeholderWithShape:shape
                                   dataType:MPSDataTypeFloat32
                                       name:@"imag"];
            MPSGraphTensor *complexIn =
                [graph complexTensorWithRealTensor:realIn
                                   imaginaryTensor:imagIn
                                              name:@"complex"];
            if (!realIn || !imagIn || !complexIn) return nil;

            MPSGraphFFTDescriptor *descriptor = [MPSGraphFFTDescriptor descriptor];
            if (!descriptor) return nil;
            descriptor.inverse = (direction == FFT_INVERSE);
            /* The CPU inverse kernels divide by the transform length. For a
             * 2-D plan MPSGraph's "size" mode divides by rows*cols. */
            descriptor.scalingMode = (direction == FFT_INVERSE)
                                         ? MPSGraphFFTScalingModeSize
                                         : MPSGraphFFTScalingModeNone;

            NSArray<NSNumber *> *axes = is2D ? @[ @0, @1 ] : @[ @1 ];
            MPSGraphTensor *spectrum =
                [graph fastFourierTransformWithTensor:complexIn
                                                 axes:axes
                                           descriptor:descriptor
                                                 name:@"fft"];
            if (!spectrum) return nil;

            MPSGraphTensor *realOut =
                [graph realPartOfTensor:spectrum name:@"realOut"];
            MPSGraphTensor *imagOut =
                [graph imaginaryPartOfTensor:spectrum name:@"imagOut"];
            if (!realOut || !imagOut) return nil;

            FFTMetalPlan *state = [[FFTMetalPlan alloc] init];
            state.graph = graph;
            state.realIn = realIn;
            state.imagIn = imagIn;
            state.realOut = realOut;
            state.imagOut = imagOut;
            state.rows = rows;
            state.cols = cols;
            state.elements = (size_t)rows * (size_t)cols;
            return state;
        }
    }
    return nil;
}

static fft_gpu_plan_t store_plan(FFTMetalPlan *state) {
    if (!state) return NULL;

    fft_gpu_plan_t plan = calloc(1, sizeof(*plan));
    if (!plan) return NULL;

    plan->state = (void *)CFBridgingRetain(state);
    return plan;
}

fft_gpu_plan_t fft_gpu_plan_1d_metal(int n, int batch, fft_direction dir) {
    if (n <= 0 || batch <= 0) return NULL;
    if (dir != FFT_FORWARD && dir != FFT_INVERSE) return NULL;
    if ((size_t)batch > SIZE_MAX / (size_t)n) return NULL;
    if (fft_gpu_init_metal() != 0) return NULL;

    return store_plan(create_plan_state(batch, n, NO, dir));
}

fft_gpu_plan_t fft_gpu_plan_2d_metal(int rows, int cols, fft_direction dir) {
    if (rows <= 0 || cols <= 0) return NULL;
    if (dir != FFT_FORWARD && dir != FFT_INVERSE) return NULL;
    if ((size_t)cols > SIZE_MAX / (size_t)rows) return NULL;
    if (fft_gpu_init_metal() != 0) return NULL;

    return store_plan(create_plan_state(rows, cols, YES, dir));
}

void fft_gpu_destroy_plan_metal(fft_gpu_plan_t plan) {
    if (!plan) return;

    if (plan->state) CFBridgingRelease(plan->state);
    free(plan);
}

static int execute_plan(FFTMetalPlan *state, fft_gpu_memory_t in,
                        fft_gpu_memory_t out) {
    if (!state || !in || !out) return -1;
    if (in->elements < state.elements || out->elements < state.elements)
        return -1;
    if (!g_metal_initialized || !g_command_queue) return -1;

    @autoreleasepool {
        MPSShape *shape = @[ @(state.rows), @(state.cols) ];
        MPSGraphTensorData *realIn =
            [[MPSGraphTensorData alloc] initWithMTLBuffer:memory_real(in)
                                                   shape:shape
                                                dataType:MPSDataTypeFloat32];
        MPSGraphTensorData *imagIn =
            [[MPSGraphTensorData alloc] initWithMTLBuffer:memory_imag(in)
                                                   shape:shape
                                                dataType:MPSDataTypeFloat32];
        if (!realIn || !imagIn) return -1;

        NSDictionary<MPSGraphTensor *, MPSGraphTensorData *> *results =
            [state.graph runWithMTLCommandQueue:g_command_queue
                                          feeds:@{state.realIn : realIn,
                                                  state.imagIn : imagIn}
                                  targetTensors:@[ state.realOut, state.imagOut ]
                               targetOperations:nil];
        if (!results) return -1;

        MPSGraphTensorData *realOut = results[state.realOut];
        MPSGraphTensorData *imagOut = results[state.imagOut];
        if (!realOut || !imagOut) return -1;

        MPSNDArray *realArray = [realOut mpsndarray];
        MPSNDArray *imagArray = [imagOut mpsndarray];
        if (!realArray || !imagArray) return -1;

        /* On a discrete GPU the result may still live in private storage.
         * The synchronize is a no-op on Apple Silicon's unified memory. */
        id<MTLCommandBuffer> commandBuffer = [g_command_queue commandBuffer];
        [realArray synchronizeOnCommandBuffer:commandBuffer];
        [imagArray synchronizeOnCommandBuffer:commandBuffer];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        [realArray readBytes:[memory_real(out) contents] strideBytes:NULL];
        [imagArray readBytes:[memory_imag(out) contents] strideBytes:NULL];
    }
    return 0;
}

int fft_gpu_execute_metal(fft_gpu_plan_t plan, fft_gpu_memory_t in,
                          fft_gpu_memory_t out, fft_direction dir) {
    /* The plan fixes the direction; the dispatch layer has no per-plan query
     * and always passes FFT_FORWARD. */
    (void)dir;
    if (!plan || !plan->state) return -1;
    return execute_plan((__bridge FFTMetalPlan *)plan->state, in, out);
}

const char *fft_gpu_get_device_name_metal(void) { return g_device_name; }

void fft_gpu_get_memory_info_metal(size_t *total, size_t *available) {
    if (!total || !available) return;

    *total = 0;
    *available = 0;
    if (!g_metal_initialized || !g_device) return;

    @autoreleasepool {
        size_t budget = (size_t)[g_device recommendedMaxWorkingSetSize];
        size_t used = (size_t)[g_device currentAllocatedSize];
        *total = budget;
        *available = budget > used ? budget - used : 0;
    }
}

int fft_gpu_set_device_metal(int device) {
    /* Every macOS build exposes one system-default device. Selecting another
     * would invalidate buffers that plans and handles already hold. */
    return device == 0 ? 0 : -1;
}

static int run_1d(complex_t *in, complex_t *out, int n, int batch,
                  fft_direction direction) {
    if (!in || !out || n <= 0 || batch <= 0) return -1;

    fft_gpu_plan_t plan = fft_gpu_plan_1d_metal(n, batch, direction);
    if (!plan) return -1;

    size_t elements = (size_t)n * (size_t)batch;
    fft_gpu_memory_t src = fft_gpu_alloc_metal(elements);
    fft_gpu_memory_t dst = fft_gpu_alloc_metal(elements);
    int rc = -1;
    if (src && dst) {
        fft_gpu_copy_h2d_metal(src, in, elements);
        rc = execute_plan((__bridge FFTMetalPlan *)plan->state, src, dst);
        if (rc == 0) fft_gpu_copy_d2h_metal(out, dst, elements);
    }

    fft_gpu_free_metal(src);
    fft_gpu_free_metal(dst);
    fft_gpu_destroy_plan_metal(plan);
    return rc;
}

int fft_gpu_dft_1d_metal(complex_t *in, complex_t *out, int n,
                         fft_direction dir) {
    return run_1d(in, out, n, 1, dir);
}

int fft_gpu_dft_1d_batch_metal(complex_t *in, complex_t *out, int n, int batch,
                               fft_direction dir) {
    return run_1d(in, out, n, batch, dir);
}

int fft_gpu_dft_2d_metal(complex_t *in, complex_t *out, int rows, int cols,
                         fft_direction dir) {
    if (!in || !out || rows <= 0 || cols <= 0) return -1;

    fft_gpu_plan_t plan = fft_gpu_plan_2d_metal(rows, cols, dir);
    if (!plan) return -1;

    size_t elements = (size_t)rows * (size_t)cols;
    fft_gpu_memory_t src = fft_gpu_alloc_metal(elements);
    fft_gpu_memory_t dst = fft_gpu_alloc_metal(elements);
    int rc = -1;
    if (src && dst) {
        fft_gpu_copy_h2d_metal(src, in, elements);
        rc = execute_plan((__bridge FFTMetalPlan *)plan->state, src, dst);
        if (rc == 0) fft_gpu_copy_d2h_metal(out, dst, elements);
    }

    fft_gpu_free_metal(src);
    fft_gpu_free_metal(dst);
    fft_gpu_destroy_plan_metal(plan);
    return rc;
}

#else /* __MAC_OS_X_VERSION_MAX_ALLOWED < 140000 */

struct fft_gpu_memory;
struct fft_gpu_plan;

int fft_gpu_init_metal(void) { return -1; }
void fft_gpu_cleanup_metal(void) {}
int fft_gpu_available_metal(void) { return 0; }
fft_gpu_memory_t fft_gpu_alloc_metal(size_t n) {
    (void)n;
    return NULL;
}
void fft_gpu_free_metal(fft_gpu_memory_t mem) { (void)mem; }
void fft_gpu_copy_h2d_metal(fft_gpu_memory_t dst, const complex_t *src,
                            size_t n) {
    (void)dst;
    (void)src;
    (void)n;
}
void fft_gpu_copy_d2h_metal(complex_t *dst, fft_gpu_memory_t src, size_t n) {
    (void)dst;
    (void)src;
    (void)n;
}
fft_gpu_plan_t fft_gpu_plan_1d_metal(int n, int batch, fft_direction dir) {
    (void)n;
    (void)batch;
    (void)dir;
    return NULL;
}
fft_gpu_plan_t fft_gpu_plan_2d_metal(int rows, int cols, fft_direction dir) {
    (void)rows;
    (void)cols;
    (void)dir;
    return NULL;
}
int fft_gpu_execute_metal(fft_gpu_plan_t plan, fft_gpu_memory_t in,
                          fft_gpu_memory_t out, fft_direction dir) {
    (void)plan;
    (void)in;
    (void)out;
    (void)dir;
    return -1;
}
void fft_gpu_destroy_plan_metal(fft_gpu_plan_t plan) { (void)plan; }
const char *fft_gpu_get_device_name_metal(void) {
    return "Metal FFT requires macOS 14 or newer";
}
void fft_gpu_get_memory_info_metal(size_t *total, size_t *available) {
    if (total) *total = 0;
    if (available) *available = 0;
}
int fft_gpu_set_device_metal(int device) {
    (void)device;
    return -1;
}
int fft_gpu_dft_1d_metal(complex_t *in, complex_t *out, int n,
                         fft_direction dir) {
    (void)in;
    (void)out;
    (void)n;
    (void)dir;
    return -1;
}
int fft_gpu_dft_1d_batch_metal(complex_t *in, complex_t *out, int n, int batch,
                               fft_direction dir) {
    (void)in;
    (void)out;
    (void)n;
    (void)batch;
    (void)dir;
    return -1;
}
int fft_gpu_dft_2d_metal(complex_t *in, complex_t *out, int rows, int cols,
                         fft_direction dir) {
    (void)in;
    (void)out;
    (void)rows;
    (void)cols;
    (void)dir;
    return -1;
}

#endif /* __MAC_OS_X_VERSION_MAX_ALLOWED >= 140000 */

#endif /* __APPLE__ */
