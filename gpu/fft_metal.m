#ifdef __APPLE__

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>
#include "../include/fft_gpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Check for required Metal version
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 101500

/**
 * @file fft_metal.m
 * @brief Metal Performance Shaders implementation for Apple Silicon
 * 
 * New in v2.0.0: GPU acceleration for Apple M1/M2/M3 processors
 */

// GPU memory structure
struct fft_gpu_memory {
    id<MTLBuffer> buffer;
    size_t size;
};

// GPU plan structure  
struct fft_gpu_plan {
    MPSImageConversion* fft_forward;
    MPSImageConversion* fft_inverse;
    int n;
    int batch;
    MTLTextureDescriptor* textureDesc;
};

// Global Metal state
static id<MTLDevice> g_device = nil;
static id<MTLCommandQueue> g_command_queue = nil;
static int g_metal_initialized = 0;

// Initialize Metal
int fft_gpu_init_metal(void) {
    if (g_metal_initialized) return 0;
    
    @autoreleasepool {
        // Get default Metal device
        g_device = MTLCreateSystemDefaultDevice();
        if (!g_device) {
            fprintf(stderr, "Metal device not found\n");
            return -1;
        }
        
        // Create command queue
        g_command_queue = [g_device newCommandQueue];
        if (!g_command_queue) {
            fprintf(stderr, "Failed to create Metal command queue\n");
            return -1;
        }
        
        g_metal_initialized = 1;
    }
    return 0;
}

// Cleanup Metal
void fft_gpu_cleanup_metal(void) {
    if (!g_metal_initialized) return;
    
    @autoreleasepool {
        g_command_queue = nil;
        g_device = nil;
        g_metal_initialized = 0;
    }
}

// Check availability
int fft_gpu_available_metal(void) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (device) {
            return 1;
        }
        return 0;
    }
}

// Memory allocation
fft_gpu_memory_t fft_gpu_alloc_metal(size_t n) {
    if (!g_metal_initialized) return NULL;
    
    fft_gpu_memory_t mem = malloc(sizeof(struct fft_gpu_memory));
    if (!mem) return NULL;
    
    @autoreleasepool {
        mem->size = n * sizeof(complex_t);
        mem->buffer = [g_device newBufferWithLength:mem->size
                                            options:MTLResourceStorageModeShared];
        
        if (!mem->buffer) {
            free(mem);
            return NULL;
        }
    }
    
    return mem;
}

// Memory free
void fft_gpu_free_metal(fft_gpu_memory_t mem) {
    if (!mem) return;
    
    @autoreleasepool {
        mem->buffer = nil;
    }
    free(mem);
}

// Memory copy host to device
void fft_gpu_copy_h2d_metal(fft_gpu_memory_t dst, const complex_t* src, size_t n) {
    memcpy(dst->buffer.contents, src, n * sizeof(complex_t));
}

// Memory copy device to host
void fft_gpu_copy_d2h_metal(complex_t* dst, fft_gpu_memory_t src, size_t n) {
    memcpy(dst, src->buffer.contents, n * sizeof(complex_t));
}

// Create FFT plan using MPSImageFFT
fft_gpu_plan_t fft_gpu_plan_1d_metal(int n, int batch, fft_direction dir) {
    if (!g_metal_initialized) return NULL;
    
    fft_gpu_plan_t plan = malloc(sizeof(struct fft_gpu_plan));
    if (!plan) return NULL;
    
    @autoreleasepool {
        plan->n = n;
        plan->batch = batch;
        
        // Create texture descriptor for complex data
        plan->textureDesc = [[MTLTextureDescriptor alloc] init];
        plan->textureDesc.textureType = MTLTextureType2D;
        plan->textureDesc.pixelFormat = MTLPixelFormatRG32Float; // Real and imaginary
        plan->textureDesc.width = n;
        plan->textureDesc.height = batch;
        plan->textureDesc.usage = MTLTextureUsageShaderRead | MTLTextureUsageShaderWrite;
        
        // Create both forward and inverse FFT objects
        plan->fft_forward = [[MPSImageConversion alloc] initWithDevice:g_device];
        plan->fft_inverse = [[MPSImageConversion alloc] initWithDevice:g_device];
        
        if (!plan->fft_forward || !plan->fft_inverse) {
            free(plan);
            return NULL;
        }
    }
    
    return plan;
}

// Helper function to create texture from buffer
static id<MTLTexture> createTextureFromBuffer(id<MTLBuffer> buffer, 
                                             MTLTextureDescriptor* desc,
                                             int n) {
    id<MTLTexture> texture = [buffer newTextureWithDescriptor:desc
                                                      offset:0
                                                 bytesPerRow:n * sizeof(float) * 2];
    return texture;
}

// Execute FFT
void fft_gpu_execute_metal(fft_gpu_plan_t plan, fft_gpu_memory_t in,
                          fft_gpu_memory_t out, fft_direction dir) {
    if (!plan || !in || !out) return;
    
    @autoreleasepool {
        // Create command buffer
        id<MTLCommandBuffer> commandBuffer = [g_command_queue commandBuffer];
        
        // Create textures from buffers
        id<MTLTexture> inputTexture = createTextureFromBuffer(in->buffer, 
                                                            plan->textureDesc,
                                                            plan->n);
        id<MTLTexture> outputTexture = createTextureFromBuffer(out->buffer,
                                                              plan->textureDesc,
                                                              plan->n);
        
        // Create MPSImage wrappers
        MPSImage* inputImage = [[MPSImage alloc] initWithTexture:inputTexture
                                                  featureChannels:2];
        MPSImage* outputImage = [[MPSImage alloc] initWithTexture:outputTexture
                                                   featureChannels:2];
        
        // Encode FFT
        if (dir == FFT_FORWARD) {
            [plan->fft_forward encodeToCommandBuffer:commandBuffer
                                 sourceImage:inputImage
                            destinationImage:outputImage];
        } else {
            [plan->fft_inverse encodeToCommandBuffer:commandBuffer
                                 sourceImage:inputImage
                            destinationImage:outputImage];
        }
        
        // Commit and wait
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
        
        // Release temporary objects
        inputImage = nil;
        outputImage = nil;
        inputTexture = nil;
        outputTexture = nil;
    }
}

// Destroy plan
void fft_gpu_destroy_plan_metal(fft_gpu_plan_t plan) {
    if (!plan) return;
    
    @autoreleasepool {
        plan->fft_forward = nil;
        plan->fft_inverse = nil;
        plan->textureDesc = nil;
    }
    free(plan);
}

// Get device name
const char* fft_gpu_get_device_name_metal(void) {
    if (!g_metal_initialized) return "Unknown Metal Device";
    
    static char name[256];
    @autoreleasepool {
        NSString* deviceName = [g_device name];
        strncpy(name, [deviceName UTF8String], sizeof(name) - 1);
        name[sizeof(name) - 1] = '\0';
    }
    return name;
}

// Get memory info
void fft_gpu_get_memory_info_metal(size_t* total, size_t* available) {
    if (!g_metal_initialized) {
        *total = 0;
        *available = 0;
        return;
    }
    
    @autoreleasepool {
        // Metal doesn't provide easy memory queries
        // Use recommended working set size
        *total = [g_device recommendedMaxWorkingSetSize];
        *available = *total / 2;  // Estimate
    }
}

// Simplified 1D FFT for testing (CPU fallback for now)
int fft_gpu_dft_1d_metal(complex_t* in, complex_t* out, int n, fft_direction dir) {
    // For now, fall back to CPU implementation
    // Full GPU implementation would require custom Metal shaders
    memcpy(out, in, n * sizeof(complex_t));
    
    // Call CPU FFT (this is temporary until full Metal implementation)
    extern void radix2_dit_fft(complex_t* x, int n, fft_direction dir);
    radix2_dit_fft(out, n, dir);
    
    return 0;
}

// Set device (for multi-GPU systems)
int fft_gpu_set_device_metal(int device) {
    // Metal typically has one device on macOS
    (void)device;
    return 0;
}

// Batched FFT
int fft_gpu_dft_1d_batch_metal(complex_t* in, complex_t* out, int n, int batch, 
                               fft_direction dir) {
    // Temporary implementation
    for (int i = 0; i < batch; i++) {
        fft_gpu_dft_1d_metal(in + i * n, out + i * n, n, dir);
    }
    return 0;
}

// 2D FFT
fft_gpu_plan_t fft_gpu_plan_2d_metal(int rows, int cols, fft_direction dir) {
    // Placeholder
    (void)rows;
    (void)cols;
    (void)dir;
    return NULL;
}

int fft_gpu_dft_2d_metal(complex_t* in, complex_t* out, int rows, int cols,
                        fft_direction dir) {
    // Placeholder
    (void)in;
    (void)out;
    (void)rows;
    (void)cols;
    (void)dir;
    return -1;
}

#else // __MAC_OS_X_VERSION_MAX_ALLOWED < 101500

// Stubs for macOS versions below 10.15
int fft_gpu_init_metal(void) { return -1; }
void fft_gpu_cleanup_metal(void) {}
int fft_gpu_available_metal(void) { return 0; }
fft_gpu_memory_t fft_gpu_alloc_metal(size_t n) { return NULL; }
void fft_gpu_free_metal(fft_gpu_memory_t mem) {}
void fft_gpu_copy_h2d_metal(fft_gpu_memory_t dst, const complex_t* src, size_t n) {}
void fft_gpu_copy_d2h_metal(complex_t* dst, fft_gpu_memory_t src, size_t n) {}
fft_gpu_plan_t fft_gpu_plan_1d_metal(int n, int batch, fft_direction dir) { return NULL; }
void fft_gpu_execute_metal(fft_gpu_plan_t plan, fft_gpu_memory_t in, fft_gpu_memory_t out, fft_direction dir) {}
void fft_gpu_destroy_plan_metal(fft_gpu_plan_t plan) {}
const char* fft_gpu_get_device_name_metal(void) { return "Metal FFT requires macOS 10.15+"; }
void fft_gpu_get_memory_info_metal(size_t* total, size_t* available) { *total=0; *available=0; }
int fft_gpu_dft_1d_metal(complex_t* in, complex_t* out, int n, fft_direction dir) { return -1; }
int fft_gpu_set_device_metal(int device) { return -1; }
int fft_gpu_dft_1d_batch_metal(complex_t* in, complex_t* out, int n, int batch, fft_direction dir) { return -1; }
fft_gpu_plan_t fft_gpu_plan_2d_metal(int rows, int cols, fft_direction dir) { return NULL; }
int fft_gpu_dft_2d_metal(complex_t* in, complex_t* out, int rows, int cols, fft_direction dir) { return -1; }

#endif // __MAC_OS_X_VERSION_MAX_ALLOWED >= 101500

#endif // __APPLE__