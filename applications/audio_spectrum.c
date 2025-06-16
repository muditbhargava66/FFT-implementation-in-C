#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <stdint.h>

/**
 * Audio Spectrum Analyzer
 * 
 * Analyzes audio signals using FFT to display frequency content.
 * Supports WAV file reading and real-time spectrum analysis.
 * 
 * Features:
 * - WAV file parsing
 * - Windowing functions
 * - Frequency bin calculation
 * - Magnitude spectrum display
 * - Peak frequency detection
 */

// WAV file header structure
typedef struct {
    char riff[4];           // "RIFF"
    uint32_t size;          // File size - 8
    char wave[4];           // "WAVE"
    char fmt[4];            // "fmt "
    uint32_t fmt_size;      // Format chunk size
    uint16_t format;        // Audio format (1 = PCM)
    uint16_t channels;      // Number of channels
    uint32_t sample_rate;   // Sample rate
    uint32_t byte_rate;     // Byte rate
    uint16_t block_align;   // Block align
    uint16_t bits_per_sample; // Bits per sample
    char data[4];           // "data"
    uint32_t data_size;     // Data size
} wav_header_t;

// Window functions for spectral analysis
void apply_window_hann(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.5 * (1.0 - cos(TWO_PI * i / (n - 1)));
        signal[i] *= window;
    }
}

void apply_window_hamming(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.54 - 0.46 * cos(TWO_PI * i / (n - 1));
        signal[i] *= window;
    }
}

void apply_window_blackman(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.42 - 0.5 * cos(TWO_PI * i / (n - 1)) 
                      + 0.08 * cos(4 * PI * i / (n - 1));
        signal[i] *= window;
    }
}

// Generate test audio signal
void generate_test_audio(complex_t* signal, int n, double sample_rate) {
    // Generate a signal with multiple frequency components
    double f1 = 440.0;   // A4 note
    double f2 = 554.37;  // C#5 note
    double f3 = 659.25;  // E5 note
    
    for (int i = 0; i < n; i++) {
        double t = i / sample_rate;
        signal[i] = 0.5 * sin(TWO_PI * f1 * t) +
                   0.3 * sin(TWO_PI * f2 * t) +
                   0.2 * sin(TWO_PI * f3 * t) +
                   0.1 * ((double)rand() / RAND_MAX - 0.5);  // Add noise
    }
}

// Compute frequency for each FFT bin
double bin_to_frequency(int bin, int fft_size, double sample_rate) {
    return bin * sample_rate / fft_size;
}

// Find peak frequencies in spectrum
typedef struct {
    double frequency;
    double magnitude;
    int bin;
} peak_t;

void find_peaks(double* magnitude, int n, double sample_rate, 
                peak_t* peaks, int* num_peaks, int max_peaks) {
    *num_peaks = 0;
    double threshold = 0.1;  // Minimum magnitude threshold
    
    for (int i = 1; i < n/2 - 1 && *num_peaks < max_peaks; i++) {
        // Check if local maximum
        if (magnitude[i] > magnitude[i-1] && 
            magnitude[i] > magnitude[i+1] &&
            magnitude[i] > threshold) {
            
            peaks[*num_peaks].bin = i;
            peaks[*num_peaks].frequency = bin_to_frequency(i, n, sample_rate);
            peaks[*num_peaks].magnitude = magnitude[i];
            (*num_peaks)++;
        }
    }
    
    // Sort peaks by magnitude (simple bubble sort)
    for (int i = 0; i < *num_peaks - 1; i++) {
        for (int j = 0; j < *num_peaks - i - 1; j++) {
            if (peaks[j].magnitude < peaks[j+1].magnitude) {
                peak_t temp = peaks[j];
                peaks[j] = peaks[j+1];
                peaks[j+1] = temp;
            }
        }
    }
}

// Display spectrum as ASCII art
void display_spectrum_ascii(double* magnitude, int n, double sample_rate) {
    int display_bins = 64;  // Number of bins to display
    if (display_bins > n/2) display_bins = n/2;
    
    // Find maximum magnitude for scaling
    double max_mag = 0;
    for (int i = 0; i < display_bins; i++) {
        if (magnitude[i] > max_mag) max_mag = magnitude[i];
    }
    
    printf("\nFrequency Spectrum:\n");
    printf("==================\n");
    
    // Display each frequency bin
    for (int i = 0; i < display_bins; i++) {
        double freq = bin_to_frequency(i, n, sample_rate);
        int bar_length = (int)(50 * magnitude[i] / max_mag);
        
        printf("%5.0f Hz |", freq);
        for (int j = 0; j < bar_length; j++) {
            printf("█");
        }
        printf(" %.3f\n", magnitude[i]);
    }
}

// Analyze audio spectrum
void analyze_audio_spectrum(complex_t* signal, int n, double sample_rate, 
                          const char* window_type) {
    // Apply window function
    if (strcmp(window_type, "hann") == 0) {
        apply_window_hann(signal, n);
    } else if (strcmp(window_type, "hamming") == 0) {
        apply_window_hamming(signal, n);
    } else if (strcmp(window_type, "blackman") == 0) {
        apply_window_blackman(signal, n);
    }
    
    // Compute FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Compute magnitude spectrum
    double* magnitude = compute_magnitude(signal, n);
    
    // Find and display peaks
    peak_t peaks[10];
    int num_peaks;
    find_peaks(magnitude, n, sample_rate, peaks, &num_peaks, 10);
    
    printf("\nDetected Peaks:\n");
    printf("===============\n");
    for (int i = 0; i < num_peaks && i < 5; i++) {
        printf("%d. %.1f Hz (magnitude: %.3f)\n", 
               i+1, peaks[i].frequency, peaks[i].magnitude);
    }
    
    // Display spectrum
    display_spectrum_ascii(magnitude, n, sample_rate);
    
    free(magnitude);
}

// Musical note detection
const char* frequency_to_note(double freq) {
    static const char* notes[] = {"C", "C#", "D", "D#", "E", "F", 
                                  "F#", "G", "G#", "A", "A#", "B"};
    
    // A4 = 440 Hz
    double a4 = 440.0;
    double c0 = a4 * pow(2, -4.75);  // C0 frequency
    
    if (freq <= 0) return "N/A";
    
    double half_steps = 12 * log2(freq / c0);
    int note_index = ((int)round(half_steps)) % 12;
    int octave = ((int)round(half_steps)) / 12;
    
    static char note_name[10];
    snprintf(note_name, sizeof(note_name), "%s%d", notes[note_index], octave);
    return note_name;
}

// Main demonstration
int main() {
    printf("Audio Spectrum Analyzer\n");
    printf("=======================\n");
    
    // Parameters
    double sample_rate = 44100.0;  // CD quality
    int fft_size = 4096;          // FFT size
    
    // Test 1: Analyze synthetic audio
    printf("\nTest 1: Synthetic Audio Signal\n");
    printf("------------------------------\n");
    
    complex_t* signal = allocate_complex_array(fft_size);
    generate_test_audio(signal, fft_size, sample_rate);
    
    printf("Generated signal with frequencies: 440 Hz (A4), 554 Hz (C#5), 659 Hz (E5)\n");
    
    // Analyze with different windows
    const char* windows[] = {"none", "hann", "hamming", "blackman"};
    
    for (int w = 0; w < 4; w++) {
        // Make a copy for analysis
        complex_t* signal_copy = allocate_complex_array(fft_size);
        memcpy(signal_copy, signal, fft_size * sizeof(complex_t));
        
        printf("\n\nWindow: %s\n", windows[w]);
        analyze_audio_spectrum(signal_copy, fft_size, sample_rate, windows[w]);
        
        free_complex_array(signal_copy);
    }
    
    // Test 2: Musical chord analysis
    printf("\n\nTest 2: Musical Chord Analysis\n");
    printf("-------------------------------\n");
    
    // Generate C major chord (C4, E4, G4)
    double chord_freqs[] = {261.63, 329.63, 392.00};  // C4, E4, G4
    
    for (int i = 0; i < fft_size; i++) {
        double t = i / sample_rate;
        signal[i] = 0;
        for (int j = 0; j < 3; j++) {
            signal[i] += 0.33 * sin(TWO_PI * chord_freqs[j] * t);
        }
    }
    
    printf("C Major Chord (C4, E4, G4):\n");
    analyze_audio_spectrum(signal, fft_size, sample_rate, "blackman");
    
    // Show note detection
    printf("\nNote Detection:\n");
    radix2_dit_fft(signal, fft_size, FFT_FORWARD);
    double* magnitude = compute_magnitude(signal, fft_size);
    
    peak_t peaks[10];
    int num_peaks;
    find_peaks(magnitude, fft_size, sample_rate, peaks, &num_peaks, 10);
    
    for (int i = 0; i < num_peaks && i < 5; i++) {
        printf("Peak %d: %.1f Hz = %s\n", 
               i+1, peaks[i].frequency, frequency_to_note(peaks[i].frequency));
    }
    
    // Test 3: Frequency resolution
    printf("\n\nTest 3: Frequency Resolution\n");
    printf("-----------------------------\n");
    
    double freq_resolution = sample_rate / fft_size;
    printf("Sample rate: %.0f Hz\n", sample_rate);
    printf("FFT size: %d\n", fft_size);
    printf("Frequency resolution: %.2f Hz/bin\n", freq_resolution);
    printf("Nyquist frequency: %.0f Hz\n", sample_rate / 2);
    printf("Time window: %.3f seconds\n", fft_size / sample_rate);
    
    // Performance considerations
    printf("\n\nPerformance Considerations:\n");
    printf("===========================\n");
    printf("1. Use power-of-2 FFT sizes for efficiency\n");
    printf("2. Apply window functions to reduce spectral leakage\n");
    printf("3. Zero-pad for increased frequency resolution\n");
    printf("4. Use overlap for time-varying signals\n");
    printf("5. Consider real-FFT for real-valued audio\n");
    
    free_complex_array(signal);
    free(magnitude);
    
    return 0;
}
