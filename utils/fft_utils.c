#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <time.h>

/**
 * FFT Utility Functions
 * 
 * Common utilities for FFT operations including:
 * - Signal generation
 * - Window functions
 * - File I/O
 * - Visualization helpers
 * - Error analysis
 */

// Signal generators with parameters
void generate_chirp_signal(complex_t* signal, int n, double f_start, double f_end, double sample_rate) {
    double chirp_rate = (f_end - f_start) / (n / sample_rate);
    
    for (int i = 0; i < n; i++) {
        double t = i / sample_rate;
        double freq = f_start + chirp_rate * t;
        signal[i] = sin(TWO_PI * freq * t);
    }
}

void generate_noise(complex_t* signal, int n, double amplitude, unsigned int seed) {
    srand(seed);
    
    for (int i = 0; i < n; i++) {
        double real = amplitude * (2.0 * rand() / RAND_MAX - 1.0);
        double imag = amplitude * (2.0 * rand() / RAND_MAX - 1.0);
        signal[i] = real + I * imag;
    }
}

void generate_multi_tone(complex_t* signal, int n, double* frequencies, 
                        double* amplitudes, int num_tones, double sample_rate) {
    memset(signal, 0, n * sizeof(complex_t));
    
    for (int i = 0; i < n; i++) {
        for (int f = 0; f < num_tones; f++) {
            signal[i] += amplitudes[f] * sin(TWO_PI * frequencies[f] * i / sample_rate);
        }
    }
}

// Window functions
void apply_window_kaiser(complex_t* signal, int n, double beta) {
    double i0_beta = 1.0;  // Simplified - would compute modified Bessel function
    
    for (int i = 0; i < n; i++) {
        double x = 2.0 * i / (n - 1) - 1.0;
        double arg = beta * sqrt(1 - x * x);
        double window = 1.0;  // Simplified Kaiser window
        signal[i] *= window;
    }
}

void apply_window_tukey(complex_t* signal, int n, double alpha) {
    int transition = (int)(alpha * n / 2);
    
    for (int i = 0; i < n; i++) {
        double window = 1.0;
        
        if (i < transition) {
            window = 0.5 * (1 + cos(PI * (2.0 * i / (alpha * n) - 1)));
        } else if (i >= n - transition) {
            window = 0.5 * (1 + cos(PI * (2.0 * (n - 1 - i) / (alpha * n) - 1)));
        }
        
        signal[i] *= window;
    }
}

// File I/O functions
int save_complex_array(const char* filename, complex_t* data, int n) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for writing\n", filename);
        return -1;
    }
    
    fprintf(fp, "# FFT Data File\n");
    fprintf(fp, "# Format: index real imag magnitude phase\n");
    fprintf(fp, "# Size: %d\n", n);
    
    for (int i = 0; i < n; i++) {
        double mag = cabs(data[i]);
        double phase = carg(data[i]);
        fprintf(fp, "%d %e %e %e %e\n", 
                i, creal(data[i]), cimag(data[i]), mag, phase);
    }
    
    fclose(fp);
    return 0;
}

int load_complex_array(const char* filename, complex_t** data, int* n) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for reading\n", filename);
        return -1;
    }
    
    char line[256];
    *n = 0;
    
    // Read header and determine size
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] != '#') {
            (*n)++;
        } else if (strstr(line, "Size:")) {
            sscanf(line, "# Size: %d", n);
            break;
        }
    }
    
    if (*n == 0) {
        fclose(fp);
        return -1;
    }
    
    // Allocate array
    *data = allocate_complex_array(*n);
    
    // Read data
    rewind(fp);
    int idx = 0;
    while (fgets(line, sizeof(line), fp) && idx < *n) {
        if (line[0] != '#') {
            int i;
            double real, imag;
            if (sscanf(line, "%d %lf %lf", &i, &real, &imag) == 3) {
                (*data)[idx++] = real + I * imag;
            }
        }
    }
    
    fclose(fp);
    return 0;
}

// Spectrum analysis utilities
void find_spectral_peaks(complex_t* spectrum, int n, double sample_rate,
                        double* peak_freqs, double* peak_mags, int* num_peaks, int max_peaks) {
    double* magnitude = compute_magnitude(spectrum, n);
    *num_peaks = 0;
    
    // Find local maxima
    for (int i = 1; i < n/2 - 1 && *num_peaks < max_peaks; i++) {
        if (magnitude[i] > magnitude[i-1] && magnitude[i] > magnitude[i+1]) {
            // Quadratic interpolation for better frequency estimate
            double y1 = magnitude[i-1];
            double y2 = magnitude[i];
            double y3 = magnitude[i+1];
            
            double delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
            double freq = (i + delta) * sample_rate / n;
            
            peak_freqs[*num_peaks] = freq;
            peak_mags[*num_peaks] = y2;
            (*num_peaks)++;
        }
    }
    
    free(magnitude);
}

double compute_snr(complex_t* signal, complex_t* noise, int n) {
    double signal_power = 0;
    double noise_power = 0;
    
    for (int i = 0; i < n; i++) {
        signal_power += cabs(signal[i]) * cabs(signal[i]);
        noise_power += cabs(noise[i]) * cabs(noise[i]);
    }
    
    signal_power /= n;
    noise_power /= n;
    
    if (noise_power > 0) {
        return 10 * log10(signal_power / noise_power);
    }
    
    return INFINITY;
}

// Visualization helpers
void print_spectrum_plot(complex_t* spectrum, int n, double sample_rate, int width) {
    double* magnitude = compute_magnitude(spectrum, n);
    
    // Find max for scaling
    double max_mag = 0;
    for (int i = 0; i < n/2; i++) {
        if (magnitude[i] > max_mag) max_mag = magnitude[i];
    }
    
    printf("\nFrequency Spectrum:\n");
    printf("Freq (Hz) |");
    for (int i = 0; i < width; i++) printf("-");
    printf("| Magnitude\n");
    
    // Plot first half (positive frequencies)
    int step = (n/2) / 20;  // Show 20 frequency bins
    if (step < 1) step = 1;
    
    for (int i = 0; i < n/2; i += step) {
        double freq = i * sample_rate / n;
        int bar_length = (int)(width * magnitude[i] / max_mag);
        
        printf("%8.1f |", freq);
        for (int j = 0; j < bar_length; j++) printf("█");
        for (int j = bar_length; j < width; j++) printf(" ");
        printf("| %.3f\n", magnitude[i]);
    }
    
    free(magnitude);
}

void save_gnuplot_script(const char* script_name, const char* data_file, 
                        const char* output_image, const char* title) {
    FILE* fp = fopen(script_name, "w");
    if (!fp) return;
    
    fprintf(fp, "#!/usr/bin/gnuplot\n");
    fprintf(fp, "set terminal png size 800,600\n");
    fprintf(fp, "set output '%s'\n", output_image);
    fprintf(fp, "set title '%s'\n", title);
    fprintf(fp, "set xlabel 'Frequency (Hz)'\n");
    fprintf(fp, "set ylabel 'Magnitude'\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "plot '%s' using 1:4 with lines title 'Spectrum'\n", data_file);
    
    fclose(fp);
}

// Zero-padding utilities
complex_t* zero_pad(complex_t* signal, int original_size, int padded_size) {
    if (padded_size < original_size) return NULL;
    
    complex_t* padded = allocate_complex_array(padded_size);
    memcpy(padded, signal, original_size * sizeof(complex_t));
    memset(padded + original_size, 0, (padded_size - original_size) * sizeof(complex_t));
    
    return padded;
}

// Frequency shift (modulation)
void frequency_shift(complex_t* signal, int n, double shift_freq, double sample_rate) {
    for (int i = 0; i < n; i++) {
        double phase = TWO_PI * shift_freq * i / sample_rate;
        signal[i] *= cexp(I * phase);
    }
}

// Main demonstration
int main() {
    printf("FFT Utility Functions Demo\n");
    printf("==========================\n");
    
    int n = 256;
    double sample_rate = 1000.0;
    
    // Test signal generation
    printf("\n1. Signal Generation:\n");
    complex_t* signal = allocate_complex_array(n);
    
    // Generate chirp
    generate_chirp_signal(signal, n, 10, 100, sample_rate);
    printf("   - Generated chirp signal (10-100 Hz)\n");
    
    // Add noise
    complex_t* noise = allocate_complex_array(n);
    generate_noise(noise, n, 0.1, 42);
    
    for (int i = 0; i < n; i++) {
        signal[i] += noise[i];
    }
    
    double snr = compute_snr(signal, noise, n);
    printf("   - Added noise (SNR = %.1f dB)\n", snr);
    
    // Test windowing
    printf("\n2. Window Functions:\n");
    apply_window_tukey(signal, n, 0.5);
    printf("   - Applied Tukey window (α = 0.5)\n");
    
    // Test FFT and analysis
    printf("\n3. Spectrum Analysis:\n");
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Find peaks
    double peak_freqs[10], peak_mags[10];
    int num_peaks;
    find_spectral_peaks(signal, n, sample_rate, peak_freqs, peak_mags, &num_peaks, 10);
    
    printf("   - Found %d spectral peaks:\n", num_peaks);
    for (int i = 0; i < num_peaks && i < 5; i++) {
        printf("     Peak %d: %.1f Hz (mag = %.2f)\n", 
               i+1, peak_freqs[i], peak_mags[i]);
    }
    
    // Test file I/O
    printf("\n4. File I/O:\n");
    const char* filename = "fft_data.txt";
    save_complex_array(filename, signal, n);
    printf("   - Saved spectrum to %s\n", filename);
    
    // Test visualization
    printf("\n5. Visualization:\n");
    print_spectrum_plot(signal, n, sample_rate, 50);
    
    // Generate gnuplot script
    save_gnuplot_script("plot_spectrum.gnu", filename, "spectrum.png", "FFT Spectrum");
    printf("\n   - Generated gnuplot script 'plot_spectrum.gnu'\n");
    printf("   - Run 'gnuplot plot_spectrum.gnu' to create spectrum.png\n");
    
    // Test zero-padding
    printf("\n6. Zero-padding:\n");
    int padded_size = 512;
    complex_t* padded = zero_pad(signal, n, padded_size);
    printf("   - Zero-padded from %d to %d samples\n", n, padded_size);
    printf("   - Frequency resolution improved from %.2f to %.2f Hz\n",
           sample_rate/n, sample_rate/padded_size);
    
    // Clean up
    free_complex_array(signal);
    free_complex_array(noise);
    free_complex_array(padded);
    
    return 0;
}
