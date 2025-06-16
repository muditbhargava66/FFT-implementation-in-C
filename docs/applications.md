# Applications Guide

This guide demonstrates practical applications of FFT in real-world scenarios.

## Audio Processing

### 1. Real-Time Spectrum Analyzer

Build a spectrum analyzer for audio signals:

```c
#include "fft_common.h"
#include "fft_algorithms.h"

typedef struct {
    int sample_rate;
    int fft_size;
    int hop_size;
    complex_t* window;
    complex_t* fft_buffer;
    double* magnitude_buffer;
} spectrum_analyzer_t;

spectrum_analyzer_t* create_spectrum_analyzer(int sample_rate, int fft_size) {
    spectrum_analyzer_t* analyzer = malloc(sizeof(spectrum_analyzer_t));
    analyzer->sample_rate = sample_rate;
    analyzer->fft_size = fft_size;
    analyzer->hop_size = fft_size / 2;  // 50% overlap
    
    // Allocate buffers
    analyzer->fft_buffer = allocate_complex_array(fft_size);
    analyzer->magnitude_buffer = malloc((fft_size/2 + 1) * sizeof(double));
    analyzer->window = allocate_complex_array(fft_size);
    
    // Precompute window
    for (int i = 0; i < fft_size; i++) {
        analyzer->window[i] = 0.5 - 0.5 * cos(2 * PI * i / (fft_size - 1));
    }
    
    return analyzer;
}

void analyze_audio_frame(spectrum_analyzer_t* analyzer, double* audio_frame) {
    // Copy and window the audio
    for (int i = 0; i < analyzer->fft_size; i++) {
        analyzer->fft_buffer[i] = audio_frame[i] * analyzer->window[i];
    }
    
    // Perform FFT
    radix2_dit_fft(analyzer->fft_buffer, analyzer->fft_size, FFT_FORWARD);
    
    // Compute magnitude spectrum
    for (int i = 0; i <= analyzer->fft_size/2; i++) {
        double mag = cabs(analyzer->fft_buffer[i]);
        analyzer->magnitude_buffer[i] = 20 * log10(mag + 1e-10);  // dB scale
    }
}

// Get frequency for a given bin
double get_bin_frequency(spectrum_analyzer_t* analyzer, int bin) {
    return bin * analyzer->sample_rate / (double)analyzer->fft_size;
}
```

### 2. Pitch Detection

Detect musical notes in audio:

```c
typedef struct {
    double frequency;
    const char* note_name;
    int octave;
} musical_note_t;

musical_note_t detect_pitch(complex_t* audio, int n, double sample_rate) {
    musical_note_t result = {0};
    
    // Apply window
    apply_window_hann(audio, n);
    
    // FFT
    radix2_dit_fft(audio, n, FFT_FORWARD);
    
    // Find peak in spectrum
    double max_magnitude = 0;
    int peak_bin = 0;
    
    // Search in musical frequency range (80 Hz to 2000 Hz)
    int start_bin = (int)(80.0 * n / sample_rate);
    int end_bin = (int)(2000.0 * n / sample_rate);
    
    for (int i = start_bin; i < end_bin && i < n/2; i++) {
        double mag = cabs(audio[i]);
        if (mag > max_magnitude) {
            max_magnitude = mag;
            peak_bin = i;
        }
    }
    
    // Refine frequency estimate using quadratic interpolation
    if (peak_bin > 0 && peak_bin < n/2 - 1) {
        double y1 = cabs(audio[peak_bin - 1]);
        double y2 = cabs(audio[peak_bin]);
        double y3 = cabs(audio[peak_bin + 1]);
        
        double delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
        result.frequency = (peak_bin + delta) * sample_rate / n;
    } else {
        result.frequency = peak_bin * sample_rate / (double)n;
    }
    
    // Convert frequency to musical note
    result = frequency_to_note(result.frequency);
    
    return result;
}

musical_note_t frequency_to_note(double freq) {
    static const char* notes[] = {"C", "C#", "D", "D#", "E", "F", 
                                  "F#", "G", "G#", "A", "A#", "B"};
    musical_note_t note;
    
    // A4 = 440 Hz reference
    double a4 = 440.0;
    double c0 = a4 * pow(2, -4.75);
    
    double half_steps = 12 * log2(freq / c0);
    int note_index = ((int)round(half_steps)) % 12;
    note.octave = ((int)round(half_steps)) / 12;
    note.note_name = notes[note_index];
    note.frequency = freq;
    
    return note;
}
```

## Digital Filtering

### 1. Frequency Domain Filtering

Implement various filters using FFT:

```c
// Filter types
typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP,
    FILTER_CUSTOM
} filter_type_t;

// Apply frequency domain filter
void fft_filter(complex_t* signal, int n, double sample_rate,
                filter_type_t type, double freq1, double freq2) {
    // Transform to frequency domain
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Apply filter
    for (int k = 0; k < n; k++) {
        double freq = k * sample_rate / n;
        if (k > n/2) freq = (n - k) * sample_rate / n;  // Negative frequencies
        
        double gain = 0.0;
        
        switch (type) {
            case FILTER_LOWPASS:
                gain = (freq <= freq1) ? 1.0 : 0.0;
                break;
                
            case FILTER_HIGHPASS:
                gain = (freq >= freq1) ? 1.0 : 0.0;
                break;
                
            case FILTER_BANDPASS:
                gain = (freq >= freq1 && freq <= freq2) ? 1.0 : 0.0;
                break;
                
            case FILTER_BANDSTOP:
                gain = (freq < freq1 || freq > freq2) ? 1.0 : 0.0;
                break;
                
            default:
                gain = 1.0;
        }
        
        signal[k] *= gain;
    }
    
    // Transform back to time domain
    radix2_dit_fft(signal, n, FFT_INVERSE);
}

// Smooth filter transitions to reduce ringing
void apply_filter_transition(complex_t* spectrum, int n, double sample_rate,
                           double cutoff_freq, double transition_width) {
    for (int k = 0; k < n; k++) {
        double freq = k * sample_rate / n;
        if (k > n/2) freq = (n - k) * sample_rate / n;
        
        double gain;
        if (freq < cutoff_freq - transition_width/2) {
            gain = 1.0;
        } else if (freq > cutoff_freq + transition_width/2) {
            gain = 0.0;
        } else {
            // Smooth transition using raised cosine
            double t = (freq - cutoff_freq + transition_width/2) / transition_width;
            gain = 0.5 * (1 + cos(PI * t));
        }
        
        spectrum[k] *= gain;
    }
}
```

### 2. Adaptive Noise Reduction

Remove noise from signals:

```c
typedef struct {
    int fft_size;
    complex_t* noise_profile;
    double noise_factor;
} noise_reducer_t;

// Learn noise profile from silent section
void learn_noise_profile(noise_reducer_t* nr, complex_t* noise_sample, int n) {
    complex_t* temp = allocate_complex_array(n);
    memcpy(temp, noise_sample, n * sizeof(complex_t));
    
    apply_window_hann(temp, n);
    radix2_dit_fft(temp, n, FFT_FORWARD);
    
    // Store magnitude spectrum as noise profile
    for (int i = 0; i < n; i++) {
        nr->noise_profile[i] = cabs(temp[i]) * nr->noise_factor;
    }
    
    free_complex_array(temp);
}

// Apply spectral subtraction
void reduce_noise(noise_reducer_t* nr, complex_t* signal, int n) {
    // Transform to frequency domain
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Spectral subtraction
    for (int k = 0; k < n; k++) {
        double mag = cabs(signal[k]);
        double phase = carg(signal[k]);
        double noise_level = creal(nr->noise_profile[k]);
        
        // Subtract noise magnitude
        mag = fmax(mag - noise_level, mag * 0.1);  // Keep 10% minimum
        
        // Reconstruct complex number
        signal[k] = mag * cexp(I * phase);
    }
    
    // Transform back
    radix2_dit_fft(signal, n, FFT_INVERSE);
}
```

## Signal Analysis

### 1. Cross-Correlation for Time Delay

Find time delay between signals:

```c
double find_time_delay(complex_t* signal1, complex_t* signal2, int n, double sample_rate) {
    int fft_size = next_power_of_two(2 * n);
    complex_t* x1 = allocate_complex_array(fft_size);
    complex_t* x2 = allocate_complex_array(fft_size);
    
    // Zero-pad signals
    memset(x1, 0, fft_size * sizeof(complex_t));
    memset(x2, 0, fft_size * sizeof(complex_t));
    memcpy(x1, signal1, n * sizeof(complex_t));
    memcpy(x2, signal2, n * sizeof(complex_t));
    
    // FFT both signals
    radix2_dit_fft(x1, fft_size, FFT_FORWARD);
    radix2_dit_fft(x2, fft_size, FFT_FORWARD);
    
    // Cross-correlation in frequency domain
    for (int k = 0; k < fft_size; k++) {
        x1[k] = conj(x1[k]) * x2[k];
    }
    
    // Inverse FFT
    radix2_dit_fft(x1, fft_size, FFT_INVERSE);
    
    // Find peak in correlation
    double max_corr = 0;
    int max_index = 0;
    
    for (int i = 0; i < fft_size; i++) {
        double corr = cabs(x1[i]);
        if (corr > max_corr) {
            max_corr = corr;
            max_index = i;
        }
    }
    
    // Convert index to time delay
    double delay;
    if (max_index <= fft_size / 2) {
        delay = max_index / sample_rate;
    } else {
        delay = (max_index - fft_size) / sample_rate;
    }
    
    free_complex_array(x1);
    free_complex_array(x2);
    
    return delay;
}
```

### 2. Spectral Analysis

Compute various spectral features:

```c
typedef struct {
    double centroid;        // Spectral center of mass
    double spread;          // Spectral spread
    double skewness;        // Asymmetry
    double kurtosis;        // Peakedness
    double flux;            // Spectral change
    double rolloff;         // 95% energy frequency
    double brightness;      // High frequency content
    double flatness;        // Tonality measure
} spectral_features_t;

spectral_features_t compute_spectral_features(complex_t* spectrum, int n, double sample_rate) {
    spectral_features_t features = {0};
    double* magnitude = compute_magnitude(spectrum, n);
    
    // Compute spectral centroid and moments
    double sum_mag = 0, sum_freq_mag = 0;
    double sum_freq2_mag = 0, sum_freq3_mag = 0, sum_freq4_mag = 0;
    
    for (int k = 1; k < n/2; k++) {  // Skip DC
        double freq = k * sample_rate / n;
        double mag = magnitude[k];
        
        sum_mag += mag;
        sum_freq_mag += freq * mag;
        sum_freq2_mag += freq * freq * mag;
        sum_freq3_mag += freq * freq * freq * mag;
        sum_freq4_mag += freq * freq * freq * freq * mag;
    }
    
    // Centroid
    if (sum_mag > 0) {
        features.centroid = sum_freq_mag / sum_mag;
        
        // Spread (standard deviation)
        double variance = sum_freq2_mag / sum_mag - features.centroid * features.centroid;
        features.spread = sqrt(variance);
        
        // Skewness
        if (features.spread > 0) {
            double mean = features.centroid;
            double m3 = sum_freq3_mag / sum_mag - 3 * mean * variance - mean * mean * mean;
            features.skewness = m3 / pow(features.spread, 3);
            
            // Kurtosis
            double m4 = sum_freq4_mag / sum_mag - 4 * mean * m3 - 6 * mean * mean * variance - mean * mean * mean * mean;
            features.kurtosis = m4 / pow(variance, 2) - 3;
        }
    }
    
    // Spectral rolloff (95% of energy)
    double cumulative_energy = 0;
    double total_energy = sum_mag;
    for (int k = 0; k < n/2; k++) {
        cumulative_energy += magnitude[k];
        if (cumulative_energy >= 0.95 * total_energy) {
            features.rolloff = k * sample_rate / n;
            break;
        }
    }
    
    // Brightness (ratio of high frequency content)
    double high_freq_energy = 0;
    int bright_threshold = (int)(1500.0 * n / sample_rate);  // Above 1500 Hz
    for (int k = bright_threshold; k < n/2; k++) {
        high_freq_energy += magnitude[k];
    }
    features.brightness = high_freq_energy / total_energy;
    
    // Spectral flatness (geometric mean / arithmetic mean)
    double log_sum = 0;
    int count = 0;
    for (int k = 1; k < n/2; k++) {
        if (magnitude[k] > 0) {
            log_sum += log(magnitude[k]);
            count++;
        }
    }
    if (count > 0 && sum_mag > 0) {
        double geometric_mean = exp(log_sum / count);
        double arithmetic_mean = sum_mag / count;
        features.flatness = geometric_mean / arithmetic_mean;
    }
    
    free(magnitude);
    return features;
}
```

## Image Processing

### 1. 2D FFT for Images

Process images in frequency domain:

```c
// 2D FFT using row-column decomposition
void fft_2d(complex_t** image, int rows, int cols, fft_direction dir) {
    // FFT each row
    for (int i = 0; i < rows; i++) {
        radix2_dit_fft(image[i], cols, dir);
    }
    
    // FFT each column
    complex_t* column = allocate_complex_array(rows);
    for (int j = 0; j < cols; j++) {
        // Extract column
        for (int i = 0; i < rows; i++) {
            column[i] = image[i][j];
        }
        
        // FFT column
        radix2_dit_fft(column, rows, dir);
        
        // Put back
        for (int i = 0; i < rows; i++) {
            image[i][j] = column[i];
        }
    }
    free_complex_array(column);
}

// Apply frequency domain filter to image
void filter_image_fft(complex_t** image, int rows, int cols, 
                     double cutoff_frequency) {
    // Forward 2D FFT
    fft_2d(image, rows, cols, FFT_FORWARD);
    
    // Apply low-pass filter
    double center_row = rows / 2.0;
    double center_col = cols / 2.0;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double dist = sqrt(pow(i - center_row, 2) + pow(j - center_col, 2));
            if (dist > cutoff_frequency) {
                image[i][j] = 0;
            }
        }
    }
    
    // Inverse 2D FFT
    fft_2d(image, rows, cols, FFT_INVERSE);
}
```

## Communications

### 1. OFDM Modulation

Implement basic OFDM for communications:

```c
typedef struct {
    int num_subcarriers;
    int fft_size;
    int cyclic_prefix_len;
    complex_t* subcarrier_data;
    complex_t* time_domain_symbol;
} ofdm_symbol_t;

// Modulate data using OFDM
void ofdm_modulate(ofdm_symbol_t* ofdm, complex_t* data_symbols) {
    // Clear FFT buffer
    memset(ofdm->time_domain_symbol, 0, ofdm->fft_size * sizeof(complex_t));
    
    // Map data to subcarriers (with guard bands)
    int guard_band = (ofdm->fft_size - ofdm->num_subcarriers) / 2;
    for (int i = 0; i < ofdm->num_subcarriers; i++) {
        ofdm->time_domain_symbol[guard_band + i] = data_symbols[i];
    }
    
    // IFFT to convert to time domain
    radix2_dit_fft(ofdm->time_domain_symbol, ofdm->fft_size, FFT_INVERSE);
    
    // Add cyclic prefix
    complex_t* output = allocate_complex_array(ofdm->fft_size + ofdm->cyclic_prefix_len);
    
    // Copy end of symbol to beginning (cyclic prefix)
    memcpy(output, 
           &ofdm->time_domain_symbol[ofdm->fft_size - ofdm->cyclic_prefix_len],
           ofdm->cyclic_prefix_len * sizeof(complex_t));
    
    // Copy full symbol
    memcpy(&output[ofdm->cyclic_prefix_len],
           ofdm->time_domain_symbol,
           ofdm->fft_size * sizeof(complex_t));
    
    // Replace with full symbol including CP
    free(ofdm->time_domain_symbol);
    ofdm->time_domain_symbol = output;
}

// Demodulate OFDM symbol
void ofdm_demodulate(ofdm_symbol_t* ofdm, complex_t* received_symbol) {
    // Remove cyclic prefix
    complex_t* symbol_no_cp = &received_symbol[ofdm->cyclic_prefix_len];
    
    // FFT to convert to frequency domain
    memcpy(ofdm->time_domain_symbol, symbol_no_cp, ofdm->fft_size * sizeof(complex_t));
    radix2_dit_fft(ofdm->time_domain_symbol, ofdm->fft_size, FFT_FORWARD);
    
    // Extract data from subcarriers
    int guard_band = (ofdm->fft_size - ofdm->num_subcarriers) / 2;
    for (int i = 0; i < ofdm->num_subcarriers; i++) {
        ofdm->subcarrier_data[i] = ofdm->time_domain_symbol[guard_band + i];
    }
}
```

## Best Practices

### 1. Window Selection

Choose appropriate windows for different applications:

| Application | Recommended Window | Reason |
|------------|-------------------|---------|
| General spectrum analysis | Hann | Good frequency resolution |
| Accurate amplitude | Flat-top | Minimal amplitude error |
| Transient detection | Rectangular | Best time resolution |
| Low sidelobe | Blackman-Harris | -92 dB sidelobes |

### 2. Zero Padding

Use zero padding for:
- Smoother spectrum visualization
- Frequency interpolation
- Preventing circular convolution

```c
// Zero pad for smoother spectrum
complex_t* zero_pad_signal(complex_t* signal, int original_size, int padded_size) {
    complex_t* padded = allocate_complex_array(padded_size);
    memset(padded, 0, padded_size * sizeof(complex_t));
    memcpy(padded, signal, original_size * sizeof(complex_t));
    return padded;
}
```

### 3. Overlap Processing

For continuous signal processing:

```c
typedef struct {
    int window_size;
    int hop_size;
    complex_t* overlap_buffer;
} overlap_processor_t;

void process_with_overlap(overlap_processor_t* proc, 
                         complex_t* input, 
                         int input_size,
                         void (*process_func)(complex_t*, int)) {
    int pos = 0;
    while (pos + proc->window_size <= input_size) {
        // Copy window of data
        complex_t* window = allocate_complex_array(proc->window_size);
        memcpy(window, &input[pos], proc->window_size * sizeof(complex_t));
        
        // Process window
        process_func(window, proc->window_size);
        
        // Handle overlap...
        
        free_complex_array(window);
        pos += proc->hop_size;
    }
}
```

## Conclusion

FFT enables efficient implementation of many signal processing applications. The key is choosing the right algorithm, parameters, and techniques for your specific use case. Always consider:

1. **Accuracy requirements** - floating vs fixed point
2. **Real-time constraints** - latency and throughput
3. **Memory limitations** - in-place vs out-of-place
4. **Power consumption** - especially for embedded systems

Experiment with the provided examples and adapt them to your needs!
