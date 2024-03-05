#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define N 8  // Size of FFT
#define PI 3.1415926535897932384626433832795  // Value of Pi 

typedef double complex cplx;  // Complex number type

// Perform in-place FFT on input array
void fft(cplx x[N]){

  cplx even[N/2], odd[N/2], twiddle;
  int i, j;
  
  // Bit reversal permutation
  for(i=0; i<N; i++) {
    int reversed = 0;
    for(j=0; j<3; j++)   // log2(N) = 3 bits needed to represent indices
      reversed = (reversed << 1) | (i >> j & 1);  
    
    if(reversed > i){
      // Swap elements at i and reversed
      cplx temp = x[i]; 
      x[i] = x[reversed];
      x[reversed] = temp;
    }
  }
  
  // Danielson-Lanczos Algorithm 
  for(int halfSize=1; halfSize<N; halfSize*=2) {

    // Combine even/odd transforms of size halfSize to get transform of size 2*halfSize 
    for(i=0; i<N; i+=2*halfSize) {
      for(j=0; j<halfSize; j++) {
        
        // Even indices 
        even[j] = x[i + j];
        // Odd indices
        odd[j]  = x[i + j + halfSize];
        
        // Twiddle factor 
        double angle = -2*PI*j/(2*halfSize);
        twiddle = cexp(I*angle);
        
        // Combine results
        x[i + j]         = even[j] + twiddle*odd[j];
        x[i + j +halfSize] = even[j] - twiddle*odd[j];
      }
    }
  }
}

// Print the FFT output
void showOutput(cplx x[N]){

  printf("FFT: "); 
  for(int i=0; i<N; i++){
    double real = creal(x[i]);
    double imag = cimag(x[i]);
    if(fabs(real) < 1e-10) real = 0; 
    if(fabs(imag) < 1e-10) imag = 0;
    printf("(%.3lf, %.3lfi) ", real, imag); 
  }
  printf("\n");
}

int main(){

  clock_t start, end;   
  double timeTaken;

  // Input sequence 
  cplx input[N] = {1, 1, 1, 1, -1, -1, -1, -1};

  printf("Input: ");
  for(int i=0; i<N; i++) 
    printf("%2.0lf ", creal(input[i]));
  printf("\n");

  // Perform FFT
  start = clock();
  fft(input);
  end = clock();

  timeTaken = ((double)(end - start)) / CLOCKS_PER_SEC;

  showOutput(input);  
  printf("Time: %.6lf seconds\n", timeTaken);

  return 0;
}