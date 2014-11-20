#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../includes/fft842.c"

// Filter Length
#define Nf 256
// Length of Signal
#define N 25600
// Sampling frequency
const double Fs = 11025;

int main(int argc, char** argv)
{
  // Input stream for filter
  std::ifstream filterIn("../data/LowPassFilter.dat");

  // filter of length Nf = 256
  double h[Nf];

  // input vairable
  double in;

  // Read in the filter data
  for(int n = 0; n < Nf; ++n)
  {
    filterIn >> in;
    h[n] = in;
  }

  // Output streams for the input x signal
  // and the output y signal
  std::ofstream x_dat("../data/x.dat");
  std::ofstream y_dat("../data/y.dat");

  // input x signal of length N = 25600
  double x[N];

  // output y signal of Length N = 25600
  double y[N];

  // f0 = f/Fs
  // Normalized frequency
  double f = atof(argv[1]);
  double f0 = f/Fs;

  // Generate input signal x[n]
  for(int n = 0; n < N; ++n)
  {
    x[n] = cos(2*M_PI*f0*n);
    x_dat << x[n] << std::endl;
  }

  double temp;
  for(int n = 0; n < N; ++n)
  {
    temp = 0;
    for(int k = 0; k < Nf; ++k)
    {
      temp += x[n-k]*h[k];
    }
    y[n] = temp;
    y_dat << y[n] << std::endl;
  }

  return 0;
}
