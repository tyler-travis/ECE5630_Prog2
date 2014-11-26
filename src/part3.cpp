#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "../includes/fft842.c"

// Filter Length
#define Nf 256
// Length of Signal
#define N 496125
// Sampling frequency
const double Fs = 11025;

complx mult(complx, complx);

int main(int argc, char** argv)
{
  // Input stream for filter
  std::ifstream filterIn("../data/LowPassFilter.dat");
  // Input stream for signal x[n]
  std::ifstream xIn("../data/flute.dat");

  // filter of length Nf = 256
  // Nf*4 for zero padding
  complx h[4*Nf];

  // input vairable
  double in;

  // Read in the filter data
  for(int n = 0; n < Nf; ++n)
  {
    filterIn >> in;
    h[n].re = in;
    h[n].im = 0;
    h[n+Nf].re = 0;
    h[n+Nf].im = 0;
    h[n+2*Nf].re = 0;
    h[n+2*Nf].im = 0;
    h[n+3*Nf].re = 0;
    h[n+3*Nf].im = 0;
  }

  // Output streams for the input x signal
  // and the output y signal
  std::ofstream x_dat("../data/x.dat");
  std::ofstream y_dat("../data/y.dat");
  std::ofstream H_dat("../data/H.dat");

  if(!strcmp(argv[1],"fft"))
  {
    // input x signal of length N = 25600
    complx x[N];

    // output y signal of Length N = 25600
    complx y[N+Nf-1];

    // Generate input signal x[n]
    for(int n = 0; n < N; ++n)
    {
      xIn >> in;
      x[n].re = in;
      x_dat << x[n].re << std::endl;
    }

    //const int nfft = 1024;
    
    int M = 256;
    int overlap = M-1;
    int nfft = 1024; 
    int stepsize = nfft - overlap;

    complx H[nfft];
    memcpy(H, h, sizeof(h));
    // generate fft
    fft842(0, nfft, H);
    for(int i = 0; i < nfft; ++i)
    {
      H_dat << H[i].re << "\t" << H[i].im <<std::endl;
    }

    complx yt[nfft];
    complx xt[nfft];

    int position = 0;
    while(position + nfft <= N)
    {
      for(int j = 0; j < nfft; ++j)
      {
        xt[j] = x[j + position];
      }

      fft842(0, nfft, xt);
      for(int k = 0; k < nfft; ++k)
      {
        yt[k] = mult(xt[k], H[k]);
      }
      fft842(1, nfft, yt);
      for(int j = M-1; j < nfft; ++j)
      {
        y[j-M+position] = yt[j];
      }
      position += stepsize;
    }
    for(int n = 0; n < N; ++n)
    {
      y_dat << y[n].re << std::endl;
    }
  }
  else
  {
    // input x signal of length N = 25600
    double x[N];

    // output y signal of Length N = 25600
    double y[N];

    // Generate input signal x[n]
    for(int n = 0; n < N; ++n)
    {
      xIn >> in;
      x[n] = in;
      x_dat << x[n] << std::endl;
    }

    double temp;
    for(int n = 0; n < N; ++n)
    {
      temp = 0;
      for(int k = 0; k < Nf; ++k)
      {
        temp += x[n-k]*h[k].re;
      }
      y[n] = temp;
      y_dat << y[n] << std::endl;
    }
  }
  return 0;
}

complx mult(complx a, complx b)
{
  complx ret;
  ret.re = a.re * b.re - a.im * b.im;
  ret.im = a.re * b.im + a.im * b.re;
  return ret;
}
