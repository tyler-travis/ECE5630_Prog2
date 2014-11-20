x = importdata('data/x.dat');
y = importdata('data/y.dat');

X = fftshift(fft(x));
Y = fftshift(fft(y));

N = length(x);
Fs = 11025;

nx = ((-N/2):(N-N/2-1))/N;
ny = nx;

figure(1); subplot(2,1,1); plot(x); xlim([0 2000])
figure(1); subplot(2,1,2); plot(nx,20*log10(abs(X)));

figure(2); subplot(2,1,1); plot(y); xlim([0 2000])
figure(2); subplot(2,1,2); plot(ny,20*log10(abs(Y))); 

Htemp = dlmread('data/H.dat');
Hf = complex(Htemp(:,1), Htemp(:,2));
figure(3); plot(20*log10(abs(Hf)));