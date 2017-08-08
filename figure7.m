%% The effects of demodulation on the discrete Fourier transform spectrum. 
%% The polynomial Fourier transform with an appropriate demodulation term 
%% can be considered as compressible for the dominant component matched
%% by the demodulation term


function figure7


close all;
clear all;
clc;
%Length of the signal
N=1024;

% An auxiliary parameter to avoid large numbers
M=32;

% time axes
t=-N/2:1:N/2-1/N;
t=t/N;

%% Signal consist of three chirps components
s=1*exp(-j*t.^2*8*M+j*t*M*16*pi)+1*exp(j*t.^2*32*M-j*t*M*16*pi)+1*exp(j*t.^2*8*M-j*t*M*16*pi);     

% Signal components demodulation for each signal component  
x1=s.*exp(j*t.^2*8*M);     
x2=s.*exp(-j*t.^2*32*M);     
x3=s.*exp(-j*t.^2*8*M);     

%% Figure 7 in the paper
figure(1), 
subplot(411),plot(abs(fft(s)),'k'),axis tight; xlabel({'Frequency samples', '(a)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(412),plot(abs(fft(x1)),'k'),axis tight; xlabel({'Frequency samples', '(b)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(413),plot(abs(fft(x2)),'k'),axis tight; xlabel({'Frequency samples', '(c)'}), ylabel('Amplitude'), colormap(1-gray) 
subplot(414),plot(abs(fft(x3)),'k'),axis tight; xlabel({'Frequency samples', '(d)'}), ylabel('Amplitude'), colormap(1-gray) 

end


