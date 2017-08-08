%% Algorithm for detection of signal components using the PFT
%% Copyright (c) 2015 Irena Orovic

function figure8

close all;
clear all;
clc;
%Length of the signal
N=1024;
% An auxilary parameter used in signal definition
M=32;

% Number of available samples / measurements
K=128;

% time axes
t=-N/2:1:N/2-1/N;
t=t/N;
T=[]; 

%% Signal consist of three chirps components
s=1*exp(-j*t.^2*8*M+j*t*M*16*pi)+1*exp(j*t.^2*32*M-j*t*M*16*pi)+1*exp(j*t.^2*8*M-j*t*M*16*pi);     

% The user needs to change demodulation parameters in the range a=-64:1:64
% in order to detect values corresponding to each signal component. 
% For simplicity we assume that we have already found values of a ; 
% As an example, try to use wrong a, e.g. a=16 and the result of plot in line (57)
% will be zero line; on the ither side for a=8, the result will be the peak at the position of one of the
% signal components;


%% Algorithm for threshold calculation 
%(applying the threshold on measurements with correct demodulations
%returns corresponding signal component); 

for  a=[8,-32,-8] 
    
  % signal demodulation - when we match the value of a corresponding to one
  % of the signal components, then the component is demodulated ; 
  
xx=s.*exp(j*a*M*t.^2);   

%% Creating CS matrix starting from the DFT matrix

B=dftmtx(N);
Binv=inv(B);
xf=B*xx';
q=randperm(N);
% random positions of measurements
qM=q(1:K);
% CS matrix corresponding to measurements positions
Ax=Binv(q(1:K),:);

% Defining the vector of measurements
y=(Ax*xf);
y=y';

%% Signal reconstruction approach based on the Single Iteration threhold based Algorithm (SIRA)
%% We are intersted in caluclated threshold value
[xp,X,Thresh]=sira(y,qM,1024,128);
T=[T Thresh];
% Reconstructed signal
sig_rec=ifft((xp)).*exp(-j*2*pi*a*M*t.^2);

%% figure,plot(abs(xp))

end

xx1=zeros(size(s));
xx2=zeros(size(s));
xx3=zeros(size(s));

%% For the illustration: x1,x2 and x3 are demodulated signal components; xx1, xx2 and xx3 are demodulated signal components with missing samples
x1=s.*exp(j*t.^2*8*M);     
x2=s.*exp(-j*t.^2*32*M);     
x3=s.*exp(-j*t.^2*8*M); 
xx1(qM)=x1(qM);
xx2(qM)=x2(qM);
xx3(qM)=x3(qM);

% Illustration (Figure 8 from the paper) - Applying the calculated
% thresholds on the DFT vectors calculated from the demodulated signal
% components with missing samples ;
figure, 
subplot(311),plot(abs(fft(xx1)),'k'), xlabel({'Frequency samples', '(a)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(1)*ones(1,1024),'r'),axis tight;
subplot(312),plot(abs(fft(xx2)),'k'), xlabel({'Frequency samples', '(b)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(2)*ones(1,1024),'r'),axis tight;
subplot(313),plot(abs(fft(xx3)),'k'), xlabel({'Frequency samples', '(c)'}), ylabel('Amplitude'), colormap(1-gray) 
hold on, plot(K*T(3)*ones(1,1024),'r'),axis tight;

end

%% Signal reconstruction based on the appropriate threshold calculation
% More details about this reconstruction algorithm can be found in "S.
% Stankovi?, I. Orovi?, and LJ. Stankovi?, “An Automated Signal
% Reconstruction Method based on Analysis of Compressive Sensed Signals in
% Noisy Environment,” Signal Processing, vol. 104, Nov 2014, pp. 43 - 50, 2014"
%
% INPUT:
%       y - measurements vector
%       qM - Positions of measurements/ available samples
%       N - total number of signal samples
%       M - number of measurements / available samples
%
% OUTPUT:
%       XX1 - Reconstructed DFT components
%       X - Initial DFT vector
%       Thresh - Threshold value
%
%% Copyright (c) 2014 Irena Orovic
function [XX1,X,Thresh]=sira(y,qM,N,M);

xr=zeros(1,N);

Amp=sum(abs(y).^2)/M;
% Variance of the noise caused by the missing samples (that are zero-valued) 
sigma_ms2=M*(N-M)/(N-1)*Amp;

% Approximative expression for probability that all noise components will
% be below certain threshold ; 

Q=1-0.99^(1/(N));

% Expression for threshold calculation ; 

Thresh=1/M*sqrt(-(sigma_ms2)*log(Q));
X=zeros(1,N);

% Calculating the initial DFT of the signal

for k=0:N-1
    X(k+1)=1/M*sum(y.*exp(-j*2*pi*k*qM/N));
end

A=Amp;    

[a,b]=find(abs(X)>Thresh);

AFULL=dftmtx(N);
Acs=AFULL(qM,b);
X1=pinv(Acs'*Acs)*(Acs'*y');

% Reconstructed and rescalled DFT
XX1=zeros(1,N);
XX1(b)=N*X1;

end