%% Algorithm for detection of signal components using the PFT
%% Illustration of misdetection when demodulation term does not match any of the components
%% Copyright (c) 2015 Irena Orovic

function figure9

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

% Thee demodualted signal versions where demodulation terms are incorrect:
% x1 for incorrect a=18, x2 for incorrect value a=-12, x3 for incorrect
% value -16; 

x1=s.*exp(j*t.^2*18*M);     
x2=s.*exp(-j*t.^2*12*M);     
x3=s.*exp(-j*t.^2*16*M);  

xx1=zeros(size(s));
xx2=zeros(size(s));
xx3=zeros(size(s));

% Defining random combination of N elements corresponding to the
% random positions of measurements
q=randperm(N);
% random positions of measurements: K out of N random elements
% corresponding to measurements positions; 

qM=q(1:K);

% Defining measurement vectors for three versions of demodulated signals
% x1, x2 and x3

xx1(qM)=x1(qM);
xx2(qM)=x2(qM);
xx3(qM)=x3(qM);


%% Algorithm for threshold calculation 
%(applying the threshold on measurements with incorrect demodulations
% returns no signal components - the threshold will be above all components); 

for  a=[18,-12,-16]; 
    %Signal demodulation 
xx=s.*exp(j*a*M*t.^2);   

%% Creating CS matrix starting from the DFT matrix
B=dftmtx(N);
Binv=inv(B);
xf=B*xx';

%creating measurement matrix
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

end

% Illustration (Figure 8 from the paper) - Applying the calculated
% thresholds on the DFT vectors in the case when demodulation terms are
% incorrect (three examples are are considered for a=18,a=-12,and a=-16); 

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
