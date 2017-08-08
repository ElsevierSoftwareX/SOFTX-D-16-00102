
% Estimated numbers of points to accurately estimate instantaneous frequency with (a) the
% polynomial distribution; (b) the complex-time distribution.

% The test signals are specificed using sig.m function:
% a) for ind=1 - the cubic phase signal for the illustration of PD
% b) for ind=2 - fast varying sinusoidal phase modulated signal for the application fo CTD

function figure11

clear all
close all
L=[];
rand('state',203)

%% For PD

% Changing the number of available samples i.e. measurements in the signal
% from 5 to 125
for M=5:5:125
% Length of the signal
N=128;

%vector specifying lag values
tau=(-1:2/N:1-2/N); 


Nk=(N)/2;
n=-Nk:Nk-1;
b=0;

ind=1;  % Set ind=1 for Polynomial distribution (PD)

%% Loading reference point fro IF calculation
load C3.mat 
REF=CSTF3;
[p1,p2]=find(abs(REF)>1);
for i=1:65
REF(p1(i),p2(i))=10;
end
%%--------------------------------------------
t1=-.5; t2=.5; 


for t=t1:2/N:t2
        b=b+1;
q=randperm(N);q1=q(1:N-M);
q=randperm(N);q2=q(1:N-M);
q=randperm(N);q3=q(1:N-M);
q=randperm(N);q4=q(1:N-M);

% Defining four terms in the Local Autocorrelation Function (LAF) for the PD
% calculation (the terms are with missing samples)
x1=sig(t+.675*tau,ind).^2;
x2=conj(sig(t-.675*tau,ind).^2);
x3=conj(sig(t+.85*tau,ind));
x4=sig(t-.85*tau,ind);

% Setting all but the M available samples to zero values (N-M missing samples in the LAF will be zero-valued);    
x1(q1)=0;x2(q2)=0;x3(q1)=0;x4(q2)=0;

 % LAF as a product of four terms in the case of PD
rxx=x1.*x2.*x3.*x4;

%creating dft matrix
B=dftmtx(N);

%Selecting random rows of the DFT matrix
qy=find(rxx~=0);

%creating measurement matrix
A=B(qy,:);
y1=(rxx(qy));

%Running one iteration kk=1 of the OMP algorithm

res = y1';
omega = [];
kk=1;
for iter=1:kk,
    
    yp = A'*res; 
    yp(omega) = 0;
    [jj,ii] = max(abs(yp));
    omega = union(omega,ii);
    a1 = A(:,omega);
    xp = pinv(a1)*y1';     
    res = y1' - a1*xp; 
    
end

xp1 = zeros(N,1);
xp1(omega) = xp;

%recovered signal in time domain
f1=(B*xp1);

X(2,1:length(f1))=fftshift(fft(flipud(f1)));
 
CSTF(b,:)=X(2,:)';  %% Compressive sensing base representation

end

l=0;
for i=1:65
if abs(CSTF(p1(i),p2(i)))>1
    l=l+1;
end
end
L=[L l];

end
figure,
SetFigureDefaults(18,4.5)
subplot(121),plot(L,'-*'),
grid on
xlabel('Number of available samples (out of 129)')
title('Number of correctly recovered points in the TF domain')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

%% For CTD


L=[];
rand('state',203)
% Changing the number of available samples i.e. measurements in the signal
% from 5 to 125
for M=5:5:125
    
% Length of the signal
N=128;

%vector specifying lag values
tau=(-1:2/N:1-2/N); 


Nk=(N)/2; n=-Nk:Nk-1;
b=0;


ind=2;  % Set ind=2 for Complex-time distribution (CTD)

%% Loading reference point fro IF calculation
load C2.mat 
REF=CSTF2;
[p1,p2]=find(abs(REF)>1);
for i=1:129
REF(p1(i),p2(i))=10;
end

t1=-1; t2=1;


for t=t1:2/N:t2
        b=b+1;
q=randperm(N);q1=q(1:N-M);
q=randperm(N);q2=q(1:N-M);
q=randperm(N);q3=q(1:N-M);
q=randperm(N);q4=q(1:N-M);

% Defining four terms in the Local Autocorrelation Function (LAF) for the PD
% calculation (the terms are with missing samples)
x1=sig(t+tau/4,ind);
x2=sig(t-tau/4,ind).^(-1);
x3=(sig(t+j*tau/4,ind).^(-j));
x4=sig(t-j*tau/4,ind).^j;

% Setting all but the M available samples to zero values (N-M missing samples in the LAF will be zero-valued);       
x1(q1)=0;x2(q2)=0;x3(q1)=0;x4(q2)=0;

% LAF as a product of four terms in the case of CTD
rxx=x1.*x2.*x3.*x4;
   
%creating dft matrix
B=dftmtx(N);

%Selecting random rows of the DFT matrix
qy=find(rxx~=0);

%creating measurement matrix
A=B(qy,:);
y1=(rxx(qy));

%Running one iteration kk=1 of the OMP algorithm

res = y1';
omega = [];
kk=1;
for iter=1:kk,
    
    yp = A'*res; 
    yp(omega) = 0;
    [jj,ii] = max(abs(yp));
    omega = union(omega,ii);
    a1 = A(:,omega);
    xp = pinv(a1)*y1';     
    res = y1' - a1*xp; 
    
end

xp1 = zeros(N,1);
xp1(omega) = xp;

%recovered signal in time domain
f1=(B*xp1);

X(2,1:length(f1))=fftshift(fft(flipud(f1)));
 
CSTF(b,:)=X(2,:)';  %% Compressive sensing base representation

end

l=0;
for i=1:129
if abs(CSTF(p1(i),p2(i)))>1
    l=l+1;
end
end

L=[L l];

end


subplot(122),plot(L,'-*'),
grid on
xlabel('Number of available samples (out of 129)')
title('Number of correctly recovered points in the TF domain')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')
end


%% Signal definition: ind=1 cubic phase signal for PD calcualtion; ind=2 - sinusoidal phase modulated sigbal for CTD calculation

function x=sig(t,ind)
  
if (ind==1)
x=exp(j*160*pi*t.^3-0*j*56*pi*t.^2-j*192*pi*t);  %% for PD

else

x=exp(j*6*sin(2.4*pi*t)+j*3*cos(1.5*pi*t)-1*j*20*pi*t.^2); %% for CTD

end

end
