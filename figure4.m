%% Algorithm for separation of stationary and nonstationary signal compoennts using the concept of Compressive sensing in STFT domain

function figure4

close all;
clear all;

%Total length of the signal
N=4096;

m=N;
M=64;
X_INIT=[];
XREZ=[];

%%-------------------------------------------------------------------------
%% Defining signal with 7 stationary components, and several nonstationary
%% components corresponding to micro-Doppler and noise; 

t=0:1/N:M/N-1/N;

mi=0;
for m=0:M:N-1
  mi=mi+1;
    tp=m/N+t;
   % Stationary signal components calculated using non-overlappign windows
    x=exp(j*tp*M*12*pi)+1.7*exp(j*tp*M*16*pi+j*pi/4)+3.5*exp(j*tp*M*18*pi+j*pi/4)+3.5*exp(j*tp*M*26*pi-j*pi/8)...
    +2*exp(-j*tp*M*20*pi-j*pi/3)+2*exp(-j*tp*M*24*pi-j*pi/3)+1.5*exp(-j*tp*M*8*pi+j*pi/8);
% Entire signal with stationary components   
XREZ=[XREZ, x]; 
   
% Adding nonstationary components o the signal of interest

     x=x+4.5*exp(j*2*pi*tp*M*round(16*sin(m/N*0.7*pi)));
    
    NE=10000;
    Dm=1/64;
    t1=Dm*6-.5/4096; x=x+5*exp(j*tp*M*(8)*pi+.3).*exp(-((tp-t1)/(3*Dm)).^NE);
    x=x+5*exp(-j*tp*M*(48)*pi+.3).*exp(-((tp-12*Dm+.5/N)/(3*Dm)).^NE);
    x=x+4.5*exp(j*tp*M*(4)*pi+.3).*exp(-((tp-52*Dm+.5/N)/(4*Dm)).^NE);
    x=x+6.5*exp(-j*tp*M*(4)*pi+.3).*exp(-((tp-22*Dm+.5/N)/(4*Dm)).^NE);
    x=x+4.5*exp(j*tp*M*16*pi+0.1).*exp(-((tp-22*Dm+.5/N)/(3*Dm)).^NE);
    x=x+4.8*exp(j*tp*M*42*pi-0.7).*exp(-((tp-26*Dm+.5/N)/(7*Dm)).^NE);
    x=x+5.5*exp(-j*tp*M*16*pi+0.9).*exp(-((tp-20*Dm+.5/N)/Dm).^NE);
    
    x=x+4*exp(j*tp*M*(10)*pi).*exp(-((tp-11*Dm+0.5/N)/Dm).^NE);
    x=x+5*exp(j*tp*M*36*pi).*exp(-((tp-57*Dm+0.5/N)/(4*Dm)).^NE);
    x=x+3*exp(-j*tp*M*42*pi).*exp(-((tp-24*Dm+0.5/N)/(3*Dm)).^NE);
    x=x+4.5*exp(-j*tp*M*52*pi).*exp(-((tp-48*Dm+0.5/N)/Dm).^NE);
    
    x=x+4*exp(-j*tp*M*4*pi).*exp(-((tp-6*Dm+0.5/N)/(Dm)).^NE);
    x=x+5.5*exp(-j*tp*M*36*pi).*exp(-((tp-25*Dm+0.5/N)/(7*Dm)).^NE);
    x=x+6*exp(j*tp*M*48*pi).*exp(-((tp-56*Dm+0.5/N)/(3*Dm)).^NE);
    x=x+4*exp(-j*tp*M*58*pi).*exp(-((tp-50*Dm+0.5/N)/(4*Dm)).^NE);
    
% Final signal with stationary and non-stationary components
    X_INIT=[X_INIT, x];
    
    %% Calculating te STFT on a window by window basis
    STFT(:,mi+1)=fftshift((fft(fftshift(x))));
end
%%-------------------------------------------------------------------------

figure(4), 

sfX=0.85;sfY=0.85; w=20; h=4.5;
set(0,'DefaultAxesPosition',[0.1,0.1,sfX,sfY])
set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontSize',10)

set(gcf,'PaperUnits','centimeters','Units','Centimeters')
p1=get(gcf,'Position');
p2=get(gcf,'PaperPosition');
p1(2)=p1(2)+p1(4)-h/sfY;
p1([3,4])=[w/sfX,h/sfY];
p2([3,4])=[w/sfX,h/sfY];
set(gcf,'Position',p1,'PaperPosition',p2)

subplot(1,2,1), plot(1/N*abs(fftshift(fft(XREZ))));
xlabel({'Frequency samples','(a)'})
ylabel('Amplitudes')
xlim([0,4096])
subplot(1,2,2), plot(1/N*abs(fftshift(fft(X_INIT))));
xlabel({'Frequency samples','(b)'})
ylabel('Amplitudes')
xlim([0,4096])
hold on
aa=1/N*abs(fftshift(fft(X_INIT)));
[m,p]=find(1/N*abs(fftshift(fft(XREZ)))>0.5);
plot(p,aa(p),'r*','markersize',4)
end

