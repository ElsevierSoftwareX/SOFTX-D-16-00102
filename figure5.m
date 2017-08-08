%% Algorithm for separation of micro-Dopler and rigid body components
%% using the concept of Compressive sensing in STFT domain

function figure5
close all;
clear all;


load rotating_reflectors.mat
x=TwoRef96r17(1:1280,:)+1*(0.05*exp(j*(1:1280)*(0.4)*pi)'+0.06*exp(j*(1:1280)*0.6*pi)'+0.08*exp(j*(1:1280)*(-0.15)*pi)'+0.0*exp(j*(1:1280)*(-0.44)*pi)'+0.045*exp(j*(1:1280)*0.15*pi)');
xpolazno=x;

x=x.';
x=[zeros(1,1348), x, zeros(1,1348)];

N=length(x);

m=N;
M=128;
XPOC=[];
XREZ=[];

x=[x zeros(1,M)];

nn=1; 
for n=0:1:N-1
     for k=0:M-1
         STFT(nn,k+1)=sum(x(1,n+1:n+M).*exp(-j*2*pi*k*(0:M-1)/M));
      end
     nn=nn+1;
end

STFT=fftshift(abs(STFT),2);
STFT1=STFT(1300:20:2600,1:2:128);
STFT1=STFT1(1:64,:);
STFT=STFT1;

figure,
SetFigureDefaults(14,4.5)
xlabel('Frequency')
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')
subplot(121),imagesc(abs(STFT)), colormap(hot)
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(a)'})
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(122),
plot(fftshift(abs(fft(xpolazno)))), xlim([1 length(xpolazno)]), ylim([1 150])
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(b)'})
ylabel('Amplitude')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')


end

