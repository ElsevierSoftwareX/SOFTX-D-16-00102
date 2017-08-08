function figure2

clear;
clc;
close all;

fs = 256;
T = 1/fs;
N = 256;
n = 0:1:(N-1);
snr = 500;

W = 0.45;
% generation of MDPS sequence dictionary
D = MDPSdictionary(N,W,5);
sx = D(:,1).'+D(:,25).'+D(:,690).'+D(:,885).';
varn = sum(sx*sx'/N)/10^(snr/10);
noisx = sqrt(varn)*randn(1,N);   
x = sx + noisx;

numdirac = 50;
diracloc = round(linspace(1,N,numdirac));
xsr = x(diracloc);  
 % basis expansion parameters
[~,gammaMSr,IND] = mp(D(diracloc,:),xsr,'mse',1e-20);
indleng = length(IND);
UMr = D(:,IND);
% MPDS expansion
xrr = UMr*gammaMSr;


x2 = x;
tw = fliplr([1:floor(N/2) -floor(N/2)+1:0])/N;
sig = 0.025;
wind = (1/sqrt(2*pi*sig^2))*exp(-0.5*tw.^2/sig^2);
W = fft(wind); %,size(WWin)
X2F = [fft(x2) fft(x2)];
STFTX = zeros(129,256);
for f=0:floor(N/2)
    STFTX(f+1,:)=ifft(X2F(f+1:f+N).*W);
end

STFTXa = zeros(129,256);
for zz = 1:length(IND)
   x3 =  (UMr(:,zz)*gammaMSr(zz));
   X3F = [fft(x3) fft(x3)];
   aX = zeros(129,256);
   for f=0:floor(N/2)
       aX(f+1,:)=ifft(X3F(f+1:f+N).*W);
   end
   STFTXa = STFTXa + aX;  
end

numdirac = 50;
diracloci = 0;
while length(diracloci)<numdirac
    diracloci = sort(unique(ceil((N-1)*abs(rand(1,ceil(1.4*numdirac))))+1));
end
diracloci = setdiff(diracloci,diracloci(round(linspace(1,length(diracloci),length(diracloci)-numdirac))));
diracloc = diracloci;
xsi = x(diracloc); 
% basis expansion parameters
[B,gammaMSi,IND] = mp(D(diracloc,:),xsi,'mse',1e-20);
indleng = length(IND);
UMi = D(:,IND);
% MPDS expansion
xMSi = UMi*gammaMSi;


STFTXb = zeros(129,256);
for zz = 1:length(IND)
   x4 =  (UMi(:,zz)*gammaMSi(zz));
   X4F = [fft(x4) fft(x4)];
   aX = zeros(129,256);
   for f=0:floor(N/2)
       aX(f+1,:)=ifft(X4F(f+1:f+N).*W);
   end
   STFTXb = STFTXb + aX;  
end

freq = 0:1:127;
subplot(2,2,1), plot(n/256,real(x),'k'), xlabel({'Time (s)', '(a)'}), ylabel('Amplitude'), axis([0 1 -0.38 0.38])
subplot(2,2,2), imagesc((0:1:255)/256,freq,abs(STFTX).^2),  axis xy, colormap(1-gray), xlabel({'Time (s)', '(b)'}), ylabel('Frequency (Hz)'), axis([0 1 0 127])
subplot(2,2,3), imagesc((0:1:255)/256,freq,abs(STFTXa).^2),  axis xy, colormap(1-gray), xlabel({'Time (s)', '(c)'}), ylabel('Frequency (Hz)'), axis([0 1 0 127])
subplot(2,2,4), imagesc((0:1:255)/256,freq,abs(STFTXb).^2),  axis xy, colormap(1-gray), xlabel({'Time (s)', '(d)'}), ylabel('Frequency (Hz)'), axis([0 1 0 127])

end


function Dic = MDPSdictionary(N,W,ad)

if nargin < 2
    error('At least N and W should be passed to the function.') 
else if nargin == 2
        ad = 5;
    end
end

if N*W >= floor(N/2)
    error('Time-bandwidth product must be less than N/2.')
end

n = 0:1:N-1;              % time vector 
dic = zeros(N,1);         % initilization of a dictionary matrix

for pp = 1:ad
    Wp = W/pp;            % splitting of the bandwidth
    NW = N*Wp;            % time-bandwidth product
    K = ceil(2*NW)+1;     % number of dpss for the bands at a particular depth
    [E,V] = dpss(N,NW,K);
    fm =W*(-1+(2*(1:pp)-1)/(pp)); % modulation frequencies for given depth
    for l = 1:length(fm)
        BM1 = (E.' .* repmat(exp(j*2*pi*fm(l)*n),size(E,2),1)).'; % modulation of dpss in a certain band at given depth
        dic = [dic BM1];
    end
end
Dic = dic(:,2:end); % the MDPS dictionary

% End of function MDPSdictionary.
end


function [B,CF,IND] = mp(D,x,cs,cond)

N = length(x);                 % length of the signal
if (N==size(D,1))< 1
    error('The duration of the signal and bases in the dictionary are not the same.')
end

if size(x,2)>size(x,1)        % making sure that the signal is a column vector
    x = x.';
end

R = x;         % initialization of the remainder              
ind = 0;       % initialization of the index vector
coef = 0;      % initialization of the coefficient vector
NB = 0;        % current number of bases in dictionary
denomdiv = sqrt(diag(D'*D));
switch lower(cs)
    case 'mse'        % matching pursuit which yield a desired MSE
        desmse = cond;
        totmse = 1;
        while totmse > desmse && NB<5*N
            tempval = abs(D'*R./denomdiv);
            gamma = find(max(tempval)==tempval,1);
            alpha = D(:,gamma)'*R/(D(:,gamma)'*D(:,gamma));
            R = R-alpha*D(:,gamma);
            coef=[coef; alpha ];
            ind = [ind gamma];
            totmse=(1/N)*sum(abs(x-D(:,ind(2:end))*coef(2:end)).^2)/var(x);
            NB = NB + 1;
        end
        ind = ind(2:end);
        coef = coef(2:end);
%         (1/N)*sum(abs(x-D(:,ind)*coef).^2)/var(x)
    case 'nob'            % matching pursuit which yields desired number of bases in a dictionary
        RNB = cond;  
        while NB<RNB
            tempval = abs(D'*R./denomdiv);
            gamma = find(max(tempval)==tempval,1);
            alpha = D(:,gamma)'*R/(D(:,gamma)'*D(:,gamma));
            R = R-alpha*D(:,gamma);
            coef=[coef; alpha ];
            ind = [ind gamma];
            NB = length(unique(ind(2:end)));
        end
        ind = ind(2:end);
        coef = coef(2:end);
end
IND = unique(ind);               % making sure there are no repeats in indexes
indleng = length(IND);

CF = zeros(indleng,1);
for k = 1:indleng                % the for loop ads the coefficients for the same basis, otherwise leaves unchanged
   if sum(IND(k)==ind)>1
       poseq = (IND(k)==ind);
       CF(k) = sum(coef(poseq));
   else
       pos = (IND(k)==ind);
       CF(k)=coef(pos);
   end
end
B = D(:,IND);
% end of function mp.
end