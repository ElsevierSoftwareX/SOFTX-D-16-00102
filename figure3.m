%% Algorithm for separation of stationary and nonstationary signal compoennts using the concept of Compressive sensing in STFT domain

function figure3

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

rand('state',1);
for ii=1:M
ISUM=randperm(M);
    GR=2*round(sum(rand(1,3)));
    STFT(ISUM(1:GR),ii)=STFT(ISUM(1:GR),ii).*(1+rand(GR,1))+6*M*(rand(GR,1)+rand(GR,1));
end

%% Rearanging STFt matrix into vector STFTvec
STFT=STFT';                
STFTvec=[];
for i=1:N/M
STFTvec=[STFTvec STFT(i,:)];
end

%% Defining matrix of DFT basis (of size MxM)
for m=0:M-1
     for k=0:M-1
         W(k+1,m+1)=exp(-j*2*pi*m*k/M);
     end
end
%% Extended NxN matrix 
W_ext = kron(eye(N/M),W);


%% Defining matrix of DFT basis (of size NxN)
 for m=0:N-1
     for k=0:N-1
         W(k+1,m+1)=exp(-j*2*pi*m*k/N);
     end
 end
%% Inverse NxN DFT matrix
WN=inv(W);

%% AFULL matrix defines the linear transform from DFT to STFT
AFULL=W_ext*WN;

STFTFULL=STFTvec';
B=(AFULL);
Binv=inv(AFULL);

S=[];
for i=0:N/M-1
S=[S; STFTFULL(i*M+1:(i+1)*M)'];
end

%%-------------------------------------------------------------------------
%% L-estimation procedure applied to sorted STFT values
ii=0.5*64;
STFT_filt=S;
for k=1:M
    STFT_sort(:,k)=sort(S(:,k));
    [aa,bb]=sort(S(:,k));
    STFT_filt(bb(ii:M),k)=0;
    STFT_filt(bb(1:10),k)=0;
end
%%-------------------------------------------------------------------------

STFTFULL=[];
for i=1:N/M
    STFTFULL=[STFTFULL STFT_filt(i,:)];
end

%%Identifying the positionas of available samples/measurements that
%%corresponds to non-zero elements in filterered STFT version
q=find(STFTFULL~=0);
A=Binv(q,:);
A1=B(q,:);

%Measurement vector
y=STFTFULL(q)';

% Initial transform vector
x0=A1'*y;

%Running the recovery Algorithm from L1-magic
xp=l1eq_pd(x0,A1,[],y,1e-3);

%recovered signal in time domain
xprec=(ifft(fftshift(xp)));

Sx=AFULL*xp;
STFT_recovered=[];
for i=0:N/M-1
STFT_recovered=[STFT_recovered; Sx(i*M+1:(i+1)*M)'];
end


figure(1), colormap(hot)

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

set(gca,'xtick',[],'ytick',[])
xlabel('frequency')
ylabel('time')
text(N,N-64,' ','HorizontalAlignment','left','VerticalAlignment','bottom')
subplot(1,4,1),imagesc(abs(STFT).^1), colormap(1-gray),
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency samples','(a)'})
ylabel('Time samples')
text(N,N-64,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(142),imagesc((abs(STFT_sort)).^1),
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency samples','(b)'})
ylabel('Time samples')
text(N,N-64,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(143),imagesc((abs(STFT_filt).^1)),
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency samples','(c)'})
ylabel('Time samples')
text(N,N-64,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(1,4,4),imagesc((abs(STFT_recovered).^1)), axis xy
set(gca,'xtick',[],'ytick',[]), 
xlabel({'Frequency samples','(d)'})
ylabel('Time samples')
text(N,N-64,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

end



% l1eq_pd.m
%
% Solve
% min_x ||x||_1  s.t.  Ax = b
%
% Recast as linear program
% min_{x,u} sum(u)  s.t.  -u <= x <= u,  Ax=b
% and use primal-dual interior point method
%
% Usage: xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if
%     the duality gap is less than pdtol).  
%     Default = 1e-3.
%
% pdmaxiter - Maximum number of primal-dual iterations.  
%     Default = 50.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)

largescale = isa(A,'function_handle');

if (nargin < 5), pdtol = 1e-3;  end
if (nargin < 6), pdmaxiter = 150;  end
if (nargin < 7), cgtol = 1e-8;  end
if (nargin < 8), cgmaxiter = 200;  end

N = length(x0);

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(N,1); ones(N,1)];

% starting point --- make sure that it is feasible
if (largescale)
  if (norm(A(x0)-b)/norm(b) > cgtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w, cgres, cgiter] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b)/norm(b) > cgtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    opts.POSDEF = true; opts.SYM = true;
%     [w, hcond] = linsolve(A*A', b, opts);
    [w, hcond] = linsolve(A*A', b);
    if (hcond < 1e-14)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = A'*w;
  end  
end
x = x0;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));

% set up for the first iteration
fu1 = x - u;
fu2 = -x - u;
lamu1 = -1./fu1;
lamu2 = -1./fu2;
if (largescale)
  v = -A(lamu1-lamu2);
  Atv = At(v);
  rpri = A(x) - b;
else
  v = -A*(lamu1-lamu2);
  Atv = A'*v;
  rpri = A*x - b;
end

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
done = (sdg < pdtol) | (pditer >= pdmaxiter);
while (~done)
  
  pditer = pditer + 1;
  
  w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
  w3 = -rpri;
  
  sig1 = -lamu1./fu1 - lamu2./fu2;
  sig2 = lamu1./fu1 - lamu2./fu2;
  sigx = sig1 - sig2.^2./sig1;
  
  if (largescale)
    w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));
    h11pfun = @(z) -A(1./sigx.*At(z));
    [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;
    Adx = A(dx);
    Atdv = At(dv);
  else
    w1p = -(w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1)));
    H11p = A*(sparse(diag(1./sigx))*A');
    opts.POSDEF = true; opts.SYM = true;
%     [dv,hcond] = linsolve(H11p, w1p, opts);
    [dv,hcond] = linsolve(H11p, w1p);
    
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;
    Adx = A*dx;
    Atdv = A'*dv;
  end
  
  du = (w2 - sig2.*dx)./sig1;
  
  dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
  dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;
  
  % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
  indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);
  s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
  indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);
  s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);
  
  % backtracking line search
  suffdec = 0;
  backiter = 0;
  while (~suffdec)
    xp = x + s*dx;  up = u + s*du; 
    vp = v + s*dv;  Atvp = Atv + s*Atdv; 
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up;  fu2p = -xp - up;  
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    suffdec = (norm([rdp; rcp; rpp]) <= (1-alpha*s)*resnorm);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
%       disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')
      xp = x;
      return
    end
  end
  
  
  % next iteration
  x = xp;  u = up;
  v = vp;  Atv = Atvp; 
  lamu1 = lamu1p;  lamu2 = lamu2p;
  fu1 = fu1p;  fu2 = fu2p;
  
  % surrogate duality gap
  sdg = -(fu1'*lamu1 + fu2'*lamu2);
  tau = mu*2*N/sdg;
  rpri = rpp;
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
  rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
  resnorm = norm([rdual; rcent; rpri]);
  
  done = (sdg < pdtol) | (pditer >= pdmaxiter);
  
 
end

end
