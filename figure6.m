%% Algorithm for separation of micro-Dopler and rigid body components
%% using the concept of Compressive sensing in STFT domain

function figure6
close all;
clear all;

%% Load micro-Dopler components 
load rotating_reflectors.mat
%% Add some rigid body components
x=TwoRef96r17(1:1280,:)+1*(0.05*exp(j*(1:1280)*(0.4)*pi)'+0.06*exp(j*(1:1280)*0.6*pi)'+0.08*exp(j*(1:1280)*(-0.15)*pi)'+0.0*exp(j*(1:1280)*(-0.44)*pi)'+0.045*exp(j*(1:1280)*0.15*pi)');
% xpolazno=x;
x=x.';
x=[zeros(1,1348), x, zeros(1,1348)];

N=length(x)

m=N;
M=128;
XPOC=[];XREZ=[];

%% Extended for STFT calculation (last window)
x=[x zeros(1,M)];

%% Calculate STFT
nn=1; 
for n=0:1:N-1
     for k=0:M-1
         STFT(nn,k+1)=sum(x(1,n+1:n+M).*exp(-j*2*pi*k*(0:M-1)/M));
      end
     nn=nn+1;
end
STFT=fftshift(abs(STFT),2);
 
%% Remove the side zeros and reduce the size of STFT matrix to ease calculation 
STFT1=STFT(1300:20:2600,1:2:128);
STFT1=STFT1(1:64,:);
STFT=STFT1;
   
  N=4096;
  M=64;
%% Rearanging STFt matrix into vector STFTvec
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
% AAA=[];
% for i=0:N/M-1
% AAA=[AAA; STFTFULL(i*M+1:(i+1)*M)'];
% end
%%-------------------------------------------------------------------------
%% L-estimation procedure applied to sorted STFT values
ii=0.5*M;
STFT_filt=S;
for k=1:M
    STFT_sort(:,k)=sort(S(:,k));
    [aa,bb]=sort(S(:,k));
    STFT_filt(bb(ii:M),k)=0;
    STFT_filt(bb(1:2),k)=0;
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
xprec=real(ifft(xp));
x=xprec';
nn=1;

% Xorig=fftshift(fft(XPOC));
% Xrecovered=xp;
% 
% xorig=XPOC;
% xrecovered=xprec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Illustration of the remaining signal parts after extracting the
% stationary rigid bodu components 
xp1=xp(1:3.2:4096,:)';
a_xp1=abs(xp1);
[poz_a,vr]=find(a_xp1'<150);

a_xp2=(a_xp1);
a_xp2(poz_a)=0;


x_FT_rec=abs(a_xp2)/3.8;
x_rem=ifft(x_FT_rec);

x_rem=[zeros(1,1348), x_rem, zeros(1,1348)];
N_rem=length(x_rem);

m_rem=N_rem;
M_rem=128;

x_rem=[x_rem zeros(1,M_rem)];
w_rem=gausswin(M_rem,2);
nn_rem=1; 


for k_rem=1:length(x_rem)-M_rem
     STFT_rem(:,k_rem)=fft(((x_rem(k_rem+(0:M_rem-1))).*w_rem'));
end

   STFT_rem=STFT_rem'; 

   STFT_rem=fftshift(abs(STFT_rem),2);
   
   
STFT_rem_1=STFT_rem(1300:20:2600,1:2:128);
STFT_rem_1=STFT_rem_1(1:64,:);
STFT_rem_2=fliplr(fftshift(abs(STFT_rem_1),2));
STFT_rem_2=[zeros(64,1) STFT_rem_2(:,1:63)];
STFT_rec=STFT_rem_2;

STFT_rem_2(:,14)=STFT(:,14);
STFT_rem_2(:,28)=STFT(:,28);
STFT_rem_2(:,33)=STFT(:,33);
STFT_rem_2(:,38)=STFT(:,38);
R=STFT-1*STFT_rem_2;


figure,
SetFigureDefaults(20,4.5)
xlabel('Frequency')
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')
subplot(141),imagesc(abs(STFT)), colormap(hot)
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(a)'})
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(142),imagesc((abs(STFT_filt))),colormap(hot)
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(b)'})
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(143),imagesc((abs(STFT_rec))),colormap(hot)
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(c)'})
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

subplot(144),imagesc((abs(R))),colormap(hot)
set(gca,'xtick',[],'ytick',[])
xlabel({'Frequency','(d)'})
ylabel('Time')
text(N,N-32,' ','HorizontalAlignment','left','VerticalAlignment','bottom')

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




