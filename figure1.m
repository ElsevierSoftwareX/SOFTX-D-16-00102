%  Compressive sensing in the ambiguity domain
%
%   The function calculates the CS-based Wigner distribution from the set
%   of random measurements taken from the corresponding ambiguity function
%
%   INPUTS: 
%   1. Standard Wigner distribution of mono-component signal
%   (WDstandard_knjiga_mono.mat); the size of WD matrix is 60x60
%   2. Transform domain matrix - 2D DFT matrix (fftmat2_60.mat, of size 3600x3600); 
%
%   The random selection of ambigutiy domain points inside the desired
%   central region in the ambiguity domain is done by using MASK
%  The reconstrcution of points of the time-frequency plane are done using
%  soft thresholding method (Copyright (c) 2010 Gabriel Peyre, see below)

% -----------------------------------------------------------------------

function figure1

clear;
clc;
close all;

% Load matrix corrresponding to the standard Wigner distribution of
% a monocomponent signal (matrix size 60x60)
load WDstandard_knjiga_mono.mat

% Define 2D mask corresponding to a small central region in the ambiguity
% domain
mask2 = zeros(size(WD));
mask2(25:35,25:35)=1;  

% ------------------------------------------------------------------------
% Illustrate the standard WD distribution and the corresponding 
% Ambiguity function obtained as 2D inverse DFT of the WD 

figure(1),
FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 800, 200]);

subplot(1,3,1)
Amb=fftshift(ifft2(WD));
imagesc(abs(Amb)); colormap(1-gray),
set(gca,'YDir','normal'); grid on
xlabel({'Time lag','(a)'}); ylabel('Frequency lag')
subplot(1,3,2),
imagesc(flipud(abs(WD'))); colormap(1-gray),
set(gca,'YDir','normal'); grid on
xlabel({'Time','(b)'}); ylabel('Frequency')
% ------------------------------------------------------------------------

% Reshape the 2D mask to 1D form in order to work with vectors instead of
% matrices
[m n] = size(mask2);
mask_vec = reshape(mask2,m*n,1);
% Positions of elements with value 1 within the mask corresponds to the
% positions of measurements in the ambiguity domain
mask_pos = find(mask_vec ==1);

% Selecting the measuremenst from the ambiguity domain that corresponds to
% the positions in vector mask_pos
tf2c = Amb;
tf2c_sparse = tf2c(mask_pos);

%-------------------------------------------------------------------------
%% PREPARING CS MATRIX
% Transform domain matrix (2D DFT of size 3600x3600)
load fftmat2_60.mat
N=60;

% Inverse 2D DFT matrix
ifftmat = conj(fftmat)/N/N;
clear fftmat

% Selecting rows in Inverse 2D DFt matrix that corresponds to the selected
% measuremenst defined by the mask_pos ;
ifftmat_sparse = ifftmat(mask_pos,:);

save cs_array_mono tf2c_sparse ifftmat_sparse

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------

%% PERFORMING THE RECONSTRUCTION OF SPARSE TIME-FREQUENCY REPRESENTATION FROM THE SELECTED MEASUREMENTS FROM THE AMBIGUITY DOMAIN ;

% load the measurement vector (tf2c_sparse) and the CS matrix
% (ifftmat_sparse) [in this example these two input vector are saved in
% cs_array_mono in the previous part of the function];
load cs_array_mono

n = size(ifftmat_sparse,2);
p = size(ifftmat_sparse,1);

ns = floor(sqrt(n));

% Regularization parameter.

lambda=0.00001;

% CS Matrix A and observation/measurements vector y.

A = ifftmat_sparse;
y = tf2c_sparse;

% List of benchmarked algorithms.
methods = {'fb', 'fista', 'nesterov'};


% operator callbacks

F = @(x)lambda*norm(x,1);
G = @(x)1/2*norm(y-A*x)^2;


% Proximal operator of F. 

ProxF = @(x,tau)perform_soft_thresholding(x,lambda*tau);


% Gradient operator of G.
GradG = @(x)A'*(A*x-y);


% Lipschitz constant.
L = norm(A)^2;

% Bench the algorithm

options.niter = 5000;
E = [];
for i=1:1 %length(methods)
    options.method = methods{i};
    [x,e] = perform_fb(zeros(n,1), ProxF, GradG, L, options);
    E(:,i) = e(:);
% Illustrate the reconstructed sparse tiem-frequency representation 
    subplot(1,3,3)
   imagesc(0:ns-1, [0:ns/2-1]/ns,flipud((abs(reshape(x,ns,ns)')))), colormap(1-gray),
    set(gca,'YDir','normal'); grid on
    xlabel({'Time','(c)'}); ylabel('Frequency')    
end

end





function y = perform_soft_thresholding(x, tau)

% perform_soft_thresholding - soft thresholding
%
%   y = perform_soft_thresholding(x, tau);
%
%   y = prox_{tau*|.|_1}(x) = max(0,1-tau/|x|)*x
%
%   Proximal operator for the scalar L1 norm.
%
%   Copyright (c) 2010 Gabriel Peyre

y = max(0,1-tau./max(abs(x),1e-10)).*x;
end

function [x,R] = perform_fb(x, ProxF, GradG, L, options)

% perform_admm - preconditionned ADMM method
%
%    [x,R] = perform_fb(x, ProxF, GradG, L, options);
%
%   Solves min_x g(x) + f(x)
%   where g is a smooth convex proper function and f is a
%   convex proper function with an easy to compute proximal operator.
%
%   Use several first order-scheme depending on options.method:
%       options.method = 'fb' : classical Foward-backward
%       options.method = 'fista' : FISTA method of Beck and Teboule
%       options.method = 'nesterov' : Nesterov scheme.
%
%   INPUTS:
%   ProxF(y,sigma) computes Prox_{sigma*F}(x)
%   GradG(x) computes \nabla f(x)
%   L is the lipschitz constant of the gradient, if g is C^2:
%       L = max_x norm( Hg(x) ) 
%       where Hg(x) is the hessian of g at x. 
%       For instance, if g(x)=1/2*|A*x-y|^2 then tau = norm(A)^2.
%   options.niter is the number of iterations.
%   options.verb is for the diaplay of iterations.
%   options.report(x) is a function to fill in R.
%
%   OUTPUTS:
%   x is the final solution.
%   R(i) = options.report(x) at iteration i.
%
%   Copyright (c) 2010 Gabriel Peyre

options.null = 0;
method = getoptions(options, 'method', 'fb');
report = getoptions(options, 'report', @(x)0);
niter = getoptions(options, 'niter', 100);
verb = getoptions(options, 'verb', 1);
fbdamping = getoptions(options, 'fbdamping', 1.8);

clear R;
t = 1;  % fista & nesterov
tt = 2/L; gg = 0; A = 0; % nesterov
y = x;
x0 = x;
for i=1:niter 
  	R(i) = report(x);
    if verb
%         progressbar(i,niter);
    end
    switch method
        case 'fb'
            x = ProxF( x-fbdamping/L*GradG(x), fbdamping/L );
        case 'fista'
            xnew = ProxF( y - 1/L*GradG(y), 1/L );
            tnew = (1+sqrt(1+4*t^2))/2;
            y = xnew + (t-1)/(tnew)*(xnew-x);
            x = xnew; t = tnew;
        case 'nesterov'
            a = (tt + sqrt(tt^2 + 4*tt*A))/2;
            v = ProxF( x0-gg, A );
            z = (A*x+a*v)/(A+a);
            x = ProxF( z - 1/L*GradG(z) , 1/L  );
            gg = gg +  a * GradG(x); % P'*(P*x-y);
            A = A + a;
        otherwise
            error('Unknown method');
            
    end      
end

end

function v = getoptions(options, name, v, mendatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0, mendatory);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<3
    error('Not enough arguments.');
end
if nargin<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mendatory
    error(['You have to provide options.' name '.']);
end 
end



