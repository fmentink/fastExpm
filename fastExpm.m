function P = fastExpm(H,varargin)
% This function efficiently implements matrix exponential matrix for sparse and full matrices. 
% This code is based on scaling, taylor series and scaling.
% Works for both CPU and GPU (with appropriate toolbox). 
% Use: P = fastExpm(H) %with default settings
%      P = fastExpm(H,convergenceCriteria)
%      P = fastExpm(H,convergenceCriteria, nonZeroTol)
%
% Two criteria are used to speed the computation and preserve the sparsity.
% [1] convergenceCriteria determine the threshold for the Taylor series (default 1e-6)
% [2] nonZeroTol strips elements smaller than nonZeroTol at each computation step to preserve sparsity, default is 1e-14.
% The code automatically switches from sparse to full if sparsity is below 15 % to maintain speed.
% If H is a gpuArray, P will be computed on the GPU.
%
% This code has been originally developed by Ilya Kuprov (http://spindynamics.org/) and has been adapted by F. Mentink-Vigier (fmentink@magnet.fsu.edu)
% If you use this code, please cite 
%  - H. J. Hogben, M. Krzystyniak, G. T. P. Charnock, P. J. Hore and I. Kuprov, Spinach – A software library for simulation of spin dynamics in large spin systems, J. Magn. Reson., 2011, 208, 179–194.
%  - I. Kuprov, Diagonalization-free implementation of spin relaxation theory for large spin systems., J. Magn. Reson., 2011, 209, 31–38.


% input check: is the user specifying its own criteria?
if nargin>1
    convergenceCriteria=varargin(2); % user specified convergenceCriteria
else
    convergenceCriteria=1e-6;
end

if nargin>2
    nonZeroTol=varargin(2); % user specified nonZeroTol
else
    nonZeroTol=1e-14;
end

%% Scaling H
mat_norm=norm(H,1); delta=1;
n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings;
H=H*scaling_factor^-1;
H=nonZeroTol*round((1/nonZeroTol)*H);

%% Run Taylor series procedure on the CPU/GPU
P=0*H; % this ensures H and P are of the same type --> CPU/GPU
P=speye(size(H)); nextTerm=P; n=1; 

% Sparsity ans size check
if nnz(H)/numel(H)>0.25 && numel(H)^0.5<=64
    H=full(H);
else
    H=sparse(H);
end

while delta>convergenceCriteria
    % Compute the next term
    if issparse(nextTerm)
        nextTerm=(1/n)*H*nextTerm; % order matters
        nextTerm=nonZeroTol*round((1/nonZeroTol)*nextTerm); % Eliminate small elements
        if nnz(nextTerm)/numel(nextTerm)>0.25
            nextTerm=full(nextTerm);
        end
    else
        nextTerm=(1/n)*nextTerm*H;
    end
    delta=norm(nextTerm,1); % check residual norm
    P=P+nextTerm; n=n+1; %Add to the total and increment the counter
end

% Reduce memory footprint: clear('H','nextTerm') could be used but may be slower. 
H=[];nextTerm=[];

%% Squaring
P=nonZeroTol*round((1/nonZeroTol)*P);
for n=1:n_squarings
    P=P*P;
    if issparse(P)
        if nnz(P)/numel(P)<0.25
            P=nonZeroTol*round((1/nonZeroTol)*P); % Eliminate small elements
        else
            P=full(P);
        end
    end
end
end
