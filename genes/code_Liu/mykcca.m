function [nalpha, nbeta, r] = mykcca(X, Y, kappa, varargin)
% mykcca(X, Y, n_comp)
%
% Kernel CCA by using Hardoon's code.
% % The value of kernel determines the used kernel. Possible values are
% 'linear', 'gauss', 'poly', 'subsets', or 'princ_angles' (default = 'gauss'). 

addpath ~/packages/kcca_package;
addpath ~/packages/drtoolbox/techniques;

[N, P1] = size(X);
[~, P2] = size(Y);

% centering.
X = X - repmat(mean(X,1), N, 1);
Y = Y - repmat(mean(Y,1), N, 1);

% compute kernel (Gram) matrix.

if nargin == 3
    kernel = 'linear';
    param1 = 0;
    param2 = 0;
elseif nargin > 3
    kernel = varargin{1};
    if length(varargin) > 1 & strcmp(class(varargin{2}), 'double'), param1 = varargin{2}; end
    if length(varargin) > 2 & strcmp(class(varargin{3}), 'double'), param2 = varargin{3}; end
end

Kx = gram(X, X, kernel, param1, param2);
tau = 1e-6;
while ~ispd(Kx)
    tau = tau*2;
    Kx = Kx + tau * eye(N);
    fprintf('Kx forced to positive-definite by tau = %f\n', tau);
end;
    
Ky = gram(Y, Y, kernel, param1, param2);
tau = 1e-6;
while ~ispd(Ky)
    tau = tau*2;
    Ky = Ky + tau * eye(N);
    fprintf('Ky forced to positive-definite by tau = %f\n', tau);
end;

[nalpha, nbeta, r, Kx, Ky] = kcanonca_reg_ver2(Kx, Ky, 10e-6, kappa, 1, 2);

