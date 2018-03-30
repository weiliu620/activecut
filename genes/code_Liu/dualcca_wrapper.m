function [vs, cors] = dualcca_wrapper(X, Y, gamma)
% dualcca_wrapper(X, Y, gamma)
% 
% Wrapper function for dualcca.m function (from the book 'Kernel methods for
% Pattern Analysis').

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
    if length(varargin) > 1 & strcmp(class(varargin{2}), 'double'), param1 = ...
            varargin{2}; 
    end
    if length(varargin) > 2 & strcmp(class(varargin{3}), 'double'), param2 = ...
            varargin{3}; 
    end
end

    
Kx = gram(X, X, kernel, param1, param2);
Ky = gram(Y, Y, kernel, param1, param2);

% construct the cell array required by dualcca
Ks{1} = Kx;
Ks{2} = Ky;

[vs, cors] = dualcca(Ks, gamma, 1000);