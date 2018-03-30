function [A, B, R, T, Xp, Yp] = mycanoncorr(X0,Y0, n_comp_X, n_comp_Y)
% [A, B, R, T] = mycanoncorr_dirty(X0, Y0, n_comp_X, n_comp_Y)
% [A, B, R, T, Xp, Yp] = mycanoncorr_dirty(X0, Y0, n_comp, n_comp_X, n_comp_Y)
%
% CCA with PCA as preprocessing. 
%
% X0, Y0: data matrix.
% n_comp_X, n_comp_Y: num of principal componnets that will be kept in preprocessing.
% A, B: projection matrix. Columns are orthogonal. Also includes the PCA's
% projection matrix.
% R: correlation between projected variates.
% T: explained variance. T = C*R.^2, where C is the variance of principal
% components, and D is the correlation between projected variates.
% Xp, Yp: projected data.
% 

[n,p1] = size(X0);
if size(Y0,1) ~= n
    error(message('stats:canoncorr:InputSizeMismatch'));
elseif n == 1
    error(message('stats:canoncorr:NotEnoughData'));
end
p2 = size(Y0,2);

% Center the variables
X0 = X0 - repmat(mean(X0,1), n, 1);
Y0 = Y0 - repmat(mean(Y0,1), n, 1);

n_comp_cca = min([n_comp_X, n_comp_Y]);
% First do a PCA preprocessing.
cache_filename = 'mycanoncorr_cache.mat';
if exist(cache_filename)
    load(cache_filename);
else
    [U0_x, S0_x, V0_x] = svd(X0, 0);
    [U0_y, S0_y, V0_y] = svd(Y0, 0);
    save(cache_filename, 'U0_x', 'S0_x', 'V0_x', 'U0_y', 'S0_y', 'V0_y', '-v7.3');
    fprintf('cache file save to %s.\n', cache_filename);
end;

% project data to principal comp.
var_x = sum(diag(S0_x).^2);
var_y = sum(diag(S0_y).^2);

X = X0 * V0_x(:, 1:n_comp_X);
Y = Y0 * V0_y(:, 1:n_comp_Y);
[U_x, S_x, V_x] = svd(X, 0);
[U_y, S_y, V_y] = svd(Y, 0);
[U, S, V] = svd(U_x' * U_y, 0);
R = diag(S);
A = V_x / S_x * U;
B = V_y / S_y * V;
    
C_x = zeros(n_comp_cca, 1);
for i = 1:n_comp_cca
    C_x(i) = sum((diag(S_x) .* U(:,i)).^2);
end;

C_y = zeros(n_comp_cca, 1);
for i = 1:n_comp_cca
    C_y(i) = sum((diag(S_y) .* V(:,i)).^2);
end;
    
% percentage of variance explained by each componenent.
var_xby = C_x .* diag(S).^2 / var_x; % X explained by Y.
var_ybx = C_y .* diag(S).^2 / var_y; % Y explained by X.

fprintf('    X: CC explains: %.3f and %.3f. Y explains: %.3f\n',...
        sum(diag(S0_x(1:n_comp_cca, 1:n_comp_cca)).^2) / var_x,...
        sum(C_x) / var_x,...
        sum(var_xby) );

fprintf('    Y: CC explains: %.3f and %.3f, X explains: %.3f\n',...
        sum(diag(S0_y(1:n_comp_cca, 1:n_comp_cca)).^2) / var_y,...
        sum(C_y) / var_y,...
        sum(var_ybx) );

% Compute the canonical variates
A = V0_x(:, 1:n_comp_X) * A;
B = V0_y(:, 1:n_comp_Y) * B;
T = var_ybx;
if nargout > 4
    Xp = X0 * A;
    Yp = Y0 * B;
end
