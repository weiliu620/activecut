function [A, B, R, T, Xp, Yp] = mycanoncorr_dirty(X0,Y0, n_comp_cca, n_comp_max)
% [A, B, R, T] = mycanoncorr_dirty(X0, Y0, n_comp)
% [A, B, R, T, Xp, Yp] = mycanoncorr_dirty(X0, Y0, n_comp)
%
% CCA with PCA as preprocessing. 
%
% X0, Y0: data matrix.
% n_comp: num of principal componnets that will be kept in preprocessing.
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

tstart = tic;
% Center the variables
X0 = X0 - repmat(mean(X0,1), n, 1);
Y0 = Y0 - repmat(mean(Y0,1), n, 1);

% First do a PCA preprocessing.
cache_filename = 'mycanoncorr_dirty_cache.mat';
if exist(cache_filename)
    load cache_filename
else
    [U0_x, S0_x, V0_x] = svd(X0, 0);
    [U0_y, S0_y, V0_y] = svd(Y0, 0);
    save(cache_filename, 'U0_x', 'S0_x', 'V0_x', 'U0_y', 'S0_y', 'V0_y');
    fprintf('cache file save to %s.\n', cache_filename);

% project data to principal comp.
cca_var_best = 0;
var_x = sum(diag(S0_x).^2);
var_y = sum(diag(S0_y).^2);
for n_comp = n_comp_cca:10:n_comp_max
    X1 = X0 * V0_x(:, 1:n_comp);
    Y1 = Y0 * V0_y(:, 1:n_comp);
    [U1_x, S1_x, V1_x] = svd(X1, 0);
    [U1_y, S1_y, V1_y] = svd(Y1, 0);
    [U, S, V] = svd(U1_x' * U1_y);
    R = diag(S);
    A = V1_x / S1_x * U;
    B = V1_y / S1_y * V;
    
    C_x = zeros(n_comp, 1);
    for i = 1:n_comp
        C_x(i) = sum((diag(S1_x) .* U(:,i)).^2);
    end;
    
    C_y = zeros(n_comp, 1);
    for i = 1:n_comp
        C_y(i) = sum((diag(S1_y) .* V(:,i)).^2);
    end;
    
    % percentage of variance explained by each componenent.
    var_xby = C_x .* diag(S).^2 / var_x; % X explained by Y.
    var_ybx = C_y .* diag(S).^2 / var_y; % Y explained by X.
    [~, idx] = sort(var_ybx, 'descend');
    var_xby = var_xby(idx);
    var_xby = var_xby(1:n_comp_cca);
    var_ybx = var_ybx(idx);
    var_ybx = var_ybx(1:n_comp_cca);
    A = A(:, idx(1:n_comp_cca));
    B = B(:, idx(1:n_comp_cca));
    R = R(idx(1:n_comp_cca));
    
    % the first n_comp_cca components explain how much variance.
    if sum(var_ybx) > cca_var_best
        fprintf('n_comp: %d, Y var explained by X: %f > %f\n', n_comp, sum(var_ybx), cca_var_best);
        cca_var_best = sum(var_ybx);
        A_best = A;
        B_best = B;
        R_best = R;
        T_best = [var_xby, var_ybx];
        n_comp_best = n_comp;
    end;
    
    fprintf('num of PC: %d. average corr: %.3f\n', n_comp, mean(R));
    fprintf('    X: %d PC explains: %.3f, %d CC explains: %.3f, %d CC explains: %.3f, Y explains: %.3f\n',...
            n_comp, sum(diag(S0_x(1:n_comp, 1:n_comp)).^2) / var_x,...
            n_comp, sum(C_x) / var_x,...
            n_comp_cca, sum(C_x(idx(1:n_comp_cca))) / var_x,...
            sum(var_xby) );
    fprintf('    Y: %d PC explains: %.3f, %d CC explains: %.3f, %d CC explains: %.3f, X explains: %.3f\n',...
            n_comp, sum(diag(S0_y(1:n_comp, 1:n_comp)).^2) / var_y,...
            n_comp, sum(C_y) / var_y,...
            n_comp_cca, sum(C_y(idx(1:n_comp_cca))) / var_y,...
            sum(var_ybx) );
end;

A = A_best;
B = B_best;
R = R_best;
T = T_best;

% Compute the canonical variates
A = V0_x(:, 1:n_comp_best) * A;
B = V0_y(:, 1:n_comp_best) * B;
if nargout > 4
    Xp = X0 * A;
    Yp = Y0 * B;
end
