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
[U0_x, S0_x, V0_x] = svd(X0, 0);
[U0_y, S0_y, V0_y] = svd(Y0, 0);

% project data to principal comp.
cca_var_best = 0;
var_x = sum(diag(S0_x).^2);
var_y = sum(diag(S0_y).^2);
for n_comp = n_comp_cca:n_comp_max
    X1 = X0(:, 1:n_comp);
    Y1 = Y0(:, 1:n_comp);
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
    var_exp = (var_xby + var_ybx) / 2;
    [var_exp_sorted, idx] = sort(var_exp, 'descend');
    
    % truncate to n_comp_cca columns.
    var_exp_sorted = var_exp_sorted(1:n_comp_cca);
    A = A(:, idx(1:n_comp_cca));
    B = B(:, idx(1:n_comp_cca));
    R = R(idx(1:n_comp_cca));
    S = S(:, idx(1:n_comp_cca));
    
    % the first n_comp_cca components explain how much variance.
    if sum(var_exp_sorted) > cca_var_best
        cca_var_best = sum(var_exp_sorted);
        A_best = A;
        B_best = B;
        R_best = R;
        T_best = var_exp_sorted;
        n_comp_best = n_comp;
        fprintf('num of PC: %d:\n', n_comp);
        fprintf('    X: PC explains: %f, %d CC explains: %f, %d CC explains: %f, Y explains: %f\n',...
            sum(diag(S0_x(1:n_comp, 1:n_comp)).^2) / var_x,...
            n_comp, sum(C_x),...
            n_comp_cca, sum((diag(S1_x).^2)) / var_x,...
            sum(diag(S)) / var_x,...
            sum(var_xby) );
        fprintf('    Y: PC explains: %f, CC explains: %f, X explains: %f\n',...
            sum(diag(S0_y(1:n_comp, 1:n_comp)).^2) / var_y,...
            sum((diag(S1_y).^2)) / var_x,...
            sum(C_y) / var_y,...
            sum(var_exp_sorted) );
    end;
end;

A = A_best;
B = B_best;
R = R_best;
T = T_best;

disp(var_xby(idx));
disp(var_ybx(idx))

% Compute the canonical variates
A = V0_x(:, 1:n_comp_best) * A;
B = V0_y(:, 1:n_comp_best) * B;
if nargout > 4
    Xp = X0 * A;
    Yp = Y0 * B;
end
