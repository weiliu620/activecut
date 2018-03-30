function [A, B, R, Xp, Yp] = mycanoncorr_clean(X0, Y0, n_comp)

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

[U0_x, S0_x, V0_x] = svd(X0, 0);
[U0_y, S0_y, V0_y] = svd(Y0, 0);
[U, S, V] = svd(U0_x(:, 1:n_comp)' * U0_y(:, 1:n_comp));
R = diag(S);
A = V0_x(:, 1:n_comp) / S0_x(1:n_comp, 1:n_comp) * U;
V0_x(:, 1:n_comp) / S0_x(1:n_comp, 1:n_comp)
B = V0_y(:, 1:n_comp) / S0_y(1:n_comp, 1:n_comp) * V;

% Compute the canonical variates
if nargout > 3
    Xp = X0 * A;
    Yp = Y0 * B;
end