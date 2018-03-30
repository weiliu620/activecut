function Xr = remove_pc(X, n_comp)
% remove the first n_comp principal components from the data. output is the
% reconstructed data with the components removed. 

[n_samples, n_features] = size(X);
col_mean = mean(X, 1);
X = X - repmat(col_mean, n_samples, 1);
[U, S, V] = svd(X);
disp(diag(S(1:5, 1:5)))
% project the data to new space. Y = XV = (USV')V = US
Y = U * S;
% remove components.
if n_comp > 0
    Y(:, 1:n_comp) = 0;
end;
% reconstruction. Xr = YV'
Xr = Y * V' + repmat(col_mean, n_samples, 1);


