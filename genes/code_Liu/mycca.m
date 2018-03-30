function mycca(X, Y, n_comp, out_file)
% mycca(X, Y, n_comp, out_file)
%
% find the n_comp components that maximize the correlation of projected X and
% Y data matrix. X and Y are samples by features.

[n_samples, n_featuresx] = size(X);
[~, n_featuresy] = size(Y);

% centering.
X = X - repmat(mean(X,1), n_samples, 1);
Y = Y - repmat(mean(Y,1), n_samples, 1);
tic;
[A, B, R, U, V] = canoncorr(X, Y);
toc;

save(out_file, 'A', 'B', 'R', 'U', 'V');


