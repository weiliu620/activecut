function [sigma2_p_hat, sigma2_t_hat, sigma2_eps_hat, ICC_u, ICC_c]...
    = iccsim(sigma2_p, sigma2_t, sigma2_eps, N, K, verbose)

p = sqrt(sigma2_p) * randn(N,1);
t = sqrt(sigma2_t) * randn(1,K);
mu = 1;

Y = repmat(mu, N, K) + repmat(p, [1 K]) + repmat(t, [N 1])...
    + sqrt(sigma2_eps) * randn(N,K); %#ok<RPMT1>
Y_id = mean(Y,2);
Y_dj = mean(Y,1);
Y_dd = sum(sum(Y)) / (N*K);

% SS = SS_p + SS_w
% SS = SS_p + SS_t + SS_e;
% SS_w = SS_t + SS_e;
SS = sum(sum((Y-repmat(Y_dd, [N K])).^2)); % total variance.
SS_p = K * sum((Y_id - Y_dd).^2); % between subjects(row)
SS_w = sum(sum((Y - repmat(Y_id, [1 K])).^2)); % within subject (row).
SS_t = N * sum((Y_dj - Y_dd).^2); % between columns.
SS_e = sum(sum((Y - repmat(Y_id, [1 K]) - repmat(Y_dj, [N 1]) + repmat(Y_dd, [N K])).^2)); % residual.
if (verbose >= 1)
    fprintf('SS=%f, SS_p=%f, SS_w=%f, SS_t=%f, SS_e=%f\n',...
        SS_p, SS_w, SS_t, SS_e);
end;

% mean square.
MS_p = SS_p / (N-1);
MS_t = SS_t / (K-1);
MS_e = SS_e / ((N-1)*(K-1));

% estimates.
sigma2_p_hat = (MS_p - MS_e) / K;
sigma2_t_hat = (MS_t - MS_e) / N;
sigma2_eps_hat = MS_e;

% intraclass correlation.
ICC_u = (MS_p - MS_e) / (MS_p + (K-1)*MS_e + (K/N)*(MS_t-MS_e));
ICC_c = (MS_p - MS_e) / (MS_p + (K-1)*MS_e);


