function res = icctest(sigma2_p, sigma2_t, sigma2_e)

% icctest(sigma2_p, sigma2_t, sigma2_e)
% sigma2_p: variance of row effect. sigma2_t is variance of column effect,
% and sigma2_e is variance of error term.

ICC_u = sigma2_p / (sigma2_p + sigma2_t + sigma2_e);
ICC_c = sigma2_p / (sigma2_p + sigma2_e);
N = 25;
K = 3;

reptime = 1000;
estmat = zeros(reptime, 5);
for m = 1:reptime
    [sigma2_p_hat, sigma2_t_hat, sigma2_e_hat, ICC_u_hat, ICC_c_hat] = iccsim(sigma2_p, sigma2_t, sigma2_e, N, K, 0);
    estmat(m,:) = [sigma2_p_hat sigma2_t_hat sigma2_e_hat ICC_u_hat ICC_c_hat];
end;

% mean squares of the estimats compared with true value.
mse = mean((estmat - repmat([sigma2_p sigma2_t sigma2_e ICC_u ICC_c], reptime,1)).^2,1);

truepar = [sigma2_p sigma2_t sigma2_e ICC_u ICC_c];
est_mean = mean(estmat,1);

truepar
est_mean
mse


    