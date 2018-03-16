function res = EmpriLL(F, fs, D, N, A, mu, sigma2)
% 
% Compute the empirical distriburtion(histogram) of the likelihood function
% p(rho|c), i.e. the conditional probability of correlation given the
% connietivity. For now if I generated two noised sequence from same sin
% wave pattern, I assume their connectivity as one, otherwise it is zero.
%
% This routine only compute p(rho | c = 1). For c = 0, we already know the
% correlation is Gaussian distribution after Pearson transform.
% F: frequency of the time series.
% mu: mean of the signal.
% fs: sampling frequency.
% D: number of time points.
% N: total number of time series, including sine signal and constant
% signal.
% A: Amplitude of sin save.
% sigma2: variance of additive noise.

ts = 1/fs; % Sampling interval.
T = (0:D-1)*ts;
X = A*repmat(sin(2*pi*F*T), N, 1) +...
    sigma2 * randn(N, D) + ones(N, D) * mu;
RHO = corr(X');
RHO = RHO - diag(diag(RHO));
B = N*N/2/200; % make sure each bin has 200 data point.
RHOS = reshape(RHO, N*N,1);
close all; figure;
hist(RHOS, B);
snr = num2str(A^2/sigma2);
title(strcat('A^2/sigma^2=', snr));
saveas(gcf, strcat('figures/empriLL', snr, '.eps'), 'eps');


Z = 0.5 * log ((1+RHOS)./(1-RHOS));
figure;
hist(Z, B);
title(strcat('A^2/sigma^2=', snr));
saveas(gcf, strcat('figures/empriLL', snr, 'z', '.eps'), 'eps');


