function [X, C] = GenTimeSeries(snr, dir, smooth, FSIZE, fwhm)
% Usage: GenTimeSeries(snr, dir, smooth, FSIZE, fwhm)
% Generate time series signal.
% F: frequency of the time series. I assume this is mu + sine wave + noise.
% mu: mean of the signal.
% fs: sampling frequency.
% D: number of time points.
% N: total number of time series, including sine signal and constant
% signal.
% label: Nx1 vector. Each element is 1 for sigmal, 0 for noise.
% snr: signal to noise ratio. Defined as A/sigma^2.
%
% output: 
% X: NxD data matrix. Each row is time course.
% C: true connectivity matrix. c_ij is the theoretical correlation
% coefficient between point i and j.

F = 0.2; % frequency of sin wave
mu = 800;
fs = 1; % sampling frequency.
D = 300;    % number of time points sampled.
N = 100;    % number of all time course.
label = zeros(N, 1);
%label(1:M/2) = 1; label(N/2+1:N/2+M/2) = 1;
label(21:30) = 1; label(71:80) = 1;
M = sum(label); % number of sin wave signals
A = 20;

ts = 1/fs;
T = (1:D)*ts;
X = zeros(N, D); Y = zeros(N, D);
sigma2 = A/snr;

X(label==1, :) = repmat(A*sin(2*pi*F*T), M, 1);
Y(label==1, :) = repmat(A*sin(2*pi*F*T), M, 1);
X = X + sigma2* randn(N, D) + ones(N, D) * mu;
Y = Y + sigma2* randn(N, D) + ones(N, D) * mu;

RNS = mycorr(X', Y');
if strcmp(smooth, 'smooth')
    X = GaussianSmooth(X, FSIZE, fwhm);
    Y = GaussianSmooth(Y, FSIZE, fwhm);
end
C = label * label'; % true labeling.
%R = corr(X', Y'); % sample correlation
R = mycorr(X', Y');
%save sinSample R snr;
save(strcat(dir, 'sinSample.mat'), 'RNS', 'R', 'snr');

% Show data
figure(1); imshow(C, 'InitialMagnification', 'fit');
%title('Ground truth Connectivity Map.');
mysave(1, strcat(dir, 'trueConn.eps'));

SF = (A + 3*sqrt(sigma2))*2;
K = 40; % num of points to draw.
Z = (X(1:K, :) - mu * ones(K,D))/SF + repmat((1:K)', 1, D);
% figure(2);  plot(Z');
% title(strcat('Simulated time series. Total time series #:', int2str(N)));
% mysave(2, strcat(dir, '1DTimeSeries.eps'));
fprintf('Signal: %d-- %d, %d -- %d\n', 1, M/2, 1+N/2, N/2+M/2);

figure(3);
imshow(R/2+0.5, 'InitialMagnification', 'fit');
%title('Sample Correlation. Scaled to (0,1)');
mysave(3, strcat(dir, 'SampleCorr.eps'));
saveas(3, strcat(dir, 'SampleCorr.png'));

B = (N*N)/200; % bin number.
figure(4);
sel = true(N,N);
hist(R(sel),B);
title('Sample correlation Histogram'); 
mysave(4, strcat(dir, 'CorrHist.eps'));

figure(5);
imshow(RNS/2+0.5, 'InitialMagnification', 'fit');
title('sample correlatoin on Non-smoothed data, scale to (0,1)');
saveas(5, strcat(dir, 'corrNoSmooth.png'));
A1 = ThresholdCorr(RNS, C);
figure(6);
imshow(A1, 'InitialMagnification', 'fit');
title('Non-Smoothed data. correlation thresholded such that false positivie > 5%');

A2 = ThresholdCorr(R, C);
figure(7);
imshow(A2, 'InitialMagnification', 'fit');
title('smoothed data. correlation thresholded such that false positivie > 5%');



function res = GaussianSmooth(X, FSIZE, fwhm)
% Gaussian smooth on 1-D signal. X is NxD vector, F is Gaussian filter.
sigma = fwhm/2.3548;
idx = -(FSIZE-1)/2:(FSIZE-1)/2;
F = exp(-idx.^2/(2*sigma^2));
F = F/sum(F); % normalize.
[N, D] = size(X);
res = zeros(N, D);
for d = 1:D
    res(:,d) = conv(X(:,d), F, 'same');
    

end;

function A = ThresholdCorr(R, trueMap)
[N1, N2] = size(R);
T = 1;

FP = 0;
while FP < 0.05
    T = T - 0.01;
    A = R > T;
    FP = sum(A(trueMap == 0))/(N1*N2);
end;

    
