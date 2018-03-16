function C = CudaSamplingGateway(CInit, R, CMask,...
    mu, sigma2, beta, maxSAIterCount,sz1, sz2)
% Usage: (C, R, mu, sigma2, beta, maxSAIterCount).
% CudaSampling interface: 
% C = CudaSampling(C, R, CMask, mu(1), mu(2), sigma2(1), sigma2(2),...
%    beta, maxSAIterCount, CRow, CCol);
%
% Gateway function for sampling in gpu.

CInit = single(CInit);
R = single(R);
CMask = uint8(CMask);
mu = single(mu);
sigma2 = single(sigma2);
beta = single(beta);
maxSAIterCount = single(maxSAIterCount);

C = CudaSampling(CInit, R, CMask, mu(1), mu(2), sigma2(1), sigma2(2),...
    beta, maxSAIterCount, sz1(1), sz1(2), sz2(1), sz2(2));
