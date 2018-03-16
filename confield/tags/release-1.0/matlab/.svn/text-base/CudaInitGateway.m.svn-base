function C = CudaInitGateway(R, mu,sigma2, CMask)
% Usage: C = CudaInitGateway(R, mu,sigma2)
% CudaInit interface: 
% C = CudaSampling(R, mu(1), mu(2), sigma2(1), sigma2(2));
%
% Gateway function for initializing matrix C by maximum likelihood.

R = single(R);
mu = single(mu);
sigma2 = single(sigma2);
CMask = uint8(CMask);
C = GPUInit(R, mu(1), mu(2), sigma2(1), sigma2(2), CMask);