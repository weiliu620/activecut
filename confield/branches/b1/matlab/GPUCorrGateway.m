function R = GPUCorrGateway(im1, im2)
% Usage: R = GPUCorrWateway(im1, im2). Gateway function for computing
% correlation of pairs elements beween columns of im1 and im2.
% im1 and im2 will be converted to single in this function.

im1 = single(im1);
im2 = single(im2);
[D1, M] = size(im1);
[D2, N] = size(im2);
if D1 ~= D2
    error('GPUCorrGateway: dimension', 'rows of im1 and im2 must be equal');
end;

% call mex function for heavy job.
R = GPUCorr(im1, im2);

