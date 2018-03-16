function [R, Cmask] = SliceAndMask(fmriVol, maskVol, s1, s2)
% Usage: ExtractSlice()
% compute mask image (same dim with C), and also correlation matrix.

[xd, yd, zd, td] = size(fmriVol);
im1 = fmriVol(:,:,s1, :);
im1 = reshape(im1, xd*yd, td);
im2 = fmriVol(:,:,s2,:);
im2 = reshape(im2, xd*yd, td);
[xdim, ydim, zdim] = size(fmriVol);
% image 1
coordFilter = coord(:,3 == s1);
coord1 = coord(coordFilter,:);
im1 = fmriVol(coordFilter,:);
% image 2
coordFilter = coord(:,3 == s2);
coord2 = coord(coordFilter,:);
im2 = fmriVol(coordFilter,:);