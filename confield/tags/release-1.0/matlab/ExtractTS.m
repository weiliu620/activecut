function [C] = ExtractTS(x, y, z, type)
% Usage: ExtractTS(x, y, z, type). x, y, z is voxel coordiante of the
% activated
% point, or the see region we're interested in.
% for SPM auditory data set, the most activated point is (62, -26, 10) -->
% (71, 44, 31).
% We take X(:, 44, 31, :) and X(:, 45, 31,:) as two one D image.

addpath('/home/sci/weiliu/packages/nifti_to_matlab'); % tool dir.
addpath('/scratch/data/auditory_spm_example'); % date dir
NII = load_nii('snrfM.nii'); % 4d nii file. already preprocessed.
xdim = NII.hdr.dime.dim(2);
ydim = NII.hdr.dime.dim(3);
zdim = NII.hdr.dime.dim(4); 
tdim = NII.hdr.dime.dim(5);
S1 = double(NII.img(:, y, z, :));
S1 = reshape(S1, xdim, tdim); % reshape to 2D.
S2 = double(NII.img(:, y+1, z, :));
S2 = reshape(S2, xdim, tdim);
S3 = double(NII.img(x, :, z, :));
S3 = reshape(S3, ydim, tdim);

S1 = S1'; S2 = S2'; S3 = S3'; % now they are tdim by xdim 
%C = corr(S1, S2);
C = mycorr(S1, S2);
save corr_fmri C;

% Visualize correlation.
figure(1); imshow((C+1)/2, 'InitialMagnification', 'fit');
title('Sample Correlation. C. Scaled to (0,1)');
figure(2);
hist(reshape(C, xdim*xdim, 1), xdim^2/100);

% try another 1d image.
%C3 = corr(S1, S3);
C3 = mycorr(S1, S3);
figure(4); imshow((C3+1)/2, 'InitialMagnification', 'fit');
title('Sample Correlation. C3. Scaled to (0,1)');
figure(5);
hist(reshape(C3, xdim*ydim, 1), xdim^2/100);
title('S1 and S2 correlation');

% visualize extracted region.
Q = double(NII.img(:,:,31, 1)); % extract one slice.
maxL = max(Q(:)); % max gray level intensity.
Q(:, y) = maxL * 1.5;
Q(:,y+1) = 0;
Q(x, :) = maxL * 1.5;
figure(3);
imshow(double(Q')/maxL, 'InitialMagnification', 'fit');



