function [R, CMask, sz1, sz2] = GetCorr(T, s1, s2)
% Obtain whole volume data, without masking.
% input: T, a scalar in (0, 1). For masking the volume.
% output: fmriVol, 3D matrix. Single precision.
%         maskVol: 3D mask matrix. uint8 data type.

rootDir = '/scratch/data/auditory_spm_example/';
greyMaskFile = strcat(rootDir, 'sM00223_seg_coreg/rc1sM00223_002.hdr');
fmriFile = strcat(rootDir,...
    'fM00223_realign_reslice/auditory_realigned.nii');


addpath('/home/sci/weiliu/packages/nifti_to_matlab'); % tool dir.
addpath('/scratch/data/auditory_spm_example'); % date dir
NII = load_nii(fmriFile); % 4d nii file. already preprocessed.
maskStruct = load_nii(greyMaskFile);

fmriVol = single(NII.img);
maskVol = single(maskStruct.img);
maskVol(maskVol > T) = 1;
maskVol(maskVol <= T) = 0;

% compute correlation R matrix, and mask matrix.
[xd, yd, zd, td] = size(fmriVol);
im1 = fmriVol(:,:,s1, :);
im1 = reshape(im1, xd*yd, td);
im2 = fmriVol(:,:,s2,:);
im2 = reshape(im2, xd*yd, td);
%R = single(mycorr(im1', im2'));
%R = single(corr(im1', im2'));
R = GPUCorrGateway(im1', im2');
sz1 = [xd; yd];
sz2 = [xd; yd];

% mask.
m1 = maskVol(:,:,s1);
m1 = reshape(m1, xd*yd, 1);
m2 = maskVol(:,:,s2);
m2 = reshape(m2, xd*yd, 1);
CMask = m1 * m2';
CMask = uint8(CMask);



