function res = segfmri(fmri_file, mask_file, label_file, n_dims, n_cls, method, varargin)
% dimension reduction on fmri BILD singal and do segmentation.

addpath ~/packages/nifti_to_matlab;
addpath ~/packages/drtoolbox;
addpath ~/packages/drtoolbox/techniques;
addpath ~/projects/old/GraphDemos;
% First normalize the data.
nifti_strt = load_untouch_nii(fmri_file);
mask_strt = load_untouch_nii(mask_file);

ts_length = nifti_strt.hdr.dime.dim(5);
vol_size = nifti_strt.hdr.dime.dim(2:4);
bold = reshape(nifti_strt.img, vol_size(1)*vol_size(2)*vol_size(3), ts_length);
mask = reshape(mask_strt.img, vol_size(1)*vol_size(2)*vol_size(3), 1);

% extract gray matter bold signal.
bold_gray = bold(mask==1, :);
n_pts = size(bold_gray, 1);

% estimate the dimension
% n_dims_est = round(intrinsic_dim(bold_gray, 'MLE'));
% fprintf('MLE estimate of intrinsic dimensions: %i.\n', n_dims_est);

% Dimensional reduction.
if isempty(varargin)
    bold_dr = compute_mapping(bold_gray, method, n_dims);
else
    bold_dr = compute_mapping(bold_gray, method, n_dims, varargin{1});
end;


% Segmentation.
try
    [clustering, centers] = KmeansPiotrDollar(bold_dr, n_cls, 'replicates', 5);
    
catch err
  warning('kmeans on the eigenvectors did not work. Spectral clustering NaN. '); 
  clustering = NaN * zeros(1,size(n_pts,1)); 
  centers = NaN;
  rethrow(err);
end

% save labels.
labels_strt = mask_strt;
labels = reshape(labels_strt.img, vol_size(1)*vol_size(2)*vol_size(3), 1);
labels(:) = 0;
labels(mask==1) = clustering;
labels_strt.img = reshape(labels, vol_size(1), vol_size(2), vol_size(3));
save_untouch_nii(labels_strt, label_file);


