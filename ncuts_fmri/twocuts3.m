% This script run normalized cuts on synthetic data for miccai 2012.


fmri_dir = ['/home/sci/weiliu/dataset/syn2/data/obs_smoothed'];
mask_file=['/home/sci/weiliu/dataset/syn2/data/mask.nii.gz'];
out_dir = ['/home/sci/weiliu/dataset/syn2/test_ncuts/'];

all_fmri = dir(fmri_dir);

[totalPts, linear2Sub] = GetTotalPts(mask_file);
fprintf('totalPts = %d\n', totalPts);
M = zeros(totalPts, totalPts);

for i = 1:length(all_fmri)
    if ~all_fmri(i).isdir
        [~, fmri_basename, ~] = fileparts(all_fmri(i).name);
        fmri_full = fullfile(fmri_dir, all_fmri(i).name);
        out_file = fullfile(out_dir, fmri_basename);
        M = M + segts(fmri_full, mask_file, 5, 0.5, out_file);
    end;
end;
    
% normalize M such that all elements are below 1.
M = M / length(all_fmri);

% Threshold M
M(M < 0.4) = 0;

fprintf('Ncuts on group...\n');
GroupSeg(M, mask_file, linear2Sub, out_dir, 'outgrp.nii', 5);



