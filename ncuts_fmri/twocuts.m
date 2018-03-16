% N-cuts script for single run. for figure 2. 

fmri_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
            'smoothed/'];

mask_file=['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'mask/avg_4mm_thr0.3_bin.nii.gz'];

out_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'ncuts/'];


all_fmri = dir(fmri_dir);

[totalPts, linear2Sub] = GetTotalPts(mask_file);
fprintf('totalPts = %d\n', totalPts);
M = zeros(totalPts, totalPts);

for i = 1:length(all_fmri)
    if ~all_fmri(i).isdir
        [~, fmri_basename, ~] = fileparts(all_fmri(i).name);
        fmri_full = fullfile(fmri_dir, all_fmri(i).name);
        out_file = fullfile(out_dir, fmri_basename);
        M = M + segts(fmri_full, mask_file, 7, 0.4, out_file);
    end;
end;
    
% normalize M such that all elements are below 1.
M = M / length(all_fmri);

% Threshold M
M(M < 0.3) = 0;

fprintf('Ncuts on group...\n');
GroupSeg(M, mask_file, linear2Sub, out_dir, 'outgrp.nii', 7);



