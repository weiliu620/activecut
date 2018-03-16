fmri_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
            'smoothed_nii/'];

mask_file=['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'mask/avg_4mm_thr0.3_bin_replicate.nii'];

out_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'ncuts4/'];

[totalPts, linear2Sub] = GetTotalPts(mask_file);

all_fmri = dir(fmri_dir);

% clean directory so all files are fmri files. 
all_fmri = all_fmri(~[all_fmri.isdir]);

% ran N-Cuts on all subjects.
for i = 1:length(all_fmri)
    [~, fmri_basename, ~] = fileparts(all_fmri(i).name);
    fmri_full = fullfile(fmri_dir, all_fmri(i).name);
    out_file = fullfile(out_dir, fmri_basename);
    M  =  segts(fmri_full, mask_file, 7, 0.4, out_file);
    save(strcat(out_dir, 'M', num2str(i, '%03d'), '.mat'), 'M', '-v7.3');
end;

