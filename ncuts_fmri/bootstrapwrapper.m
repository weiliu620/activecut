function bootstrapwrapper(bs_start, bs_end, bs_sample_size)
% wrapper function to run boostrap sampling and then N-cuts on each BS
% sample.


fmri_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
            'smoothed_nii'];

mask_file=['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'mask/avg_4mm_thr0.3_bin_replicate.nii'];

out_dir = ['/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/' ...
           'ncuts4/'];

subcls_dir = ['/home/sci/weiliu/dataset/ncuts_subs/'];

[totalPts, linear2Sub] = GetTotalPts(mask_file);


all_fmri = dir(fmri_dir);

% clean directory so all files are fmri files. 
all_fmri = all_fmri(~[all_fmri.isdir]);

% load all subjects clustering .mat file once and for all.
for subIdx = 1: length(all_fmri)
    load(strcat(subcls_dir, 'M', num2str(subIdx, '%03d'), '.mat'));
    allM{subIdx} = M;
    fprintf('load subject %d done\n', subIdx);
end;    

for bsid = bs_start:bs_end
    fprintf('----------> Boostrap %d begin:\n', bsid);
    % Generate boostrap samples.
    bs_sample = randi(length(all_fmri), bs_sample_size, 1)
    grpM = zeros(totalPts, totalPts);
    for i = 1:length(bs_sample)
        grpM = grpM + allM{bs_sample(i)};
        fprintf('    bsid = %d, load sub %d done.\n', bsid, i);
    end;
    
    % normalize M such that all elements are below 1.
    grpM = grpM / bs_sample_size;

    % Threshold M
    grpM(grpM < 0.3) = 0;
    
    outfile = strcat('outgrp', num2str(bsid, '%03d'), '.nii');
    GroupSeg(grpM, mask_file, linear2Sub, out_dir, outfile, 7);
end;

        