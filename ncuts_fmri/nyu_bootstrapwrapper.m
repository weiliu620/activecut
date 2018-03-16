function nyu_bootstrapwrapper(bs_start, bs_end, bs_sample_size)
addpath ~/projects/ncuts_fmri
addpath ~/projects/ncuts_fmri/Ncut_9

mask_file=['/usr/sci/scratch/weiliu/NYU_test_retest/ncuts_bootstrap/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii'];

subcls_dir=['/usr/sci/scratch/weiliu/NYU_test_retest/ncuts_bootstrap/ncuts_allmat'];
out_dir = ['/usr/sci/scratch/weiliu/NYU_test_retest/ncuts_bootstrap/out'];

[totalPts, linear2Sub] = GetTotalPts(mask_file);
allsubs = dir(subcls_dir);
allsubs = allsubs(~[allsubs.isdir]);
maskStruct = load_untouch_nii(mask_file);

for subid = 1:length(allsubs)
    load(fullfile(subcls_dir, allsubs(subid).name), 'M');
    allM{subid} = M;
    fprintf('Loading %d: %s done.\n', subid, allsubs(subid).name);
end

% change seed to generate different random number. 
for bsid = bs_start:bs_end
    fprintf('----------> Boostrap %d begin:\n', bsid);
    % Generate boostrap samples.
    bs_sample = randi(length(allsubs), bs_sample_size, 1);
    grpM = zeros(totalPts, totalPts);
    for i = 1:length(bs_sample)
        grpM = grpM + double(allM{bs_sample(i)});
        fprintf('    bsid = %d, adding sub %d: %s done.\n', bsid, i, allsubs(bs_sample(i)).name);
    end;
    
    % normalize M such that all elements are below 1.
    grpM = grpM / bs_sample_size;

    % Threshold M
    grpM(grpM < 0.3) = 0;
    [NcutDiscrete, ~, ~] = ncutW(grpM, 7);
    
    % Save group label map. 
    outfile = fullfile(out_dir, strcat('grp', num2str(bsid, '%03d'), '.nii') );
    grpLabelStruct = maskStruct;
    for n = 1:totalPts
        x = linear2Sub(n,1);
        y = linear2Sub(n,2);
        z = linear2Sub(n,3);
        grpLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
    end;
    save_untouch_nii(grpLabelStruct, outfile);
    fprintf('output label file saved to %s.\n', outfile);
end;
