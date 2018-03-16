function S = nyu_groupcuts(affdir, thre, maskFile, numClusters, outFile)
% Sum over all affinity matrix in files in 'affdir', and run Ncuts to get
% group segmentation. The aff. matrix is named 'M' in the mat files.

addpath ~/packages/nifti_to_matlab;
addpath ~/projects/ncuts_fmri/Ncut_9;


% compute linear2sub and totalPts.
[totalPts, linear2Sub] = GetTotalPts(maskFile);
fprintf('totalPts = %d\n', totalPts);
S = zeros(totalPts, totalPts);

% Read in aff. matrix from all mat files and sum over them.
numSubs = 0;
all_aff_files = dir(strcat(affdir, '/*.mat'));
for j = 1: length(all_aff_files)
    if ~all_aff_files(j).isdir
        numSubs = numSubs + 1; % number of subjects.
        load(fullfile(affdir, all_aff_files(j).name));
        S = S + double(M);
        fprintf('Reading aff. files from %s done. numSubs = %i\n', ...
                all_aff_files(j).name, numSubs);
    end;
end;

clear M;
pause(5);

S = S / numSubs;

% Threshold aff. matrix.
S(S < thre) = 0;

fprintf('Ncuts on group...\n');

maskStruct = load_untouch_nii(maskFile);
[NcutDiscrete, ~, ~] = ncutW(S, numClusters);

clear S;
% Save group label map.
grpLabelStruct = maskStruct;
[totalPts, ~] = size(linear2Sub);
for n = 1:totalPts
    x = linear2Sub(n,1);
    y = linear2Sub(n,2);
    z = linear2Sub(n,3);
    grpLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
end;

save_untouch_nii(grpLabelStruct, outFile);
fprintf('output label file saved to %s.\n', outFile);



