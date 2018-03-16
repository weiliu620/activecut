function M = two_cuts(fmriDir, maskFile, numClusters, weightThresh, destDir, bsfix)

% Segment both subject and group level in this function. This fun is a merge
% of nyu_segts and nyu_groupcuts
%
% M = two_cuts(fmriDir, maskFile, numClusters, weightThresh, destDir, bsfix)



addpath ~/packages/nifti_to_matlab;
addpath ~/projects/ncuts_fmri/Ncut_9/;

% Load mask file.
maskStruct = load_untouch_nii(maskFile);

totalPts = 0;
linear2Sub = zeros(totalPts,3, 'uint8');
sub2Linear = zeros(maskStruct.hdr.dime.dim(2), maskStruct.hdr.dime.dim(3), ...
                   maskStruct.hdr.dime.dim(4), 'uint32');
for x = 1:maskStruct.hdr.dime.dim(2)
    for y = 1:maskStruct.hdr.dime.dim(3)
        for z = 1:maskStruct.hdr.dime.dim(4)
            if maskStruct.img(x,y,z) > 0
                % save coordinate mapping data structure.
                totalPts = totalPts + 1;
                linear2Sub(totalPts,:) = [x,y,z];
                sub2Linear(x,y,z) = totalPts;
            end; % if
        end;
    end;
end;

all_aff_files = dir(strcat(fmriDir, '/*.nii.gz'));

S = zeros(totalPts, totalPts);
numSubs = 0;
for j = 1: length(all_aff_files)
    if ~all_aff_files(j).isdir
        numSubs = numSubs + 1; % number of subjects.
        [~, sub_name, ~] = fileparts(all_aff_files(j).name); % remove .gz
        [~, sub_name, ~] = fileparts(sub_name); % remove .nii
        destFile = fullfile(destDir, strcat(sub_name, bsfix, '.nii.gz') );
        M = sub_seg(fullfile(fmriDir, all_aff_files(j).name), linear2Sub, sub2Linear, ...
                    maskStruct, destFile, weightThresh, numClusters);
        S = S + M;
    end;
end;

S = S / numSubs;
S(S < weightThresh) = 0;
fprintf('Ncuts on group...\n');
[NcutDiscrete, ~, ~] = ncutW(S, numClusters);

% Save group label map.
grpLabelStruct = maskStruct;
for n = 1:totalPts
    x = linear2Sub(n,1);
    y = linear2Sub(n,2);
    z = linear2Sub(n,3);
    grpLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
end;

grpFile = fullfile(destDir, strcat('grp', bsfix, 'nii.gz'));
save_untouch_nii(grpLabelStruct, outFile);


        

function M = sub_seg(subjectFile, linear2Sub, sub2Linear, maskStruct, destFile, ...
                     weightThresh, numClusters)
% Rearrange nii data into a NxT. N is number of voxels within mask. T is
% time points. Init D with a N' larger than N to save time for allocatign
% memeory.

niiStruct = load_untouch_nii(subjectFile);
T = niiStruct.hdr.dime.dim(5);
[totalPts, ~] = size(linear2Sub);
fprintf('sub_seg(): Filling data matrix D...\n');
D = zeros(totalPts, T);
for n = 1:totalPts
    D(n,:) = reshape(niiStruct.img(linear2Sub(n,1), linear2Sub(n,2), ...
                                   linear2Sub(n,3),:), 1, T);
end;

% Noramalized D so D(n,:) has zero mean and unit variance.
D = D - repmat(mean(D,2), 1, T);
D = D ./ repmat(std(D,0,2), 1, T);

% Compute weights metrix. For normalized vectors, correlation is just the
% dot product divided by T.
fprintf('Compute weight matrix. Patience...');

W = D * D' / (T-1);
W(W < weightThresh) = 0;

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W, numClusters);

outLabelStruct = maskStruct;

for n = 1:totalPts
    x = linear2Sub(n,1);
    y = linear2Sub(n,2);
    z = linear2Sub(n,3);
    outLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
end;

save_untouch_nii(outLabelStruct, destFile);
fprintf('output label file saved to %s.\n', destFile);

NcutDiscrete = single(NcutDiscrete);
M = NcutDiscrete * NcutDiscrete';
