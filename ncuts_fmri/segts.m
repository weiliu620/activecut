function M = segts(srcFile, maskFile, numClusters, weightThresh, destFile)
% res = segts(srcFile, maskFile, numCluster, destFile)
%
% Segment fmri (either real or simulated) image by Normalized Cuts.

% add library so matlab can read nifti file.
addpath ~/packages/nifti_to_matlab;
addpath ~/projects/ncuts_fmri/Ncut_9/;

% load fmri file. Need to tell if it is .gz file.
[srcPath, tmp, srcFileExt] = fileparts(srcFile);
if strcmp(srcFileExt, '.gz')
    srcFile_nii = gunzip(srcFile);
    niiStruct = load_untouch_nii(srcFile_nii{1});
    delete(srcFile_nii{1});
else
    niiStruct = load_untouch_nii(srcFile);
end;

[maskPath, tmp, maskFileExt] = fileparts(maskFile);
if strcmp(maskFileExt, '.gz')
    maskFile_nii = gunzip(maskFile);
    maskStruct = load_untouch_nii(maskFile_nii{1});
    delete(maskFile_nii{1});
else
    maskStruct = load_untouch_nii(maskFile);
end;

% Rearrange nii data into a NxT. N is number of voxels within mask. T is
% time points. Init D with a N' larger than N to save time for allocatign
% memeory.

% Maximal number of points.
N = maskStruct.hdr.dime.dim(2) * maskStruct.hdr.dime.dim(3) ...
    * maskStruct.hdr.dime.dim(4);
T = niiStruct.hdr.dime.dim(5);

D = zeros(N, T);
linear2Sub = zeros(N,3);
sub2Linear = zeros(maskStruct.hdr.dime.dim(2), maskStruct.hdr.dime.dim(3), ...
    maskStruct.hdr.dime.dim(4));
totalPts = 0;
fprintf('begin filling data matrix D...\n');
startic = tic;
for x = 1:maskStruct.hdr.dime.dim(2)
    for y = 1:maskStruct.hdr.dime.dim(3)
        for z = 1:maskStruct.hdr.dime.dim(4)
            if maskStruct.img(x,y,z) > 0
                % save coordinate mapping data structure.
                totalPts = totalPts + 1;
                linear2Sub(totalPts,:) = [x,y,z];
                sub2Linear(x,y,z) = totalPts;
                D(totalPts,:) = reshape(niiStruct.img(x,y,z,:), ...
                        1, niiStruct.hdr.dime.dim(5));
            end; % if
        end;
    end;
end;
fprintf('Done with filling data matrix at %f seconds.\n', toc(startic));

D = D(1:totalPts,:);
fprintf('total number of gray matter data pts: %d\n', totalPts);

% Noramalized D so D(n,:) has zero mean and unit variance.
D = D - repmat(mean(D,2), 1, T);
D = D ./ repmat(std(D,0,2), 1, T);

% Compute weights metrix. For normalized vectors, correlation is just the
% dot product divided by T.
fprintf('Compute weight matrix. Patience...');

W = D * D' / (T-1);
W(W < weightThresh) = 0;
fprintf('Done with %f seconds.\n', toc(startic));

fprintf('running Ncuts...\n');
[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W, numClusters);
fprintf('Done with N-cuts at  %f seconds.\n', toc(startic));

% use mask image as initial image of output label map.
outLabelStruct = maskStruct;

for n = 1:totalPts
    x = linear2Sub(n,1);
    y = linear2Sub(n,2);
    z = linear2Sub(n,3);
    outLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
end;
fprintf('Done with outLabelStruct at  %f seconds.\n', toc(startic));

% Fill small holes if 4 neighbors have same labels.
% outLabelStruct = fill_hole(outLabelStruct, maskStruct, 4);
fprintf('Done with fill_hole at  %f seconds.\n', toc(startic));

save_untouch_nii(outLabelStruct, destFile);
fprintf('output label file saved to %s.\n', destFile);

% Convert label map into NxN matrix.
M = NcutDiscrete * NcutDiscrete';
fprintf('Done with affinity matrix M  at  %f seconds.\n', toc(startic));



% Fill small hole in label map by majority voting, just like Heuvel 2008
% paper did.
function outStr = fill_hole(inStr, maskStr, thresh)

for x = 1:maskStr.hdr.dime.dim(2)
    for y = 1:maskStr.hdr.dime.dim(3)
        for z = 1:maskStr.hdr.dime.dim(4)
            if maskStr.img(x,y,z) > 0            
                [new_label, freq] = mode([
                    inStr.img(x-1,y,z),
                    inStr.img(x+1,y,z),
                    inStr.img(x,y-1,z),
                    inStr.img(x,y+1,z),
                    inStr.img(x,y,z-1),                    
                    inStr.img(x,y,z+1)]);
                if (freq >= thresh) && (new_label ~= 0)
                    inStr.img(x,y,z) = new_label;
                end;
            end;
        end;
    end;
end;

outStr = inStr;
             