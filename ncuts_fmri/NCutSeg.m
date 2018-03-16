function res = NCutSeg(srcFile, maskFile, numClusters, disKernel, destFile)
% res = NCutSeg(srcFile, numClusters)

% segment image by Normalized Cuts.
% input: srcFile = image to segment.
%        numClusters = number of clusters.

% add library so matlab can read nifti file.
addpath ~/packages/nifti_to_matlab;
addpath ~/projects/ncuts_fmri/Ncut_9/;

% use load_untouch_nii instead of load_nii since I don't want to do any
% transform on the data.
niiStruct = load_untouch_nii(srcFile);
maskStruct = load_untouch_nii(maskFile);

% Rearrange nii data into a cell array. Each cell might be a scalar
% intensity value, or a time series, depending on what data we work on.
D = cell(maskStruct.hdr.dime.dim(2) * maskStruct.hdr.dime.dim(3) ...
    * maskStruct.hdr.dime.dim(4),1);
totalPts = 0;
linear2Sub = zeros(maskStruct.hdr.dime.dim(2) * maskStruct.hdr.dime.dim(3)...
    * maskStruct.hdr.dime.dim(4), 3, 'uint16');
sub2Linear = zeros(maskStruct.hdr.dime.dim(2:4), 'uint32');

for x = 1:maskStruct.hdr.dime.dim(2)
    for y = 1:maskStruct.hdr.dime.dim(3)
        for z = 1:maskStruct.hdr.dime.dim(4)
            if maskStruct.img(x,y,z) > 0
                % save coordinate mapping data structure.
                totalPts = totalPts + 1;
                linear2Sub(totalPts,:) = [x,y,z];
                sub2Linear(x,y,z) = totalPts;
                if niiStruct.hdr.dime.dim(1) <= 3
                    % We are using scalar data.
                    D{totalPts} = niiStruct.img(x,y,z);
                elseif niiStruct.hdr.dime.dim(1) == 4
                    % We are using fmri time series.
                    D{totalPts} = reshape(niiStruct.img(x,y,z,:), ...
                        niiStruct.hdr.dime.dim(5),1);
                end;
            end;
        end;
    end;
end;
D = D(1:totalPts);

% Define kernel function for computing graph weights.
if isequal(disKernel, 'recsqr')
    kernelFun = @(x,a) exp( - (x - a)^2);
elseif isequal(disKernel, 'corr')
    kernelFun = @(x,a) dot(x,a) / length(x);
end;


% Init a Weight matrix.
W = zeros(totalPts, totalPts);
% Each iteration, compute weights between D{i} and each of D{i+1}...D{end}.
% for i = 1:totalPts
%     W(i,i+1:end) = cellfun(@(x) kernelFun(x, D{i}), D(i+1:end) );
%     fprintf('kernelFun done for data point %d.\n', i);
% end;
tic;
for i = 1:totalPts
    fprintf('working on idx %d.\n', i);
    for j = (i+1):totalPts
        W(i,j) = exp(-(D{i} - D{j})^2);
    end;
end;
fprintf('workin on W takes %f seconds.\n', toc);
        

% fill lower triangular matrix and diagonal part.
W = W + W' + diag(ones(totalPts, 1));

% Run Ncuts.
fprintf('running Ncuts...\n');
[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(W,numClusters);

% Write labels to image.
fprintf('write labels to image.\n');
% use mask image as initial image of output label map.
outLabelStruct = maskStruct;

for n = 1:totalPts
    x = linear2Sub(n,1);
    y = linear2Sub(n,2);
    z = linear2Sub(n,3);
    outLabelStruct.img(x,y,z) = find(NcutDiscrete(n,:));
end;

save_untouch_nii(outLabelStruct, destFile);

