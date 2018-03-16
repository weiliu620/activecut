function GroupSeg(W, maskFile, linear2Sub, outFile, K)

[maskPath, tmp, maskFileExt] = fileparts(maskFile);

maskStruct = load_untouch_nii(maskFile);
[NcutDiscrete, ~, ~] = ncutW(W, K);

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