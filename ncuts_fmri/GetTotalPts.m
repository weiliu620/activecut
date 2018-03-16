function [totalPts, linear2Sub] = GetTotalPts(maskFile)

addpath ~/packages/nifti_to_matlab;
maskStruct = load_untouch_nii(maskFile);

% Maximal number of points.
N = maskStruct.hdr.dime.dim(2) * maskStruct.hdr.dime.dim(3) ...
    * maskStruct.hdr.dime.dim(4);
linear2Sub = zeros(N,3);

totalPts = 0;
for x = 1:maskStruct.hdr.dime.dim(2)
    for y = 1:maskStruct.hdr.dime.dim(3)
        for z = 1:maskStruct.hdr.dime.dim(4)
            if maskStruct.img(x,y,z) > 0

                totalPts = totalPts + 1;
                linear2Sub(totalPts,:) = [x,y,z];

            end; % if
        end;
    end;
end;

% truncate at actual number of pts.
linear2Sub = linear2Sub(1:totalPts, :);