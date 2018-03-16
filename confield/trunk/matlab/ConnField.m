function res = ConnField(dataset, method, dir)
% Usage: ConnField(dataset, method, dir)
% main function for posterior connectivity.

s1 = 30; 
s2 = 30; % slice number to extract.
maskThresh = 0.5; % threshold to mask gray matter.

switch(dataset)
    case 'sin'
        if exist(strcat(dir, 'sinSample.mat'), 'file') == 2
            dataSaved = load(strcat(dir, 'sinSample.mat'));
            %X = dataSaved.X;
            RNS = dataSaved.RNS;
            RNS = RNS - diag(diag(RNS));
            N = size(RNS, 1);
            RNS = RNS + diag(max(RNS(:)) * ones(N,1));
            fprintf('load synthetic "sin" data file.\n');
        else
            error(ConnField: datafile, 'Data file not found.')
        end;
        if strcmp(method, 'cutoff')
            CutCorr(RNS, 0.12, dir);
            return;
        else
            EMPost(RNS, method, dir);
        end;
    case 'auditory'
        if exist(strcat(dir, 'auditoryData.mat'), 'file') == 2
            dataSaved = load(strcat(dir, 'auditoryData.mat'));
        else
        %[fmriVol, coord] = ExtractVol(maskThresh);
        %[im1, coord1, im2, coord2] = ExtractSlice(fmriVol, coord, s1, s2);
        %R = mycorr(im1', im2');
        %R = single(R);
        %EMPostReal(R, coord1, corod2, dir);
        end;
        [R, CMask, sz1, sz2] = GetCorr(maskThresh, s1, s2);
        C = EMPostReal(R, CMask, dir, sz1, sz2);
        
    case 'other'
end;


function CutCorr(R, T, dir)
R(R < T) = 0;
R(R >= T) = 1;
figure(1);
imshow(R, 'InitialMagnification', 'fit');
pause(1);
saveas(1, strcat(dir, 'RCutoff.eps'), 'epsc');
