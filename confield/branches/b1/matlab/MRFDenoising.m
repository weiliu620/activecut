function MRFDenoising(savedState)
%Usage: MRFDenoising(savedState)
eta = 1;
I = imread('happy4.tif', 'tiff');
I = im2double(I);

[maxRow, maxCol] = size(I);
defaultStream.State = savedState;
N = double(rand(maxRow, maxCol) > 0.9); % noise
IN = xor(I,N); % multiply noise (not additie)
subplot(2,3,1); imshow(I);
title('Original image.');
subplot(2,3,2); imshow(IN);
title('add 10% noise');

Y = IN * 2 - 1; % Scale to [-1, 1].
XInit = (rand(maxRow, maxCol) > 0.5) * 2 -1; %init X


beta = 0 % add prior weights.
X = XInit;
subplot(2,3,3);
for k = 1 : 5
    allpos = randperm(maxRow*maxCol);
    for m = 1:maxRow*maxCol
        pos = allpos(m);
        [i j] = ind2sub([maxRow, maxCol], pos);
        % Select a pixel at random
        neighborhood = pos + [-1 1 -maxRow maxRow]; % Find indicies of neighbours
        neighborhood( find( [i == 1, i == maxRow, j == 1, j == maxCol] ) ) = [];
        % pesky boundaries...thrash them
        nagree = sum( X(i,j) == X(neighborhood) );
        ndisagree = sum( X(i,j) ~= X(neighborhood) );
        change = nagree - ndisagree;
        deltaE = beta*(nagree - ndisagree) + eta * (Y(i,j)/X(i,j));
        if deltaE < 0
            X(i,j) = - X(i,j);
        end;
    end;
    imshow((X+1)/2); pause(0.1);
    title('no prior. beta == 0');
    fprintf('beta = %d, iteration %d\n', beta, k);
end;


beta = 0.3
X = XInit;
subplot(2,3,4);
for k = 1 : 5
    allpos = randperm(maxRow*maxCol);
    for m = 1:maxRow*maxCol
        pos = allpos(m);
        [i j] = ind2sub([maxRow, maxCol], pos);
        % Select a pixel at random
        neighborhood = pos + [-1 1 -maxRow maxRow]; % Find indicies of neighbours
        neighborhood( find( [i == 1, i == maxRow, j == 1, j == maxCol] ) ) = [];
        % pesky boundaries...thrash them
        nagree = sum( X(i,j) == X(neighborhood) );
        ndisagree = sum( X(i,j) ~= X(neighborhood) );
        change = nagree - ndisagree;
        deltaE = beta*(nagree - ndisagree) + eta * (Y(i,j)/X(i,j));
        if deltaE < 0
            X(i,j) = - X(i,j);
        end;
    end;
    imshow((X+1)/2);pause(0.1);
    title('beta = 0.3')
    fprintf('beta = %d, iteration %d\n', beta, k);
end;

beta = 0.6 % add prior weights.
X = XInit;
subplot(2,3,5);
for k = 1 : 10
    allpos = randperm(maxRow*maxCol);
    for m = 1:maxRow*maxCol
        pos = allpos(m);
        [i j] = ind2sub([maxRow, maxCol], pos);
        % Select a pixel at random
        neighborhood = pos + [-1 1 -maxRow maxRow]; % Find indicies of neighbours
        neighborhood( find( [i == 1, i == maxRow, j == 1, j == maxCol] ) ) = [];
        % pesky boundaries...thrash them
        nagree = sum( X(i,j) == X(neighborhood) );
        ndisagree = sum( X(i,j) ~= X(neighborhood) );
        change = nagree - ndisagree;
        deltaE = beta*(nagree - ndisagree) + eta * (Y(i,j)/X(i,j));
        if deltaE < 0
            X(i,j) = - X(i,j);
        end;
    end;
    imshow((X+1)/2); pause(0.1);
    title('beta = 0.6');
    fprintf('beta = %d, iteration %d\n', beta, k);
end;






    
