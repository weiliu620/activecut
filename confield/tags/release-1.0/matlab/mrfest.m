function res =  mrfest(alpha, beta)
N = 100;
I = double(rand(N,N) > 0.5);
iterCount = 1;
figure(1);
imshow(I); title(strcat('iteration:', int2str(iterCount)));
if isoctave()
  print -depsc mrfest01.eps;
else
  saveas(gcf, strcat('figures/mrfest', int2str(iterCount), '.eps'), 'eps');
end;


[maxRow, maxCol] = size(I);
p = zeros(2,1);
while(iterCount < 100)
    allpos = randperm(maxRow*maxCol);
    for n = 1:maxRow*maxCol
        pos = allpos(n);
        [i,j] = ind2sub([maxRow, maxCol], pos);
        neighbors = pos + [-1 1 -maxRow maxRow];
        neighbors([i == 1, i == maxRow, j == 1, j == maxCol]) = [];
        U0 = sum(I(neighbors) ~= 0) - sum(I(neighbors) == 0);
        U1 = -U0;
        p(1) = exp(-beta*U0);
        p(2) = exp(-alpha - beta*U1);
        p = p/sum(p); %normalize.
        if rand > p(1)
            I(i,j) = 1;
        else
            I(i,j) = 0;
        end;
    end;
    figure(1);
    if isoctave()
        imshow(I);
    else
        imshow(I, 'InitialMagnification', 'fit');
    end;
    title(strcat('iteration:', int2str(iterCount)));
    iterCount = iterCount +1;
    pause(0.01);
end;


% Now compute the log-likelihood given beta. suppose alpha is known.
allBeta = (beta-2:0.1:beta+2);
LL = zeros(size(allBeta));
figure(2);
for k = 1:length(allBeta)
    beta = allBeta(k);
    for n = 1:maxRow*maxCol
        [i,j] = ind2sub([maxRow, maxCol], n);
        neighbors = n + [-1 1 -maxRow maxRow];
        neighbors([i == 1, i == maxRow, j == 1, j == maxCol]) = [];
        U0 = sum(I(neighbors) ~= 0) - sum(I(neighbors) == 0);
        U1 = -U0;
        if (I(i,j) == 1)
            U = U1;
        else
            U = U0;
        end;
        Z = exp(-alpha*0 - beta*U0) + exp(-alpha*1 - beta * U1);
        LL(k) = LL(k) + (-alpha * I(i,j) - beta * U - log(Z));
    end;
    fprintf('log-likelihood for beta: %d is: %d\n', allBeta(k), LL(k));
    plot(allBeta(1:k), LL(1:k), 'bo-'); pause(0.01);
end;

plot(allBeta, LL, 'bo-');
title('Log-likelihood over beta');
llfile = 'figures/mrfest02.eps';
if isoctave()
    print -depsc llfile;
else
    saveas(2, llfile, 'epsc');
end;

    
        
        
