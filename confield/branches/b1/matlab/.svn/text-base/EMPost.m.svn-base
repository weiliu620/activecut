function res = EMPost(R, method, dir)
% Given sample correlation matrix R,
% compute the posterior matrix C, where c_ij = p(c_ij|r_ij).
% EMethod:
% 1) 'sampling': 
% 2) 'Annealing': change T during iteration of sampling.
% 3) 'smf': first sampling, then iteratively compute mean field.
% 4) 'amf': first annealing, then iteratively compute mean field.
        
% init model parameters.
mu = [0 0.3]; % Bold guess. mu(1) is mean for noconnectivity.
sigma2 = [0.1 0.1]; % when A^2/sigma^2 is 10. Roughly guess.

% initialize
[N, M] = size(R); % number of voxels.
C = zeros(N, M); % Connectvity matri. Cij = {-1, 1}
coor = randperm(N*M);
gc = 1./(sqrt(2*pi*sigma2)); % Gaussian constant.
plotInterval = 5;
for n = coor
    [i,j] = ind2sub([N N], n);
    lhood(2) = gc(2) * exp(-(R(i,j)-mu(2))^2/(2*sigma2(2)));
    lhood(1) = gc(1) * exp(-(R(i,j)-mu(1))^2/(2*sigma2(1)));
    if lhood(2) > lhood(1)
        C(i,j) = 1;
    else
        C(i,j) = -1;
    end;
end;

T = 5; % Annealing temperature.
beta = 1/T;
%maxIterCount = 20;
%maxSAIterCount = 40;
iterCount = 1;
para = [mu, sigma2]; % all parameters.
allPost = zeros(N, M, 2);
allBeta = beta;
change = inf;
maxMFIter = 5;
maxSAIterCount = 10;
Nk = 0.5 * N * M * ones(2,1);
while(change > 10^(-3))
    % E step. Compute posterior probability.
    SAIterCount = 2;
    switch lower(method)
        case {'mrf'}
            while(SAIterCount <= maxSAIterCount)
                
                C = ESampling(C, R, mu, sigma2, beta);
                
                % show connectivity map.
                figure(1); imshow((C+1)*0.5, 'InitialMagnification', 'fit');
                title(strcat('Posterior Connectivity: Sampled'));
                
                fprintf('EM iteration: %d. SA iteration: %d\n'...
                    , iterCount, SAIterCount);
                SAIterCount = SAIterCount + 1;
                if iterCount == 1
                    imwrite(uint8(255*0.5*(C+1)), strcat(dir,'postConn.gif'), 'gif', 'LoopCount', 100);
                else
                    imwrite(uint8(255*0.5*(C+1)), strcat(dir,'postConn.gif'), 'gif',...
                        'WriteMode', 'append');
                end;
                
                C = EMeanField(C, R, mu, sigma2, beta, maxMFIter);
                % show connectivity map.
                figure(1); imshow((C+1)*0.5, 'InitialMagnification', 'fit');
                title(strcat('Posterior Connectivity: Expected value'));
                
                % Now the field is stable (equililbrium), so  we save posterior for
                % estimating parameters.
                allPost = ELocalProb(C, R, mu, sigma2, beta);
                figure(12); imshow(allPost(:,:,2), 'InitialMagnification', 'fit');
                %title('Posterior map for C = 1 (have connectivity)');
                pause(0.01);
            end;
        case{'ml'}
            allPost = MLEstep(R, mu, sigma2, Nk(1), Nk(2));
            figure(12); imshow(allPost(:,:,2) > 0.5);
            %title('Posterior map for C = 1 (have connectivity)');
    end
    
    % M step. estimate parameters.
    % # of c = 0 in triangular matrix.
    Nk(1) = sum(sum(allPost(:,:,1)));
    Nk(2) = sum(sum(allPost(:,:,2))); % c = 1
    % mu(1) = 0. No change. 
    CR = allPost(:,:,2) .* R;
    CR = CR - diag(diag(CR));
    mu(2) = sum(CR(:))/Nk(2);
    
    RM1 = allPost(:,:,1) .* (R - mu(1)).^2; 
    RM2 = allPost(:,:,2) .* (R - mu(2)).^2; 
    sigma2(1) = sum(RM1(:))/Nk(1);
    sigma2(2) = sum(RM2(:))/Nk(2); 
    if strcmp(method, 'mrf')
        % Estimate beta. This need a lot of work. We begin with previous beta
        % as a init value for Newton's method.
        epsilon = 10^(-6);
        maxBeta =10;
        step = inf;
        newBeta = beta;
        while (step > epsilon && newBeta < maxBeta)
            J = 0; DELL =0 ;
            for n = 1:N*M
                [i,j] = ind2sub([N M], n);
                neighbors = n + [-1 1 -N N];
                neighbors([i==1, i==N, j==1, j==M]) = [];
                nsum = sum(C(neighbors));
                J = J + nsum^2 * (tanh(beta*nsum)^2 - 1); %2nd deriv
                DELL = DELL + nsum*(C(i,j) - tanh(beta*nsum));
            end;
            step = -DELL/J;
            newBeta = beta + step;
            if newBeta > maxBeta || isnan(newBeta)
                %fprintf('beta is already very large. Unchanged.\n');
                beta = maxBeta;
                fprintf('beta set to 100\n');
            else
                beta = newBeta;
                fprintf('step = %d, new beta = %2.2f\n', step, beta);
            end;
            
            allBeta = [allBeta; beta];
            figure(2);
            plot(allBeta, 'ro-');
            title('Estimation of beta');
        end;
    end;
    % show the histogram.
    sel = logical(allPost(:,:,2) < 0.5);
    R0 = R(sel);
    sel = logical(allPost(:,:,2) > 0.5);
    R1 = R(sel);

    figure(3);
    hist(R0, Nk(1)/200);
    title('sample correlation given c = 0');
    figure (4); hist(R1, Nk(2)/20);
    title('sample correlation given c = 1');
    
    % save parameters and dispaly the trend.
    if mod(iterCount, plotInterval) == 0
        para = [para; [mu, sigma2]];
        change = sum(abs(para(end,:) - para(end-1,:))); % if change big, continue.
        Xcoord = 0:(plotInterval):iterCount;
        figure(5);
        plot(Xcoord, para(:,2), 'ro-',...
            Xcoord, para(:,3), 'b*--',...
            Xcoord, para(:,4),  'r*--');
        title('trend of parameters in EM iteration.');
        xlabel('EM iteration');
        legend('mu1', 'sgima^2 of c=0', 'sigma^2 of c=1');
    end;
    
    fprintf('iteration: %d\n', iterCount);
    iterCount = iterCount +1;        
end;


mysave(12, strcat(dir, 'postP1.eps'));
if strcmp(method, 'mrf')
    mysave(1, strcat(dir, 'postConnMap.eps'));
    mysave(2, strcat(dir, 'beta.eps'));
end;
mysave(3, strcat(dir, 'dond_corr0.eps'));
mysave(4, strcat(dir, 'dond_corr1.eps'));
mysave(5, strcat(dir, 'para_est.eps'));

function C = ESampling(C, R, mu, sigma2, beta)

[N, M] = size(C);
gc = 1./(sqrt(2*pi*sigma2)); % Gaussian constant.
coor = randperm(N*M);
for n = coor
    [i,j] = ind2sub([N M], n);
    neighbors = n + [-1 1 -N N];
    neighbors([i==1, i==N, j==1, j==M]) = [];
    nsum = sum(C(neighbors));
    denom = 2*cosh(beta*nsum);
    pr(1) = exp(beta*(-1)*nsum)/denom;
    pr(2) = exp(beta*(+1)*nsum)/denom;
    % Likelihood.
    lhood(2) = gc(2) * exp(-(R(i,j)-mu(2))^2/(2*sigma2(2)));
    lhood(1) = gc(1) * exp(-(R(i,j)-mu(1))^2/(2*sigma2(1)));
    % Posterior.
    post(1) = pr(1) * lhood(1);
    post(2) = pr(2) * lhood(2);
    % Normalize so they sum to one.
    post = post/sum(post);
    % Gibbs Sampling.
    if rand < post(1)
        C(i,j) = -1;
    else
        C(i,j) = 1;
    end;    
end;

function allPost = MLEstep(R, mu, sigma2, Nk1, Nk2)

[N, M] = size(R);
gc = 1./(sqrt(2*pi*sigma2)); % Gaussian constant.
allPost = zeros(N,M,2);
for n = 1:N*M
    [i,j] = ind2sub([N M], n);
    pr(1) = Nk1;
    pr(2) = Nk2;
    pr = pr/sum(pr);
    % Likelihood.
    lhood(2) = gc(2) * exp(-(R(i,j)-mu(2))^2/(2*sigma2(2)));
    lhood(1) = gc(1) * exp(-(R(i,j)-mu(1))^2/(2*sigma2(1)));
    % Posterior.
    allPost(i,j,1) = pr(1) * lhood(1);
    allPost(i,j,2) = pr(2) * lhood(2);
    % Normalize so they sum to one.
    allPost(i,j,:) = allPost(i,j,:)/sum(allPost(i,j,:));
end;

function allPost = ELocalProb(C, R, mu, sigma2, beta)
% compute local posterior probability p(x | nbor(x))
[N, M] = size(C);
gc = 1./(sqrt(2*pi*sigma2)); % Gaussian constant.
allPost = zeros(N,M,2);
for n = 1:N*M
    [i,j] = ind2sub([N M], n);
    
    neighbors = n + [-1 1 -N N];
    neighbors([i==1, i==N|i==j-1, j==1|i==j-1, j==N]) = [];
    nsum = sum(C(neighbors));
    denom = 2*cosh(beta*nsum);
    pr(1) = exp(beta*(-1)*nsum)/denom;
    pr(2) = exp(beta*(+1)*nsum)/denom;
    % Likelihood.
    lhood(2) = gc(2) * exp(-(R(i,j)-mu(2))^2/(2*sigma2(2)));
    lhood(1) = gc(1) * exp(-(R(i,j)-mu(1))^2/(2*sigma2(1)));
    % Posterior.
    allPost(i,j, 1) = pr(1) * lhood(1);
    allPost(i,j, 2) = pr(2) * lhood(2);
    % Normalize so they sum to one.
    allPost(i,j,:) = allPost(i,j,:)/sum(allPost(i,j,:));    
end;

function C = EMeanField(C, R, mu, sigma2, beta, maxMFIter)
% mean field 
[N, M] = size(C);
gc = 1./(sqrt(2*pi*sigma2)); % Gaussian constant.
for k = 1:maxMFIter
    coor = randperm(N*M);
    for n = coor
        [i,j] = ind2sub([N M], n);
            neighbors = n + [-1 1 -N N];
            neighbors([i==1, i==N, j==1, j==M]) = [];
            nsum = sum(C(neighbors));
            denom = 2*cosh(beta*nsum);
            pr(1) = exp(beta*(-1)*nsum)/denom;
            pr(2) = exp(beta*(+1)*nsum)/denom;
            % Likelihood.
            lhood(2) = gc(2) * exp(-(R(i,j)-mu(2))^2/(2*sigma2(2)));
            lhood(1) = gc(1) * exp(-(R(i,j)-mu(1))^2/(2*sigma2(1)));
            % Posterior.
            post(1) = pr(1) * lhood(1);
            post(2) = pr(2) * lhood(2);
            % Normalize so they sum to one.
            post = post/sum(post);
            % we dont' sample from posterior. Instead we compute the
            % expected value. Mean field theory.
            C(i,j) = post(2) - post(1) ;
    end;
end;



