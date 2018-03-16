function res = GibbsPriorSim(N, L, T)
% Usage:
% Use Gibbs Sampling to simulate a Gibbs distribution, i.e. Markov random
% field prior, on a 2D lattice. 

% Setup parameter
interval = 200; % interval for saving snapshot.
I = ceil(rand(N, N)*L);
close all; figure;
imshow(I/L);
saveas(gcf, strcat('Gibbs01', '.eps'), 'eps');
for k = 1:1000
    U = zeros(L,1); % Energy function.
    CP = zeros(L,1); % Conditional probability p(x | other data).
    for i = 2:N-1
        for j = 2:N-1
            for l = 1:L
                U(l) = sum([I(i-1,j), I(i+1,j), I(i,j-1), I(i,j+1)] ~= l);
                CP(l) = exp(-U(l)/T);
            end;
            CP = CP/sum(CP);
            F = CP(1);
            l = 1;
            U = rand;
            while(U > F)
                l = l+1;
                F = F + CP(l);
            end;
            I(i,j) = l;
        end;
    end;
    fprintf('Iteration: %d\n', k);
    if (mod(k, interval/10) == 0)
        imshow(I/L); pause(0.01);
        title(strcat('Iteration', int2str(k)));
    end;
    if (mod(k, interval) == 0)
        saveas(gcf, strcat('Gibbs', int2str(k), '.eps'), 'eps');
    end;
end;
            
            