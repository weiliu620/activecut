
D = csvread(gene_samples);
D = D(:, 2:end);
[numRep, numVar] = size(D);
D = D - repmat(mean(D, 1), numRep, 1);
C = D' * D;

    
    

