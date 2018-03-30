function relmat = make_probe_gene_relmat(gene_file)
% return a matrix with (i,j) set to 1 if i and j test same genes.

gene_list = dlmread(gene_file);
gene_list = int32(gene_list);

[n_genes, ~] = size(gene_list);
relmat = zeros(n_genes, 'uint8');
for gidx = 1:n_genes
    relmat(gidx, gidx:end) = (gene_list(gidx) == gene_list(gidx:end) );
    if (gidx % 100 == 0)
        fprintf('gidx: %i\n');
    end;
end;
