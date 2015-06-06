function [geneinfo] = extract_gene_info(gene_table, ids)
    % Extract additional info about a gene id from a platform gene table.
    [~,types] = size(gene_table);
    n_id = numel(ids);
    geneinfo = cell(n_id, types);
    table_ids = num2str(cell2mat(gene_table(:, 1)));
    for i = 1:numel(ids)
        geneinfo(i, :) = gene_table(strcmp(table_ids, ids(i)), :);
    end
end
