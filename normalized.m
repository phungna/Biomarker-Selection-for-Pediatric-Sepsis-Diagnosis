function n_data = normalized(gene_data)
n_data = (gene_data - min(gene_data))/(max(gene_data) - min(gene_data))
end