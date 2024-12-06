##selecionar os top50 genes mais diferencialmente expressos no dataframe de resultados pelo padj

top50_genes <- dds_res[order(dds_res$padj),][1:50,]

##selecionar os nomes desses genes

top50_genes.names <- dds_res[order(dds_res$padj),][1:50,]
rownames(top50_genes.names) <- sub("\\..*", "", rownames(top50_genes.names))
top50_genes.names$symbol <- mapIds(org.Hs.eg.db, keys = rownames(top50_genes.names), keytype = "ENSEMBL", column = "SYMBOL")

##selecinar as contagens normalizadas de acordo com os 50 genes mais diferencialmente expressos

mat <- counts(dds, normalized = TRUE)[rownames(top50_genes),]

##ajustar mat em z-score

mat_z <- t(apply(mat, 1, scale))

##nomear as colunas de mat_z

colnames(mat_z) <- colnames(mat)

##substituir o file_id de colnames(mat) pelo seu valor correspondente de submitter_id no dataframe de metadata

colnames(mat) <- metadata$submitter_id[match(colnames(mat), metadata$file_id)]

##gerar o heatmap

Heatmap(mat_z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat), name = "Z-SCORE", row_labels = top50_genes.names$symbol, row_names_gp = gpar(fontsize = 7))
