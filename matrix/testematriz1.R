# Suponha que o manifesto está em um arquivo ou já carregado
manifesto <- read.delim('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/manifesto.txt', header = TRUE) # Substitua pelo seu arquivo
head(manifesto)

# Obtenha os IDs das colunas na ordem desejada a partir do manifesto
manifest_ids <- manifesto$Sample_ID # Substitua 'Sample_ID' pela coluna do arquivo com os IDs

# Verifique os IDs presentes no raw_counts e no manifesto
colnames(raw_counts) # IDs atuais da matriz
head(manifest_ids)    # IDs desejados (do manifesto)

# Reorganizar a matriz
ordered_cols <- match(manifest_ids, colnames(raw_counts))
raw_counts <- raw_counts[, ordered_cols, drop = FALSE] # Organiza conforme o manifesto

# Verifique a matriz reorganizada
head(raw_counts)

saveRDS(raw_counts, file = 'C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/TCGA_COAD_RNASEQ_REORDERED.rds')
