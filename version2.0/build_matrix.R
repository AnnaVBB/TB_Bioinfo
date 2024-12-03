library(dplyr)

# Carregar manifesto
manifesto <- read.delim('C:/Users/Usuário/Downloads/TCGA/manifest_coad_filtered.txt', header = TRUE)
manifest_ids <- manifesto$file_id  # Ajuste para a coluna que contém os IDs no manifesto
print(manifest_ids)  # Verifique os IDs do manifesto

# Listar diretórios de RNA-Seq disponíveis
rnaseq_coad.dirs <- list.files('C:/Users/Usuário/Downloads/TCGA/COAD/')
print(rnaseq_coad.dirs)

# Filtrar diretórios com base nos IDs do manifesto
filtered_dirs <- rnaseq_coad.dirs[rnaseq_coad.dirs %in% manifest_ids]
missing_ids <- setdiff(manifest_ids, filtered_dirs)

# Aviso sobre IDs ausentes
if (length(missing_ids) > 0) {
  warning("Os seguintes IDs estão ausentes nos diretórios: ", paste(missing_ids, collapse = ", "))
}

# Ordenar diretórios conforme a ordem no manifesto
ordered_dirs <- filtered_dirs[match(manifest_ids, filtered_dirs)]
print(ordered_dirs)  # Verifique a ordem

# Verificar se há IDs ausentes no manifesto
if (any(is.na(ordered_dirs))) {
  stop("Erro: Há IDs no manifesto que não correspondem a nenhum diretório disponível.")
}

# Inicializar matriz com o primeiro diretório
first_file <- list.files(paste0('C:/Users/Usuário/Downloads/TCGA/COAD/', ordered_dirs[1]), pattern = '.tsv')
raw_counts <- read.delim(paste0('C:/Users/Usuário/Downloads/TCGA/COAD/', ordered_dirs[1], '/', first_file), skip = 1)
raw_counts <- raw_counts[, c("gene_id", "unstranded")]
colnames(raw_counts)[2] <- ordered_dirs[1]

# Loop para adicionar os demais diretórios
for (i in ordered_dirs[-1]) {
  file_name <- list.files(paste0('C:/Users/Usuário/Downloads/TCGA/COAD/', i), pattern = '.tsv')
  tmp <- read.delim(paste0('C:/Users/Usuário/Downloads/TCGA/COAD/', i, '/', file_name), skip = 1)
  tmp <- tmp[, c("gene_id", "unstranded")]
  colnames(tmp)[2] <- i
  raw_counts <- merge(raw_counts, tmp, by = "gene_id", all = TRUE)  # Manter todos os genes, mesmo se estiverem ausentes em algum arquivo
}

# Garantir que a coluna gene_id seja a primeira
raw_counts <- raw_counts %>%
  select(gene_id, everything())  # Mover gene_id para a primeira coluna

# Transformar a matriz para numérica, excluindo a coluna gene_id
raw_counts_matrix <- as.matrix(raw_counts[, -1])  # Remove a coluna gene_id para a matriz
rownames(raw_counts_matrix) <- raw_counts$gene_id  # Adiciona os IDs dos genes como nomes das linhas

# Verificar se a matriz é numérica
if (!all(sapply(raw_counts_matrix, is.numeric))) {
  stop("Erro: Valores na matriz 'raw_counts_matrix' não são numéricos.")
}

# Salvar matriz gerada
saveRDS(raw_counts_matrix, file = 'C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_RNASEQ_ORDERED.rds')
