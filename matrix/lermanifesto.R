# Suponha que o manifesto está em um arquivo ou já carregado
# Ler o manifesto
manifesto <- read.delim('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD_manifesto_clinico.tsv', header = TRUE)
manifest_ids <- manifesto$Sample_ID # Ajuste para a coluna correta
print(manifest_ids) # Verifique os IDs do manifesto

# Listar diretórios disponíveis
rnaseq_coad.dirs <- list.files('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/')
print(rnaseq_coad.dirs)

# Identificar quais IDs do manifesto estão disponíveis
filtered_dirs <- rnaseq_coad.dirs[rnaseq_coad.dirs %in% manifest_ids]
missing_ids <- setdiff(manifest_ids, filtered_dirs)

# Aviso sobre IDs ausentes
if (length(missing_ids) > 0) {
  warning("Os seguintes IDs estão ausentes nos diretórios: ", paste(missing_ids, collapse = ", "))
}

# Ordenar os diretórios conforme o manifesto
ordered_dirs <- filtered_dirs[match(manifest_ids, filtered_dirs)]
print(ordered_dirs) # Verifique a ordem

# Verificar se há NA em ordered_dirs
if (any(is.na(ordered_dirs))) {
  stop("Há IDs no manifesto que não correspondem a nenhum diretório disponível.")
}

# Inicializar a matriz usando o primeiro diretório
raw_counts <- read.delim(
  paste0('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/TCGA_COAD_RNASEQ/', ordered_dirs[1], '/', 
         list.files(paste0('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/TCGA_COAD_RNASEQ', ordered_dirs[1]), pattern = '.tsv')), 
  skip = 1
)[-c(1:4), ]
raw_counts <- raw_counts[, c(1, 4)]
colnames(raw_counts)[2] <- ordered_dirs[1]

# Loop para adicionar os demais diretórios na ordem do manifesto
for (i in ordered_dirs[-1]) {
  tmp <- read.delim(
    paste0('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/', i, '/', 
           list.files(paste0('C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD/', i), pattern = '.tsv')), 
    skip = 1
  )[-c(1:4), ]
  tmp <- tmp[, 4]
  raw_counts <- cbind(raw_counts, tmp)
  colnames(raw_counts)[ncol(raw_counts)] <- i
}

# Ajustar nomes das linhas e converter para matriz numérica
rownames(raw_counts) <- raw_counts$gene_id
raw_counts <- as.matrix(raw_counts[, -1])

# Certifique-se de que raw_counts é um data frame ao atribuir os nomes das linhas
rownames(raw_counts) <- raw_counts[, 1]  # Assume-se que a primeira coluna contém os IDs dos genes
raw_counts <- raw_counts[, -1]           # Remover a coluna usada como nomes das linhas

# Transformar para matriz numérica
raw_counts <- as.matrix(raw_counts)

# Verificar se todos os valores são numéricos
if (!is.numeric(raw_counts)) {
  stop("Erro: Valores na matriz 'raw_counts' não são numéricos.")
}

# Salvar a matriz final
saveRDS(raw_counts, file = 'C:/Users/Usuário/Downloads/Bioinformatica/TCGA/TCGA_COAD_RNASEQ_ORDERED.rds')

# Verificar resultados
head(raw_counts)
dim(raw_counts)

# Salvar a matriz gerada
saveRDS(raw_counts, file = 'C:/Users/Usuário/Downloads/Bioinformatica/TCGA/TCGA_COAD_RNASEQ_ORDERED.rds')

# Verificar o resultado final
head(raw_counts)
dim(raw_counts)
