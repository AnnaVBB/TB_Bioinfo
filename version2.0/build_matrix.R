library(dplyr)
#DEseq e análises

# Carregar metadados
metadata <- read.delim("C:/Users/Usuário/Downloads/TCGA/manifest_coad_filtered.txt", header = TRUE)

# Obter IDs das colunas da matriz de contagem (excluindo gene_id)
matrix_ids <- setdiff(colnames(raw_counts), "gene_id")
metadata_ids <- uuid_to_tcga$file_id

# Diagnóstico: IDs não correspondentes
missing_ids <- setdiff(matrix_ids, metadata_ids)

# Exibir IDs que estão na matriz, mas não nos metadados
if (length(missing_ids) > 0) {
  warning("IDs ausentes nos metadados: ", paste(missing_ids, collapse = ", "))
  stop("Erro: Existem IDs na matriz de contagem que não têm correspondência nos metadados.")
}

# Criar o dataframe de condições somente com IDs correspondentes
sample_conditions <- data.frame(
  row.names = matrix_ids,
  condition = uuid_to_tcga$condition[match(matrix_ids, uuid_to_tcga$file_id)]
)

# Verificar se todas as amostras têm condição associada
if (any(is.na(sample_conditions$condition))) {
  stop("Erro: Algumas amostras não têm condição associada. Verifique o mapeamento.")
}

print(missing_ids)

# Verificar a matriz de condições
head(sample_conditions)

# --- Processar Matriz de Contagem ---
# Remover coluna de IDs dos genes, se presente
if ("gene_id" %in% colnames(raw_counts)) {
  row.names(raw_counts) <- raw_counts$gene_id
  raw_counts <- raw_counts[, -1]
}

# Garantir que a matriz é numérica
raw_counts <- as.matrix(raw_counts)
mode(raw_counts) <- "numeric"

# Verificar se há valores negativos
if (any(raw_counts < 0)) {
  warning("A matriz contém valores negativos. Substituindo por 0.")
  raw_counts[raw_counts < 0] <- 0
}

# Remover genes com valores ausentes
raw_counts <- raw_counts[complete.cases(raw_counts), ]

# Verificar estrutura final da matriz de contagem
summary(as.vector(raw_counts))
all(raw_counts == floor(raw_counts))  # Deve retornar TRUE

# --- Verificar Dados ---
print(dim(raw_counts))
print(table(sample_conditions$condition))

# Converter condição para fator
sample_conditions$condition <- as.factor(sample_conditions$condition)


# --- Criar Objeto DESeqDataSet ---
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = sample_conditions,
  design = ~ condition
)

# --- Análise com DESeq ---
dds <- DESeq(dds)

# Extrair resultados
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Ordenar por p-valor ajustado
res <- res[order(res$padj), ]

# Visualizar os primeiros genes
head(res)

# Salvar o objeto DESeqDataSet
saveRDS(dds, file = 'C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_DESeq.rds')

##### Gráfico MA
plotMA(res, main = "MA Plot - Tumor vs Normal")

# Filtrar Genes Significativos
significant_genes <- res[which(res$padj < 0.05), ]
write.csv(significant_genes, file = "C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_DESeq_significant_genes.csv", row.names = TRUE)

library(ggplot2)

# Adicionar coluna para significância
res$significant <- ifelse(res$padj < 0.05, "Yes", "No")

# Criar MA plot com ggplot2
ggplot(res, aes(x = log2(baseMean), y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "MA Plot", x = "log2(Mean Normalized Counts)", y = "log2(Fold Change)") +
  theme(legend.position = "top")
