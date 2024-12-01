library(DESeq2)
library(dplyr)

# Preparar os dados clínicos para DESeq2
coldata <- coad_manifesto %>%
  select(file_id, cases.0.samples.0.tissue_type) %>%
  rename(sample_type = cases.0.samples.0.tissue_type)

# Converter 'sample_type' para fator explicitamente
coldata$sample_type <- as.factor(coldata$sample_type)

if (any(is.na(coldata$sample_type))) {
  stop("Existem valores NA em 'sample_type'. Verifique os dados clínicos.")
}

# Converter 'sample_type' para fator explicitamente
coldata$sample_type <- as.factor(coldata$sample_type)

# Certificar-se de que os IDs do manifesto e da matriz estão alinhados
# Verifica se os IDs em coldata$file_id estão em aligned_RNASEQ_COAD
missing_ids <- setdiff(coldata$file_id, colnames(aligned_RNASEQ_COAD))
if (length(missing_ids) > 0) {
  stop("Os seguintes IDs do manifesto não estão na matriz RNA-seq: ", paste(missing_ids, collapse = ", "))
}

aligned_RNASEQ_COAD <- aligned_RNASEQ_COAD[, coldata$file_id]

# Criar o objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = aligned_RNASEQ_COAD,
                              colData = coldata,
                              design = ~ cases.0.samples.0.tissue_type)

# Normalizar os dados e realizar análise de expressão diferencial
dds <- DESeq(dds)
results <- results(dds)

# Ver os genes diferencialmente expressos
head(results[order(results$padj), ])

# Salvar os resultados em um arquivo
write.csv(results, "genes_diferencialmente_expressos_coad.csv", row.names = TRUE)
