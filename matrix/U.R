library(dplyr)
# Filtrar amostras de tumor primário e tecido normal
manifesto <- manifesto %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))

# Verificar as amostras disponíveis no manifesto
table(manifesto$sample_type)

# Vincular manifesto e dados clínicos pelo `bcr_patient_barcode`
clinical_manifest_combined <- clinical_data_combined %>%
  inner_join(manifesto, by = "bcr_patient_barcode")

# Garantir que temos ambos os tipos de amostra para cada paciente
paired_samples <- clinical_manifest_combined %>%
  group_by(bcr_patient_barcode) %>%
  filter(n() > 1)  # Inclui apenas os pacientes com tumor e tecido normal

# Separar as amostras de tumor e normal
tumor_samples <- paired_samples %>% filter(sample_type == "Primary Tumor")
normal_samples <- paired_samples %>% filter(sample_type == "Solid Tissue Normal")


# Função para carregar dados de expressão gênica
load_expression_data <- function(file_path) {
  data <- read.delim(file_path, header = TRUE, row.names = 1)
  return(data)
}

# Caminho base para os arquivos de expressão
base_path <- "C:/Users/Usuário/Downloads/Bioinformatica/LUAD/"

# Carregar os arquivos do manifesto
expression_data_list <- lapply(manifesto$file_name, function(file_name) {
  load_expression_data(paste0(base_path, file_name))
})

# Combinar os arquivos de expressão em uma matriz
expression_matrix <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE), expression_data_list)

# Ajustar nomes das colunas e linhas
rownames(expression_matrix) <- expression_matrix$Row.names
expression_matrix <- expression_matrix[, -1]

# Garantir que as colunas (amostras) da matriz de expressão correspondam ao manifesto
colnames(expression_matrix) <- sub("-[0-9]{2}[A-Z]$", "", colnames(expression_matrix))

# Filtrar para manter apenas amostras emparelhadas
common_samples <- intersect(colnames(expression_matrix), paired_samples$submitter_id_full)
expression_matrix <- expression_matrix[, common_samples]

# Vincular os dados clínicos a essas amostras
final_clinical_data <- paired_samples %>%
  filter(submitter_id_full %in% common_samples)

# Verificar dimensões finais
dim(expression_matrix)
dim(final_clinical_data)

library(DESeq2)

# Preparar os dados clínicos para DESeq2
coldata <- final_clinical_data %>%
  select(submitter_id_full, sample_type) %>%
  mutate(sample_type = ifelse(sample_type == "Primary Tumor", "tumor", "normal"))

# Criar o objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = expression_matrix,
                              colData = coldata,
                              design = ~ sample_type)

# Normalizar os dados e realizar análise de expressão diferencial
dds <- DESeq(dds)
results <- results(dds)

# Ver os genes diferencialmente expressos
head(results[order(results$padj), ])
