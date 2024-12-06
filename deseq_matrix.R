# Instalar o BiocManager se não estiver instalado
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar DESeq2
BiocManager::install("DESeq2")

library(DESeq2)

''' Se você tem os UUIDs de arquivos na matriz de contagem e 
os IDs TCGA (amostra/submissor) no manifesto, o objetivo é mapear 
os IDs do manifesto para os IDs da matriz de contagem e criar a 
matriz sample_conditions para análise diferencial.'''

# Carregar metadados
metadata <- read.delim("C:/Users/Usuário/Downloads/TCGA/manifest_coad_filtered.txt", header = TRUE)

# Inspecionar colunas para ver correspondência
head(metadata)

# Criar um mapeamento entre UUID e Submitter ID
uuid_to_tcga <- metadata[, c("file_id", "submitter_sample")]

# Inspecionar o mapeamento
head(uuid_to_tcga)

# Adicionar condição com base no código de tipo de tecido no submitter_sample
# Sufixos "-11A", "-11B" indicam amostras normais; tudo o mais é considerado tumoral
uuid_to_tcga$condition <- ifelse(
  grepl("-11[A-Z]?$", uuid_to_tcga$submitter_sample), "Normal", "Tumor"
)

# Verificar o mapeamento com condições
head(uuid_to_tcga)

# Selecionar apenas as amostras presentes na matriz de contagem
sample_conditions <- data.frame(
  row.names = colnames(raw_counts),  # UUIDs das colunas da matriz
  condition = uuid_to_tcga$condition[match(colnames(raw_counts), uuid_to_tcga$file_id)]
)

# Verificar a matriz de condições
head(sample_conditions)

###Verificação da matriz
# Remover coluna de IDs (caso necessário)
row.names(raw_counts) <- raw_counts$gene_id
raw_counts <- raw_counts[, -1]

# Garantir que é uma matriz numérica
raw_counts <- as.matrix(raw_counts)
mode(raw_counts) <- "numeric"

# Verificar se há valores negativos
summary(as.vector(raw_counts))

# Garantir que todos os valores são inteiros
all(raw_counts == floor(raw_counts))

####

# Remover genes com valores ausentes
raw_counts <- raw_counts[complete.cases(raw_counts), ]

# Verificar se há valores negativos na matriz de contagem
summary(as.vector(raw_counts))

# Substituir valores negativos por 0
raw_counts[raw_counts < 0] <- 0

# Criar o objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = sample_conditions,
                              design = ~ condition)


dds <- DESeq(dds)
# Extrair resultados
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Ordenar por p-valor ajustado
res <- res[order(res$padj), ]

# Visualizar os primeiros genes
head(res)


# Salvar o objeto DESeqDataSet se necessário
saveRDS(dds, file = 'C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_DESeq.rds')
dds <- DESeq(dds)
results_dds <- results(dds)
