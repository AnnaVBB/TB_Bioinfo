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

