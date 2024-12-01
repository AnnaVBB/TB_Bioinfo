# Carregar o manifesto clínico
coad_manifesto <- read.delim("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD_manifesto_clinico.tsv", 
                             header = TRUE, sep = "\t")

# Visualizar o manifesto (opcional)
View(coad_manifesto)

# Verificar se a matriz de RNA-seq já está carregada no ambiente
if (!exists("TCGA_COAD_RNASEQ")) {
  stop("A matriz de RNA-seq (TCGA_COAD_RNASEQ) não foi encontrada no ambiente.")
}

# Garantir que os IDs do manifesto estejam na coluna correta
if (!"file_id" %in% colnames(coad_manifesto)) {
  stop("A coluna 'file_id' não foi encontrada no manifesto clínico.")
}

# Alinhar as colunas da matriz de RNA-seq com os IDs do manifesto
aligned_RNASEQ_COAD <- TCGA_COAD_RNASEQ[, coad_manifesto$file_id, drop = FALSE]

# Visualizar a matriz alinhada (opcional)
View(aligned_RNASEQ_COAD)

#Verificar estão na mesma ordem
all(colnames(aligned_RNASEQ_COAD) == coad_manifesto$cases.0.samples.0.tissue_type)

saveRDS(aligned_RNASEQ_COAD, file = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/aligned_RNASEQ_COAD.rds")

