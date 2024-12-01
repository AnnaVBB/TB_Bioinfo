#Genes diferencialmente expressos
#Vias de sinalização do câncer (Gene Set Enrichment Analysis)
#LUAD, COAD e BRCA 
#Tecido tumoral x Normal 

#Visualização e filtragem dos dados clínicos 
####Instalação dos pacotes
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("data.table", quietly = TRUE))
  install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# Instalar pacotes Bioconductor
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))

# Instalar pacotes CRAN
install.packages(c("survminer", "survival", "tidyverse"))

###bibliotecas 
library(TCGAbiolinks)
library(survminer)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)
library(DESeq2)
library(dplyr)

# getting clinical data for TCGA-BRCA, TCGA-LUAD and TCGA-COAD cohort
clinical_brca <- GDCquery_clinic("TCGA-BRCA")
clinical_luad <- GDCquery_clinic("TCGA-LUAD")
clinical_coad <- GDCquery_clinic("TCGA-COAD")
#filtrar as colunas (selecionar as que forem relevantes para o trabalho)
columns_to_keep <- c("project", 
                     "submitter_id", 
                     "tissue_or_organ_of_origin", 
                     "tumor_grade", 
                     "ajcc_pathologic_t", 
                     "classification_of_tumor", 
                     "diagnosis_id", 
                     "gender", 
                     "ethnicity", 
                     "age_at_index", 
                     "vital_status", 
                     "bcr_patient_barcode")

clinical_brca_filtered <- clinical_brca %>%
  select(all_of(columns_to_keep))

clinical_luad_filtered <- clinical_luad %>%
  select(all_of(columns_to_keep))

clinical_coad_filtered <- clinical_coad %>%
  select(all_of(columns_to_keep))

# Combinar os dados em uma única tabela (opcional, pois iremos trabalhar com os dados separadamente)
clinical_data_combined <- bind_rows(
  clinical_brca_filtered,
  clinical_luad_filtered,
  clinical_coad_filtered
)

#Transformando a tabela em arquivo .tsv para manipula-la no python (vamos unir os dados do manifesto com os dados clínicos)
# Caminho para salvar o arquivo
output_path <- "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/clinical_data_combined.tsv"

# Salvar a tabela como um arquivo .tsv
write.table(clinical_data_combined, 
            file = output_path, 
            sep = "\t",        # Define o separador como tabulação
            row.names = FALSE, # Não incluir os números das linhas
            quote = FALSE)     # Não adicionar aspas ao redor dos valores

cat("Tabela salva em:", output_path, "\n")

#Depois de feita a combinação das tabelas, se quiser visualiza-las aqui no R, utilizar:
brca_manf_clinc <- read.delim("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/BRCA_manifesto_clinico.tsv", 
                              header = FALSE, sep = "\t")
luad_manf_clinc <- read.delim("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/LUAD_manifesto_clinico.tsv", 
                              header = FALSE, sep = "\t")
coad_manf_clinc <- read.delim("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD_manifesto_clinico.tsv", 
                              header = FALSE, sep = "\t")
View(brca_manf_clinc)
View(luad_manf_clinc)
View(coad_manf_clinc)
'''













