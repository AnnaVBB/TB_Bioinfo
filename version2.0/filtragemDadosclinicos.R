#para puxar e filtrar a matriz de dados clínicos referente ao TCGA (coad nesse caso), considerando o manifesto com os submitters IDS já prontos:
# Carregar TCGAbiolinks
library(TCGAbiolinks)
library(SummarizedExperiment)

#project- submitter_id--ajcc_pathologic_stage- gender --- vital_status -- age_at_index-

manifest_luad <- read.delim("C:/Users/Home/Documents/USP2024/bioinfo/tb1/manifest_luad_filtered.txt")
  
  ##puxar os dados clínicos
  clinical_luad <- GDCquery_clinic("TCGA-LUAD")
  
  ##filtrar as linhas que correspondem ao submitter ID do manifesto
  clinical_luad_filtered <- clinical_luad[clinical_luad$submitter_id %in% manifest_luad$submitter_id, ]
  
  ##selecionar as colunas interessantes, no caso submitter ID sempre, gender e age_at_index
  clinical_luad_filtered <- clinical_luad_filtered[, c(1,2,4,38, 40,41)]
 -------------------------------------------------------------------------------------------- 
  
  ##inserir as informacoes clinicas de clinical_coad_filtered em manifest coad
  manifest_luad$gender <- clinical_luad$gender[match(manifest_luad$submitter_id, clinical_luad$submitter_id)]
 manifest_luad$age_at_index <- clinical_luad$age_at_index[match(manifest_luad$submitter_id, clinical_luad$submitter_id)]
  
  ##transformar o manifest coad no dataframe de dados clinicos final, selecionando as colunas file_id sempre e as outras
  clinical.luad <- manifest_luad[, c(1, 3, 5, 6)]
  
 # Transformar os nomes das linhas em file id e remover a coluna file_id
rownames(clinical.coad) <- clinical.coad$file_id
View(clinical.coad)
clinical.coad <- clinical.coad[, -2]