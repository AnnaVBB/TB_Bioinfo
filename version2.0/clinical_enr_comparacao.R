#Corrigir o código, não está funcionando

'''
Integração com dados clínicos ou experimentais
Relacione as vias enriquecidas com características clínicas específicas (ex.: sobrevida, estágio tumoral).
Analise a expressão dos genes associados às vias exclusivas ou compartilhadas usando dados experimentais (ex.: RNA-seq).
Identifique padrões de expressão gênica ou atividade de vias que diferem entre os cânceres.
'''
View(clinical_coad_final)
##### Relacionar os DEGs ou vias enriquecidas (GO/KEGG) com dados clínicos, como idade, estágio patológico e vitalidade.
library(survival)
library(survminer)
library(DESeq2)

#Selecionar os top 10 genes 
# Converter os resultados do DESeq2 para um data frame
dds <- DESeq(dds)
res_df <- res.df
normalized_counts <- counts(dds, normalized = TRUE)


# Filtrar os 10 genes mais significativos ordenados pelo p-valor ajustado
top_genes <- res_df %>%
  filter(!is.na(padj)) %>%    # Remover linhas com NA no p-valor ajustado
  arrange(padj) %>%           # Ordenar por p-valor ajustado
  head(10)                    # Selecionar os 10 primeiros genes

print(top_genes)

rownames(top_genes) <- gsub("\\..*$", "", rownames(top_genes))  # Remove sufixos, se existirem
gene_ids <- intersect(rownames(normalized_counts), rownames(top_genes))
if (length(gene_ids) == 0) {
  stop("Nenhum gene correspondente encontrado entre normalized_counts e top_genes.")
}

gene_expression <- normalized_counts[gene_ids, ]

# Kaplan-Meier para um gene DEG específico
# Verifica e limpa os dados de expressão
gene_expression <- normalized_counts[rownames(normalized_counts) %in% rownames(top_genes), ]
aggregate_expression <- rowMeans(gene_expression)

if (nrow(gene_expression) == 0) {
  stop("Nenhum gene correspondente encontrado em normalized_counts.")
}

# Agrupa a expressão
expression_group <- ifelse(aggregate_expression > median(aggregate_expression), "High", "Low")

# Atualiza os dados clínicos
survival_data <- clinical_coad_final
survival_data$expression_group <- expression_group[match(survival_data$submitter_id, names(expression_group))]
table(survival_data$expression_group, useNA = "ifany")

# Remove linhas com valores ausentes em expression_group
survival_data <- survival_data[!is.na(survival_data$expression_group), ]

# Verifica e corrige a variável vital_status
survival_data$vital_status <- trimws(tolower(survival_data$vital_status))  # Normaliza para minúsculas e remove espaços
survival_data$vital_status <- ifelse(survival_data$vital_status == "dead", "Dead", "Alive")  # Uniformiza para "Dead" e "Alive"

# Remove linhas com valores ausentes em vital_status
survival_data <- survival_data[!is.na(survival_data$vital_status), ]

# Verifica o número de observações restantes
print(nrow(survival_data))
print(table(survival_data$expression_group, useNA = "ifany"))
print(table(survival_data$vital_status, useNA = "ifany"))

# Ajusta o modelo Kaplan-Meier
fit <- survfit(Surv(age_at_index, vital_status == "Dead") ~ expression_group, data = survival_data)

# Plota o gráfico de Kaplan-Meier
library(survminer)
ggsurvplot(fit, data = survival_data, pval = TRUE, 
           risk.table = TRUE, conf.int = TRUE, 
           palette = c("#E41A1C", "#377EB8"),
           legend.title = "Expression Group",
           legend.labs = c("High", "Low"))

#Investigar se pacientes com diferentes estágios ou subtipos tumorais têm vias enriquecidas distintas.
##Subdivida os pacientes pelos valores de ajcc_pathologic_stage (por exemplo, estágio I vs estágio IV).
##Compare o enriquecimento de vias GO/KEGG entre os subgrupos, priorizando os termos compartilhados e exclusivos.
##Use gráficos como heatmaps para visualizar a associação entre os estágios e as vias mais significativas.

library(pheatmap)

# Criar uma matriz com -log10(p.adjust) para diferentes estágios
go_matrix <- data.frame(
  StageI = -log10(ego_stageI$p.adjust),
  StageIV = -log10(ego_stageIV$p.adjust)
)
rownames(go_matrix) <- ego_stageI$Description

pheatmap(go_matrix, cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE)


#Avaliar como variáveis clínicas afetam a sobrevivência ou a expressão de DEGs
##Use modelos de regressão logística para predizer vital_status (morto/vivo) com base em age_at_index, ajcc_pathologic_stage, e outros fatores.
##Realize uma regressão multivariada de Cox para avaliar o impacto combinado de fatores clínicos na sobrevivência.

cox_model <- coxph(Surv(age_at_index, vital_status == "dead") ~ 
                     gender + ajcc_pathologic_stage + tissue_type, 
                   data = clinical_data)
summary(cox_model)


#Relacionar enriquecimento funcional de vias (KEGG/GO) com características clínicas.
##Combine os resultados de enriquecimento funcional (GO/KEGG) com a tabela clínica, associando pacientes a vias enriquecidas.
# 1. Preparar os dados
# Certifique-se de que clinical_coad_final e enrichment_results têm a coluna "submitter_id"
# e que os IDs estão consistentes.

if (!"submitter_id" %in% colnames(clinical_coad_final)) {
  stop("A coluna 'submitter_id' não está presente nos dados clínicos.")
}
enrichment_results$submitter_id <- paste0("Sample_", 1:nrow(enrichment_results))
merged_data <- merge(clinical_coad_final, enrichment_results, by = "submitter_id")

print(head(merged_data))

# Verificar se as colunas necessárias existem no merged_data
if (!"age_at_index" %in% colnames(merged_data) || !"pathway_score" %in% colnames(merged_data)) {
  stop("As colunas 'age_at_index' ou 'pathway_score' estão ausentes em merged_data.")
}
# Correlação
correlation_result <- cor.test(merged_data$age_at_index, merged_data$pathway_score)

# Exibir resultados
print(correlation_result)
