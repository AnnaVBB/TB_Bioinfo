# Carregar pacotes necessários
library(ggplot2)

# Carregar os resultados de expressão diferencial do CSV
results_df <- read.csv("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/TCGA_bioinfo/genes_diferencialmente_expressos_coad.csv", 
                       row.names = 1)  # Use row.names = 1 se a primeira coluna for o gene_id

# Adicionar uma coluna para indicar se o gene é significativo
results_df$significant <- ifelse(results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1, "Significant", "Not Significant")

# Criar o gráfico de Volcano
ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "green")) +
  labs(title = "Volcano Plot", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

