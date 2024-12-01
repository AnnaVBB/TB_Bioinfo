library(ggplot2)

# Carregar os dados
results <- read.csv("C:/Users/Usuário/Downloads/Bioinformatica/TCGA/TCGA_bioinfo/genes_diferencialmente_expressos_coad.csv")

#Primeiro, você precisa carregar os dados e identificar os genes com maior log2FoldChange e menor padj.

# Filtrar genes significativos (padj < 0.05 é um critério comum)
significant_genes <- results[results$padj < 0.05, ]

# Identificar os genes com maior log2FoldChange
max_fc_gene <- significant_genes[which.max(significant_genes$log2FoldChange), ]
min_fc_gene <- significant_genes[which.min(significant_genes$padj), ]

# Exibir os resultados
print(max_fc_gene)
print(min_fc_gene)

results_clean <- results[!is.na(results$log2FoldChange) & !is.na(results$padj), ]


# Criar gráfico de volcano
ggplot(results_clean, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("violet", "darkviolet"), 
                     labels = c("not significant", "significant")) +  # Mudar os rótulos
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal() +
  # Adicionar rótulos apenas para genes significativos
  geom_text(data = significant_genes[significant_genes$log2FoldChange > 2 | significant_genes$log2FoldChange < -2, ],
            aes(label = ""), 
            check_overlap = TRUE, size = 3) +
  guides(color = guide_legend(title = "Significance (padj<0.05"))  # Título da legenda



