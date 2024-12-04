#Volcano plot2
# Loading relevant libraries 
library(tidyverse)  # Includes ggplot2 for visualization, dplyr for data manipulation
library(RColorBrewer)  # For custom color palettes
library(ggrepel)  # For better annotations

# Load the data
results <- read.csv("C:/Users/Usuário/Downloads/TCGA/TCGA_COAD_DESeq_significant_genes.csv")
results <- results %>% drop_na(log2FoldChange, pvalue)

# Criar a coluna diffexpressed com base nos limiares de significância
results <- results %>%
  mutate(
    diffexpressed = case_when(
      log2FoldChange > 2 & pvalue < 0.5 ~ "Superexpresso",
      log2FoldChange < -2 & pvalue < 0.5 ~ "Subexpresso",
      TRUE ~ "Não significativo"
    )
  )


# Biostatsquid theme
theme_set(
  theme_classic(base_size = 20) +
    theme(
      axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1.0), color = 'black'),
      axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20, 0, 0, 0), size = rel(1.0), color = 'black'),
      plot.title = element_text(hjust = 0.5)
    )
)

# Create a basic volcano plot
ggplot(data = results, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("#00AFBB", "grey", "red"),  # Custom colors
    labels = c("Subexpresso", "Não significativo", "Superexpresso")  # Labels
  )
