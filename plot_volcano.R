#Volcano plot (excluindo variáveis)
# Carregar os pacotes necessários
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnhancedVolcano")
BiocManager::install("AnnotationDbi")

library(org.Hs.eg.db)
library(EnhancedVolcano)
library(AnnotationDbi)  

# --- Processar Dados ---

# Remover versões dos IDs ENSEMBL, caso presentes
rownames(res.df) <- gsub("\\..*", "", rownames(res.df))

# Filtrar IDs inválidos ou não relacionados a genes (e.g., N_ambiguous, N_unmapped)
valid_ids <- grep("^ENSG", rownames(res.df), value = TRUE)
res.df <- res.df[valid_ids, ]

# Mapear símbolos gene -> ENSEMBL com tratamento para múltiplos símbolos
mapped_genes <- mapIds(
  org.Hs.eg.db,
  keys = rownames(res.df),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = function(x) paste(unique(x), collapse = ", ")
)

# Adicionar os símbolos ao dataframe
res.df$symbol <- mapped_genes
res.df$log2FoldChange <- as.numeric(res.df$log2FoldChange)
res.df$padj <- as.numeric(res.df$padj)

# Verificar se há NAs no mapeamento
na_genes <- rownames(res.df)[is.na(res.df$symbol)]
if (length(na_genes) > 0) {
  warning("Os seguintes IDs não foram mapeados para nenhum símbolo: ", paste(na_genes, collapse = ", "))
}

head(res_sig)

library(EnhancedVolcano)

# --- Criar Volcano Plot ---
EnhancedVolcano(
  res.df,
  x = "log2FoldChange",     # Coluna do log2 Fold Change
  y = "padj",               # Coluna do p-value ajustado
  lab = res.df$symbol,                      # Rótulos (nomes dos genes ou símbolos)
  xlab = "Log2 Fold Change",  
  ylab = "-Log10 Adjusted P-value",
  xlim = c(-10, 10),                        # Define o intervalo do eixo X
  ylim = c(0, max(-log10(res.df$padj), na.rm = TRUE) + 5),     # Ajuste para o eixo Y (opcional)
  pCutoff = 0.05,           # Valor de p-value ajustado de corte
  FCcutoff = 2,             # Fold Change de corte
  pointSize = 1,            # Tamanho dos pontos
  labSize = 3.0,
  colAlpha = 0.6,
  legendPosition = "right"
)


#Verificando a consistencia
Counts <- counts(dds)
total_genes <- nrow(res.df)
print(paste("Genes iniciais:", nrow(raw_counts)))          # Total inicial
print(paste("Genes pós-remover zeros:", nrow(Counts)))     # Após remoção de linhas com zero
print(paste("Genes válidos:", length(valid_ids)))          # IDs ENSG válidos
print(paste("Genes mapeados:", sum(!is.na(mapped_genes)))) # Após mapeamento de IDs
print(paste("Genes analisados:", nrow(res.df)))            # Final para análise


