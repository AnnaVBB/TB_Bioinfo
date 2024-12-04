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

# Verificar se há NAs no mapeamento
na_genes <- rownames(res.df)[is.na(res.df$symbol)]
if (length(na_genes) > 0) {
  warning("Os seguintes IDs não foram mapeados para nenhum símbolo: ", paste(na_genes, collapse = ", "))
}

# --- Criar Volcano Plot ---
EnhancedVolcano(
  res.df,
  x = "log2FoldChange",
  y = "padj",
  lab = res.df$symbol,
  xlab = "Log2 Fold Change",  
  ylab = "-Log10 Adjusted P-value",
  xlim = c(-5, 5),
  ylim = c(0, max(-log10(res.df$padj), na.rm = TRUE) + 5),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 0.5,
  labSize = 3.0,
  colAlpha = 0.6,
  legendPosition = "right"
)
