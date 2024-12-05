# Carregar pacotes necessários
library(org.Hs.eg.db)
library(AnnotationDbi)

# Copiar a tabela results_df para res_df2
res_df2<- results_df

# Verificar os rownames originais
print("IDs originais nos rownames:")
head(rownames(res_df2))

# Remover versões (ponto seguido de número) dos IDs Ensembl
cleaned_ids <- sub("\\..*", "", rownames(res_df2))  # Remover sufixos após o ponto

# Garantir que os IDs sejam únicos
if (anyDuplicated(cleaned_ids) > 0) {
  print("IDs duplicados encontrados após remover versões. Ajustando para garantir unicidade.")
  cleaned_ids <- make.unique(cleaned_ids)  # Tornar os IDs únicos
}

# Atualizar os rownames de res_df2 com os IDs corrigidos
rownames(res_df2) <- cleaned_ids

# Verificar os rownames corrigidos
print("IDs corrigidos (sem versões e únicos):")
head(rownames(res_df2))

# Adicionar a coluna SYMBOL mapeando ENSEMBL -> SYMBOL
res_df2$symbol <- mapIds(
  org.Hs.eg.db,                # Banco de dados org.Hs.eg.db
  keys = rownames(res_df2),     # IDs Ensembl corrigidos (rownames de res_df2)
  keytype = "ENSEMBL",         # Tipo de ID de entrada
  column = "SYMBOL",           # Tipo de ID de saída
  multiVals = "first"          # Caso haja múltiplos valores, usar o primeiro
)

# Verificar se a coluna foi adicionada corretamente
print("Tabela com coluna SYMBOL adicionada:")
head(res_df2)

# Salvar a tabela final em um arquivo CSV (opcional)
write.csv(res_df2, "res_df2_with_symbol.csv", row.names = TRUE)
