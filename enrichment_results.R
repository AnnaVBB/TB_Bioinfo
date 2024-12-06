
# Carregar dados de enriquecimento
ego_coad <- read.csv("GO_enrichment_COAD.csv", row.names = 1)
ekegg_coad <- read.csv("KEGG_enrichment_COAD.csv", row.names = 1)

ego_luad <- read.csv("GO_enrichment_LUAD.csv", row.names = 1)
ekegg_luad <- read.csv("KEGG_enrichment_LUAD.csv", row.names = 1)

# Consolidar GO e KEGG para COAD e LUAD
ego_coad$Cancer <- "COAD"
ekegg_coad$Cancer <- "COAD"
ego_luad$Cancer <- "LUAD"
ekegg_luad$Cancer <- "LUAD"

# Adicionar um identificador único para as vias
ego_coad$ID <- paste("GO", ego_coad$Description, sep = ":")
ekegg_coad$ID <- paste("KEGG", ekegg_coad$Description, sep = ":")
ego_luad$ID <- paste("GO", ego_luad$Description, sep = ":")
ekegg_luad$ID <- paste("KEGG", ekegg_luad$Description, sep = ":")

# Combinar os resultados em um único data frame
enrichment_results <- rbind(
  ego_coad[, c("ID", "Cancer", "p.adjust")],
  ekegg_coad[, c("ID", "Cancer", "p.adjust")],
  ego_luad[, c("ID", "Cancer", "p.adjust")],
  ekegg_luad[, c("ID", "Cancer", "p.adjust")]
)

# Renomear colunas para clareza
colnames(enrichment_results) <- c("submitter_id", "Cancer", "p.adjust")

# Adicionar escore de enriquecimento (-log10(p.adjust))
enrichment_results$pathway_score <- -log10(enrichment_results$p.adjust)

# Verificar a estrutura
print(head(enrichment_results))
