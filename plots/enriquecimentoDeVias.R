#Instalar os pacotes
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db") 

#definir os clustters primeiro
#PCA (deseq2)

library(clusterProfiler)
library(org.Hs.eg.db)

# Converter os IDs dos genes para IDs Entrez (se necessário)
gene_list <- significant_genes$gene_id  # Supondo que 'gene_id' seja a coluna com os IDs dos genes
entrez_ids <- bitr(gene_list, fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Realizar a análise de enriquecimento de vias KEGG
kegg_enrichment <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                              organism = 'hsa', 
                              pvalueCutoff = 0.05)

# Visualizar os resultados
dotplot(kegg_enrichment)


'''
Análise de Enriquecimento de Gene Ontology (GO): Além das vias KEGG, você pode querer realizar uma análise de enriquecimento GO para entender as funções biológicas associadas aos genes diferencialmente expressos.

Visualização de Resultados: Utilize gráficos adicionais como gráficos de barras para as vias de sinalização ou GO que estão significativamente enriquecidas.

Validação de Resultados: Considere validar seus resultados usando dados de expressão gênica em outros conjuntos de dados ou experimentos de validação.

Análise Integrativa: Explore possíveis associações entre as vias identificadas e as características clínicas das amostras (normal vs. tumoral) para obter insights adicionais.


