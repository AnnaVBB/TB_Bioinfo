import pandas as pd

# Caminhos dos arquivos
manifesto_filtrado_path = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/LUAD_manifesto_clinico.tsv" #Mudar LUAD por BRCA, COAD ou outro que queira
manifesto_original_path = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/gdc_manifest.BRCA_LUAD_COAD.txt"
output_manifesto_path = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/LUAD_manifesto_completo.txt" #Mudar LUAD por BRCA, COAD ou outro que queira

# Carregar o manifesto filtrado e o manifesto original
manifesto_filtrado = pd.read_csv(manifesto_filtrado_path, sep='\t')
manifesto_original = pd.read_csv(manifesto_original_path, sep='\t')

# Verificar se a coluna 'id' existe em ambos os DataFrames
if 'id' not in manifesto_filtrado.columns or 'id' not in manifesto_original.columns:
    raise ValueError("A coluna 'id' deve estar presente em ambos os manifestos.")

# Realizar a união (merge) dos dados usando 'id' como chave
merged_manifesto = pd.merge(manifesto_filtrado, manifesto_original, on='id', how='left')

# Verificar se as colunas obrigatórias agora estão presentes
colunas_obrigatorias = ['id', 'file_name', 'md5', 'size', 'state']
for coluna in colunas_obrigatorias:
    if coluna not in merged_manifesto.columns:
        raise ValueError(f"A coluna obrigatória '{coluna}' não foi encontrada no manifesto unido.")

# Selecionar apenas as colunas obrigatórias
merged_manifesto = merged_manifesto[colunas_obrigatorias]

# Salvar o novo manifesto
merged_manifesto.to_csv(output_manifesto_path, sep='\t', index=False)

print(f"Manifesto unido salvo em: {output_manifesto_path}")
