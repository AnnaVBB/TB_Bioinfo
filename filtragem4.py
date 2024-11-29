import pandas as pd

#Adicionar uma coluna com valores numéricos para classificar os barcodes e colocar em ordem crescente, 
#de modo que barcodes iguais recebam o mesmo valor numérico
  
# Caminho para o arquivo com os dados
file_path = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\gdc_manifest.BRCA_LUAD_COAD.txt.map2submitterID"

# Carregar o arquivo em um DataFrame
df = pd.read_csv(file_path, sep='\t')

# Verificar se as colunas necessárias existem
if "cases.0.submitter_id" not in df.columns or "cases.0.samples.0.sample_type" not in df.columns:
    raise ValueError("As colunas 'cases.0.submitter_id' e/ou 'cases.0.samples.0.sample_type' não foram encontradas no arquivo.")

# Criar um mapeamento único para cada barcode
unique_barcodes = df['cases.0.submitter_id'].unique()
barcode_map = {barcode: idx + 1 for idx, barcode in enumerate(unique_barcodes)}

# Adicionar a nova coluna com os números sequenciais
df['barcode_numeric_id'] = df['cases.0.submitter_id'].map(barcode_map)

# Filtrar para manter apenas barcodes duplicados
barcode_counts = df['barcode_numeric_id'].value_counts()
duplicated_barcodes = barcode_counts[barcode_counts > 1].index

# Filtrar as linhas com barcodes duplicados
filtered_df = df[df['barcode_numeric_id'].isin(duplicated_barcodes)]

# Ordenar os dados pelo 'barcode_numeric_id' em ordem crescente
filtered_df_sorted = filtered_df.sort_values(by='barcode_numeric_id')

# Salvar o resultado em um novo arquivo
output_path = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\filtragem4.tsv"
filtered_df_sorted.to_csv(output_path, sep='\t', index=False)

print(f"Arquivo filtrado e organizado com IDs numéricos salvo em: {output_path}")
