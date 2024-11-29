import pandas as pd

#Selecionar apenas as linhas que possuam mais de um barcode  

# Caminho para o arquivo com os dados
file_path = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\gdc_manifest.BRCA_LUAD_COAD.txt.map2submitterID"

# Carregar o arquivo em um DataFrame
# Substitua '\t' por ',' se o arquivo for CSV
df = pd.read_csv(file_path, sep='\t')

# Verificar se a coluna existe no DataFrame
if "cases.0.submitter_id" not in df.columns:
    raise ValueError("A coluna 'cases.0.submitter_id' não foi encontrada no arquivo.")

# Contar as ocorrências de cada submitter_id
duplicates = df[df['cases.0.submitter_id'].duplicated(keep=False)]

# Agrupar os duplicados e organizar sequencialmente
result = duplicates.groupby('cases.0.submitter_id').apply(
    lambda group: group.reset_index(drop=True)
).reset_index(drop=True)

# Salvar o resultado em um novo arquivo
output_path = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\gdc_manifest.BRCA_LUAD_COAD.txt.map2submitterID.tsv"
result.to_csv(output_path, sep='\t', index=False)

print(f"Arquivo filtrado salvo em: {output_path}")
