import pandas as pd 
import requests 

# Carregar o arquivo .map2submitterID em um DataFrame
file_path = 'C:\\Users\\Usuário\\Downloads\\Bioinformatica\\Projeto\\gdc_manifest_TCGA_LG_maf.txt.map2submitterID'
data_df = pd.read_csv(file_path, sep='\t', low_memory=False)

# Verificar as primeiras linhas do DataFrame para entender a estrutura dos dados
print(data_df.columns)

'''
# Verificar a estrutura do DataFrame para identificar a coluna de tipo de amostra
print(data_df.columns)  # Isso ajuda a confirmar o nome da coluna relevante

# Filtrar amostras tumorais
tumor_samples = data_df[data_df['cases.samples.sample_type'].str.contains('Primary Tumor', na=False)]

# Filtrar amostras normais/benignas
normal_samples = data_df[data_df['cases.samples.sample_type'].str.contains('Solid Tissue Normal', na=False)]

# Exibir um resumo dos resultados
print("Amostras tumorais:")
print(tumor_samples[['file_id', 'file_name', 'cases.samples.submitter_id']].head())

print("\nAmostras normais:")
print(normal_samples[['file_id', 'file_name', 'cases.samples.submitter_id']].head())

tumor_samples.to_csv('tumor_samples.csv', index=False)
normal_samples.to_csv('normal_samples.csv', index=False)
'''