import pandas as pd

#Selecionar uma amostra normal e uma tumoral para cada barcode (portanto, para barcodes duplicados no documento)
#Excluindo as linhas que possuam somente um tipo de amostra 

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

# Filtrar barcodes duplicados
duplicates = df[df['cases.0.submitter_id'].duplicated(keep=False)]

# Criar um DataFrame para armazenar as amostras selecionadas
filtered_samples = pd.DataFrame()

# Iterar pelos barcodes únicos duplicados
for barcode in duplicates['cases.0.submitter_id'].unique():
    # Selecionar todas as amostras para o barcode atual
    barcode_samples = duplicates[duplicates['cases.0.submitter_id'] == barcode]
    
    # Selecionar uma amostra tumoral (Primary Tumor)
    tumor_sample = barcode_samples[barcode_samples['cases.0.samples.0.sample_type'] == "Primary Tumor"].head(1)
    
    # Selecionar uma amostra normal (Solid Tissue Normal)
    normal_sample = barcode_samples[barcode_samples['cases.0.samples.0.sample_type'] == "Solid Tissue Normal"].head(1)
    
    # Concatenar as amostras selecionadas ao DataFrame final
    filtered_samples = pd.concat([filtered_samples, tumor_sample, normal_sample], ignore_index=True)

# Salvar o resultado em um novo arquivo
output_path = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\filtragem3.tsv"
filtered_samples.to_csv(output_path, sep='\t', index=False)

print(f"Arquivo filtrado com IDs numéricos salvo em: {output_path}")
