import pandas as pd
import tarfile
import os

# Adicionar a qual projeto/tipo de cancer cada id do manifesto corresponde através dos dados clínicos. 
# Combina os dataframes com os dados processados em um arquivo csv


# Lista de arquivos TAR e os tipos de câncer correspondentes
tar_files = {
    "LUAD": r'C:\Users\Usuário\Downloads\Bioinformatica\TCGA\LUAD_clinical.gz', 
    "BRCA": r'C:\Users\Usuário\Downloads\Bioinformatica\TCGA\BRCA_clinical.gz',
    "COAD": r'C:\Users\Usuário\Downloads\Bioinformatica\TCGA\COAD_clinical.gz'
}

# Caminho para extrair os arquivos
base_extract_path = r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\Extraido'

# Colunas que queremos manter
columns_to_remain = [
    'case_id', 'case_submitter_id', 'age_at_index', 'gender', 
    'race', 'ethnicity', 'tumor_grade', 'primary_site', 
    'tissue_or_organ_of_origin', 'tumor_stage', 'year_of_diagnosis'
]

# Função para processar e limpar os arquivos TAR
def process_tar_file(tar_path, cancer_type):
    extract_path = os.path.join(base_extract_path, cancer_type)
    os.makedirs(extract_path, exist_ok=True)  # Criar pasta se não existir

    # Extrair os arquivos
    with tarfile.open(tar_path, 'r:gz') as tar_ref:
        tar_ref.extractall(extract_path)

    # Procurar pelo arquivo clínico específico
    clinical_file = os.path.join(extract_path, 'clinical.tsv')
    if os.path.exists(clinical_file):
        # Processar o arquivo clínico
        df = pd.read_csv(clinical_file, sep='\t')
        # Selecionar apenas as colunas desejadas
        df_cleaned = df[columns_to_remain].copy()
        df_cleaned['cancer_type'] = cancer_type  # Adicionar coluna identificadora
        return df_cleaned
    else:
        print(f"Arquivo 'clinical.tsv' não encontrado em {extract_path}")
        return pd.DataFrame()  # Retornar DataFrame vazio se não encontrar o arquivo

# Processar cada arquivo TAR e combinar os resultados
dataframes = []
for cancer_type, tar_path in tar_files.items():
    print(f"Processando {cancer_type}...")
    df_cleaned = process_tar_file(tar_path, cancer_type)
    if not df_cleaned.empty:
        dataframes.append(df_cleaned)

# Combinar todos os DataFrames em um único
combined_df = pd.concat(dataframes, ignore_index=True)

# Exibir as colunas finais e uma amostra dos dados
print("Colunas combinadas:", combined_df.columns.tolist())
print(combined_df.head())

# Salvar o resultado em um arquivo CSV (opcional)
output_path = r'C:\Users\Usuário\Downloads\Bioinformatica\Projeto\clinical_combined.csv'
combined_df.to_csv(output_path, index=False)
print(f"Dados combinados salvos em: {output_path}")
