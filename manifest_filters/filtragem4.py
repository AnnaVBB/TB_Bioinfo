import pandas as pd

# Combinar os dados clíncos com o arquivo de manifesto filtrado. 
# Separar em três arquivos de acordo com o tipo de câncer 

# Caminhos dos arquivos
manifesto_path = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/filtragem4.tsv"  
clinical_path = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/clinical_data_combined.tsv" 

# Carregar os dados do manifesto
manifesto = pd.read_csv(manifesto_path, sep='\t')

# Carregar os dados clínicos
clinical_data = pd.read_csv(clinical_path, sep='\t')

# Verificar se as colunas necessárias existem
if "cases.0.submitter_id" not in manifesto.columns or "submitter_id" not in clinical_data.columns:
    raise ValueError("As colunas 'cases.0.submitter_id' no manifesto ou 'submitter_id' nos dados clínicos não foram encontradas.")

# Renomear colunas para facilitar o merge
manifesto = manifesto.rename(columns={"cases.0.submitter_id": "submitter_id"})

# Realizar a união (merge) dos dados
combined_data = pd.merge(manifesto, clinical_data, on="submitter_id", how="inner")

# Verificar se a coluna do tipo de câncer existe
if "project" not in combined_data.columns:
    raise ValueError("A coluna 'project' com o tipo de câncer não foi encontrada nos dados clínicos.")

# Exibir a quantidade de amostras por tipo de câncer
cancer_types = combined_data["project"].value_counts()

print("Quantidade de amostras por tipo de câncer:")
print(cancer_types)

# Exibir os primeiros registros das tabelas combinadas
print("\nPré-visualização dos dados combinados:")
print(combined_data.head())


# Caminhos de saída para os arquivos
output_dir = "C:/Users/Usuário/Downloads/Bioinformatica/"
brca_output = output_dir + "BRCA_manifesto_clinico.tsv"
luad_output = output_dir + "LUAD_manifesto_clinico.tsv"
coad_output = output_dir + "COAD_manifesto_clinico.tsv"

# Verificar se a coluna 'project' está presente nos dados combinados
if "project" not in combined_data.columns:
    raise ValueError("A coluna 'project' com o tipo de câncer não foi encontrada nos dados combinados.")

# Separar os dados por tipo de câncer
brca_data = combined_data[combined_data["project"].str.contains("BRCA", case=False, na=False)]
luad_data = combined_data[combined_data["project"].str.contains("LUAD", case=False, na=False)]
coad_data = combined_data[combined_data["project"].str.contains("COAD", case=False, na=False)]

# Salvar os dados filtrados em arquivos separados
brca_data.to_csv(brca_output, sep='\t', index=False)
luad_data.to_csv(luad_output, sep='\t', index=False)
coad_data.to_csv(coad_output, sep='\t', index=False)

# Mensagens de confirmação
print(f"Dados do projeto BRCA salvos em: {brca_output}")
print(f"Dados do projeto LUAD salvos em: {luad_output}")
print(f"Dados do projeto COAD salvos em: {coad_output}")
