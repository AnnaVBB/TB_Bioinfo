import os
import pandas as pd
#CÓDIGO PARA EXTRAIR OS ARQUIVOS TSV DAS PASTAS E MONTAR A MATRIZ DESEQ (NÃO FUNCIONOU)

# Caminho da pasta onde os arquivos estão localizados
data_dir = "C:/Users/Usuário/Downloads/Bioinformatica/TCGA/COAD"  # Altere para o caminho correto
# Inicializa um dicionário para armazenar os dados de contagem
count_data = {}

# Percorrer todas as pastas no diretório de dados
for folder in os.listdir(data_dir):
    folder_path = os.path.join(data_dir, folder)
    
    if os.path.isdir(folder_path):
        # Cada pasta deve conter arquivos .tsv
        for file in os.listdir(folder_path):
            if file.endswith(".tsv"):
                file_path = os.path.join(folder_path, file)
                
                # Carregar o arquivo .tsv
                print(f"Lendo arquivo: {file_path}")  # Adiciona uma mensagem de depuração
                df = pd.read_csv(file_path, sep='\t')
                
                # Verifica se as colunas necessárias estão presentes
                if 'gene_id' in df.columns and 'unstranded' in df.columns:
                    # Use 'unstranded' ou 'fpkm_unstranded' ou o que for adequado para seu caso
                    count_data[folder] = df.set_index('gene_id')['unstranded'].to_dict()
                else:
                    print(f"Colunas ausentes no arquivo: {file_path}")  # Adiciona uma mensagem de depuração

# Criar um DataFrame a partir do dicionário
count_matrix = pd.DataFrame.from_dict(count_data, orient='index').fillna(0)

# Verifica se a matriz de contagem foi criada corretamente
if count_matrix.empty:
    print("A matriz de contagem está vazia. Verifique se os arquivos .tsv contêm dados válidos.")
else:
    # Transpor a matriz para que genes sejam linhas e amostras sejam colunas
    count_matrix = count_matrix.T

    # Salvar a matriz como um arquivo .tsv
    count_matrix.to_csv("matriz_deseq.tsv", sep='\t')

    print("Matriz DESeq salva em: matriz_deseq.tsv")


