import pandas as pd

# Caminho do arquivo TSV
arquivo_tsv = "C:\\Users\\Usuário\\Downloads\\Bioinformatica\\TCGA\\filtragem4.tsv"

# Ler o arquivo TSV
df = pd.read_csv(arquivo_tsv, sep="\t")

# Extrair os barcodes únicos da coluna "case.0.submitter_id"
barcodes_unicos = df["cases.0.submitter_id"].unique()

# Converter o vetor para uma lista e exibir os valores separados por vírgulas
print(", ".join(barcodes_unicos)) 