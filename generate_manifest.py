import argparse
import requests
import csv

def usage():
    parser = argparse.ArgumentParser(description="Download GDC Manifest with Correct Headers")
    parser.add_argument("-i", "--input", dest="manifest", type=str, required=True, 
                        help="Input manifest file downloaded from GDC portal or a file containing file IDs in the first column.")
    return parser.parse_args()

def loadManifest(manifest_path):
    file_ids = []
    with open(manifest_path, 'r') as manifest_file:
        for line in manifest_file:
            columns = line.strip().split("\t")
            if len(columns[0]) == 36:  # Verifica IDs válidos
                file_ids.append(columns[0])
    return file_ids

def makeParams(file_ids):
    return {
        "filters": {
            "op": "in",
            "content": {
                "field": "files.file_id",
                "value": file_ids
            }
        },
        "format": "TSV",
        "fields": "file_id,file_name,cases.submitter_id,cases.samples.submitter_id,cases.samples.tissue_type",
        "size": len(file_ids)
    }

def gdcAPI(file_ids, manifest_path):
    url = "https://api.gdc.cancer.gov/files/"
    params = makeParams(file_ids)

    response = requests.post(url, json=params)
    if response.status_code != 200:
        raise Exception(f"API Error: {response.status_code} - {response.text}")
    
    # Processa o retorno da API
    output_path = manifest_path + ".map2submitterID.tsv"
    with open(output_path, 'w', newline='') as output_file:
        reader = csv.reader(response.text.splitlines(), delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')

        # Ajusta os cabeçalhos para refletir os campos desejados
        original_headers = next(reader)
        renamed_headers = [
            "submitter_sample","tissue_type","submitter_id", "file_id","file_name"
        ]
        writer.writerow(renamed_headers)

        for row in reader:
            # Cria um novo row com a ordem correta
            new_row = [
                row[0],  # submitter_sample
                row[1],  # tissue_type
                row[2],  # submitter_id
                row[3],  # file_id
                row[4],  # file_name
            ]
            writer.writerow(new_row)

def main():
    args = usage()
    manifest_path = args.manifest

    # Carrega os IDs do manifesto original
    file_ids = loadManifest(manifest_path)

    # Faz a solicitação para a API e gera o novo manifesto
    gdcAPI(file_ids, manifest_path)

if __name__ == "__main__":
    main()
