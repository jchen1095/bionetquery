import requests


KEGG_URL = 'https://rest.kegg.jp'

def get_kegg_image_with_pathway(path_names): #MUST BE A KEGG PATHWAY NAME
    
    try:
        for name in path_names:
            response = requests.get(f"{KEGG_URL}/get/{name}")
            response.raise_for_status()
            with open(f'./outputs/images/image_{name}.jpg', 'wb') as f:
                f.write(response.content)
            print('successfully retrieved image')
    except requests.exceptions.RequestException as e:
        print('no image for pathway', e)


def get_kegg_pathways(gene_id): #MUST BE KEGG ID
    path_ids = []
    try:
        response = requests.get(f"{KEGG_URL}/link/pathway/{gene_id}")
        response.raise_for_status()
        content = response.text.splitlines()  
        for line in content: 
            print(line)
            path = line.split('\t')
            for part in path:
                if part.startswith('path:'):
                    path_id = part.split('path:')[1]
                    path_id = path_id.strip()
                    path_ids.append(path_id)
    except requests.exceptions.RequestException as e:
        print('no pathways found', e)
    print(path_ids)
    return path_ids


def get_kegg_id(gene_id): #assumes uniprot id, can also use ncbi ids but you need to change the code
    try:
        response = requests.get(f"{KEGG_URL}/conv/genes/uniprot:{gene_id}")
        response.raise_for_status()
        content = response.text
        kegg_id = content.split('\t')
        print(kegg_id)
        for part in kegg_id:
            if part.startswith('hsa:'):
                kegg_id = part.strip()
                print(kegg_id)
                return kegg_id
        else:
            raise ValueError("KEGG ID not found in response")

    except requests.exceptions.RequestException as e:
        print(f"reqiest to KEGG API failed: {e}")
    except ValueError as ve:
        print(f"Did not retrieve KEGG ID: {ve}")

practice = ["P31213"]

def main():
    for gene in practice: 
        id = get_kegg_id(gene)
        pathways = get_kegg_pathways(id)
        get_kegg_image_with_pathway(pathways)

main()