import pandas as pd
import sys
import os
import csv
import requests

# import obonet
# G = obonet.read_obo("http://purl.obolibrary.org/obo/cl/cl-basic.obo")

nodes_file = "/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_nodes.csv"
nodes_df = pd.read_csv(nodes_file)
edges_file = "/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_edges.csv"
edges_df = pd.read_csv(edges_file)
cpdb_genes = "/Users/JenChen/Desktop/SIBMI/bionetquery/cellphonedb_data/gene_input.csv"
cpdb_gene_df = pd.read_csv(edges_file)
cpdb_interactions = "/Users/JenChen/Desktop/SIBMI/bionetquery/cellphonedb_data/interaction_input.csv"
cpdb_interaction_df = pd.read_csv(edges_file)


def get_files():
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/out.csv"
    # return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"

#given a cell type ID, find interacting cell types and return their biomarkers
def get_related_cells_and_biomarkers(cell_type_id):
    cell_ids = nodes_df.loc[nodes_df['cell_id'] == cell_type_id, 'id'].tolist()
    #there should only be one id per cell id 
    relevant_ids = []
    for id in cell_ids:
        filtered_edges = edges_df[(edges_df['source'] == id) | (edges_df['target'] == id)]
        relevant_ids = filtered_edges['source'].tolist() + filtered_edges['target'].tolist()
        relevant_ids = list(set(relevant_ids))
    relevant_biomarkers = []
    print(relevant_ids)
    for cellt in relevant_ids:
        print(cellt)
        df = pd.read_csv(get_files(sys.argv[1]))
        mask = df.apply(lambda row: check_ct_id_col(row, cellt), axis=1) #filter out the rows where the token is found in the cell type col
        filtered_rows = df[mask]
        biomarker_col = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID")]
        filtered_rows = filtered_rows[biomarker_col]
        relevant_biomarkers.append(filtered_rows)
    print(relevant_biomarkers)
    return relevant_biomarkers

def check_ct_id_col(row, token):
    ct_columns = [col for col in row.index if col.endswith('ID') and col.startswith('CT')]
    for col in ct_columns:
        if token in str(row[col]):
            return True
    return False

def search_cpdb(list_biomarkers) : #must be a list of their hgnc ids
    ###returns a list of the hgnc symbols of related biomarkers in cpdb
    found_symbols = []
    found_uniprot = []
    all_bm_symbols = []
    #get the hgnc symbols for all the biomarkers
    for bm in list_biomarkers:
        hgnc_symbol = get_symbol_from_hgnc_id(bm)
        all_bm_symbols.append(hgnc_symbol)
    
    #check if the symbols are in the cpdb gene csv. if it is get the uniprot ids
    with open(cpdb_genes, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            symbol_in_csv = row['hgnc_symbol']
            if symbol_in_csv in all_bm_symbols:
                uniprot_id = row['uniprot']
                found_uniprot.append(uniprot_id)
                found_symbols.append(symbol_in_csv)
    #find interaction partners if any
    related_bm_uniprot = []
    with open(cpdb_interactions, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            partner_a = row['partner_a']
            partner_b = row['partner_b']
            if partner_a == found_uniprot:
                related_bm_uniprot.append(partner_b)
            if partner_b == found_uniprot:
                related_bm_uniprot.append(partner_a)

    related_hgnc_symbol = []
    with open(cpdb_genes, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            uniprot_id_to_check = row['uniprot']
            if uniprot_id_to_check == related_bm_uniprot:
                related_hgnc_symbol.append(row['hgnc_symbol'])
    return related_hgnc_symbol


practice_list = ['130', '2155']



    
        


            

    
def get_symbol_from_hgnc_id(hgnc_id):
    url = f"https://rest.genenames.org/search/hgnc_id/{hgnc_id}"
    headers = {
        'Accept': 'application/xml'
    }
    
    try:
        response = requests.get(url, headers=headers)
        if response.status_code == 200:
           
            xml_response = response.text
          
            import xml.etree.ElementTree as ET
            root = ET.fromstring(xml_response)
            
            symbol = root.find(".//str[@name='symbol']").text
            print(symbol)
            return symbol
        else:
            print(f"Error {response.status_code}: {response.reason}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return None

print(search_cpdb(practice_list))