import pandas as pd
import sys
from itertools import chain
import json
import obonet

import os
import requests

# =====================================
# This block handles terminal input and file retrieval
# =====================================

#terminal input should be python function_one.py anatomical_system 
print("Arguments:", sys.argv[1:])


def get_files():
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/out.csv"
    # return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"

common_words = ["cell", "neuron"]
#if cell type has more than 1 word, filter out the tokens 


# =====================================
# This block handles cellxgene data
# =====================================

G = obonet.read_obo("http://purl.obolibrary.org/obo/cl/cl-basic.obo")

def get_cell_ontology_id(cell_type):
    #looking thru obonet first
    for node_id, data in G.nodes(data=True):
        if 'name' in data and data['name'] == cell_type:
            return node_id
        if 'synonym' in data:
            for syn in data['synonym']:
                if syn == cell_type:
                    return node_id

    req = requests.get(f"http://www.ebi.ac.uk/ols4/api/search?q={cell_type}&exact=false&ontology=cl")
    if req.status_code == 200:
        data = req.json()
        docs = data["response"]["docs"]
       
        if docs:
            first_id = docs[0]["obo_id"]
            return first_id
        #could possibly for loop to try other ids

        

# def find_cl_id(gene_name):
#     obo_id = get_cell_ontology_id(gene_name)
#     if obo_id == None:
#         pass
# The cellxgene API URLs contain a "latest snapshot identifier" (presumably points to a specific version of the data)
curr_id = requests.get("https://cellguide.cellxgene.cziscience.com/latest_snapshot_identifier").text
     

# cl_id = "CL:0000653" # Podocyte
cl_id = sys.argv[1]
     

canonical_info = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/canonical_marker_genes/{cl_id.replace(':', '_')}.json").json()
data_driven_info = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/computational_marker_genes/{cl_id.replace(':', '_')}.json").json()

canonical_load = json.loads(canonical_info)
data_driven = json.loads(data_driven_info)

print(canonical_info)
# def get_biomarkers_from_cl(string){
    
#     symbol = data["symbol"]
#     name = data["symbol"]
    
# }

# =====================================
# This block handles HuBMAP data
# =====================================

def get_biomarkers(cell_type, file):
    results = []
    df = pd.read_csv(file)
    tokens = cell_type.lower().split('_') #tokenize query
    # print(tokens)
    counter = 0
    for token in tokens: 
        counter+=1
        mask = df.apply(lambda row: only_check_ct_col(row, token), axis=1) #filter out the rows where the token is found in the cell type col
        filtered_rows = df[mask]
        if len(sys.argv) == 2:
            biomarker_col = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID")]
            filtered_rows = filtered_rows[biomarker_col]
        elif len(sys.argv) > 2:
            if sys.argv[2]=="id":
                ids_only = [col for col in df.columns if col.startswith('BGene') and col.endswith("ID")]
                filtered_rows = filtered_rows[ids_only]
            elif sys.argv[2]=="name":
                names_only = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID") and not col.endswith("LABEL")]
                filtered_rows = filtered_rows[names_only]
        non_empty_cols = filtered_rows.columns[filtered_rows.notna().any()]
        current_results = list((filtered_rows[non_empty_cols].values.tolist()))
        results = [item for sublist in current_results for item in sublist if pd.notna(item)]
        
    
    return results
    

def only_check_ct_col(row, token):
    ct_columns = [col for col in row.index if col.startswith('CT')]
    for col in ct_columns:
       
        if token.lower() in str(row[col]).lower():
            return True
    return False

# =====================================
# Search function
# =====================================

def search(): #method to call full search
    file_to_use = get_files()
    search_results = get_cell_ontology_id(sys.argv[1])
    print(search_results)
    # if(sys.argv[1].startswith("CL:")):
    #     search_results = get_cell_ontology_id(sys.argv[1])
    # else:
    #     search_results = get_biomarkers(sys.argv[1], file_to_use)
    # for result in search_results:
    #     pass
    #     print(result)
    #get_related_biomarkers("CL:0000895")

search()

