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
    return f"./bionetquery/outputs/out.csv" #this is the combined file
    # return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"

common_words = ["cell", "neuron"]
#if cell type has more than 1 word, filter out the tokens 


# =====================================
# This block handles getting the CL id given a cell type input
# =====================================

G = obonet.read_obo("http://purl.obolibrary.org/obo/cl/cl-basic.obo")
#TODO: cl-full instead of cl-basic (takes longer to retrieve)

def get_cell_ontology_id(cell_type):
    #looking thru obonet first
    for node_id, data in G.nodes(data=True):
        if 'name' in data and data['name'] == cell_type:
            print("found in name")
            return node_id
        if 'synonym' in data:
            for syn in data['synonym']:
                if syn == cell_type:
                    print("found in synonym")
                    return node_id

    req = requests.get(f"http://www.ebi.ac.uk/ols4/api/search?q={cell_type}&exact=false&ontology=cl")
    if req.status_code == 200:
        data = req.json()
        docs = data["response"]["docs"]
        if docs:
            first_id = docs[0]["obo_id"]
            print(docs[0])
            print("found in api")
            return first_id
        
        #could possibly for loop to try other ids

        
# =====================================
# This block handles retrieving markers from cellxgene 
# =====================================


# The cellxgene API URLs contain a "latest snapshot identifier" (presumably points to a specific version of the data)
curr_id = requests.get("https://cellguide.cellxgene.cziscience.com/latest_snapshot_identifier").text
     
print(curr_id)
# cl_id = "CL:0000653" # Podocyte

     
# NOTE: "Marker genes are not available for blood or small populations of cells"  https://cellxgene.cziscience.com/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_5__Find%20Marker%20Genes

def get_biomarkers_from_cl(string):
    cl_id = str(get_cell_ontology_id(string))
    print(cl_id)
    canonical_info_req = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/canonical_marker_genes/{cl_id.replace(':', '_')}.json")
    data_driven_info_req = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/computational_marker_genes/{cl_id.replace(':', '_')}.json")
    data_ok = data_driven_info_req.text
    canonical_ok = canonical_info_req.text
    filtered_data_driven_info = []
    filtered_canonical_info = []
    if data_ok:
        if data_driven_info_req.status_code == 200:
            data_driven_info = data_driven_info_req.json()
            for entry in data_driven_info:
                filtered_data_driven_info.append({'name': entry['name'], 'symbol': entry['symbol'], 'gene_ontology_term_id': entry['gene_ontology_term_id']})
        else:
            print("did not find data driven markers")
    else:
            print("did not find data driven markers")
    if canonical_ok:
        if canonical_info_req.status_code == 200:
            canonical_info = canonical_info_req.json()
            for entry in canonical_info:
                filtered_canonical_info.append({'name': entry['name'], 'symbol': entry['symbol']})
        else:
            print("did not find canonical markers")
    else:
        print("did not find canonical markers")
    # print(filtered_canonical_info)
    # unique_canonical_info = list(set(filtered_canonical_info))
    # unique_data_driven_info = list(set(filtered_data_driven_info))
    # print(unique_canonical_info)
    return filtered_canonical_info, filtered_data_driven_info

    
    
# =====================================
# This block handles HuBMAP data
# =====================================

def get_biomarkers(cell_type, file):
    results = []
    df = pd.read_csv(file)
    print(cell_type)
    mask = df.apply(lambda row: row.str.contains(cell_type, case=False).any(), axis=1) #filter out the rows where the token is found in the cell type col
    filtered_rows = df[mask]
    if len(sys.argv) == 2: ### TODO: Maybe better organize the names/ids/labels that result from this block because sometimes when you do unique/set you lose the connections between label and name/id
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
    unique_results = list(set(results)) #TODO: do this later/downstream not at the level ur creating
    return unique_results
    

def only_check_ct_col(row, token):
    ct_columns = [col for col in row.index if col.startswith('CT')]
    for col in ct_columns:
        if token.lower() in str(row[col]).lower():
            return True
    return False

# =====================================
# Search function
# =====================================

def process_input():
    input=sys.argv[1]
    processed_string = input.replace('_', ' ')
    processed_string = str(processed_string)
    return processed_string


def search(): #method to call full search
    file_to_use = get_files()
    hubmap_results = get_biomarkers(process_input(), file_to_use)
    canonical_markers, data_driven_markers = get_biomarkers_from_cl(process_input())
    print(data_driven_markers)
    print(canonical_markers)
    print(hubmap_results)
    return hubmap_results
    #print(data_driven_markers)
    # if(sys.argv[1].startswith("CL:")):
    #     search_results = get_cell_ontology_id(sys.argv[1])
    # else:
    #     search_results = get_biomarkers(sys.argv[1], file_to_use)
    # for result in search_results:
    #     pass
    #     print(result)
    #get_related_biomarkers("CL:0000895")

search()

#TODO: Combine biomarkers but also provide data source/reasoning ("why is this biomarker included in list of results") [hubmap, cellxgene]
#would allow a user to filter downstream 
#unify/standardize schema of all the database returns 
#users don't *have* to know where data comes from but if they want to they can find it 
#if users want to add another database make it easier for them to process all the data into this one standardize return schema 
#Another property: biomarker type (canonical or data driven) ("method") because then you can include everything in the same list and user can downstream filter
#right now hubmap is only canonical 
#description, symbol, hgnc id, cl id, synonyms
#TODO: all the downstream filtering can be applied later 


