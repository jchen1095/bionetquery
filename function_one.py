import pandas as pd
import sys
import os
import requests

# The cellxgene API URLs contain a "latest snapshot identifier" (presumably points to a specific version of the data)
curr_id = requests.get("https://cellguide.cellxgene.cziscience.com/latest_snapshot_identifier").text
     

cl_id = "CL:0000653" # Podocyte
# cl_id = sys.argv[1]
     

canonical_info = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/canonical_marker_genes/{cl_id.replace(':', '_')}.json").json()
data_driven_info = requests.get(f"https://cellguide.cellxgene.cziscience.com/{curr_id}/computational_marker_genes/{cl_id.replace(':', '_')}.json").json()
     

print(canonical_info)


#terminal input should be python function_one.py anatomical_system 
print("Arguments:", sys.argv[1:])

def get_files():
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/out.csv"
    # return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"

common_words = ["cell", "neuron"]
#if cell type has more than 1 word, filter out the tokens 


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
        if sys.argv[2]== None:
            biomarker_col = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID")]
            filtered_rows = filtered_rows[biomarker_col]
        if sys.argv[2]=="id":
            ids_only = [col for col in df.columns if col.startswith('BGene') and col.endswith("ID")]
            filtered_rows = filtered_rows[ids_only]
        elif sys.argv[2]=="name":
            names_only = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID") and not col.endswith("LABEL")]
            filtered_rows = filtered_rows[names_only]
        non_empty_cols = filtered_rows.columns[filtered_rows.notna().any()]
        results.append(filtered_rows[non_empty_cols])
    
    return results
    

def only_check_ct_col(row, token):
    ct_columns = [col for col in row.index if col.startswith('CT')]
    for col in ct_columns:
       
        if token.lower() in str(row[col]).lower():
            return True
    return False

def search(): #method to call full search
    file_to_use = get_files()
    search_results = get_biomarkers(sys.argv[1], file_to_use)
    for result in search_results:
        print(result)
    #get_related_biomarkers("CL:0000895")

search()

