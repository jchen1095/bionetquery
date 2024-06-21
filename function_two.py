import pandas as pd
import sys
import os
import requests

# import obonet
# G = obonet.read_obo("http://purl.obolibrary.org/obo/cl/cl-basic.obo")

nodes_file = "/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_nodes.csv"
nodes_df = pd.read_csv(nodes_file)
edges_file = "/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_edges.csv"
edges_df = pd.read_csv(edges_file)


def get_files():
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/out.csv"
    # return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"

#given a cell type ID, find interacting cell types and return their biomarkers
def get_related_biomarkers(cell_type_id):
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
