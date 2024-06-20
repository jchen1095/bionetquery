import pandas as pd
import sys
import os

#terminal input should be python function_one.py anatomical_system 
print("Arguments:", sys.argv[1:])

def get_files(string):
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/hubmap_asct_data/hubmap_{string}.csv"


def get_biomarkers(cell_type, file):
    results = []
    df = pd.read_csv(file)
    tokens = [cell_type.lower()] + cell_type.lower().split('_') #tokenize query
    print(tokens)
    counter = 0
    for token in tokens: 
        counter+=1
        mask = df.apply(lambda row: only_check_ct_col(row, token), axis=1) #filter out the rows where the token is found in the cell type col
        filtered_rows = df[mask]
        biomarker_col = [col for col in df.columns if col.startswith('BGene')]
        filtered_rows = filtered_rows[biomarker_col]
        results.append(filtered_rows)
    
    print(counter)
    return results
    

def only_check_ct_col(row, token):
    ct_columns = [col for col in row.index if col.startswith('CT')]
    
    for col in ct_columns:
       
        if token.lower() in str(row[col]).lower():
            return True
    return False

def search(): #method to call full search
    file_to_use = get_files(sys.argv[1])
    search_results = get_biomarkers(sys.argv[2], file_to_use)
    for result in search_results:
        print(result)

search()
