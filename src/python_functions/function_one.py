import pandas as pd
import sys
from itertools import chain
import json
import obonet
import os
import math
import numpy as np
from pydantic import BaseModel, TypeAdapter, field_validator
import requests
from typing import List,Optional, Dict, Union
import requests


API_URL = "https://cellguide.cellxgene.cziscience.com"

class AdditionalMetadata(BaseModel):
    """More specific properties from different sources"""
    cellxgene_canonical: Optional[list]
    cellxgene_computational: Optional[list]
    hubmap: Optional[list]

class BioQueryBiomarker(BaseModel):
    """A biomarker object for bioquery"""
    source: list
    type: str #should this be a list?
    symbol: str
    label: list
    additionalMetadata: AdditionalMetadata


    # TODO: More fields which should be supported by bioquery responses

class CellXGeneCanonicalMarkerGene(BaseModel):
    tissue: str
    symbol: str
    label: str
    publication: str
    publication_titles: str

class CellXGeneComputationalMarkerGene(BaseModel):
    me: float
    pc: float
    marker_score: float
    specificity: float
    gene_ontology_term_id: str
    symbol: str
    label: str
    groupby_dims: dict


class HubmapBiomarker(BaseModel):
    label:str
    hgnc_id:int
    symbol:str
    anatomical_structures:list
    cell_types:list
    




# # These are the response types from the cellguide API
# CellGuideCanonicalMarkerGenes = TypeAdapter(List[CellGuideCanonicalMarkerGene])
# CellGuideComputationalMarkerGenes = TypeAdapter(List[CellGuideComputationalMarkerGene])


# =====================================
# This block handles terminal input and file retrieval
# =====================================

#terminal input should be python function_one.py anatomical_system 
print("Arguments:", sys.argv[1:])


def get_files():
    return f"/Users/JenChen/Desktop/SIBMI/bionetquery/outputs/out.csv" #this is the combined file
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
    if data_ok:
        if data_driven_info_req.status_code == 200:
            data_driven_info = data_driven_info_req.json()
        else:
            print("did not find data driven markers")
    else:
            print("did not find data driven markers")
    if canonical_ok:
        if canonical_info_req.status_code == 200:
            canonical_info = canonical_info_req.json()
        else:
            print("did not find canonical markers")
    else:
        print("did not find canonical markers")
    # print(filtered_canonical_info)
    # unique_canonical_info = list(set(filtered_canonical_info))
    # unique_data_driven_info = list(set(filtered_data_driven_info))
    # print(unique_canonical_info)
    return canonical_info, data_driven_info

    
    
# =====================================
# This block handles HuBMAP data
# =====================================

def get_biomarkers(cell_type, file):
    
    df = pd.read_csv(file)
    print(cell_type)
    mask = df.apply(lambda row: row.str.contains(cell_type, case=False).any(), axis=1) #filter out the rows where the token is found in the cell type col
    filtered_rows = df[mask]
    # if len(sys.argv) == 2: ### TODO: Maybe better organize the names/ids/labels that result from this block because sometimes when you do unique/set you lose the connections between label and name/id
    #     biomarker_col = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID")]
    #     filtered_rows = filtered_rows[biomarker_col]
    # elif len(sys.argv) > 2: 
    #     if sys.argv[2]=="id":
    #         ids_only = [col for col in df.columns if col.startswith('BGene') and col.endswith("ID")]
    #         filtered_rows = filtered_rows[ids_only]
    #     elif sys.argv[2]=="name":
    #         names_only = [col for col in df.columns if col.startswith('BGene') and not col.endswith("ID") and not col.endswith("LABEL")]
    #         filtered_rows = filtered_rows[names_only]
   
    return filtered_rows
    

def only_check_ct_col(row, token):
    ct_columns = [col for col in row.index if col.startswith('CT')]
    for col in ct_columns:
        if token.lower() in str(row[col]).lower():
            return True
    return False

# =====================================
# Processing HuBMAP data
# =====================================
#add hgnc id and cl id 
# def isNaN(num): #i can't figure out how to handle nan
#     if num is not num:
#         if isinstance(num, str):
#             print("is nan")
#             return ""
#         if isinstance(num, int):
#             print("is nan")
#             return -1
#     else:
#         print("did not identify nan")
#         return num

def get_hubmap_biomarker_objects(filtered_rows):    
    hubmap_biomarkers = []
    
    for index, row in filtered_rows.iterrows():  # Use iterrows() to iterate over rows
        as_names = [name for name in row.index if name.startswith('AS') and not name.endswith('LABEL') and not name.endswith('ID')]
        ct_names = [name for name in row.index if name.startswith('CT') and not name.endswith('LABEL') and not name.endswith('ID')]
       
        for x in range(1, 15):  # Assuming this is the maximum number of biomarkers returned by Hubmap
            label_col = f'BGene/{x}/LABEL'
            id_col = f'BGene/{x}/ID'
            symbol_col = f'BGene/{x}'
            
            if label_col in row.index and id_col in row.index and symbol_col in row.index:
                gene_id = row[id_col]
                symbol_result = row[symbol_col]
                label_result = row[label_col]
                if pd.isna(gene_id):
                    gene_id = -1  
                
                if pd.isna(symbol_result):
                    symbol_result = "" 
                
                if pd.isna(label_result):
                    label_result = "" 
                
                print(gene_id, symbol_result, label_result)
                if gene_id is not -1: 
                    split = str(gene_id).split(':')
                    id_num = split[1]
                    gene_id = int(id_num)
                
                hubmap_biomarkers.append(HubmapBiomarker(
                    label=label_result,
                    hgnc_id=gene_id,
                    symbol=symbol_result,
                    anatomical_structures=as_names,
                    cell_types=ct_names
                ))
    
    print(hubmap_biomarkers)
    return hubmap_biomarkers

    # {
    #    anatomical structures: []
    #    cell types []

    # }
# =====================================
# Combining functions 
# =====================================

#TODO: making sure no label duplicates


def combine_biomarkers(canonical_marker_genes,computational_marker_genes, hubmap_biomarkers):
    merged_biomarkers = []

    merged_by_symbol = {}
    for marker in canonical_marker_genes:
        symbol = marker.symbol

        if symbol not in merged_by_symbol:
            #  new merged biomarker entry
            merged_by_symbol[symbol] = {
                'symbol': symbol,
                'source': ["cellxgene"],
                'type': "marker_gene",
                'label': [marker.label],
                'additionalMetadata': {
                    'cellxgene_canonical': [{
                        'tissue': marker.tissue,
                        'publication': marker.publication,
                        'publication_titles': marker.publication_titles
                    }],
                    'cellxgene_computational': None,
                    'hubmap': None  # Add hubmap data if available
                }
            }
        else:
            # update metadata
            merged_by_symbol[symbol]['additionalMetadata']['cellxgene_canonical'].append({
                'tissue': marker.tissue,
                'publication': marker.publication,
                'publication_titles': marker.publication_titles
            })

    for marker in computational_marker_genes:
        symbol = marker.symbol

        if symbol not in merged_by_symbol:
            # Create new biomarker entry
            merged_by_symbol[symbol] = {
                'symbol': symbol,
                'source': ["cellxgene"],
                'type': "marker_gene",
                'label': [marker.label],
                'additionalMetadata': {
                    'cellxgene_canonical': None,
                    'cellxgene_computational': {
                        'me': marker.me,
                        'pc': marker.pc,
                        'marker_score': marker.marker_score,
                        'specificity': marker.specificity,
                        'gene_ontology_term_id': marker.gene_ontology_term_id,
                        'groupby_dims': marker.groupby_dims
                    },
                    'hubmap': None  # Add hubmap data if available
                }
            }
        else:
            # update biomarker entry
            if merged_by_symbol[symbol]['additionalMetadata']['cellxgene_computational']:
                existing_computational = merged_by_symbol[symbol]['additionalMetadata']['cellxgene_computational']
                existing_computational.update({
                    'me': marker.me,
                    'pc': marker.pc,
                    'marker_score': marker.marker_score,
                    'specificity': marker.specificity,
                    'gene_ontology_term_id': marker.gene_ontology_term_id,
                    'groupby_dims': marker.groupby_dims
                })
            else:
                merged_by_symbol[symbol]['additionalMetadata']['cellxgene_computational'] = {
                    'me': marker.me,
                    'pc': marker.pc,
                    'marker_score': marker.marker_score,
                    'specificity': marker.specificity,
                    'gene_ontology_term_id': marker.gene_ontology_term_id,
                    'groupby_dims': marker.groupby_dims
                }

    for marker in hubmap_biomarkers:
        symbol = marker.symbol
        if symbol not in merged_by_symbol:
            # Create new biomarker entry
            merged_by_symbol[symbol] = {
                'symbol': symbol,
                'source': ["hubmap"],
                'type': "marker_gene",
                'label': [marker.label],
                'additionalMetadata': {
                    'cellxgene_canonical': None,
                    'cellxgene_computational': None,
                    'hubmap': {
                        'anatomical_structures': marker.anatomical_structures,
                        'cell_types': marker.cell_types
                    } 
                }
            }
        else:
            merged_by_symbol[symbol]['additionalMetadata']['hubmap'] = {
                'anatomical_structures': marker.anatomical_structures,
                'cell_types': marker.cell_types
                }
            merged_by_symbol[symbol]['source'].append("hubmap")

                                     
    # convert to list of biomarker dictionaries
    merged_biomarkers = list(merged_by_symbol.values())
        
    return merged_biomarkers
    

def process_biomarkers(canonical_info, data_driven_info, hubmap_biomarkers):
    print(data_driven_info)
    canonical_marker_genes = [
        CellXGeneCanonicalMarkerGene(
            tissue=entry['tissue'],
            symbol=entry['symbol'],
            label=entry['name'],
            publication=entry['publication'],
            publication_titles=entry['publication_titles']
        )
        for entry in canonical_info]
    computational_marker_genes = [
        CellXGeneComputationalMarkerGene(
            me=entry['me'],
            pc=entry['pc'],
            marker_score=entry['marker_score'],
            specificity=entry['specificity'],
            gene_ontology_term_id=entry['gene_ontology_term_id'],
            symbol=entry['symbol'],
            label=entry['name'],
            groupby_dims=entry['groupby_dims'],
        )
        for entry in data_driven_info]
    combined_biomarkers = combine_biomarkers(canonical_marker_genes, computational_marker_genes, hubmap_biomarkers)
    return combined_biomarkers
    # for entry in data_driven_info:
    #     filtered_data_driven_info.append({'sources': 'cellxgene', 'name': entry['name'], 'symbol': entry['symbol'], 'method': 'data_driven', 
    #             'gene_ontology_term_id': entry['gene_ontology_term_id'], 'me':entry['me'], 'pc': entry['pc'], })
    # for entry in canonical_info:
    #     filtered_canonical_info.append({'sources': 'cellxgene', 'name': entry['name'], 'symbol': entry['symbol']})

# def combine_hubmap_data(list_biomarkers):

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
    hubmap_markers = get_hubmap_biomarker_objects(hubmap_results)
    canonical_markers, data_driven_markers = get_biomarkers_from_cl(process_input())
    # print(data_driven_markers)
    # print(canonical_markers)
   
    processed_data = process_biomarkers(canonical_markers, data_driven_markers, hubmap_markers)
    json_string = json.dumps(processed_data, indent=4)
    output_file = "cellxgene2.json"

    with open(output_file, 'w') as f:
        f.write(json_string)

    print(f"JSON data has been written to {output_file}")
    # print(hubmap_results)
    # return hubmap_results
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

