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
from pydantic_classes import AdditionalMetadata,BioQueryBiomarker,CellXGeneCanonicalMarkerGene,CellXGeneComputationalMarkerGene,HubmapBiomarker
from enum import Enum

class MarkerSource(Enum):
    CELLXGENE_CANONICAL = 'cellxgene_canonical'
    CELLXGENE_COMPUTATIONAL = 'cellxgene_computational'
    HUBMAP = 'hubmap'
    CELLPHONEDB = 'cellphonedb'

class OutputType(Enum):
    DICTIONARY = 'dictionary'
    LIST_SYMBOLS = 'list_symbols'
    LIST_UNIPROT = 'list_uniprot'


API_URL = "https://cellguide.cellxgene.cziscience.com"
data_file = "./outputs/out.csv"


# TODO: Unit testing/examples
# Utils/helper functions
# TODO: More fields which should be supported by bioquery responses

# =====================================
# This block handles terminal input and file retrieval
# =====================================

#terminal input should be python function_one.py anatomical_system 
print("Arguments:", sys.argv[1:])

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
            print("found in api")
            return first_id
        
        #could possibly for loop to try other ids

        
# =====================================
# This block handles retrieving markers from cellxgene 
# =====================================


# The cellxgene API URLs contain a "latest snapshot identifier" (presumably points to a specific version of the data)
curr_id = requests.get("https://cellguide.cellxgene.cziscience.com/latest_snapshot_identifier").text
     
# cl_id = "CL:0000653" # Podocyte

     
# NOTE: "Marker genes are not available for blood or small populations of cells"  https://cellxgene.cziscience.com/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_5__Find%20Marker%20Genes

def get_biomarkers_from_cl(string):
    cl_id = str(get_cell_ontology_id(string))

    canonical_info_req = requests.get(f"{API_URL}/{curr_id}/canonical_marker_genes/{cl_id.replace(':', '_')}.json")
    data_driven_info_req = requests.get(f"{API_URL}/{curr_id}/computational_marker_genes/{cl_id.replace(':', '_')}.json")
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

#TODO: how will you handle duplicates of data?
    
    
# =====================================
# This block handles HuBMAP data
# =====================================

def get_hubmap_biomarkers(cell_type, file):
    
    df = pd.read_csv(file)
    mask = df.apply(lambda row: row.str.contains(cell_type, case=False).any(), axis=1) #filter out the rows where the token is found in the cell type col
    filtered_rows = df[mask]
   
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

def get_hubmap_biomarker_objects(cell_type, filtered_rows):    
    hubmap_biomarkers = []
    ct_names_all = []
    for index, row in filtered_rows.iterrows():  # Use iterrows() to iterate over rows

        as_names = [row[name] for name in row.index if name.startswith('AS') and not name.endswith('LABEL') and not name.endswith('ID') and not pd.isna(row[name])]
        ct_names = [row[name] for name in row.index if name.startswith('CT') and not name.endswith('LABEL') and not name.endswith('ID') and not pd.isna(row[name])]
        if not any(
            str(ct_name).strip().lower() == cell_type.lower() or 
            f" {cell_type.lower()}" in str(ct_name).strip().lower()
            for ct_name in ct_names
        ):
            continue

        ct_names_all.append((index, ct_names))

        for x in range(1, 15):  # Assuming this is the maximum number of biomarkers returned by Hubmap
            label_col = f'BGene/{x}/LABEL'
            id_col = f'BGene/{x}/ID'
            symbol_col = f'BGene/{x}'
            
            if label_col in row.index and id_col in row.index and symbol_col in row.index:
                gene_id = row[id_col]
                symbol_result = row[symbol_col]
                label_result = row[label_col]
                if pd.isna(gene_id):
                    gene_id = ""  
                
                if pd.isna(symbol_result):
                    symbol_result = "" 
                
                if pd.isna(label_result):
                    label_result = "" 
                
                hubmap_biomarkers.append(HubmapBiomarker(
                    label=label_result,
                    id=gene_id,
                    symbol=symbol_result,
                    anatomical_structures=as_names,
                    cell_types=ct_names
                ))
    
        output_file = 'ct_names'
        # Write CT names to file
        with open(output_file, 'w') as f:
            for index, ct_names in ct_names_all:
                f.write(f'Row {index}: {ct_names}\n')
    return hubmap_biomarkers

# =====================================
# Combining functions 
# =====================================

#TODO: making sure no label duplicates
def make_cellxgene_canonical_metadata(biomarker):
    return {
        'tissue': biomarker.tissue,
        'publication': biomarker.publication,
        'publication_titles': biomarker.publication_titles
    }

def make_cellxgene_computational_metadata(biomarker):
    return {
        'me': biomarker.me,
        'pc': biomarker.pc,
        'marker_score': biomarker.marker_score,
        'specificity': biomarker.specificity,
        'gene_ontology_term_id': biomarker.gene_ontology_term_id,
        'groupby_dims': biomarker.groupby_dims
    }

def make_hubmap_metadata(biomarker):
    return {
        'anatomical_structures': biomarker.anatomical_structures,
        'cell_types': biomarker.cell_types
    }

def make_cellphonedb_metadata(biomarker):
    return {
        'related_biomarkers': [biomarker.partner_a]
    }


def update_labels(symbol_dict, biomarker):
    if biomarker.label not in symbol_dict[biomarker.symbol].label:
        symbol_dict[biomarker.symbol].label.append(biomarker.label) #TODO: make label updating/checking more comprehensive
    return symbol_dict

def update_metadata(symbol_dict, biomarker, marker_source):
    metadata_func_map = {
        MarkerSource.CELLXGENE_CANONICAL: make_cellxgene_canonical_metadata,
        MarkerSource.CELLXGENE_COMPUTATIONAL: make_cellxgene_computational_metadata,
        MarkerSource.HUBMAP: make_hubmap_metadata,
        MarkerSource.CELLPHONEDB: make_cellphonedb_metadata
    }

    metadata_func = metadata_func_map.get(marker_source)
    if metadata_func:
        possible_metadata = metadata_func(biomarker)

        if marker_source is not MarkerSource.CELLPHONEDB:
            if possible_metadata not in getattr(symbol_dict[biomarker.symbol].additionalMetadata, marker_source.value.lower(), []):
                getattr(symbol_dict[biomarker.symbol].additionalMetadata, marker_source.value.lower()).append(possible_metadata)
        else:
            if possible_metadata not in getattr(symbol_dict[biomarker.symbol].additionalMetadata, marker_source.value.lower(), {}):
                symbol_dict[biomarker.symbol].additionalMetadata.marker_source.value = (possible_metadata)


    return symbol_dict

def combine_biomarkers(canonical_marker_genes,computational_marker_genes, hubmap_biomarkers):
    merged_by_symbol = {}

    for marker in canonical_marker_genes:
        symbol = marker.symbol
        if symbol not in merged_by_symbol:
            #  new merged biomarker entry
            merged_by_symbol[symbol] = BioQueryBiomarker(
                symbol=symbol,
                source=["cellxgene"],
                type="marker_gene",
                uniprot_id="",
                label=[marker.label],
                additionalMetadata=AdditionalMetadata(
                    cellxgene_canonical=[{
                        'tissue': marker.tissue,
                        'publication': marker.publication,
                        'publication_titles': marker.publication_titles
                    }],
                    cellxgene_computational=[],
                    hubmap=[],  # Add hubmap data if available
                    cellphonedb={},
                ),
            )
        else:
            # update metadata
            merged_by_symbol = update_labels(merged_by_symbol, marker)
            merged_by_symbol = update_metadata(merged_by_symbol, marker, MarkerSource.CELLXGENE_CANONICAL)
            
    for marker in computational_marker_genes:
        symbol = marker.symbol

        if symbol not in merged_by_symbol:
            # Create new biomarker entry
            merged_by_symbol[symbol] = BioQueryBiomarker(
                symbol=symbol,
                source=["cellxgene"],
                type="marker_gene",
                uniprot_id="",
                label=[marker.label],
                additionalMetadata=AdditionalMetadata(
                    cellxgene_canonical=None,
                    cellxgene_computational=[{
                        'me': marker.me,
                        'pc': marker.pc,
                        'marker_score': marker.marker_score,
                        'specificity': marker.specificity,
                        'gene_ontology_term_id': marker.gene_ontology_term_id,
                        'groupby_dims': marker.groupby_dims
                    }],
                    hubmap=[],  # Add hubmap data if available
                    cellphonedb={},
                ),
            )
        else:
            # update biomarker entry
            merged_by_symbol = update_labels(merged_by_symbol, marker)
            merged_by_symbol = update_metadata(merged_by_symbol, marker, MarkerSource.CELLXGENE_COMPUTATIONAL)

    for marker in hubmap_biomarkers:
        symbol = marker.symbol
        if symbol not in merged_by_symbol:
            # Create new biomarker entry
            merged_by_symbol[symbol] =  BioQueryBiomarker(
                symbol=symbol,
                source=["hubmap"],
                type="marker_gene",
                uniprot_id="",
                label=[marker.label],
                additionalMetadata=AdditionalMetadata(
                    cellxgene_canonical= [],
                    cellxgene_computational=[],
                    hubmap=[{
                        'anatomical_structures': marker.anatomical_structures,
                        'cell_types': marker.cell_types
                    }],
                    cellphonedb={}, 
                )
            )
        else:
            merged_by_symbol = update_labels(merged_by_symbol, marker)
            merged_by_symbol = update_metadata(merged_by_symbol, marker, MarkerSource.HUBMAP)
            if "hubmap" not in merged_by_symbol[symbol].source:
                merged_by_symbol[symbol].source.append("hubmap")

                                     
    # convert to list of biomarker dictionaries
   
    merged_biomarkers = list(merged_by_symbol.values())
    biomarker_list = []
    for bm in merged_biomarkers:
        biomarker_list.append(bm.dict())
        
    return biomarker_list
    

def process_biomarkers(canonical_info, data_driven_info, hubmap_biomarkers):
    
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
    

# =====================================
# Search function
# =====================================

def process_output_into_list_symbols(final_biomarkers_list):
    symbols_list = []
    for biomarker in final_biomarkers_list:
        symbols_list.append(biomarker['symbol'])
    return symbols_list



def process_input(input):
    print(input)
    processed_string = input.replace('_', ' ')
    processed_string = str(processed_string)
    print(processed_string)
    return processed_string


def search(input_string, output_type): #method to call full search
    file_to_use = data_file
    cell_type = process_input(input_string)
    print("did this")
    hubmap_results = get_hubmap_biomarkers(cell_type, file_to_use)
    hubmap_markers = get_hubmap_biomarker_objects(cell_type, hubmap_results)
    canonical_markers, data_driven_markers = get_biomarkers_from_cl(process_input(input_string))
    
    processed_data = process_biomarkers(canonical_markers, data_driven_markers, hubmap_markers)
    if output_type == OutputType.LIST_SYMBOLS:
        final_list = process_output_into_list_symbols(processed_data)
        output_file = './outputs/av_phage.txt'
        with open(output_file, 'w') as file:
            for gene in final_list:
                file.write(gene + '\n')
        return final_list
    else:
        json_string = json.dumps(processed_data, indent=4)
        output_file = "./outputs/av_macrophage.json"

        with open(output_file, 'w') as f:
            f.write(json_string)

        print(f"c2bm {output_file}")
   

search(sys.argv[1], OutputType.LIST_SYMBOLS)