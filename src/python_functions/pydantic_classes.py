from pydantic import BaseModel, TypeAdapter
from typing import List,Optional, Dict, Union


class AdditionalMetadata(BaseModel):
    """More specific properties from different sources"""
    cellphonedb: Optional[list]
    cellxgene_canonical: Optional[list]
    cellxgene_computational: Optional[list]
    hubmap: Optional[list]

class BioQueryBiomarker(BaseModel):
    """A biomarker object for bioquery"""
    source: list
    symbol: str
    uniprot_id: str
    type:str
    label: list
    additionalMetadata: AdditionalMetadata


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

class CellPhoneDBBiomarker(BaseModel):
    source: str
    partner_a: str
    partner_b: str
    uniprot_id_a: str
    uniprot_id_b: str
    