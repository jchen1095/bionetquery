# import pytest

import sys
sys.path.append('src')

from python_functions.c2bm import search, OutputType
from python_functions.related_biomarkers import main, MarkerProperty
from python_functions.network import search_kegg

def test_podocyte_c2bm():
    getting_biomarkers = search("podocyte", OutputType.LIST_SYMBOLS)
    print(getting_biomarkers)
    related_biomarkers = main(getting_biomarkers, MarkerProperty.HGNC_SYMBOL, OutputType.LIST_UNIPROT)
    print(related_biomarkers)
    search_kegg(related_biomarkers)


#I'm not sure how you would be able to do ground truth for a lot of examples it would require manually looking through all the datasets
#TODO: Make a test search/process biomarkers so that you can try different cell types
