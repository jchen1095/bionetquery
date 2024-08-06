from c2bm import search, OutputType
from related_biomarkers import main, MarkerProperty
from network import search_kegg

def test_podocyte_c2bm():
    getting_biomarkers = search("T cell", OutputType.LIST_SYMBOLS)
    print(getting_biomarkers)
    related_biomarkers = main(getting_biomarkers, MarkerProperty.HGNC_SYMBOL, OutputType.LIST_UNIPROT)
    print(related_biomarkers)
    # search_kegg(related_biomarkers)


test_podocyte_c2bm()