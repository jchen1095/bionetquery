import anndata as an 
import cellphonedb 
from cellphonedb.src.core.methods import cpdb_analysis_method


raw_anndata = an.read_h5ad("bionetquery/secondary_analysis.h5ad")
df_anndata = an.AnnData.to_df(raw_anndata)
print(raw_anndata)
print(df_anndata) #where obs_names initializes the index, and var_names the columns.


