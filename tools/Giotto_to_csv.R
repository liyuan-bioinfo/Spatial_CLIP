# transform RData or RDS files to csv
library(Giotto)
library(Matrix)
setwd("temp")

# my_giotto@raw_exprs            my_giotto@spatial_locs         my_giotto@gene_ID              my_giotto@dimension_reduction  my_giotto@instructions
# my_giotto@norm_expr            my_giotto@cell_metadata        my_giotto@spatial_network      my_giotto@nn_network           my_giotto@offset_file
# my_giotto@norm_scaled_expr     my_giotto@gene_metadata        my_giotto@spatial_grid         my_giotto@images               my_giotto@OS_platform
# my_giotto@custom_expr          my_giotto@cell_ID              my_giotto@spatial_enrichment   my_giotto@parameters

file1_name = "bsw_validation_deepCell.RDS"
my_giotto = readRDS(file1_name)

# Save expression matrix
expr_matrix <- as.matrix(my_giotto@raw_exprs)
write.csv(expr_matrix, "bsw_validation_deepCell.csv")

# Save metadata
write.csv(as.data.frame(my_giotto@cell_metadata), "bsw_validation_deepCell_cell_metadata.csv")
write.csv(as.data.frame(my_giotto@gene_metadata), "bsw_validation_deepCell_gene_metadata.csv")
write.csv(as.data.frame(my_giotto@spatial_locs), "bsw_validation_deepCell_spatial_locs.csv")

# under PYTHON env.
# transform .csv to .h5ad
{
  # import scanpy as sc
  # import pandas as pd
  # import os
  # import anndata
  
  # dataset = "bsw_validation_deepCell"
  
  # # Load matrix market file
  # expr = pd.read_csv(f"{dataset}.csv", index_col=0)
  
  # # Load metadata
  # cell_meta = pd.read_csv(f"{dataset}_cell_metadata.csv", index_col=1)
  # gene_meta = pd.read_csv(f"{dataset}_gene_metadata.csv", index_col=1)
  # spatial_loc = pd.read_csv(f"{dataset}_spatial_locs.csv", index_col=3)
  
  # # Create AnnData
  # adata = anndata.AnnData(X=expr.T, obs=cell_meta, var=gene_meta)
  
  # adata.obsm['spatial'] = spatial_loc[['sdimx', 'sdimy']].values
  # adata.write(f"{dataset}_giotto_converted.h5ad")
}





