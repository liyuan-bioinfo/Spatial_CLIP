# Deal with the non-unique features
sp_adata = sc.read_h5ad(f"./features/{slide_name}/{slide_name}_filter.h5ad")
sp_adata.var.index = sp_adata.var["gene_name"]
sp_adata.var_names = sp_adata.var_names.astype(str)

sp_adata.var_names_make_unique()  
