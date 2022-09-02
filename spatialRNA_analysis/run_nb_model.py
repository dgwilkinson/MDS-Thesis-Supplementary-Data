#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Import relevant modules
import sys, scanpy as sc, anndata, numpy as np, pandas as pd, os
import cell2location, scvi, matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import torch
print(torch.get_num_threads())

# In[ ]:


# Read Command Line Arguments
# Format: ref_obj, output_folder, sample_id
ref_obj, out_path, sample_id = sys.argv[1:]
ref_obj = os.path.expanduser(ref_obj)
out_path = os.path.expanduser(out_path)


# In[18]:


# Read annotated integrated snRNA-seq h5ad into scanpy AnnData object
rna_data = sc.read(ref_obj)
# Raw slot redundant to X
del rna_data.raw


# In[19]:


# Filter reference data (scanpy PreProcessing)
print(rna_data.shape)
# Filter cells and genes for entirely empty entries
sc.pp.filter_cells(rna_data, min_genes = 1)
sc.pp.filter_genes(rna_data, min_cells = 1)

# Separate MT genes from GCM
rna_data.var["MT_feat"] = [feat.startswith("MT-") for feat in rna_data.var_names]
rna_data.obsm["MT"] = rna_data[:, rna_data.var["MT_feat"].values].X.toarray()
rna_data = rna_data[:, ~rna_data.var["MT_feat"].values]

# Calculate feature-wise qc statistics
rna_data.var["n_cells"] = (rna_data.X.toarray() > 0).sum(axis = 0)
rna_data.var["mean_nzero"] = (rna_data.X.toarray().sum(axis = 0))/rna_data.var["n_cells"]
rna_data.var["SYMBOL"] = rna_data.var_names

# Perform "permisive gene selection" for widely/highly expressed genes
mean_gene_thresh = np.exp(5e-2)
cell_count_thresh1 = np.log10(rna_data.n_obs * 5e-4)
cell_perc_thresh2 = 3e-2

from cell2location.utils.filtering import filter_genes
permitted_genes = filter_genes(rna_data,
                               cell_count_cutoff = cell_count_thresh1,
                               cell_percentage_cutoff2 = cell_perc_thresh2,
                               nonz_mean_cutoff = mean_gene_thresh)

rna_data = rna_data[:, permitted_genes]

print(rna_data.shape)


# In[22]:


# Prepare reference data for regression model
cell2location.models.RegressionModel.setup_anndata(adata = rna_data,
                                                  batch_key = "orig.ident",
                                                  labels_key = "cell_type")

NB_model = cell2location.models.RegressionModel(rna_data)

# View model parameters to verify correct setup
NB_model.view_anndata_setup()


# In[45]:


NB_model.train(max_epochs = 100, batch_size = 1024, lr = 0.01)


# In[49]:


# Sample the Bayesian Posterior Distribution of the NB Regression to 
# generate cell type specific feature statistics
rna_data = NB_model.export_posterior(rna_data, 
                                     sample_kwargs = {"num_samples":1000,
                                                      "batch_size":1024})

# Results are stored in varm
# rna_data.varm["means_per_cluster_mu_fg"]


# In[60]:


# Save NB model and AnnData object containing results
NB_model.save(out_path+"NB_models",prefix = sample_id + "_NB_", overwrite = True)
rna_data.write(out_path+sample_id+"_NB_results.h5ad")

