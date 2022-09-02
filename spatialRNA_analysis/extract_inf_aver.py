#!/usr/bin/env python
# coding: utf-8

# In[15]:


import scanpy as sc, anndata, pandas as pd, numpy as np
import os, sys, re


# In[ ]:


# Read in Command line arguments 
# Format: results_path, sample_id, out_path
results_path, sample_id, out_path = sys.argv[1:]


# In[3]:


results = os.path.expanduser(results_path + sample_id + "_NB_results.h5ad")


# In[5]:


# Read in NB results AnnData object
rna_data = sc.read_h5ad(results)


# In[23]:


# Extract inferred feature averages across cell types 
metric_name = "means_per_cluster_mu_fg"
cell_types = rna_data.uns["mod"]["factor_names"]
inf_aver = rna_data.varm[metric_name].copy()

# Rename columns with cell type labels 
inf_aver.rename(columns = lambda name: re.sub(metric_name + "_",
                                              "", name),
               inplace = True)


# In[31]:


# Verify the retention of all cell types
print(f"All Cell Types Retained?: {inf_aver.columns.symmetric_difference(cell_types).empty}")


# In[ ]:


# Write Results to csv for rapid loading in next step
inf_aver.to_csv(out_path + sample_id + "_inf_aver.csv")

