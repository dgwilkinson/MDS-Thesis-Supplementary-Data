#!/usr/bin/env python
# coding: utf-8

# In[91]:


import scanpy as sc, scipy, os, sys, anndata, pandas as pd, re


# In[ ]:


# Read Command Line arguments
# Format: res_path, out_path, sample_id
sample_id = sys.argv[-1]
res_path, out_path = map(os.path.expanduser, sys.argv[1:-1])


# In[ ]:


# I accidentally erased the results from ACH0010 so I skip processing this sample 
if sample_id == "ACH0010":
    quit()


# In[69]:


# Load in Visium Results object
visium_res = sc.read_h5ad(res_path)


# In[93]:


# Extract the means model output from the metadata
ct_mean_abd = visium_res.obsm["means_cell_abundance_w_sf"].copy()
ct_mean_abd.rename(columns = lambda name: re.sub("meanscell_abundance_w_sf_", "", name),
                  inplace = True)


# In[ ]:


# Save mean abundances to csv
ct_mean_abd.to_csv(out_path + sample_id + "_mean_abd.csv")

