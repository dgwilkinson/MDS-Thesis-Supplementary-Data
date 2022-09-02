#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc, cell2location, numpy as np, pandas as pd
import sys, os, anndata


# In[ ]:


# Read Command Line Arguments
# Format: sample_path, inf_aver_path, fig_path
sample_id = sys.argv[-1]
sample_path, inf_aver_path, out_path = list(map(os.path.expanduser, sys.argv[1:-1]))


# In[21]:


# Read in inf_aver csv and processed visium sample
inf_aver = pd.read_csv(inf_aver_path, index_col = 0)
# Note Ribosomal and MT genes have already been filtered
visium_obj = sc.read_h5ad(sample_path)


# In[28]:


# Identify common features between inf_aver and visium_obj
common_feats = inf_aver.index.intersection(visium_obj.var_names)
# Subset to intersection
visium_obj = visium_obj[:, common_feats].copy()
inf_aver = inf_aver.loc[common_feats, :].copy()


# In[31]:


# Initialise cell2location spatial mapping model
cell2location.models.Cell2location.setup_anndata(adata = visium_obj,
                                                batch_key = "orig.ident",
                                                labels_key = "opt_clust")
mapping_model = cell2location.models.Cell2location(visium_obj,
                                                  cell_state_df = inf_aver,
                                                  N_cells_per_location = 8,
                                                  detection_alpha = 20)
# Verify correct model setup
mapping_model.view_anndata_setup()


# In[36]:


import torch
torch.set_num_threads(1)
print("Number of Threads:", torch.get_num_threads())


# In[37]:


# Run cell2location spatial mapping model
mapping_model.train(max_epochs = 500,
                   batch_size = None,
                   train_size = 1,
                   lr = 0.01)


# In[48]:


# Export Bayesian model posterior estimates on cell abundance to visium_obj
visium_obj = mapping_model.export_posterior(visium_obj,
                                           sample_kwargs = {"num_samples":1000})

# Save modified visium_obj results
visium_obj.write(out_path + sample_id + "_abundances.h5ad")

