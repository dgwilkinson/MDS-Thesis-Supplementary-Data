#!/usr/bin/env python
# coding: utf-8

# In[3]:


import torch, os
print(torch.get_num_threads())
print(os.environ.get("OMP_NUM_THREADS"))
print(os.environ.get("MKL_NUM_THREADS"))
print(os.environ.get("OPENBLAS_NUM_THREADS"))
print(os.environ.get("BLIS_NUM_THREADS"))
