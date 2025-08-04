#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# ---

# ## Code to load the data and detect doublets

# In[1]:


import scrublet as scr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy



def MovePlots(plotpattern, subplotdir):
    os.system('mkdir -p '+str(sc.settings.figdir)+'/'+subplotdir)
    os.system('mv '+str(sc.settings.figdir)+'/*'+plotpattern+'** '+str(sc.settings.figdir)+'/'+subplotdir)


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.settings.figdir = '../results/images/preprocessing/'
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures

sys.executable


# In[3]:


# Benjamini-Hochberg and Bonferroni FDR helper functions.

def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.
    
    Input:
        * pvals - vector of p-values to correct
    """
    pvalues = np.array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def bonf(pvalues):
    """
    Computes the Bonferroni FDR correction.
    
    Input:
        * pvals - vector of p-values to correct
    """
    new_pvalues = np.array(pvalues) * len(pvalues)
    new_pvalues[new_pvalues>1] = 1
    return new_pvalues


# ## Scrublet

# In[4]:



def runScrublet(samples, data_dir):

    for sample in reversed(list(samples)):
        print(sample)
        #import data
        adata_sample = sc.read_10x_mtx(data_dir+sample,var_names='gene_symbols',cache=True) #reading the data
        adata_sample.var_names_make_unique()
        
        #rename cells to SAMPLE_BARCODE
        adata_sample.obs_names = [sample+'_'+i for i in adata_sample.obs_names]
        
        #do some early filtering to retain meaningful cells for doublet inspection
        sc.pp.filter_cells(adata_sample, min_genes=200)
        sc.pp.filter_genes(adata_sample, min_cells=3)

        #convert to lower to be species agnostic: human mito start with MT-, mouse with mt-
        mito_genes = [name for name in adata_sample.var_names if name.lower().startswith('mt-')]
        # for each cell compute fraction of counts in mito genes vs. all genes
        # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
        #adata_sample.obs['percent_mito'] = np.sum(
        #    adata_sample[:, mito_genes].X, axis=1).A1 / np.sum(adata_sample.X, axis=1).A1
        #adata_sample = adata_sample[adata_sample.obs['percent_mito'] < 0.2, :]

        #set up and run Scrublet, seeding for replicability
        np.random.seed(0)
        scrub = scr.Scrublet(adata_sample.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
        adata_sample.obs['scrublet_score'] = doublet_scores

        #overcluster prep. run turbo basic scanpy pipeline
        sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
        sc.pp.log1p(adata_sample)
        sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
        sc.pp.scale(adata_sample, max_value=10)
        sc.tl.pca(adata_sample, svd_solver='arpack')
        sc.pp.neighbors(adata_sample)
        #overclustering proper - do basic clustering first, then cluster each cluster
        sc.tl.leiden(adata_sample)
        adata_sample.obs['leiden'] = [str(i) for i in adata_sample.obs['leiden']]
        for clus in np.unique(adata_sample.obs['leiden']):
            adata_sub = adata_sample[adata_sample.obs['leiden']==clus].copy()
            sc.tl.leiden(adata_sub)
            adata_sub.obs['leiden'] = [clus+','+i for i in adata_sub.obs['leiden']]
            adata_sample.obs.loc[adata_sub.obs_names,'leiden'] = adata_sub.obs['leiden']

        #compute the cluster scores - the median of Scrublet scores per overclustered cluster
        for clus in np.unique(adata_sample.obs['leiden']):
            adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_cluster_score'] =                 np.median(adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_score'])
        #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
        med = np.median(adata_sample.obs['scrublet_cluster_score'])
        mask = adata_sample.obs['scrublet_cluster_score']>med
        mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
        #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
        zscores = (adata_sample.obs['scrublet_cluster_score'].values - med) / (1.4826 * mad)
        adata_sample.obs['zscore'] = zscores
        pvals = 1-scipy.stats.norm.cdf(zscores)
        adata_sample.obs['bh_pval'] = bh(pvals)
        adata_sample.obs['bonf_pval'] = bonf(pvals)

        #create results data frame for single sample and copy stuff over from the adata object
        scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
        for score in scorenames:
            scrublet_sample[score] = adata_sample.obs[score]

        #write out complete sample scores
        scrublet_sample.to_csv(outputDir+sample+'.csv')


# In[ ]:





# ### Reading Metadata

# In[5]:


#Toxoplasma infected metadata
metaTg = pd.read_csv('../metadata/meta_exp_infection_Tg_scell.csv',index_col=0)
metaTg['donor'] = metaTg['donor'].astype('str')
print('Number of samples Tg: ', metaTg.index.size)


# In[6]:


metaTg


# In[ ]:





# ## Reading and saving the files as .h5 for the pipeline
# Please remember that STARsolo v2.7.9a returns the file genes.tsv.gz, however scanpy now only reads the file features.tsv.gz.
# Just change the name of the file

# In[7]:


#there's loads of clustering going on, so set verbosity low unless you enjoy walls of text
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)

scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']







# #### Toxoplasma

# In[8]:


#output directory
out_dir = '../results/'


if not os.path.exists(out_dir+'scrublet-scores'):
    os.makedirs(out_dir+'scrublet-scores')
    #loop over the subfolders of the rawdata folder

outputDir=out_dir+'scrublet-scores/'


# In[15]:


runScrublet(metaTg.index.to_list(), data_dir='../data/')


# In[ ]:




