import sys, os, re, warnings, getpass
from pathlib import Path
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

import matplotlib
import matplotlib.pyplot as plt


# For sequencing runs SLX12978 and SLX14831

# Reminder: Always verify paths for semantic intent.

home = str(Path.home())
user = getpass.getuser()
basedir = os.path.join(home, 'datafloor/users', user, 'AEH/')

sc.settings.writedir = os.path.join(basedir, 'data_intermediate/SLX_14831_12978_COMBO/')

data_in_1 = os.path.join(basedir, 'data_in/SLX_12978')
data_in_2 = os.path.join(basedir, 'data_in/SLX_14831')

# DOD1
sigad9 = sc.read_10x_h5(os.path.join(data_in_1, 'D9_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigae9 = sc.read_10x_h5(os.path.join(data_in_1, 'E9_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigaf9 = sc.read_10x_h5(os.path.join(data_in_1, 'F9_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigag9 = sc.read_10x_h5(os.path.join(data_in_1, 'G9_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigah9 = sc.read_10x_h5(os.path.join(data_in_1, 'H9_filtered_gene_bc_matrices_h5.h5'), 'hg19')

# DOD2
sigab10 = sc.read_10x_h5(os.path.join(data_in_2, 'SIGAB10_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigac10 = sc.read_10x_h5(os.path.join(data_in_2, 'SIGAC10_filtered_gene_bc_matrices_h5.h5'), 'hg19')
sigad10 = sc.read_10x_h5(os.path.join(data_in_2, 'SIGAD10_filtered_gene_bc_matrices_h5.h5'), 'hg19')


sigad9.var_names_make_unique()
sigae9.var_names_make_unique()
sigaf9.var_names_make_unique()
sigag9.var_names_make_unique()
sigah9.var_names_make_unique()

sigab10.var_names_make_unique()
sigac10.var_names_make_unique()
sigad10.var_names_make_unique()


sigad9.obs['donor'] = 'DOD1'
sigae9.obs['donor'] = 'DOD1'
sigaf9.obs['donor'] = 'DOD1'
sigag9.obs['donor'] = 'DOD1'
sigah9.obs['donor'] = 'DOD1'

sigab10.obs['donor'] = 'DOD2'
sigac10.obs['donor'] = 'DOD2'
sigad10.obs['donor'] = 'DOD2'


sigad9.obs['library'] = 'SIGAD9'
sigae9.obs['library'] = 'SIGAE9'
sigaf9.obs['library'] = 'SIGAF9'
sigag9.obs['library'] = 'SIGAG9'
sigah9.obs['library'] = 'SIGAH9'

sigab10.obs['library'] = 'SIGAB10'
sigac10.obs['library'] = 'SIGAC10'
sigad10.obs['library'] = 'SIGAD10'


sigad9.obs['organ'] = 'PB'
sigae9.obs['organ'] = 'BM'
sigaf9.obs['organ'] = 'BM'
sigag9.obs['organ'] = 'SPL'
sigah9.obs['organ'] = 'SPL'

sigab10.obs['organ'] = 'BM'
sigac10.obs['organ'] = 'SPL'
sigad10.obs['organ'] = 'PB'


data = sigad9.concatenate([sigae9, sigaf9, sigag9, sigah9, sigab10, sigac10, sigad10])


# Predict and remove putative doublets

samples = list( data.obs['library'].unique() )
samples

scrub = {}
doublet_scores = {}
predicted_doublets = {}

# approximate expected multiplet rate
# from page 6 by '# of Cells Recovered'
emr = {'SIGAD9': 0.020,
       'SIGAE9': 0.022,
       'SIGAF9': 0.021,
       'SIGAG9': 0.019,
       'SIGAH9': 0.021,
       'SIGAB10': 0.048, 
       'SIGAC10': 0.039, 
       'SIGAD10': 0.050}


for sample in samples:
    print(sample)
    scrub[sample] = scr.Scrublet( data[np.array(data.obs['library'] == sample), :].X, 
                                  expected_doublet_rate=emr[sample] 
                                )
    doublet_scores[sample], predicted_doublets[sample] = scrub[sample].scrub_doublets()
    print('\n\n')


# adjust cutoff manually based on the histogram

sample = 'SIGAD9'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.19)

sample = 'SIGAE9'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.14)

sample = 'SIGAF9'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.15)

sample = 'SIGAG9'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.15)

sample = 'SIGAH9'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.13)

sample = 'SIGAB10'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.20)

sample = 'SIGAC10'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.20)

sample = 'SIGAD10'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.25)



for sample in samples:
    scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))
    scrub[sample].plot_embedding('UMAP', order_points=True)


data_doublets = os.path.join( sc.settings.writedir , 'doublets')

if not os.path.exists( data_doublets ):
    os.makedirs( data_doublets )
    
for key in doublet_scores:
    np.savetxt(os.path.join(data_doublets, key+'_doublet_scores.txt'), doublet_scores[key])
    
    
doublet_scores_list = []
for key in doublet_scores:
    doublet_scores_list += list(doublet_scores[key])
    
data.obs['doublet_scores'] = doublet_scores_list

# create the boolean mask to filter out predicted doublets
predicted_doublets_mask = []

for key in predicted_doublets:
    predicted_doublets_mask += list(predicted_doublets[key])
    
predicted_singletons_mask = [not i for i in predicted_doublets_mask]


data = data[np.array(predicted_singletons_mask), :].copy()
print('Removing %d cells due to doublet scoring' %(len(predicted_singletons_mask) - sum(predicted_singletons_mask)))
    

sc.write('SLX14831_12978_filtered_gene_bc_expression_minus_putative_doublets', data)
