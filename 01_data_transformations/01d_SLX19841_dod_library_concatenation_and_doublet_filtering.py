import sys, os, re, warnings, getpass
from pathlib import Path
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr

import matplotlib
import matplotlib.pyplot as plt


# Reminder: Always verify paths for semantic intent.

home = str(Path.home())
user = getpass.getuser()
basedir = os.path.join(home, 'databoard/users', user, 'AEH/')

sc.settings.writedir = os.path.join(basedir, 'SLX19841/analysis/h5ad/')

# For sequencing run SLX19841 DOD samples

samples = ['SIGAA12', 'SIGAB12', 'SIGAC12', 'SIGAD12']
libraries = {}

processed_path = os.path.join( basedir, 'Nicole/2020/SLX19841/processed/' )

for sample in samples:
    libraries[sample] = sc.read_10x_h5(os.path.join(processed_path, sample, 'outs', 
                                                    'filtered_feature_bc_matrix.h5'), 'hg19')
    libraries[sample].var_names_make_unique()
    libraries[sample].obs['library'] = sample


libraries['SIGAA12'].obs['donor'] = 'KSP29'
libraries['SIGAB12'].obs['donor'] = 'KSP32'
libraries['SIGAC12'].obs['donor'] = 'KSP29'
libraries['SIGAD12'].obs['donor'] = 'KSP32'


libraries['SIGAA12'].obs['organ'] = 'BM'
libraries['SIGAB12'].obs['organ'] = 'BM'
libraries['SIGAC12'].obs['organ'] = 'SPL'
libraries['SIGAD12'].obs['organ'] = 'SPL'
    

data = libraries['SIGAA12'].concatenate( [ libraries['SIGAB12'], libraries['SIGAC12'], libraries['SIGAD12'] ] )


# Predict and remove putative doublets

scrub = {}
doublet_scores = {}
predicted_doublets = {}

for sample in samples:
    print(sample)
    scrub[sample] = scr.Scrublet( data[np.array(data.obs['library'] == sample), :].X )

    doublet_scores[sample], predicted_doublets[sample] = scrub[sample].scrub_doublets()
    print('\n\n')
    
for sample in samples:
    print(sample, ':', sum(predicted_doublets[sample]))    
    
    
sample = 'SIGAA12'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.23)
sum(predicted_doublets[sample])
    
    
sample = 'SIGAB12'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.23)
sum(predicted_doublets[sample])


sample = 'SIGAC12'
scrub[sample].plot_histogram()
#
#

sample = 'SIGAD12'
scrub[sample].plot_histogram()
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=0.4)
sum(predicted_doublets[sample])



data_doublets = os.path.join( sc.settings.writedir , '..', 'doublets')
if not os.path.exists( data_doublets ):
    os.makedirs( data_doublets )

for key in doublet_scores:
    np.savetxt(os.path.join(data_doublets, key+'_doublet_scores.txt'), doublet_scores[key])

doublet_scores_list = []
for key in doublet_scores:
    #print(key)
    doublet_scores_list += list(doublet_scores[key])
    
data.obs['doublet_scores'] = doublet_scores_list
len(doublet_scores_list)


predicted_doublets_mask = []
for key in predicted_doublets:
    #print(key)
    predicted_doublets_mask += list(predicted_doublets[key])
len(predicted_doublets_mask)


predicted_singletons_mask = [not i for i in predicted_doublets_mask]
data = data[np.array(predicted_singletons_mask), :].copy()

print('Removing %d cells due to doublet scoring' %(len(predicted_singletons_mask) - sum(predicted_singletons_mask)))


sc.write('SLX19841_DOD_filtered_gene_bc_expression_minus_putative_doublets', data)
