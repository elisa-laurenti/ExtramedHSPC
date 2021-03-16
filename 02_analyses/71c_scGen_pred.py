import os, sys, json, operator, getpass, pickle
from pathlib import Path
from datetime import datetime

import scgen

import numpy as np
import pandas as pd
import scanpy as sc


home = str(Path.home())
user = getpass.getuser()

basedir = os.path.join(home, 'datafloor/users', user, '2020/SLX19841/')

sc.settings.writedir = os.path.join(basedir, 'analysis/h5ad/')


data = sc.read('COMBO10_NO_SPL3_Seurat3_lognorm')

data = data[:, data.var.highly_variable].copy()


train = data[data.obs.organ.isin(['PB']), : ].copy()

train.obs['condition'] = train.obs.donor.apply(lambda x: 'DOD' if x in ['DOD1', 'DOD2'] else 'LD')


scg = scgen.VAEArithKeras(x_dimension=train.shape[1], 
                          model_path=os.path.join(sc.settings.writedir, '..', 'models/lognorm_dod2ld_pred_hvg') )

scg.train(train_data=train, n_epochs=100)


predict = data[data.obs.organ == 'SPL', :].copy()
predict.obs['condition'] = 'DOD'


pred, delta = scg.predict(adata=train, 
                          adata_to_predict=predict,
                          conditions={"ctrl": "DOD", "stim": "LD"}, 
                          cell_type_key="leiden.1.2", 
                          condition_key="condition")


pred_adata = sc.AnnData(pred, obs={"condition":["pred"]*len(pred)}, var={"var_names":train.var_names})

sc.write('COMBO10_SPL_lognorm_LD_prediction_hvg', pred_adata)

pickle.dump( delta,  open(  os.path.join(sc.settings.writedir, '..', 'models/dod2ld_pred_hvg.delta'), 'wb'  ) )
