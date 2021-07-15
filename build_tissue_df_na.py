import pandas as pd
import numpy as np
import sys

tissue_type = sys.argv[1]

with open('rm_probe.txt', 'r') as f:
    line = f.readline()
probes_to_rm = line.split("|")
probes_to_rm = [int(i)-2 for i in probes_to_rm]

def build_data_matrix(tissue_type):
    with open('data_TCGA/dir_'+tissue_type+'.txt', 'r') as f:
        first = True
        for file in f:
            data = pd.read_csv('data_TCGA/'+file.strip(), sep="\t")
            data_rm = data.drop(probes_to_rm)
            beta = data_rm['Beta_value'].values
            if first:
                first = False
                data_matrix = beta
            else:
                data_matrix = np.vstack((data_matrix, beta))
                
    probe_names = data_rm['Composite Element REF'].values
    data_df = pd.DataFrame(data_matrix.T, index=probe_names)
    data_df_na = data_df.dropna(thresh=int(data_df.shape[1]*0.95), axis=0)
    data_df_na.to_csv(tissue_type+"_beta.zip", index=True, header=False, compression="zip")

build_data_matrix(tissue_type)
