###
### This is for parsing results from blast or hmmr search.
### Author: Huiluo Cao
### Unversity of Hong Kong
### October, 2017

import pandas as pd
import csv, re
import os
import numpy as np
##import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

####clean hmmr gene2type list###
hmmr_df = pd.read_table('C:/Users/Huiluo/Desktop/HKU-microbiology/microbiome/ref_seqs/HMP/hmmr_gff3/hmmr-db_list.txt', sep = '\t', skiprows = 1, engine='python', header = None)
hmmr_ndf = hmmr_df[[0,1,4,8,10]]
hmmr_ndf.columns = ['ref_gene', 'short_name', 'genesymbol', 'product_name','resistance_type']
print(pd.DataFrame(hmmr_ndf['resistance_type'].unique(), columns = ['resistance_type']))
hmmr_arg_type = pd.DataFrame(hmmr_ndf['resistance_type'].unique(), columns = ['resistance_type'])

files = os.listdir(os.curdir)
for file in files:
    if file.endswith("hmmr_arg_type"):
        sample_id = file.split("_")[0]
        print(sample_id)
        if sample_id.split(".")[1] == "an":
            df = pd.read_csv(file, sep = '\t',  header = 0)
            hmmr_ndf2_type = df.groupby('resistance_type')['resistance_type'].count().reset_index(name = '%s'%sample_id)
            print(hmmr_ndf2_type)
            hmmr_arg_type = pd.merge(hmmr_arg_type, hmmr_ndf2_type, how = 'left', on = ['resistance_type'])
hmmr_arg_type = hmmr_arg_type.fillna(0)
hmmr_arg_type.iloc[:,1:] = hmmr_arg_type.iloc[:,1:].astype(int)
print(hmmr_arg_type, hmmr_arg_type['resistance_type'], hmmr_arg_type.ix[:,1])
with open(os.path.join(os.curdir,"an_hmmr_arg_type"),'w') as g:
        hmmr_arg_type.to_csv(g, header = True, index = False, sep = "\t")

plt.bar(np.arange(len(hmmr_arg_type['resistance_type'])), hmmr_arg_type.iloc[:,2])
plt.title("Abundance of ARG types of human microbiota")
plt.xlabel("Sample_ID")
plt.ylabel("Number of ARGs")
plt.show()
