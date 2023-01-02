#%%
import numpy as np
import pandas as pd
import os
import sys

# get the directory of the file
dir_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(dir_path)

# list the files in the data folder
data_path = os.path.join(dir_path, 'batches_210')

# get a list of the files in the data folder
data_files = os.listdir(data_path)

# sort the list by the number after the last underscore but before the .csv
data_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

#%%
# get a list of df's from each csv
dataframes = []
for file in data_files:
    if file.endswith('.csv'):
        new_df = pd.read_csv(f'{data_path}/{file}', index_col=0)
        # drop the rows with all 0s
        new_df = new_df.drop(new_df[new_df['e_tot'] == 0].index)
        dataframes.append(new_df)

#%%
# combine all the dataframes into one
df = pd.concat(dataframes)

# remove the .pdb in the index column and format as int
df.index = df.index.str.replace('.pdb', '').astype(int)

# order the rows by the number in the index
df = df.sort_index()
#%%

# delete the duplicate indices
df = df[~df.index.duplicated(keep='first')]

# %% Copy the dataframe to a csv file
df.to_csv('energies_training.csv')

# %%
# determine which indices were skipped
# get the indices of the dataframe
indices = df.index
# get the range of indices
indices_range = np.arange(0, 31389+1)

# which indices are missing from incies_range?
indices_diff = np.setdiff1d(indices_range, indices)


# %%
