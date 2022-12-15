#%%
import pandas as pd
import numpy as np
import os

# get the path to the directory that this file is in
dir_path = os.path.dirname(os.path.realpath(__file__))
# get a path to the without_energies directory
without_energies_dir = os.path.join(dir_path, "without_energies")
# list all the csv files in the without_energies directory
cluster_csvs = [f for f in os.listdir(without_energies_dir) if f.endswith('.csv')]
# list the numbers of the csv files
cluster_nums = [int(f.split('_')[1]) for f in cluster_csvs]

#%%
# get a path to the energies_from_ambertools directory
ambertools_dir = os.path.join(dir_path, "../energies_from_ambertools")
# open the energies_from_ambertools.csv file
energies_df = pd.read_csv(os.path.join(ambertools_dir, "energies_training.csv"),index_col=0)
# remove the '.pdb' from the index of the energies_df and turn to int
energies_df.index = [int(f.split('.')[0]) for f in energies_df.index]

# %%
# Let's look at cluster_972_features.csv first
for cluster_csv in cluster_csvs:
    cluster_df = pd.read_csv(os.path.join(without_energies_dir, cluster_csv),index_col=0)

    # concatenate the energies_df and cluster_df horizontally if they have the same index
    # Create a list of the two DataFrames that you want to concatenate
    dfs = [cluster_df, energies_df]

    # Find the intersection of the indices of the two DataFrames
    intersection = pd.Index.intersection(cluster_df.index, energies_df.index)

    # Subset both DataFrames so that they only include rows with indices that are in the intersection
    dfs_subsetted = [df[df.index.isin(intersection)] for df in dfs]

    # Horizontally concatenate the subsetted DataFrames
    cluster_w_energies_df = pd.concat(dfs_subsetted, axis=1)

    # add .pdb to the end of the index again
    cluster_w_energies_df.index = [str(f) + '.pdb' for f in cluster_w_energies_df.index]

    # Ouptut csv file to the with_energies directory but format the last 8 columns in scientific notation with 3 decimal points
    cluster_w_energies_df.to_csv(os.path.join(dir_path, "with_energies", cluster_csv))
