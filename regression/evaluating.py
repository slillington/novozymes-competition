#%%
import pandas as pd
import numpy as np
import os
os.getcwd()

# get the directory of the file
dir_path = os.path.dirname(os.path.abspath(__file__))
os.chdir(dir_path)

# list the files in the data folder
cleaned_data_path = os.path.join(dir_path, 'datasets_cleaned')

# get a list of the files in the data folder
cleaned_data_files = os.listdir(cleaned_data_path)

# sort the feature files list by the centroid number
cleaned_data_files.sort(key=lambda x: int(x.split('_')[2]))

# list the centroid numbers in the feature files list
cluster_centroids = [int(x.split('_')[2]) for x in cleaned_data_files]

# get a list of df's from each cluster
dataframes = []
for file in cleaned_data_files:
    if file.endswith('.csv'):
        dataframes.append(pd.read_csv(f'{cleaned_data_path}/{file}', index_col=0))
    
#%%
# Make a 3x5 grid of plots
import matplotlib.pyplot as plt

for column in dataframes[0].columns:
    fig, axes = plt.subplots(3, 5, figsize=(15, 7),sharex=True,sharey=True)
    axes = axes.ravel()
    # Loop through each cluster and plot histograms of pH
    for i, df in enumerate(dataframes):
        # plot in the correct subplot
        axes[i].hist(df[column], bins=15)
        axes[i].set_title(f'Cluster {cluster_centroids[i]}')
        axes[i].set_xlabel(column)
    plt.show()

#%% Plot the 


# %%
