"""
Cleaning the features (from the `features` directory)
Sam Lobo
"""
#%%
import pandas as pd
import numpy as np
import os

# Note the path to this file
dir_path = os.path.dirname(os.path.abspath(__file__))
# Output directory
output_dir = os.path.join(dir_path, 'datasets_cleaned')
# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)    

# Note the path of the data folder
raw_data_path = os.path.join(dir_path, '../features/with_energies')
# List the files in the data folder
feature_files = os.listdir(raw_data_path)
# Sort the feature files list by the centroid number
feature_files.sort(key=lambda x: int(x.split('_')[1]))
# Get a list of just the centroid numbers in the feature files list
cluster_centroids = [int(x.split('_')[1]) for x in feature_files]

# Load all the csvs from raw_data_path into a list of dataframes
dataframes = []
for file in feature_files:
    if file.endswith('.csv'):
        dataframes.append(pd.read_csv(f'{raw_data_path}/{file}', index_col=0))

#%%
### SORT BY MELTING POINT, CLEAN THE DATA, AND REOUTPUT THE CSV FILES
for i, df in enumerate(dataframes):
    # Rearrange the rows so that the melting points are in ascending order
    df = df.sort_values(by='tm')

    # If Tm is 25 exactly then delete the row
    df = df[df['tm'] != 25]
    # If pH > 14 then delete the row (novozymes or kaggle or scientists mislabelled it)
    df = df[df['ph'] <= 14]

    # Clean by Sally's criteria
    # if the cluster centroid is 12642, skip to the next iteration of for loop
    if cluster_centroids[i] == 12642:
        continue
    if cluster_centroids[i] == 16540:
        # if Tm=20, delete the row
        df = df.drop(df[df['tm'] == 20].index)

    # resave the csv file in the output directory
    df.to_csv(f'{output_dir}/cleaned_cluster_{cluster_centroids[i]}_features.csv')
