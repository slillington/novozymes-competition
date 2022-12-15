"""
Combining the cleaned datasets into one big dataset and scaling the features.
Sam Lobo
"""
#%%
import numpy as np
import pandas as pd
import os

def scale_and_combine():
    # Note the path to this file
    dir_path = os.path.dirname(os.path.abspath(__file__))
    # Output directory
    output_dir = os.path.join(dir_path, 'datasets_scaled_and_cleaned')
    # Note the directory of the cleaned data with energies
    cleaned_data_path = os.path.join(dir_path, 'datasets_cleaned')
    # List the files in the data folder
    cleaned_data_files = os.listdir(cleaned_data_path)
    # Sort the feature files list by the centroid number
    cleaned_data_files.sort(key=lambda x: int(x.split('_')[2]))
    # Get a list of just the centroid numbers in the feature files list
    cluster_centroids = [int(x.split('_')[2]) for x in cleaned_data_files]

    # # turns out we deleted 13268 because it had a fake melting point of 25
    # # so we need to delete it from the list of centroids and from the list of files
    # cluster_centroids.remove(13268)
    # cleaned_data_files.remove('cleaned_cluster_13268_features.csv')

    # Make a list of dataframes from each cluster
    dataframes = []
    for file in cleaned_data_files:
        if file.endswith('.csv'):
            dataframes.append(pd.read_csv(f'{cleaned_data_path}/{file}', index_col=0))

    # Standardize (or normalize) each dataframe
    # Further clean each dataframe if you'd like
    scaled_dataframes = []
    for i, df in enumerate(dataframes):
        # Get rid of the .pdb at the end of the index and turn it to int
        df.index = [int(f.split('.')[0]) for f in df.index]

        # # Filtering by pH is surely important! Play with this later.
        # # delete pH's that are smaller than 5
        # df = df[df['ph'] >= 5]

        # If you want to ignore the energy features, uncomment this:
        # Filtering out the energy features (beginning with e_)
        # df = df.filter(regex='^(?!e_)') # This is a regex made by Copilot, test it before using

        # Make every column (except the Tm & pH columns!) have a mean of 0
        df.iloc[:, 2:] -= df.iloc[:, 1:].mean()
        # and a standard deviation of 1 (this is called standardizing)
        df.iloc[:, 2:] /= df.iloc[:, 1:].std()
        # Later play with NORMALIZING instead of standardizing ^

        # Subtract the Tm column by the Tm with index = cluster_centroids[i]
        try: df['tm'] = df['tm'] - df.loc[cluster_centroids[i], 'tm']
        except KeyError:
            # Subtract by the mean of the Tm column
            df['tm'] = df['tm'] - df['tm'].mean()
            # print(f'Cluster centroid {cluster_centroids[i]} was not found in the dataset. Subtracting by the mean of the Tm column.')
            # print('This is probably because we deleted it - we probably thought its melting point was fake since it was 25C.\n')
        
        # normalize the melting point column to fit between -1 and 1
        df['tm'] /= df['tm'].abs().max()

        # Add the dataframe to the list of cleaned dataframes
        scaled_dataframes.append(df)

        # Add '.pdb' to the end of the index
        df.index = df.index.astype(str) + '.pdb'

        # Save the dataframe to a csv file in the output directory
        df.to_csv(f'{output_dir}/scaled_cluster_{cluster_centroids[i]}_features.csv')

    # Combine all the dataframes into one big dataframe (use the cluster centroid as a key)
    combined_df = pd.concat(scaled_dataframes, keys=[f'cluster{cluster_centroid}' for cluster_centroid in cluster_centroids], names=['cluster', 'pdb'])
    
    # replace nan values with 0 (exclude index)
    combined_df.iloc[:, 1:] = combined_df.iloc[:, 1:].fillna(0)

    # Save the combined dataframe to a csv file and ignore the index
    combined_df.to_csv('dataset_combined.csv', index=True)
    
    # Split the combined dataframe into X (features) and y (melting point)
    X = combined_df.iloc[:, 1:].values
    y = combined_df['tm'].values

    return combined_df, X, y # X is the features, y is ΔTm

# make sure the input argument cluster is an int
def scale_and_combine_without_one_cluster(drop_cluster):
    # Note the path to this file
    dir_path = os.path.dirname(os.path.abspath(__file__))
    # Output directory
    output_dir = os.path.join(dir_path, 'datasets_scaled_and_cleaned')
    # Note the directory of the cleaned data with energies
    cleaned_data_path = os.path.join(dir_path, 'datasets_cleaned')
    # List the files in the data folder
    cleaned_data_files = os.listdir(cleaned_data_path)
    # Sort the feature files list by the centroid number
    cleaned_data_files.sort(key=lambda x: int(x.split('_')[2]))
    # Get a list of just the centroid numbers in the feature files list
    cluster_centroids = [int(x.split('_')[2]) for x in cleaned_data_files]

    # Make a list of dataframes from each cluster
    dataframes = []
    for file in cleaned_data_files:
        if file.endswith('.csv'):
            dataframes.append(pd.read_csv(f'{cleaned_data_path}/{file}', index_col=0))

    # Standardize (or normalize) each dataframe
    # Further clean each dataframe if you'd like
    scaled_dataframes = []
    for i, df in enumerate(dataframes):
        # Get rid of the .pdb at the end of the index and turn it to int
        df.index = [int(f.split('.')[0]) for f in df.index]

        # # Filtering by pH is surely important! Play with this later.
        # # delete pH's that are smaller than 5
        # df = df[df['ph'] >= 5]

        # If you want to ignore the energy features, uncomment this:
        # Filtering out the energy features (beginning with e_)
        # df = df.filter(regex='^(?!e_)') # This is a regex made by Copilot, test it before using

        # Drop the e_vdw and e_vdw14 columns
        # df = df.drop(columns=['e_vdw', 'e_vdw14'])


        # Make every column (except the Tm & pH columns!) have a mean of 0
        df.iloc[:, 2:] -= df.iloc[:, 1:].mean()
        # and a standard deviation of 1 (this is called standardizing)
        df.iloc[:, 2:] /= df.iloc[:, 1:].std()
        # # This is an alternative block to normalize instead of standardizing
        # # normalize every column (except the Tm & pH columns!)
        # df.iloc[:, 2:] /= df.iloc[:, 1:].abs().max()

        # Subtract the Tm column by the Tm with index = cluster_centroids[i]
        try: df['tm'] = df['tm'] - df.loc[cluster_centroids[i], 'tm']
        except KeyError:
            # Subtract by the mean of the Tm column
            df['tm'] = df['tm'] - df['tm'].mean()
            # print(f'Cluster centroid {cluster_centroids[i]} was not found in the dataset. Subtracting by the mean of the Tm column.')
            # print('This is probably because we deleted it - we probably thought its melting point was fake since it was 25C.\n')
        
        # normalize the melting point column to fit between -1 and 1
        df['tm'] /= df['tm'].abs().max()

        # Add the dataframe to the list of cleaned dataframes
        scaled_dataframes.append(df)

        # Add '.pdb' to the end of the index
        df.index = df.index.astype(str) + '.pdb'

        # Save the dataframe to a csv file in the output directory
        df.to_csv(f'{output_dir}/scaled_cluster_{cluster_centroids[i]}_features.csv')

    # Combine all the dataframes into one big dataframe (use the cluster centroid as a key)
    combined_df = pd.concat(scaled_dataframes, keys=[f'cluster{cluster_centroid}' for cluster_centroid in cluster_centroids], names=['cluster', 'pdb'])
    
    ### THIS ONE LINE IS THE ONLY DIFFERENCE BETWEEN THIS FUNCTION AND THE ONE ABOVE
    # delete the rows with key cluster{drop_cluster}
    combined_df = combined_df.drop(index=f'cluster{drop_cluster}')

    # replace nan values with 0 (exclude index)
    combined_df.iloc[:, 1:] = combined_df.iloc[:, 1:].fillna(0)

    # Save the combined dataframe to a csv file and ignore the index
    combined_df.to_csv('dataset_combined.csv', index=True)
    
    # Split the combined dataframe into X (features) and y (melting point)
    X = combined_df.iloc[:, 1:].values
    y = combined_df['tm'].values

    return combined_df, X, y # X is the features, y is ΔTm

#%%
if __name__ == '__main__':
    combined_df, X, y = scale_and_combine()
    print(combined_df)

