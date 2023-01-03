"""
Predict delta melting temperatures
Sally Jiao and Sam Lobo
"""
#%%
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from scale_and_combine import scale_and_combine
from scale_and_combine import scale_and_combine_without_one_cluster
import os
from sklearn.utils import shuffle

# Make sure the current working directory is the directory of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

## CHANGE OPTIONS HERE TO MODIFY DATASETS USED, TYPE OF FIT ##################
# load the dataframe, features (X), and metric (y = melting points) for the combined dataset
# adjust the scale_and_combine function in scale_and_combine.py to see if you can get better results
validation_cluster = 24305
df_combined, X, y = scale_and_combine_without_one_cluster(validation_cluster)
X, y = shuffle(X, y, random_state=0)


# df_combined, X, y = scale_and_combine()
model_type = 'LASSO' # options: linear, RF, ANN, LASSO
print(f'MODEL TYPE: {model_type}')

# print notes
Notes = 'Dropping no columns'
print(f'Notes: {Notes}')

# Load the validation dataset and store its features and melting points
fname_dataset_validation = f'datasets_scaled_and_cleaned/scaled_cluster_{validation_cluster}_features.csv' # use None if don't want this
df_validation_cluster = pd.read_csv(fname_dataset_validation, index_col=0)
df_validation_cluster = df_validation_cluster.fillna(0) # replacing nan with 0
X_validation = df_validation_cluster.iloc[:, 1:].values
y_validation = df_validation_cluster['tm'].values
X_validation, y_validation = shuffle(X_validation, y_validation, random_state=0)

#%%

K = 10 # number of folds
kf = KFold(n_splits=K, shuffle=True)
train_rmse = np.zeros(K)
test_rmse = np.zeros(K)
corr_spearman_train = np.zeros(K)
corr_spearman_test = np.zeros(K)
if fname_dataset_validation is not None:
    validation_rmse = np.zeros(K)
    validation_spearman = np.zeros(K)

#%%
i_split = 0
# for each fold
for train_index, test_index in kf.split(X):

    #print('\nSplit number', i_split)

    # get training and testing set
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]

    if model_type == 'RF': # default RF settings
        model = RandomForestRegressor()
        model.fit(X_train, y_train.flatten())
        print('coefficients', model.feature_importances_)
        # Find the coefficients of the features in the model
        coefficients = model.feature_importances_
        # Find the names of the features that have been removed by LASSO
        removed_features = [feature for feature, coef in zip(df_combined.columns[1:], coefficients) if coef < 0.01]
        print("Highly reduced features:", removed_features)

        # Find the names and coefficients of the remaining features
        remaining_features = [(feature, coef) for feature, coef in zip(df_combined.columns[1:], coefficients) if coef >= 0.01]
    elif model_type == 'linear':
        model = LinearRegression()
        model.fit(X_train, y_train)        
        # print('coefficients')
        # print(model.coef_)
        # print('intercept')
        # print(model.intercept_)
    elif model_type == 'ANN': # 2 hidden layers, 35 nodes each
        model = MLPRegressor(hidden_layer_sizes=(35,35,), max_iter=10000)
        model.fit(X_train, y_train.flatten())
    elif model_type == 'LASSO': # alpha parameter 0.01
        model = Lasso(alpha=0.005)
        model.fit(X_train, y_train)
        print('coefficients')
        # Find the coefficients of the features in the model
        coefficients = model.coef_
        # Find the names of the features that have been removed by LASSO
        removed_features = [feature for feature, coef in zip(df_combined.columns[1:], coefficients) if coef == 0]
        print("Removed features:", removed_features)

        # Find the names and coefficients of the remaining features
        remaining_features = [(feature, coef) for feature, coef in zip(df_combined.columns[1:], coefficients) if coef != 0]
        print("Remaining features:", remaining_features)
        # print(model.coef_)
        # print('intercept')
        # print(model.intercept_)
        
    # compute predicted metrics
    y_train_predict = model.predict(X_train)
    y_test_predict = model.predict(X_test)

    # compute fit quality metrics
    train_rmse[i_split] = np.sqrt(mean_squared_error(y_train, y_train_predict))
    test_rmse[i_split] = np.sqrt(mean_squared_error(y_test, y_test_predict))

    corr_spearman_train[i_split] = spearmanr(y_train, y_train_predict)[0]
    corr_spearman_test[i_split] = spearmanr(y_test, y_test_predict)[0]

    if fname_dataset_validation is not None:
        y_validation_predict = model.predict(X_validation)
        validation_rmse[i_split] = np.sqrt(mean_squared_error(y_validation, y_validation_predict))
        validation_spearman[i_split] = spearmanr(y_validation, y_validation_predict)[0]
    
    i_split += 1

    # plot fits
    if K == 5:
        plt.subplot(2, 3, i_split)
        plt.plot(y_train, y_train_predict, 'o', label='Training')
        plt.plot(y_test, y_test_predict, 'o', label='Testing')
        if fname_dataset_validation is not None:
            plt.plot(y_validation, y_validation_predict, 'o', label='Validation')

        # put label 'Actual' on bottom row of plots
        if i_split == 4 or i_split == 5 or i_split == 6: plt.xlabel('Actual')
        # put label 'Predicted' on left column of plots
        if i_split == 1 or i_split == 4: plt.ylabel('Predicted')
        if i_split == 1: plt.legend(edgecolor='None', fontsize=6)


    if K == 10:
        plt.subplot(3, 4, i_split)
        plt.plot(y_train, y_train_predict, 'o', label='Training')
        plt.plot(y_test, y_test_predict, 'o', label='Testing')
        if fname_dataset_validation is not None:
            plt.plot(y_validation, y_validation_predict, 'o', label='Validation')
        
        # put label 'Actual' on bottom row of plots
        if i_split == 7 or i_split == 8 or i_split == 9 or i_split == 10: plt.xlabel('Actual')
        # put label 'Predicted' on left column of plots
        if i_split == 1 or i_split == 5 or i_split == 9: plt.ylabel('Predicted')
        if i_split == 1: plt.legend(edgecolor='None', fontsize=6)

#%%
# print fit quality metrics
# print(f'training rmse:      {np.around(train_rmse,3)}')
# print(f'testing rmse:       {np.around(test_rmse,3)}')
# print(f'\ntraining spearman:  {np.around(corr_spearman_train,3)}')
# print(f'testing spearman:   {np.around(corr_spearman_test,3)}')
if fname_dataset_validation is not None:
    # print(f'\nvalidation rmse:    {np.around(validation_rmse,3)}')
    # print(f'validation spearman:{np.around(validation_spearman,3)}\n')
    # print the mean +/- std dev of the validation spearman correlation
    print(f'validation spearman mean: {np.around(np.mean(validation_spearman),3)} +/- {np.around(np.std(validation_spearman),3)}')
plt.show()
