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
df_combined, X_train, y_train, X_test = scale_and_combine()
X_train, y_train = shuffle(X_train, y_train, random_state=0)

# df_combined, X, y = scale_and_combine()
model_type = 'LASSO' # options: linear, RF, ANN, LASSO
print(f'MODEL TYPE: {model_type}')

# print notes
Notes = 'None'
# Notes = 'ignoring features starting with n_'
print(f'Notes: {Notes}')

# Load the validation dataset and store its features and melting points
fname_dataset_validation = None
# fname_dataset_validation = f'datasets_scaled_and_cleaned/scaled_cluster_{validation_cluster}_features.csv' # use None if don't want this
# df_validation_cluster = pd.read_csv(fname_dataset_validation, index_col=0)
# df_validation_cluster = df_validation_cluster.fillna(0) # replacing nan with 0
# X_validation = df_validation_cluster.iloc[:, 1:].values
# y_validation = df_validation_cluster['tm'].values
# X_validation, y_validation = shuffle(X_validation, y_validation, random_state=0)


#%%
# for each fold

    #print('\nSplit number', i_split)

    # # get training and testing set
    # X_train, X_test = X[train_index], X[test_index]
    # y_train, y_test = y[train_index], y[test_index]

if model_type == 'RF': # default RF settings
    model = RandomForestRegressor()
    model.fit(X_train, y_train.flatten())
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
    # print('coefficients')
    # print(model.coef_)
    # print('intercept')
    # print(model.intercept_)

# compute predicted metrics
y_train_predict = model.predict(X_train)

# compute fit quality metrics
train_rmse = np.sqrt(mean_squared_error(y_train, y_train_predict))
corr_spearman_train = spearmanr(y_train, y_train_predict)[0]

# compute Tm for test set
y_validation_predict = model.predict(X_test)
# seqence ids
seq_ids = np.arange(31390,33802+1)
# create dataframe
df = pd.DataFrame({'seq_id':seq_ids, 'tm':y_validation_predict})
df['tm'] = np.around(df['tm'],6)
df.to_csv('submission5.csv',index=False)

# np.savetxt('submission.csv', np.hstack((seq_id.reshape(-1, 1), metrics_validation_predict.reshape(-1, 1))), delimiter=',', header='seq_id,tm', fmt=['%d', '%6.2f'],
#            comments='')

#%%
# plot results for training set
plt.plot(y_train, y_train_predict, 'o', label='Training')
plt.xlabel('Actual')
plt.ylabel('Predicted')

# print fit quality metrics
print('training rmse', train_rmse)
print('training spearman', corr_spearman_train)
plt.savefig('fit.png',dpi=500,transparent=False,bbox_inches="tight")
plt.show() 
# %%
