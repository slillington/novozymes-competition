"""
Predict delta melting temperatures
Sally Jiao

"""
import numpy as np
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

## CHANGE OPTIONS HERE TO MODIFY DATASETS USED, TYPE OF FIT ##################

exclude_set = '972'
#fname_dataset = 'dataset_combined_exclude{}_ph7.5.csv'.format(exclude_set)
fname_dataset = 'dataset_combined_exclude{}_ph7.5_tol1.0.csv'.format(exclude_set)
i_wt = -1 # index of line in cleaned dataset with WT
model_type = 'RF' # options: linear, RF, ANN, LASSO
fname_dataset_validation = 'datasets_ph/{}_ph7.5_tol1.0.csv'.format(exclude_set) # use None if don't want this

##############################################################################

# load data
dataset = np.loadtxt(fname_dataset, delimiter=',', skiprows=1)
features = dataset[:, 1:]
metrics = dataset[:, 0]

# compute delta features and metrics
if not 'dataset_combined' in fname_dataset:
    features = features - features[i_wt, :]
    metrics = metrics - metrics[i_wt]

# load validation data and compute delta features and metrics
if fname_dataset_validation is not None:
    dataset_validation = np.loadtxt(fname_dataset_validation, delimiter=',', skiprows=1)
    features_validation = dataset_validation[:, 1:] 
    metrics_validation = dataset_validation[:, 0] 

print('shape of features', np.shape(features))
print('length of metrics', len(metrics))

print('std dev of metrics', np.std(metrics))

kf = KFold(n_splits=5, shuffle=True)
train_rmse = np.zeros(5)
test_rmse = np.zeros(5)
corr_spearman_train = np.zeros(5)
corr_spearman_test = np.zeros(5)
if fname_dataset_validation is not None:
    validation_rmse = np.zeros(5)
    validation_spearman = np.zeros(5)
    
i_split = 0
# for each fold
for train_index, test_index in kf.split(features):

    print('\nSplit number', i_split)

    # get training and testing set
    features_train, features_test = features[train_index], features[test_index]
    metrics_train, metrics_test = metrics[train_index], metrics[test_index]

    # scale features and metrics
    scaler_features = StandardScaler()
    features_train_scaled = scaler_features.fit_transform(features_train)
    features_test_scaled = scaler_features.transform(features_test)
    
    scaler_metrics = StandardScaler()
    metrics_train_scaled = scaler_metrics.fit_transform(metrics_train.reshape(-1, 1))

    if model_type == 'RF': # default RF settings
        model = RandomForestRegressor()
        model.fit(features_train_scaled, metrics_train_scaled.flatten())
        print('coefficients', model.feature_importances_)
    elif model_type == 'linear':
        model = LinearRegression()
        model.fit(features_train_scaled, metrics_train_scaled)        
        print('coefficients')
        print(model.coef_)
        print('intercept')
        print(model.intercept_)
    elif model_type == 'ANN': # 2 hidden layers, 35 nodes each
        model = MLPRegressor(hidden_layer_sizes=(35,35,), max_iter=10000)
        model.fit(features_train_scaled, metrics_train_scaled.flatten())
    elif model_type == 'LASSO': # alpha parameter 0.01
        model = Lasso(alpha=0.01)
        model.fit(features_train_scaled, metrics_train_scaled)
        print('coefficients')
        print(model.coef_)
        print('intercept')
        print(model.intercept_)
    
    # compute predicted metrics
    metrics_train_predict = scaler_metrics.inverse_transform(model.predict(features_train_scaled).reshape(-1, 1))
    metrics_test_predict = scaler_metrics.inverse_transform(model.predict(features_test_scaled).reshape(-1, 1))

    # compute fit quality metrics
    train_rmse[i_split] = np.sqrt(mean_squared_error(metrics_train, metrics_train_predict))
    test_rmse[i_split] = np.sqrt(mean_squared_error(metrics_test, metrics_test_predict))

    corr_spearman_train[i_split] = spearmanr(metrics_train, metrics_train_predict)[0]
    corr_spearman_test[i_split] = spearmanr(metrics_test, metrics_test_predict)[0]

    if fname_dataset_validation is not None:
        metrics_validation_predict = scaler_metrics.inverse_transform(model.predict(scaler_features.transform(features_validation)).reshape(-1, 1))
        validation_rmse[i_split] = np.sqrt(mean_squared_error(metrics_validation, metrics_validation_predict))
        validation_spearman[i_split] = spearmanr(metrics_validation, metrics_validation_predict)[0]
    
    i_split += 1

    # plot fits
    plt.subplot(2, 3, i_split)
    plt.plot(metrics_train, metrics_train_predict, 'o', label='Training')
    plt.plot(metrics_test, metrics_test_predict, 'o', label='Testing')
    if fname_dataset_validation is not None:
        plt.plot(metrics_validation, metrics_validation_predict, 'o', label='Validation')
    if i_split == 4:
        plt.xlabel('Actual')
        plt.ylabel('Predicted')
    if i_split == 1:
        plt.legend(edgecolor='None', fontsize=6)
    
# print fit quality metrics
print('training rmse', train_rmse)
print('testing rmse', test_rmse)
print('training spearman', corr_spearman_train)
print('testing spearman', corr_spearman_test)
if fname_dataset_validation is not None:
    print('validation rmse', validation_rmse)
    print('validation spearman', validation_spearman)
plt.savefig('fit.png',dpi=500,transparent=False,bbox_inches="tight")
plt.show()
