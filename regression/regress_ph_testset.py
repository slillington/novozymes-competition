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

fname_dataset = 'dataset_combined_all_ph7.5.csv'
model_type = 'RF' # options: linear, RF, ANN, LASSO
fname_dataset_validation = '../tools/cluster_test_set_features.csv' # use None if don't want this
i_wt_validation = 0

##############################################################################

# load data
dataset = np.loadtxt(fname_dataset, delimiter=',', skiprows=1)
features = dataset[:, 1:]
metrics = dataset[:, 0]

# load validation data and compute delta features and metrics
dataset_validation = np.loadtxt(fname_dataset_validation, delimiter=',', skiprows=1)
tm_wt = 0
features_validation = dataset_validation[:, 3:] - dataset_validation[i_wt_validation, 3:] # do not use pH as feature
seq_id = dataset_validation[:, 0]

print('shape of features', np.shape(features))
print('length of metrics', len(metrics))

print('std dev of metrics', np.std(metrics))

# scale features and metrics
scaler_features = StandardScaler()
features_train_scaled = scaler_features.fit_transform(features)    
scaler_metrics = StandardScaler()
metrics_train_scaled = scaler_metrics.fit_transform(metrics.reshape(-1, 1))

# fit model
if model_type == 'RF': # default RF settings
    model = RandomForestRegressor()
    model.fit(features_train_scaled, metrics_train_scaled.flatten())
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

# compute fit quality metrics
train_rmse = np.sqrt(mean_squared_error(metrics, metrics_train_predict))
corr_spearman_train = spearmanr(metrics, metrics_train_predict)[0]

# compute Tm for test set
metrics_validation_predict = scaler_metrics.inverse_transform(model.predict(scaler_features.transform(features_validation)).reshape(-1, 1)) + tm_wt
np.savetxt('submission.csv', np.hstack((seq_id.reshape(-1, 1), metrics_validation_predict.reshape(-1, 1))), delimiter=',', header='seq_id,tm', fmt=['%d', '%6.2f'],
           comments='')

# plot results for training set
plt.plot(metrics, metrics_train_predict, 'o', label='Training')
plt.xlabel('Actual')
plt.ylabel('Predicted')

# print fit quality metrics
print('training rmse', train_rmse)
print('training spearman', corr_spearman_train)
plt.savefig('fit.png',dpi=500,transparent=False,bbox_inches="tight")
plt.show()
