"""
Predict melting temperatures
Sally Jiao
"""
import numpy as np
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error

dataset = np.loadtxt('dataset.csv', delimiter=',')

features = dataset[:, 1:]
metrics = dataset[:, 0]

kf = KFold(n_splits=5, shuffle=True)
train_rmse = np.zeros(5)
test_rmse = np.zeros(5)

i_split = 0
for train_index, test_index in kf.split(patterns_all):

    print('\nSplit number', i_split)
    
    features_train, features_test = features[train_index], features[test_index]
    metrics_train, metrics_test = metrics[train_index], metrics[test_index]

    # scale features and metrics
    scaler_features = StandardScaler()
    features_train_scaled = scaler_features.fit_transform(features_train)
    features_test_scaled = scaler_features.transform(features_test)
    
    scaler_metrics = StandardScaler()
    metrics_scaled = scaler_metrics.fit_transform(metrics_train)

    model = LinearRegression()
    model.fit(features_train_scaled, metrics_train_scaled)

    print('coefficients')
    print(model.coef_)
    print('intercept')
    print(model.intercept_)
    
    metrics_train_predict = scaler_metrics.inverse_transform(model.predict(features_train_scaled))
    metrics_test_predict = scaler_metrics.inverse_transform(model.predict(features_test_scaled))

    train_rmse[i_split] = np.sqrt(mean_squared_error(metrics_train, metrics_train_predict))
    test_rmse[i_split] = np.sqrt(mean_squared_error(metrics_test, metrics_test_predict))
    
    i_split += 1

print('mean train_rmse', np.mean(train_rmse))
print('mean test_rmse', np.mean(test_rmse))
print('training rmse', train_rmse)
print('testing rmse', test_rmse)

