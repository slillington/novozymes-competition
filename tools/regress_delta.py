"""
Predict melting temperatures
Sally Jiao
"""
import numpy as np
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

dataset = np.loadtxt('972_data.csv', delimiter=',', skiprows=1)
#dataset = np.loadtxt('972_data_cp.csv', delimiter=',', skiprows=1)

features = dataset[:, 1:]
metrics = dataset[:, 0]

i_wt = 0
features = features - features[i_wt, :]
metrics = metrics - metrics[i_wt]

print('shape of features', np.shape(features))
print('length of metrics', len(metrics))

print('std dev of metrics', np.std(metrics))

kf = KFold(n_splits=5, shuffle=True)
train_rmse = np.zeros(5)
test_rmse = np.zeros(5)

i_split = 0
for train_index, test_index in kf.split(features):

    print('\nSplit number', i_split)
    
    features_train, features_test = features[train_index], features[test_index]
    metrics_train, metrics_test = metrics[train_index], metrics[test_index]

    # scale features and metrics
    scaler_features = StandardScaler()
    features_train_scaled = scaler_features.fit_transform(features_train)
    features_test_scaled = scaler_features.transform(features_test)
    
    scaler_metrics = StandardScaler()
    metrics_train_scaled = scaler_metrics.fit_transform(metrics_train.reshape(-1, 1))

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

    plt.subplot(2, 3, i_split)
    plt.plot(metrics_train, metrics_train_predict, 'o')
    plt.plot(metrics_test, metrics_test_predict, 'o')
    if i_split == 4:
        plt.xlabel('Actual')
        plt.ylabel('Predicted')

print('mean train_rmse', np.mean(train_rmse))
print('mean test_rmse', np.mean(test_rmse))
print('training rmse', train_rmse)
print('testing rmse', test_rmse)

plt.show()
