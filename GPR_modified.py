# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 00:44:58 2022

@author: lshenaf@connect.ust.hk
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor 
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel


## T=40℃
train_X = np.array([0, 3.7, 7.4, 11.1]).reshape(-1, 1)
#train_y = np.array([8.39E-3, 1.37E-2, 1.64E-2, 1.27E-2]).reshape(-1, 1) #40C
#train_y = np.array([1.93E-2, 2.70E-2, 3.31E-2, 2.46E-2]).reshape(-1, 1) #50C
#train_y = np.array([3.73E-2, 4.85E-2, 6.33E-2, 4.59E-2]).reshape(-1, 1) #60C
#train_y = np.array([6.90E-2, 8.70E-2, 1.17E-1, 7.92E-2]).reshape(-1, 1) #70C
#train_y = np.array([1.10E-1, 1.43E-1, 1.84E-1, 1.44E-1]).reshape(-1, 1) #80C
#train_y = np.array([1.79E-1, 2.24E-1, 3.30E-1, 2.28E-1]).reshape(-1, 1) #90C
#train_y = np.array([2.86E-1, 3.82E-1, 4.71E-1, 3.46E-1]).reshape(-1, 1) #100C
train_y = np.array([5.11E-1, 5.90E-1, 8.09E-1, 5.56E-1]).reshape(-1, 1) #110C
test_X = np.arange(0, 1, 0.001).reshape(-1, 1)


#Normalizing 
train_X_norm = np.array([0,0.33333333,0.6666667,1]).reshape(-1,1)
train_y_norm = np.array([i / np.max(train_y) for i in train_y]).reshape(-1, 1)


# fit GPR
kernel = ConstantKernel(constant_value=0.2, constant_value_bounds=(1e-4, 1e4)) * RBF(length_scale=0.5, length_scale_bounds=(1e-4, 1e4)) + WhiteKernel(noise_level=2e-5, noise_level_bounds=(2e-5,2e-5))
gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=2)
gpr.fit(train_X_norm, train_y_norm)
print(gpr.get_params())
mu, cov = gpr.predict(test_X, return_cov=True)
test_y = mu.ravel()
uncertainty = 1.96 * np.sqrt(np.diag(cov))

# plotting
plt.figure()

plt.fill_between(test_X.ravel() * max(train_X), (test_y + uncertainty) * max(train_y), (test_y - uncertainty) * max(train_y), alpha=0.4)
plt.plot(test_X * max(train_X), test_y * max(train_y), label="predict")
plt.scatter(train_X, train_y, label="train", c="red", marker="x")
plt.legend()

test_y_max = max(test_y + uncertainty) * max(train_y)
test_y_max_index = np.argmax(test_y + uncertainty)
print('The predited maximum ionic conductivity is:', test_y_max, 'mS/cm')
print('The index of the predited maximum ionic conductivity is:', test_y_max_index * max(train_X) / 1000)


f = open('110.txt', 'w')
for i in range(len(test_X)):
    f.write(str(test_X[i] * max(train_X)) + '      ')
    f.write(str(test_y[i] * max(train_y)) + '      ')
    f.write(str((test_y[i] + uncertainty[i]) * max(train_y)) + '      ')
    f.write(str((test_y[i] - uncertainty[i]) * max(train_y)) + '      ' + '\n')

f.write('The predited maximum ionic conductivity is:' + str(test_y_max) + 'mS/cm' + '\n')
f.write('The index of the predited maximum ionic conductivity is:' + str(test_y_max_index))
f.close()