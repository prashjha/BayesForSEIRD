import numpy as np

samples = np.load('validation_posterior_samples.npy')
#print(samples.shape)

mean = np.mean(samples, axis=0)
std = np.cov(samples, rowvar=False)
print('mean: {}'.format(mean))
print('std: {}'.format(std))
