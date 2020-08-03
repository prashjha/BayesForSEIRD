import sys

from SALib.analyze import morris
from SALib.analyze import delta
from SALib.sample.morris import sample
import matplotlib.pyplot as plt

import numpy as np

# import forward model
# from model import model

# define parameters and range
PARAM_DIM = 9
problem = {
 'num_vars': PARAM_DIM,
 'names': ['A', 'beta_e', 'nu_s', 'nu_i', 'gamma_e', 'gamma_r', 'gamma_d', 'sigma', 'tr'],
 'groups': None,
 'bounds': [[10., 1200.],
            [5.e-5, 1.e-3],
            [1.e-7, 1.e-4],
            [5.e-9, 5.e-6],
            [0.8 * 1./6., 1.25 * 1./6.],
            [0.75 * 1./24., 1.33 * 1./24.],
            [0.5 * 1./160., 2. * 1./160.],
            [0.8 * 1./7., 1.25 * 1./7.],
            [0.2, 20.]]
}

# Generate samples. Number of samples are N * (num_paramsa+1)
param_values = sample(problem, N=200, num_levels=4,
                      optimal_trajectories=None)

print('N samples: {}'.format(len(param_values)))
np.savetxt('samples_small_param_dim_' + str(PARAM_DIM) + '.txt', param_values)
