import dolfin as dl
import numpy as np
import math
from datetime import datetime

STATE = 0
PARAMETER = 1

def validate_date(date_text):
    try:
        datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        raise ValueError("Incorrect data format, should be YYYY-MM-DD")

class seird_misfit:
    def __init__(self, infected_cases, deceased_cases, date, simulation_time, truncated_gaussian = False, data_start_date="2020-03-06"):
        
        self.noise_variance = None
        self.truncated_gaussian = truncated_gaussian
        
        validate_date(date)
        d0 = datetime.strptime(data_start_date, "%Y-%m-%d")
        d1 = datetime.strptime(date, "%Y-%m-%d")
        day_index = abs((d1-d0)).days
        
        if not infected_cases.size == deceased_cases.size:
            raise IndexError("The total decease cases data and the total infected cases data does not match in size. Please re-check the data")
        
        self.data = np.empty((round(simulation_time), 2))
        self.data[:, 0] = infected_cases[day_index+1:(day_index + round(simulation_time)+1):1]
        self.data[:, 1] = deceased_cases[day_index+1:(day_index + round(simulation_time)+1):1]

    def set_noise_variance(self, percentages):
        """Setting the noise variance based on the percentages (numpy array of 2 values) at one standard deviation"""
        self.noise_variance = np.empty_like(self.data)
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                self.noise_variance[i,j] = (self.data[i,j]*percentages[j])**2
        
    def cost(self, x):
        
        if self.noise_variance is None:
            raise ValueError("Noise Variance must be specified")
        
        if not x[STATE].shape == self.data.shape:
            raise IndexError("The state output is of shape ", x[STATE].shape, ", while data is of shape ", self.data.shape, ".")
            
        if self.truncated_gaussian:
            misfit = 0.0
            for i in range(self.data.shape[0]):
                for j in range(self.data.shape[1]):
                    if x[STATE][i, j] < self.data[i,j]:
                        misfit = math.inf
                        break
                    else:
                        misfit += 0.5*(x[STATE][i, j] - self.data[i,j])**2/(self.noise_variance[i,j]) + math.log(2.)
        
            return misfit
        
        else:
            return np.sum(np.divide(np.power(self.data - x[STATE], 2), 2*self.noise_variance))
            
