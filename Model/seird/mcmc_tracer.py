import numpy as np
import dolfin as dl

STATE = 0
PARAMETER = 1

class FullTracer:
    def __init__(self, Vh, simulation_time, dir_name, save = True, print = True):
        self.posterior = None
        self.qoi = None
        self.acceptance_rate = None
        self.cost = None
        self._dir_name = dir_name
        self._time = round(simulation_time)
        self._parameter = None
        self._param_dim = Vh[PARAMETER].dim()
        self.accepted = 0
        self.rejected = 0
        self._save = save
        self._print = print
        
    def append(self, current, q):
        
        if self._parameter is None:
            self._parameter = current.m.get_local()
            acceptance_rate = 0.0
        elif not np.array_equal(self._parameter, current.m.get_local()):
            self.accepted += 1
            acceptance_rate = self.accepted/float(self.accepted + self.rejected+1)
            self.write(current, acceptance_rate)
        else:
            self.rejected += 1
            acceptance_rate = self.accepted/float(self.accepted + self.rejected+1)
        if self._print:
            print("Acceptance ratio: ", acceptance_rate*100, "%")
    
    def write(self, current, acceptance_rate):
        self._parameter = current.m.get_local()
        m_vector = current.m.get_local().reshape((1,self._param_dim))
        q_vector = current.u.reshape((1, self._time, 2))
        ar_vector = np.array([acceptance_rate])
        cost_vector = np.array([current.cost])
        if self.posterior is None:
            self.posterior = m_vector
        else:
            self.posterior = np.append(self.posterior, m_vector, axis = 0)
        if self.qoi is None:
            self.qoi = q_vector
        else:
            self.qoi = np.append(self.qoi, q_vector, axis = 0)
        if self.acceptance_rate is None:
            self.acceptance_rate = ar_vector
        else:
            self.acceptance_rate = np.append(self.acceptance_rate, ar_vector)
        if self.cost is None:
            self.cost = cost_vector
        else:
            self.cost = np.append(self.cost, cost_vector)
        if self._save and self.accepted % 10 == 0:
            self.save()
            
    def save(self):
        np.save(self._dir_name +'/param_samples.npy', self.posterior)
        np.save(self._dir_name +'/qoi.npy', self.qoi)
        np.save(self._dir_name +'/acceptance_rate.npy', self.acceptance_rate)
        np.save(self._dir_name +'/cost.npy', self.cost)
        
