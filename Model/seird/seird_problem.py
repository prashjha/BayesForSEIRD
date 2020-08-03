from __future__ import absolute_import, division, print_function

import dolfin as dl
import numpy as np
from picard_solver import picard_solver
from seird_forms import seird_forms

STATE = 0
PARAMETER = 1

class seird_fwd_problem:
    def __init__(self, Vh, simulation_time, dt, init_vectors, subdmn_path = './subdomains/', mesh_tag = 'mesh_5h', qoi_type='state', reset_exposed = True, save = False, save_days = True, out_path = './fwd_result/'):
        
        self.Vh = Vh
        self.init_vectors = init_vectors
        self.dt = dt
        self._save = save
        self._reset_exposed = reset_exposed
        self._qoi_type = qoi_type

        self._save_days = save_days
        
        if not simulation_time.is_integer():
            raise ValueError("The total simulation time is should be a whole number.")
        else:
            self.T = round(simulation_time)
        
        self.nt = round(self.T/self.dt)
        if not self.nt*self.dt == self.T:
            raise ValueError("t = 1, 2, ... cannot be reached with the given time step.")
            
        self.out_freq = round(1./dt)
    
        _linear_solver = dl.PETScKrylovSolver("gmres", "ilu")
        _linear_solver.parameters["nonzero_initial_guess"] = True
        
        self._solver = picard_solver(_linear_solver)
        self.problem = seird_forms(Vh, dt, save, out_path, subdmn_path, mesh_tag, qoi_type)
        self.u_0 = self.generate_pde_state()
        self.u = self.generate_pde_state()
        self.u_ic = self.generate_pde_state()
        
        self._set_initial_conditions(self.u_ic)
            
    def generate_state(self):
        """ Return a vector in the shape of the fwd model output. """
        return np.empty((self.T,2))
    
    def generate_state_district(self):
        """ Return a vector in the shape of the fwd model output. """
        return np.empty((self.T, 26, 2))

    def generate_pde_state(self):
        """ Return a list of vectors that correspons to the SEIRD state vectors. """
        return [dl.Function(self.Vh[STATE]).vector() for i in range(5)]

    def generate_parameter(self):
        """ Return a vector in the shape of the parameter. """
        return dl.Function(self.Vh[PARAMETER]).vector()

    def init_parameter(self, m):
        """ Initialize the parameter. """
        dummy = self.generate_parameter()
        m.init( dummy.mpi_comm(), dummy.local_range() )
        
    def _assign_vectors(self, u1, u2):
        for i in range(5):
            u1[i].zero()
            u1[i].axpy(1, u2[i])

    def _set_initial_conditions(self,u):
        for i in range(5):
            if not (self._reset_exposed and i == 1):
                u[i].zero()
                u[i].set_local(self.init_vectors[i])
        
    def solveFwd(self, out, x):
        
#        print("Solve with parameters ", np.exp(x[PARAMETER].get_local()))
        
        self._assign_vectors(self.u, self.u_ic)
        if self._reset_exposed:
            self.u[1].zero()
            self.u[1].axpy(np.exp(x[PARAMETER].get_local())[-1], self.u[2])
        
        self.problem.set_parameters(x[PARAMETER])
        store_index = 0
        time = 0.0
        for time_index in range (self.nt):
            self._assign_vectors(self.u_0, self.u)
            self._solver.solve(self.problem, self.u, self.u_0)
            if time_index % self.out_freq == 0:
                if self._qoi_type == 'state'
                    out[store_index, :] = self.problem.evaluate(self.u)*10000.
                else:
                    out[store_index] = self.problem.evaluate_district(self.u)*10000.
                store_index +=1
            time += self.dt
            if self._save:
                if self._save_days and time_index % self.out_freq == 0:
                    self.problem.save(self.u, time)

                if self._save_days == False:
                    self.problem.save(self.u, time)
        if not store_index == self.T:
            raise Exception("The forwad solve output does not match the data.")

