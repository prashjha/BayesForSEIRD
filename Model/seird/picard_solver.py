SU = 0
EX = 1
IN = 2
RE = 3
DE = 4

class picard_solver:
    def __init__(self, linear_solver, tol = 1.0e-6, max_iter = 100):
    
        self._max_iter = max_iter
        self._tol = tol
        self._solver = linear_solver
        self.u_k = None
        
    def _assign_vectors(self, u1, u2):
        for i in range(5):
            u1[i].zero()
            u1[i].axpy(1., u2[i])
        
    def solve(self, problem, u, u_0):
        if self.u_k is None:
            self.u_k = [u_0[i].copy() for i in range(5)]
        else:
            self._assign_vectors(self.u_k, u_0)
        converged = False
        iter = 0
        while not converged and iter < self._max_iter:
            H, b = problem.assemble_systems(self.u_k, u_0)
            for i in range(5):
                self._solver.solve(H[i], u[i], b[i])
            error = problem.error_norm(u, self.u_k)
            if error <= self._tol:
                converged = True
            else:
                self._assign_vectors(self.u_k, u)
            iter += 1
