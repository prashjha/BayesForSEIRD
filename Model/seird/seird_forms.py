import dolfin as dl
import numpy as np
import math

STATE = 0
PARAMETER = 1

SU = 0
EX = 1
IN = 2
RE = 3
DE = 4

dl.parameters["form_compiler"]["optimize"]     = True
dl.parameters["form_compiler"]["cpp_optimize"] = True
dl.parameters["form_compiler"]["representation"] = "uflacs"
dl.parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"

class DistrictSubDomain(dl.SubDomain):
    
    def set(self, district, dist_geom):
        self.district = district
        self.dist_geom = dist_geom
    
    def inside(self, x, on_boundary):
        
        dist = []
        return self.dist_geom[self.district].contains(Point(x[0],x[1]))
    

class seird_forms:
    def __init__(self, Vh, dt, save, out_path, subdmn_path, mesh_tag, qoi_type):
    
        self.Vh = Vh
        self._param_dim = Vh[PARAMETER].dim()
        self.dt =dt
        self._out_path = out_path
        self._file = None
        self._qoi_type = qoi_type
        if save:
            self._file= [None]*5
            for i in range(5):
                names = ["susceptible", "exposed", "infected", "recovered", "deceased"]
                self._file[i] = dl.File(out_path + names[i] + '.pvd', "compressed")
        
        self.m = dl.Function(Vh[PARAMETER])
        self.gamma = [None]*5
        self.sigma   = None
        self.beta = [None]*5
        self.nu = [None]*5
        self.A = None
        
        self.help = dl.Function(Vh[STATE])
        
        self.p = dl.TestFunction(self.Vh[STATE])
        self.u_trial = dl.TrialFunction(self.Vh[STATE])
        self.M = dl.assemble(self.u_trial*self.p*dl.dx)
        
        self.z = dl.assemble(dl.Constant(1.0)*self.p*dl.dx)
        
        self.u_0 = [dl.Function(self.Vh[STATE]) for i in range(5)]
        self.u_k = [dl.Function(self.Vh[STATE]) for i in range(5)]
        
        self.H = [None]*5
        self.b = [None]*5
        
        self.a = [None]*5
        self.L = [None]*5
        
        # load subdomains
        self.subdmn_path = subdmn_path
        if qoi_type == 'district':
            subdomain = dl.MeshFunction("size_t", self.Vh[STATE].mesh(), 2)
            dl.File(self.subdmn_path + 'subdomain_' + mesh_tag + '.xml.gz') >> subdomain
            dx_dist = dl.Measure('dx', domain=self.Vh[STATE].mesh(), subdomain_data=subdomain)
            
            # for integration over districts
            self.z_dist = [dl.assemble(dl.Constant(1.0)*self.p*dx_dist(i)) for i in range(25)]
            
            # read area factor (integration of unit function over district subdomain is 
            # smaller than the district area)
            self.dist_area_factor = np.load(self.subdmn_path + 'area_factor_' + mesh_tag + '.npy')
        else:
            self.z_dist = []
            self.dist_area_factor = []
        
        
    def set_parameters(self, param):
    
        self.m.vector().zero()
        self.m.vector().axpy(1., param)
            
        if self._param_dim == 7:
            (self.A, self.beta[EX], self.nu[SU], self.gamma[RE], self.gamma[DE], self.sigma, tr) = dl.split(self.m)
            self.beta[IN] = self.beta[EX]
            self.nu[EX] = self.nu[SU]
            self.nu[RE] = self.nu[SU]
            self.nu[IN] = dl.Constant(-8*math.log(10.))
            self.gamma[EX] = dl.Constant(math.log(1./6.))
        elif self._param_dim == 8:
            (self.A, self.beta[EX], self.beta[IN], self.nu[SU], self.nu[EX], self.nu[IN], self.nu[RE], tr) = dl.split(self.m)
            self.gamma = [None, dl.Constant(math.log(1./6.)), None, dl.Constant(math.log(1./24.)), dl.Constant(math.log(1./160.))]
            self.sigma   = dl.Constant(math.log(1./7.))
        elif self._param_dim == 5:
            (self.A, self.beta[EX], self.nu[SU], self.nu[IN], tr) = dl.split(self.m)
            self.beta[IN] = self.beta[EX]
            self.nu[EX] = self.nu[SU]
            self.nu[RE] = self.nu[SU]
            self.gamma = [None, dl.Constant(math.log(1./6.)), None, dl.Constant(math.log(1./24.)), dl.Constant(math.log(1./160.))]
            self.sigma   = dl.Constant(math.log(1./7.))
        elif self._param_dim == 4:
            (self.A, self.beta[EX], self.nu[SU], self.nu[IN]) = dl.split(self.m)
            self.beta[IN] = self.beta[EX]
            self.nu[EX] = self.nu[SU]
            self.nu[RE] = self.nu[SU]
            self.gamma = [None, dl.Constant(math.log(1./6.)), None, dl.Constant(math.log(1./24.)), dl.Constant(math.log(1./160.))]
            self.sigma   = dl.Constant(math.log(1./7.))
        elif self._param_dim == 9:
            (self.A, self.beta[EX], self.nu[SU], self.nu[IN], self.gamma[EX], self.gamma[RE], self.gamma[DE], self.sigma, tr) = dl.split(self.m)
            self.beta[IN] = self.beta[EX]
            self.nu[EX] = self.nu[SU]
            self.nu[RE] = self.nu[SU]
        else:
            raise IndexError("The parameter dimension should be 4, 5, 7, 8 or 9.")
        
        self._set_varf_forms()
        
    def vector2Function(self, u_func, u_vec):
    
        u_func.vector().zero()
        u_func.vector().axpy(1., u_vec)
        
    def error_norm(self, u_k, u_0):
        
        error = 0.0
        for i in range(5):
            diff = u_k[i].copy()
            diff.axpy(-1., u_0[i])
            error +=math.sqrt(diff.inner(self.M*diff))

        return error
    
    def _set_varf_forms(self):
    
        F = [None]*5
    
        un_k = sum(self.u_k)
            
        # Define weak form for phi_s evolution
        f_s = - (1.- dl.exp(self.A)/un_k)*dl.exp(self.beta[IN])*self.u_trial*self.u_k[IN] - (1.- dl.exp(self.A)/un_k)*dl.exp(self.beta[EX])*self.u_trial*self.u_k[EX]
        F[SU]     = (self.u_trial-self.u_0[SU])*self.p*dl.dx - self.dt*f_s*self.p*dl.dx +\
                            self.dt*(un_k*dl.exp(self.nu[SU])*dl.inner(dl.grad(self.u_trial), dl.grad(self.p)))*dl.dx
                            
        # Define the weak form for phi_e evolution
        f_e = - dl.exp(self.sigma)*self.u_trial - dl.exp(self.gamma[EX])*self.u_trial + (1.- dl.exp(self.A)/un_k)*dl.exp(self.beta[IN])*self.u_k[SU]*self.u_k[IN] + (1.- dl.exp(self.A)/un_k)*dl.exp(self.beta[EX])*self.u_k[SU]*self.u_trial
        F[EX]     = (self.u_trial-self.u_0[EX])*self.p*dl.dx - self.dt*f_e*self.p*dl.dx +\
                            self.dt*(un_k*dl.exp(self.nu[EX])*dl.inner(dl.grad(self.u_trial), dl.grad(self.p)))*dl.dx

        # Define the weak form for phi_i evolution
        f_i = - dl.exp(self.gamma[DE])*self.u_trial - dl.exp(self.gamma[RE])*self.u_trial + dl.exp(self.sigma)*self.u_trial
        F[IN]     = (self.u_trial-self.u_0[IN])*self.p*dl.dx - self.dt*f_i*self.p*dl.dx +\
                            self.dt*(un_k*dl.exp(self.nu[IN])*dl.inner(dl.grad(self.u_trial), dl.grad(self.p)))*dl.dx
                            
        # Define the weak form for phi_r evolution
        f_r = dl.exp(self.gamma[RE])*self.u_k[IN] + dl.exp(self.gamma[EX])*self.u_k[EX]
        F[RE]     = (self.u_trial-self.u_0[RE])*self.p*dl.dx - self.dt*f_r*self.p*dl.dx +\
                            self.dt*(un_k*dl.exp(self.nu[RE])*dl.inner(dl.grad(self.u_trial), dl.grad(self.p)))*dl.dx
                            
        # Define the weak form for phi_d evolution
        F[DE]     = (self.u_trial-self.u_0[DE])*self.p*dl.dx - self.dt*dl.exp(self.gamma[DE])*self.u_k[IN]*self.p*dl.dx
            
        for i in range(5):
            self.a[i], self.L[i] = dl.lhs(F[i]), dl.rhs(F[i])
        
    def assemble_systems(self, u_k_vec, u_0_vec):
        
        for i in range(5):
            self.vector2Function(self.u_0[i], u_0_vec[i])
            self.vector2Function(self.u_k[i], u_k_vec[i])

        for i in range(5):
            if not self.H[i] is None:
                if not i == 4:
                    dl.assemble(self.a[i], tensor = self.H[i])
            else:
                self.H[i] = dl.assemble(self.a[i])
            if not self.b[i] is None:
                dl.assemble(self.L[i], tensor = self.b[i])
            else:
                self.b[i] = dl.assemble(self.L[i])
                
        return self.H, self.b
    
    def integrate_over_dist(self, u, i, j):
        
        # need to perform integration over district subdomain
        return self.dist_area_factor[i] * u[j].inner(self.z_dist[i])
        
    def evaluate(self, u):
        int_values = np.empty(3)
        for i in range (3):
            int_values[i] = u[i+2].inner(self.z)
            
        return np.array([np.sum(int_values), int_values[2]])
    
    def evaluate_district(self, u):
        
        # data info:
        # vector of qoi where each row correspond to district
        # last row correspond to total cases
        # 
        # Each row has two elements, one for total infected and other for deceased
        
        int_values = np.zeros((26, 3))
        
        # compute total qoi
        for i in range(3):
            int_values[25, i] = u[i+2].inner(self.z)
           
        # compute district qoi
        for i in range (25):
            for j in range(3):
                # compute integral of species j over district subdomain i
                int_values[i, j] = self.integrate_over_dist(u, i, j+2)
                
        # collect total infected and deceased qoi
        int_dist_qoi = np.zeros((26, 2))
        for i in range(26):
            int_dist_qoi[i, 1] = int_values[i, 2]
            for j in range(3):
                int_dist_qoi[i, 0] += int_values[i, j]
                
        return int_dist_qoi
        
        
    def save(self, u, time):
        if self._file is None:
            raise Exception("Files for saving are not created.")
        for i in range (5):
            self.vector2Function(self.help, u[i])
            self._file[i] << (self.help, time)
        
