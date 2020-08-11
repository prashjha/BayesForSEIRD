import dolfin as dl
import numpy as np
import math
from datetime import datetime

SU = 0
EX = 1
IN = 2
RE = 3
DE = 4

def validate_date(date_text):
    try:
        datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        raise ValueError("Incorrect data format, should be YYYY-MM-DD")

class seird_ic(dl.UserExpression):
    def __init__(self, geom, pop, **kwargs):
        super().__init__(**kwargs)
        if not geom.shape[0] == pop.size:
            raise Exception("The size of the arrays are not matching.")
        self._centroid = geom[:, 0:2]
        self._area = geom[:, 2]
        self._population = pop
        
    def compute_value(self, centroid, area, pop, x):
    
        return 2*pop/area*math.exp(-2*math.pi/area*((x[0]-centroid[0])**2 + (x[1]-centroid[1])**2))

    def eval(self, value, x):
        value[0] = 0.
        for i in range(self._area.size):
            if round(self._population[i]) > 0:
                value[0] += self.compute_value(self._centroid[i,:], self._area[i], self._population[i]/10000., x)
                    
    def value_shape(self):
        return ()
        
def generate_initial_condition(mesh, county_geom, total_pop, infected_cases, recovered_cases, deceased_cases, infected_total_state, deceased_state, save_path = './', transmit_rate = 1.16, date = None, FE_polynomial = 1, save_numpy = True, save_vtu = True):
    
    # Define function space and all functions
    Vu = dl.FunctionSpace(mesh, "Lagrange", FE_polynomial)
    u = [dl.Function(Vu) for i in range (5)]

    
    # Determine the day index based on the given date
    if date is None:
        day_index = 0
    else:
        validate_date(date)
        d0 = datetime.strptime("2020-03-06", "%Y-%m-%d")
        d1 = datetime.strptime(date, "%Y-%m-%d")
        day_index = abs((d1-d0)).days
    
    # Determine the infected and deseased data on the given date
    infected_pop = infected_cases[day_index, :]
    deceased_pop = deceased_cases[day_index, :]
    recovered_pop = recovered_cases[day_index,:]
    active_infected_pop = infected_pop - deceased_pop - recovered_pop
    susceptible_pop = total_pop - transmit_rate*active_infected_pop - deceased_pop - recovered_pop- active_infected_pop
    
    # Interpolate initial conditions of the deceased/recovered/infected density based on the total cases
    deceased_pop_ic = seird_ic(county_geom, deceased_pop)
    u[DE] = dl.interpolate(deceased_pop_ic, Vu)
    recovered_pop_ic = seird_ic(county_geom, recovered_pop)
    u[RE] = dl.interpolate(recovered_pop_ic, Vu)
    active_infected_pop_ic = seird_ic(county_geom, active_infected_pop)
    u[IN] = dl.interpolate(active_infected_pop_ic, Vu)
    exposed_pop_ic = seird_ic(county_geom, transmit_rate*active_infected_pop)
    u[EX] = dl.interpolate(exposed_pop_ic, Vu)
    susceptible_pop_ic = seird_ic(county_geom, susceptible_pop)
    u[SU] = dl.interpolate(susceptible_pop_ic, Vu)
    
    r1 = deceased_state[day_index]/(10000*dl.assemble(u[DE]*dl.dx))
    print(r1)
    u[DE] = dl.project(r1*u[DE], Vu)
    
    r2 = (infected_total_state[day_index] - deceased_state[day_index])/(10000*dl.assemble((u[IN]+u[RE])*dl.dx))
    print(r2)
    u[IN] = dl.project(r2*u[IN], Vu)
    u[RE] = dl.project(r2*u[RE], Vu)
    
    names = ["susceptible", "exposed", "infected", "recovered", "deceased"]
    
    for i in range(5):
        print(names[i] + ": ", dl.assemble(10000*u[i]*dl.dx))
    if save_numpy:
        for i in range(5):
            np.save(save_path + names[i] + '_ic.npy', u[i].vector().get_local())
    if save_vtu:
        for i in range(5):
            file = dl.File(save_path + names[i] + '_ic.pvd', "compressed")
            file << u[i]
