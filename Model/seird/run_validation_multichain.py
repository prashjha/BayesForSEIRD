import dolfin as dl
import numpy as np
import hippylib as hl
from datetime import datetime
from seird_problem import seird_fwd_problem
from seird_initial_with_recovered import generate_initial_condition
from seird_misfit_validation import seird_misfit
from mcmc_tracer import FullTracer

import os
import argparse

STATE = 0
PARAMETER = 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Model SEIRD')
    parser.add_argument('--nsamples',
                        default=5000,
                        type=int,
                        help="Number of MCMC samples")
    # AMAL: start
    parser.add_argument('--chain_id',
                        default=1,
                        type=int,
                        help="ID of MCMC chain")

    parser.add_argument('--sim_time',
                        default=30,
                        type=float,
                        help="Total simulation time")
    
    parser.add_argument('--time_step',
                        default=0.1,
                        type=float,
                        help="Size of time step")
                        
    parser.add_argument('--data_path',
                        default='../../Data/',
                        type=str,
                        help="Relative path to data")
                        
    parser.add_argument('--mesh_path',
                        default='../../Data/mesh/',
                        type=str,
                        help="Relative path to mesh")

    parser.add_argument('--mesh_file',
                        default='mesh_5h',
                        type=str,
                        help="Mesh filename such as mesh_5h or mesh_4h")
    
    parser.add_argument('--ic_path',
                        default='../ic_result/',
                        type=str,
                        help="Relative path to initial condition files")

    parser.add_argument('--generate_ic',
                        default=0,
                        type=int,
                        help="Generate initial condition 0 - false, 1 - true")

    parser.add_argument('--subdmn_path',
                        default='../subdomains/',
                        type=str,
                        help="Relative path to subdomain files")
    
    parser.add_argument('--out_path',
                        default='./',
                        type=str,
                        help="Relative path to store output files")

    parser.add_argument('--samples_path',
                        default='../../Results/calibration_results/',
                        type=str,
                        help="Relative path to parameter samples")

    parser.add_argument('--samples_file',
                        default='calibration_posterior_samples.npy',
                        type=str,
                        help="File name containing samples")

    parser.add_argument('--noise_inf',
                        default=0.05,
                        type=float,
                        help="Noise in the infected data")

    parser.add_argument('--noise_dec',
                        default=0.02,
                        type=float,
                        help="Noise in the deceased data")

    parser.add_argument('--pcn_s',
                        default=0.3,
                        type=float,
                        help="PCN s parameter")
    
    args = parser.parse_args()
    print('Arguments passed: ')
    print(args)
    try:
        dl.set_log_active(False)
    except:
        pass

    # set initial day of simulation
    date = "2020-06-01"

    # load mesh
    mesh_path = args.mesh_path
    mesh_fname = args.mesh_file
    mesh = dl.Mesh(mesh_path + mesh_fname +  ".xml")

    # FE space
    FE_polynomial = 1
    Vu = dl.FunctionSpace(mesh, "Lagrange", FE_polynomial)

    # read COVID-19 data
    data_path = args.data_path
    infected_total_state = np.loadtxt(data_path + 'covid_7July2020/infected_total_state.txt')
    deceased_state = np.loadtxt(data_path + 'covid_7July2020/deceased_state.txt')
        
    # Define elements: P1 and real number
    P1  = dl.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    R   = dl.FiniteElement("R", mesh.ufl_cell(), 0)

    # Define differen scenario
    param_dim = 7

    # Define the parameter space: R^7
    Vh_PARAMETER = dl.VectorFunctionSpace(mesh, "R", degree=0, dim=param_dim)

    # Define the state space for each species
    Vh = [Vu, Vh_PARAMETER, None]

    # simulation time and time step
    simulation_time = args.sim_time
    dt = args.time_step

    # initial condition
    ic_path = args.ic_path + mesh_fname + '/'
    generate_ic = args.generate_ic
    if generate_ic == 1:
        print("Generating initial conditions...")
        geom_path = data_path + 'geography/'

        county_geom_details = np.loadtxt(geom_path + 'county_geom_details.txt')
        county_population = np.loadtxt(geom_path + 'county_population.txt')
        infected_cases = np.loadtxt(data_path + 'infected_total_county.txt')
        deceased_cases = np.loadtxt(data_path + 'deceased_county.txt')
        recovered_cases = np.loadtxt(data_path + 'recovered_county.txt')

        #Generate initial conditions
        generate_initial_condition(mesh, county_geom_details, county_population, infected_cases, recovered_cases, deceased_cases, infected_total_state, deceased_state, save_path=ic_path, FE_polynomial = FE_polynomial, date = date)

    # read initial condition
    u_init = [None]*5
    names = ["susceptible", "exposed", "infected", "recovered", "deceased"]
    for i in range(5):
        u_init[i] = np.load(ic_path + names[i] + '_ic.npy')

    # describe qoi type: 
    # 'state', 'district'
    qoi_type = 'state'

    # create pde problem
    pde = seird_fwd_problem(Vh, simulation_time, dt, u_init, subdmn_path = args.subdmn_path, mesh_tag=mesh_fname, qoi_type=qoi_type)

    # create misfit function
    misfit = seird_misfit(infected_total_state, deceased_state, date, simulation_time, validation = True, validation_day = 10)
    noise_inf = args.noise_inf
    noise_dec = args.noise_dec
    misfit.set_noise_variance([noise_inf, noise_dec])
    
    # read calibration posterior samples and compute mean and std for Gaussian
    # approximation of posterior
    prior_path = args.samples_path
    prior_file = args.samples_file
    prior_samples = np.load(prior_path + prior_file)
    if prior_samples.shape[0] == param_dim:
        p_mean = np.mean(prior_samples, axis = 1)
        p_cov = np.cov(prior_samples)
    elif prior_samples.shape[1] == param_dim:
        p_mean = np.mean(prior_samples, axis = 0)
        p_cov = np.cov(prior_samples, rowvar = False)
    else:
        raise IndexError('Wrong dimension for prior samples')
    
    print(p_mean)
    print(p_cov)
    
    # base prior for validation on calibration posterior
    mean = pde.generate_parameter()
    mean.set_local(p_mean)
    prior = hl.GaussianRealPrior(Vh[PARAMETER], p_cov, mean=mean)

    # create model
    model = hl.Model(pde, prior, misfit)
    
    chain_id = args.chain_id
    
    kernel = hl.pCNKernel(model)
    pcn_s = args.pcn_s
    kernel.parameters["s"] = pcn_s
    chain = hl.MCMC(kernel)
    
    # output path
    out_path = args.out_path
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    
    # save data
    validation_result_path = out_path + '/run_' + str(chain_id) + '/'
    try:
        np.save(validation_result_path + 'data.npy', misfit.data)
    except:
        os.mkdir(validation_result_path)
        np.save(validation_result_path + 'data.npy', misfit.data)
    
    tracer = FullTracer(Vh, simulation_time, validation_result_path, print = False)

    chain.parameters["number_of_samples"]     = args.nsamples
    chain.parameters["burn_in"]               = 0
    chain.parameters["print_progress"]        = 10
    chain.parameters["print_level"]           = -1
    
    for idx in range(chain_id-1):
        chain.consume_random() # exhaust randoms used in previous chains
        
    noise = dl.Vector()
    prior.init_vector(noise,"noise")
    hl.parRandom.normal(1., noise)
    m0 = dl.Vector()
    prior.init_vector(m0, 0)
    prior.sample(noise,m0)

    n_accept = chain.run(m0, qoi = None, tracer = tracer)
    print("Number accepted = {0}".format(n_accept))
    tracer.save()
