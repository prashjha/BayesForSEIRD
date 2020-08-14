import dolfin as dl
import numpy as np
from datetime import datetime
from seird_problem import seird_fwd_problem
from seird_initial_with_recovered import generate_initial_condition

import os
import argparse

import sys
import timeit

STATE = 0
PARAMETER = 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' Model SEIRD')

    parser.add_argument('--num_sim',
                        default=1,
                        type=int,
                        help="Total number of sims")

    parser.add_argument('--sim_id',
                        default=0,
                        type=int,
                        help="Simulation id")
    
    parser.add_argument('--sim_time',
                        default=80,
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
                        default='../ic_result_15July2020/',
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
                        default='../../Results/calibration_11August2020/',
                        type=str,
                        help="Relative path to parameter samples")

    parser.add_argument('--samples_file',
                        default='param_samples_mcmc_10.npy',
                        type=str,
                        help="File name containing samples")

     
    args = parser.parse_args()
    print('Arguments passed: ')
    print(args)
    try:
        dl.set_log_active(False)
    except:
        pass

    # set initial day of simulation
    date = "2020-07-15"

    # load mesh
    mesh_path = args.mesh_path
    mesh_fname = args.mesh_file
    mesh = dl.Mesh(mesh_path + mesh_fname +  ".xml")

    # FE space
    FE_polynomial = 1
    Vu = dl.FunctionSpace(mesh, "Lagrange", FE_polynomial)

    # read COVID-19 data
    data_path = args.data_path
    infected_total_state = np.loadtxt(data_path + 'covid_11August2020/infected_total_state.txt')
    deceased_state = np.loadtxt(data_path + 'covid_11August2020/deceased_state.txt')

    # save data 
    d0 = datetime.strptime("2020-07-15", "%Y-%m-%d")
    d1 = datetime.strptime(date, "%Y-%m-%d")
    day_index = abs((d1-d0)).days

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
        infected_cases = np.loadtxt(data_path + 'covid_11August2020/infected_total_county.txt')
        deceased_cases = np.loadtxt(data_path + 'covid_11August2020/deceased_county.txt')
        recovered_cases = np.loadtxt(data_path + 'covid_11August2020/recovered_county.txt')

        #Generate initial conditions
        generate_initial_condition(mesh, county_geom_details, county_population, infected_cases, recovered_cases, deceased_cases, infected_total_state, deceased_state, save_path=ic_path, FE_polynomial = FE_polynomial)

    # read initial condition
    u_init = [None]*5
    names = ["susceptible", "exposed", "infected", "recovered", "deceased"]
    for i in range(5):
        u_init[i] = np.load(ic_path + names[i] + '_ic.npy')

    # describe qoi type: 
    # 'state', 'district'
    qoi_type = 'district'

    # create pde problem
    pde = seird_fwd_problem(Vh, simulation_time, dt, u_init, subdmn_path = args.subdmn_path, mesh_tag=mesh_fname, qoi_type=qoi_type)

    # load parameter samples
    samples_path = args.samples_path
    samples_file = args.samples_file
    param_samples = np.load(samples_path + samples_file)
    #param_samples = param_samples[:5]

    # get sims and sim id
    num_sim = args.num_sim+1
    sim_id = args.sim_id
    if num_sim == 0 or num_sim == sim_id:
        raise IndexError('Number of sim or sim id is invalid')

    # output and log paths
    out_path = args.out_path
    qoi_path = out_path + '/qoi/'
    log_path = out_path + '/log/'
    sys.stdout = open(log_path + 'sim_' + str(sim_id) + '.log', 'w')
    print('\n\nsim id: {0:<8} num sim: {1}\n\n'.format(sim_id, num_sim))

    # save real data
    if sim_id == 0:
        data_last_day = 25
        data = np.empty((data_last_day, 2))
        data[:, 0] = infected_total_state[day_index+1:(day_index + data_last_day+1):1]
        data[:, 1] = deceased_state[day_index+1:(day_index + data_last_day+1):1]
    
        np.save(qoi_path + 'data.npy', data)

    # solve
    sim_time = []
    qoi = []
    counter = 0
    current_run = 0

    for param in param_samples:

        # divide the tasks into num_sim and let sim_id perform task 
        # associated to it
        if counter % num_sim == sim_id:

            print('\n\nsimulation counter: {}\n'.format(counter))

            start = timeit.default_timer()

            #### solve for current sample
            x = [pde.generate_state_district(), pde.generate_parameter(), None]
            m_test = np.array(param)

            x[PARAMETER].set_local(m_test)
            pde.solveFwd(x[STATE], x)

            # save QoI
            qoi.append(x[STATE])
            if current_run % 5 == 0:
                np.save(qoi_path + 'qoi_' + str(sim_id) + '.npy', qoi)

            stop = timeit.default_timer()
            sim_time.append(stop - start)

            current_run += 1
            
        counter += 1


    # save QoI
    np.save(qoi_path + 'qoi_' + str(sim_id) + '.npy', qoi)

    # output mean time
    print('mean sim time: {}'.format(np.mean(sim_time)))


