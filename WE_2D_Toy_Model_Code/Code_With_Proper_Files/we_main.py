import sys
main_directory = '/scratch/users/sahn1/WE_2D_Toy_Model'  # TODO: set main directory for WE simulation
sys.path.append(main_directory)
import os
os.chdir(main_directory)
import numpy as np
import we_global_variables as gv
import we_functions
from time import time


def weighted_ensemble_simulation(input_parameters_file, input_initial_values_file):
    # set simulation parameters
    we_functions.set_parameters(input_parameters_file)

    # create python objects for walkers and balls
    walker_list = [None]*(gv.max_num_balls*gv.num_walkers)
    temp_walker_list = [None]*(gv.max_num_balls*gv.num_walkers)
    vacant_walker_indices = []
    balls = np.zeros((1, gv.num_cvs+2))
    ball_to_walkers = {}

    # create walkers and their directories
    we_functions.initialize(input_initial_values_file, walker_list, temp_walker_list, ball_to_walkers,
                            vacant_walker_indices)

    for step_num in range(gv.initial_step_num, gv.initial_step_num + gv.max_num_steps):
        # reset ball objects so that balls are newly created at every step
        if gv.balls_flag == 0 and step_num != gv.initial_step_num:
            balls = np.zeros((1, gv.num_cvs+2))
            ball_to_walkers = {}
            gv.current_num_balls = 0

        if gv.flag == 1 and step_num == gv.initial_step_num:
            pass
        else:
            gv.first_walker = 0
            gv.last_walker = gv.num_occupied_balls*gv.num_walkers-1
        print 'running   ' + str(step_num+1) + '-th step'

        # first, run simulation
        t0 = time()
        if (gv.flag == 2 or gv.flag == 3) and step_num == gv.initial_step_num:
            pass
        else:
            we_functions.m_simulation(walker_list)

        # second, create balls and assign walkers to balls
        t1 = time()
        new_balls = we_functions.binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers)

        # third, perform spectral clustering
        if step_num == gv.sc_start:
            gv.transition_matrix = np.zeros((new_balls.shape[0], new_balls.shape[0]))
            gv.ref_balls = new_balls
        if gv.sc_start < step_num <= gv.sc_start+gv.sc_steps:
            we_functions.spectral_clustering(step_num, temp_walker_list, new_balls)

        # fourth, resample walkers for every ball
        we_functions.resampling(walker_list, temp_walker_list, new_balls, ball_to_walkers, vacant_walker_indices)

        # finally, output the results in text files
        we_functions.print_status(step_num, walker_list, new_balls, ball_to_walkers)
        balls = new_balls
        t2 = time()

        os.chdir(gv.main_directory+'/WE')
        f = open('time_record.txt', 'a')
        f.write(str(step_num+1) + '-th step: ' + 'simulation time: ' + str(t1-t0) + ' ' + 'post-processing time: ' +
                str(t2-t1) + '\n')
        f.close()
        
weighted_ensemble_simulation('we_parameters.txt', 'we_initial_values.txt')
