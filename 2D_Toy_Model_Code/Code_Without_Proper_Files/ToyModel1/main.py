import numpy as np
from time import time
import os
import sys
main_directory = '/scratch/users/sahn1/2D_Toy_Model'  # TODO: set main directory for CAS simulation
sys.path.append(main_directory)
os.chdir(main_directory)
import global_variables as gv
import functions


def CAS_simulation(input_initial_values_file):
    # set simulation parameters.
    functions.set_parameters()

    # create python objects for walkers and macrostates.
    # walker_list keeps track of previous information whereas temp_walker_list keeps track of current/new information.
    if gv.enhanced_sampling_flag == 2:
        walker_list = [None]*(gv.num_balls_for_sc*gv.num_walkers*10)
        temp_walker_list = [None]*(gv.num_balls_for_sc*gv.num_walkers*10)
    else:
        walker_list = [None]*(gv.num_balls_limit*gv.num_walkers*2)
        temp_walker_list = [None]*(gv.num_balls_limit*gv.num_walkers*2)

    # balls is recorded in the following order: ball coordinates / ball radius / ball key / # of walkers
    balls = np.zeros((1, gv.num_cvs+3))
    balls_array = np.zeros((1, gv.num_cvs))
    # maps ball coordinates to walkers
    ball_to_walkers = {}
    # list of ball clusters (used only when gv.enhanced_sampling_flag = 2)
    ball_clusters_list = {}

    # create walkers and their directories.
    functions.initialize(input_initial_values_file, walker_list)

    for step_num in range(gv.max_num_steps):
        # reset ball objects so that balls are newly created at every step.
        if gv.balls_flag == 0:
            balls = np.zeros((1, gv.num_cvs+3))
            balls_array = np.zeros((1, gv.num_cvs))
            ball_to_walkers = {}
            ball_clusters_list = {}
            gv.current_num_balls = 0

        print 'running ' + str(step_num+1) + '-th step'

        # first, run simulation.
        t0 = time()
        functions.m_simulation(walker_list)

        # second, create macrostates and assign or bin walkers to macrostates.
        t1 = time()
        if gv.enhanced_sampling_flag == 1:
            balls, balls_array = functions.threshold_binning(step_num, walker_list, temp_walker_list, balls,
                                                             balls_array, ball_to_walkers)
        else:
            balls, balls_array = functions.binning(step_num, walker_list, temp_walker_list, balls, balls_array,
                                                   ball_to_walkers)

        t2 = time()
        # third, perform spectral clustering if enhanced_sampling_flag = 2.
        if gv.enhanced_sampling_flag == 2 and gv.num_balls_for_sc <= gv.num_occupied_balls and gv.sc_performed == 0 \
                and gv.sc_start == -1:
            # start fixing macrostates from this point on until we finish calculating the transition matrix
            gv.balls_flag = 1
            gv.sc_start = step_num
        if gv.enhanced_sampling_flag == 2 and gv.sc_performed == 0 and gv.sc_start != -1:
            functions.calculate_trans_mat_for_sc(step_num, temp_walker_list, balls, balls_array)
        if gv.enhanced_sampling_flag == 2 and gv.sc_performed == 1:
            functions.spectral_clustering(step_num, balls, ball_clusters_list)
            # fourth, resample walkers for every macrostate.
            if gv.sc_performed == 1:
                balls = functions.resampling_for_sc(walker_list, temp_walker_list, balls, ball_to_walkers,
                                                    ball_clusters_list)
            else:
                balls = functions.resampling(walker_list, temp_walker_list, balls, ball_to_walkers)
        else:
            balls = functions.resampling(walker_list, temp_walker_list, balls, ball_to_walkers)

        # finally, output the results as text files.
        balls = functions.print_status(step_num, walker_list, balls, ball_to_walkers, ball_clusters_list)
        t3 = time()

        os.chdir(gv.main_directory+'/CAS')
        f = open('time_record.txt', 'a')
        f.write(str(step_num+1) + '-th step: simulation time: ' + str(t1-t0) + ' binning time: ' + str(t2-t1) +
                ' resampling time: ' + str(t3-t2) + '\n')
        f.close()

CAS_simulation('initial_values.txt')
