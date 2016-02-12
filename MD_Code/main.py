import numpy as np
from time import time
import os
import sys
main_directory = '/scratch/users/sahn1/Triazine'  # TODO: set main directory for CAS simulation
sys.path.append(main_directory)
os.chdir(main_directory)
import global_variables as gv
import functions


def CAS_simulation(input_initial_values_file):
    # set simulation parameters
    functions.set_parameters()

    # create python objects for walkers and balls
    if gv.enhanced_sampling_flag == 3:
        walker_list = [None]*(gv.num_balls_limit*gv.num_walkers_for_sc*2)
        temp_walker_list = [None]*(gv.num_balls_limit*gv.num_walkers_for_sc*2)
    else:
        walker_list = [None]*(gv.num_balls_limit*gv.num_walkers*2)
        temp_walker_list = [None]*(gv.num_balls_limit*gv.num_walkers*2)
    vacant_walker_indices = []
    balls = np.zeros((1, gv.num_cvs+3))  # ball coordinates / ball radius / ball key / # of walkers
    ball_to_walkers = {}
    key_to_ball = {}
    ball_clusters_list = {}

    # create walkers and their directories
    new_balls = functions.initialize(input_initial_values_file, walker_list, temp_walker_list, balls, ball_to_walkers,
                                     vacant_walker_indices)
    balls = new_balls

    for step_num in range(gv.initial_step_num, gv.initial_step_num + gv.max_num_steps):
        # reset ball objects so that balls are newly created at every step
        if gv.balls_flag == 0 and step_num != gv.initial_step_num:
            balls = np.zeros((1, gv.num_cvs+3))
            ball_to_walkers = {}
            key_to_ball = {}
            ball_clusters_list = {}
            gv.current_num_balls = 0

        if gv.simulation_flag != 0 and step_num == gv.initial_step_num:
            pass
        else:
            gv.first_walker = 0
            gv.last_walker = gv.total_num_walkers-1
        print 'running   ' + str(step_num+1) + '-th step'
        os.chdir(gv.main_directory)
        f = open('bash_script_input_file.txt', 'w')
        f.write(str(gv.first_walker))
        f.write(' first_' + str(gv.last_walker) + '_last')
        f.close()

        # first, run simulation or clean up with bash script
        t0 = time()
        if (gv.simulation_flag == 2 or gv.simulation_flag == 3) and step_num == gv.initial_step_num:
            pass
        #elif gv.simulation_flag == 1 and step_num == gv.initial_step_num:
            #os.system('./clean_up.sh')
        else:
            os.system('./simulations.sh')

        # second, create balls and assign walkers to balls
        t1 = time()
        if gv.enhanced_sampling_flag == 2:
            new_balls = functions.threshold_binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers,
                                                    key_to_ball)
        else:
            new_balls = functions.binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers, key_to_ball)

        # third, perform spectral clustering if enhanced_sampling_flag = 3
        if gv.enhanced_sampling_flag == 3 and gv.num_balls_for_sc <= gv.num_occupied_balls \
                and step_num != gv.initial_step_num and gv.sc_performed == 0:
            functions.spectral_clustering(step_num, temp_walker_list, new_balls,  ball_clusters_list)
            # fourth, resample walkers for every ball
            if gv.sc_performed == 1:
                functions.resampling_for_sc(walker_list, temp_walker_list, new_balls, ball_to_walkers,
                                               ball_clusters_list, vacant_walker_indices)
            else:
                functions.resampling(walker_list, temp_walker_list, new_balls, ball_to_walkers, vacant_walker_indices)
        else:
            functions.resampling(walker_list, temp_walker_list, new_balls, ball_to_walkers, vacant_walker_indices)

        # finally, output the results in text files
        functions.print_status(step_num, walker_list, new_balls, ball_to_walkers, ball_clusters_list, key_to_ball)
        balls = new_balls
        t2 = time()

        os.chdir(gv.main_directory+'/CAS')
        f = open('time_record.txt', 'a')
        f.write(str(step_num+1) + '-th step: ' + 'simulation time: ' + str(t1-t0) + ' ' + 'post-processing time: ' +
                str(t2-t1) + '\n')
        f.close()

CAS_simulation('initial_values.txt')
