import copy
import itertools
import os
import shutil
import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError
# from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.covariance import EllipticEnvelope
import check_state_function
import global_variables as gv
import parameters as p
import walker


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(gv.num_cvs):
        if gv.angle_cvs[i] == 0:
            distance += (values[i] - center[i])**2
        else:
            distance += min(360.0 - abs(values[i] - center[i]), abs(values[i] - center[i]))**2
    if abs(distance) < 1.0e-10:
        distance = 0.0
    return np.sqrt(distance)


def closest_ball(walker_coordinates, balls_array):
    distance = np.zeros((balls_array.shape[0],))
    inside = np.zeros((balls_array.shape[0],))
    for i in range(gv.num_cvs):
        if gv.separate_radii_flag == 1:
            radius = gv.radius[i]
        else:
            radius = gv.radius
        if gv.angle_cvs[i] == 0:
            distance_from_center = (balls_array[:, i] - walker_coordinates[i])**2
            inside[distance_from_center <= radius**2] += 1
            distance += distance_from_center
        else:
            distance_from_center = np.minimum(360.0 - np.abs(balls_array[:, i] - walker_coordinates[i]),
                                   np.abs(balls_array[:, i] - walker_coordinates[i]))**2
            inside[distance_from_center <= radius**2] += 1
            distance += distance_from_center
    if gv.separate_radii_flag == 1:
        inside_max = np.max(inside)
        indices = np.where(inside == inside_max)[0]
        index = np.argmin(distance[indices])
        closest_ball_key = indices[index]
    else:
        closest_ball_key = np.argmin(distance)
        if distance[closest_ball_key] <= radius**2 or abs(distance[closest_ball_key] - radius**2) < 1.0e-10:
            inside_max = gv.num_cvs
        else:
            inside_max = 0
    return closest_ball_key, inside_max


def set_parameters():
    gv.main_directory = p.main_directory
    gv.initial_configuration_directory = p.initial_configuration_directory
    gv.simulation_flag = p.simulation_flag
    gv.balls_flag = p.balls_flag
    gv.flux_flag = p.flux_flag
    gv.num_states = p.num_states
    gv.enhanced_sampling_flag = p.enhanced_sampling_flag
    # macrostates can't be fixed if we use threshold binning or spectral clustering
    if (gv.enhanced_sampling_flag == 1 or gv.enhanced_sampling_flag == 2) and gv.balls_flag == 1:
        gv.balls_flag = 0
    gv.num_balls_limit = p.num_balls_limit
    gv.separate_radii_flag = p.separate_radii_flag
    gv.radius = p.radius
    gv.num_walkers = p.num_walkers
    gv.num_cvs = p.num_cvs
    gv.angle_cvs = p.angle_cvs
    gv.initial_step_num = p.initial_step_num
    gv.max_num_steps = p.max_num_steps
    gv.num_occupied_balls = p.num_occupied_balls
    gv.first_walker = p.first_walker
    gv.last_walker = p.last_walker
    if gv.enhanced_sampling_flag == 1:
        gv.less_or_greater_flag = p.less_or_greater_flag
        gv.static_threshold_flag = p.static_threshold_flag
        gv.threshold_values = p.threshold_values
        gv.properties_to_keep_track = p.properties_to_keep_track
    elif gv.enhanced_sampling_flag == 2:
        gv.num_occupied_clusters = p.num_occupied_clusters
        gv.num_balls_for_sc = p.num_balls_for_sc
        gv.num_clusters = p.num_clusters
        gv.num_walkers_for_sc = p.num_walkers_for_sc
        gv.num_steps_for_sc = p.num_steps_for_sc
    elif gv.enhanced_sampling_flag == 3:
        gv.initial_step_num_for_eq = p.initial_step_num_for_eq
        gv.num_steps_for_eq = p.num_steps_for_eq
        gv.eq_frequency = p.eq_frequency
        gv.num_steps_in_bw = p.num_steps_in_bw

    gv.current_num_balls = 0
    # when restarting a simulation, the total number of walkers depends on whether it stopped right after
    # spectral clustering or not
    if gv.enhanced_sampling_flag == 2 and gv.simulation_flag != 0 and gv.num_occupied_clusters != 0:
        gv.total_num_walkers = gv.num_occupied_clusters*gv.num_walkers_for_sc
        gv.num_occupied_clusters = 0
    else:
        if gv.simulation_flag == 0:
            gv.total_num_walkers = gv.num_occupied_balls
        else:
            gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers
    gv.sc_performed = 0
    gv.sc_start = -1
    gv.prev_balls_flag = gv.balls_flag
    gv.resampling_performed = 0


def initialize(input_initial_values_file, walker_list, temp_walker_list, balls, balls_array, ball_to_walkers):
    # first populate walker_list with walker objects. initial values are random.
    for i in range(len(walker_list)):
        walker_list[i] = walker.Walker([-1000.0]*gv.num_cvs, [-1000.0]*gv.num_cvs, i, [-1000.0]*gv.num_cvs,
                                       [-1000.0]*gv.num_cvs, 0, 0, 0, 0.0, -1)

    # new simulation
    if gv.simulation_flag == 0:
        # all walkers have equally divided weights
        initial_weight = 1.0/gv.total_num_walkers
        f = open(input_initial_values_file, 'r')
        if gv.flux_flag == 1:
            flux_f = open('initial_states.txt', 'r')
        # for each occupied ball (usually 1 because one initial state is provided but multiple can be provided)
        for n in range(gv.num_occupied_balls):
            # read initial values from file
            line = f.readline().strip().split()
            initial_values = [float(entry) for entry in line]
            # if fluxes are calculated, obtain initial state
            if gv.flux_flag == 1:
                line = flux_f.readline().strip()
                initial_state = int(line)
            for i in range(n, n+1):
                walker_list[i].set(initial_values, initial_weight)
                if gv.flux_flag == 1:
                    walker_list[i].state = initial_state
        f.close()
        if gv.flux_flag == 1:
            flux_f.close()

        # make walker directories
        os.system('mkdir CAS')
        os.chdir(gv.main_directory + '/CAS')
        for n in range(gv.num_occupied_balls):
            for i in range(n, n+1):
                walker_directory = gv.main_directory + '/CAS/walker' + str(i)
                os.mkdir(walker_directory)
                # TODO: change initial configurations' name
                os.system('cp '+str(gv.initial_configuration_directory)+'/minim_'+str(n)+'.gro '+str(walker_directory)+
                          '/minim.gro')

        # in case created balls are kept throughout simulation
        if gv.balls_flag == 1:
            os.chdir(gv.main_directory)
            created_balls = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')
            balls = np.zeros((created_balls.shape[0], gv.num_cvs+2))
            balls[:, 0:gv.num_cvs] = created_balls[:, 0:gv.num_cvs]
            balls_array = balls[:, 0:gv.num_cvs]
            gv.current_num_balls = balls.shape[0]
            for j in range(balls.shape[0]):
                current_ball_center = balls[j][0:gv.num_cvs].tolist()
                ball_to_walkers[tuple(current_ball_center)] = []
                balls[j][gv.num_cvs] = j

    # restarting simulation in the middle of simulation
    elif gv.simulation_flag == 1 or gv.simulation_flag == 2:
        # loop through all walkers and obtain their weights, coordinates, ball coordinates, state, etc.
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            f = open('trajectory.txt', 'r')
            if gv.initial_step_num > 1:
                lines = f.readlines()[-2:]
                previous_line = lines[0].strip().split()
                previous_coordinates = [float(entry) for entry in previous_line]
                current_line = lines[1].strip().split()
                current_coordinates = [float(entry) for entry in current_line]
            else:
                lines = f.readlines()[-1:]
                previous_line = lines[0].strip().split()
                previous_coordinates = [float(entry) for entry in previous_line]
                current_coordinates = previous_coordinates
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            f.close()
            f = open('ball_trajectory.txt', 'r')
            if gv.initial_step_num > 1:
                lines = f.readlines()[-2:]
                previous_line = lines[0].strip().split()
                previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                current_line = lines[1].strip().split()
                current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
            else:
                lines = f.readlines()[-1:]
                previous_line = lines[0].strip().split()
                previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                current_line = previous_line
                current_ball_center = previous_ball_center
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].previous_ball_key = int(previous_line[gv.num_cvs+1])
            walker_list[i].current_ball_key = int(current_line[gv.num_cvs+1])
            f.close()
            if gv.flux_flag == 1:
                walker_list[i].state = int(current_line[-1])
            # in case created balls are kept throughout simulation
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                created_balls = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')
                balls = np.zeros((created_balls.shape[0], gv.num_cvs+2))
                balls[:, 0:gv.num_cvs] = created_balls[:, 0:gv.num_cvs]
                balls_array = balls[:, 0:gv.num_cvs]
                gv.current_num_balls = balls.shape[0]
                for j in range(balls.shape[0]):
                    current_ball_center = balls[j][0:gv.num_cvs].tolist()
                    ball_to_walkers[tuple(current_ball_center)] = []
                    balls[j][gv.num_cvs] = j

    # restarting simulation in the middle of binning
    elif gv.simulation_flag == 3:
        # loop through all walkers and obtain their weights, coordinates, ball coordinates, state, etc.
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            f = open('trajectory.txt', 'r')
            if gv.initial_step_num > 1:
                lines = f.readlines()[-2:]
                previous_line = lines[0].strip().split()
                previous_coordinates = [float(entry) for entry in previous_line]
                current_line = lines[1].strip().split()
                current_coordinates = [float(entry) for entry in current_line]
            else:
                lines = f.readlines()[-1:]
                previous_line = lines[0].strip().split()
                previous_coordinates = [float(entry) for entry in previous_line]
                current_coordinates = previous_coordinates
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            f.close()
            # if walker is already binned to a ball, delete the binning and have binning start from scratch
            num_lines = sum(1 for line in open('ball_trajectory.txt'))
            if num_lines > gv.initial_step_num:
                os.system('sed -i \'$d\' ball_trajectory.txt')
            f = open('ball_trajectory.txt', 'r')
            if gv.initial_step_num > 1:
                lines = f.readlines()[-2:]
                previous_line = lines[0].strip().split()
                previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                current_line = lines[1].strip().split()
                current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
            else:
                lines = f.readlines()[-1:]
                previous_line = lines[0].strip().split()
                previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                current_line = previous_line
                current_ball_center = previous_ball_center
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].previous_ball_key = int(previous_line[gv.num_cvs+1])
            walker_list[i].current_ball_key = int(current_line[gv.num_cvs+1])
            f.close()
            if gv.flux_flag == 1:
                walker_list[i].state = int(current_line[-1])
            # in case created balls are kept throughout simulation
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                created_balls = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')
                balls = np.zeros((created_balls.shape[0], gv.num_cvs+2))
                balls[:, 0:gv.num_cvs] = created_balls[:, 0:gv.num_cvs]
                balls_array = balls[:, 0:gv.num_cvs]
                gv.current_num_balls = balls.shape[0]
                for j in range(balls.shape[0]):
                    current_ball_center = balls[j][0:gv.num_cvs].tolist()
                    ball_to_walkers[tuple(current_ball_center)] = []
                    balls[j][gv.num_cvs] = j

    # restarting simulation in the middle of resampling
    # this is more tricky than the previous cases, since some walker directories are entirely gone or partially gone
    # and so we have to get the remaining walkers and make new walkers that will carry the remaining weights
    # that have disappeared from this process. the simulation will restart from the binning stage and minor differences
    # caused from this temporary fix will disappear after running a few subsequent steps of the simulation.
    elif gv.simulation_flag == 4:
        total_weight = 0.0
        total_num_walkers = 0
        if gv.enhanced_sampling_flag == 2:
            occupied_indices = np.zeros(gv.num_balls_for_sc*gv.num_walkers*100, int)
        else:
            occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers*2, int)
        vacant_walker_indices = []
        # first, obtain previous macrostates' weights to keep track of which weights are missing and keep track of
        # walkers to the previous macrostates.
        previous_ball_to_walkers = {}
        os.chdir(gv.main_directory + '/CAS')
        previous_balls_weights = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')

        # second, loop through all of the remaining walkers and check to see if all of the files exist or not.
        for i in range(gv.last_walker+1):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            # TODO: change simulation files' names
            if os.path.isfile(walker_directory + '/weight_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/ball_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/trajectory.txt') and \
                    os.path.isfile(walker_directory + '/traj.xtc') and \
                    os.path.isfile(walker_directory + '/minim.gro') and \
                    os.path.isfile(walker_directory + '/minim.tpr'):
                os.chdir(walker_directory)
                f = open('weight_trajectory.txt', 'r')
                weight = float(f.readlines()[-1].strip())
                f.close()
                f = open('ball_trajectory.txt', 'r')
                if gv.initial_step_num > 1:
                    lines = f.readlines()[-2:]
                    previous_line = lines[0].strip().split()
                    previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                    previous_ball_key = int(previous_line[gv.num_cvs])
                    current_line = lines[1].strip().split()
                    current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
                    current_ball_key = int(current_line[gv.num_cvs])
                    current_state = int(current_line[gv.num_cvs+1])
                else:
                    lines = f.readlines()[-1:]
                    previous_line = lines[0].strip().split()
                    previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
                    previous_ball_key = int(previous_line[gv.num_cvs])
                    current_line = previous_line
                    current_ball_center = previous_ball_center
                    current_ball_key = previous_ball_key
                    current_state = int(current_line[gv.num_cvs+1])
                f.close()
                # if there's still weight left over in the macrostate of interest,
                # obtain the walker's weights, coordinates, ball coordinates, state, etc.
                if previous_balls_weights[previous_ball_key][-1]-weight >= 0.0:
                    walker_list[i].weight = weight
                    total_weight += weight
                    previous_balls_weights[previous_ball_key][-1] -= weight
                    if tuple(previous_ball_center) in previous_ball_to_walkers:
                        previous_ball_to_walkers[tuple(previous_ball_center)].append(i)
                    else:
                        previous_ball_to_walkers[tuple(previous_ball_center)] = [i]
                    if tuple(current_ball_center) in ball_to_walkers:
                        ball_to_walkers[tuple(current_ball_center)].append(i)
                    else:
                        ball_to_walkers[tuple(current_ball_center)] = [i]
                    gv.current_num_balls += 1
                    walker_list[i].previous_ball_center = previous_ball_center
                    walker_list[i].current_ball_center = current_ball_center
                    f = open('trajectory.txt', 'r')
                    if gv.initial_step_num > 1:
                        lines = f.readlines()[-2:]
                        previous_line = lines[0].strip().split()
                        previous_coordinates = [float(entry) for entry in previous_line]
                        current_line = lines[1].strip().split()
                        current_coordinates = [float(entry) for entry in current_line]
                    else:
                        lines = f.readlines()[-1:]
                        previous_line = lines[0].strip().split()
                        previous_coordinates = [float(entry) for entry in previous_line]
                        current_coordinates = previous_coordinates
                    f.close()
                    walker_list[i].previous_coordinates = previous_coordinates
                    walker_list[i].current_coordinates = current_coordinates
                    walker_list[i].state = current_state
                    walker_list[i].previous_ball_key = previous_ball_key
                    walker_list[i].current_ball_key = current_ball_key
                    temp_walker_list[i] = walker.Walker(previous_coordinates, current_coordinates, i,
                                                        previous_ball_center, current_ball_center, previous_ball_key,
                                                        current_ball_key, 0, weight, current_state)
                    total_num_walkers += 1
                    occupied_indices[i] = 1
                # if not, then delete the walker and potentially use the walker for other remaining weights.
                else:
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('rm -rf walker' + str(i))
                    vacant_walker_indices.append(i)
            # otherwise, then delete the incomplete walker and potentially use the walker for other remaining weights.
            else:
                if os.path.isdir(walker_directory):
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('rm -rf walker' + str(i))
                vacant_walker_indices.append(i)

        # third, create new walkers for the remaining weights by looping through the remaining weights.
        excess_index = gv.last_walker+1
        for i in range(previous_balls_weights.shape[0]):
            if previous_balls_weights[i][-1] > 0.0:
                previous_ball_center = previous_balls_weights[i][0:gv.num_cvs].tolist()
                # if a reference walker exists, copy that walker and put the remaining weight.
                if tuple(previous_ball_center) in previous_ball_to_walkers:
                    reference_walker = previous_ball_to_walkers[tuple(previous_ball_center)][0]
                # if not, then use the walker that's closest to the macrostate of interest as the reference walker.
                else:
                    distance = 0.0
                    reference_ball_center = previous_ball_center
                    first = 0
                    for ball_center in previous_ball_to_walkers:
                        distance_from_center = calculate_distance_from_center(ball_center, previous_ball_center)
                        if first == 0:
                            distance = distance_from_center
                            reference_ball_center = ball_center
                            first += 1
                        else:
                            if distance_from_center < distance:
                                distance = distance_from_center
                                reference_ball_center = ball_center
                    reference_walker = previous_ball_to_walkers[reference_ball_center][0]
                reference_walker_directory = gv.main_directory + '/CAS/walker' + str(reference_walker)
                if len(vacant_walker_indices) > 0:
                    walker_index = vacant_walker_indices.pop(0)
                else:
                    walker_index = excess_index
                    excess_index += 1
                walker_directory = gv.main_directory + '/CAS/walker' + str(walker_index)
                shutil.copytree(reference_walker_directory, walker_directory)
                weight = previous_balls_weights[i][-1]
                previous_balls_weights[i][-1] -= weight
                os.chdir(walker_directory)
                f = open('weight_trajectory.txt', 'w')
                f.write(str(weight) + '\n')
                walker_list[walker_index].weight = weight
                total_weight += weight
                f.close()
                f = open('ball_trajectory.txt', 'r')
                lines = f.readlines()[-1:]
                current_line = lines[0].strip().split()
                current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
                current_ball_key = int(current_line[gv.num_cvs])
                current_state = int(current_line[gv.num_cvs+1])
                f.close()
                walker_list[walker_index].previous_ball_center = previous_ball_center
                walker_list[walker_index].current_ball_center = current_ball_center
                ball_to_walkers[tuple(current_ball_center)].append(walker_index)
                os.system('sed -i \'$d\' ball_trajectory.txt')
                f = open('trajectory.txt', 'r')
                if gv.initial_step_num > 1:
                    lines = f.readlines()[-2:]
                    previous_line = lines[0].strip().split()
                    previous_coordinates = [float(entry) for entry in previous_line]
                    current_line = lines[1].strip().split()
                    current_coordinates = [float(entry) for entry in current_line]
                else:
                    lines = f.readlines()[-1:]
                    previous_line = lines[0].strip().split()
                    previous_coordinates = [float(entry) for entry in previous_line]
                    current_coordinates = previous_coordinates
                f.close()
                walker_list[walker_index].previous_coordinates = previous_coordinates
                walker_list[walker_index].current_coordinates = current_coordinates
                walker_list[walker_index].state = current_state
                walker_list[walker_index].current_ball_key = current_ball_key
                temp_walker_list[walker_index] = walker.Walker(previous_coordinates, current_coordinates, walker_index,
                                                               previous_ball_center, current_ball_center,
                                                               previous_ball_key, current_ball_key, 0, weight,
                                                               current_state)
                total_num_walkers += 1
                occupied_indices[walker_index] = 1

        # finally, re-index the walkers so that the walkers have indices in order from 0 to total_num_walkers-1
        gv.total_num_walkers = excess_index-1
        if total_num_walkers >= gv.total_num_walkers:
            for i in range(total_num_walkers, excess_index):
                new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/CAS')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))
        else:
            for i in range(gv.total_num_walkers, excess_index):
                new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/CAS')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))
            for i in range(total_num_walkers, gv.total_num_walkers):
                if occupied_indices[i] == 1:
                    new_index = vacant_walker_indices.pop()
                    while new_index >= total_num_walkers:
                        new_index = vacant_walker_indices.pop()
                    occupied_indices[new_index] = 1
                    walker_list[new_index].copy_walker(walker_list[i])
                    # rename the directory with name 'i' to 'new_index'
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('mv walker' + str(i) + ' walker' + str(new_index))
        # reset the number of balls since the simulation will restart from the binning step
        gv.total_num_walkers = total_num_walkers
        gv.num_occupied_balls = 0
        gv.current_num_balls = 0
        gv.first_walker = 0
        gv.last_walker = total_num_walkers-1
        # in case created balls are kept throughout simulation
        if gv.balls_flag == 1:
            os.chdir(gv.main_directory + '/CAS')
            created_balls = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')
            balls = np.zeros((created_balls.shape[0], gv.num_cvs+2))
            balls[:, 0:gv.num_cvs] = created_balls[:, 0:gv.num_cvs]
            balls_array = balls[:, 0:gv.num_cvs]
            gv.current_num_balls = balls.shape[0]
            for j in range(balls.shape[0]):
                current_ball_center = balls[j][0:gv.num_cvs].tolist()
                ball_to_walkers[tuple(current_ball_center)] = []
                balls[j][gv.num_cvs] = j
    return balls, balls_array


def binning(step_num, walker_list, temp_walker_list, balls, balls_array, ball_to_walkers):
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)  # convert from list to array
    walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights
    flux = np.zeros((gv.num_states, gv.num_states))
    start = 0  # indicates whether we are dealing with the very first walker or not

    # loop through all of the walkers in descending order based on their weights
    for i in walker_indices:
        # first, go to walker directory i.
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)

        # second, obtain new coordinates' values.
        if os.path.exists(walker_directory + '/coordinates.out'):
            coordinates = np.loadtxt('coordinates.out')
            if gv.num_cvs > 1:
                new_coordinates = coordinates.tolist()
            else:
                new_coordinates = [float(coordinates)]
            rm_command = 'rm -rf *.out'
            os.system(rm_command)

            # also, record the new coordinates' values on the trajectory file
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()
        # if new coordinates' values have been already recorded, then read the last line of the trajectory file.
        else:
            f = open('trajectory.txt', 'r')
            new_coordinates = f.readlines()[-1].strip().split()
            new_coordinates = [float(coordinate) for coordinate in new_coordinates]
            f.close()

        # third, obtain previous information from walker_list[i].
        previous_coordinates = walker_list[i].current_coordinates
        previous_ball_center = walker_list[i].current_ball_center
        previous_ball_key = walker_list[i].current_ball_key
        previous_state = walker_list[i].state
        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight

        # calculate fluxes if needed.
        if gv.flux_flag == 1:
            current_state = check_state_function.check_state_function(new_coordinates)
            if previous_state != -1 and current_state == -1:
                current_state = previous_state
            if previous_state != -1 and current_state != -1:
                flux[previous_state, current_state] += weight
        else:
            current_state = -1

        inside = 0  # indicates whether walker is inside an existing macrostate or not, i.e., binned to a macrostate
        # fourth, bin walker to a macrostate.
        # if we're dealing with the very first walker, create the very first macrostate for the walker.
        if start == 0 and gv.balls_flag == 0:
            start += 1
            inside += 1
            current_ball_center = [float(coordinate) for coordinate in new_coordinates]
            center_key_num = copy.deepcopy(current_ball_center)
            balls_array[gv.current_num_balls] = np.asarray(center_key_num)
            center_key_num.append(gv.current_num_balls)
            center_key_num.append(1)
            balls[gv.current_num_balls] = np.asarray(center_key_num)
            ball_to_walkers[tuple(current_ball_center)] = [i]
            temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                current_ball_center, previous_ball_key, gv.current_num_balls,
                                                initial_step_num, weight, current_state)
            gv.current_num_balls += 1

        # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
        if inside == 0:
            current_ball_key, inside = closest_ball(new_coordinates, balls_array)

            # case 1: walker is inside some macrostate or is not but needs to be binned to the nearest macrostate due to
            # reaching the maximum number of macrostates limit and/or balls_flag = 1.
            if inside == gv.num_cvs or (inside != gv.num_cvs and (gv.current_num_balls == gv.num_balls_limit or
                                                                          gv.balls_flag == 1)):
                balls[current_ball_key][gv.num_cvs+1] += 1
                current_ball_center = balls[current_ball_key][0:gv.num_cvs].tolist()
                ball_to_walkers[tuple(current_ball_center)].append(i)
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                    current_ball_center, previous_ball_key, current_ball_key,
                                                    initial_step_num, weight, current_state)

            # case 2: walker is not inside any macrostate and the maximum number of macrostates limit has not been
            # reached, so create a new macrostate centered around the walker.
            else:
                current_ball_center = [float(coordinate) for coordinate in new_coordinates]
                center_key_num = copy.deepcopy(current_ball_center)
                balls_array = np.append(balls_array, [np.asarray(center_key_num)], axis=0)
                center_key_num.append(gv.current_num_balls)
                center_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_key_num)], axis=0)
                ball_to_walkers[tuple(current_ball_center)] = [i]
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                    current_ball_center, previous_ball_key, gv.current_num_balls,
                                                    initial_step_num, weight, current_state)
                gv.current_num_balls += 1

    # fifth, loop through all of the walkers once more to assign them to their true nearest macrostates
    if gv.balls_flag == 0:
        for i in walker_indices:
            current_coordinates = temp_walker_list[i].current_coordinates
            new_ball_key, inside = closest_ball(current_coordinates, balls_array)
            if inside == gv.num_cvs or (inside != gv.num_cvs and (gv.current_num_balls == gv.num_balls_limit or gv.balls_flag == 1)):
                new_ball_center = balls[new_ball_key][0:gv.num_cvs].tolist()
                old_ball_key = temp_walker_list[i].current_ball_key
                old_ball_center = temp_walker_list[i].current_ball_center
                balls[old_ball_key][gv.num_cvs+1] -= 1
                balls[new_ball_key][gv.num_cvs+1] += 1
                ball_to_walkers[tuple(old_ball_center)].remove(i)
                ball_to_walkers[tuple(new_ball_center)].append(i)
                temp_walker_list[i].current_ball_key = new_ball_key
                temp_walker_list[i].current_ball_center = new_ball_center

    # sixth, record the new macrostate on the ball trajectory file.
    for i in walker_indices:
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)
        current_ball_center = temp_walker_list[i].current_ball_center
        current_ball_key = temp_walker_list[i].current_ball_key
        current_state = temp_walker_list[i].state
        center_key_state = copy.deepcopy(current_ball_center)
        center_key_state.append(current_ball_key)
        center_key_state.append(current_state)
        f = open('ball_trajectory.txt', 'a')
        f.write(' '.join(map(lambda coordinate: str(coordinate), center_key_state)))
        f.write('\n')
        f.close()

    # seventh, delete empty macrostates
    if gv.balls_flag == 0:
        delete_list = []
        for i in range(balls.shape[0]):
            if balls[i][gv.num_cvs+1] == 0:
                delete_list.append(i)
        balls = np.delete(balls, delete_list, 0)
        balls_array = np.delete(balls_array, delete_list, 0)

    # output the total flux for this particular step to a text file, if needed.
    if gv.flux_flag == 1:
        os.chdir(gv.main_directory + '/CAS')
        np.savetxt('flux_' + str(step_num+1) + '.txt', flux, fmt=' %1.5e')

    if gv.balls_flag == 1 and gv.enhanced_sampling_flag == 0:
        # output the transition matrix for this particular step.
        transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
        for i in range(gv.total_num_walkers):
            previous_coordinates = temp_walker_list[i].previous_coordinates
            previous_ball_key, inside = closest_ball(previous_coordinates, balls_array)
            transition_matrix[previous_ball_key][temp_walker_list[i].current_ball_key] += temp_walker_list[i].weight
        os.chdir(gv.main_directory + '/CAS')
        np.savetxt('transition_matrix_' + str(step_num+1) + '.txt', transition_matrix, fmt=' %1.10e')
    return balls, balls_array


def threshold_binning(step_num, walker_list, temp_walker_list, balls, balls_array, ball_to_walkers):
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)  # convert from list to array
    walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights
    flux = np.zeros((gv.num_states, gv.num_states))

    # if threshold values change throughout the simulation, the following objects are needed.
    if gv.static_threshold_flag == 0:
        new_threshold_values = gv.threshold_values
        ref_walker_binning_value = 0
        ref_walker_properties_value = 0.0
        ref_walker_index = 0
    # otherwise, the list of "leftover" walkers, i.e., walkers that did not meet the threshold requirement, is needed.
    else:
        leftover_walker_indices = []

    # loop through all of the walkers in descending order based on their weights.
    for i in walker_indices:
        # go to walker directory i.
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)

        # obtain new coordinates' values.
        if os.path.exists(walker_directory + '/coordinates.out'):
            coordinates = np.loadtxt('coordinates.out')
            if gv.num_cvs > 1:
                new_coordinates = coordinates.tolist()
            else:
                new_coordinates = [float(coordinates)]
            rm_command = 'rm -rf *.out'
            os.system(rm_command)

            # also, record the new coordinates' values on the trajectory file
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()
        # if new coordinates' values have been already recorded, then read the last line of the trajectory file.
        else:
            f = open('trajectory.txt', 'r')
            new_coordinates = f.readlines()[-1].strip().split()
            new_coordinates = [float(coordinate) for coordinate in new_coordinates]
            f.close()

        # third, obtain previous information from walker_list[i].
        previous_coordinates = walker_list[i].current_coordinates
        previous_ball_center = walker_list[i].current_ball_center
        previous_ball_key = walker_list[i].current_ball_key
        previous_state = walker_list[i].state
        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight

        # calculate fluxes if needed.
        if gv.flux_flag == 1:
            current_state = check_state_function.check_state_function(new_coordinates)
            if previous_state != -1 and current_state == -1:
                current_state = previous_state
            if previous_state != -1 and current_state != -1:
                flux[previous_state, current_state] += weight
        else:
            current_state = -1

        temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                            previous_ball_center, previous_ball_key, gv.current_num_balls,
                                            initial_step_num, weight, current_state)

        # if threshold values change throughout the simulation, the walker with the lowest or highest value
        # needs to be found.
        if gv.static_threshold_flag == 0:
            # first, find out which properties need to be kept track of.
            properties_to_keep_track = []
            for k in range(len(gv.properties_to_keep_track)):
                if gv.properties_to_keep_track[k] < 0:
                    properties_to_keep_track.append(weight)
                else:
                    properties_to_keep_track.append(new_coordinates[gv.properties_to_keep_track[k]])
            # second, find out if walker needs to be binned to one designated "leftover" macrostate
            # (which in this case, the walker_binning_value goes up)
            # and the degree to which this walker needs to be binned to it
            # (the larger the degree, the larger the walker_properties_value).
            walker_binning_value = 0
            walker_properties_value = 0.0
            if gv.less_or_greater_flag == 0:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] < gv.threshold_values[m]:
                        walker_binning_value += 1
                    walker_properties_value += (gv.threshold_values[m]-properties_to_keep_track[m])
            else:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        walker_binning_value += 1
                    walker_properties_value += (properties_to_keep_track[m]-gv.threshold_values[m])
            if i == 0:  # the very first walker
                ref_walker_index = i
                ref_walker_binning_value = walker_binning_value
                ref_walker_properties_value = walker_properties_value
            else:  # otherwise, check to see if this walker can replace the reference walker
                if walker_binning_value == ref_walker_binning_value \
                        and walker_properties_value <= ref_walker_properties_value:
                    ref_walker_index = i
                    ref_walker_binning_value = walker_binning_value
                    ref_walker_properties_value = walker_properties_value
                elif walker_binning_value < ref_walker_binning_value:
                    ref_walker_index = i
                    ref_walker_binning_value = walker_binning_value
                    ref_walker_properties_value = walker_properties_value
        # if threshold values do not change throughout the simulation,
        # simply check and see if the walker needs to be binned to one designated "leftover" macrostate
        # (there's no need to check the degree to which the walker needs to be binned,
        # since that's only relevant for getting the walker with the lowest or highest value).
        else:
            properties_to_keep_track = []
            for k in range(len(gv.properties_to_keep_track)):
                if gv.properties_to_keep_track[k] < 0:
                    properties_to_keep_track.append(weight)
                else:
                    properties_to_keep_track.append(new_coordinates[gv.properties_to_keep_track[k]])
            walker_binning_value = 0
            if gv.less_or_greater_flag == 0:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] < gv.threshold_values[m]:
                        walker_binning_value += 1
            else:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        walker_binning_value += 1
            if walker_binning_value > 0:
                leftover_walker_indices.append(i)

    # if threshold values change throughout the simulation, replace all walkers with the "reference walker,"
    # i.e., the walker that had the lowest or highest values.
    if gv.static_threshold_flag == 0:
        new_coordinates = temp_walker_list[ref_walker_index].current_coordinates
        current_state = temp_walker_list[ref_walker_index].state
        center_key_num = copy.deepcopy(new_coordinates)
        balls_array[gv.current_num_balls] = np.asarray(center_key_num)
        center_key_num.append(gv.current_num_balls)
        center_key_num.append(gv.num_walkers)
        balls[gv.current_num_balls] = np.asarray(center_key_num)
        ball_to_walkers[tuple(new_coordinates)] = []
        center_key_state = copy.deepcopy(new_coordinates)
        center_key_state.append(gv.current_num_balls)
        center_key_state.append(current_state)
        ref_walker_directory = gv.main_directory + '/CAS/walker' + str(ref_walker_index)
        os.chdir(ref_walker_directory)
        f = open('ball_trajectory.txt', 'a')
        f.write(' '.join(map(lambda coordinate: str(coordinate), center_key_state)))
        f.write('\n')
        f.close()
        f = open('trajectory.txt', 'a')
        f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
        f.write('\n')
        f.close()

        for i in walker_indices:
            temp_walker_list[i].current_coordinates = new_coordinates
            temp_walker_list[i].current_ball_center = new_coordinates
            temp_walker_list[i].current_ball_key = gv.current_num_balls
            ball_to_walkers[tuple(new_coordinates)].append(i)
            if i != ref_walker_index:
                os.chdir(gv.main_directory + '/CAS')
                os.system('rm -rf walker' + str(i))
                old_directory = gv.main_directory + '/CAS/walker' + str(ref_walker_index)
                new_directory = gv.main_directory + '/CAS/walker' + str(i)
                shutil.copytree(old_directory, new_directory)
        gv.current_num_balls += 1

    # if threshold values do not change throughout the simulation,
    # bin the walkers that do not meet the threshold requirement to one designated "leftover" macrostate
    # and bin the rest of walkers normally as done in the function binning, i.e., create macrostates for the walkers.
    else:
        walker_indices_list = walker_indices.tolist()
        # if there are any walkers that did not meet the threshold requirement,
        if len(leftover_walker_indices) > 0:
            # choose the first walker in the leftover list to be the center of the "leftover" macrostate
            ref_walker_index = leftover_walker_indices[0]
            current_ball_center = temp_walker_list[ref_walker_index].current_coordinates
            center_key_num = copy.deepcopy(current_ball_center)
            balls_array[gv.current_num_balls] = np.asarray(center_key_num)
            center_key_num.append(gv.current_num_balls)
            center_key_num.append(0)
            balls[gv.current_num_balls] = np.asarray(center_key_num)
            ball_to_walkers[tuple(current_ball_center)] = []

            # and bin the rest of the leftover walkers to the leftover macrostate.
            for i in leftover_walker_indices:
                walker_indices_list.remove(i)
                current_state = temp_walker_list[i].state
                new_coordinates = temp_walker_list[i].current_coordinates
                temp_walker_list[i].current_ball_center = current_ball_center
                temp_walker_list[i].current_ball_key = gv.current_num_balls
                center_key_state = copy.deepcopy(current_ball_center)
                center_key_state.append(gv.current_num_balls)
                center_key_state.append(current_state)
                walker_directory = gv.main_directory + '/CAS/walker' + str(i)
                os.chdir(walker_directory)
                f = open('ball_trajectory.txt', 'a')
                f.write(' '.join(map(lambda coordinate: str(coordinate), center_key_state)))
                f.write('\n')
                f.close()
                f = open('trajectory.txt', 'a')
                f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
                f.write('\n')
                f.close()
                balls[gv.current_num_balls][gv.num_cvs+1] += 1
                ball_to_walkers[tuple(current_ball_center)].append(i)
            gv.current_num_balls += 1

        # bin the rest of the walkers normally as done in the function binning.
        start = 0  # indicates whether we are dealing with the very first walker or not for regular binning
        for i in walker_indices_list:
            new_coordinates = temp_walker_list[i].current_coordinates
            inside = 0  # indicates whether walker is inside an exisiting macrostate or not, i.e., binned to a macrostate
            # if we're dealing with the very first walker, create the very first ball for the walker.
            if start == 0:
                start += 1
                inside += 1
                current_ball_center = new_coordinates
                center_key_num = copy.deepcopy(current_ball_center)
                if gv.current_num_balls == 0:
                    balls_array[gv.current_num_balls] = np.asarray(center_key_num)
                else:
                    balls_array = np.append(balls_array, [np.asarray(center_key_num)], axis=0)
                center_key_num.append(gv.current_num_balls)
                center_key_num.append(1)
                if gv.current_num_balls == 0:
                    balls[gv.current_num_balls] = np.asarray(center_key_num)
                else:
                    balls = np.append(balls, [np.asarray(center_key_num)], axis=0)
                ball_to_walkers[tuple(current_ball_center)] = [i]
                temp_walker_list[i].current_ball_center = current_ball_center
                temp_walker_list[i].current_ball_key = gv.current_num_balls
                gv.current_num_balls += 1

            # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
            if inside == 0:
                current_ball_key, inside = closest_ball(new_coordinates, balls_array)

                # case 1: walker is inside some macrostate or is not but needs to be binned to the nearest macrostate
                # due to reaching the maximum number of macrostates limit.
                if inside == gv.num_cvs or (inside != gv.num_cvs and (gv.current_num_balls == gv.num_balls_limit
                                                                      or gv.balls_flag == 1)):
                    balls[current_ball_key][gv.num_cvs+1] += 1
                    current_ball_center = balls[current_ball_key][0:gv.num_cvs].tolist()
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                    temp_walker_list[i].current_ball_center = current_ball_center
                    temp_walker_list[i].current_ball_key = current_ball_key

                # case 2: walker is not inside any macrostate and the maximum number of macrostates limit
                # has not been reached, so create a new macrostate centered around the walker.
                else:
                    current_ball_center = [coordinate for coordinate in new_coordinates]
                    center_key_num = copy.deepcopy(current_ball_center)
                    balls_array = np.append(balls_array, [np.asarray(center_key_num)], axis=0)
                    center_key_num.append(gv.current_num_balls)
                    center_key_num.append(1)
                    balls = np.append(balls, [np.asarray(center_key_num)], axis=0)
                    ball_to_walkers[tuple(current_ball_center)] = [i]
                    temp_walker_list[i].current_ball_center = current_ball_center
                    temp_walker_list[i].current_ball_key = gv.current_num_balls
                    gv.current_num_balls += 1

        # loop through all of the walkers once more to assign them to their true nearest macrostates
        for i in walker_indices_list:
            current_coordinates = temp_walker_list[i].current_coordinates
            new_ball_key, inside = closest_ball(current_coordinates, balls_array)
            if inside == gv.num_cvs or (inside != gv.num_cvs and (gv.current_num_balls == gv.num_balls_limit or
                                                                          gv.balls_flag == 1)):
                new_ball_center = balls[new_ball_key][0:gv.num_cvs].tolist()
                old_ball_key = temp_walker_list[i].current_ball_key
                old_ball_center = temp_walker_list[i].current_ball_center
                balls[old_ball_key][gv.num_cvs+1] -= 1
                balls[new_ball_key][gv.num_cvs+1] += 1
                ball_to_walkers[tuple(old_ball_center)].remove(i)
                ball_to_walkers[tuple(new_ball_center)].append(i)
                temp_walker_list[i].current_ball_key = new_ball_key
                temp_walker_list[i].current_ball_center = new_ball_center

        # record the new coordinates' values on the trajectory file and the new macrostate on the ball trajectory file.
        for i in walker_indices_list:
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            new_coordinates = temp_walker_list[i].current_coordinates
            current_ball_center = temp_walker_list[i].current_ball_center
            current_ball_key = temp_walker_list[i].current_ball_key
            current_state = temp_walker_list[i].state
            center_key_state = copy.deepcopy(current_ball_center)
            center_key_state.append(current_ball_key)
            center_key_state.append(current_state)
            f = open('ball_trajectory.txt', 'a')
            f.write(' '.join(map(lambda coordinate: str(coordinate), center_key_state)))
            f.write('\n')
            f.close()
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()

        # delete empty macrostates
        delete_list = []
        for i in range(balls.shape[0]):
            if balls[i][gv.num_cvs+1] == 0:
                delete_list.append(i)
        balls = np.delete(balls, delete_list, 0)
        balls_array = np.delete(balls_array, delete_list, 0)

    # output the total flux for this particular step to a text file, if needed.
    if gv.flux_flag == 1:
        os.chdir(gv.main_directory + '/CAS')
        np.savetxt('flux_' + str(step_num+1) + '.txt', flux, fmt=' %1.5e')
    # update threshold values if they are better
    if gv.static_threshold_flag == 0:
        threshold_replace_value = 0
        if gv.less_or_greater_flag == 0:
            for m in range(len(gv.properties_to_keep_track)):
                if new_threshold_values[m] > gv.threshold_values[m]:
                    threshold_replace_value += 1
                else:
                    threshold_replace_value -= 1
        else:
            for m in range(len(gv.properties_to_keep_track)):
                if new_threshold_values[m] < gv.threshold_values[m]:
                    threshold_replace_value += 1
                else:
                    threshold_replace_value -= 1
        if threshold_replace_value > 0:
            gv.threshold_values = new_threshold_values
    return balls, balls_array


def delta2(c1, c2):
    min_dist = np.inf
    for i in xrange(0, len(c1)):
        for j in xrange(0, len(c2)):
            p1 = c1[i, :]
            p2 = c2[j, :]
            dist = np.sqrt(np.sum(np.square(p2 - p1)))
            if dist < min_dist:
                min_dist = dist
    return min_dist


def delta1(c):
    max_dist = 0
    for i in xrange(0, len(c)):
        for j in xrange(0, len(c)):
            if i == j:
                continue
            p1 = c[i, :]
            p2 = c[j, :]
            dist = np.sqrt(np.sum(np.square(p2 - p1)))
            if dist > max_dist:
                max_dist = dist
    return max_dist


def minDelta2(ball_coords):
    column = ball_coords.shape[1]-1
    num_clusters = int(np.max(ball_coords[:, column])+1)
    min_delta2 = np.inf
    for i in xrange(0, num_clusters):
        for j in xrange(0, num_clusters):
            if i == j:
                continue
            i = float(i)
            j = float(j)
            c1 = ball_coords[ball_coords[:, column] == i, :-1]
            c2 = ball_coords[ball_coords[:, column] == j, :-1]
            d2 = delta2(c1, c2)
            if d2 < min_delta2:
                min_delta2 = d2
    return min_delta2


def maxDelta1(ball_coords):
    column = ball_coords.shape[1]-1
    num_clusters = int(np.max(ball_coords[:, column])+1)
    max_delta1 = 0
    for i in xrange(0, num_clusters):
        i = float(i)
        c1 = ball_coords[ball_coords[:, column] == i, :-1]
        d1 = delta1(c1)
        if d1 > max_delta1:
            max_delta1 = d1
    return max_delta1


def dunn(ball_coords):
    num = minDelta2(ball_coords)
    den = maxDelta1(ball_coords)
    if den == 0:
        return -1
    else:
        return num/den


def create_outlier_labels(outlier_labels, new_outlier_label, matrix):
    clf = EllipticEnvelope(contamination=0.05)
    try:
        clf.fit(matrix)
        inliers = clf.predict(matrix) == 1
        i = 0
        assert len(matrix) == len(outlier_labels[outlier_labels == -1])
        for label in clf.predict(matrix):
            while outlier_labels[i] != -1:
                i += 1
            if label == -1:
                outlier_labels[i] = new_outlier_label
            i += 1
        return outlier_labels, inliers
    except ValueError:  # singular cov matrix
        return outlier_labels, [True] * len(matrix)


def merge_with_outliers(outlier_labels, labels):
    #assert len(labels) == len(outlier_labels[outlier_labels == -1]), '%d, %d, %s, %s' % (len(labels), len(outlier_labels[outlier_labels == -1]), str(labels), str(outlier_labels))
    assert len(labels) == len(outlier_labels), '%d, %d, %s, %s' % (len(labels), len(outlier_labels), str(labels), str(outlier_labels))
    rv = []
    i = 0
    #j = 0
    while True:
        while i < len(outlier_labels) and outlier_labels[i] != -1:
            rv.append(outlier_labels[i])
            i += 1
        while i < len(outlier_labels) and i < len(labels) and outlier_labels[i] == -1:
        #while i < len(outlier_labels) and j < len(labels) and outlier_labels[i] == -1:
            rv.append(labels[i])  #rv.append(labels[j])
            i += 1
            #j += 1
        if i == len(outlier_labels):
            break
    return np.array(rv)


def calculate_trans_mat(step_num, temp_walker_list, balls, balls_array):
    if (gv.enhanced_sampling_flag == 2 and step_num == gv.sc_start) or \
            (gv.enhanced_sampling_flag == 3 and step_num == gv.initial_step_num_for_eq):
        gv.trans_mat = np.zeros((balls.shape[0], balls.shape[0]))
    if gv.enhanced_sampling_flag == 2 and step_num == gv.sc_start+gv.num_steps_for_sc:
        gv.sc_performed = 1  # indicate spectral clustering can be started after this step

    if (gv.enhanced_sampling_flag == 2 and (step_num == gv.sc_start or gv.num_steps_for_sc == 0)) or \
            (gv.enhanced_sampling_flag == 3 and (step_num == gv.initial_step_num_for_eq or gv.num_steps_for_eq == 0)):
        for i in range(gv.total_num_walkers):
            previous_coordinates = temp_walker_list[i].previous_coordinates
            previous_ball_key, inside = closest_ball(previous_coordinates, balls_array)
            gv.trans_mat[previous_ball_key][temp_walker_list[i].current_ball_key] += temp_walker_list[i].weight
    else:
        for i in range(gv.total_num_walkers):
            gv.trans_mat[temp_walker_list[i].previous_ball_key][temp_walker_list[i].current_ball_key] \
                += temp_walker_list[i].weight


def reweighting(step_num, balls):
    # transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
    # canonical ensemble. equation is from Prinz, et al JCP (2011).
    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = (gv.trans_mat[i][j]+gv.trans_mat[j][i])/\
                                          (2.0*(step_num-gv.initial_step_num_for_eq-gv.num_steps_for_eq+1))

    zero_row_indices = [i for i in range(new_transition_matrix.shape[0]) if abs(np.sum(new_transition_matrix[i, :])) < 1.0e-20]
    for i in reversed(zero_row_indices):
        new_transition_matrix = np.delete(new_transition_matrix, i, 0)
        new_transition_matrix = np.delete(new_transition_matrix, i, 1)
    row_sum = np.sum(new_transition_matrix, axis=1)
    for i in range(new_transition_matrix.shape[0]):
        new_transition_matrix[i, :] /= row_sum[i]
    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('transition_matrix_' + str(step_num+1) + '.txt', new_transition_matrix, fmt=' %1.10e')

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evalues = evalues[idx]
    final_evalues = np.real(evalues)
    evectors = evectors[:, idx]
    final_evectors = np.real(evectors)
    np.savetxt('evalues_' + str(step_num+1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num+1) + '.txt', final_evectors, fmt=' %1.10e')
    eq_weights = np.zeros((balls.shape[0],))
    eq_weights_index = 0
    x_axis = np.zeros((new_transition_matrix.shape[0],))
    for i in range(eq_weights.shape[0]):
        if i in zero_row_indices:
            eq_weights[i] = 0.0
        else:
            eq_weights[i] = abs(final_evectors[eq_weights_index, 0])
            x_axis[eq_weights_index] = i
            eq_weights_index += 1
    eq_weights /= np.sum(eq_weights)  # normalize
    np.savetxt('x_axis_' + str(step_num+1) + '.txt', x_axis,fmt=' %1.5f')
    return eq_weights


def spectral_clustering(step_num, balls):
    ball_clusters_list = {}
    # transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
    # canonical ensemble. equation is from Prinz, et al JCP (2011).
    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = (gv.trans_mat[i][j]+gv.trans_mat[j][i])/\
                                          (2.0*(step_num-gv.sc_start-gv.num_steps_for_sc+1))

    row_sum = np.sum(new_transition_matrix, axis=1)
    for i in range(new_transition_matrix.shape[0]):
        if row_sum[i] > 0.0:
            new_transition_matrix[i, :] /= row_sum[i]
    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('transition_matrix_' + str(step_num+1) + '.txt', new_transition_matrix, fmt=' %1.10e')

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evalues = evalues[idx]
    final_evalues = np.real(evalues)
    evectors = evectors[:, idx]
    final_evectors = np.real(evectors)
    np.savetxt('evalues_' + str(step_num+1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num+1) + '.txt', final_evectors, fmt=' %1.10e')

    # second, normalize the second evector by the first evector values -> good approximation to committor functions.
    num_clusters = gv.num_clusters
    normalized_second_evector = np.zeros((final_evectors.shape[0], 1))
    for i in range(final_evectors.shape[0]):
        if abs(final_evectors[i, 0]) > 0.0:
            normalized_second_evector[i] = final_evectors[i, 1] / abs(final_evectors[i, 0])
        else:
            normalized_second_evector[i] = 0.0

    if np.min(normalized_second_evector) != 0.0:
        normalized_second_evector /= np.min(normalized_second_evector)  # normalize

    # third, use the normalized second evector to cluster macrostates using k-means.
    clustering_matrix = normalized_second_evector  #np.hstack((balls, normalized_second_evector))
    if abs(evalues[0]-evalues[1]) > 1.0e-14:  # avoid having eigenvalue 1 with multiplicity more than 1
        cont = True
    else:
        cont = False
        gv.sc_performed = 0
    #outlier_labels = np.ones(len(matrix))*-1
    #outliers_exist = 0
    while cont:
        try:
            centroids, labels = kmeans2(clustering_matrix, num_clusters, minit='points', iter=200, missing='raise')
            #labels = merge_with_outliers(outlier_labels, labels)
            break
        except ClusterError:
            num_clusters -= 1
        if num_clusters <= 1:
            break
    # if the number of clusters is less than or equal to 1, spectral clustering is canceled entirely.
    if num_clusters <= 1:
        gv.sc_performed = 0

        """
        # otherwise, silhouette scores are calculated and macrostates are labeled as outliers or not.
        else:
            unique = np.unique(labels)
            if len(unique) > 1:
                try:
                    silhouette_avg = silhouette_score(matrix, labels)
                    sample_silhouette_values = silhouette_samples(matrix, labels)
                except ValueError:
                    silhouette_avg = -1
                    sample_silhouette_values = [-2] * num_clusters
            else:
                silhouette_avg = 0
                sample_silhouette_values = [-1] * num_clusters

            cont = False
            if silhouette_avg > 0.8 and num_clusters >= 2:
                outliers_exist = 1
                outlier_labels, inliers = create_outlier_labels(outlier_labels, num_clusters, clustering_matrix)
                num_clusters += 1
                labels = merge_with_outliers(outlier_labels, labels)
                '''
                if len(clustering_matrix[inliers]) == len(clustering_matrix):
                    # couldn't remove any outliers; singular cov matrix (?)
                    cont = False
                    with open('outlier_removal_' + str(step_num + 1) + '.txt', 'a') as outlier_f:
                        print >>outlier_f, "Couldn't remove any outliers; just continuing"
                else:
                    cont = True
                    num_clusters -= 1
                    clustering_matrix = clustering_matrix[inliers]
                    with open('outlier_removal_' + str(step_num + 1) + '.txt', 'a') as outlier_f:
                        print >>outlier_f, 'Removing %d outliers from data as cluster %d' % (len(inliers[inliers == False]), num_clusters - 1)
                '''
            # dunn index is calculated for the entire clustering result.
            if not cont:
                with open('dunn_index_' + str(step_num + 1) + '.txt', 'w') as dunn_index_f:
                    labeled_matrix = np.zeros((matrix.shape[0], matrix.shape[1] + 1))
                    labeled_matrix[:, 0:matrix.shape[1]] = matrix
                    labeled_matrix[:, matrix.shape[1]] = labels
                    print >>dunn_index_f, dunn(labeled_matrix)
                    print >>dunn_index_f, "The average silhouette_score is: %f" % silhouette_avg
                    for i in xrange(int(max(labels))+1):
                        print >>dunn_index_f, "The average silhouette score for cluster %d is: %f" % (i, np.mean(sample_silhouette_values[labels == i]))
        """

    # finally, if clustering using k-means was successful, the results are output into text files
    # and python objects for subsequent resampling.
    num_balls = 0
    gv.balls_flag = p.balls_flag  # reset balls flag to original option
    gv.sc_start = -1
    if gv.sc_performed == 1:
        f = open('ball_clustering_'+str(step_num+1)+'.txt', 'w')
        """
        if outliers_exist == 1:
            # if outliers exist, first loop through the big main clusters
            for i in range(num_clusters):
                first = 0  # used for picking out the reference macrostate that will represent the center of the cluster
                for j in range(balls.shape[0]):
                    if labels[j] == i and first == 0:
                        first += 1
                        ref_ball_center = balls[j, 0:gv.num_cvs].tolist()
                        ball_cluster = copy.deepcopy(ref_ball_center)
                        ball_cluster.append(i)
                        ball_cluster.append(abs(final_evectors[j, 0]))
                        ball_cluster.append(final_evectors[j, 1])
                        ball_cluster.append(final_evectors[j, 2])
                        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                        f.write('\n')
                        ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
                        balls[j][gv.num_cvs+2] -= 1
                        num_balls += 1
                    elif labels[j] == i and first != 0:
                        ball_center = balls[j, 0:gv.num_cvs].tolist()
                        ball_cluster = copy.deepcopy(ball_center)
                        ball_cluster.append(i)
                        ball_cluster.append(abs(final_evectors[j, 0]))
                        ball_cluster.append(final_evectors[j, 1])
                        ball_cluster.append(final_evectors[j, 2])
                        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                        f.write('\n')
                        ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
                        balls[j][gv.num_cvs+2] -= 1
                        num_balls += 1
            # then, loop through the small, individual clusters that were labeled as outliers
            for j in range(balls.shape[0]):
                if labels[j] >= num_clusters:
                    ball_center = balls[j, 0:gv.num_cvs].tolist()
                    ball_cluster = copy.deepcopy(ball_center)
                    ball_cluster.append(labels[j])
                    ball_cluster.append(abs(final_evectors[j, 0]))
                    ball_cluster.append(final_evectors[j, 1])
                    ball_cluster.append(final_evectors[j, 2])
                    f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                    f.write('\n')
                    ball_clusters_list[tuple(ball_center)] = [tuple(ball_center)]
                    balls[j][gv.num_cvs+2] -= 1
                    num_balls += 1
        """
        for i in range(num_clusters):
            first = 0  # used for picking out the reference macrostate that will represent the center of the cluster
            for j in range(balls.shape[0]):
                if labels[j] == i and first == 0:
                    first += 1
                    ref_ball_center = balls[j, 0:gv.num_cvs].tolist()
                    ball_cluster = copy.deepcopy(ref_ball_center)
                    ball_cluster.append(i)
                    ball_cluster.append(abs(final_evectors[j, 0]))
                    ball_cluster.append(final_evectors[j, 1])
                    ball_cluster.append(final_evectors[j, 2])
                    f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                    f.write('\n')
                    ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
                    balls[j][gv.num_cvs+1] -= 1
                    num_balls += 1
                elif labels[j] == i and first != 0:
                    ball_center = balls[j, 0:gv.num_cvs].tolist()
                    ball_cluster = copy.deepcopy(ball_center)
                    ball_cluster.append(i)
                    ball_cluster.append(abs(final_evectors[j, 0]))
                    ball_cluster.append(final_evectors[j, 1])
                    ball_cluster.append(final_evectors[j, 2])
                    f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                    f.write('\n')
                    ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
                    balls[j][gv.num_cvs+1] -= 1
                    num_balls += 1
        f.close()
    if num_balls != balls.shape[0]:
        gv.sc_performed = 0
    return ball_clusters_list


def resampling_with_eq(walker_list, temp_walker_list, balls, ball_to_walkers, eq_weights):
    gv.resampling_performed = 1
    num_occupied_balls = 0
    occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers*2, int)
    excess_index = gv.total_num_walkers
    vacant_walker_indices = []
    initial_total_weight = 0.0
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs+1]) > 0:
            initial_total_weight += eq_weights[current_ball]
    if initial_total_weight != 1.0:
        factor = 1.0/initial_total_weight
        eq_weights *= factor

    # loop through each macrostate and perform resampling within each macrostate
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs+1]) > 0:
            num_occupied_balls += 1
            current_ball_center = balls[current_ball][0:gv.num_cvs].tolist()
            initial_weights = [temp_walker_list[i].weight for i in ball_to_walkers[tuple(current_ball_center)]]
            initial_indices = [temp_walker_list[i].global_index for i in ball_to_walkers[tuple(current_ball_center)]]
            # reset ball_to_walkers and balls
            ball_to_walkers[tuple(current_ball_center)] = []
            balls[current_ball][gv.num_cvs+1] = 0

            num_states = 1
            states = [-1]
            num_walkers_for_each_state = [len(initial_indices)]

            # if fluxes are calculated, we need to resample separately for each state,
            # so check to see if more than one state exists in the macrostate.
            if gv.flux_flag == 1:
                num_states = 0
                states = []
                num_walkers_for_each_state = []
                states_list = range(gv.num_states)
                for state in states_list:
                    num_walkers = 0
                    for i in initial_indices:
                        walker_state = temp_walker_list[i].state
                        if walker_state == state:
                            num_walkers += 1
                    if num_walkers != 0:
                        num_states += 1
                        states.append(state)
                        num_walkers_for_each_state.append(num_walkers)

            target_num_walkers = int(np.floor(float(gv.num_walkers)/num_states))
            remainder = gv.num_walkers-target_num_walkers*num_states
            # resample separately for each state in the macrostate
            for state_num, state in enumerate(states):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last state if there are any
                if remainder != 0 and state_num == num_states-1:
                    target_num_walkers += remainder

                weights = [float]*num_walkers_for_each_state[state_num]
                indices = [int]*num_walkers_for_each_state[state_num]

                # if the macrostate only consists of one state
                if num_states == 1:
                    weights = initial_weights
                    indices = initial_indices
                # otherwise, need to pick out the walkers that are in the particular state of interest
                else:
                    i = 0
                    for j in initial_indices:
                        walker_state = temp_walker_list[j].state
                        if state == walker_state:
                            weights[i] = temp_walker_list[j].weight
                            indices[i] = temp_walker_list[j].global_index
                            i += 1

                indices_copy = [i for i in indices]
                neg_weights = np.array([-i for i in weights])  # convert from list to array
                sorted_list = list(np.argsort(neg_weights))  # sort walkers in descending order based on their weights

                total_weight = eq_weights[current_ball]/num_states
                factor = total_weight/np.sum(weights)
                weights_copy = [i*factor for i in weights]
                weights = weights_copy
                target_weight = total_weight/target_num_walkers
                x = sorted_list.pop()
                while True:
                    x_weight = weights[x]
                    current_walker = indices[x]
                    if x_weight >= target_weight or len(sorted_list) == 0:
                        r = max(1, int(np.floor(x_weight/target_weight)))
                        r = min(r, target_num_walkers-new_num_walkers)
                        new_num_walkers += r
                        for _ in itertools.repeat(x, r):
                            new_indices.append(current_walker)
                            new_weights.append(target_weight)
                        if new_num_walkers < target_num_walkers and x_weight-r*target_weight > 0.0:
                            sorted_list.append(x)
                            weights[x] = x_weight-r*target_weight
                        if len(sorted_list) > 0:
                            x = sorted_list.pop()
                        else:
                            break
                    else:
                        y = sorted_list.pop()
                        y_weight = weights[y]
                        xy_weight = x_weight+y_weight
                        p = np.random.random()
                        if p < y_weight/xy_weight:
                            x = y
                        weights[x] = xy_weight

                for x in indices_copy:
                    if x not in new_indices:
                        vacant_walker_indices.append(x)
                        # remove walker y directory
                        os.chdir(gv.main_directory + '/CAS')
                        os.system('rm -rf walker' + str(x))

                # assign the resampled walkers to particular indices
                for index_num, global_index in enumerate(new_indices):
                    # if the global index is not used up, use it
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[index_num]
                        ball_to_walkers[tuple(current_ball_center)].append(global_index)
                        directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        os.chdir(directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()
                    # otherwise, use one of the vacant walker indices or the next smallest index available
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(temp_walker_list[global_index])
                        walker_list[new_index].weight = new_weights[index_num]
                        ball_to_walkers[tuple(current_ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        new_directory = gv.main_directory + '/CAS/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        os.system('sed -i \'$ d\' weight_trajectory.txt')
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()
                    balls[current_ball][gv.num_cvs+1] += 1

    total_num_walkers = num_occupied_balls*gv.num_walkers
    if excess_index-total_num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    # finally, re-index the walkers so that the walkers have indices in order from 0 to total_num_walkers-1
    if total_num_walkers >= gv.total_num_walkers:
        for i in range(total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(total_num_walkers, gv.total_num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= total_num_walkers:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/CAS')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))

    gv.total_num_walkers = total_num_walkers
    gv.num_occupied_balls = num_occupied_balls
    return balls


def resampling_for_sc(walker_list, temp_walker_list, balls, ball_to_walkers, ball_clusters_list):
    gv.resampling_performed = 1
    num_occupied_clusters = 0
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    occupied_indices = np.zeros(gv.num_balls_for_sc*gv.num_walkers*100, int)
    excess_index = gv.total_num_walkers
    vacant_walker_indices = []
    # loop through each cluster and perform resampling within each cluster
    for current_cluster in ball_clusters_list:
        if len(ball_clusters_list[current_cluster]) > 0:
            num_occupied_clusters += 1
            initial_target_num_walkers = gv.num_walkers_for_sc

            initial_weights = []
            initial_indices = []
            for ball_center in ball_clusters_list[current_cluster]:
                if len(ball_to_walkers[ball_center]) > 0:
                    for walker_index in ball_to_walkers[ball_center]:
                        initial_weights.append(temp_walker_list[walker_index].weight)
                        initial_indices.append(temp_walker_list[walker_index].global_index)
                    # reset ball_to_walkers and balls
                    ball_to_walkers[ball_center] = []
                    ball_key = temp_walker_list[walker_index].current_ball_key
                    balls[ball_key][gv.num_cvs+1] = 0

            initial_weights_array = np.array(initial_weights)  # convert from list to array
            walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights
            temp_initial_indices = initial_indices  # sorted indices based on descending order of weights
            initial_indices = [temp_initial_indices[i] for i in walker_indices]

            num_states = 1
            states = [0]
            num_walkers_for_each_state = [len(initial_indices)]

            # if fluxes are calculated, we need to resample separately for each state,
            # so check to see if more than one state exists in the macrostate/cluster.
            if gv.flux_flag == 1:
                num_states = 0
                states = []
                num_walkers_for_each_state = []
                states_list = range(gv.num_states)
                states_list.append(-1)
                for i in states_list:
                    num_walkers = 0
                    for j in initial_indices:
                        state = temp_walker_list[j].state
                        if state == i:
                            num_walkers += 1
                    if num_walkers != 0:
                        num_states += 1
                        states.append(i)
                        num_walkers_for_each_state.append(num_walkers)

            target_num_walkers = int(np.floor(float(initial_target_num_walkers)/num_states))
            remainder = initial_target_num_walkers-target_num_walkers*num_states
            # resample separately for each state in the macrostate/cluster
            for state_num, state in enumerate(states):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last state if there are any
                if remainder != 0 and state_num == num_states-1:
                    target_num_walkers += remainder

                weights_bin = [float]*num_walkers_for_each_state[state_num]
                indices_bin = [int]*num_walkers_for_each_state[state_num]

                # if the macrostate only consists of one state
                if num_states == 1:
                    weights_bin = initial_weights
                    indices_bin = initial_indices
                # otherwise, need to pick out the walkers that are in the particular state of interest
                else:
                    k = 0
                    for j in initial_indices:
                        walker_state = temp_walker_list[j].state
                        if state == walker_state:
                            weights_bin[k] = temp_walker_list[j].weight
                            indices_bin[k] = temp_walker_list[j].global_index
                            k += 1

                total_weight = np.sum(weights_bin)
                target_weight = total_weight/target_num_walkers
                x = indices_bin.pop()
                while True:
                    x_weight = weights[x]
                    if x_weight >= target_weight or len(indices_bin) == 0:
                        r = max(1, int(np.floor(x_weight/target_weight)))
                        r = min(r, target_num_walkers-new_num_walkers)
                        new_num_walkers += r
                        for item in np.repeat(x, r):
                            new_indices.append(item)
                            new_weights.append(target_weight)
                        if new_num_walkers < target_num_walkers and x_weight-r*target_weight > 0.0:
                            indices_bin.append(x)
                            weights[x] = x_weight-r*target_weight
                        if len(indices_bin) > 0:
                            x = indices_bin.pop()
                        else:
                            break
                    else:
                        y = indices_bin.pop()
                        y_weight = weights[y]
                        xy_weight = x_weight + y_weight
                        p = np.random.random()
                        # swap x and y
                        if p < y_weight / xy_weight:
                            temp = x
                            x = y
                            y = temp
                        weights[x] = xy_weight
                        if y not in new_indices:
                            vacant_walker_indices.append(y)
                            # remove walker y directory
                            os.chdir(gv.main_directory + '/CAS')
                            os.system('rm -rf walker' + str(y))

                # assign the resampled walkers to particular indices
                for index_num, global_index in enumerate(new_indices):
                    # if the global index is not used up, use it
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[index_num]
                        ball_key = walker_list[global_index].current_ball_key
                        if balls[ball_key][gv.num_cvs+1] == 0:
                            num_occupied_balls += 1
                        balls[ball_key][gv.num_cvs+1] += 1
                        ball_center = walker_list[global_index].current_ball_center
                        ball_to_walkers[tuple(ball_center)].append(global_index)
                        directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        os.chdir(directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()
                    # otherwise, use one of the vacant walker indices or the next smallest index available
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(temp_walker_list[global_index])
                        walker_list[new_index].weight = new_weights[index_num]
                        ball_key = walker_list[new_index].current_ball_key
                        if balls[ball_key][gv.num_cvs+1] == 0:
                            num_occupied_balls += 1
                        balls[ball_key][gv.num_cvs+1] += 1
                        ball_center = walker_list[new_index].current_ball_center
                        ball_to_walkers[tuple(ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        new_directory = gv.main_directory + '/CAS/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        os.system('sed -i \'$ d\' weight_trajectory.txt')
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()

    total_num_walkers = num_occupied_clusters*gv.num_walkers_for_sc
    if excess_index - total_num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    # finally, re-index the walkers so that the walkers have indices in order from 0 to total_num_walkers-1
    if total_num_walkers >= gv.total_num_walkers:
        for i in range(total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(total_num_walkers, gv.total_num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= total_num_walkers:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/CAS')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))

    gv.total_num_walkers = total_num_walkers
    gv.num_occupied_balls = num_occupied_balls
    gv.num_occupied_clusters = num_occupied_clusters
    return balls


def resampling(step_num, walker_list, temp_walker_list, balls, ball_to_walkers):
    gv.resampling_performed = 1
    if gv.enhanced_sampling_flag == 2 and gv.sc_performed == 1:
        gv.sc_performed = 0
    elif gv.enhanced_sampling_flag == 2 and gv.sc_performed == 0 and gv.sc_start != -1 and \
                    step_num == gv.sc_start + gv.num_steps_for_sc:  # in case spectral clustering failed
        gv.sc_performed = 1
        gv.sc_start = -1
    num_occupied_balls = 0
    if gv.enhanced_sampling_flag == 2:
        occupied_indices = np.zeros(gv.num_balls_for_sc*gv.num_walkers*100, int)
    else:
        occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers*2, int)
    excess_index = gv.total_num_walkers
    vacant_walker_indices = []
    # loop through each macrostate and perform resampling within each macrostate
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs+1]) > 0:
            num_occupied_balls += 1
            current_ball_center = balls[current_ball][0:gv.num_cvs].tolist()
            initial_weights = [temp_walker_list[i].weight for i in ball_to_walkers[tuple(current_ball_center)]]
            initial_indices = [temp_walker_list[i].global_index for i in ball_to_walkers[tuple(current_ball_center)]]
            # reset ball_to_walkers and balls
            ball_to_walkers[tuple(current_ball_center)] = []
            balls[current_ball][gv.num_cvs+1] = 0

            num_states = 1
            states = [-1]
            num_walkers_for_each_state = [len(initial_indices)]

            # if fluxes are calculated, we need to resample separately for each state,
            # so check to see if more than one state exists in the macrostate.
            if gv.flux_flag == 1:
                num_states = 0
                states = []
                num_walkers_for_each_state = []
                states_list = range(gv.num_states)
                for state in states_list:
                    num_walkers = 0
                    for i in initial_indices:
                        walker_state = temp_walker_list[i].state
                        if walker_state == state:
                            num_walkers += 1
                    if num_walkers != 0:
                        num_states += 1
                        states.append(state)
                        num_walkers_for_each_state.append(num_walkers)

            target_num_walkers = int(np.floor(float(gv.num_walkers)/num_states))
            remainder = gv.num_walkers-target_num_walkers*num_states
            # resample separately for each state in the macrostate
            for state_num, state in enumerate(states):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last state if there are any
                if remainder != 0 and state_num == num_states-1:
                    target_num_walkers += remainder

                weights = [float]*num_walkers_for_each_state[state_num]
                indices = [int]*num_walkers_for_each_state[state_num]

                # if the macrostate only consists of one state
                if num_states == 1:
                    weights = initial_weights
                    indices = initial_indices
                # otherwise, need to pick out the walkers that are in the particular state of interest
                else:
                    i = 0
                    for j in initial_indices:
                        walker_state = temp_walker_list[j].state
                        if state == walker_state:
                            weights[i] = temp_walker_list[j].weight
                            indices[i] = temp_walker_list[j].global_index
                            i += 1

                indices_copy = [i for i in indices]
                neg_weights = np.array([-i for i in weights])  # convert from list to array
                sorted_list = list(np.argsort(neg_weights))  # sort walkers in descending order based on their weights

                total_weight = np.sum(weights)
                target_weight = total_weight/target_num_walkers
                x = sorted_list.pop()
                while True:
                    x_weight = weights[x]
                    current_walker = indices[x]
                    if x_weight >= target_weight or len(sorted_list) == 0:
                        r = max(1, int(np.floor(x_weight/target_weight)))
                        r = min(r, target_num_walkers-new_num_walkers)
                        new_num_walkers += r
                        for _ in itertools.repeat(x, r):
                            new_indices.append(current_walker)
                            new_weights.append(target_weight)
                        if new_num_walkers < target_num_walkers and x_weight-r*target_weight > 0.0:
                            sorted_list.append(x)
                            weights[x] = x_weight-r*target_weight
                        if len(sorted_list) > 0:
                            x = sorted_list.pop()
                        else:
                            break
                    else:
                        y = sorted_list.pop()
                        y_weight = weights[y]
                        xy_weight = x_weight+y_weight
                        p = np.random.random()
                        if p < y_weight/xy_weight:
                            x = y
                        weights[x] = xy_weight

                for x in indices_copy:
                    if x not in new_indices:
                        vacant_walker_indices.append(x)
                        # remove walker y directory
                        os.chdir(gv.main_directory + '/CAS')
                        os.system('rm -rf walker' + str(x))

                # assign the resampled walkers to particular indices
                for index_num, global_index in enumerate(new_indices):
                    # if the global index is not used up, use it
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[index_num]
                        ball_to_walkers[tuple(current_ball_center)].append(global_index)
                        directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        os.chdir(directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()
                    # otherwise, use one of the vacant walker indices or the next smallest index available
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(temp_walker_list[global_index])
                        walker_list[new_index].weight = new_weights[index_num]
                        ball_to_walkers[tuple(current_ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        new_directory = gv.main_directory + '/CAS/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        os.system('sed -i \'$ d\' weight_trajectory.txt')
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[index_num] + '\n')
                        f.close()
                    balls[current_ball][gv.num_cvs+1] += 1

    total_num_walkers = num_occupied_balls*gv.num_walkers
    if excess_index-total_num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    # finally, re-index the walkers so that the walkers have indices in order from 0 to total_num_walkers-1
    if total_num_walkers >= gv.total_num_walkers:
        for i in range(total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/CAS')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(total_num_walkers, gv.total_num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= total_num_walkers:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/CAS')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))

    gv.total_num_walkers = total_num_walkers
    gv.num_occupied_balls = num_occupied_balls
    return balls


def print_status(step_num, walker_list, balls, ball_to_walkers):
    gv.resampling_performed = 0
    os.chdir(gv.main_directory + '/CAS')
    total_weight = 0.0
    f = open('total_weight_on_each_ball_' + str(step_num+1) + '.txt', 'w')
    for current_ball in range(balls.shape[0]):
        ball_center = balls[current_ball][0:gv.num_cvs].tolist()
        weights = [walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
        total_weight += np.sum(weights)
        ball_center_weights = copy.deepcopy(ball_center)
        ball_center_weights.append(np.sum(weights))
        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_center_weights)))
        f.write('\n')

        # reset walkers and number of walkers that belong in each ball
        balls[current_ball][gv.num_cvs+1] = 0
        ball_to_walkers[tuple(ball_center)] = []
    f.close()

    # verify that total weight of all balls is 1.0
    f = open('total_weight.txt', 'a')
    if gv.enhanced_sampling_flag == 2:
        f.write(str(step_num+1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + ' '
                + str(gv.num_occupied_clusters) + ' ' + str(gv.total_num_walkers) + '\n')
        gv.num_occupied_clusters = 0
    else:
        f.write(str(step_num+1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + ' '
                + str(gv.total_num_walkers) + '\n')
    return balls
