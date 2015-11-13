import numpy as np
import os
import shutil
import copy
from scipy import special
from scipy.cluster.vq import kmeans2, ClusterError
import walker
import we_global_variables as gv
import we_check_state_function
import we_parameters as p


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += (values[i] - center[i]) ** 2
    if abs(distance) < 1.0e-10:
        distance = 0.0
    return np.sqrt(distance)


def set_parameters():
    gv.main_directory = p.main_directory
    gv.initial_configuration_directory = p.initial_configuration_directory
    gv.simulation_flag = p.simulation_flag
    gv.balls_flag = p.balls_flag
    gv.sorting_flag = p.sorting_flag
    gv.rate_flag = p.rate_flag
    gv.num_states = p.num_states
    gv.enhanced_sampling_flag = p.enhanced_sampling_flag
    gv.num_balls_limit = p.num_balls_limit
    gv.radius = p.radius
    gv.num_walkers = p.num_walkers
    gv.num_cvs = p.num_cvs
    gv.lower_bound = p.lower_bound
    gv.upper_bound = p.upper_bound
    gv.initial_step_num = p.initial_step_num
    gv.max_num_steps = p.max_num_steps
    gv.num_occupied_balls = p.num_occupied_balls
    gv.first_walker = p.first_walker
    gv.last_walker = p.last_walker
    if gv.enhanced_sampling_flag == 2:
        gv.less_or_greater_flag = p.less_or_greater_flag
        gv.static_threshold_flag = p.static_threshold_flag
        gv.threshold_values = p.threshold_values
        gv.properties_to_keep_track = p.properties_to_keep_track
    elif gv.enhanced_sampling_flag == 3:
        gv.num_balls_for_sc = p.num_balls_for_sc
        gv.num_clusters = p.num_clusters
        gv.num_walkers_for_sc = p.num_walkers_for_sc
        gv.timestep = p.timestep

    ball_volume = (np.pi ** (gv.num_cvs / 2) * gv.radius ** gv.num_cvs) / special.gamma((gv.num_cvs / 2) + 1)
    gv.max_num_balls = 0
    if ball_volume != 0.0:
        gv.max_num_balls = int(np.floor((gv.upper_bound - gv.lower_bound) ** gv.num_cvs / ball_volume))
    if gv.max_num_balls > gv.num_balls_limit or gv.max_num_balls < gv.num_balls_limit * 1e-2:
        gv.max_num_balls = gv.num_balls_limit
    print 'max # of balls (n_b) = ' + str(gv.max_num_balls)
    gv.current_num_balls = 0
    gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers


def initialize(input_initial_values_file, walker_list, temp_walker_list, ball_to_walkers, vacant_walker_indices):
    for i in range(len(walker_list)):
        walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, i, 0.0, [-1000.0] * gv.num_cvs,
                                       [-1000.0] * gv.num_cvs, 0, 0.0, 0.0, 0, 0.0, -1)

    if gv.simulation_flag == 0:  # new simulation
        initial_weight = 1.0/gv.total_num_walkers
        f = open(input_initial_values_file, 'r')
        for n in range(gv.num_occupied_balls):
            initial_values = [None] * gv.num_cvs
            for i in range(gv.num_cvs):
                initial_values[i] = float(f.readline())
            if gv.rate_flag == 1:
                initial_state = we_check_state_function.check_state_function(initial_values)
                print initial_state
            for i in range(n * gv.num_walkers, (n + 1) * gv.num_walkers):
                walker_list[i].set(initial_values, initial_weight)
                if gv.rate_flag == 1:
                    walker_list[i].state = initial_state
        f.close()

        os.system('mkdir WE')
        os.chdir(gv.main_directory + '/WE')
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            shutil.copytree(gv.initial_configuration_directory, walker_directory)

    elif gv.simulation_flag == 1:  # restarting simulation in the middle of simulation
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            trajectory = np.loadtxt('trajectory.txt')
            previous_coordinates = trajectory[-2].tolist()
            current_coordinates = trajectory[-1].tolist()
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            ball_trajectory = np.loadtxt('ball_trajectory.txt')
            previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
            current_ball_center = ball_trajectory[-1][0:gv.num_cvs].tolist()
            current_ball_radius = ball_trajectory[-1][gv.num_cvs]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = current_ball_radius
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                         current_ball_center)
            if gv.rate_flag == 1:
                walker_list[i].state = int(ball_trajectory[-1][-1])

    elif gv.simulation_flag == 2:  # restarting simulation in the middle of binning
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            trajectory = np.loadtxt('trajectory.txt')
            previous_coordinates = trajectory[-2].tolist()
            current_coordinates = trajectory[-1].tolist()
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            num_lines = sum(1 for line in open('ball_trajectory.txt'))
            # if walker is already binned to a ball, delete the binning and have binning start from scratch
            if num_lines > gv.initial_step_num:
                os.system('sed -i \'$d\' ball_trajectory.txt')
            ball_trajectory = np.loadtxt('ball_trajectory.txt')
            previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
            current_ball_center = ball_trajectory[-1][0:gv.num_cvs].tolist()
            current_ball_radius = ball_trajectory[-1][gv.num_cvs]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = current_ball_radius
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                         current_ball_center)
            if gv.rate_flag == 1:
                walker_list[i].state = int(ball_trajectory[-1][-1])

    elif gv.simulation_flag == 3:  # restarting simulation in the middle of resampling
        total_weight = 0.0
        previous_ball_to_walkers = {}
        previous_balls_weights = np.loadtxt('total_weight_of_each_ball_' + str(gv.initial_step_num) + '.txt')
        previous_balls_walker_count = np.zeros((previous_balls_weights.shape[0], previous_balls_weights.shape[1]))
        for i in range(previous_balls_weights.shape[0]):
            previous_balls_walker_count[i] = previous_balls_weights[i]
            previous_balls_walker_count[i][-1] = gv.num_walkers

        # TODO: make sure that gv.num_occupied_balls is equal to the highest walker number inside the WE folder

        for i in range(gv.num_occupied_balls + 1):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            # if all of the files exist in the walker folder, it is a complete walker
            if os.path.isfile(walker_directory + '/weight_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/ball_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/trajectory.txt') and \
                    os.path.isfile(walker_directory + '/traj.xtc') and os.path.isfile(walker_directory + '/minim.gro'):
                os.chdir(walker_directory)
                f = open('weight_trajectory.txt', 'r')
                weight = float(f.readlines()[-1].strip())
                walker_list[i].weight = weight
                total_weight += weight
                f.close()

                ball_trajectory = np.loadtxt('ball_trajectory.txt')
                previous_ball = ball_trajectory[-2].tolist()
                previous_ball_key = previous_ball[gv.num_cvs+1]
                previous_ball_center = previous_ball[0:gv.num_cvs]
                previous_balls_weights[previous_ball_key][-1] -= weight
                if previous_balls_weights[previous_ball_key][-1] < 0.0:
                    print 'ERROR: weight is ' + str(previous_balls_weights[previous_ball_key][-1]) + ' for walker ' + \
                          str(i) + ' with ball_key ' + str(previous_ball_key)
                previous_balls_walker_count[previous_ball_key][-1] -= 1
                if previous_balls_walker_count[previous_ball_key][-1] < 0:
                    print 'ERROR: walker count is ' + str(previous_balls_walker_count[previous_ball_key][-1]) + \
                          ' for walker ' + str(i) + ' with ball key ' + str(previous_ball_key)
                if tuple(previous_ball_center) in previous_ball_to_walkers:
                    previous_ball_to_walkers[tuple(previous_ball_center)].append(i)
                else:
                    previous_ball_to_walkers[tuple(previous_ball_center)] = [i]
                current_ball = ball_trajectory[-1].tolist()
                current_ball_radius = current_ball[gv.num_cvs]
                current_ball_key = current_ball[gv.num_cvs+1]
                current_ball_center = current_ball[0:gv.num_cvs]
                if tuple(current_ball_center) in ball_to_walkers:
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                else:
                    ball_to_walkers[tuple(current_ball_center)] = [i]
                    gv.current_num_balls += 1
                walker_list[i].radius = current_ball_radius
                walker_list[i].previous_ball_center = previous_ball_center
                walker_list[i].current_ball_center = current_ball_center

                trajectory = np.loadtxt('trajectory.txt')
                previous_coordinates = trajectory[-2].tolist()
                current_coordinates = trajectory[-1].tolist()
                walker_list[i].previous_coordinates = previous_coordinates
                walker_list[i].current_coordinates = current_coordinates
                if gv.rate_flag == 1:
                    current_state = int(current_ball[-1])
                else:
                    current_state = -1
                walker_list[i].state = current_state
                walker_list[i].ball_key = current_ball_key
                previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                               previous_ball_center)
                current_distance_from_center = calculate_distance_from_center(current_coordinates, current_ball_center)
                walker_list[i].previous_distance_from_center = previous_distance_from_center
                walker_list[i].current_distance_from_center = current_distance_from_center

                temp_walker_list[i] = walker.Walker(previous_coordinates, current_coordinates, i, current_ball_radius,
                                                    previous_ball_center, current_ball_center, current_ball_key,
                                                    previous_distance_from_center, current_distance_from_center, 0,
                                                    weight, current_state)

            # otherwise, it is an incomplete walker that needs missing files
            else:
                if os.path.isdir(walker_directory):
                    os.chdir(gv.main_directory + '/WE')
                    os.system('rm -rf walker' + str(i))
                vacant_walker_indices.append(i)

        # create new walkers for the remaining weights
        excess_index = gv.num_occupied_balls + 1
        for i in range(previous_balls_weights.shape[0]):
            if previous_balls_weights[i][-1] > 0.0:
                if previous_balls_walker_count[i][-1] <= 0:
                    print 'ERROR: at least one walker should exist if there is a weight of ' + \
                          str(previous_balls_weights[i][-1]) + ' for walker ' + str(i)
                else:
                    current_ball_center = previous_balls_weights[i][0:gv.num_cvs].tolist()
                    reference_walker = ball_to_walkers[tuple(current_ball_center)][0]
                    reference_walker_directory = gv.main_directory + '/WE/walker/' + str(reference_walker)
                    if len(vacant_walker_indices) > 0:
                        walker_index = vacant_walker_indices.pop(0)
                    else:
                        walker_index = excess_index
                        excess_index += 1
                    walker_directory = gv.main_directory + '/WE/walker' + str(walker_index)
                    shutil.copytree(reference_walker_directory, walker_directory)

                    weight = previous_balls_weights[i][-1]
                    previous_balls_weights[i][-1] -= weight

                    os.chdir(walker_directory)
                    f = open('weight_trajectory.txt', 'w')
                    f.write(str(weight) + '\n')
                    walker_list[walker_index].weight = weight
                    total_weight += weight
                    f.close()

                    ball_to_walkers[tuple(current_ball_center)].append(walker_index)
                    walker_list[walker_index].current_ball_center = current_ball_center

                    trajectory = np.loadtxt('trajectory.txt')
                    previous_coordinates = trajectory[-2].tolist()
                    current_coordinates = trajectory[-1].tolist()
                    walker_list[walker_index].previous_coordinates = previous_coordinates
                    walker_list[walker_index].current_coordinates = current_coordinates
                    ball_trajectory = np.loadtxt('ball_trajectory.txt')
                    previous_ball_center = ball_trajectory[-2][0:gv.num_cvs].tolist()
                    walker_list[i].previous_ball_center = previous_ball_center
                    current_state = ball_trajectory[-1][-1]
                    current_ball_key = ball_trajectory[-1][gv.num_cvs+1]
                    current_ball_radius = ball_trajectory[-1][gv.num_cvs]
                    walker_list[walker_index].state = current_state
                    walker_list[walker_index].ball_key = current_ball_key
                    walker_list[walker_index].radius = current_ball_radius
                    previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                   previous_ball_center)
                    current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                  current_ball_center)
                    walker_list[i].previous_distance_from_center = previous_distance_from_center
                    walker_list[i].current_distance_from_center = current_distance_from_center

                    temp_walker_list[walker_index] = walker.Walker(previous_coordinates, current_coordinates,
                                                                   walker_index, current_ball_radius,
                                                                   previous_ball_center, current_ball_center,
                                                                   current_ball_key, previous_distance_from_center,
                                                                   current_distance_from_center, 0, weight,
                                                                   current_state)

        # check if total weight is 1.0
        if total_weight != 1.0:
            print 'ERROR: total weight is ' + str(total_weight)


def binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers, key_to_ball):
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)
    flux = np.zeros((gv.num_states, gv.num_states))
    if gv.sorting_flag == 1:
        walker_indices = np.argsort(initial_weights_array)  # sort walkers in ascending order based on their weights
    else:
        walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    start = 0  # indicates whether we are dealing with the very first walker or not
    if gv.enhanced_sampling_flag == 2:
        ref_walker = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, 0, 0.0, [-1000.0] * gv.num_cvs,
                                   [-1000.0] * gv.num_cvs, 0, 0.0, 0.0, 0, 0.0, -1)
        ref_walker_binning_value = len(gv.properties_to_keep_track)
        ref_walker_properties_value = 0.0
        if gv.static_threshold_flag == 0:
            new_threshold_values = gv.threshold_values

    for i in walker_indices:
        # first, go to walker directory i
        walker_directory = gv.main_directory + '/WE/walker' + str(i)
        os.chdir(walker_directory)

        # then, obtain new coordinates' values
        if os.path.exists(walker_directory + '/coordinates.out'):
            new_coordinates = np.loadtxt('coordinates.out')
            new_coordinates = new_coordinates.tolist()
            rm_command = 'rm *.out'
            os.system(rm_command)

            # also, write the new coordinates' values  on the trajectory file
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()
        else:
            f = open('trajectory.txt', 'r')
            new_coordinates = f.readlines()[-1].strip().split()
            new_coordinates = [float(coordinate) for coordinate in new_coordinates]
            f.close()

        previous_coordinates = walker_list[i].current_coordinates
        previous_ball_center = walker_list[i].current_ball_center
        previous_distance_from_center = walker_list[i].current_distance_from_center
        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight

        if gv.rate_flag == 1:
            state = we_check_state_function.check_state_function(new_coordinates)
            if walker_list[i].state != -1 and state == -1:
                state = walker_list[i].state
            if walker_list[i].state != -1 and state != -1:
                flux[walker_list[i].state, state] += walker_list[i].weight
        else:
            state = -1

        if gv.enhanced_sampling_flag == 2:
            properties_to_keep_track = []
            for k in range(len(gv.properties_to_keep_track)):
                if gv.properties_to_keep_track[k] < 0:
                    properties_to_keep_track.append(weight)
                else:
                    properties_to_keep_track.append(new_coordinates[gv.properties_to_keep_track[k]])
            walker_binning_value = 0
            walker_properties_value = 0.0
            if gv.less_or_greater_flag == 0:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] < gv.threshold_values[m]:
                        walker_binning_value += 1
                        walker_properties_value += (gv.threshold_values[m]-properties_to_keep_track[m])
            elif gv.less_or_greater_flag == 1:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        walker_binning_value += 1
                        walker_properties_value += (properties_to_keep_track[m]-gv.threshold_values[m])

        inside = 0  # indicates whether we are dealing with the very first walker or not
        # if balls_flag = 0 and if we're dealing with the very first walker, create the very first ball for the walker
        if (gv.balls_flag == 0 and start == 0) or (gv.balls_flag == 1 and start == 0 and step_num == 0):
            start += 1
            inside += 1
            current_ball_center = [coordinate for coordinate in new_coordinates]
            ball_to_walkers[tuple(current_ball_center)] = [i]
            temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                previous_ball_center, current_ball_center, gv.current_num_balls,
                                                previous_distance_from_center, 0.0, initial_step_num, weight, state)
            if gv.enhanced_sampling_flag == 2:
                ref_walker = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius, previous_ball_center,
                                           current_ball_center, gv.current_num_balls,
                                           previous_distance_from_center, 0.0, initial_step_num, weight, state)
                ref_walker_binning_value = walker_binning_value
                ref_walker_properties_value = walker_properties_value
            center_r_key_num = copy.deepcopy(current_ball_center)
            center_r_key_num.append(gv.radius)
            center_r_key_num.append(gv.current_num_balls)
            center_r_key_num.append(1)
            balls[gv.current_num_balls] = np.asarray(center_r_key_num)
            key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
            gv.current_num_balls += 1

        distance = 0.0
        ball_key = 0
        # otherwise, loop through all of the balls and find the ball that has a center nearest the walker
        if inside == 0:
            for j in range(balls.shape[0]):
                current_ball_center = balls[j][0:gv.num_cvs].tolist()
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                if distance_from_center <= gv.radius or abs(distance_from_center - gv.radius) < 1.0e-10:
                    inside += 1
                if distance == 0.0:
                    distance = distance_from_center
                    ball_key = j
                else:
                    if distance_from_center < distance:
                        distance = distance_from_center
                        ball_key = j

            # walker is inside some ball
            if inside != 0:
                balls[ball_key][gv.num_cvs+2] += 1
                current_ball_center = balls[ball_key][0:gv.num_cvs].tolist()
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, ball_key,
                                                    previous_distance_from_center, distance_from_center,
                                                    initial_step_num, weight, state)
                if gv.enhanced_sampling_flag == 2 and ((gv.balls_flag == 1 and start == 0) or (walker_binning_value ==
                    ref_walker_binning_value and walker_properties_value < ref_walker_properties_value) or
                    walker_binning_value < ref_walker_binning_value):
                    ref_walker = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                               previous_ball_center, current_ball_center, ball_key,
                                               previous_distance_from_center, distance_from_center, initial_step_num,
                                               weight, state)
                    ref_walker_binning_value = walker_binning_value
                    ref_walker_properties_value = walker_properties_value
                    new_threshold_values = properties_to_keep_track
                ball_to_walkers[tuple(current_ball_center)].append(i)

            # or walker does not belong in any ball -> create a new ball
            elif gv.enhanced_sampling_flag != 2:
                current_ball_center = [coordinate for coordinate in new_coordinates]
                ball_to_walkers[tuple(current_ball_center)] = [i]
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, gv.current_num_balls,
                                                    previous_distance_from_center, 0.0, initial_step_num, weight, state)
                center_r_key_num = copy.deepcopy(current_ball_center)
                center_r_key_num.append(gv.radius)
                center_r_key_num.append(gv.current_num_balls)
                center_r_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                gv.current_num_balls += 1

            # or if enhanced_sampling_flag = 2 and ref_walker is a "better" walker in terms of its values -> replace
            # walker with ref_walker
            elif gv.enhanced_sampling_flag == 2 and ((walker_binning_value == ref_walker_binning_value and
                walker_properties_value > ref_walker_properties_value) or walker_binning_value >
                ref_walker_binning_value):
                balls[ref_walker.ball_key][gv.num_cvs+2] += 1
                temp_walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, 0, 0.0,
                                                    [-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, 0, 0.0, 0.0, 0, 0.0,
                                                    -1)
                temp_walker_list[i].copy_walker(ref_walker)
                temp_walker_list[i].global_index = i
                temp_walker_list[i].weight = weight
                current_ball_center = temp_walker_list[i].current_ball_center
                ball_to_walkers[tuple(current_ball_center)].append(i)

            # or if enhanced_sampling_flag = 2 and walker is a "better" or "equivalent" walker in terms of its values
            # -> create a new ball
            elif gv.enhanced_sampling_flag == 2 and ((gv.balls_flag == 1 and start == 0) or (walker_binning_value ==
                ref_walker_binning_value and walker_properties_value <= ref_walker_properties_value) or
                walker_binning_value < ref_walker_binning_value):
                current_ball_center = [coordinate for coordinate in new_coordinates]
                ball_to_walkers[tuple(current_ball_center)] = [i]
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, gv.current_num_balls,
                                                    previous_distance_from_center, 0.0, initial_step_num, weight, state)
                ref_walker = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius, previous_ball_center,
                                           current_ball_center, gv.current_num_balls, previous_distance_from_center,
                                           0.0, initial_step_num, weight, state)
                ref_walker_binning_value = walker_binning_value
                ref_walker_properties_value = walker_properties_value
                new_threshold_values = properties_to_keep_track
                center_r_key_num = copy.deepcopy(current_ball_center)
                center_r_key_num.append(gv.radius)
                center_r_key_num.append(gv.current_num_balls)
                center_r_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                gv.current_num_balls += 1

        # finally, write the new ball on the trajectory file
        if gv.enhanced_sampling_flag != 2:
            current_ball_center = temp_walker_list[i].current_ball_center
            ball_key = temp_walker_list[i].ball_key
            center_r_key_state = copy.deepcopy(current_ball_center)
            center_r_key_state.append(gv.radius)
            center_r_key_state.append(ball_key)
            center_r_key_state.append(state)
            f = open('ball_trajectory.txt', 'a')
            f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
            f.write('\n')
            f.close()

    # if enhanced_sampling_flag = 2, replace "inadequate" walkers with ref_walker
    if gv.enhanced_sampling_flag == 2:
        for i in walker_indices:
            new_coordinates = walker_list[i].current_coordinates
            weight = walker_list[i].weight
            properties_to_keep_track = []
            for k in range(len(gv.properties_to_keep_track)):
                if gv.properties_to_keep_track[k] < 0:
                    properties_to_keep_track.append(weight)
                else:
                    properties_to_keep_track.append(new_coordinates[gv.properties_to_keep_track[k]])
            walker_binning_value = 0
            walker_properties_value = 0.0
            if gv.less_or_greater_flag == 0:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] < gv.threshold_values[m]:
                        walker_binning_value += 1
                        walker_properties_value += (gv.threshold_values[m]-properties_to_keep_track[m])
            elif gv.less_or_greater_flag == 1:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        walker_binning_value += 1
                        walker_properties_value += (properties_to_keep_track[m]-gv.threshold_values[m])
            if (walker_binning_value > ref_walker_binning_value or (walker_binning_value == ref_walker_binning_value and
                walker_properties_value > ref_walker_properties_value)):
                previous_ball_center = temp_walker_list[i].current_ball_center
                previous_ball_key = temp_walker_list[i].ball_key
                balls[previous_ball_key][gv.num_cvs+2] -= 1
                balls[ref_walker.ball_key][gv.num_cvs+2] += 1
                temp_walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, 0, 0.0,
                                                    [-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, 0, 0.0, 0.0, 0, 0.0,
                                                    -1)
                temp_walker_list[i].copy_walker(ref_walker)
                temp_walker_list[i].global_index = i
                temp_walker_list[i].weight = weight
                ball_to_walkers[tuple(previous_ball_center)].remove(i)
                current_ball_center = temp_walker_list[i].current_ball_center
                ball_to_walkers[tuple(current_ball_center)].append(i)
            current_ball_center = temp_walker_list[i].current_ball_center
            ball_key = temp_walker_list[i].ball_key
            center_r_key_state = copy.deepcopy(current_ball_center)
            center_r_key_state.append(gv.radius)
            center_r_key_state.append(ball_key)
            center_r_key_state.append(state)
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            os.chdir(walker_directory)
            f = open('ball_trajectory.txt', 'a')
            f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
            f.write('\n')
            f.close()
            num_lines = sum(1 for line in open('trajectory.txt'))
            if num_lines > step_num:
                os.system('sed -i \'$d\' trajectory.txt')
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()

    os.chdir(gv.main_directory + '/WE')
    np.savetxt('balls_' + str(step_num + 1) + '.txt', balls, fmt=' %+1.5f')
    if gv.rate_flag == 1:
        np.savetxt('flux_' + str(step_num + 1) + '.txt', flux, fmt=' %1.5e')
    # update threshold values
    if gv.enhanced_sampling_flag == 2 and gv.static_threshold_flag == 0:
        gv.threshold_values = new_threshold_values
    return balls
def delta2(c1, c2):
  minDist = np.inf
  for i in xrange(0, len(c1)):
    for j in xrange(0, len(c2)):
      p1 = c1[i,:]
      p2 = c2[j,:]
      dist = np.sqrt(np.sum(np.square(p2 - p1)))
      if dist < minDist:
        minDist = dist
  return minDist

def delta1(c):
  maxDist = 0
  for i in xrange(0, len(c)):
    for j in xrange(0, len(c)):
      if i == j:
        continue
      p1 = c[i,:]
      p2 = c[j,:]
      dist = np.sqrt(np.sum(np.square(p2 - p1)))
      if dist > maxDist:
        maxDist = dist
  return maxDist

def minDelta2(ball_coords):
  column = ball_coords.shape[1]-1
  num_clusters = int(np.max(ball_coords[:,column])+1)
  min_delta2 = np.inf
  for i in xrange(0,num_clusters):
    for j in xrange(0,num_clusters):
      if i == j:
        continue
      i = float(i)
      j = float(j)
      c1 = ball_coords[ball_coords[:,column] == i,:-1]
      c2 = ball_coords[ball_coords[:,column] == j,:-1]
      d2 = delta2(c1, c2)
      if d2 < min_delta2:
        min_delta2 = d2
  return min_delta2

def maxDelta1(ball_coords):
  column = ball_coords.shape[1]-1
  num_clusters = int(np.max(ball_coords[:,column])+1)
  max_delta1 = 0
  for i in xrange(0,num_clusters):
    i = float(i)
    c1 = ball_coords[ball_coords[:,column] == i,:-1]
    d1 = delta1(c1)
    if d1 > max_delta1:
      max_delta1 = d1
  return max_delta1

def dunn(ball_coords):
  return minDelta2(ball_coords) / maxDelta1(ball_coords)

def spectral_clustering(step_num, temp_walker_list, balls, ball_clusters_list):
    transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(gv.total_num_walkers):
        if temp_walker_list[i].previous_distance_from_center <= gv.radius:
            previous_coordinates = temp_walker_list[i].previous_coordinates
        else:
            previous_coordinates = temp_walker_list[i].previous_ball_center
        if temp_walker_list[i].current_distance_from_center <= gv.radius:
            current_coordinates = temp_walker_list[i].current_coordinates
        else:
            current_coordinates = temp_walker_list[i].current_ball_center

        previous_distance = 0.0
        previous_ball_key = 0
        current_distance = 0.0
        current_ball_key = 0
        for j in range(balls.shape[0]):
            ball_center = balls[j][0:gv.num_cvs].tolist()
            previous_distance_from_center = calculate_distance_from_center(ball_center, previous_coordinates)
            current_distance_from_center = calculate_distance_from_center(ball_center, current_coordinates)
            if previous_distance == 0.0:
                previous_distance = previous_distance_from_center
                previous_ball_key = j
            else:
                if previous_distance_from_center < previous_distance:
                    previous_distance = previous_distance_from_center
                    previous_ball_key = j
            if current_distance == 0.0:
                current_distance = current_distance_from_center
                current_ball_key = j
            else:
                if current_distance_from_center < current_distance:
                    current_distance = current_distance_from_center
                    current_ball_key = j
        transition_matrix[previous_ball_key][current_ball_key] += temp_walker_list[i].weight

    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = (transition_matrix[i][j] + transition_matrix[j][i]) / 2.0
    for i in range(new_transition_matrix.shape[0]):
        row_sum = np.sum(new_transition_matrix[i])
        for j in range(new_transition_matrix.shape[1]):
            if row_sum != 0.0:
                new_transition_matrix[i][j] /= row_sum

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evectors = evectors[:, idx]
    eq_vector = abs(np.real(evectors[:, 0]))
    eq_vec_diag_matrix = np.diag(eq_vector)
    inv_eq_vec_diag_matrix = np.zeros((eq_vec_diag_matrix.shape[0], eq_vec_diag_matrix.shape[0]))
    for i in range(inv_eq_vec_diag_matrix.shape[0]):
        if eq_vec_diag_matrix[i][i] != 0.0:
            inv_eq_vec_diag_matrix[i][i] = 1.0 / eq_vec_diag_matrix[i][i]
    symmetric_transition_matrix = np.dot(np.sqrt(eq_vec_diag_matrix),
                                         np.dot(new_transition_matrix, np.sqrt(inv_eq_vec_diag_matrix)))

    final_evalues, final_evectors = np.linalg.eig(symmetric_transition_matrix)
    idx = abs(final_evalues).argsort()[::-1]
    final_evalues = np.real(final_evalues[idx])
    final_evectors = np.real(final_evectors[:, idx])

    num_clusters = gv.num_clusters
    for i in range(len(final_evalues)):
        if abs(final_evalues[i]) < 1.0e-10:
            final_evalues[i] = 0.0
    second_evector = final_evectors[:, 1]
    second_evector = second_evector.reshape(second_evector.shape[0], 1)
    log_second_evector = np.zeros((second_evector.shape[0], 1))
    for i in range(second_evector.shape[0]):
        if second_evector[i] < 0.0:
            log_second_evector[i] = -np.log(-second_evector[i])
        elif second_evector[i] == 0.0 or second_evector[i] == 1.0:
            log_second_evector[i] = 0.0
        else:
            log_second_evector[i] = np.log(second_evector[i])


    '''
    sorted_second_evector = np.sort(second_evector, axis=0)
    second_evector_order = np.ndarray.argsort(second_evector)
    num_balls = int(np.ceil(len(sorted_second_evector) / num_clusters))
    array_of_clusters = [sorted_second_evector[i:i + num_balls] for i in
                         range(0, len(sorted_second_evector), num_balls)]
    array_of_orderings = [second_evector_order[i:i + num_balls] for i in range(0, len(second_evector_order), num_balls)]
    num_clusters = len(array_of_clusters)
    '''

    matrix = np.hstack((balls, log_second_evector))

    while True:
        try:
            centroids, labels = kmeans2(matrix, num_clusters, minit='points', iter=100, missing='raise')
            break
        except ClusterError:
            num_clusters -= 1

    os.chdir(gv.main_directory + '/WE')
    with open('dunn_index_' + str(step_num + 1) + '.txt', 'w') as dunn_index_f:
        labeled_matrix = np.zeros((matrix.shape[0], matrix.shape[1] + 1))
        labeled_matrix[:,0:matrix.shape[1]] = matrix
        labeled_matrix[:,matrix.shape[1]] = labels
        print >>dunn_index_f, dunn(labeled_matrix)
    f = open('ball_clustering_' + str(step_num + 1) + '.txt', 'w')

    '''
    for i in range(num_clusters):
        first = 0
        cluster = array_of_clusters[i]
        ordering = array_of_orderings[i]
        for j in range(cluster.shape[0]):
            if first == 0:
                first += 1
                ref_ball_center = balls[ordering[j], 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ref_ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[ordering[j], 0]))
                ball_cluster.append(second_evector[ordering[j]])
                ball_cluster.append(final_evectors[ordering[j], 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
            else:
                ball_center = balls[ordering[j], 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[ordering[j], 0]))
                ball_cluster.append(second_evector[ordering[j]])
                ball_cluster.append(final_evectors[ordering[j], 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
    '''

    for i in range(num_clusters):
        first = 0
        for j in range(balls.shape[0]):
            if labels[j] == i and first == 0:
                first += 1
                ref_ball_center = balls[j, 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ref_ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[j, 0]))
                ball_cluster.append(log_second_evector[j, 0])
                ball_cluster.append(final_evectors[j, 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
                balls[j][gv.num_cvs+2] -= 1
            elif labels[j] == i and first != 0:
                ball_center = balls[j, 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[j, 0]))
                ball_cluster.append(log_second_evector[j, 0])
                ball_cluster.append(final_evectors[j, 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
                balls[j][gv.num_cvs+2] -= 1
    f.close()

    np.savetxt('evalues_' + str(step_num + 1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num + 1) + '.txt', final_evectors, fmt=' %1.10e')
    np.savetxt('transition_matrix_' + str(step_num + 1) + '.txt', symmetric_transition_matrix, fmt=' %1.10e')


def resampling_for_sc(walker_list, temp_walker_list, balls, ball_to_walkers, ball_clusters_list):
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    occupied_indices = np.zeros(gv.max_num_balls*gv.num_walkers_for_sc, int)
    excess_index = gv.total_num_walkers
    vacant_walker_indices = []
    for current_cluster in ball_clusters_list:
        if len(ball_clusters_list[current_cluster]) > 0:
            num_occupied_balls += 1

            num_bins = len(ball_clusters_list[current_cluster])
            if num_bins > gv.num_balls_for_sc:
                num_bins = gv.num_balls_for_sc
            bins = ball_clusters_list[current_cluster][0:num_bins]

            target_num_walkers = int(np.floor(float(gv.num_walkers_for_sc)/num_bins))
            remainder = gv.num_walkers_for_sc-target_num_walkers*num_bins

            for b, ball_center in enumerate(bins):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last bin if there are any
                if remainder != 0 and b == num_bins-1:
                    target_num_walkers += remainder

                weights_bin = []
                indices_bin = []
                for walker_index in ball_to_walkers[ball_center]:
                    weights_bin.append(temp_walker_list[walker_index].weight)
                    indices_bin.append(temp_walker_list[walker_index].global_index)
                # reset ball_to_walkers
                ball_to_walkers[ball_center] = []
                if b == num_bins-1:
                    for ball in ball_clusters_list[current_cluster][num_bins:]:
                        for walker_index in ball_to_walkers[ball]:
                            weights_bin.append(temp_walker_list[walker_index].weight)
                            indices_bin.append(temp_walker_list[walker_index].global_index)
                        # reset ball_to_walkers
                        ball_to_walkers[ball] = []

                weights_array = np.array(weights_bin)
                walker_indices = np.argsort(-weights_array)
                temp_indices = indices_bin
                # sorted indices based on descending order of weights
                indices_bin = [temp_indices[i] for i in walker_indices]

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
                        xy_weight = x_weight+y_weight
                        p = np.random.random()
                        # swap x and y
                        if p < y_weight/xy_weight:
                            temp = x
                            x = y
                            y = temp
                        weights[x] = xy_weight
                        if y not in new_indices:
                            vacant_walker_indices.append(y)

                for ni, global_index in enumerate(new_indices):
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[ni]
                        ball_to_walkers[ball_center].append(global_index)
                        ball_key = walker_list[global_index].ball_key
                        balls[ball_key][gv.num_cvs+2] += 1
                        directory = gv.main_directory + '/WE/walker' + str(global_index)
                        os.chdir(directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[ni] + '\n')
                        f.close()
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(walker_list[global_index])
                        ball_to_walkers[ball_center].append(new_index)
                        ball_key = walker_list[global_index].ball_key
                        balls[ball_key][gv.num_cvs+2] += 1
                        old_directory = gv.main_directory + '/WE/walker' + str(global_index)
                        new_directory = gv.main_directory + '/WE/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % walker_list[new_index].weight + '\n')
                        f.close()

    if excess_index - num_occupied_balls*gv.num_walkers_for_sc != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    if num_occupied_balls*gv.num_walkers_for_sc >= gv.total_num_walkers:
        for i in range(num_occupied_balls*gv.num_walkers_for_sc, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(num_occupied_balls*gv.num_walkers_for_sc, gv.total_num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= num_occupied_balls*gv.num_walkers_for_sc:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/WE')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))

    while len(vacant_walker_indices) > 0:
        vacant_walker_indices.pop()
    gv.num_occupied_balls = num_occupied_balls
    gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers_for_sc


def resampling(walker_list, temp_walker_list, balls, ball_to_walkers, vacant_walker_indices):
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    occupied_indices = np.zeros(gv.max_num_balls*gv.num_walkers, int)
    excess_index = gv.total_num_walkers
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs+2]) > 0:
            num_occupied_balls += 1
            current_ball_center = balls[current_ball][0:gv.num_cvs].tolist()
            initial_weights = [temp_walker_list[i].weight for i in ball_to_walkers[tuple(current_ball_center)]]
            initial_weights_array = np.array(initial_weights)
            walker_indices = np.argsort(-initial_weights_array)
            initial_indices = [temp_walker_list[i].global_index for i in ball_to_walkers[tuple(current_ball_center)]]
            temp_initial_indices = initial_indices
            # sorted indices based on descending order of weights
            initial_indices = [temp_initial_indices[i] for i in walker_indices]

            num_bins = 1
            true_num_bins = 1
            bins = [0]
            num_walkers_bin = [len(initial_indices)]

            if gv.enhanced_sampling_flag == 1:
                distance_from_center_list = [temp_walker_list[i].current_distance_from_center for i in
                                             ball_to_walkers[tuple(current_ball_center)]]
                std = np.sqrt(np.var(distance_from_center_list))
                if std != 0.0:
                    num_bins = int(np.ceil(gv.radius/std))+2
                    true_num_bins = 0
                    bins = []
                    num_walkers_bin = []
                    for nb in range(num_bins):
                        num_walkers = 0
                        for ii in initial_indices:
                            distance = temp_walker_list[ii].current_distance_from_center
                            if nb == 0:
                                if distance <= std:
                                    num_walkers += 1
                            elif nb == num_bins-1:
                                if nb*std < distance:
                                    num_walkers += 1
                            else:
                                if nb*std < distance <= (nb+1)*std:
                                    num_walkers += 1
                        if num_walkers != 0:
                            true_num_bins += 1
                            bins.append(nb)
                            num_walkers_bin.append(num_walkers)
                    if true_num_bins <= 1:
                        num_bins = 1
                        true_num_bins = 1
                        bins = [0]
                        num_walkers_bin = [len(initial_indices)]
                        std = 0.0

            target_num_walkers = int(np.floor(float(gv.num_walkers)/true_num_bins))
            remainder = gv.num_walkers-target_num_walkers*true_num_bins
            # reset ball_to_walkers
            ball_to_walkers[tuple(current_ball_center)] = []

            for b, bin_index in enumerate(bins):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last bin if there are any
                if remainder != 0 and b == (true_num_bins-1):
                    target_num_walkers += remainder

                weights_bin = [float] * num_walkers_bin[b]
                indices_bin = [int] * num_walkers_bin[b]

                if gv.enhanced_sampling_flag != 1 or (gv.enhanced_sampling_flag == 1 and std == 0.0):
                    weights_bin = initial_weights
                    indices_bin = initial_indices

                elif gv.enhanced_sampling_flag == 1 and std != 0.0:
                    k = 0
                    for j in initial_indices:
                        distance = temp_walker_list[j].current_distance_from_center
                        if bin_index == 0:
                            if distance <= std:
                                weights_bin[k] = temp_walker_list[j].weight
                                indices_bin[k] = temp_walker_list[j].global_index
                                k += 1
                        elif bin_index == num_bins - 1:
                            if bin_index * std < distance:
                                weights_bin[k] = temp_walker_list[j].weight
                                indices_bin[k] = temp_walker_list[j].global_index
                                k += 1
                        else:
                            if bin_index * std < distance <= (bin_index + 1) * std:
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
                        xy_weight = x_weight+y_weight
                        p = np.random.random()
                        # swap x and y
                        if p < y_weight/xy_weight:
                            temp = x
                            x = y
                            y = temp
                        weights[x] = xy_weight
                        if y not in new_indices:
                            vacant_walker_indices.append(y)
                            # remove walker y directory
                            os.chdir(gv.main_directory + '/WE')
                            os.system('rm -rf walker' + str(y))

                if b == 0:  # reset balls
                    balls[current_ball][gv.num_cvs+2] = 0
                for ni, global_index in enumerate(new_indices):
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[ni]
                        ball_to_walkers[tuple(current_ball_center)].append(global_index)
                        directory = gv.main_directory + '/WE/walker' + str(global_index)
                        os.chdir(directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % new_weights[ni] + '\n')
                        f.close()
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(walker_list[global_index])
                        ball_to_walkers[tuple(current_ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/WE/walker' + str(global_index)
                        new_directory = gv.main_directory + '/WE/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % walker_list[new_index].weight + '\n')
                        f.close()
                    balls[current_ball][gv.num_cvs+2] += 1

    if excess_index-num_occupied_balls*gv.num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    if num_occupied_balls*gv.num_walkers >= gv.total_num_walkers:
        for i in range(num_occupied_balls*gv.num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.total_num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(num_occupied_balls*gv.num_walkers, gv.total_num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= num_occupied_balls*gv.num_walkers:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])
                # rename the directory with name 'i' to 'new_index'
                os.chdir(gv.main_directory + '/WE')
                os.system('mv walker' + str(i) + ' walker' + str(new_index))

    while len(vacant_walker_indices) > 0:
        vacant_walker_indices.pop()
    gv.num_occupied_balls = num_occupied_balls
    gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers


def print_status(step_num, walker_list, balls, ball_to_walkers, ball_clusters_list, key_to_ball):
    os.chdir(gv.main_directory + '/WE')
    total_weight = 0.0
    f = open('total_weight_on_each_ball_' + str(step_num + 1) + '.txt', 'w')
    for current_ball in range(balls.shape[0]):
        ball_center = balls[current_ball][0:gv.num_cvs].tolist()
        weights = [walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
        total_weight += np.sum(weights)
        ball_center_weights = copy.deepcopy(ball_center)
        ball_center_weights.append(np.sum(weights))
        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_center_weights)))
        f.write('\n')
        # reset walkers and number of walkers that belong in each ball
        balls[current_ball][gv.num_cvs+2] = 0
        ball_to_walkers[tuple(ball_center)] = []
        key_to_ball[tuple(ball_center)] = []
        ball_clusters_list[tuple(ball_center)] = []
    f.close()

    # verify that total weight of all balls is 1.0
    f = open('total_weight.txt', 'a')
    f.write(str(step_num + 1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + '\n')
