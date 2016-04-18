import numpy as np
import os
import shutil
import copy
from scipy import special
from scipy.cluster.vq import kmeans2, ClusterError
import walker
import global_variables as gv
import check_state_function
import parameters as p
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.covariance import EllipticEnvelope


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        if gv.angle_cvs[i] == 0:
            distance += (values[i] - center[i]) ** 2
        else:
            if values[i] - center[i] > 180.0:
                distance += (values[i] - center[i] - 360.0) ** 2
            elif values[i] - center[i] < -180.0:
                distance += (values[i] - center[i] + 360.0) ** 2
            else:
                distance += (values[i] - center[i]) ** 2
    if abs(distance) < 1.0e-10:
        distance = 0.0
    return np.sqrt(distance)


def set_parameters():
    gv.main_directory = p.main_directory
    gv.initial_configuration_directory = p.initial_configuration_directory
    gv.simulation_flag = p.simulation_flag
    gv.balls_flag = p.balls_flag
    gv.rate_flag = p.rate_flag
    gv.num_states = p.num_states
    gv.enhanced_sampling_flag = p.enhanced_sampling_flag
    gv.num_balls_limit = p.num_balls_limit
    gv.radius = p.radius
    gv.num_walkers = p.num_walkers
    gv.num_cvs = p.num_cvs
    gv.lower_bound = p.lower_bound
    gv.upper_bound = p.upper_bound
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
        gv.num_occupied_big_clusters = p.num_occupied_big_clusters
        gv.num_occupied_small_clusters = p.num_occupied_small_clusters
        gv.num_balls_for_sc = p.num_balls_for_sc
        gv.num_clusters = p.num_clusters
        gv.num_walkers_for_sc = p.num_walkers_for_sc

    ball_volume = (np.pi ** (gv.num_cvs / 2) * gv.radius ** gv.num_cvs) / special.gamma((gv.num_cvs / 2) + 1)
    if ball_volume != 0.0:
        max_num_balls = int(np.floor((gv.upper_bound - gv.lower_bound) ** gv.num_cvs / ball_volume))*2
    if max_num_balls < gv.num_balls_limit:
        gv.num_balls_limit = max_num_balls
    print 'max # of balls (n_b) = ' + str(gv.num_balls_limit)
    gv.current_num_balls = 0
    if gv.enhanced_sampling_flag == 2 and gv.simulation_flag != 0 and \
                            gv.num_occupied_small_clusters+gv.num_occupied_big_clusters != 0:
        gv.total_num_walkers = gv.num_occupied_small_clusters*gv.num_walkers \
                               + gv.num_occupied_big_clusters*gv.num_walkers_for_sc
        gv.num_occupied_small_clusters = 0
        gv.num_occupied_big_clusters = 0
    else:
        gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers
    gv.sc_performed = 0


def initialize(input_initial_values_file, walker_list, temp_walker_list, balls, ball_to_walkers, vacant_walker_indices):
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
                initial_state = check_state_function.check_state_function(initial_values)
                print initial_state
            for i in range(n * gv.num_walkers, (n + 1) * gv.num_walkers):
                walker_list[i].set(initial_values, initial_weight)
                if gv.rate_flag == 1:
                    walker_list[i].state = initial_state
        f.close()

        os.system('mkdir CAS')
        os.chdir(gv.main_directory + '/CAS')
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            shutil.copytree(gv.initial_configuration_directory, walker_directory)

    elif gv.simulation_flag == 1 or gv.simulation_flag == 2:  # restarting simulation in the middle of simulation
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            f = open('trajectory.txt', 'r')
            lines = f.readlines()[-2:]
            previous_line = lines[0].strip().split()
            previous_coordinates = [float(entry) for entry in previous_line]
            current_line = lines[1].strip().split()
            current_coordinates = [float(entry) for entry in current_line]
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            f.close()
            f = open('ball_trajectory.txt', 'r')
            lines = f.readlines()[-2:]
            previous_line = lines[0].strip().split()
            previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
            current_line = lines[1].strip().split()
            current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = float(current_line[gv.num_cvs])
            f.close()
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                         current_ball_center)
            if gv.rate_flag == 1:
                walker_list[i].state = int(current_line[-1])
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
                gv.current_num_balls = balls.shape[0]

    elif gv.simulation_flag == 3:  # restarting simulation in the middle of binning
        for i in range(gv.total_num_walkers):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            f = open('trajectory.txt', 'r')
            lines = f.readlines()[-2:]
            previous_line = lines[0].strip().split()
            previous_coordinates = [float(entry) for entry in previous_line]
            current_line = lines[1].strip().split()
            current_coordinates = [float(entry) for entry in current_line]
            walker_list[i].previous_coordinates = previous_coordinates
            walker_list[i].current_coordinates = current_coordinates
            f.close()
            num_lines = sum(1 for line in open('ball_trajectory.txt'))
            # if walker is already binned to a ball, delete the binning and have binning start from scratch
            if num_lines > gv.initial_step_num:
                os.system('sed -i \'$d\' ball_trajectory.txt')
            f = open('ball_trajectory.txt', 'r')
            lines = f.readlines()[-2:]
            previous_line = lines[0].strip().split()
            previous_ball_center = [float(entry) for entry in previous_line[0:gv.num_cvs]]
            current_line = lines[1].strip().split()
            current_ball_center = [float(entry) for entry in current_line[0:gv.num_cvs]]
            walker_list[i].previous_ball_center = previous_ball_center
            walker_list[i].current_ball_center = current_ball_center
            walker_list[i].radius = float(current_line[gv.num_cvs])
            f.close()
            walker_list[i].previous_distance_from_center = calculate_distance_from_center(previous_coordinates,
                                                                                          previous_ball_center)
            walker_list[i].current_distance_from_center = calculate_distance_from_center(current_coordinates,
                                                                                         current_ball_center)
            if gv.rate_flag == 1:
                walker_list[i].state = int(current_line[-1])
            if gv.balls_flag == 1:
                os.chdir(gv.main_directory + '/CAS')
                balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
                gv.current_num_balls = balls.shape[0]

    elif gv.simulation_flag == 4:  # restarting simulation in the middle of resampling
        total_weight = 0.0
        total_num_walkers = 0
        occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers*2, int)
        previous_ball_to_walkers = {}
        os.chdir(gv.main_directory + '/CAS')
        previous_balls_weights = np.loadtxt('total_weight_on_each_ball_' + str(gv.initial_step_num) + '.txt')

        for i in range(gv.last_walker+1):
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            # if all of the files exist in the walker folder, it is a complete walker
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

                ball_trajectory = np.loadtxt('ball_trajectory.txt')
                previous_ball = ball_trajectory[-2].tolist()
                previous_ball_key = previous_ball[gv.num_cvs+1]
                previous_ball_center = previous_ball[0:gv.num_cvs]
                if previous_balls_weights[previous_ball_key][-1]-weight >= 0.0:
                    walker_list[i].weight = weight
                    total_weight += weight
                    previous_balls_weights[previous_ball_key][-1] -= weight
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
                    previous_distance_from_center = calculate_distance_from_center(previous_coordinates, previous_ball_center)
                    current_distance_from_center = calculate_distance_from_center(current_coordinates, current_ball_center)
                    walker_list[i].previous_distance_from_center = previous_distance_from_center
                    walker_list[i].current_distance_from_center = current_distance_from_center

                    temp_walker_list[i] = walker.Walker(previous_coordinates, current_coordinates, i,
                                                        current_ball_radius, previous_ball_center, current_ball_center,
                                                        current_ball_key, previous_distance_from_center,
                                                        current_distance_from_center, 0, weight, current_state)
                    total_num_walkers += 1
                    occupied_indices[i] = 1
                else:
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('rm -rf walker' + str(i))
                    vacant_walker_indices.append(i)
            # otherwise, it is an incomplete walker that needs missing files
            else:
                if os.path.isdir(walker_directory):
                    os.chdir(gv.main_directory + '/CAS')
                    os.system('rm -rf walker' + str(i))
                vacant_walker_indices.append(i)

        # create new walkers for the remaining weights
        excess_index = gv.last_walker+1
        for i in range(previous_balls_weights.shape[0]):
            if previous_balls_weights[i][-1] > 0.0:
                previous_ball_center = previous_balls_weights[i][0:gv.num_cvs].tolist()
                if tuple(previous_ball_center) in previous_ball_to_walkers:
                    reference_walker = previous_ball_to_walkers[tuple(previous_ball_center)][0]
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

                ball_trajectory = np.loadtxt('ball_trajectory.txt')
                current_ball_center = ball_trajectory[-1][0:gv.num_cvs].tolist()
                walker_list[walker_index].previous_ball_center = previous_ball_center
                walker_list[walker_index].current_ball_center = current_ball_center
                ball_to_walkers[tuple(current_ball_center)].append(walker_index)
                os.system('sed -i \'$d\' ball_trajectory.txt')
                trajectory = np.loadtxt('trajectory.txt')
                previous_coordinates = trajectory[-2].tolist()
                current_coordinates = trajectory[-1].tolist()
                walker_list[walker_index].previous_coordinates = previous_coordinates
                walker_list[walker_index].current_coordinates = current_coordinates
                current_state = ball_trajectory[-1][-1]
                current_ball_key = ball_trajectory[-1][gv.num_cvs+1]
                current_ball_radius = ball_trajectory[-1][gv.num_cvs]
                walker_list[walker_index].state = current_state
                walker_list[walker_index].ball_key = current_ball_key
                walker_list[walker_index].radius = current_ball_radius
                previous_distance_from_center = calculate_distance_from_center(previous_coordinates, previous_ball_center)
                current_distance_from_center = calculate_distance_from_center(current_coordinates, current_ball_center)
                walker_list[walker_index].previous_distance_from_center = previous_distance_from_center
                walker_list[walker_index].current_distance_from_center = current_distance_from_center

                temp_walker_list[walker_index] = walker.Walker(previous_coordinates, current_coordinates, walker_index,
                                                               current_ball_radius, previous_ball_center,
                                                               current_ball_center, current_ball_key,
                                                               previous_distance_from_center,
                                                               current_distance_from_center, 0, weight, current_state)
                total_num_walkers += 1
                occupied_indices[walker_index] = 1

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
        gv.total_num_walkers = total_num_walkers
        gv.num_occupied_balls = 0
        gv.current_num_balls = 0
        gv.first_walker = 0
        gv.last_walker = total_num_walkers-1
        if gv.balls_flag == 1:
            os.chdir(gv.main_directory + '/CAS')
            balls = np.loadtxt('balls_' + str(gv.initial_step_num) + '.txt')
            gv.current_num_balls = balls.shape[0]
    return balls


def binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers, key_to_ball):
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)
    flux = np.zeros((gv.num_states, gv.num_states))
    flux_num_walkers = np.zeros((gv.num_states, gv.num_states))
    walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    start = 0  # indicates whether we are dealing with the very first walker or not

    for i in walker_indices:
        # first, go to walker directory i
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)

        # then, obtain new coordinates' values
        if os.path.exists(walker_directory + '/coordinates.out'):
            coordinates = np.loadtxt('coordinates.out')
            if gv.num_cvs > 1:
                new_coordinates = coordinates.tolist()
            else:
                new_coordinates = [float(coordinates)]
            rm_command = 'rm -rf *.out'
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
            state = check_state_function.check_state_function(new_coordinates)
            if walker_list[i].state != -1 and state == -1:
                state = walker_list[i].state
            if walker_list[i].state != -1 and state != -1:
                flux[walker_list[i].state, state] += walker_list[i].weight
                flux_num_walkers[walker_list[i].state, state] += 1
        else:
            state = -1

        inside = 0  # indicates whether we are dealing with the very first walker or not
        # if we're dealing with the very first walker, create the very first ball for the walker
        if (gv.balls_flag == 0 and start == 0) or (gv.balls_flag == 1 and start == 0 and step_num == 0):
            start += 1
            inside += 1
            current_ball_center = [coordinate for coordinate in new_coordinates]
            center_r_key_num = copy.deepcopy(current_ball_center)
            center_r_key_num.append(gv.radius)
            center_r_key_num.append(gv.current_num_balls)
            center_r_key_num.append(1)
            balls[gv.current_num_balls] = np.asarray(center_r_key_num)
            ball_to_walkers[tuple(current_ball_center)] = [i]
            key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
            temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                previous_ball_center, current_ball_center, gv.current_num_balls,
                                                previous_distance_from_center, 0.0, initial_step_num, weight, state)
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
                if j == 0:
                    distance = distance_from_center
                    ball_key = j
                else:
                    if distance_from_center < distance:
                        distance = distance_from_center
                        ball_key = j

            # walker is inside some ball or is not but needs to be binned to the nearest ball
            if inside != 0 or (inside == 0 and gv.current_num_balls == gv.num_balls_limit):
                balls[ball_key][gv.num_cvs+2] += 1
                current_ball_center = balls[ball_key][0:gv.num_cvs].tolist()
                ball_to_walkers[tuple(current_ball_center)].append(i)
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, ball_key,
                                                    previous_distance_from_center, distance_from_center,
                                                    initial_step_num, weight, state)

            # or walker does not belong in any ball -> create a new ball
            elif inside == 0 and gv.current_num_balls < gv.num_balls_limit:
                current_ball_center = [coordinate for coordinate in new_coordinates]
                center_r_key_num = copy.deepcopy(current_ball_center)
                center_r_key_num.append(gv.radius)
                center_r_key_num.append(gv.current_num_balls)
                center_r_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                ball_to_walkers[tuple(current_ball_center)] = [i]
                key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                                    previous_ball_center, current_ball_center, gv.current_num_balls,
                                                    previous_distance_from_center, 0.0, initial_step_num, weight, state)
                gv.current_num_balls += 1

        # finally, write the new ball on the trajectory file
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

    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('balls_' + str(step_num + 1) + '.txt', balls, fmt=' %+1.5f')
    if gv.rate_flag == 1:
        np.savetxt('flux_' + str(step_num + 1) + '.txt', flux, fmt=' %1.5e')
        np.savetxt('flux_num_walkers_' + str(step_num + 1) + '.txt', flux_num_walkers, fmt=' %d')

    transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(gv.total_num_walkers):
        previous_coordinates = temp_walker_list[i].previous_coordinates
        previous_distance = 0.0
        previous_ball_key = 0
        for j in range(balls.shape[0]):
            ball_center = balls[j][0:gv.num_cvs].tolist()
            previous_distance_from_center = calculate_distance_from_center(ball_center, previous_coordinates)
            if j == 0:
                previous_distance = previous_distance_from_center
                previous_ball_key = j
            else:
                if previous_distance_from_center < previous_distance:
                    previous_distance = previous_distance_from_center
                    previous_ball_key = j
        transition_matrix[previous_ball_key][temp_walker_list[i].ball_key] += temp_walker_list[i].weight

    # transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
    # canonical ensemble. equation is from Prinz, et al JCP (2011).
    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = transition_matrix[i][j] + transition_matrix[j][i]

    row_sum = np.sum(new_transition_matrix, axis=1)
    for i in range(new_transition_matrix.shape[0]):
        if row_sum[i] != 0.0:
            new_transition_matrix[i, :] /= row_sum[i]
    np.savetxt('transition_matrix_' + str(step_num + 1) + '.txt', new_transition_matrix, fmt=' %1.10e')

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evalues = evalues[idx]
    final_evalues = np.real(evalues)
    evectors = evectors[:, idx]
    final_evectors = np.real(evectors)
    np.savetxt('evalues_' + str(step_num + 1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num + 1) + '.txt', final_evectors, fmt=' %1.10e')
    return balls


def threshold_binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers, key_to_ball):
    initial_weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    initial_weights_array = np.array(initial_weights)
    flux = np.zeros((gv.num_states, gv.num_states))
    flux_num_walkers = np.zeros((gv.num_states, gv.num_states))
    walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    if gv.static_threshold_flag == 0:
        new_threshold_values = gv.threshold_values
        ref_walker_binning_value = 0
        ref_walker_properties_value = 0.0
        ref_walker_index = 0
    else:
        ref_walker_indices = []

    for i in walker_indices:
        # first, go to walker directory i
        walker_directory = gv.main_directory + '/CAS/walker' + str(i)
        os.chdir(walker_directory)

        # then, obtain new coordinates' values
        if os.path.exists(walker_directory + '/coordinates.out'):
            coordinates = np.loadtxt('coordinates.out')
            if gv.num_cvs > 1:
                new_coordinates = coordinates.tolist()
            else:
                new_coordinates = [float(coordinates)]
            rm_command = 'rm -rf *.out'
            os.system(rm_command)
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
            state = check_state_function.check_state_function(new_coordinates)
            if walker_list[i].state != -1 and state == -1:
                state = walker_list[i].state
            if walker_list[i].state != -1 and state != -1:
                flux[walker_list[i].state, state] += walker_list[i].weight
                flux_num_walkers[walker_list[i].state, state] += 1
        else:
            state = -1

        temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, gv.radius,
                                            previous_ball_center, previous_ball_center, gv.current_num_balls,
                                            previous_distance_from_center, previous_distance_from_center,
                                            initial_step_num, weight, state)

        if gv.static_threshold_flag == 0:
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
            else:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        walker_binning_value += 1
                    walker_properties_value += (properties_to_keep_track[m]-gv.threshold_values[m])
            if i == 0:
                ref_walker_index = i
                ref_walker_binning_value = walker_binning_value
                ref_walker_properties_value = walker_properties_value
            else:
                if walker_binning_value == ref_walker_binning_value \
                        and walker_properties_value <= ref_walker_properties_value:
                    ref_walker_index = i
                    ref_walker_binning_value = walker_binning_value
                    ref_walker_properties_value = walker_properties_value
                elif walker_binning_value < ref_walker_binning_value:
                    ref_walker_index = i
                    ref_walker_binning_value = walker_binning_value
                    ref_walker_properties_value = walker_properties_value
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
                ref_walker_indices.append(i)

    if gv.static_threshold_flag == 0:
        new_coordinates = temp_walker_list[ref_walker_index].current_coordinates
        state = temp_walker_list[ref_walker_index].state
        center_r_key_num = copy.deepcopy(new_coordinates)
        center_r_key_num.append(gv.radius)
        center_r_key_num.append(gv.current_num_balls)
        center_r_key_num.append(gv.num_walkers)
        balls[gv.current_num_balls] = np.asarray(center_r_key_num)
        ball_to_walkers[tuple(new_coordinates)] = []
        key_to_ball[tuple(new_coordinates)] = gv.current_num_balls

        center_r_key_state = copy.deepcopy(new_coordinates)
        center_r_key_state.append(gv.radius)
        center_r_key_state.append(gv.current_num_balls)
        center_r_key_state.append(state)
        ref_walker_directory = gv.main_directory + '/CAS/walker' + str(ref_walker_index)
        os.chdir(ref_walker_directory)
        f = open('ball_trajectory.txt', 'a')
        f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
        f.write('\n')
        f.close()
        f = open('trajectory.txt', 'a')
        f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
        f.write('\n')
        f.close()

        for i in walker_indices:
            temp_walker_list[i].current_coordinates = new_coordinates
            temp_walker_list[i].current_ball_center = new_coordinates
            temp_walker_list[i].ball_key = gv.current_num_balls
            temp_walker_list[i].current_distance_from_center = 0.0
            ball_to_walkers[tuple(new_coordinates)].append(i)
            if i != ref_walker_index:
                os.chdir(gv.main_directory + '/CAS')
                os.system('rm -rf walker' + str(i))
                old_directory = gv.main_directory + '/CAS/walker' + str(ref_walker_index)
                new_directory = gv.main_directory + '/CAS/walker' + str(i)
                shutil.copytree(old_directory, new_directory)

        gv.current_num_balls += 1

    else:
        walker_indices_list = walker_indices.tolist()
        if len(ref_walker_indices) > 0:
            ref_walker_index = ref_walker_indices[0]
            current_ball_center = temp_walker_list[ref_walker_index].current_coordinates
            center_r_key_num = copy.deepcopy(current_ball_center)
            center_r_key_num.append(gv.radius)
            center_r_key_num.append(gv.current_num_balls)
            center_r_key_num.append(0)
            balls[gv.current_num_balls] = np.asarray(center_r_key_num)
            ball_to_walkers[tuple(current_ball_center)] = []
            key_to_ball[tuple(current_ball_center)] = gv.current_num_balls

            for i in ref_walker_indices:
                walker_indices_list.remove(i)
                state = temp_walker_list[i].state
                new_coordinates = temp_walker_list[i].current_coordinates
                temp_walker_list[i].current_ball_center = current_ball_center
                temp_walker_list[i].ball_key = gv.current_num_balls
                temp_walker_list[i].current_distance_from_center = \
                    calculate_distance_from_center(current_ball_center, temp_walker_list[i].current_coordinates)
                center_r_key_state = copy.deepcopy(current_ball_center)
                center_r_key_state.append(gv.radius)
                center_r_key_state.append(gv.current_num_balls)
                center_r_key_state.append(state)
                walker_directory = gv.main_directory + '/CAS/walker' + str(i)
                os.chdir(walker_directory)
                f = open('ball_trajectory.txt', 'a')
                f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
                f.write('\n')
                f.close()
                f = open('trajectory.txt', 'a')
                f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
                f.write('\n')
                f.close()
                balls[gv.current_num_balls][gv.num_cvs+2] += 1
                ball_to_walkers[tuple(current_ball_center)].append(i)

            gv.current_num_balls += 1

        start = 0

        for i in walker_indices_list:
            state = temp_walker_list[i].state
            new_coordinates = temp_walker_list[i].current_coordinates
            inside = 0
            if start == 0:
                start += 1
                inside += 1
                current_ball_center = new_coordinates
                center_r_key_num = copy.deepcopy(current_ball_center)
                center_r_key_num.append(gv.radius)
                center_r_key_num.append(gv.current_num_balls)
                center_r_key_num.append(1)
                if gv.current_num_balls == 0:
                    balls[gv.current_num_balls] = np.asarray(center_r_key_num)
                else:
                    balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                ball_to_walkers[tuple(current_ball_center)] = [i]
                key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                temp_walker_list[i].current_ball_center = current_ball_center
                temp_walker_list[i].ball_key = gv.current_num_balls
                temp_walker_list[i].current_distance_from_center = 0.0
                gv.current_num_balls += 1

            distance = 0.0
            ball_key = 0
            if inside == 0:
                for j in range(balls.shape[0]):
                    current_ball_center = balls[j][0:gv.num_cvs].tolist()
                    distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                    if distance_from_center <= gv.radius or abs(distance_from_center - gv.radius) < 1.0e-10:
                        inside += 1
                    if j == 0:
                        distance = distance_from_center
                        ball_key = j
                    else:
                        if distance_from_center < distance:
                            distance = distance_from_center
                            ball_key = j

                if inside != 0 or (inside == 0 and gv.current_num_balls == gv.num_balls_limit):
                    balls[ball_key][gv.num_cvs+2] += 1
                    current_ball_center = balls[ball_key][0:gv.num_cvs].tolist()
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                    temp_walker_list[i].current_ball_center = current_ball_center
                    temp_walker_list[i].ball_key = gv.current_num_balls
                    temp_walker_list[i].current_distance_from_center = \
                        calculate_distance_from_center(current_ball_center, new_coordinates)

                elif inside == 0 and gv.current_num_balls < gv.num_balls_limit:
                    current_ball_center = [coordinate for coordinate in new_coordinates]
                    center_r_key_num = copy.deepcopy(current_ball_center)
                    center_r_key_num.append(gv.radius)
                    center_r_key_num.append(gv.current_num_balls)
                    center_r_key_num.append(1)
                    balls = np.append(balls, [np.asarray(center_r_key_num)], axis=0)
                    ball_to_walkers[tuple(current_ball_center)] = [i]
                    key_to_ball[tuple(current_ball_center)] = gv.current_num_balls
                    temp_walker_list[i].current_ball_center = current_ball_center
                    temp_walker_list[i].ball_key = gv.current_num_balls
                    temp_walker_list[i].current_distance_from_center = 0.0
                    gv.current_num_balls += 1

            current_ball_center = temp_walker_list[i].current_ball_center
            ball_key = temp_walker_list[i].ball_key
            center_r_key_state = copy.deepcopy(current_ball_center)
            center_r_key_state.append(gv.radius)
            center_r_key_state.append(ball_key)
            center_r_key_state.append(state)
            walker_directory = gv.main_directory + '/CAS/walker' + str(i)
            os.chdir(walker_directory)
            f = open('ball_trajectory.txt', 'a')
            f.write(' '.join(map(lambda coordinate: str(coordinate), center_r_key_state)))
            f.write('\n')
            f.close()
            f = open('trajectory.txt', 'a')
            f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
            f.write('\n')
            f.close()

    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('balls_' + str(step_num + 1) + '.txt', balls, fmt=' %+1.5f')
    if gv.rate_flag == 1:
        np.savetxt('flux_' + str(step_num + 1) + '.txt', flux, fmt=' %1.5e')
        np.savetxt('flux_num_walkers_' + str(step_num + 1) + '.txt', flux_num_walkers)
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

    transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(gv.total_num_walkers):
        previous_coordinates = temp_walker_list[i].previous_coordinates
        previous_distance = 0.0
        previous_ball_key = 0
        for j in range(balls.shape[0]):
            ball_center = balls[j][0:gv.num_cvs].tolist()
            previous_distance_from_center = calculate_distance_from_center(ball_center, previous_coordinates)
            if j == 0:
                previous_distance = previous_distance_from_center
                previous_ball_key = j
            else:
                if previous_distance_from_center < previous_distance:
                    previous_distance = previous_distance_from_center
                    previous_ball_key = j
        transition_matrix[previous_ball_key][temp_walker_list[i].ball_key] += temp_walker_list[i].weight

    # transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
    # canonical ensemble. equation is from Prinz, et al JCP (2011).
    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = transition_matrix[i][j] + transition_matrix[j][i]

    row_sum = np.sum(new_transition_matrix, axis=1)
    for i in range(new_transition_matrix.shape[0]):
        if row_sum[i] != 0.0:
            new_transition_matrix[i, :] /= row_sum[i]
    np.savetxt('transition_matrix_' + str(step_num + 1) + '.txt', new_transition_matrix, fmt=' %1.10e')

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evalues = evalues[idx]
    final_evalues = np.real(evalues)
    evectors = evectors[:, idx]
    final_evectors = np.real(evectors)
    np.savetxt('evalues_' + str(step_num + 1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num + 1) + '.txt', final_evectors, fmt=' %1.10e')
    return balls


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


def spectral_clustering(step_num, temp_walker_list, balls, ball_clusters_list):
    gv.sc_performed = 1
    transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(gv.total_num_walkers):
        previous_coordinates = temp_walker_list[i].previous_coordinates
        previous_distance = 0.0
        previous_ball_key = 0
        for j in range(balls.shape[0]):
            ball_center = balls[j][0:gv.num_cvs].tolist()
            previous_distance_from_center = calculate_distance_from_center(ball_center, previous_coordinates)
            if j == 0:
                previous_distance = previous_distance_from_center
                previous_ball_key = j
            else:
                if previous_distance_from_center < previous_distance:
                    previous_distance = previous_distance_from_center
                    previous_ball_key = j
        transition_matrix[previous_ball_key][temp_walker_list[i].ball_key] += temp_walker_list[i].weight

    # transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
    # canonical ensemble. equation is from Prinz, et al JCP (2011).
    new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(new_transition_matrix.shape[0]):
        for j in range(new_transition_matrix.shape[1]):
            new_transition_matrix[i][j] = transition_matrix[i][j] + transition_matrix[j][i]

    row_sum = np.sum(new_transition_matrix, axis=1)
    for i in range(new_transition_matrix.shape[0]):
        if row_sum[i] != 0.0:
            new_transition_matrix[i, :] /= row_sum[i]
    os.chdir(gv.main_directory + '/CAS')
    np.savetxt('transition_matrix_' + str(step_num + 1) + '.txt', new_transition_matrix, fmt=' %1.10e')

    '''
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
    '''

    evalues, evectors = np.linalg.eig(new_transition_matrix.T)
    idx = abs(evalues).argsort()[::-1]
    evalues = evalues[idx]
    final_evalues = np.real(evalues)
    evectors = evectors[:, idx]
    final_evectors = np.real(evectors)
    np.savetxt('evalues_' + str(step_num + 1) + '.txt', final_evalues, fmt=' %1.10e')
    np.savetxt('evectors_' + str(step_num + 1) + '.txt', final_evectors, fmt=' %1.10e')

    num_clusters = gv.num_clusters
    normalized_second_evector = np.zeros((final_evectors.shape[0], 1))
    for i in range(final_evectors.shape[0]):
        if final_evectors[i, 0] != 0.0:
            normalized_second_evector[i] = final_evectors[i, 1] / abs(final_evectors[i, 0])
        else:
            normalized_second_evector[i] = 0.0

    '''
    sorted_second_evector = np.sort(second_evector, axis=0)
    second_evector_order = np.ndarray.argsort(second_evector)
    num_balls = int(np.ceil(len(sorted_second_evector) / num_clusters))
    array_of_clusters = [sorted_second_evector[i:i + num_balls] for i in
                         range(0, len(sorted_second_evector), num_balls)]
    array_of_orderings = [second_evector_order[i:i + num_balls] for i in range(0, len(second_evector_order), num_balls)]
    num_clusters = len(array_of_clusters)
    '''

    matrix = normalized_second_evector  #np.hstack((balls, normalized_second_evector))
    clustering_matrix = matrix
    cont = True
    outlier_labels = np.ones(len(matrix)) * -1
    outliers_exist = 0
    while cont:
        while True:
            try:
                centroids, labels = kmeans2(clustering_matrix, num_clusters, minit='points', iter=100, missing='raise')
                labels = merge_with_outliers(outlier_labels, labels)
                break
            except ClusterError:
                num_clusters -= 1

        if num_clusters <= 1:
            gv.sc_performed = 0
            break
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
            if not cont:
                with open('dunn_index_' + str(step_num + 1) + '.txt', 'w') as dunn_index_f:
                    labeled_matrix = np.zeros((matrix.shape[0], matrix.shape[1] + 1))
                    labeled_matrix[:, 0:matrix.shape[1]] = matrix
                    labeled_matrix[:, matrix.shape[1]] = labels
                    print >>dunn_index_f, dunn(labeled_matrix)
                    print >>dunn_index_f, "The average silhouette_score is: %f" % silhouette_avg
                    for i in xrange(int(max(labels))+1):
                        print >>dunn_index_f, "The average silhouette score for cluster %d is: %f" % (i, np.mean(sample_silhouette_values[labels == i]))

    if gv.sc_performed == 1:
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

        if outliers_exist == 1:
            for i in range(num_clusters-1):
                first = 0
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

            cluster_num = num_clusters-1
            for j in range(balls.shape[0]):
                if labels[j] == num_clusters-1:
                    ball_center = balls[j, 0:gv.num_cvs].tolist()
                    ball_cluster = copy.deepcopy(ball_center)
                    ball_cluster.append(cluster_num)
                    ball_cluster.append(abs(final_evectors[j, 0]))
                    ball_cluster.append(final_evectors[j, 1])
                    ball_cluster.append(final_evectors[j, 2])
                    f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                    f.write('\n')
                    ball_clusters_list[tuple(ball_center)] = [tuple(ball_center)]
                    balls[j][gv.num_cvs + 2] -= 1
                    cluster_num += 1
        else:
            for i in range(num_clusters):
                first = 0
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
        f.close()


def resampling_for_sc(walker_list, temp_walker_list, balls, ball_to_walkers, ball_clusters_list, vacant_walker_indices):
    gv.sc_performed = 1
    num_occupied_big_clusters = 0
    num_occupied_small_clusters = 0
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers_for_sc*2, int)
    excess_index = gv.total_num_walkers
    for current_cluster in ball_clusters_list:
        if len(ball_clusters_list[current_cluster]) > 0:
            if len(ball_clusters_list[current_cluster]) > 1:
                num_occupied_big_clusters += 1
                initial_target_num_walkers = gv.num_walkers_for_sc
            else:
                num_occupied_small_clusters += 1
                initial_target_num_walkers = gv.num_walkers

            initial_weights = []
            initial_indices = []
            for ball_center in ball_clusters_list[current_cluster]:
                for walker_index in ball_to_walkers[ball_center]:
                    initial_weights.append(temp_walker_list[walker_index].weight)
                    initial_indices.append(temp_walker_list[walker_index].global_index)
                # reset ball_to_walkers and balls
                ball_to_walkers[ball_center] = []
                ball_key = temp_walker_list[walker_index].ball_key
                balls[ball_key][gv.num_cvs+2] = 0

            initial_weights_array = np.array(initial_weights)
            walker_indices = np.argsort(-initial_weights_array)
            temp_initial_indices = initial_indices
            # sorted indices based on descending order of weights
            initial_indices = [temp_initial_indices[i] for i in walker_indices]

            num_bins = 1
            true_num_bins = 1
            bins = [0]
            num_walkers_bin = [len(initial_indices)]

            if gv.rate_flag == 1:
                num_bins = gv.num_states
                true_num_bins = 0
                bins = []
                num_walkers_bin = []
                for i in range(gv.num_states):
                    num_walkers = 0
                    for j in initial_indices:
                        state = temp_walker_list[j].state
                        if state == i:
                            num_walkers += 1
                    if num_walkers != 0:
                        true_num_bins += 1
                        bins.append(i)
                        num_walkers_bin.append(num_walkers)
                if true_num_bins <= 1:
                    num_bins = 1
                    true_num_bins = 1
                    bins = [0]
                    num_walkers_bin = [len(initial_indices)]

            target_num_walkers = int(np.floor(float(initial_target_num_walkers)/true_num_bins))
            remainder = initial_target_num_walkers-target_num_walkers*true_num_bins

            for b, bin_index in enumerate(bins):
                new_weights = []
                new_indices = []
                new_num_walkers = 0
                # add the remaining walkers to the very last bin if there are any
                if remainder != 0 and b == (true_num_bins - 1):
                    target_num_walkers += remainder

                weights_bin = [float] * num_walkers_bin[b]
                indices_bin = [int] * num_walkers_bin[b]

                if num_bins == 1:
                    weights_bin = initial_weights
                    indices_bin = initial_indices
                else:
                    k = 0
                    for j in initial_indices:
                        state = temp_walker_list[j].state
                        if bin_index == state:
                            weights_bin[k] = temp_walker_list[j].weight
                            indices_bin[k] = temp_walker_list[j].global_index
                            k += 1

                total_weight = np.sum(weights_bin)
                target_weight = total_weight / target_num_walkers

                x = indices_bin.pop()
                while True:
                    x_weight = weights[x]
                    if x_weight >= target_weight or len(indices_bin) == 0:
                        r = max(1, int(np.floor(x_weight / target_weight)))
                        r = min(r, target_num_walkers - new_num_walkers)
                        new_num_walkers += r
                        for item in np.repeat(x, r):
                            new_indices.append(item)
                            new_weights.append(target_weight)
                        if new_num_walkers < target_num_walkers and x_weight - r * target_weight > 0.0:
                            indices_bin.append(x)
                            weights[x] = x_weight - r * target_weight
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

                for ni, global_index in enumerate(new_indices):
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[ni]
                        ball_key = walker_list[global_index].ball_key
                        if balls[ball_key][gv.num_cvs+2] == 0:
                            num_occupied_balls += 1
                        balls[ball_key][gv.num_cvs+2] += 1
                        ball_center = walker_list[global_index].current_ball_center
                        ball_to_walkers[tuple(ball_center)].append(global_index)
                        directory = gv.main_directory + '/CAS/walker' + str(global_index)
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
                        ball_key = walker_list[new_index].ball_key
                        if balls[ball_key][gv.num_cvs+2] == 0:
                            num_occupied_balls += 1
                        balls[ball_key][gv.num_cvs+2] += 1
                        ball_center = walker_list[new_index].current_ball_center
                        ball_to_walkers[tuple(ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        new_directory = gv.main_directory + '/CAS/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        os.system('sed -i \'$ d\' weight_trajectory.txt')
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % walker_list[new_index].weight + '\n')
                        f.close()

    total_num_walkers = num_occupied_big_clusters*gv.num_walkers_for_sc + num_occupied_small_clusters*gv.num_walkers
    if excess_index - total_num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

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

    while len(vacant_walker_indices) > 0:
        vacant_walker_indices.pop()
    gv.num_occupied_balls = num_occupied_balls
    gv.total_num_walkers = total_num_walkers
    gv.num_occupied_big_clusters = num_occupied_big_clusters
    gv.num_occupied_small_clusters = num_occupied_small_clusters


def resampling(walker_list, temp_walker_list, balls, ball_to_walkers, vacant_walker_indices):
    gv.sc_performed = 0
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.total_num_walkers)]
    occupied_indices = np.zeros(gv.num_balls_limit*gv.num_walkers*2, int)
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

            if gv.rate_flag == 1:
                num_bins = gv.num_states
                true_num_bins = 0
                bins = []
                num_walkers_bin = []
                for i in range(gv.num_states):
                    num_walkers = 0
                    for j in initial_indices:
                        state = temp_walker_list[j].state
                        if state == i:
                            num_walkers += 1
                    if num_walkers != 0:
                        true_num_bins += 1
                        bins.append(i)
                        num_walkers_bin.append(num_walkers)
                if true_num_bins <= 1:
                    num_bins = 1
                    true_num_bins = 1
                    bins = [0]
                    num_walkers_bin = [len(initial_indices)]

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

                if num_bins == 1:
                    weights_bin = initial_weights
                    indices_bin = initial_indices
                else:
                    k = 0
                    for j in initial_indices:
                        state = temp_walker_list[j].state
                        if bin_index == state:
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
                            os.chdir(gv.main_directory + '/CAS')
                            os.system('rm -rf walker' + str(y))

                if b == 0:  # reset balls
                    balls[current_ball][gv.num_cvs+2] = 0
                for ni, global_index in enumerate(new_indices):
                    if occupied_indices[global_index] == 0:
                        occupied_indices[global_index] = 1
                        walker_list[global_index].copy_walker(temp_walker_list[global_index])
                        walker_list[global_index].weight = new_weights[ni]
                        ball_to_walkers[tuple(current_ball_center)].append(global_index)
                        directory = gv.main_directory + '/CAS/walker' + str(global_index)
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
                        old_directory = gv.main_directory + '/CAS/walker' + str(global_index)
                        new_directory = gv.main_directory + '/CAS/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        os.system('sed -i \'$ d\' weight_trajectory.txt')
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % walker_list[new_index].weight + '\n')
                        f.close()
                    balls[current_ball][gv.num_cvs+2] += 1

    total_num_walkers = num_occupied_balls*gv.num_walkers
    if excess_index-total_num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

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

    while len(vacant_walker_indices) > 0:
        vacant_walker_indices.pop()
    gv.num_occupied_balls = num_occupied_balls
    gv.total_num_walkers = gv.num_occupied_balls*gv.num_walkers


def print_status(step_num, walker_list, balls, ball_to_walkers, ball_clusters_list, key_to_ball):
    os.chdir(gv.main_directory + '/CAS')
    total_weight = 0.0
    f = open('total_weight_on_each_ball_' + str(step_num + 1) + '.txt', 'w')
    for current_ball in range(balls.shape[0]):
        ball_center = balls[current_ball][0:gv.num_cvs].tolist()
        weights = [walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
        if np.sum(weights) > 0.0:
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
    if gv.enhanced_sampling_flag == 2:
        f.write(str(step_num + 1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + ' '
                + str(gv.num_occupied_big_clusters) + ' ' + str(gv.num_occupied_small_clusters) + '\n')
        gv.num_occupied_big_clusters = 0
        gv.num_occupied_small_clusters = 0
    else:
        f.write(str(step_num + 1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + '\n')
