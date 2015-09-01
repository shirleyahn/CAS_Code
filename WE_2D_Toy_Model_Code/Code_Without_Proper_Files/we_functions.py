import numpy as np
import walker
import we_global_variables as gv
import we_energy_function as ef
from scipy import special
import os
import copy
import we_check_state_function


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += (values[i]-center[i])**2
    if abs(distance) < 1.0e-12:
        distance = 0.0
    return np.sqrt(distance)


def set_parameters(input_parameter_file):
    with open(input_parameter_file, 'r') as f:
        gv.main_directory = f.readline().strip()
        f.readline()
        gv.balls_flag = int(f.readline())
        gv.resample_less_flag = int(f.readline())
        gv.rate_flag = int(f.readline())
        gv.num_states = int(f.readline())
        f.readline()
        gv.less_or_greater_flag = int(f.readline())
        gv.static_threshold_flag = int(f.readline())
        temp_list = f.readline().strip().split()
        gv.threshold_values = [float(i) for i in temp_list]
        temp_list = f.readline().strip().split()
        gv.properties_to_keep_track = [int(i) for i in temp_list]
        f.readline()
        gv.enhanced_sampling_flag = int(f.readline())
        gv.num_balls_limit = int(f.readline())
        gv.radius = float(f.readline())
        gv.num_walkers = int(f.readline())
        gv.num_cvs = int(f.readline())
        temp_list = f.readline().strip().split()
        gv.grid_dimensions = [float(i) for i in temp_list]
        f.readline()
        gv.m_steps_per_step = int(f.readline())
        gv.step_size = float(f.readline())
        gv.beta = float(f.readline())
        f.readline()
        gv.max_num_steps = int(f.readline())
        gv.num_occupied_balls = int(f.readline())

    ball_volume = (np.pi**(gv.num_cvs/2)*gv.radius**gv.num_cvs)/special.gamma((gv.num_cvs/2)+1)
    grid_volume = 1.0
    ii = 0
    for i in range(gv.num_cvs):
        grid_volume *= (gv.grid_dimensions[ii+1]-gv.grid_dimensions[ii])
        ii += 2
    gv.max_num_balls = 0
    if ball_volume != 0.0:
        gv.max_num_balls = int(np.floor(grid_volume/ball_volume))
    if gv.max_num_balls > gv.num_balls_limit or gv.max_num_balls < gv.num_balls_limit*1e-2:
        gv.max_num_balls = gv.num_balls_limit
    print 'max # of balls (n_b) = ' + str(gv.max_num_balls)
    gv.current_num_balls = 0


def initialize(input_initial_values_file, walker_list):
    for i in range(len(walker_list)):
        walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, [-1000.0] * gv.num_cvs, i, [-1000.0] * gv.num_cvs,
                                       [-1000.0] * gv.num_cvs, 0, 0.0, 0, 0.0, -1)

    initial_weight = 1.0 / (gv.num_walkers * gv.num_occupied_balls)
    f = open(input_initial_values_file, 'r')
    for n in range(gv.num_occupied_balls):
        initial_values = [None] * gv.num_cvs
        for i in range(gv.num_cvs):
            initial_values[i] = float(f.readline())
        if gv.rate_flag == 1:
            initial_state = we_check_state_function.check_state_function(initial_values)
        for i in range(n * gv.num_walkers, (n + 1) * gv.num_walkers):
            walker_list[i].set(initial_values, initial_weight)
            if gv.rate_flag == 1:
                walker_list[i].state = initial_state
    f.close()

    os.system('mkdir WE')


def m_simulation(walker_list):
    for i in range(gv.num_occupied_balls*gv.num_walkers):
        previous_coordinates = walker_list[i].current_coordinates
        temp_x = walker_list[i].current_coordinates[0]
        temp_y = walker_list[i].current_coordinates[1]
        for j in range(gv.m_steps_per_step):
            direction = np.random.randint(0, 4)
            if direction == 0:  # move to left
                new_x = temp_x - gv.step_size
                new_y = temp_y
            elif direction == 1:  # move to right
                new_x = temp_x + gv.step_size
                new_y = temp_y
            elif direction == 2:  # move to top
                new_x = temp_x
                new_y = temp_y + gv.step_size
            else:  # move to bottom
                new_x = temp_x
                new_y = temp_y - gv.step_size
            old_energy = ef.energy_function(temp_x, temp_y)
            new_energy = ef.energy_function(new_x, new_y)
            if new_energy - old_energy <= 0.0:  # accept move
                temp_x = new_x
                temp_y = new_y
            else:
                random_number = np.random.uniform(0.0, 1.0)
                if random_number < np.exp(-(new_energy-old_energy)*gv.beta):  # accept move
                    temp_x = new_x
                    temp_y = new_y
        if abs(temp_x) < 1.0e-12:
            temp_x = 0.0
        if abs(temp_y) < 1.0e-12:
            temp_y = 0.0
        walker_list[i].set([temp_x, temp_y])
        walker_list[i].previous_coordinates = previous_coordinates


def binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers):
    initial_weights = [walker_list[i].weight for i in range(gv.num_occupied_balls*gv.num_walkers)]
    initial_weights_array = np.array(initial_weights)
    flux = np.zeros((gv.num_states, gv.num_states))
    if gv.resample_less_flag == 1:
        walker_indices = np.argsort(initial_weights_array)  # sort walkers in ascending order based on their weights
    else:
        walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    start = 0  # indicates whether we are dealing with the very first walker or not
    new_threshold_values = gv.threshold_values  # only needed if gv.static_threshold_flag = 0
    for i in walker_indices:
        new_coordinates = walker_list[i].current_coordinates
        if gv.rate_flag == 1:
            state = we_check_state_function.check_state_function(new_coordinates)
            if walker_list[i].state != -1 and state != -1 and walker_list[i].state != state:
                flux[walker_list[i].state, state] += walker_list[i].weight
            elif walker_list[i].state != -1 and state == -1:
                state = walker_list[i].state
        else:
            state = -1
        previous_coordinates = walker_list[i].previous_coordinates
        previous_ball_center = walker_list[i].current_ball_center
        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight
        inside = 0  # indicates whether we are dealing with the very first walker or not
        # when we're  dealing with the very first walker, create the very first ball for the walker and we're done
        if start == 0:
            start += 1
            inside += 1
            current_ball_center = [coordinate for coordinate in new_coordinates]
            ball_to_walkers[tuple(current_ball_center)] = [i]
            temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                current_ball_center, gv.current_num_balls, 0.0, initial_step_num,
                                                weight, state)
            center_key_num = copy.deepcopy(current_ball_center)
            center_key_num.append(gv.current_num_balls)
            center_key_num.append(1)
            balls[gv.current_num_balls] = np.asarray(center_key_num)
            gv.current_num_balls += 1

        distance = 0.0
        ball_key = 0
        # otherwise, loop through all of the balls and find the ball that has a center nearest the walker
        if inside == 0:
            for j in range(balls.shape[0]):
                current_ball_center = balls[j][:-2].tolist()
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                if distance_from_center <= gv.radius:
                    inside += 1
                if distance == 0.0:
                    distance = distance_from_center
                    ball_key = j
                else:
                    if distance_from_center < distance:
                        distance = distance_from_center
                        ball_key = j

            properties_to_keep_track = []
            for k in range(len(gv.properties_to_keep_track)):
                if gv.properties_to_keep_track[k] < 0:
                    properties_to_keep_track.append(weight)
                else:
                    properties_to_keep_track.append(new_coordinates[gv.properties_to_keep_track[k]])
            bin_walker = 0
            if gv.resample_less_flag == 1 and gv.less_or_greater_flag == 0:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] < gv.threshold_values[m]:
                        bin_walker += 1
            elif gv.resample_less_flag == 1 and gv.less_or_greater_flag == 1:
                for m in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[m] > gv.threshold_values[m]:
                        bin_walker += 1
            # walker is inside some ball or if we need to resample less and the properties to keep track of are all less
            # or greater than the threshold values
            if inside != 0 or (gv.resample_less_flag == 1 and gv.less_or_greater_flag == 0 and bin_walker != 0) \
                    or (gv.resample_less_flag == 1 and gv.less_or_greater_flag == 1 and bin_walker != 0):
                balls[ball_key][gv.num_cvs+1] += 1
                current_ball_center = balls[ball_key][:-2].tolist()
                distance_from_center = calculate_distance_from_center(current_ball_center, new_coordinates)
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                    current_ball_center, balls[ball_key][-2], distance_from_center,
                                                    initial_step_num, weight, state)
                if tuple(current_ball_center) in ball_to_walkers:
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                else:
                    ball_to_walkers[tuple(current_ball_center)] = [i]
            # walker is not inside any existing ball, so create a new ball
            else:
                current_ball_center = [coordinate for coordinate in new_coordinates]
                ball_to_walkers[tuple(current_ball_center)] = [i]
                temp_walker_list[i] = walker.Walker(previous_coordinates, new_coordinates, i, previous_ball_center,
                                                    current_ball_center, gv.current_num_balls, 0.0, initial_step_num,
                                                    weight, state)
                center_key_num = copy.deepcopy(current_ball_center)
                center_key_num.append(gv.current_num_balls)
                center_key_num.append(1)
                balls = np.append(balls, [np.asarray(center_key_num)], axis=0)
                gv.current_num_balls += 1
            # update threshold values if gv.static_threshold_flag = 0
            if gv.static_threshold_flag == 0 and gv.less_or_greater_flag == 0:
                for n in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[n] > new_threshold_values[n]:
                        new_threshold_values[n] = properties_to_keep_track[n]
            elif gv.static_threshold_flag == 0 and gv.less_or_greater_flag == 1:
                for n in range(len(gv.properties_to_keep_track)):
                    if properties_to_keep_track[n] < new_threshold_values[n]:
                        new_threshold_values[n] = properties_to_keep_track[n]

    os.chdir(gv.main_directory + '/WE')
    np.savetxt('balls_' + str(step_num+1) + '.txt', balls, fmt=' %+1.5f')
    if gv.rate_flag == 1:
        np.savetxt('flux_' + str(step_num+1) + '.txt', flux, fmt=' %1.5e')
    if gv.static_threshold_flag == 0:
        gv.threshold_values = new_threshold_values
    return balls


def calculating_transition(step_num,  temp_walker_list, balls):
    transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
    for i in range(gv.num_occupied_balls*gv.num_walkers):
        previous_coordinates = temp_walker_list[i].previous_ball_center
        inside = 0
        distance = 0.0
        ball_key = 0
        for j in range(balls.shape[0]):
            ball_center = balls[j][:-2].tolist()
            distance_from_center = calculate_distance_from_center(ball_center, previous_coordinates)
            if distance_from_center <= gv.radius:
                inside += 1
            if distance == 0.0:
                distance = distance_from_center
                ball_key = j
            else:
                if distance_from_center < distance:
                    distance = distance_from_center
                    ball_key = j
        if inside != 0:
            transition_matrix[ball_key][temp_walker_list[i].ball_key] += temp_walker_list[i].weight
    os.chdir(gv.main_directory + '/WE')
    np.savetxt('transition_matrix_' + str(step_num+1) + '.txt', transition_matrix, fmt=' %1.5e')


def resampling(walker_list, temp_walker_list, balls, ball_to_walkers):
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.num_occupied_balls*gv.num_walkers)]
    occupied_indices = np.zeros(gv.max_num_balls*gv.num_walkers, int)
    excess_index = gv.num_occupied_balls*gv.num_walkers
    vacant_walker_indices = []
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs+1]) > 0:
            num_occupied_balls += 1
            current_ball_center = balls[current_ball][:-2].tolist()
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
                distance_from_center_list = [temp_walker_list[i].distance_from_center for i in
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
                            distance = temp_walker_list[ii].distance_from_center
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

                weights_bin = [float]*num_walkers_bin[b]
                indices_bin = [int]*num_walkers_bin[b]

                if gv.enhanced_sampling_flag == 0 or (gv.enhanced_sampling_flag == 1 and std == 0.0):
                    weights_bin = initial_weights
                    indices_bin = initial_indices

                elif gv.enhanced_sampling_flag == 1 and std != 0.0:
                    k = 0
                    for j in initial_indices:
                        distance = temp_walker_list[j].distance_from_center
                        if bin_index == 0:
                            if distance <= std:
                                weights_bin[k] = temp_walker_list[j].weight
                                indices_bin[k] = temp_walker_list[j].global_index
                                k += 1
                        elif bin_index == num_bins-1:
                            if bin_index*std < distance:
                                weights_bin[k] = temp_walker_list[j].weight
                                indices_bin[k] = temp_walker_list[j].global_index
                                k += 1
                        else:
                            if bin_index*std < distance <= (bin_index+1)*std:
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

                if b == 0:  # reset balls
                    balls[current_ball][gv.num_cvs+1] = 0
                for ni, global_index in enumerate(new_indices):
                    previous_coordinates = temp_walker_list[global_index].previous_coordinates
                    current_coordinates = temp_walker_list[global_index].current_coordinates
                    previous_ball_center = temp_walker_list[global_index].previous_ball_center
                    ball_key = temp_walker_list[global_index].ball_key
                    if occupied_indices[global_index] == 0:
                        walker_list[global_index].set(current_coordinates, new_weights[ni])
                        walker_list[global_index].previous_coordinates = previous_coordinates
                        walker_list[global_index].previous_ball_center = previous_ball_center
                        walker_list[global_index].current_ball_center = current_ball_center
                        walker_list[global_index].distance_from_center = \
                            calculate_distance_from_center(current_ball_center, current_coordinates)
                        walker_list[global_index].ball_key = ball_key
                        occupied_indices[global_index] = 1
                        ball_to_walkers[tuple(current_ball_center)].append(global_index)
                    else:
                        if len(vacant_walker_indices) > 0:
                            new_index = vacant_walker_indices.pop()
                        else:
                            new_index = excess_index
                            excess_index += 1
                        occupied_indices[new_index] = 1
                        walker_list[new_index].copy_walker(walker_list[global_index])
                        ball_to_walkers[tuple(current_ball_center)].append(new_index)
                    balls[current_ball][gv.num_cvs] += 1

    if excess_index-num_occupied_balls*gv.num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    if num_occupied_balls >= gv.num_occupied_balls:
        for i in range(num_occupied_balls*gv.num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
    else:
        for i in range(gv.num_occupied_balls*gv.num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
        for i in range(num_occupied_balls*gv.num_walkers, gv.num_occupied_balls*gv.num_walkers):
            if occupied_indices[i] == 1:
                new_index = vacant_walker_indices.pop()
                while new_index >= num_occupied_balls*gv.num_walkers:
                    new_index = vacant_walker_indices.pop()
                occupied_indices[new_index] = 1
                walker_list[new_index].copy_walker(walker_list[i])

    gv.num_occupied_balls = num_occupied_balls


def print_status(step_num, walker_list, balls, ball_to_walkers):
    os.chdir(gv.main_directory + '/WE')
    total_weight = 0.0
    f = open('total_weight_on_each_ball_' + str(step_num+1) + '.txt', 'w')
    for current_ball in range(balls.shape[0]):
        ball_center = balls[current_ball][:-2].tolist()
        weights = [walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
        total_weight += np.sum(weights)
        ball_center.append(np.sum(weights))
        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_center)))
        f.write('\n')
        # reset walkers and number of walkers that belong in each ball
        balls[current_ball][gv.num_cvs+1] = 0
        ball_to_walkers[tuple(ball_center)] = []
    f.close()

    # verify that total weight of all balls is 1.0
    f = open('total_weight.txt', 'a')
    f.write(str(step_num+1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + '\n')
