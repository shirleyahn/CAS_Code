import numpy as np
import walker
import we_global_variables as gv
import we_energy_function as ef
from scipy import special
import os
import shutil
import copy


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += (values[i]-center[i])**2
    return np.sqrt(distance)


def set_parameters(input_parameter_file):
    with open(input_parameter_file, 'r') as f:
        gv.main_directory = f.readline().strip()
        f.readline()
        gv.flag = int(f.readline())
        gv.balls_flag = int(f.readline())
        gv.resample_less_flag = int(f.readline())
        gv.less_or_greater_flag = int(f.readline())
        gv.static_threshold_flag = int(f.readline())
        f.readline()
        gv.threshold = float(f.readline())
        gv.property_to_keep_track = int(f.readline())
        gv.enhanced_sampling_flag = int(f.readline())
        gv.num_balls_limit = int(f.readline())
        f.readline()
        gv.radius = float(f.readline())
        gv.num_walkers = int(f.readline())
        gv.num_cvs = int(f.readline())
        gv.grid_dimensions = f.readline().strip().split()
        f.readline()
        gv.m_steps_per_step = int(f.readline())
        gv.step_size = float(f.readline())
        gv.beta = float(f.readline())
        f.readline()
        gv.initial_step_num = int(f.readline())
        gv.max_num_steps = int(f.readline())
        gv.num_occupied_balls = int(f.readline())
        gv.first_walker = int(f.readline())
        gv.last_walker = int(f.readline())

    ball_volume = (np.pi**(gv.num_cvs/2)*gv.radius**gv.num_cvs)/special.gamma((gv.num_cvs/2)+1)
    grid_volume = 1.0
    ii = 0
    for i in range(gv.num_cvs):
        grid_volume *= (float(gv.grid_dimensions[ii+1])-float(gv.grid_dimensions[ii]))
        ii += 2
    gv.max_num_balls = 0
    if ball_volume != 0.0:
        gv.max_num_balls = int(np.floor(grid_volume/ball_volume))
    if gv.max_num_balls > gv.num_balls_limit or gv.max_num_balls < gv.num_balls_limit*1e-2:
        gv.max_num_balls = gv.num_balls_limit
    print 'max # of balls (n_b) = ' + str(gv.max_num_balls)
    gv.current_num_balls = 0


def initialize(input_initial_values_file, walker_list, temp_walker_list, ball_to_walkers, vacant_walker_indices):
    for i in range(len(walker_list)):
        walker_list[i] = walker.Walker([-1000.0] * gv.num_cvs, i, [-1000.0] * gv.num_cvs, 0.0, 0, 0.0)

    if gv.flag == 0:  # new simulation
        initial_weight = 1.0/(gv.num_walkers*gv.num_occupied_balls)
        f = open(input_initial_values_file, 'r')
        for n in range(gv.num_occupied_balls):
            initial_values = [None] * gv.num_cvs
            for i in range(gv.num_cvs):
                initial_values[i] = float(f.readline())
            for i in range(n*gv.num_walkers, (n+1)*gv.num_walkers):
                walker_list[i].set(initial_values, initial_weight)
        f.close()

        os.system('mkdir WE')
        os.chdir(gv.main_directory+'/WE')
        for i in range(gv.num_walkers*gv.num_occupied_balls):
            os.system('mkdir walker' + str(i))

    elif gv.flag == 1:  # restarting simulation in the middle of simulation
        for i in range(gv.num_walkers*gv.num_occupied_balls):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()

    elif gv.flag == 2:  # restarting simulation in the middle of binning
        for i in range(gv.num_walkers*gv.num_occupied_balls):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            os.chdir(walker_directory)
            f = open('weight_trajectory.txt', 'r')
            weight = float(f.readlines()[-1].strip())
            walker_list[i].weight = weight
            f.close()
            num_lines = sum(1 for line in open('ball_trajectory.txt'))
            # if walker is already binned to a ball, delete the binning and have binning start from scratch
            if num_lines > gv.initial_step_num:
                os.system('sed -i \'$d\' ball_trajectory.txt')

    elif gv.flag == 3:  # restarting simulation in the middle of resampling
        total_weight = 0.0
        previous_ball_to_walkers = {}
        previous_balls_order = {}
        previous_balls_weights = np.loadtxt('total_weight_of_each_ball_' + str(gv.initial_step_num) + '.txt')
        previous_balls_walker_count = np.zeros((previous_balls_weights.shape[0], previous_balls_weights.shape[1]))
        for i in range(previous_balls_weights.shape[0]):
            previous_ball_center = previous_balls_weights[i, 0:-1].tolist()
            previous_balls_order[tuple(previous_ball_center)] = i
            previous_balls_walker_count[i] = previous_balls_weights[i]
            previous_balls_walker_count[i, -1] = gv.num_walkers

        # TODO: make sure that gv.num_occupied_balls is equal to the highest walker number inside the WE folder
        for i in range(gv.num_occupied_balls+1):
            walker_directory = gv.main_directory + '/WE/walker' + str(i)
            # if all of the files exist in the walker folder, it is a complete walker
            if os.path.isfile(walker_directory + '/weight_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/ball_trajectory.txt') and \
                    os.path.isfile(walker_directory + '/trajectory.txt'):
                os.chdir(walker_directory)
                f = open('weight_trajectory.txt', 'r')
                weight = float(f.readlines()[-1].strip())
                walker_list[i].weight = weight
                total_weight += weight
                f.close()

                ball_trajectory = np.loadtxt('ball_trajectory.txt')
                previous_ball_center = ball_trajectory[-2].tolist()
                previous_ball_key = previous_balls_order[tuple(previous_ball_center)]
                previous_balls_weights[previous_ball_key, -1] -= weight
                if previous_balls_weights[previous_ball_key, -1] < 0.0:
                    print 'ERROR: weight is ' + str(previous_balls_weights[previous_ball_key, -1]) + ' for walker ' + \
                          str(i) + ' with ball_key ' + str(previous_ball_key)
                previous_balls_walker_count[previous_ball_key, -1] -= 1
                if previous_balls_walker_count[previous_ball_key, -1] < 0:
                    print 'ERROR: walker count is ' + str(previous_balls_walker_count[previous_ball_key, -1]) + \
                          ' for walker ' + str(i) + ' with ball key ' + str(previous_ball_key)
                if tuple(previous_ball_center) in previous_ball_to_walkers:
                    previous_ball_to_walkers[tuple(previous_ball_center)].append(i)
                else:
                    previous_ball_to_walkers[tuple(previous_ball_center)] = [i]
                current_ball_center = ball_trajectory[-1].tolist()
                if tuple(current_ball_center) in ball_to_walkers:
                    ball_to_walkers[tuple(current_ball_center)].append(i)
                else:
                    ball_to_walkers[tuple(current_ball_center)] = [i]
                    gv.current_num_balls += 1
                walker_list[i].ball_center = current_ball_center

                f = open('trajectory.txt', 'r')
                current_coordinates = f.readlines()[-1].strip().split()
                current_coordinates = [float(coordinate) for coordinate in current_coordinates]
                walker_list[i].coordinates = current_coordinates
                distance_from_center = calculate_distance_from_center(current_ball_center, current_coordinates)
                walker_list[i].distance_from_center = distance_from_center

                temp_walker_list[i] = walker.Walker(current_coordinates, i, current_ball_center, distance_from_center,
                                                    0, weight)

            # otherwise, it is an incomplete walker that needs missing files
            else:
                if os.path.isdir(walker_directory):
                    os.chdir(gv.main_directory + '/WE')
                    os.system('rm -rf walker' + str(i))
                vacant_walker_indices.append(i)

        # create new walkers for the remaining weights
        walker_index = 0
        excess_index = gv.num_occupied_balls+1
        for i in range(previous_balls_weights.shape[0]):
            if previous_balls_weights[i, -1] > 0.0:
                if previous_balls_walker_count[i, -1] <= 0:
                    print 'ERROR: at least one walker should exist if there is a weight of ' + \
                          str(previous_balls_weights[i, -1]) + ' for walker ' + str(i)
                else:
                    ball_center = previous_balls_weights[i, 0:-1].tolist()
                    reference_walker = ball_to_walkers[tuple(ball_center)][0]
                    reference_walker_directory = gv.main_directory + '/WE/walker/' + str(reference_walker)
                    if len(vacant_walker_indices) > 0:
                        walker_index = vacant_walker_indices.pop(0)
                    else:
                        walker_index = excess_index
                        excess_index += 1
                    walker_directory = gv.main_directory + '/WE/walker' + str(walker_index)
                    shutil.copytree(reference_walker_directory, walker_directory)

                    weight = previous_balls_weights[i, -1]
                    previous_balls_weights[i, -1] -= weight

                    os.chdir(walker_directory)
                    f = open('weight_trajectory.txt', 'w')
                    f.write(str(weight) + '\n')
                    walker_list[walker_index].weight = weight
                    total_weight += weight
                    f.close()

                    ball_to_walkers[tuple(ball_center)].append(walker_index)
                    walker_list[walker_index].ball_center = ball_center

                    f = open('trajectory.txt', 'r')
                    current_coordinates = f.readlines()[-1].strip().split()
                    current_coordinates = [float(coordinate) for coordinate in current_coordinates]
                    walker_list[walker_index].coordinates = current_coordinates
                    distance_from_center = calculate_distance_from_center(ball_center, current_coordinates)
                    walker_list[walker_index].distance_from_center = distance_from_center

                    temp_walker_list[walker_index] = walker.Walker(current_coordinates, walker_index, ball_center,
                                                                   distance_from_center, 0, weight)

        # check if total weight is 1.0
        if total_weight != 1.0:
            print 'ERROR: total weight is ' + str(total_weight)


def m_simulation(walker_list):
    for i in range(gv.num_occupied_balls*gv.num_walkers):
        temp_x = walker_list[i].coordinates[0]
        temp_y = walker_list[i].coordinates[1]
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
        walker_list[i].set([temp_x, temp_y])
        new_coordinates = walker_list[i].coordinates
        walker_directory = gv.main_directory + '/WE/walker' +str(i)
        os.chdir(walker_directory)
        f = open('trajectory.txt', 'a')
        f.write(' '.join(str(coordinate) for coordinate in new_coordinates))
        f.write('\n')
        f.close()


def binning(step_num, walker_list, temp_walker_list, balls, ball_to_walkers):
    initial_weights = [walker_list[i].weight for i in range(gv.num_occupied_balls*gv.num_walkers)]
    initial_weights_array = np.array(initial_weights)
    if gv.resample_less_flag == 1:
        walker_indices = np.argsort(initial_weights_array)  # sort walkers in ascending order based on their weights
    else:
        walker_indices = np.argsort(-initial_weights_array)  # sort walkers in descending order based on their weights

    start = 0  # indicates whether we are dealing with the very first walker or not
    new_threshold = gv.threshold  # only needed if gv.static_threshold_flag = 0
    for i in walker_indices:
        # first, go to walker directory i
        walker_directory = gv.main_directory + '/WE/walker' + str(i)
        os.chdir(walker_directory)

        # then, obtain new coordinates' values
        f = open('trajectory.txt', 'r')
        new_coordinates = f.readlines()[-1].strip().split()
        new_coordinates = [float(coordinate) for coordinate in new_coordinates]
        f.close()

        initial_step_num = walker_list[i].initial_step_num
        weight = walker_list[i].weight
        inside = 0  # indicates whether we are dealing with the very first walker or not
        # when we're  dealing with the very first walker, create the very first ball for the walker and we're done
        if start == 0:
            start += 1
            inside += 1
            ball_center = [coordinate for coordinate in new_coordinates]
            ball_to_walkers[tuple(ball_center)] = [i]
            temp_walker_list[i] = walker.Walker(new_coordinates, i, ball_center, 0.0, initial_step_num, weight)
            center_num_balls = copy.deepcopy(ball_center)
            center_num_balls.append(1)
            balls[gv.current_num_balls] = np.asarray(center_num_balls)
            gv.current_num_balls += 1

        distance = 0.0
        balls_key = 0
        # otherwise, loop through all of the balls and find the ball that has a center nearest the walker
        if inside == 0:
            for j in range(balls.shape[0]):
                ball_center = balls[j][:-1].tolist()
                distance_from_center = calculate_distance_from_center(ball_center, new_coordinates)
                if distance_from_center <= gv.radius:
                    inside += 1
                if distance == 0.0:
                    distance = distance_from_center
                    balls_key = j
                else:
                    if distance_from_center < distance:
                        distance = distance_from_center
                        balls_key = j

            if gv.property_to_keep_track < 0:
                property_to_keep_track = weight
            else:
                property_to_keep_track = new_coordinates[gv.property_to_keep_track]
            # walker is inside some ball or if we need to resample less and the weight is less or greater than threshold
            if inside != 0 or (gv.resample_less_flag == 1 and gv.less_or_greater_flag == 0 and property_to_keep_track < gv.threshold) \
                    or (gv.resample_less_flag == 1 and gv.less_or_greater_flag == 1 and property_to_keep_track > gv.threshold):
                balls[balls_key, gv.num_cvs] += 1
                ball_center = balls[balls_key][:-1].tolist()
                distance_from_center = calculate_distance_from_center(ball_center, new_coordinates)
                temp_walker_list[i] = walker.Walker(new_coordinates, i, ball_center, distance_from_center,
                                                    initial_step_num, weight)
                if tuple(ball_center) in ball_to_walkers:
                    ball_to_walkers[tuple(ball_center)].append(i)
                else:
                    ball_to_walkers[tuple(ball_center)] = [i]
            # walker is not inside any existing ball, so create a new ball
            else:
                ball_center = [coordinate for coordinate in new_coordinates]
                ball_to_walkers[tuple(ball_center)] = [i]
                temp_walker_list[i] = walker.Walker(new_coordinates, i, ball_center, 0.0, initial_step_num, weight)
                center_num_balls = copy.deepcopy(ball_center)
                center_num_balls.append(1)
                balls = np.append(balls, [np.asarray(center_num_balls)], axis=0)
                gv.current_num_balls += 1
            if gv.static_threshold_flag == 0 and gv.less_or_greater_flag == 0 and property_to_keep_track > new_threshold:
                new_threshold = property_to_keep_track
            elif gv.static_threshold_flag == 0 and gv.less_or_greater_flag == 1 and property_to_keep_track < new_threshold:
                new_threshold = property_to_keep_track

        # finally, write the new ball on the trajectory file
        ball_center = temp_walker_list[i].ball_center
        f = open('ball_trajectory.txt', 'a')
        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_center)))
        f.write('\n')
        f.close()

    os.chdir(gv.main_directory + '/WE')
    np.savetxt('balls_' + str(step_num+1) + '.txt', balls, fmt=' %+1.5f')
    if gv.static_threshold_flag == 0:
        gv.threshold = new_threshold
    return balls


def resampling(walker_list, temp_walker_list, balls, ball_to_walkers, vacant_walker_indices):
    num_occupied_balls = 0
    weights = [walker_list[i].weight for i in range(gv.num_occupied_balls*gv.num_walkers)]
    occupied_indices = np.zeros(gv.max_num_balls*gv.num_walkers, int)
    excess_index = gv.num_occupied_balls*gv.num_walkers
    for current_ball in range(balls.shape[0]):
        if int(balls[current_ball][gv.num_cvs]) > 0:
            num_occupied_balls += 1
            ball_center = balls[current_ball][:-1].tolist()
            initial_weights = [temp_walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
            initial_weights_array = np.array(initial_weights)
            walker_indices = np.argsort(-initial_weights_array)
            initial_indices = [temp_walker_list[i].global_index for i in ball_to_walkers[tuple(ball_center)]]
            temp_initial_indices = initial_indices
            # sorted indices based on descending order of weights
            initial_indices = [temp_initial_indices[i] for i in walker_indices]

            num_bins = 1
            true_num_bins = 1
            bins = [0]
            num_walkers_bin = [len(initial_indices)]

            if gv.enhanced_sampling_flag == 1:
                distance_from_center_list = [temp_walker_list[i].distance_from_center for i in
                                             ball_to_walkers[tuple(ball_center)]]
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
            ball_to_walkers[tuple(ball_center)] = []

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
                            # remove walker y directory
                            os.chdir(gv.main_directory + '/WE')
                            os.system('rm -rf walker' + str(y))

                if b == 0:  # reset balls
                    balls[current_ball][gv.num_cvs] = 0
                for ni, global_index in enumerate(new_indices):
                    coordinates = temp_walker_list[global_index].coordinates
                    if occupied_indices[global_index] == 0:
                        walker_list[global_index].set(coordinates, new_weights[ni])
                        walker_list[global_index].ball_center = ball_center
                        walker_list[global_index].distance_from_center = \
                            calculate_distance_from_center(ball_center, coordinates)
                        occupied_indices[global_index] = 1
                        ball_to_walkers[tuple(ball_center)].append(global_index)
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
                        ball_to_walkers[tuple(ball_center)].append(new_index)
                        old_directory = gv.main_directory + '/WE/walker' + str(global_index)
                        new_directory = gv.main_directory + '/WE/walker' + str(new_index)
                        shutil.copytree(old_directory, new_directory)
                        os.chdir(new_directory)
                        # write new weights on the trajectory file
                        f = open('weight_trajectory.txt', 'a')
                        f.write('% 1.20e' % walker_list[new_index].weight + '\n')
                        f.close()
                    balls[current_ball][gv.num_cvs] += 1

    if excess_index-num_occupied_balls*gv.num_walkers != len(vacant_walker_indices):
        print 'Something wrong with resampling'

    if num_occupied_balls >= gv.num_occupied_balls:
        for i in range(num_occupied_balls*gv.num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
    else:
        for i in range(gv.num_occupied_balls*gv.num_walkers, excess_index):
            new_index = vacant_walker_indices.pop()
            occupied_indices[new_index] = 1
            walker_list[new_index].copy_walker(walker_list[i])
            # rename the directory with name 'i' to 'new_index'
            os.chdir(gv.main_directory + '/WE')
            os.system('mv walker' + str(i) + ' walker' + str(new_index))
        for i in range(num_occupied_balls*gv.num_walkers, gv.num_occupied_balls*gv.num_walkers):
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


def print_status(step_num, walker_list, balls, ball_to_walkers):
    os.chdir(gv.main_directory + '/WE')
    total_weight = 0.0
    f = open('total_weight_on_each_ball_' + str(step_num+1) + '.txt', 'w')
    for current_ball in range(balls.shape[0]):
        ball_center = balls[current_ball][:-1].tolist()
        weights = [walker_list[i].weight for i in ball_to_walkers[tuple(ball_center)]]
        total_weight += np.sum(weights)
        ball_center.append(np.sum(weights))
        f.write(' '.join(map(lambda coordinate: str(coordinate), ball_center)))
        f.write('\n')
        # reset walkers and number of walkers that belong in each ball
        balls[current_ball][gv.num_cvs] = 0
        ball_to_walkers[tuple(ball_center[:-1])] = []
    f.close()

    # verify that total weight of all balls is 1.0
    f = open('total_weight.txt', 'a')
    f.write(str(step_num+1) + ' ' + str(total_weight) + ' ' + str(gv.num_occupied_balls) + '\n')