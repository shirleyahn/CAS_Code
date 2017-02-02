import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError

num_states = 2
num_cvs = 6
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values_first')
#traj_file = traj_file[0:6000]
state_file = np.loadtxt('initial_states_first', dtype=int)
#state_file = state_file[0:6000]
radius = 80.0
num_clusters = 10
boundary = 1.3


def closest_ball(coordinates, balls, num_cvs):
    distance = np.zeros((balls.shape[0],))
    for i in range(num_cvs):
        distance += np.minimum(360.0 - np.abs(balls[:, i] - coordinates[i]), np.abs(balls[:, i] - coordinates[i])) ** 2
    return np.argmin(distance)


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += min(360.0 - abs(values[i] - center[i]), abs(values[i] - center[i])) ** 2
    return np.sqrt(distance)


def check_state_function(coordinates):
    if -100.0 <= coordinates[0] <= -30.0 and -90.0 <= coordinates[1] <= -10.0 and -100.0 <= coordinates[2] <= -30.0 and \
       -90.0 <= coordinates[3] <= -10.0 and -100.0 <= coordinates[4] <= -30.0 and -90.0 <= coordinates[5] <= -10.0:
        return 0
    elif -180.0 <= coordinates[0] <= -55.0 and (105.0 <= coordinates[1] <= 180.0 or -180.0 <= coordinates[1] <= -155.0) and \
         -180.0 <= coordinates[2] <= -55.0 and (105.0 <= coordinates[3] <= 180.0 or -180.0 <= coordinates[3] <= -155.0) and \
         -180.0 <= coordinates[4] <= -55.0 and (105.0 <= coordinates[5] <= 180.0 or -180.0 <= coordinates[5] <= -155.0):
        return 1
    else:
        return -1


def convert_angles_to_cos_sin(balls):
    new_balls = np.zeros((balls.shape[0], 2*balls.shape[1]))
    for i in range(balls.shape[0]):
        for j in range(balls.shape[1]):
            new_balls[i, 2*j] = np.cos(balls[i, j]*np.pi/180)
            new_balls[i, 2*j+1] = np.sin(balls[i, j]*np.pi/180)
    return new_balls


# first, separate out folded and unfolded states from the rest for partitioning
folded_balls = np.zeros((1, num_cvs))
folded_ball_to_traj = np.zeros((1,))
unfolded_balls = np.zeros((1, num_cvs))
unfolded_ball_to_traj = np.zeros((1,))
balls = np.zeros((1, num_cvs))
ball_to_traj = np.zeros((1,))
states = np.zeros((1,), dtype=int)
folded_present = 0
unfolded_present = 0
neither_present = 0
for i in range(traj_file.shape[0]):
    current_state = check_state_function(traj_file[i])
    if current_state == 0:
        if folded_present == 0:
            folded_balls[0] = traj_file[i]
            folded_ball_to_traj[0] = i*time_step
            folded_present += 1
        else:
            folded_balls = np.append(folded_balls, [np.asarray(traj_file[i])], axis=0)
            folded_ball_to_traj = np.append(folded_ball_to_traj, [np.asarray(i*time_step)], axis=0)
    elif current_state == 1:
        if unfolded_present == 0:
            unfolded_balls[0] = traj_file[i]
            unfolded_ball_to_traj[0] = i*time_step
            unfolded_present += 1
        else:
            unfolded_balls = np.append(unfolded_balls, [np.asarray(traj_file[i])], axis=0)
            unfolded_ball_to_traj = np.append(unfolded_ball_to_traj, [np.asarray(i*time_step)], axis=0)
    else:
        if neither_present == 0:
            balls[0] = traj_file[i]
            ball_to_traj[0] = i*time_step
            states[0] = state_file[i]
            neither_present += 1
        else:
            balls = np.append(balls, [np.asarray(traj_file[i])], axis=0)
            ball_to_traj = np.append(ball_to_traj, [np.asarray(i*time_step)], axis=0)
            states = np.append(states, [np.asarray(state_file[i])], axis=0)
unique, counts = np.unique(states, return_counts=True)
print folded_balls.shape, unfolded_balls.shape, balls.shape
print np.asarray((unique,counts)).T
folded_centroid, label = kmeans2(convert_angles_to_cos_sin(folded_balls), 1, minit='points', iter=200, missing='raise')
unfolded_centroid, label = kmeans2(convert_angles_to_cos_sin(unfolded_balls), 1, minit='points', iter=200, missing='raise')
folded_centroid = folded_centroid.reshape((num_cvs*2,))
unfolded_centroid = unfolded_centroid.reshape((num_cvs*2,))
closest_folded = closest_ball(folded_centroid, convert_angles_to_cos_sin(folded_balls), num_cvs*2)
folded_centroid = folded_balls[closest_folded]
closest_unfolded = closest_ball(unfolded_centroid, convert_angles_to_cos_sin(unfolded_balls), num_cvs*2)
unfolded_centroid = unfolded_balls[closest_unfolded]

# second, cover the folded states with macrostates
new_folded_balls = np.zeros((1, num_cvs))
new_folded_ball_to_traj = np.zeros((1,))
balls_count = np.zeros((1, 2), dtype=int)
balls_assignment = np.zeros((traj_file.shape[0],), dtype=int)
current_num_balls = 0
index = 0
for i in range(folded_balls.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        new_folded_balls[index] = folded_balls[i]
        new_folded_ball_to_traj[index] = folded_ball_to_traj[i]
        balls_count[current_num_balls, 0] += 1
        ball_index = int(folded_ball_to_traj[i]/time_step)
        balls_assignment[ball_index] = current_num_balls
        current_num_balls += 1
        index += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(folded_balls[i], new_folded_balls, num_cvs)
        current_ball_center = new_folded_balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, folded_balls[i])
        if distance_from_center <= radius:
            inside += 1
            balls_count[ball_key, 0] += 1
            ball_index = int(folded_ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            new_folded_balls = np.append(new_folded_balls, [np.asarray(folded_balls[i])], axis=0)
            new_folded_ball_to_traj = np.append(new_folded_ball_to_traj, [np.asarray(folded_ball_to_traj[i])], axis=0)
            balls_count = np.append(balls_count, [np.asarray((1, 0))], axis=0)
            ball_index = int(folded_ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = current_num_balls
            current_num_balls += 1
            index += 1

# third, loop through all of the walkers once more to assign them to their true nearest macrostates
for i in range(folded_balls.shape[0]):
    ball_index = int(folded_ball_to_traj[i]/time_step)
    old_ball_key = int(balls_assignment[ball_index])
    new_ball_key = closest_ball(folded_balls[i], new_folded_balls, num_cvs)
    balls_count[old_ball_key, 0] -= 1
    balls_count[new_ball_key, 0] += 1
    balls_assignment[ball_index] = new_ball_key

# fourth, delete empty macrostates
delete_list = []
index = 0
for i in range(new_folded_balls.shape[0]):
    if balls_count[i, 0] != 0:
        for j in range(traj_file.shape[0]):
            if balls_assignment[j] == i:
                balls_assignment[j] = index
        index += 1
    else:
        delete_list.append(i)
new_folded_balls = np.delete(new_folded_balls, delete_list, 0)
new_folded_ball_to_traj = np.delete(new_folded_ball_to_traj, delete_list, 0)
balls_count = np.delete(balls_count, delete_list, 0)
current_num_balls = new_folded_balls.shape[0]

# second, cover the unfolded states with macrostates
new_unfolded_balls = np.zeros((1, num_cvs))
new_unfolded_ball_to_traj = np.zeros((1,))
index = 0
for i in range(unfolded_balls.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        new_unfolded_balls[index] = unfolded_balls[i]
        new_unfolded_ball_to_traj[index] = unfolded_ball_to_traj[i]
        balls_count = np.append(balls_count, [np.asarray((0, 1))], axis=0)
        ball_index = int(unfolded_ball_to_traj[i]/time_step)
        balls_assignment[ball_index] = current_num_balls
        current_num_balls += 1
        index += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(unfolded_balls[i], new_unfolded_balls, num_cvs)
        current_ball_center = new_unfolded_balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, unfolded_balls[i])
        if distance_from_center <= radius:
            inside += 1
            ball_key += new_folded_balls.shape[0]
            balls_count[ball_key, 1] += 1
            ball_index = int(unfolded_ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            new_unfolded_balls = np.append(new_unfolded_balls, [np.asarray(unfolded_balls[i])], axis=0)
            new_unfolded_ball_to_traj = np.append(new_unfolded_ball_to_traj, [np.asarray(unfolded_ball_to_traj[i])], axis=0)
            balls_count = np.append(balls_count, [np.asarray((0, 1))], axis=0)
            ball_index = int(unfolded_ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = current_num_balls
            current_num_balls += 1
            index += 1

# third, loop through all of the walkers once more to assign them to their true nearest macrostates
for i in range(unfolded_balls.shape[0]):
    ball_index = int(unfolded_ball_to_traj[i]/time_step)
    old_ball_key = int(balls_assignment[ball_index])
    new_ball_key = closest_ball(unfolded_balls[i], new_unfolded_balls, num_cvs)+new_folded_balls.shape[0]
    balls_count[old_ball_key, 1] -= 1
    balls_count[new_ball_key, 1] += 1
    balls_assignment[ball_index] = new_ball_key

# fourth, delete empty macrostates
delete_list = []
index = new_folded_balls.shape[0]
for i in range(new_unfolded_balls.shape[0]):
    ball_index = i+new_folded_balls.shape[0]
    if balls_count[ball_index, 1] != 0:
        for j in range(traj_file.shape[0]):
            if balls_assignment[j] == ball_index:
                balls_assignment[j] = index
        index += 1
    elif balls_count[ball_index, 0] == 0 and balls_count[ball_index, 1] == 0:
        delete_list.append(i)
new_unfolded_balls = np.delete(new_unfolded_balls, delete_list, 0)
new_unfolded_ball_to_traj = np.delete(new_unfolded_ball_to_traj, delete_list, 0)
if len(delete_list) != 0:
    delete_list += new_folded_balls.shape[0]
    balls_count = np.delete(balls_count, delete_list, 0)

# second, cover the rest of the states with macrostates
new_balls = np.zeros((1, num_cvs))
new_ball_to_traj = np.zeros((1,))
new_states = np.zeros((1,), dtype=int)
new_distances = np.zeros((1, 2))
index = 0
for i in range(balls.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        new_balls[index] = balls[i]
        new_ball_to_traj[index] = ball_to_traj[i]
        new_states[index] = states[i]
        new_distances[index, 0] = calculate_distance_from_center(folded_centroid, balls[i])
        new_distances[index, 1] = calculate_distance_from_center(unfolded_centroid, balls[i])
        if states[i] == 0:
            balls_count = np.append(balls_count, [np.asarray((1, 0))], axis=0)
        else:
            balls_count = np.append(balls_count, [np.asarray((0, 1))], axis=0)
        ball_index = int(ball_to_traj[i]/time_step)
        balls_assignment[ball_index] = current_num_balls
        current_num_balls += 1
        index += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(balls[i], new_balls, num_cvs)
        current_ball_center = new_balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, balls[i])
        if distance_from_center <= radius:
            inside += 1
            new_distances[ball_key, 0] += calculate_distance_from_center(folded_centroid, balls[i])
            new_distances[ball_key, 1] += calculate_distance_from_center(unfolded_centroid, balls[i])
            ball_key += new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
            balls_count[ball_key, states[i]] += 1
            ball_index = int(ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            new_balls = np.append(new_balls, [np.asarray(balls[i])], axis=0)
            new_ball_to_traj = np.append(new_ball_to_traj, [np.asarray(ball_to_traj[i])], axis=0)
            new_states = np.append(new_states, [np.asarray(states[i])], axis=0)
            folded_distance = calculate_distance_from_center(folded_centroid, balls[i])
            unfolded_distance = calculate_distance_from_center(unfolded_centroid, balls[i])
            new_distances = np.append(new_distances, [np.asarray((folded_distance, unfolded_distance))], axis=0)
            if states[i] == 0:
                balls_count = np.append(balls_count, [np.asarray((1, 0))], axis=0)
            else:
                balls_count = np.append(balls_count, [np.asarray((0, 1))], axis=0)
            ball_index = int(ball_to_traj[i]/time_step)
            balls_assignment[ball_index] = current_num_balls
            current_num_balls += 1
            index += 1

# third, loop through all of the walkers once more to assign them to their true nearest macrostates
for i in range(balls.shape[0]):
    ball_index = int(ball_to_traj[i]/time_step)
    old_ball_key = int(balls_assignment[ball_index])
    new_ball_key = closest_ball(balls[i], new_balls, num_cvs)+new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
    balls_count[old_ball_key, states[i]] -= 1
    balls_count[new_ball_key, states[i]] += 1
    balls_assignment[ball_index] = new_ball_key

# fourth, delete empty macrostates
delete_list = []
index = new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
for i in range(new_balls.shape[0]):
    ball_index = i+new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
    if balls_count[ball_index, 0] != 0 or balls_count[ball_index, 1] != 0:
        for j in range(traj_file.shape[0]):
            if balls_assignment[j] == ball_index:
                balls_assignment[j] = index
        index += 1
    elif balls_count[ball_index, 0] == 0 and balls_count[ball_index, 1] == 0:
        delete_list.append(i)
new_balls = np.delete(new_balls, delete_list, 0)
new_ball_to_traj = np.delete(new_ball_to_traj, delete_list, 0)
new_states = np.delete(new_states, delete_list, 0)
new_distances = np.delete(new_distances, delete_list, 0)
if len(delete_list) != 0:
    delete_list += new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
    balls_count = np.delete(balls_count, delete_list, 0)
current_num_balls = new_folded_balls.shape[0]+new_unfolded_balls.shape[0]

for i in range(new_distances.shape[0]):
    ball_index = i+new_folded_balls.shape[0]+new_unfolded_balls.shape[0]
    if balls_count[ball_index, 0] != 0:
        new_distances[i, 0] /= balls_count[ball_index, 0]
        #new_distances[i, 0] = calculate_distance_from_center(folded_centroid, new_balls[i])
    if balls_count[ball_index, 1] != 0:
        new_distances[i, 1] /= balls_count[ball_index, 1]
        #new_distances[i, 1] = calculate_distance_from_center(unfolded_centroid, new_balls[i])

# fifth, calculate the transition matrix and its eigenvalues and eigenvectors
transition_matrix_size = new_folded_balls.shape[0]+new_unfolded_balls.shape[0]+new_balls.shape[0]
transition_matrix = np.zeros((transition_matrix_size, transition_matrix_size))
for i in range(traj_file.shape[0]-1):
    previous_ball_key = balls_assignment[i]
    current_ball_key = balls_assignment[i+1]
    transition_matrix[previous_ball_key][current_ball_key] += 1.0

new_transition_matrix = np.zeros((transition_matrix_size, transition_matrix_size))
for i in range(new_transition_matrix.shape[0]):
    for j in range(new_transition_matrix.shape[1]):
        new_transition_matrix[i][j] = (transition_matrix[i][j]+transition_matrix[j][i])/2.0

row_sum = np.sum(new_transition_matrix, axis=1)
for i in range(new_transition_matrix.shape[0]):
    if row_sum[i] != 0.0:
        new_transition_matrix[i, :] /= row_sum[i]

evalues, evectors = np.linalg.eig(new_transition_matrix.T)
idx = abs(evalues).argsort()[::-1]
final_evalues = evalues[idx]
final_evectors = evectors[:, idx]
np.savetxt('evalues.txt', final_evalues, fmt=' %1.10e')
np.savetxt('evectors.txt', final_evectors, fmt=' %1.10e')

# sixth, normalize the second evector by the first evector values -> good approximation to committor functions.
normalized_second_evector = np.zeros((final_evectors.shape[0], 1))
for i in range(final_evectors.shape[0]):
    if final_evectors[i, 0] != 0.0:
        normalized_second_evector[i] = final_evectors[i, 1] / abs(final_evectors[i, 0])
    else:
        normalized_second_evector[i] = 0.0

#if np.min(normalized_second_evector) != 0.0:
    #normalized_second_evector /= np.min(normalized_second_evector)  # normalize
np.savetxt('committor_function.txt', normalized_second_evector, fmt=' %1.10e')
#print normalized_second_evector[0], normalized_second_evector[1]

# seventh, cluster states with committor function
folded_centroid, label = kmeans2(normalized_second_evector[0:new_folded_balls.shape[0]], 1, minit='points', iter=200, missing='raise')
unfolded_centroid, label = kmeans2(normalized_second_evector[new_folded_balls.shape[0]:new_folded_balls.shape[0]+new_unfolded_balls.shape[0]], 1, minit='points', iter=200, missing='raise')
print folded_centroid, unfolded_centroid
closest_folded = closest_ball(folded_centroid, normalized_second_evector[0:new_folded_balls.shape[0]], 1)
folded_centroid = new_folded_balls[closest_folded]
closest_unfolded = closest_ball(unfolded_centroid, normalized_second_evector[new_folded_balls.shape[0]:new_folded_balls.shape[0]+new_unfolded_balls.shape[0]], 1)
unfolded_centroid = new_unfolded_balls[closest_unfolded]

neither_committor = normalized_second_evector[new_folded_balls.shape[0]+new_unfolded_balls.shape[0]:]
if num_clusters > 1:
    interval = boundary*2/num_clusters
    clusters = np.arange(-boundary,boundary,interval)
    print clusters

if num_clusters > 1:
    balls = np.zeros((num_states+new_balls.shape[0], num_cvs))
    ball_to_traj = np.zeros((num_states+new_balls.shape[0],))
    states = np.zeros((num_states+new_balls.shape[0],))
else:
    balls = np.zeros((num_states+1, num_cvs))
    ball_to_traj = np.zeros((num_states+1,))
    states = np.zeros((num_states+1,))

balls[0] = folded_balls[0]
ball_to_traj[0] = folded_ball_to_traj[0]
states[0] = 0

balls[1] = unfolded_balls[0]
ball_to_traj[1] = unfolded_ball_to_traj[0]
states[1] = 1

if num_clusters > 1:
    balls_file = np.zeros((num_states+new_balls.shape[0], num_cvs+3))
else:
    balls_file = np.zeros((num_states+1, num_cvs+3))
balls_file[0,0:num_cvs] = balls[0]
balls_file[1,0:num_cvs] = balls[1]
balls_file[0,num_cvs] = 0
balls_file[1,num_cvs] = 1
balls_file[0,num_cvs+1] = 0
balls_file[1,num_cvs+1] = 1
balls_file[:,num_cvs+2] = 0

if num_clusters == 1:
    balls[num_states] = new_balls[0]
    ball_to_traj[num_states] = new_ball_to_traj[0]
    states[num_states] = new_states[0]
    balls_file[num_states,0:num_cvs] = new_balls[0]
    balls_file[num_states,num_cvs] = num_states
    balls_file[num_states,num_cvs+1] = num_states
else:
    starting_index = num_states
    new_balls_count = balls_count[new_folded_balls.shape[0]+new_unfolded_balls.shape[0]:,:]
    for i in range(num_clusters):
        first = 0
        neither_balls = np.zeros((1, num_cvs))
        neither_ball_to_traj = np.zeros((1,))
        neither_states = np.zeros((1,))
        states_count = np.zeros((1, 2))
        distance = np.zeros((1, 2))
        for j in range(new_balls.shape[0]):
            if clusters[i] <= neither_committor[j] and ((i < num_clusters-1 and neither_committor[j] < clusters[i+1]) or (i == num_clusters-1 and neither_committor[j] <= boundary)) and first == 0:
                first += 1
                neither_balls[0] = new_balls[j]
                neither_ball_to_traj[0] = new_ball_to_traj[j]
                neither_states[0] = new_states[j]
                states_count += new_balls_count[j]
                folded_distance = calculate_distance_from_center(folded_centroid, new_balls[j])
                unfolded_distance = calculate_distance_from_center(unfolded_centroid, new_balls[j])
                distance += np.asarray((folded_distance, unfolded_distance))
                #print neither_committor[j], new_balls[j], new_balls_count[j], folded_distance, unfolded_distance
            elif clusters[i] <= neither_committor[j] and ((i < num_clusters-1 and neither_committor[j] < clusters[i+1]) or (i == num_clusters-1 and neither_committor[j] <= boundary)) and first != 0:
                neither_balls = np.append(neither_balls, [np.asarray(new_balls[j])], axis=0)
                neither_ball_to_traj = np.append(neither_ball_to_traj, [np.asarray(new_ball_to_traj[j])], axis=0)
                neither_states = np.append(neither_states, [np.asarray(new_states[j])], axis=0)
                states_count += new_balls_count[j]
                folded_distance = calculate_distance_from_center(folded_centroid, new_balls[j])
                unfolded_distance = calculate_distance_from_center(unfolded_centroid, new_balls[j])
                distance += np.asarray((folded_distance, unfolded_distance))
                #print neither_committor[j], new_balls[j], new_balls_count[j], folded_distance, unfolded_distance
        print clusters[i], states_count, neither_balls[0]
        #if states_count[0, 0] != 0:
            #print distance[0, 0]/states_count[0, 0]
        #else:
            #print distance[0, 0]
        #if states_count[0, 1] != 0:
            #print distance[0, 1]/states_count[0, 1]
        #else:
            #print distance[0, 1]
        for j in range(neither_balls.shape[0]):
            index = starting_index+j
            balls[index] = neither_balls[j]
            balls_file[index,0:num_cvs] = neither_balls[j]
            balls_file[index,num_cvs] = starting_index
            balls_file[index,num_cvs+1] = num_states+i
            ball_to_traj[index] = neither_ball_to_traj[j]
            states[index] = neither_states[j]
        starting_index += neither_balls.shape[0]

np.savetxt('balls_0.txt', balls_file, fmt=' %1.5e')
