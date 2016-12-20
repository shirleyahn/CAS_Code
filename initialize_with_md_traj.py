import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError

num_states = 2
num_cvs = 6
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values')
traj_file = traj_file[0:6002,:]
state_file = np.loadtxt('initial_states')
state_file = state_file[0:6002]
radius = 80.0
num_clusters = 1


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
    elif -180.0 <= coordinates[0] <= -55.0 and (105.0 <= coordinates[1] <= 180.0 or -180.0 <= coordinates[1] <= -155.0) \
         -180.0 <= coordinates[2] <= -55.0 and (105.0 <= coordinates[3] <= 180.0 or -180.0 <= coordinates[3] <= -155.0) \
         -180.0 <= coordinates[4] <= -55.0 and (105.0 <= coordinates[5] <= 180.0 or -180.0 <= coordinates[5] <= -155.0):
        return 1
    else:
        return -1


# first, separate out folded and unfolded states from the rest for partitioning
folded_balls = np.zeros((1, num_cvs))
folded_ball_to_traj = np.zeros((1,))
unfolded_balls = np.zeros((1, num_cvs))
unfolded_ball_to_traj = np.zeros((1,))
clusters = np.zeros((traj_file.shape[0],))
balls = np.zeros((1, num_cvs))
ball_to_traj = np.zeros((1,))
states = np.zeros((1,))
folded_present = 0
unfolded_present = 0
neither_present = 0
for i in range(traj_file.shape[0]):
    clusters[i] = check_state_function(traj_file[i])
    if clusters[i] == 0:
        if folded_present == 0:
            folded_balls[0] = traj_file[i]
            folded_ball_to_traj[0] = i*time_step
            folded_present += 1
        else:
            folded_balls = np.append(folded_balls, [np.asarray(traj_file[i])], axis=0)
            folded_ball_to_traj = np.append(folded_ball_to_traj, [np.asarray(i*time_step)], axis=0)
    elif clusters[i] == 1:
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

# second, cover the rest of the states with macrostates
new_balls = np.zeros((1, num_cvs))
new_ball_to_traj = np.zeros((1,))
new_balls_count = np.zeros((1,))
new_states = np.zeros((1,))
new_balls_assignment = np.zeros((balls.shape[0],))
current_num_balls = 0
for i in range(balls.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        new_balls[current_num_balls] = balls[i]
        new_ball_to_traj[current_num_balls] = ball_to_traj[i]
        new_balls_count[current_num_balls] += 1
        new_states[current_num_balls] = states[i]
        new_balls_assignment[i] = current_num_balls
        current_num_balls += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(balls[i], new_balls, num_cvs)
        current_ball_center = new_balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, balls[i])
        if distance_from_center <= radius:
            inside += 1
            new_balls_count[ball_key] += 1
            new_balls_assignment[i] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            new_balls = np.append(new_balls, [np.asarray(balls[i])], axis=0)
            new_ball_to_traj = np.append(new_ball_to_traj, [np.asarray(ball_to_traj[i])], axis=0)
            new_balls_count = np.append(new_balls_count, [np.asarray(1)], axis=0)
            new_states = np.append(new_states, [np.asarray(states[i])], axis=0)
            new_balls_assignment[i] = current_num_balls
            current_num_balls += 1

# third, loop through all of the walkers once more to assign them to their true nearest macrostates
for i in range(balls.shape[0]):
    old_ball_key = int(new_balls_assignment[i])
    new_ball_key = closest_ball(balls[i], new_balls, num_cvs)
    new_balls_count[old_ball_key] -= 1
    new_balls_count[new_ball_key] += 1
    new_balls_assignment[i] = new_ball_key

# fourth, delete empty macrostates
delete_list = []
for i in range(new_balls.shape[0]):
    if new_balls_count[i] == 0:
        delete_list.append(i)
new_balls = np.delete(new_balls, delete_list, 0)
new_ball_to_traj = np.delete(new_ball_to_traj, delete_list, 0)
new_balls_count = np.delete(new_balls_count, delete_list, 0)
new_states = np.delete(new_states, delete_list, 0)
new_balls_assignment += num_states

# fifth, calculate the transition matrix and its eigenvalues and eigenvectors
transition_matrix = np.zeros((new_balls.shape[0]+num_states, new_balls.shape[0]+num_states))
current_index = 0
for i in range(traj_file.shape[0]-1):
    if clusters[i] == 0:
        previous_ball_key = 0
    elif clusters[i] == 1:
        previous_ball_key = 1
    else:
        previous_ball_key = int(new_balls_assignment[current_index])

    if clusters[i+1] == 0:
        current_ball_key = 0
    elif clusters[i+1] == 1:
        current_ball_key = 1
    else:
        if clusters[i] == -1:
            current_index += 1
        current_ball_key = int(new_balls_assignment[current_index])

    transition_matrix[previous_ball_key][current_ball_key] += 1.0

new_transition_matrix = np.zeros((transition_matrix.shape[0], transition_matrix.shape[0]))
for i in range(new_transition_matrix.shape[0]):
    for j in range(new_transition_matrix.shape[1]):
        new_transition_matrix[i][j] = (transition_matrix[i][j] + transition_matrix[j][i]) / 2.0

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
normalized_second_evector = np.zeros((final_evectors.shape[0]-num_states, 1))
for i in range(final_evectors.shape[0]-num_states):
    if final_evectors[i+num_states, 0] != 0.0:
        normalized_second_evector[i] = final_evectors[i+num_states, 1] / abs(final_evectors[i+num_states, 0])
    else:
        normalized_second_evector[i] = 0.0

if np.min(normalized_second_evector) != 0.0:
    normalized_second_evector /= np.min(normalized_second_evector)  # normalize
np.savetxt('committor_function.txt', normalized_second_evector, fmt=' %1.10e')

# seventh, cluster states with committor function
if num_clusters > 1:
    while True:
        try:
            centroids, labels = kmeans2(normalized_second_evector, num_clusters, minit='points', iter=200, missing='raise')
            break
        except ClusterError:
            num_clusters -= 1
        if num_clusters <= 1:
            break

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

balls[1] = unfolded_balls[1]
ball_to_traj[1] = unfolded_ball_to_traj[1]
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
    for i in range(num_clusters):
        first = 0
        neither_balls = np.zeros((1, num_cvs))
        neither_ball_to_traj = np.zeros((1,))
        neither_states = np.zeros((1,))
        for j in range(new_balls.shape[0]):
            if labels[j] == i and first == 0:
                first += 1
                neither_balls[0] = new_balls[j]
                neither_ball_to_traj[0] = new_ball_to_traj[j]
                neither_states[0] = new_states[j]
            elif labels[j] == i and first != 0:
                neither_balls = np.append(neither_balls, [np.asarray(new_balls[j])], axis=0)
                neither_ball_to_traj = np.append(neither_ball_to_traj, [np.asarray(new_ball_to_traj[j])], axis=0)
                neither_states = np.append(neither_states, [np.asarray(new_states[j])], axis=0)
        for j in range(neither_balls.shape[0]):
            index = starting_index+j
            balls[index] = neither_balls[j]
            balls_file[index,0:num_cvs] = neither_balls[j]
            balls_file[index,num_cvs] = starting_index
            balls_file[index,num_cvs+1] = num_states+i
            ball_to_traj[index] = neither_ball_to_traj[j]
            states[index] = neither_states[j]
        starting_index += neither_balls.shape[0]

np.savetxt('balls_0.txt', balls_file, fmt=' %1.10e')

for i in range(balls.shape[0]):
    print check_state_function(balls[i]), balls[i], states[i]
