import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError, whiten

num_states = 2
num_cvs = 6
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values')
traj_file = traj_file[0:6002,:]
state_file = np.loadtxt('initial_states')
state_file = state_file[0:6002]
radius = 160.0
num_clusters = 10


def convert_angles_to_cos_sin(balls):
    new_balls = np.zeros((balls.shape[0], 2*balls.shape[1]))
    for i in range(balls.shape[0]):
        for j in range(balls.shape[1]):
            new_balls[i, 2*j] = np.cos(balls[i, j]*np.pi/180)
            new_balls[i, 2*j+1] = np.sin(balls[i, j]*np.pi/180)
    return new_balls


def closest_ball(coordinates, balls, num_cvs):
    distance = np.zeros((balls.shape[0],))
    for i in range(num_cvs):
        distance += np.minimum(360.0 - np.abs(balls[:, i] - coordinates[i]), np.abs(balls[:, i] - coordinates[i])) ** 2
    return np.argmin(distance)


def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += min(360.0 - abs(values[i] - center[i]), abs(values[i] - center[i])) ** 2
    if abs(distance) < 1.0e-10:
        distance = 0.0
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


# first, cover the states with macrostates
balls = np.zeros((1, num_cvs))
ball_to_traj = np.zeros((1,))
balls_count = np.zeros((1,))
states = np.zeros((1,))
balls_assignment = np.zeros((traj_file.shape[0],))
current_num_balls = 0
for i in range(traj_file.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        balls[current_num_balls] = np.asarray(traj_file[i])
        ball_to_traj[current_num_balls] = i*time_step
        balls_count[current_num_balls] += 1
        states[current_num_balls] = state_file[i]
        balls_assignment[i] = current_num_balls
        current_num_balls += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(traj_file[i], balls, num_cvs)
        current_ball_center = balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, traj_file[i])
        if distance_from_center <= radius or abs(distance_from_center - radius) < 1.0e-10:
            inside += 1
            balls_count[ball_key] += 1
            balls_assignment[i] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            balls = np.append(balls, [np.asarray(traj_file[i])], axis=0)
            ball_to_traj = np.append(ball_to_traj, [np.asarray(i*time_step)], axis=0)
            balls_count = np.append(balls_count, [np.asarray(1)], axis=0)
            states = np.append(states, [np.asarray(state_file[i])], axis=0)
            balls_assignment[i] = current_num_balls
            current_num_balls += 1

# second, loop through all of the walkers once more to assign them to their true nearest macrostates
for i in range(traj_file.shape[0]):
    old_ball_key = int(balls_assignment[i])
    new_ball_key = closest_ball(traj_file[i], balls, num_cvs)
    balls_count[old_ball_key] -= 1
    balls_count[new_ball_key] += 1
    balls_assignment[i] = new_ball_key

# third, delete empty macrostates
delete_list = []
for i in range(balls.shape[0]):
    if balls_count[i] == 0:
        delete_list.append(i)
balls = np.delete(balls, delete_list, 0)
ball_to_traj = np.delete(ball_to_traj, delete_list, 0)
balls_count = np.delete(balls_count, delete_list, 0)
states = np.delete(states, delete_list, 0)
balls_assignment += num_states

# fourth, calculate the transition matrix and its eigenvalues and eigenvectors
transition_matrix = np.zeros((balls.shape[0]+num_states, balls.shape[0]+num_states))
current_index = 0
for i in range(traj_file.shape[0]-1):
    previous_ball_key = int(balls_assignment[i])
    current_ball_key = int(balls_assignment[i+1])
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

# fifth, normalize the second evector by the first evector values -> good approximation to committor functions.
normalized_second_evector = np.zeros((final_evectors.shape[0], 1))
for i in range(final_evectors.shape[0]):
    if final_evectors[i, 0] != 0.0:
        normalized_second_evector[i] = final_evectors[i, 1] / abs(final_evectors[i, 0])
    else:
        normalized_second_evector[i] = 0.0

if np.min(normalized_second_evector) != 0.0:
    normalized_second_evector /= np.min(normalized_second_evector)  # normalize
np.savetxt('committor_function.txt', normalized_second_evector, fmt=' %1.10e')

# sixth, group macrostates into folded, unfolded, and neither
folded_balls = np.zeros((1, num_cvs))
folded_ball_to_traj = np.zeros((1,))
folded_committor = np.zeros((1, 1))
unfolded_balls = np.zeros((1, num_cvs))
unfolded_ball_to_traj = np.zeros((1,))
unfolded_committor = np.zeros((1, 1))
neither_balls = np.zeros((1, num_cvs))
neither_ball_to_traj = np.zeros((1,))
neither_committor = np.zeros((1, 1))
neither_states = np.zeros((1,))
folded_present = 0
unfolded_present = 0
neither_present = 0
for i in range(balls.shape[0]):
    current_state = check_state_function(balls[i])
    if current_state == 0:
        if folded_present == 0:
            folded_balls[0] = balls[i]
            folded_ball_to_traj[0] = ball_to_traj[i]
            folded_committor[0] = normalized_second_evector[i]
            folded_present += 1
        else:
            folded_balls = np.append(folded_balls, [np.asarray(balls[i])], axis=0)
            folded_ball_to_traj = np.append(folded_ball_to_traj, [np.asarray(ball_to_traj[i])], axis=0)
            folded_committor = np.append(folded_committor, [np.asarray(normalized_second_evector[i])], axis=0)
    elif current_state == 1:
        if unfolded_present == 0:
            unfolded_balls[0] = balls[i]
            unfolded_ball_to_traj[0] = ball_to_traj[i]
            unfolded_committor[0] = normalized_second_evector[i]
            unfolded_present += 1
        else:
            unfolded_balls = np.append(unfolded_balls, [np.asarray(balls[i])], axis=0)
            unfolded_ball_to_traj = np.append(unfolded_ball_to_traj, [np.asarray(ball_to_traj[i])], axis=0)
            unfolded_committor = np.append(unfolded_committor, [np.asarray(normalized_second_evector[i])], axis=0)
    else:
        if neither_present == 0:
            neither_balls[0] = balls[i]
            neither_ball_to_traj[0] = ball_to_traj[i]
            neither_committor[0] = normalized_second_evector[i]
            neither_states[0] = states[i]
            neither_present += 1
        else:
            neither_balls = np.append(neither_balls, [np.asarray(balls[i])], axis=0)
            neither_ball_to_traj = np.append(neither_ball_to_traj, [np.asarray(ball_to_traj[i])], axis=0)
            neither_committor = np.append(neither_committor, [np.asarray(normalized_second_evector[i])], axis=0)
            neither_states = np.append(neither_states, [np.asarray(states[i])], axis=0)

# seventh, cluster states with committor function
if folded_committor.shape[0] > 1:
    while True:
        try:
            folded_centroids, folded_labels = kmeans2(folded_committor, 1, minit='points', iter=200, missing='raise')
            break
        except ClusterError:
            break
else:
    folded_centroids = folded_committor

if unfolded_committor.shape[0] > 1:
    while True:
        try:
            unfolded_centroids, unfolded_labels = kmeans2(unfolded_committor, 1, minit='points', iter=200, missing='raise')
            break
        except ClusterError:
            break
else:
    unfolded_centroids = unfolded_committor

while True:
    try:
        neither_centroids, neither_labels = kmeans2(neither_committor, num_clusters, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        num_clusters -= 1
    if num_clusters <= 1:
        break

new_balls = np.zeros((num_states+num_clusters, num_cvs))
new_ball_to_traj = np.zeros((num_states+num_clusters,))
new_states = np.zeros((num_states+num_clusters,))

ball_key = closest_ball(folded_centroids, folded_committor, 1)
new_balls[0] = folded_balls[ball_key]
new_ball_to_traj[0] = folded_ball_to_traj[ball_key]
new_states[0] = 0

ball_key = closest_ball(unfolded_centroids, unfolded_committor, 1)
new_balls[1] = unfolded_balls[ball_key]
new_ball_to_traj[1] = unfolded_ball_to_traj[ball_key]
new_states[1] = 1

for i in range(neither_centroids.shape[0]):
    ball_key = closest_ball(neither_centroids[i], neither_committor, 1)
    new_balls[num_states+i] = neither_balls[ball_key]
    new_ball_to_traj[num_states+i] = neither_ball_to_traj[ball_key]
    new_states[num_states+i] = neither_states[ball_key]

np.savetxt('initial_values.txt', new_balls, fmt=' %1.10e')
np.savetxt('ball_to_traj.txt', new_ball_to_traj, fmt=' %1.10e')
np.savetxt('initial_states.txt', new_states, fmt=' %d')

balls_file = np.zeros((new_balls.shape[0], num_cvs+3))
balls_file[:,0:num_cvs] = new_balls
balls_file[:,num_cvs] = radius
balls_file[:,num_cvs+1] = np.arange(new_balls.shape[0])
balls_file[:,num_cvs+2] = 0
np.savetxt('balls_0.txt', balls_file, fmt=' %1.10e')

for i in range(new_balls.shape[0]):
    print check_state_function(new_balls[i]), new_balls[i], new_states[i]
