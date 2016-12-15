import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError

num_states = 2
num_cvs = 6
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values')
traj_file = traj_file[0:6002,:]
state_file = np.loadtxt('initial_states')
state_file = state_file[0:6002]
num_clusters_kmeans = 200
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
while True:
    try:
        centroids, labels = kmeans2(balls, num_clusters_kmeans, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        num_clusters_kmeans -= 1
    if num_clusters_kmeans <= 1:
        break

labels += num_states
new_balls = np.zeros((1, num_cvs))
new_ball_to_traj = np.zeros((1,))
new_states = np.zeros((1,))
neither_present = 0
for i in range(centroids.shape[0]):
    ball_key = closest_ball(centroids[i], balls, num_cvs)
    if neither_present == 0:
        new_balls[0] = balls[ball_key]
        new_ball_to_traj[0] = ball_to_traj[ball_key]
        new_states[0] = states[ball_key]
        neither_present += 1
    else:
        new_balls = np.append(new_balls, [np.asarray(balls[ball_key])], axis=0)
        new_ball_to_traj = np.append(new_ball_to_traj, [np.asarray(ball_to_traj[ball_key])], axis=0)
        new_states = np.append(new_states, [np.asarray(states[ball_key])], axis=0)

# third, calculate the transition matrix and its eigenvalues and eigenvectors
transition_matrix = np.zeros((new_balls.shape[0]+num_states, new_balls.shape[0]+num_states))
current_index = 0
for i in range(traj_file.shape[0]-1):
    if clusters[i] == 0:
        previous_ball_key = 0
    elif clusters[i] == 1:
        previous_ball_key = 1
    else:
        previous_ball_key = int(labels[current_index])

    if clusters[i+1] == 0:
        current_ball_key = 0
    elif clusters[i+1] == 1:
        current_ball_key = 1
    else:
        if clusters[i] == -1:
            current_index += 1
        current_ball_key = int(labels[current_index])

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

# fourth, normalize the second evector by the first evector values -> good approximation to committor functions.
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
while True:
    try:
        folded_centroids, folded_labels = kmeans2(folded_balls, 1, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        break

while True:
    try:
        unfolded_centroids, unfolded_labels = kmeans2(unfolded_balls, 1, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        break

while True:
    try:
        centroids, labels = kmeans2(normalized_second_evector, num_clusters, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        num_clusters -= 1
    if num_clusters <= 1:
        break

balls = np.zeros((num_states, num_cvs))
ball_to_traj = np.zeros((num_states,))
states = np.zeros((num_states,))

folded_centroids = folded_centroids.reshape((num_cvs,))
ball_key = closest_ball(folded_centroids, folded_balls, num_cvs)
balls[0] = folded_balls[ball_key]
ball_to_traj[0] = folded_ball_to_traj[ball_key]
states[0] = 0

unfolded_centroids = unfolded_centroids.reshape((num_cvs,))
ball_key = closest_ball(unfolded_centroids, unfolded_balls, num_cvs)
balls[1] = unfolded_balls[ball_key]
ball_to_traj[1] = unfolded_ball_to_traj[ball_key]
states[1] = 1

for i in range(centroids.shape[0]):
    ball_key = closest_ball(centroids[i], normalized_second_evector, 1)
    balls = np.append(balls, [np.asarray(new_balls[ball_key])], axis=0)
    ball_to_traj = np.append(ball_to_traj, [np.asarray(new_ball_to_traj[ball_key])], axis=0)
    states = np.append(states, [np.asarray(new_states[ball_key])], axis=0)

np.savetxt('initial_values.txt', balls, fmt=' %1.10e')
np.savetxt('ball_to_traj.txt', ball_to_traj, fmt=' %1.10e')
np.savetxt('initial_states.txt', states, fmt=' %d')

balls_file = np.zeros((balls.shape[0], num_cvs+3))
balls_file[:,0:num_cvs] = balls
balls_file[:,num_cvs] = num_clusters_kmeans
balls_file[:,num_cvs+1] = np.arange(balls.shape[0])
balls_file[:,num_cvs+2] = 0
np.savetxt('balls_0.txt', balls_file, fmt=' %1.10e')

for i in range(balls.shape[0]):
    print check_state_function(balls[i]), balls[i], states[i]
