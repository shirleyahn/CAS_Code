import numpy as np
from scipy.cluster.vq import kmeans2, ClusterError, whiten

num_states = 2
num_cvs = 6
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values')
state_file = np.loadtxt('initial_states')
num_clusters = 200
radius = 180.0


def convert_angles_to_cos_sin(balls):
    new_balls = np.zeros((balls.shape[0], 2*balls.shape[1]))
    for i in range(balls.shape[0]):
        for j in range(balls.shape[1]):
            new_balls[i, 2*j] = np.cos(balls[i, j]*np.pi/180)
            new_balls[i, 2*j+1] = np.sin(balls[i, j]*np.pi/180)
    return new_balls


def closest_ball(coordinates, balls):
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
    #elif -180.0 <= coordinates[0] <= -55.0 and (-180.0 <= coordinates[1] <= -155.0) \
         #-180.0 <= coordinates[2] <= -55.0 and (-180.0 <= coordinates[3] <= -155.0) \
         #-180.0 <= coordinates[4] <= -55.0 and (-180.0 <= coordinates[5] <= -155.0):
        #return 2
    else:
        return -1

state_balls = np.zeros((num_states, num_cvs))
state_ball_to_traj = np.zeros((num_states,))
state_states = np.zeros((num_states,))
clusters = np.zeros((traj_file.shape[0],))
new_traj_file = np.zeros((1, num_cvs))
new_ball_to_traj = np.zeros((1,))
new_states = np.zeros((1,))
folded_present = 0
unfolded_present_1 = 0
unfolded_present_2 = 0
first = 0
for i in range(traj_file.shape[0]):
    clusters[i] = check_state_function(traj_file[i])
    if folded_present == 0 and clusters[i] == 0:
        state_balls[0] = traj_file[i]
        state_ball_to_traj[0] = i*time_step
        state_states[0] = 0
        folded_present += 1
    elif unfolded_present_1 == 0 and clusters[i] == 1:
        state_balls[1] = traj_file[i]
        state_ball_to_traj[1] = i*time_step
        state_states[1] = 1
        unfolded_present_1 += 1
    elif unfolded_present_2 == 0 and clusters[i] == 2:
        state_balls[2] = traj_file[i]
        state_ball_to_traj[2] = i*time_step
        state_states[2] = 2
        unfolded_present_2 += 1
    elif clusters[i] == -1:
        if first == 0:
            new_traj_file[0] = traj_file[i]
            new_ball_to_traj[0] = i*time_step
            new_states[0] = state_file[i]
            first += 1
        else:
            new_traj_file = np.append(new_traj_file, [np.asarray(traj_file[i])], axis=0)
            new_ball_to_traj = np.append(new_ball_to_traj, [np.asarray(i*time_step)], axis=0)
            new_states = np.append(new_states, [np.asarray(state_file[i])], axis=0)

clustering_matrix = convert_angles_to_cos_sin(new_traj_file)
while True:
    try:
        centroids, labels = kmeans2(clustering_matrix, num_clusters, minit='points', iter=200, missing='raise')
        break
    except ClusterError:
        num_clusters -= 1
    if num_clusters <= 1:
        break

balls = np.zeros((num_clusters, num_cvs))
ball_to_traj = np.zeros((num_clusters,))
states = np.zeros((num_clusters,))
for i in range(balls.shape[0]):
     ball_key = closest_ball(centroids[i], clustering_matrix)
     balls[i] = new_traj_file[ball_key]
     ball_to_traj[i] = new_ball_to_traj[ball_key]
     states[i] = new_states[ball_key]

new_balls = np.append(state_balls, balls, axis=0)
new_ball_to_traj = np.append(state_ball_to_traj, ball_to_traj, axis=0)
new_states = np.append(state_states, states, axis=0)
labels += num_states

transition_matrix = np.zeros((new_balls.shape[0], new_balls.shape[0]))
current_index = 0
for i in range(traj_file.shape[0]-1):
    if clusters[i] == 0:
        previous_ball_key = 0
    elif clusters[i] == 1:
        previous_ball_key = 1
    elif clusters[i] == 2:
        previous_ball_key = 2
    elif clusters[i] == -1:
        previous_ball_key = labels[current_index]

    if clusters[i+1] == 0:
        current_ball_key = 0
    elif clusters[i+1] == 1:
        current_ball_key = 1
    elif clusters[i+1] == 2:
        current_ball_key = 2
    elif clusters[i+1] == -1:
        if clusters[i] == -1:
            current_index += 1
        current_ball_key = labels[current_index]
   
    transition_matrix[previous_ball_key][current_ball_key] += 1.0

new_transition_matrix = np.zeros((new_balls.shape[0], new_balls.shape[0]))
for i in range(new_transition_matrix.shape[0]):
    for j in range(new_transition_matrix.shape[1]):
        new_transition_matrix[i][j] = (transition_matrix[i][j] + transition_matrix[j][i])/2.0

row_sum = np.sum(new_transition_matrix, axis=1)
for i in range(new_transition_matrix.shape[0]):
    if row_sum[i] != 0.0:
        new_transition_matrix[i, :] /= row_sum[i]

evalues, evectors = np.linalg.eig(new_transition_matrix.T)
idx = abs(evalues).argsort()[::-1]
final_evalues = evalues[idx]
final_evectors = evectors[:, idx]

idx = abs(final_evectors[:, 0]).argsort()[::-1]
final_balls = new_balls[idx, :]
final_ball_to_traj = new_ball_to_traj[idx]
final_states = new_states[idx]

"""
balls = np.zeros((1, num_cvs))
ball_to_traj = np.zeros((1,))
balls_count = np.zeros((1,))
balls_assignment = np.zeros((new_traj_file.shape[0],))
current_num_balls = 0
for i in range(new_traj_file.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        current_ball_center = new_traj_file[i]
        balls[current_num_balls] = np.asarray(current_ball_center)
        ball_to_traj[current_num_balls] = new_ball_to_traj[i]
        balls_count[current_num_balls] += 1
        balls_assignment[i] = current_num_balls
        current_num_balls += 1
 
    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(new_traj_file[i], balls)
        current_ball_center = balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, new_traj_file[i])
        if distance_from_center <= radius or abs(distance_from_center - radius) < 1.0e-10:
            inside += 1
            balls_count[ball_key] += 1
            balls_assignment[i] = ball_key
 
        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            current_ball_center = new_traj_file[i]
            balls = np.append(balls, [np.asarray(current_ball_center)], axis=0)
            ball_to_traj = np.append(ball_to_traj, [np.asarray(new_ball_to_traj[i])], axis=0)
            balls_count = np.append(balls_count, [np.asarray(1)], axis=0)
            balls_assignment[i] = current_num_balls
            current_num_balls += 1

for i in range(new_traj_file[i].shape[0]):
    old_ball_key = int(balls_assignment[i])
    new_ball_key = closest_ball(new_traj_file[i], balls)
    balls_count[old_ball_key] -= 1
    balls_count[new_ball_key] += 1
    balls_assignment[i] = new_ball_key

for i in range(new_traj_file[i].shape[0]):
    old_ball_key = int(balls_assignment[i])
    if balls_count[old_ball_key] == 1:
        mask = np.ones(balls.shape[0], dtype=bool)
        mask[old_ball_key] = 0
        new_ball_key = closest_ball(new_traj_file[i], balls[mask])
        new_ball_center = balls[new_ball_key].tolist()
        distance_from_center = calculate_distance_from_center(new_ball_center, new_traj_file[i])
        if distance_from_center <= radius or abs(distance_from_center - radius) < 1.0e-10:
            balls_count[old_ball_key] -= 1
            balls_count[new_ball_key] += 1
            balls_assignment[i] = new_ball_key

delete_list = []
for i in range(balls.shape[0]):
    if balls_count[i] == 0:
        delete_list.append(i)
np.delete(balls, delete_list, 0)
np.delete(ball_to_traj, delete_list, 0)
np.delete(balls_count, delete_list, 0)

idx = balls_count.argsort()[::-1]
new_balls = balls[idx, :]
new_ball_to_traj = ball_to_traj[idx]
new_balls_count = balls_count[idx]

final_balls = np.append(state_balls, new_balls, axis=0)
final_ball_to_traj = np.append(state_ball_to_traj, new_ball_to_traj, axis=0)
"""

np.savetxt('initial_values.txt', final_balls, fmt=' %1.10e')
np.savetxt('ball_to_traj.txt', final_ball_to_traj, fmt=' %1.10e')
np.savetxt('initial_states.txt', final_states, fmt=' %d')

balls_file = np.zeros((final_balls.shape[0], num_cvs+3))
balls_file[:,0:num_cvs] = final_balls
balls_file[:,num_cvs] = radius
balls_file[:,num_cvs+1] = np.arange(final_balls.shape[0])
balls_file[:,num_cvs+2] = 0

np.savetxt('balls_0.txt', balls_file, fmt=' %1.10e')

for i in range(final_balls.shape[0]):
    print check_state_function(final_balls[i]), final_balls[i], final_states[i]
