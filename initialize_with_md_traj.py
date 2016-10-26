import numpy as np

num_cvs = 6
radius = 100.0
time_step = 500.0  # in ps
traj_file = np.loadtxt('initial_values.txt')

def calculate_distance_from_center(center, values):
    distance = 0.0
    for i in range(len(center)):
        distance += min(360.0 - abs(values[i] - center[i]), abs(values[i] - center[i])) ** 2
    if abs(distance) < 1.0e-10:
        distance = 0.0
    return np.sqrt(distance)


def closest_ball(walker_coordinates, balls):
    distance = np.zeros((balls.shape[0],))
    for i in range(num_cvs):
        distance += np.minimum(360.0 - np.abs(balls[:, i] - walker_coordinates[i]), np.abs(balls[:, i] - walker_coordinates[i])) ** 2
    #distance = np.sum((balls_array - walker_coordinates)**2, axis=1)
    return np.argmin(distance)


# first, partition the free energy landscape by creating macrostates
balls = np.zeros((1, num_cvs))
balls_array = np.zeros((traj_file.shape[0],))
ball_to_traj = np.zeros((1,))
current_num_balls = 0
for i in range(traj_file.shape[0]):
    inside = 0
    # if we're dealing with the very first walker, create the very first macrostate for the walker.
    if i == 0:
        inside += 1
        current_ball_center = traj_file[i]
        balls[current_num_balls] = np.asarray(current_ball_center)
        balls_array[i] = current_num_balls
        ball_to_traj[current_num_balls] = np.asarray(i*time_step)
        current_num_balls += 1

    # otherwise, loop through the existing macrostates and find the macrostate with a center nearest to the walker.
    if inside == 0:
        ball_key = closest_ball(traj_file[i], balls)
        current_ball_center = balls[ball_key].tolist()
        distance_from_center = calculate_distance_from_center(current_ball_center, traj_file[i])
        if distance_from_center <= radius or abs(distance_from_center - radius) < 1.0e-10:
            inside += 1
            balls_array[i] = ball_key

        # walker is not inside any macrostate, so create a new macrostate centered around the walker.
        if inside == 0:
            current_ball_center = traj_file[i]
            balls = np.append(balls, [np.asarray(current_ball_center)], axis=0)
            balls_array[i] = current_num_balls
            ball_to_traj = np.append(ball_to_traj, [np.asarray(i*time_step)], axis=0)
            current_num_balls += 1

# second, calculate the transition matrix of the macrostates
transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
for i in range(traj_file.shape[0]-1):
    previous_ball_key = balls_array[i]
    current_ball_key = balls_array[i+1]
    transition_matrix[previous_ball_key][current_ball_key] += 1

# transition matrix should fulfill detailed balance if simulation is run under Hamiltonian dynamics in the
# canonical ensemble. equation is from Prinz, et al JCP (2011).
new_transition_matrix = np.zeros((balls.shape[0], balls.shape[0]))
for i in range(new_transition_matrix.shape[0]):
    for j in range(new_transition_matrix.shape[1]):
        new_transition_matrix[i][j] = (transition_matrix[i][j] + transition_matrix[j][i])/2.0

row_sum = np.sum(new_transition_matrix, axis=1)
for i in range(new_transition_matrix.shape[0]):
    if row_sum[i] != 0.0:
        new_transition_matrix[i, :] /= row_sum[i]
np.savetxt('transition_matrix.txt', new_transition_matrix, fmt=' %1.10e')

evalues, evectors = np.linalg.eig(new_transition_matrix.T)
idx = abs(evalues).argsort()[::-1]
evalues = evalues[idx]
final_evalues = np.real(evalues)
evectors = evectors[:, idx]
final_evectors = np.real(evectors)
np.savetxt('evalues.txt', final_evalues, fmt=' %1.10e')
np.savetxt('evectors.txt', final_evectors, fmt=' %1.10e')
np.savetxt('eq_weights.txt', abs(final_evectors[:, 0]), fmt=' %1.10e')
np.savetxt('ball_to_traj.txt', ball_to_traj, fmt=' %1.10e')
