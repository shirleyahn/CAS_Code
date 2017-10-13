import numpy as np
from scipy.sparse import lil_matrix, linalg


def one_dim_energy_function(x):
    return ((x+5.0)**2*(x-5.0)**2)/1000.0 + 3*np.exp(-x**2/10.0) - x/10.0


def two_dim_energy_function(x, y):
    return np.exp(-x**2)+y**2


def x_coord(x_step, x_index):
    return x_step*x_index+x_left


def y_coord(y_step, y_index):
    return y_step*y_index+y_down


def check_state_function(x, y):
    if -1.0 <= x <= 0.0 and -1.0 <= y <= 1.0 and np.sqrt((x+1.0)**2 + y**2) <= 0.4:
        return 0
    elif 0.0 <= x <= 1.0 and -1.0 <= y <= 1.0 and np.sqrt((x-1.0)**2 + y**2) <= 0.4:
        return 1
    else:
        return -1


beta = 10.0
delta_x = 0.05
delta_y = 0.05
delta_t = 1.0
x_right = 1.0
x_left = -1.0
x_total = x_right-x_left
y_up = 1.0
y_down = -1.0
y_total = y_up-y_down
trans_mat = lil_matrix(((int(x_total/delta_x)+1)*(int(y_total/delta_y)+1),
                        (int(x_total/delta_x)+1)*(int(y_total/delta_y)+1)), dtype='float64')
mod_num_x = int(x_total/delta_x)+1
mod_num_y = int(y_total/delta_y)+1
states_index = np.zeros(trans_mat.shape[0])
states_eq_weights = np.zeros(2)

for i in range(trans_mat.shape[0]):
    x = x_coord(delta_x, i / mod_num_y)
    y = y_coord(delta_y, i % mod_num_y)
    states_index[i] = check_state_function(x, y)
    for j in range(trans_mat.shape[1]):
        if i/mod_num_y == j/mod_num_y and abs(j%mod_num_y-i%mod_num_y) == 1:  # up and down
            old_x = x_coord(delta_x, i/mod_num_y)
            old_y = y_coord(delta_y, i%mod_num_y)
            new_x = x_coord(delta_x, j/mod_num_y)
            new_y = y_coord(delta_y, j%mod_num_y)            
            old_energy = two_dim_energy_function(old_x, old_y)
            new_energy = two_dim_energy_function(new_x, new_y)
            if new_energy - old_energy <= 0.0:
                trans_mat[i, j] = 0.25
            else:
                trans_mat[i, j] = np.exp(-beta*(new_energy-old_energy))*0.25
        elif abs(j-i) == mod_num_y:  # right and left
            old_x = x_coord(delta_x, i/mod_num_y)
            old_y = y_coord(delta_y, i%mod_num_y)
            new_x = x_coord(delta_x, j/mod_num_y)
            new_y = y_coord(delta_y, j%mod_num_y)
            old_energy = two_dim_energy_function(old_x, old_y)
            new_energy = two_dim_energy_function(new_x, new_y)
            if new_energy - old_energy <= 0.0:
                trans_mat[i, j] = 0.25
            else:
                trans_mat[i, j] = np.exp(-beta*(new_energy-old_energy))*0.25

row_sum = trans_mat.sum(axis=1)
for i in range(trans_mat.shape[0]):
    trans_mat[i, i] = 1.0 - row_sum[i]
trans_mat = trans_mat.tocsr()

evalues, evectors = linalg.eigs(trans_mat.T, k=10)
np.savetxt('evalues.txt', evalues, fmt=' %1.10e')
np.savetxt('evectors.txt', evectors, fmt=' %1.10e')

eq_weights = np.zeros(trans_mat.shape[0])
for i in range(trans_mat.shape[0]):
    eq_weights[i] = abs(np.real(evectors[i, 0]))
eq_weights /= np.sum(eq_weights)

for i in range(trans_mat.shape[0]):
    if states_index[i] > -1:
        states_eq_weights[int(states_index[i])] += eq_weights[i]
np.savetxt('states_eq_weights.txt', states_eq_weights, fmt=' %1.10e')

rates = np.zeros(evalues.shape[0])
for i in range(evalues.shape[0]):
    if evalues[i] != 0.0:
        rates[i] = (-np.log(np.real(evalues[i]))/delta_t)*0.5
np.savetxt('exact_rates.txt', rates, fmt=' %1.10e')

