import numpy as np


def one_dim_energy_function(x):
    return ((x+5.0)**2*(x-5.0)**2)/1000.0 + 3*np.exp(-x**2/10.0) - x/10.0


def two_dim_energy_function(x, y):
    return np.exp(-x**2)+y**2


def x_coord(x_step, x_index):
    return x_step*x_index+x_left


def y_coord(y_step, y_index):
    return y_step*y_index+y_down


beta = 1.0/0.1
delta_x = 0.02
delta_y = 0.02
delta_t = 1.0
x_right = 1.0
x_left = -1.0
x_total = x_right-x_left
y_up = 1.0
y_down = -1.0
y_total = y_up-y_down
trans_mat = np.zeros((int(x_total/delta_x)*int(y_total/delta_y), int(x_total/delta_x)*int(y_total/delta_y)))

for i in range(0, trans_mat.shape[0]):
    for j in range(0, trans_mat.shape[1]):
        if i != j:
            old_x_index = i / int(y_total/delta_y)
            old_y_index = i % int(y_total/delta_y)
            new_x_index = j / int(y_total/delta_y)
            new_y_index = j % int(y_total/delta_y)
            if (old_x_index == new_x_index and abs(new_y_index-old_y_index) == 1) \
                    or (abs(new_x_index-old_x_index) == 1 and old_y_index == new_y_index):
                old_x = x_coord(delta_x, old_x_index)
                old_y = y_coord(delta_y, old_y_index)
                new_x = x_coord(delta_x, new_x_index)
                new_y = y_coord(delta_y, new_y_index)
                old_energy = np.exp(-beta*two_dim_energy_function(old_x, old_y))
                new_energy = np.exp(-beta*two_dim_energy_function(new_x, new_y))
                trans_mat[i, j] = np.minimum(1.0, new_energy/old_energy)*0.25

row_sum = np.sum(trans_mat, axis=1)
for i in range(trans_mat.shape[0]):
    trans_mat[i, i] = 1.0 - row_sum[i]

evalues, evectors = np.linalg.eig(trans_mat.T)
idx = abs(evalues).argsort()[::-1]
evalues = evalues[idx]
final_evalues = np.real(evalues)
evectors = evectors[:, idx]
final_evectors = np.real(evectors)
np.savetxt('trans_mat.txt', trans_mat, fmt=' %1.10e')
np.savetxt('evalues.txt', final_evalues, fmt=' %1.10e')
np.savetxt('evectors.txt', final_evectors, fmt=' %1.10e')
np.savetxt('rates.txt', -np.log(final_evalues)/delta_t, fmt=' %1.10e')
