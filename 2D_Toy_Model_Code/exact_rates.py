import numpy as np
from time import time
from scipy.sparse import lil_matrix, linalg, identity


def one_dim_energy_function(x):
    return ((x+5.0)**2*(x-5.0)**2)/1000.0 + 3*np.exp(-x**2/10.0) - x/10.0


def two_dim_energy_function(x, y):
    return np.exp(-x**2)+y**2
    #return 2.25*(x-0.5*np.sin(2*np.pi*y))**2+2.25*np.cos(2*np.pi*y)


def x_coord(x_step, x_index):
    return x_step*x_index+x_left


def y_coord(y_step, y_index):
    return y_step*y_index+y_down


beta = 15.0
delta_x = 0.05
delta_y = 0.05
delta_t = 1.0
x_right = 1.5
x_left = -1.5
x_total = x_right-x_left
y_up = 1.0
y_down = -1.0
y_total = y_up-y_down
trans_mat = lil_matrix(((int(x_total/delta_x)+1)*(int(y_total/delta_y)+1),
                        (int(x_total/delta_x)+1)*(int(y_total/delta_y)+1)), dtype='float64')
mod_num_x = int(x_total/delta_x)+1
mod_num_y = int(y_total/delta_y)+1

t0 = time()
for i in range(trans_mat.shape[0]):
    for j in range(trans_mat.shape[1]):
        if (i/mod_num_y == j/mod_num_y and abs(j%mod_num_y-i%mod_num_y) == 1):  # up and down
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
        elif (abs(j-i) == mod_num_y):  # right and left
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
t1 = time()

evalues, evectors = linalg.eigs(trans_mat.T, k=2)
np.savetxt('evalues.txt', evalues, fmt=' %1.10e')
np.savetxt('evectors.txt', evectors, fmt=' %1.10e')

'''
## inverse iteration
v_0 = lil_matrix(((int(x_total/delta_x)+1)*(int(y_total/delta_y)+1),1), dtype='float64')
v_0[0] = 1.0
eigen = 0.0
for k in range(100):
    w=linalg.inv(trans_mat-0.999*identity((int(x_total/delta_x)+1)*(int(y_total/delta_y)+1))).dot(v_0)
    v_0 = w/np.linalg.norm(w.toarray())
    eigen = v_0.T.dot(trans_mat.dot(v_0))
    print eigen
'''

rates = np.zeros(evalues.shape[0])
for i in range(evalues.shape[0]):
    if evalues[i] != 0.0:
        rates[i] = (-np.log(np.real(evalues[i]))/delta_t)*0.5

np.savetxt('rates.txt', rates, fmt=' %1.10e')
t2 = time()
print t2-t1, t1-t0
