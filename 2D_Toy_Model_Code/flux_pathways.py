import numpy as np
import os

num_steps = 70
simulation_time = 1.0  # in ns

#version = 'test'
os.chdir('CAS')
flux_file = open('flux_total.txt', 'w')
a_to_b_via_c = 0.0
a_to_b_via_d = 0.0
b_to_a_via_c = 0.0
b_to_a_via_d = 0.0

for i in range(1, num_steps+1):
    flux_matrix = np.loadtxt('flux_' + str(i) + '.txt')
    if np.sum(flux_matrix[0, :]) != 0.0:
        a_to_b_via_c = a_to_b_via_c + (flux_matrix[0, 2]+flux_matrix[0, 3])/(np.sum(flux_matrix[0, :])*simulation_time)
    if np.sum(flux_matrix[1, :]) != 0.0:
        a_to_b_via_d = a_to_b_via_d + (flux_matrix[1, 2]+flux_matrix[1, 3])/(np.sum(flux_matrix[1, :])*simulation_time)
    if np.sum(flux_matrix[2, :]) != 0.0:
        b_to_a_via_c = b_to_a_via_c + (flux_matrix[2, 0]+flux_matrix[2, 1])/(np.sum(flux_matrix[2, :])*simulation_time)
    if np.sum(flux_matrix[3, :]) != 0.0:
        b_to_a_via_d = b_to_a_via_d + (flux_matrix[3, 0]+flux_matrix[3, 1])/(np.sum(flux_matrix[3, :])*simulation_time)
    rate_a_to_b_via_c = a_to_b_via_c/i
    rate_a_to_b_via_d = a_to_b_via_d/i
    rate_b_to_a_via_c = b_to_a_via_c/i
    rate_b_to_a_via_d = b_to_a_via_d/i
    flux_file.write(str(i*simulation_time) + ' ' + str(rate_a_to_b_via_c) + ' ' + str(rate_a_to_b_via_d) + ' '+ str(rate_b_to_a_via_c) + ' ' + str(rate_b_to_a_via_d))
    flux_file.write('\n')

flux_file.close()
os.system('mv flux_total.txt ..')
