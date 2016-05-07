import numpy as np
import os

num_steps = 5000
simulation_time = 1.0  # in ns

version = '2balls_10steps_beta_10_3'
os.chdir('CAS_'+str(version))
flux_file = open('flux_total_'+str(version)+'.txt', 'w')
a_to_b_via_c = 0.0
a_to_b_via_d = 0.0
b_to_a_via_c = 0.0
b_to_a_via_d = 0.0

for i in range(1, num_steps+1):
    flux_matrix = np.loadtxt('flux_' + str(i) + '.txt')
    if np.sum(flux_matrix[0:2, :]) != 0.0:
        a_to_b_via_c = a_to_b_via_c + np.sum(flux_matrix[0, 2:4])/(np.sum(flux_matrix[0:2, :])*simulation_time)
        a_to_b_via_d = a_to_b_via_d + np.sum(flux_matrix[1, 2:4])/(np.sum(flux_matrix[0:2, :])*simulation_time)
    if np.sum(flux_matrix[2:4, :]) != 0.0:
        b_to_a_via_c = b_to_a_via_c + np.sum(flux_matrix[2, 0:2])/(np.sum(flux_matrix[2:4, :])*simulation_time)
        b_to_a_via_d = b_to_a_via_d + np.sum(flux_matrix[3, 0:2])/(np.sum(flux_matrix[2:4, :])*simulation_time)
    rate_a_to_b_via_c = a_to_b_via_c/i
    rate_a_to_b_via_d = a_to_b_via_d/i
    rate_b_to_a_via_c = b_to_a_via_c/i
    rate_b_to_a_via_d = b_to_a_via_d/i
    flux_file.write(str(i*simulation_time) + ' ' + str(rate_a_to_b_via_c) + ' ' + str(rate_a_to_b_via_d) + ' '+ str(rate_b_to_a_via_c) + ' ' + str(rate_b_to_a_via_d))
    flux_file.write('\n')

flux_file.close()
os.system('mv flux_total_'+str(version)+'.txt ..')
