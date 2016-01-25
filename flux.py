import numpy as np
import os

num_steps = 1000
simulation_time = 1.0  # in ns

os.chdir('CAS')
flux_file = open('flux_total.txt', 'w')
a_to_b = 0.0
b_to_a = 0.0

for i in range(1, num_steps):
    flux_matrix = np.loadtxt('flux_' + str(i) + '.txt')
    if flux_matrix[0, 0]+flux_matrix[0, 1] != 0.0:
        a_to_b = a_to_b + flux_matrix[0, 1]/((flux_matrix[0, 0]+flux_matrix[0, 1])*simulation_time)
    if flux_matrix[1, 0]+flux_matrix[1, 1] != 0.0:
        b_to_a = b_to_a + flux_matrix[1, 0]/((flux_matrix[1, 0]+flux_matrix[1, 1])*simulation_time)
    rate_a_to_b = a_to_b/i
    rate_b_to_a = b_to_a/i
    flux_file.write(str(i*simulation_time) + ' ' + str(rate_a_to_b) + ' ' + str(rate_b_to_a))
    flux_file.write('\n')

flux_file.close()
#os.system('mv flux_total.txt ..')
