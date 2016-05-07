import numpy as np
import os
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

version = '1ball_10steps_sc_200_40_1'
simulation_time = 1.0  # in ns
num_steps = 1
flux_file = open('flux_total_'+str(version)+'.txt', 'w')
os.chdir('CAS_'+str(version))
a_to_b = 0.0
b_to_a = 0.0

for subdir, dirs, files in os.walk(os.getcwd()):
    for file in sorted(files, key=numericalSort):
        if file[0:5] == "flux_" and int(file[5:-4]) == num_steps:
            flux_matrix = np.loadtxt(file)
            if np.sum(flux_matrix[0, :]) != 0.0:
                a_to_b = a_to_b + (np.sum(flux_matrix[0, 1])/(np.sum(flux_matrix[0, :])*simulation_time))
            if np.sum(flux_matrix[1, :]) != 0.0:
                b_to_a = b_to_a + (np.sum(flux_matrix[1, 0])/(np.sum(flux_matrix[1, :])*simulation_time))
            rate_a_to_b = a_to_b/num_steps
            rate_b_to_a = b_to_a/num_steps
            flux_file.write(str(num_steps*simulation_time) + ' ' + str(rate_a_to_b) + ' ' + str(rate_b_to_a))
            flux_file.write('\n')
            num_steps += 1

flux_file.close()
