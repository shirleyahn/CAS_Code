import numpy as np
import os
import re
numbers = re.compile(r'(\d+)')


def numerical_sort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


version = 13
new_version = 13
flux_file = open('flux_total_180_v'+str(new_version)+'.txt', 'w')
mfpt_forward_file = open('mfpt_forward_180_v'+str(new_version)+'.txt', 'w')
mfpt_backward_file = open('mfpt_backward_180_v'+str(new_version)+'.txt', 'w')
simulation_time = 0.5  # in ns
initial_step = 1 
os.chdir('a180_'+str(version)+'/CAS')
total_balls_file = np.loadtxt('total_weight.txt')
total_balls = 0 
num_min1_trans = 1 
num_min2_trans = 1 
a_to_b = 0.0 
b_to_a = 0.0 
rate_a_to_b = 0.0 
rate_b_to_a = 0.0 


for subdir, dirs, files in os.walk(os.getcwd()):
    for file in sorted(files, key=numerical_sort):
        if file[0:5] == "flux_" and int(file[5:-4]) == initial_step: #and int(file[5:-4]) <= 20: 
            flux_matrix = np.loadtxt(file)
            if np.sum(flux_matrix[0, :]) != 0.0 and np.sum(flux_matrix[1, :]) == 0.0:  # just min1
                new_a_to_b = flux_matrix[0, 1]/(np.sum(flux_matrix[0, :])*simulation_time)
                a_to_b = a_to_b + new_a_to_b
                mfpt_forward_file.write(str(new_a_to_b)+'\n')
                rate_a_to_b = a_to_b/num_min1_trans
                num_min1_trans += 1
            elif np.sum(flux_matrix[0, :]) == 0.0 and np.sum(flux_matrix[1, :]) != 0.0:  # just min2
                new_b_to_a = flux_matrix[1, 0]/(np.sum(flux_matrix[1, :])*simulation_time)
                b_to_a = b_to_a + new_b_to_a
                mfpt_backward_file.write(str(new_b_to_a)+'\n')
                rate_b_to_a = b_to_a/num_min2_trans
                num_min2_trans += 1
            elif np.sum(flux_matrix[0, :]) != 0.0 and np.sum(flux_matrix[1, :]) != 0.0:  # both
                new_a_to_b = flux_matrix[0, 1]/(np.sum(flux_matrix[0, :])*simulation_time)
                a_to_b = a_to_b + new_a_to_b
                mfpt_forward_file.write(str(new_a_to_b)+'\n')
                rate_a_to_b = a_to_b/num_min1_trans
                num_min1_trans += 1
                new_b_to_a = flux_matrix[1, 0]/(np.sum(flux_matrix[1, :])*simulation_time)
                b_to_a = b_to_a + new_b_to_a
                mfpt_backward_file.write(str(new_b_to_a)+'\n')
                rate_b_to_a = b_to_a/num_min2_trans
                num_min2_trans += 1
            #if total_balls_file[initial_step-1, 3] == 0:
            total_balls += total_balls_file[initial_step-1, -1] 
            #else:
                #total_balls += total_balls_file[initial_step-1, 3]*800/10
            flux_file.write(str(initial_step*simulation_time) + ' ' + str(total_balls) + ' ' + str(rate_a_to_b) + ' ' + str(rate_b_to_a) + "\n")
            initial_step += 1

flux_file.close()
mfpt_forward_file.close()
mfpt_backward_file.close()
