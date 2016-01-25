import numpy as np
import os

os.chdir('CAS')
flux_file = open('flux_total.txt', 'w')
a_to_b = 0.0
b_to_a = 0.0

for i in range(1, 1000):
    flux_matrix = np.loadtxt('flux_' + str(i) + '.txt')
    if flux_matrix[0, 0]+flux_matrix[0, 1] != 0.0:
        a_to_b = a_to_b + flux_matrix[0, 1]/(flux_matrix[0, 0]+flux_matrix[0, 1])
    if flux_matrix[1, 0]+flux_matrix[1, 1] != 0.0:
        b_to_a = b_to_a + flux_matrix[1, 0]/(flux_matrix[1, 0]+flux_matrix[1, 1])
    rate_a_to_b = a_to_b/i
    rate_b_to_a = b_to_a/i
    flux_file.write(str(i) + ' ' + str(rate_a_to_b) + ' ' + str(rate_b_to_a))
    flux_file.write('\n')

flux_file.close()
os.system('mv flux_total.txt ..')
