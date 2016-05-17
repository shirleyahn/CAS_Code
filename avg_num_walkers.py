import numpy as np

num_walkers = 100
num_walkers_for_sc = 0
version = '1ball_10steps_188_1'
total_weight = np.loadtxt('/Users/Ahn/Documents/CAS/CAS_'+version+'/total_weight.txt')
total_walkers = 0
for i in range(total_weight.shape[0]):
    #if total_weight[i, 3] == 0.0 and total_weight[i, 4] == 0.0:
    total_walkers += total_weight[i, 2]*num_walkers
    #else:
        #total_walkers += total_weight[i, 3]*num_walkers_for_sc + total_weight[i, 4]*num_walkers
print total_walkers/total_weight[-1, 0]
