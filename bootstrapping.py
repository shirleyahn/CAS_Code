import numpy as np

mfpt_file = np.loadtxt('mfpt_1.txt')
total_sim_length = 24  # in microseconds
fluxes = np.zeros((100, ))
flux_final = np.zeros((total_sim_length, 3))

for i in range(total_sim_length):
    sample_size = int(((i+1)/float(total_sim_length))*mfpt_file.shape[0])
    if sample_size != 0:
        for j in range(fluxes.shape[0]):
            mfpt_array = np.zeros((sample_size,))
            indices = np.random.choice(mfpt_file.shape[0], sample_size, replace=True)
            for k in range(sample_size):
                mfpt_array[k] = mfpt_file[indices[k]]
            fluxes[j] = 1.0/np.mean(mfpt_array)
        flux_final[i, 0] = (i+1)*1000.0  # in nanoseconds
        flux_final[i, 1] = np.mean(fluxes)
        flux_final[i, 2] = np.std(fluxes)

np.savetxt('flux_1.txt', flux_final, fmt=' %1.5e')
