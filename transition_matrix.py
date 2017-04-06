import numpy as np
import os


def compute_transition_matrix(file_directory, size, num_matrices):
    os.chdir(file_directory)
    root_directory = os.getcwd()
    new_indices = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1])

    for num_samples in range(2, num_matrices+1):
        avg_transition_matrix = np.zeros((size, size))
        matrices_array = np.zeros((1, size, size))
        for subdir, dirs, files in os.walk(root_directory):
            for file in files:
                if file[0:18] == 'transition_matrix_' and int(file[18:-4]) <= num_samples:
                    os.chdir(subdir)
                    transition_matrix = np.loadtxt(file)
                    new_transition_matrix = np.zeros((size, size))
                    for i in range(size):
                        for j in range(size):
                            new_transition_matrix[i, j] = transition_matrix[new_indices[i], new_indices[j]]
                    row_sum = np.sum(new_transition_matrix, axis=1)
                    for i in range(new_transition_matrix.shape[0]):
                        if row_sum[i] != 0.0:
                            new_transition_matrix[i, :] /= row_sum[i]
                    avg_transition_matrix += new_transition_matrix
                    new_transition_matrix = new_transition_matrix.reshape((1, size, size))
                    if int(file[18:-4]) == 1:
                        matrices_array = new_transition_matrix
                    else:
                        matrices_array = np.concatenate((matrices_array, new_transition_matrix))

        avg_transition_matrix /= num_samples
        np.savetxt('avg_transition_matrix_'+str(num_samples)+'.txt', avg_transition_matrix, fmt=' %1.10e')

        std_transition_matrix = np.std(matrices_array, axis=0, ddof=0)/np.sqrt(num_samples)
        np.savetxt('std_transition_matrix_'+str(num_samples)+'.txt', std_transition_matrix, fmt=' %1.10e')


if __name__ == '__main__':
    compute_transition_matrix('CAS_180_v4', 12, 8)
