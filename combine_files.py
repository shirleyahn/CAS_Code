import numpy as np
import os


def combine_files(first_file_num, second_file_num):
    file1 = np.loadtxt('ball_clustering_'+str(first_file_num)+'.txt')
    file2 = np.loadtxt('total_weight_on_each_ball_'+str(second_file_num)+'.txt')
    combined_array = np.zeros((file1.shape[0]+file2.shape[0], file2.shape[1]+3))
    combined_array[0:file1.shape[0], 0:file2.shape[1]] = file1[:, 0:file2.shape[1]]  # coordinates
    combined_array[0:file1.shape[0], file2.shape[1]] = file1[:, file2.shape[1]+1]  # eq weights
    combined_array[0:file1.shape[0], file2.shape[1]+1] = file1[:, file2.shape[1]]  # clusters
    combined_array[0:file1.shape[0], file2.shape[1]+2] = file1[:, file2.shape[1]+2]  # committor
    combined_array[0:file1.shape[0], file2.shape[1]+3] = 0  # indicating before
    combined_array[file1.shape[0]:, 0:file2.shape[1]] = file2[:,0:file2.shape[1]]  # coordinates
    combined_array[file1.shape[0]:, file2.shape[1]] = file2[:, file2.shape[1]]  # weights
    combined_array[file1.shape[0]:, file2.shape[1]+1] = -1  # filler
    combined_array[file1.shape[0]:, file2.shape[1]+2] = -1  # filler
    combined_array[file1.shape[0]:, file2.shape[1]+3] = 1  # indicating after
    np.savetxt('total_weight_on_each_ball_'+str(first_file_num)+'_'+str(second_file_num)+'.txt', combined_array, fmt='%1.10e')


if __name__ == '__main__':
    os.chdir('CAS_60_sc1_1')
    combine_files(133, 134)  # first file number is when clustering happens, second file number is right after clustering
