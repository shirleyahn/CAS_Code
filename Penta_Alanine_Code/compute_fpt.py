import numpy as np
from scipy import stats
import os


def check_state_function(coordinates):
    if -100.0 <= coordinates[0] <= -30.0 and -90.0 <= coordinates[1] <= -10.0 and -100.0 <= coordinates[2] <= -30.0 and \
       -90.0 <= coordinates[3] <= -10.0 and -100.0 <= coordinates[4] <= -30.0 and -90.0 <= coordinates[5] <= -10.0:
        return 0
    elif -180.0 <= coordinates[0] <= -55.0 and (105.0 <= coordinates[1] <= 180.0 or -180.0 <= coordinates[1] <= -155.0) and \
         -180.0 <= coordinates[2] <= -55.0 and (105.0 <= coordinates[3] <= 180.0 or -180.0 <= coordinates[3] <= -155.0) and \
         -180.0 <= coordinates[4] <= -55.0 and (105.0 <= coordinates[5] <= 180.0 or -180.0 <= coordinates[5] <= -155.0):
        return 1
    else:
        return -1


def compute_fpt(del_t, file_directory, file_name, num_lines):
    os.chdir(file_directory)
    traj_file = open(file_name, 'r')
    data = np.zeros(num_lines, int)
    for i in range(num_lines):
        line = traj_file.readline().strip().split()
        coordinates = [float(coord) for coord in line]
        data[i] = int(check_state_function(coordinates[2:]))
    traj_file.close()

    output_file = open('mfpt.txt', 'w')
    conformation = -1
    prev_time = 0.0
    total_min1_time = 0.0 
    total_min2_time = 0.0 
    num_min1_trans = 0 
    num_min2_trans = 0 
    max_min1_time = 0.0 
    min_min1_time = 10**10
    max_min2_time = 0.0 
    min_min2_time = 10**10
    for i in range(num_lines):
        if data[i] != -1 and conformation == -1:  # very first transition from "none" to either min1 or min2
            conformation = data[i]
            prev_time = i*del_t
        elif data[i] != -1 and conformation != -1 and data[i] != conformation:  # the rest of the transitions (min1 to min2 or min2 to min1)
            if conformation == 0:
                num_min1_trans += 1
                time = i*del_t - prev_time
                total_min1_time += time
                if time > max_min1_time:
                    max_min1_time = time
                if time < min_min1_time: 
                    min_min1_time = time
                output_file.write("min1 -> min2 " + str(time) + ' ' + str(i*del_t) + "\n")
            else:
                num_min2_trans += 1
                time = i*del_t - prev_time
                total_min2_time += time
                if time > max_min2_time: 
                    max_min2_time = time
                if time < min_min2_time: 
                    min_min2_time = time
                output_file.write("min2 -> min1 " + str(time) + ' ' + str(i*del_t) + "\n")
            conformation = data[i]
            prev_time = i*del_t

    avg_min1_time = 0.0 
    if num_min1_trans != 0:
        avg_min1_time = total_min1_time/num_min1_trans
    avg_min2_time = 0.0 
    if num_min2_trans != 0: 
        avg_min2_time = total_min2_time/num_min2_trans
    ratio = 0.0 
    if avg_min1_time != 0.0 and avg_min2_time != 0.0: 
        ratio = avg_min2_time/avg_min1_time
    if num_min1_trans != 0:
        output_file.write("min min1 -> min2 " + str(min_min1_time) + "\n")
        output_file.write("max min1 -> min2 " + str(max_min1_time) + "\n")
        output_file.write("avg min1 -> min2 " + str(avg_min1_time) + "\n")
        output_file.write("avg flux min1 -> min2 " + str(1.0/avg_min1_time) + "\n")
        output_file.write("# of min1 -> min2 " + str(num_min1_trans) + "\n")
    if num_min2_trans != 0:
        output_file.write("min min2 -> min1 " + str(min_min2_time) + "\n")
        output_file.write("max min2 -> min1 " + str(max_min2_time) + "\n")
        output_file.write("avg min2 -> min1 " + str(avg_min2_time) + "\n")
        output_file.write("avg flux min2 -> min1 " + str(1.0/avg_min2_time) + "\n")
        output_file.write("# of min2 -> min1 " + str(num_min2_trans) + "\n")
    if num_min1_trans != 0 and num_min2_trans != 0: 
        output_file.write("ratio " + str(ratio) + "\n")
    output_file.close()
    
    """
    os.chdir("..")
    root_directory = os.getcwd()
    output_file = open('mfpt_total.txt','w')
    total_min1_time = 0.0
    total_min2_time = 0.0
    num_min1_trans = 0
    num_min2_trans = 0
    max_min1_time = 0.0
    min_min1_time = 10**10
    max_min2_time = 0.0
    min_min2_time = 10**10
    first_mfpt_array = []
    second_mfpt_array = []
    for subdir, dirs, files in os.walk(root_directory):
        for file in files:
            if file == 'mfpt.txt':
                os.chdir(subdir)
                input_file = open(file, 'r')
                while True:
                    line = input_file.readline()
                    line = line.split()
                    if not line:
                        break
                    if line[0] == "min1":
                        time = float(line[3])
                        current_time = float(line[4])
                        num_min1_trans += 1
                        total_min1_time += time
                        if time > max_min1_time:
                            max_min1_time = time
                        if time < min_min1_time:
                            min_min1_time = time
                        output_file.write("min1 -> min2 " + str(time) + ' ' + str(current_time) + "\n")
                        first_mfpt_array.append(float(time))
                    if line[0] == "min2":
                        time = float(line[3])
                        current_time = float(line[4])
                        num_min2_trans += 1
                        total_min2_time += time
                        if time > max_min2_time:
                            max_min2_time = time
                        if time < min_min2_time:
                            min_min2_time = time
                        output_file.write("min2 -> min1 " + str(time) + ' ' + str(current_time) + "\n")
                        second_mfpt_array.append(float(time))
                input_file.close()
    avg_min1_time = 0.0 
    if num_min1_trans != 0:
        avg_min1_time = total_min1_time/num_min1_trans
    avg_min2_time = 0.0 
    if num_min2_trans != 0: 
        avg_min2_time = total_min2_time/num_min2_trans
    ratio = 0.0 
    if avg_min1_time != 0.0 and avg_min2_time != 0.0: 
        ratio = avg_min2_time/avg_min1_time
    if num_min1_trans != 0:
        output_file.write("min min1 -> min2 " + str(min_min1_time) + "\n")
        output_file.write("max min1 -> min2 " + str(max_min1_time) + "\n")
        output_file.write("avg min1 -> min2 " + str(avg_min1_time) + "\n")
        output_file.write("avg flux min1 -> min2 " + str(1.0/avg_min1_time) + "\n")
        output_file.write("std flux min1 -> min2 " + str(1.0/np.std(np.asarray(first_mfpt_array))) + "\n")
        output_file.write("# of min1 -> min2 " + str(num_min1_trans) + "\n")
    if num_min2_trans != 0:
        output_file.write("min min2 -> min1 " + str(min_min2_time) + "\n")
        output_file.write("max min2 -> min1 " + str(max_min2_time) + "\n")
        output_file.write("avg min2 -> min1 " + str(avg_min2_time) + "\n")
        output_file.write("avg flux min2 -> min1 " + str(1.0/avg_min2_time) + "\n")
        output_file.write("std flux min1 -> min2 " + str(1.0/np.std(np.asarray(second_mfpt_array))) + "\n")
        output_file.write("# of min2 -> min1 " + str(num_min2_trans) + "\n")
    if num_min1_trans != 0 and num_min2_trans != 0: 
        output_file.write("ratio " + str(ratio) + "\n")
    output_file.close() 
    """
if __name__ == '__main__':
    compute_fpt(0.001, 'RUN00', 'traj_angles.txt', 3000001)
