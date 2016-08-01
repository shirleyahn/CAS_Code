import numpy as np


def check_state_function(coordinates):
    if -100.0 <= float(coordinates[0]) <= -30.0 and -90.0 <= float(coordinates[1]) <= -10.0:
        return 1
    elif -180.0 <= float(coordinates[0]) <= -55.0 and (105.0 <= float(coordinates[1]) <= 180.0 or -180.0 <= float(coordinates[1]) <= -155.0):
        return 2
    else:
        return 0


def check_folded_unfolded(file_name, num_lines, num_cvs, num_cvs_for_conformation):
    input_file = open(file_name, 'r')
    output_file = open('result.txt', 'w')
    folded_unfolded = np.zeros(num_cvs/num_cvs_for_conformation, dtype='int64')
    result = np.zeros((1, num_cvs/num_cvs_for_conformation), dtype='int64')
    for i in range(num_lines):
        values = input_file.readline()
        values = (values.strip()).split()
        for j in range(num_cvs/num_cvs_for_conformation):
            folded_unfolded[j] = check_state_function(values[j*num_cvs_for_conformation:j*num_cvs_for_conformation+2])
        for item in folded_unfolded.tolist():
            output_file.write("%s " % item) 
        output_file.write("\n")
        if i == 0:
            result = np.sort(folded_unfolded.reshape(1, num_cvs/num_cvs_for_conformation))
        else:
            result = np.append(result, np.sort(folded_unfolded.reshape(1, num_cvs/num_cvs_for_conformation)), axis=0)
    input_file.close()
    output_file.close()

    result = result.tolist()
    setOfItems = []
    newListOfItems = []
    for item in result:
        if item in setOfItems:
            continue
        setOfItems.append(item)
        temp = list(item)
        occurence = result.count(item)
        temp.append(occurence)
        newListOfItems.append(temp)

    for i in range(len(newListOfItems)):
        print newListOfItems[i]

if __name__ == '__main__':
    check_folded_unfolded('CAS_60_sc2/total_weight_on_each_ball_170.txt', 520, 6, 2)
