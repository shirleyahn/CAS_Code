import numpy as np


def check_state_function(coordinates):
    if -100.0 <= float(coordinates[0]) <= -30.0 and -90.0 <= float(coordinates[1]) <= -10.0:
        return 1
    elif -180.0 <= float(coordinates[0]) <= -55.0 and (105.0 <= float(coordinates[1]) <= 180.0 or -180.0 <= float(coordinates[1]) <= -155.0):
        return 2
    else:
        return 0


def check_folded_unfolded(file_name, num_lines, num_cvs, num_cvs_for_conformation, num_clusters):
    input_file = open(file_name, 'r')
    output_file = open('folded_unfolded_result.txt', 'w')
    folded_unfolded = np.zeros(num_cvs/num_cvs_for_conformation, dtype='int64')
    result = np.zeros((num_lines, num_cvs/num_cvs_for_conformation, num_clusters), dtype='int64')

    for i in range(num_lines):
        values = input_file.readline()
        values = (values.strip()).split()
        for j in range(num_cvs/num_cvs_for_conformation):
            folded_unfolded[j] = check_state_function(values[j*num_cvs_for_conformation:j*num_cvs_for_conformation+2])
        for item in folded_unfolded.tolist():
            output_file.write("%s" % item)
        cluster_num = int(values[num_cvs])
        output_file.write(' ' + str(cluster_num))
        output_file.write('\n')
        # NOTE: get rid of sort if order is important
        result[i, :, cluster_num] = np.sort(folded_unfolded.reshape(1, num_cvs/num_cvs_for_conformation))
    input_file.close()
    output_file.close()

    for i in range(num_clusters):
        cluster_result = result[:,:,i].tolist()
        setOfItems = []
        newListOfItems = []
        for item in cluster_result:
            if item in setOfItems:
                continue
            setOfItems.append(item)
            temp = list(item)
            occurence = cluster_result.count(item)
            temp.append(occurence)
            newListOfItems.append(temp)

        print 'cluster number: ' + str(i)
        for j in range(len(newListOfItems)):
            print newListOfItems[j]

if __name__ == '__main__':
    check_folded_unfolded('CAS_60_sc/ball_clustering_170.txt', 520, 6, 2, 5)
