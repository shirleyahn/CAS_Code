import numpy as np


def check_state_function(coordinates):
    if -100.0 <= float(coordinates[0]) <= -30.0 and -90.0 <= float(coordinates[1]) <= -10.0:
        return 1
    elif -180.0 <= float(coordinates[0]) <= -55.0 and (105.0 <= float(coordinates[1]) <= 180.0 or -180.0 <= float(coordinates[1]) <= -155.0):
        return 2
    else:
        return 0


def check_pathways(conformations_file_name, paths_file_name, fluxes_file_name, num_lines, num_cvs, num_cvs_for_conformation):
    conformations_file = np.loadtxt(conformations_file_name)
    paths_file = open(paths_file_name, 'r')
    fluxes_file = np.loadtxt(fluxes_file_name)
    output_file = open('top_paths_v2.txt', 'w')
    folded_unfolded = np.zeros(num_cvs/num_cvs_for_conformation, dtype='int64')
    pathways_dict = {}
    pathways_dict2 = {}
    for i in range(fluxes_file.shape[0]):
        values = paths_file.readline()
        values = (values.strip()).split()
        key = ''
        key2 = ''
        for j in range(len(values)):
            for k in range(num_cvs/num_cvs_for_conformation):
                folded_unfolded[k] = check_state_function(conformations_file[values[j], k*num_cvs_for_conformation:k*num_cvs_for_conformation+2])
                if 0 < j < len(values)-1 and k == 0:
                    key += str(values[j])
                elif j == 0 and k == 0:
                    key += 'folded'
                elif j == len(values)-1 and k == 0:
                    key += 'unfolded'
                key2 += str(folded_unfolded[k])
            for item in folded_unfolded.tolist():
                output_file.write("%s" % item)
            key += ' '
            key2 += ' '
            output_file.write(" ")
        output_file.write("\n")
        if key in pathways_dict:
            pathways_dict[key] += fluxes_file[i]
        else:
            pathways_dict[key] = fluxes_file[i]
        if key2 in pathways_dict2:
            pathways_dict2[key2] += fluxes_file[i]
        else:
            pathways_dict2[key2] = fluxes_file[i]
    paths_file.close()
    output_file.close()
    output_file = open('top_paths_v3.txt', 'w')
    for i in sorted(pathways_dict, key=pathways_dict.get, reverse=True):
        output_file.write(str(i) + ' ' + str(pathways_dict[i]) + '\n')
    output_file.close()
    output_file = open('top_paths_v4.txt', 'w')
    for i in sorted(pathways_dict2, key=pathways_dict2.get, reverse=True):
        output_file.write(str(i) + ' ' + str(pathways_dict2[i]) + '\n')
    output_file.close()


if __name__ == '__main__':
    check_pathways('CAS_60_sc1_1/total_weight_on_each_ball_133.txt', 'CAS_60_sc1_1/top_paths.txt', 'CAS_60_sc1_1/relative_fluxes.txt', 20360, 6, 2)
