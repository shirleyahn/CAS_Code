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
    for i in range(fluxes_file.shape[0]):
        values = paths_file.readline()
        values = (values.strip()).split()
        key = ''
        for j in range(len(values)):
            for k in range(num_cvs/num_cvs_for_conformation):
                folded_unfolded[k] = check_state_function(conformations_file[values[j], k*num_cvs_for_conformation:k*num_cvs_for_conformation+2])
                key += str(folded_unfolded[k])
            for item in folded_unfolded.tolist():
                output_file.write("%s" % item)
            key += ' '
            output_file.write(" ")
        output_file.write("\n")
        if key in pathways_dict:
            pathways_dict[key] += fluxes_file[i]
        else:
            pathways_dict[key] = fluxes_file[i]
    paths_file.close()
    output_file.close()
    for i in sorted(pathways_dict, key=pathways_dict.get, reverse=True):
        print i, pathways_dict[i]


if __name__ == '__main__':
    check_pathways('CAS_60_sc1_1/total_weight_on_each_ball_133.txt', 'CAS_60_sc1_1/top_paths.txt', 'CAS_60_sc1_1/relative_fluxes.txt', 20360, 6, 2)
