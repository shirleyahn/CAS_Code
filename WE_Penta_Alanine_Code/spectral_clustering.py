import numpy as np
import copy
import os
from scipy.cluster.vq import kmeans2, ClusterError


def spectral_clustering(balls_file, evectors_file, num_clusters, new_ball_clustering_file_name):
    evectors = np.loadtxt(evectors_file)
    balls = np.loadtxt(balls_file)

    second_evector = evectors[:, 1]
    second_evector = second_evector.reshape(second_evector.shape[0], 1)
    log_second_evector = np.zeros((second_evector.shape[0], 1))
    for i in range(second_evector.shape[0]):
        if second_evector[i] < 0.0:
            log_second_evector[i] = -np.log(-second_evector[i])
        elif second_evector[i] == 0.0 or second_evector[i] == 1.0:
            log_second_evector[i] = 0.0
        else:
            log_second_evector[i] = np.log(second_evector[i])

    '''
    sorted_second_evector = np.sort(second_evector, axis=0)
    second_evector_order = np.ndarray.argsort(second_evector)
    num_balls = int(np.ceil(len(sorted_second_evector) / num_clusters))
    array_of_clusters = [sorted_second_evector[i:i + num_balls] for i in
                         range(0, len(sorted_second_evector), num_balls)]
    array_of_orderings = [second_evector_order[i:i + num_balls] for i in range(0, len(second_evector_order), num_balls)]
    num_clusters = len(array_of_clusters)
    '''

    while True:
        try:
            centroids, labels = kmeans2(log_second_evector, num_clusters, minit='points', iter=30, missing='raise')
            break
        except ClusterError:
            num_clusters -= 1

    f = open(new_ball_clustering_file_name, 'w')

    '''
    for i in range(num_clusters):
        first = 0
        cluster = array_of_clusters[i]
        ordering = array_of_orderings[i]
        for j in range(cluster.shape[0]):
            if first == 0:
                first += 1
                ref_ball_center = balls[ordering[j], 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ref_ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[ordering[j], 0]))
                ball_cluster.append(second_evector[ordering[j]])
                ball_cluster.append(final_evectors[ordering[j], 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
            else:
                ball_center = balls[ordering[j], 0:gv.num_cvs].tolist()
                ball_cluster = copy.deepcopy(ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(final_evectors[ordering[j], 0]))
                ball_cluster.append(second_evector[ordering[j]])
                ball_cluster.append(final_evectors[ordering[j], 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
    '''

    for i in range(num_clusters):
        first = 0
        for j in range(balls.shape[0]):
            if labels[j] == i and first == 0:
                first += 1
                ref_ball_center = balls[j, 0:6].tolist()
                ball_cluster = copy.deepcopy(ref_ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(evectors[j, 0]))
                ball_cluster.append(log_second_evector[j, 0])
                ball_cluster.append(evectors[j, 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                #ball_clusters_list[tuple(ref_ball_center)] = [tuple(ref_ball_center)]
                #balls[j][gv.num_cvs+2] -= 1
            elif labels[j] == i and first != 0:
                ball_center = balls[j, 0:6].tolist()
                ball_cluster = copy.deepcopy(ball_center)
                ball_cluster.append(i)
                ball_cluster.append(abs(evectors[j, 0]))
                ball_cluster.append(log_second_evector[j, 0])
                ball_cluster.append(evectors[j, 2])
                f.write(' '.join(map(lambda coordinate: str(coordinate), ball_cluster)))
                f.write('\n')
                #ball_clusters_list[tuple(ref_ball_center)].append(tuple(ball_center))
                #balls[j][gv.num_cvs+2] -= 1
    f.close()

    #np.savetxt('evalues_' + str(step_num + 1) + '.txt', final_evalues, fmt=' %1.10e')
    #np.savetxt('evectors_' + str(step_num + 1) + '.txt', final_evectors, fmt=' %1.10e')
    #np.savetxt('transition_matrix_' + str(step_num + 1) + '.txt', symmetric_transition_matrix, fmt=' %1.10e')

os.chdir('/Users/Ahn/Dropbox (Stanford Mechanics)/Hee Sun Shirley/WE_Enhanced_Sampling/Penta_Alanine/WE_v2')
spectral_clustering('balls_62.txt', 'evectors_62.txt', 5, 'ball_clustering_62_new.txt')
