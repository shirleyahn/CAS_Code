import numpy as np
import sys
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, Isomap


def calculate_distance_from_center(point1, point2):
    distance = 0.
    for i in range(point1.shape[0]):
            if point1[i] - point2[i] > 180.0:
                distance += (point1[i] - point2[i] - 360.0) ** 2
            elif point1[i] - point2[i] < -180.0:
                distance += (point1[i] - point2[i] + 360.0) ** 2
            else:
                distance += (point1[i] - point2[i]) ** 2
    return np.sqrt(distance)


def compute_angular_distances(matrix):
    similarities = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i in range(similarities.shape[0]):
        for j in range(similarities.shape[1]):
            similarities[i, j] = calculate_distance_from_center(matrix[i, :], matrix[j, :])
    return similarities


def principal_component_analysis(file_name, dimension):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    pca = PCA(n_components=2)
    transformed_matrix = pca.fit_transform(matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for j in xrange(balls.shape[0]):
        ball_coords[j, 0:dimension] = balls[j, 0:dimension].tolist()
        ball_coords[j, dimension:dimension+2] = transformed_matrix[j]
        ball_coords[j, dimension+2] = balls[j, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[j, :]])


def multidimensioanl_scaling(file_name, dimension):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    similarities = compute_angular_distances(matrix)
    mds = MDS(n_components=2, metric=True, n_init=4, max_iter=300, verbose=0, eps=1e-6, n_jobs=1, random_state=None,
              dissimilarity='precomputed')
    #transformed_matrix = mds.fit_transform(similarities)
    transformed_matrix = mds.fit(similarities).embedding_
    # Rescale the data
    #transformed_matrix *= np.sqrt((matrix ** 2).sum()) / np.sqrt((transformed_matrix ** 2).sum())
    # Rotate the data
    #clf = PCA(n_components=2)
    #transformed_matrix = clf.fit_transform(transformed_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for j in xrange(balls.shape[0]):
        ball_coords[j, 0:dimension] = balls[j, 0:dimension].tolist()
        ball_coords[j, dimension:dimension+2] = transformed_matrix[j]
        ball_coords[j, dimension+2] = balls[j, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[j, :]])


def isomap(file_name, dimension):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    imap = Isomap(n_neighbors=10, n_components=2, eigen_solver='auto', tol=0, max_iter=None, path_method='auto',
                  neighbors_algorithm='auto')
    transformed_matrix = imap.fit_transform(matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for j in xrange(balls.shape[0]):
        ball_coords[j, 0:dimension] = balls[j, 0:dimension].tolist()
        ball_coords[j, dimension:dimension+2] = transformed_matrix[j]
        ball_coords[j, dimension+2] = balls[j, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[j, :]])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print >>sys.stderr, "Need 2 args: file name, dimension"
        sys.exit(1)
    file_name = sys.argv[1]
    dimension = int(sys.argv[2])
    #principal_component_analysis(file_name, dimension)
    multidimensioanl_scaling(file_name, dimension)
    #isomap(file_name, dimension)
