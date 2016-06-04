import numpy as np
import sys
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, Isomap


def calculate_distance(point1, point2):
    distance = 0.0
    for i in range(point1.shape[0]):
        distance += (point1[i]-point2[i])**2
    return np.sqrt(distance)


def calculate_angular_distance(point1, point2):
    distance = 0.0
    for i in range(point1.shape[0]):
        distance += min(360.0-abs(point1[i]-point2[i]), abs(point1[i]-point2[i]))**2
    return np.sqrt(distance)


def convert_matrix_to_angular(matrix):
    new_matrix = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i in range(new_matrix.shape[0]):
        for j in range(new_matrix.shape[1]):
            new_matrix[i, j] = calculate_angular_distance(matrix[i, :], matrix[j, :])
    np.savetxt('matrix.txt', new_matrix, fmt=' %1.10e')
    return new_matrix


def convert_angles_to_cos_sin(balls):
    new_balls = np.zeros((balls.shape[0], 2*balls.shape[1]))
    for i in range(new_balls.shape[0]):
        for j in range(new_balls.shape[1]):
            new_balls[i, 2*j] = np.cos(balls[i, j]*np.pi/180)
            new_balls[i, 2*j+1] = np.sin(balls[i, j]*np.pi/180)
    return new_balls


def principal_component_analysis(file_name, dimension):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    pca = PCA(n_components=2)
    transformed_matrix = pca.fit_transform(matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])


def multidimensioanl_scaling(file_name, dimension):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_matrix_to_angular(matrix)
    mds = MDS(n_components=2, metric=True, n_init=4, max_iter=300, verbose=0, eps=1e-6, n_jobs=1, random_state=None,
              dissimilarity='precomputed')
    #transformed_matrix = mds.fit_transform(new_matrix)
    transformed_matrix = mds.fit(new_matrix).embedding_
    # Rescale the data
    #transformed_matrix *= np.sqrt((matrix ** 2).sum()) / np.sqrt((transformed_matrix ** 2).sum())
    # Rotate the data
    #clf = PCA(n_components=2)
    #transformed_matrix = clf.fit_transform(transformed_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])


def isomap(file_name, dimension, num_neighbors):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    convert_matrix_to_angular(matrix)
    new_matrix = convert_angles_to_cos_sin(matrix)
    imap = Isomap(n_neighbors=num_neighbors, n_components=2, eigen_solver='auto', tol=0, max_iter=None,
                  path_method='auto', neighbors_algorithm='auto')
    transformed_matrix = imap.fit_transform(new_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    np.savetxt('labels.txt', balls[:, dimension], fmt=' %1.10e')


def diffusion_map(file_name, dimension, epsilon):
    # Given N data points x_n, n=1,...,N where each x_n\in R^{p},
    # the distance (similarity) between any two points is given by the formula
    # L_{i,j}=K(x_i,x_j)=exp(-||x_i-x_j||^2/(2*epsilon)) with Gaussian kernel of width epsilon
    # and a diagonal normalization matrix D with D_i=\sum_{j=1}^{N}L_{i,j}.
    # Then we solve the normalized eigenvalue problem L*v=lambda*D*v or M*v=lambda*v where M=D^{-1}L
    # and use the first few eigenvectors of M for low-dimensional representation of data.

    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    convert_matrix_to_angular(matrix)
    kernel_matrix = np.zeros([balls.shape[0], balls.shape[0]])
    diag_normalization_matrix = np.zeros([balls.shape[0], balls.shape[0]])
    for i in range(kernel_matrix.shape[0]):
        for j in range(kernel_matrix.shape[1]):
            kernel_matrix[i, j] = np.exp(-calculate_angular_distance(balls[i][0:dimension], balls[j][0:dimension])**2/epsilon)
    for i in range(diag_normalization_matrix.shape[0]):
        diag_normalization_matrix[i, i] = np.sum(kernel_matrix[i, :])  # normalizes rows of kernel matrix
    diag_normalization_matrix_inv = np.linalg.inv(diag_normalization_matrix)
    new_kernel_matrix = np.dot(np.dot(diag_normalization_matrix_inv, kernel_matrix), diag_normalization_matrix_inv)  # L^{(\alpha)}=D^{-\alpha}LD^{-\alpha}
    new_diag_normalization_matrix = np.zeros([balls.shape[0], balls.shape[0]])
    for i in range(new_diag_normalization_matrix.shape[0]):
        new_diag_normalization_matrix[i, i] = np.sum(new_kernel_matrix[i, :])  # normalizes rows of new kernel matrix, D^{(\alpha)}
    evals, evecs = np.linalg.eig(np.dot(np.linalg.inv(new_diag_normalization_matrix), new_kernel_matrix))  # M = (D^{(\alpha)})^{-1}L^{(\alpha)}
    col = np.argsort(-evals)  # sort eigenvalues in descending order
    first_eval = evals[col[0]]
    second_eval = evals[col[1]]
    first_evec = evecs[:, col[0]]
    second_evec = evecs[:, col[1]]

    ball_coords = np.zeros([balls.shape[0], dimension+3])
    for i in range(ball_coords.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension] = first_eval*first_evec[i]
        ball_coords[i, dimension+1] = second_eval*second_evec[i]
        ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    np.savetxt('labels.txt', balls[:, dimension], fmt=' %1.10e')


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print >>sys.stderr, "Need 4 args: file name, original dimension, # of neighbors (iso) / epsilon (diff), method (pca, mds, iso, diff)"
        sys.exit(1)
    file_name = sys.argv[1]
    dimension = int(sys.argv[2])
    if sys.argv[4] == 'pca':
        principal_component_analysis(file_name, dimension)
    elif sys.argv[4] == 'mds':
        multidimensioanl_scaling(file_name, dimension)
    elif sys.argv[4] == 'iso':
        num_neighbors = int(sys.argv[3])
        isomap(file_name, dimension, num_neighbors)
    elif sys.argv[4] == 'diff':
        epsilon = float(sys.argv[3])
        diffusion_map(file_name, dimension, epsilon)
