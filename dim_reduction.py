import numpy as np
import sys
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, Isomap, SpectralEmbedding


def calculate_distance(point1, point2):
    distance = 0.0
    for i in range(point1.shape[0]):
        distance += (point1[i]-point2[i])**2
    return np.sqrt(distance)


def calculate_distance_matrix(matrix):
    distance_matrix = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i in range(distance_matrix.shape[0]):
        for j in range(distance_matrix.shape[1]):
            distance_matrix[i, j] = calculate_distance(matrix[i, :], matrix[j, :])
    return distance_matrix


def calculate_angular_distance(point1, point2):
    distance = 0.0
    for i in range(point1.shape[0]):
        distance += min(360.0-abs(point1[i]-point2[i]), abs(point1[i]-point2[i]))**2
    return np.sqrt(distance)


def calculate_angular_distance_matrix(matrix):
    angular_distance_matrix = np.zeros((matrix.shape[0], matrix.shape[0]))
    for i in range(angular_distance_matrix.shape[0]):
        for j in range(angular_distance_matrix.shape[1]):
            angular_distance_matrix[i, j] = calculate_angular_distance(matrix[i, :], matrix[j, :])
    return angular_distance_matrix


def convert_angles_to_cos_sin(balls):
    new_balls = np.zeros((balls.shape[0], 2*balls.shape[1]))
    for i in range(balls.shape[0]):
        for j in range(balls.shape[1]):
            new_balls[i, 2*j] = np.cos(balls[i, j]*np.pi/180)
            new_balls[i, 2*j+1] = np.sin(balls[i, j]*np.pi/180)
    return new_balls


def principal_component_analysis(file_name, dimension, label):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_angles_to_cos_sin(matrix)
    pca = PCA(n_components=2)
    transformed_matrix = pca.fit_transform(new_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        if label == 'cluster':
            ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        elif label == 'eq':
            ball_coords[i, dimension+2] = (-0.0019872041*300*np.log(abs(balls[i, dimension+1]))).tolist()
        elif label == 'committor':
            ball_coords[i, dimension+2] = (balls[i, dimension+2]/abs(balls[i, dimension+1])).tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    #np.savetxt('labels.txt', ball_coords[:, dimension+2], fmt=' %1.10e')


def multidimensioanl_scaling(file_name, dimension, label):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_angles_to_cos_sin(matrix)
    mds = MDS(n_components=2, metric=True, n_init=4, max_iter=300, verbose=0, eps=1e-6, n_jobs=1, random_state=None,
              dissimilarity='euclidean')
    transformed_matrix = mds.fit_transform(new_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        if label == 'cluster':
            ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        elif label == 'eq':
            ball_coords[i, dimension+2] = (-0.0019872041*300*np.log(abs(balls[i, dimension+1]))).tolist()
        elif label == 'committor':
            ball_coords[i, dimension+2] = (balls[i, dimension+2]/abs(balls[i, dimension+1])).tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    #np.savetxt('labels.txt', ball_coords[:, dimension+2], fmt=' %1.10e')


def isomap(file_name, dimension, num_neighbors, label):
    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_angles_to_cos_sin(matrix)
    imap = Isomap(n_neighbors=num_neighbors, n_components=2, eigen_solver='auto', tol=0, max_iter=None,
                  path_method='auto', neighbors_algorithm='auto')
    transformed_matrix = imap.fit_transform(new_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        if label == 'cluster':
            ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        elif label == 'eq':
            ball_coords[i, dimension+2] = (-0.0019872041*300*np.log(abs(balls[i, dimension+1]))).tolist()
        elif label == 'committor':
            ball_coords[i, dimension+2] = (balls[i, dimension+2]/abs(balls[i, dimension+1])).tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    #np.savetxt('labels.txt', ball_coords[:, dimension+2], fmt=' %1.10e')


def spectral_embedding(file_name, dimension, num_neighbors, label):
    # aka Laplacian eigenmaps
    # finds a low dimensional representation of the data using a spectral decomposition of the graph Laplacian
    # graph = discrete approximation of the low dimensional manifold in the high dimensional space
    # minimization of the cost function based on the graph preserves local distances

    # 1. weighted graph construction
    # 2. graph Laplacian construction: unnormalized: L = D-A, normalized: L=D^{-1/2}(D-A)D^{-1/2}
    # 3. partial eigenvalue decomposition

    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_angles_to_cos_sin(matrix)
    s_embedding = SpectralEmbedding(n_components=2, affinity='nearest_neighbors', gamma=None, random_state=None,
                                    eigen_solver=None, n_neighbors=num_neighbors)
    transformed_matrix = s_embedding.fit_transform(new_matrix)
    ball_coords = np.zeros((balls.shape[0], dimension+3))
    for i in xrange(balls.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension:dimension+2] = transformed_matrix[i]
        if label == 'cluster':
            ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        elif label == 'eq':
            ball_coords[i, dimension+2] = (-0.0019872041*300*np.log(abs(balls[i, dimension+1]))).tolist()
        elif label == 'committor':
            ball_coords[i, dimension+2] = (balls[i, dimension+2]/abs(balls[i, dimension+1])).tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    #np.savetxt('labels.txt', ball_coords[:, dimension+2], fmt=' %1.10e')


def diffusion_map(file_name, dimension, epsilon, label):
    # Given N data points x_n, n=1,...,N where each x_n\in R^{p},
    # the distance (similarity) between any two points is given by the formula
    # L_{i,j}=K(x_i,x_j)=exp(-||x_i-x_j||^2/epsilon) with Gaussian kernel of width epsilon
    # and a diagonal normalization matrix D with D_i=\sum_{j=1}^{N}L_{i,j}.
    # Then we solve the normalized eigenvalue problem L*v=lambda*D*v or M*v=lambda*v where M=D^{-1}L
    # and use the first few eigenvectors of M for low-dimensional representation of data.
    # Note: set the parameter alpha = 0, reducing it to classical graph Laplacian normalization

    balls = np.loadtxt(file_name)
    matrix = balls[:, 0:dimension]
    new_matrix = convert_angles_to_cos_sin(matrix)
    kernel_matrix = np.zeros([new_matrix.shape[0], new_matrix.shape[0]])
    for i in range(kernel_matrix.shape[0]):
        for j in range(kernel_matrix.shape[1]):
            kernel_matrix[i, j] = np.exp(-calculate_distance(new_matrix[i, :], new_matrix[j, :])**2/epsilon)
    v = np.diag(np.sum(kernel_matrix, axis=1))
    v_inv = np.linalg.inv(v)
    A = np.dot(v_inv, kernel_matrix)
    evals, evecs = np.linalg.eig(A)
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    first_eval = evals[1]
    second_eval = evals[2]
    first_evec = evecs[:, 1]
    second_evec = evecs[:, 2]
    ball_coords = np.zeros([balls.shape[0], dimension+3])
    for i in range(ball_coords.shape[0]):
        ball_coords[i, 0:dimension] = balls[i, 0:dimension].tolist()
        ball_coords[i, dimension] = first_eval*first_evec[i]
        ball_coords[i, dimension+1] = second_eval*second_evec[i]
        if label == 'cluster':
            ball_coords[i, dimension+2] = balls[i, dimension].tolist()
        elif label == 'eq':
            ball_coords[i, dimension+2] = (-0.0019872041*300*np.log(abs(balls[i, dimension+1]))).tolist()
        elif label == 'committor':
            ball_coords[i, dimension+2] = (balls[i, dimension+2]/abs(balls[i, dimension+1])).tolist()
        print ' '.join([str(x) for x in ball_coords[i, :]])
    #np.savetxt('labels.txt', ball_coords[:, dimension+2], fmt=' %1.10e')


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print >>sys.stderr, "Need 5 args: file name, original dimension, # of neighbors (iso, sembed) / epsilon (diff), method (pca, mds, iso, sembed, diff), label (cluster, eq, committor)"
        sys.exit(1)
    file_name = sys.argv[1]
    dimension = int(sys.argv[2])
    label = sys.argv[5]
    if sys.argv[4] == 'pca':
        principal_component_analysis(file_name, dimension, label)
    elif sys.argv[4] == 'mds':
        multidimensioanl_scaling(file_name, dimension, label)
    elif sys.argv[4] == 'iso':
        num_neighbors = int(sys.argv[3])
        isomap(file_name, dimension, num_neighbors, label)
    elif sys.argv[4] == 'sembed':
        num_neighbors = int(sys.argv[3])
        spectral_embedding(file_name, dimension, num_neighbors, label)
    elif sys.argv[4] == 'diff':
        epsilon = float(sys.argv[3])
        diffusion_map(file_name, dimension, epsilon, label)
