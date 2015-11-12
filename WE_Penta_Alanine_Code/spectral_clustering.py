#! /usr/bin/env python

import numpy as np
import copy
from scipy.cluster.vq import kmeans2, ClusterError
import sys

def delta2(c1, c2):
  minDist = np.inf
  for i in xrange(0, len(c1)):
    for j in xrange(0, len(c2)):
      p1 = c1[i,:]
      p2 = c2[j,:]
      dist = np.sqrt(np.sum(np.square(p2 - p1)))
      if dist < minDist:
        minDist = dist
  return minDist

def delta1(c):
  maxDist = 0
  for i in xrange(0, len(c)):
    for j in xrange(0, len(c)):
      if i == j:
        continue
      p1 = c[i,:]
      p2 = c[j,:]
      dist = np.sqrt(np.sum(np.square(p2 - p1)))
      if dist > maxDist:
        maxDist = dist
  return maxDist

def minDelta2(ball_coords):
  column = ball_coords.shape[1]-1
  num_clusters = int(np.max(ball_coords[:,column])+1)
  min_delta2 = np.inf
  for i in xrange(0,num_clusters):
    for j in xrange(0,num_clusters):
      if i == j:
        continue
      i = float(i)
      j = float(j)
      c1 = ball_coords[ball_coords[:,column] == i,:-1]
      c2 = ball_coords[ball_coords[:,column] == j,:-1]
      d2 = delta2(c1, c2)
      if d2 < min_delta2:
        min_delta2 = d2
  return min_delta2

def maxDelta1(ball_coords):
  column = ball_coords.shape[1]-1
  num_clusters = int(np.max(ball_coords[:,column])+1)
  max_delta1 = 0
  for i in xrange(0,num_clusters):
    i = float(i)
    c1 = ball_coords[ball_coords[:,column] == i,:-1]
    d1 = delta1(c1)
    if d1 > max_delta1:
      max_delta1 = d1
  return max_delta1

def dunn(ball_coords):
  return minDelta2(ball_coords) / maxDelta1(ball_coords)

def spectral_clustering(balls_file, evectors_file, num_clusters):
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
    a = balls
    b = log_second_evector
    matrix = np.hstack((balls, log_second_evector))

    while True:
        try:
            centroids, labels = kmeans2(matrix, num_clusters, minit='points', iter=200, missing='raise')
            break
        except ClusterError:
            num_clusters -= 1

    ball_coords = np.zeros((balls.shape[0], 10))
    for j in xrange(balls.shape[0]):
        ball_coords[j,0:6] = balls[j, 0:6].tolist()
        ball_coords[j, 6] = abs(evectors[j, 0])
        ball_coords[j, 7] = log_second_evector[j, 0]
        ball_coords[j, 8] = evectors[j, 2]
        ball_coords[j, 9] = labels[j]
        print ' '.join([str(x) for x in ball_coords[j, :]])
    dunnIndex = dunn(ball_coords)
    print >>sys.stderr, "Dunn index: %f" % dunnIndex

if __name__ == "__main__":
  if len(sys.argv) != 4:
    print >>sys.stderr, "Need 3 args: balls file, evectors file, number of clusters"
    sys.exit(1)
  balls_file = sys.argv[1]
  evectors_file = sys.argv[2]
  num_clusters = int(sys.argv[3])
  spectral_clustering(balls_file, evectors_file, num_clusters)
