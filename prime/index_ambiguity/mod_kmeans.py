from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
'''
Author      : Uervirojnangkoorn, M.
Created     : 8/17/2016
Description : Handy class to do K-means clustering.
'''
import random
import numpy as np
import matplotlib.pyplot as plt
from random import randint

class kmeans_handler(object):
  """
  From a given N-dimensional array, find K clusters.
  Algorithms taken from Stanford CS221 pseudo-python code.
  See http://stanford.edu/~cpiech/cs221/handouts/kmeans.html
  """
  def __init__(self):
    """
    Constructor
    """
    self.MAX_ITERS = 1000

  def run(self, dataset, k, flag_plot=False):
    #define color for plotting
    colors = ['#%06X' % randint(0, 0xFFFFFF) for i in range(k)]
    #initialize controids randomly
    n_data, n_features = dataset.shape
    centroids = np.array([random.random() for i in range(k*n_features)]).reshape((k,n_features))
    scale = np.max(dataset, axis=0)
    centroids = centroids * scale
    if flag_plot:
      plt.scatter(dataset[:,0], dataset[:,1], s=10, marker='x', c='b')
      plt.scatter(centroids[:,0], centroids[:,1], s=20, marker='o', c='k')
      plt.title('Initial step (k=%3.0f)'%(k))
      plt.show()
    #initialize book keeping vars
    n_iters = 0
    old_centroids = None
    #run the k-means algorithm
    while not self.should_stop(old_centroids, centroids, n_iters):
      #save old centroids for convergenence test. book keeping.
      old_centroids = centroids[:]
      n_iters += 1
      #assign labels to each datapoint based on centroids
      labels = self.get_labels(dataset, centroids)
      #assign centroids based on datapoint labels
      centroids = self.get_centroids(dataset, labels, k)
      if flag_plot:
        for i in range(k):
          s = dataset[labels==i]
          plt.scatter(s[:,0], s[:,1], s=10, marker='x', c=colors[i])
          plt.scatter(centroids[i,0], centroids[i,1], s=20, marker='o', c='k')
        plt.title('Cycle %3.0f (k=%3.0f)'%(n_iters, k))
        plt.show()
    return centroids, labels

  def should_stop(self, old_centroids, centroids, n_iters):
    if n_iters > self.MAX_ITERS: return True
    return np.all(old_centroids == centroids)

  def get_labels(self, dataset, centroids):
    #For each element in the dataset, choose the closest centroid.
    #Make that centroid the element's label.
    n_data, n_features = dataset.shape
    n_centroids = len(centroids)
    dist = np.zeros((n_data, n_centroids))
    for c,i in zip(centroids, range(n_centroids)):
      dist[:,i] = np.sqrt(((dataset-c)**2).sum(axis=1))
    labels = np.argmin(dist, axis=1)
    return labels

  def get_centroids(self, dataset, labels, k):
    #For each label group, calculate a new geometric mean as a new centroid.
    #For empty group, assign a new random centroid.
    n_data, n_features = dataset.shape
    centroids = np.zeros((k, n_features))
    for i in range(k):
      s = dataset[labels==i]
      if len(s) > 0:
        centroids[i] = np.power(np.sqrt(np.prod(s, axis=0)**2), 1/len(s))
      else:
        scale = np.max(dataset, axis=0)
        centroids[i] = np.array([random.random() for i in range(n_features)]) * scale
    return centroids

if __name__=="__main__":
  #unit test
  points = np.vstack([np.random.multivariate_normal(mean, \
                 0.03 * np.diag([1,1]), 20) \
                 for mean in [(1, 1), (2, 4), (3, 2)]])
  k = 3
  kmh = kmeans_handler()
  centroids, labels = kmh.run(points, k, flag_plot=True)
