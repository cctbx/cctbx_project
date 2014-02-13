from __future__ import division
import os
import re
import sys
from libtbx import easy_pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
#from scipy.cluster.vq import kmeans, kmeans2, vq
import  scipy.cluster.hierarchy as hcluster
#import libtbx.load_env
from cctbx import sgtbx

class target:
  
  def __init__(self, path_to_integration_dir):
    """ Creates a list of (point group, unit cell) tuples, and a list of niggli cells from the recursively walked 
        paths. Can take more than one argument for multiple folders."""
    self.all_uc     = []
    self.niggli_ucs = []
    for arg in path_to_integration_dir:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          path = os.path.join(dirpath, filename)
          try:
            # Warn on error, but continue directory traversal.
            d = easy_pickle.load(path)
            pg = d['pointgroup']
            uc = d['current_orientation'][0].unit_cell()
          except KeyError:
              sys.stderr.write(
               "Could not extract point group and unit cell from %s\n" % path)
          except Exception:
              sys.stderr.write(
                    "Could not read %s\n" % path)
          else:
            self.all_uc.append((pg, uc))
            self.niggli_ucs.append(uc.niggli_cell().parameters())
    
  def find_distance(self, G6a, G6b):
    """ Retursn the distance between two cells, already in the G6 convention. Curently trivial Euclidian """
    val = abs((np.sum(G6a - G6b)))**0.5
    return val
       
  def make_G6(self, uc):
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    a = uc[0]**2 
    b = uc[1]**2 
    c = uc[2]**2 
    d = 2*uc[1]*uc[2]*math.sin(3)
    e = 2*uc[0]*uc[2]*math.sin(4)
    f = 2*uc[0]*uc[1]*math.sin(5)
    return [a,b,c,d,e,f]
    
  def cluster(self, threshold):
    """ Do basic hierarchical clustering on the Niggli cells created at the start """
    # 1. Create a numpy array of G6 cells
    G6_cells = []
    for n_cell in self.niggli_ucs: 
      G6_cells.append(self.make_G6(n_cell))
    self.G6_cells = np.array(G6_cells)
    # 2. Do hierarchichal clustering on this, using the find_distance method above. 
    self.clusters = hcluster.fclusterdata(self.G6_cells, 
                                          threshold, 
                                          criterion='distance', 
                                          metric=lambda a, b: self.find_distance(a,b))
    # 3. print out some information that is useful.
    print "{} clusters have been identified.".format(max(self.clusters))
    print "{:^14} {:<11} {:<11} {:<11} {:<11} {:<11} {:<11}".format(
                             "Num in cluster", 
                             "Med_a", "Med_b", "Med_c", 
                             "Med_alpha", "Med_beta", "Med_gamma")
    cluster_idx = 1
    #import pdb; pdb.set_trace()
    for cluster in range(max(self.clusters)):
      this_cluster = np.array([self.niggli_ucs[i] for i in range(len(self.niggli_ucs)) if self.clusters[i]==cluster+1])
      print "{:<14} {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f}) {:<5.2f}({:<4.2f}) {:<5.2f}({:<4.2f}) {:<5.2f}({:<4.2f})".format(
        len(this_cluster),
        np.median(this_cluster[:,0]), np.std(this_cluster[:,0]),
        np.median(this_cluster[:,1]), np.std(this_cluster[:,1]),
        np.median(this_cluster[:,2]), np.std(this_cluster[:,2]),
        np.median(this_cluster[:,3]), np.std(this_cluster[:,3]),
        np.median(this_cluster[:,4]), np.std(this_cluster[:,4]),
        np.median(this_cluster[:,5]), np.std(this_cluster[:,5]))
    print "Standard deviations are in brackets."           

    
  def plot_clusters(self, clusters):
    """ Plot Niggli cells -- one plot for (a,b,c) and one plot for (alpha, beta, gamma) -- colour 
    coded by cluster index.  """
    colors = ['b', 'y', 'g', 'c', 'm', 'r', 'k']
    fig = plt.figure('unit_cells_dimensions')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('a [A]')
    ax.set_ylabel('b [A]')
    ax.set_zlabel('c [A]')
    for i in range(len(self.niggli_ucs)):
      ax.scatter(np.array(self.niggli_ucs)[i,0],
                 np.array(self.niggli_ucs)[i,1],
                 np.array(self.niggli_ucs)[i,2],
                 c=colors[clusters[i]], marker='o')
    fig = plt.figure('unit_cells_angles')
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(self.niggli_ucs)):
      ax.scatter(np.array(self.niggli_ucs)[i,3],
                 np.array(self.niggli_ucs)[i,4],
                 np.array(self.niggli_ucs)[i,5],
                 c=colors[clusters[i]], marker='o')
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    ax.set_zlabel('gamma')
    plt.show()

#    centroids,_ = kmeans(reduced_data[:,:3], num_clusters)
#    idx,_ = vq(reduced_data[:,:3],centroids)
