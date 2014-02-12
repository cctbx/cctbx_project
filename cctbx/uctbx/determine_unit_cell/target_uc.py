from __future__ import division
import os
import re
import sys
from libtbx import easy_pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.cluster.vq import kmeans, kmeans2, vq
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
    a = G6a[0]**2 - G6b[0]**2
    b = G6a[1]**2 - G6b[1]**2
    c = G6a[2]**2 - G6b[2]**2
    d = G6a[3]**2 - G6b[3]**2
    e = G6a[4]**2 - G6b[4]**2
    f = G6a[5]**2 - G6b[5]**2
    return sqrt(a + b + c + d + e + f)

  def make_G6(self, uc):
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    

#  def cluster(num_clusters)
#    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
#    reduced_data = np.array(data)
#    # These next lines, where I actually cluster, are curently b******. 
#    #Need to put something meaningfull here
#    centroids,_ = kmeans(reduced_data[:,:3], num_clusters)
#    idx,_ = vq(reduced_data[:,:3],centroids)
#    
#    fig = plt.figure('unit_cells_dimensions')
#    ax = fig.add_subplot(111, projection='3d')
#    ax.set_xlabel('a [um]')
#    ax.set_ylabel('b [um]')
#    ax.set_zlabel('c [um]')
#    ax.scatter(centroids[:,0],centroids[:,1],centroids[:,2],'k')
#    for i in range(num_clusters):
#      ax.scatter(reduced_data[idx==i,0],reduced_data[idx==i,1],reduced_data[idx==i,2],c=colors[i], marker='o')
#    fig = plt.figure('unit_cells_angles')
#    ax = fig.add_subplot(111, projection='3d')
#    for i in range(num_clusters):
#      ax.scatter(reduced_data[idx==i,3],reduced_data[idx==i,4],reduced_data[idx==i,5],c=colors[i], marker='o')
#    
#    ax.set_xlabel('alpha')
#    ax.set_ylabel('beta')
#    ax.set_zlabel('gamma')
#    
#    
#    plt.show()
