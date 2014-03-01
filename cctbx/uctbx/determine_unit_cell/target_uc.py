from __future__ import division
import os
import sys
from libtbx import easy_pickle
import numpy as np
import math
#from scipy.cluster.vq import kmeans, kmeans2, vq
#import libtbx.load_env
from cctbx.uctbx.determine_unit_cell import NCDist

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

  def find_distance(self, G6a, G6b, key):
    """ Retursn the distance between two cells, already in the G6 convention.
        Optionally (key =1)  trivial Euclidian or NCDist (key =2 )"""
    if key is 1: # trivial euclidian
      return abs((np.sum(G6a - G6b)))**0.5
    if key is 2: # Andrews-Bernstein distance
      return  NCDist(G6a , G6b)


  def make_G6(self, uc):
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    a = uc[0]**2
    b = uc[1]**2
    c = uc[2]**2
    d = 2*uc[1]*uc[2]*math.cos(uc[3])
    e = 2*uc[0]*uc[2]*math.cos(uc[4])
    f = 2*uc[0]*uc[1]*math.cos(uc[5])
    return [a,b,c,d,e,f]

  def cluster(self, threshold, method, linkage):
    """ Do basic hierarchical clustering on the Niggli
    cells created at the start """
    dist_method = 2
    if dist_method is 2:
      print ("Using Andrews-Bernstein Distance from " +
      "Andrews & Bernstein J Appl Cryst 47:346 (2014).")
    # 1. Create a numpy array of G6 cells
    G6_cells = []
    for n_cell in self.niggli_ucs:
      G6_cells.append(self.make_G6(n_cell))
    self.G6_cells = np.array(G6_cells)
    # 2. Do hierarchichal clustering on this, using the find_distance method above.
    import  scipy.cluster.hierarchy as hcluster
    self.clusters = hcluster.fclusterdata(self.G6_cells,
                                          threshold**2,
                                          criterion='distance',
                                          metric=lambda a, b: self.find_distance(a,b,dist_method))
    # 3. print out some information that is useful.
    print "{} clusters have been identified.".format(max(self.clusters))
    print "{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
                             "Num in cluster",
                             "Med_a", "Med_b", "Med_c",
                             "Med_alpha", "Med_beta", "Med_gamma")
    cluster_idx = 1
    #import pdb; pdb.set_trace()
    singletons = []
    for cluster in range(max(self.clusters)):
      this_cluster = np.array([self.niggli_ucs[i]
                               for i in range(len(self.niggli_ucs))
                               if self.clusters[i]==cluster+1])
      if len(this_cluster) is not 1:
        print "".join([("{:<14} {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f})" +
                        " {:<5.1f}({:<4.1f})").format(
          len(this_cluster),
          np.median(this_cluster[:,0]), np.std(this_cluster[:,0]),
          np.median(this_cluster[:,1]), np.std(this_cluster[:,1]),
          np.median(this_cluster[:,2]), np.std(this_cluster[:,2])),
          " {:<6.2f}({:<4.2f}) {:<6.2f}({:<4.2f}) {:<6.2f}({:<4.2f})".format(
          np.median(this_cluster[:,3]), np.std(this_cluster[:,3]),
          np.median(this_cluster[:,4]), np.std(this_cluster[:,4]),
          np.median(this_cluster[:,5]), np.std(this_cluster[:,5]))])
      else:
        singletons.append("".join([("{:<14} {:<11.1f} {:<11.1f}" +
                          " {:<11.1f}").format(
          '',
          np.median(this_cluster[:,0]),
          np.median(this_cluster[:,1]),
          np.median(this_cluster[:,2])),
          " {:<12.1f} {:<12.1f} {:<12.1f}".format(
          np.median(this_cluster[:,3]),
          np.median(this_cluster[:,4]),
          np.median(this_cluster[:,5])), '\n']))
    print "{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
                             "Num in cluster",
                             "Med_a", "Med_b", "Med_c",
                             "Med_alpha", "Med_beta", "Med_gamma")
    print "Standard deviations are in brackets."
    print  str(len(singletons)) + " singletons:  \n"
    print "".join(singletons)


  def plot_clusters(self, clusters):
    """ Plot Niggli cells -- one plot for (a,b,c) and one plot for (alpha, beta, gamma) -- colour
    coded by cluster index.  """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.cm as mpl_cmaps
    cmap = mpl_cmaps.Paired
    fig = plt.figure('unit_cells_dimensions')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('a [A]')
    ax.set_ylabel('b [A]')
    ax.set_zlabel('c [A]')
    for i in range(len(self.niggli_ucs)):
      ax.scatter(np.array(self.niggli_ucs)[i,0],
                 np.array(self.niggli_ucs)[i,1],
                 np.array(self.niggli_ucs)[i,2],
                 c=cmap(clusters[i-1]/max(clusters)), marker='o', s=20)
    fig = plt.figure('unit_cells_angles')
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(self.niggli_ucs)):
      ax.scatter(np.array(self.niggli_ucs)[i,3],
                 np.array(self.niggli_ucs)[i,4],
                 np.array(self.niggli_ucs)[i,5],
                 c=cmap(clusters[i-1]/max(clusters)), marker='o', s=20)
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    ax.set_zlabel('gamma')
    plt.show()

#    centroids,_ = kmeans(reduced_data[:,:3], num_clusters)
#    idx,_ = vq(reduced_data[:,:3],centroids)
