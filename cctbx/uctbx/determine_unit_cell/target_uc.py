from __future__ import division
import os
import sys
from libtbx import easy_pickle
import numpy as np
import math
#from scipy.cluster.vq import kmeans, kmeans2, vq
#import libtbx.load_env
from cctbx.uctbx.determine_unit_cell import NCDist
import  scipy.cluster.hierarchy as hcluster

class target:

  def __init__(self, path_to_integration_dir):
    """ Creates a list of (point group, unit cell) tuples, and a list of niggli cells from the recursively walked
        paths. Can take more than one argument for multiple folders."""
    self.pgs        = []
    self.niggli_ucs = []
    self.names      = []
    for arg in path_to_integration_dir:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          path = os.path.join(dirpath, filename)
          try:
            # Warn on error, but continue directory traversal.
            d = easy_pickle.load(path)
            pg   = d['pointgroup']
            uc   = d['current_orientation'][0].unit_cell()
            name = filename
          except KeyError:
              sys.stderr.write(
               "Could not extract point group and unit cell from %s\n" % path)
          except Exception:
              sys.stderr.write(
                    "Could not read %s\n" % path)
          else:
            self.pgs.append(pg)
            self.niggli_ucs.append(uc.niggli_cell().parameters())
            self.names.append(name)

  def print_ucs(self, outfile="niggli_ucs.csv"):
    out_str = ["File name, Point group, a, b, c, alpha, beta, gamma"]
    for i in range(len(self.names)):
       out_str.append("{}, {}, {}, {}, {}, {}, {}, {}".format(
                  self.names[i], self.pgs[i], 
                  self.niggli_ucs[i][0], self.niggli_ucs[i][1],
                  self.niggli_ucs[i][2], self.niggli_ucs[i][3],
                  self.niggli_ucs[i][4], self.niggli_ucs[i][5]))
    with open("{}.csv".format(outfile), 'w') as outfile:             
        outfile.write("\n".join(out_str))

  def find_distance(self, G6a, G6b, key):
    """ Retursn the distance between two cells, already in the G6 convention.
        Optionally (key =1)  trivial Euclidian or NCDist (key =2 )"""
    if key is 1: # trivial euclidian
      return math.sqrt(np.sum((G6a - G6b)**2))
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

  def cluster(self, threshold, method, linkage_method):

    """ Do basic hierarchical clustering on the Niggli
    cells created at the start """
    self.threshold=threshold
    dist_method = 2
    if dist_method is 2:
      print ("Using Andrews-Bernstein Distance from " +
      "Andrews & Bernstein J Appl Cryst 47:346 (2014).")
    # 1. Create a numpy array of G6 cells
    self.G6_cells = np.array([self.make_G6(n_cell) for n_cell in self.niggli_ucs])
    # 2. Do hierarchichal clustering on this, using the find_distance method above.
    import  scipy.spatial.distance as dist
    pair_distances = dist.pdist(self.G6_cells,
                                metric=lambda a, b: self.find_distance(a,b,dist_method))
    print "Distances have been calculated"
    self.this_linkage  = hcluster.linkage(pair_distances,
                                     method=linkage_method,
                                     metric=lambda a, b: self.find_distance(a,b,dist_method))
    self.clusters = hcluster.fcluster(self.this_linkage,
                                          threshold,
                                          criterion='distance')

    # 3. print out some information that is useful.
    print "{} clusters have been identified.".format(max(self.clusters))
    print "{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
                             "Num in cluster",
                             "Med_a", "Med_b", "Med_c",
                             "Med_alpha", "Med_beta", "Med_gamma")
    singletons = []
    for cluster in range(max(self.clusters)):
      this_cluster = np.array([self.niggli_ucs[i]
                               for i in range(len(self.niggli_ucs))
                               if self.clusters[i]==cluster+1])
      this_cluster_pg = ([self.pgs[i] for i in range(len(self.pgs))
                               if self.clusters[i]==cluster+1])
      assert len(this_cluster_pg) == len(this_cluster)
      all_pgs={}
      for pg in this_cluster_pg:
          if pg in all_pgs.keys():
              all_pgs[pg] += 1
          else:
              all_pgs[pg] = 1

      if len(this_cluster) != 1:
        pg_strings = ["{} images in {}".format(all_pgs[pg], pg) for pg in  all_pgs]
        point_group_string = ", ".join(pg_strings)+"."
        print "".join([("{:<14} {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f})" \
                        " {:<5.1f}({:<4.1f})").format(
          len(this_cluster),
          np.median(this_cluster[:,0]), np.std(this_cluster[:,0]),
          np.median(this_cluster[:,1]), np.std(this_cluster[:,1]),
          np.median(this_cluster[:,2]), np.std(this_cluster[:,2])),
          " {:<6.2f}({:<4.2f}) {:<6.2f}({:<4.2f}) {:<6.2f}({:<4.2f})".format(
          np.median(this_cluster[:,3]), np.std(this_cluster[:,3]),
          np.median(this_cluster[:,4]), np.std(this_cluster[:,4]),
          np.median(this_cluster[:,5]), np.std(this_cluster[:,5]))])
        print "  --> " +  point_group_string
      else:
        singletons.append("".join([("{:<14} {:<11.1f} {:<11.1f}" \
                          " {:<11.1f}").format(
          this_cluster_pg[0],
          this_cluster[0,0],
          this_cluster[0,1],
          this_cluster[0,2]),
          " {:<12.1f} {:<12.1f} {:<12.1f}".format(
          this_cluster[0,3],
          this_cluster[0,4],
          this_cluster[0,5]), '\n']))
    print "Standard deviations are in brackets."
    print  str(len(singletons)) + " singletons:  \n"
    print "{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
                             "Point group",
                             "Med_a", "Med_b", "Med_c",
                             "Med_alpha", "Med_beta", "Med_gamma")
    print "".join(singletons)

def plot_clusters(ucs, log=False, outname='clustering'):
    """ Plot Niggli cells -- one plot for (a,b,c) and one plot for
    (alpha, beta, gamma) -- colour coded by cluster index.  """
    
    import pylab
    fig = pylab.figure()
    hcluster.dendrogram(ucs.this_linkage,
                      labels=["{:<4.1f}, {:<4.1f}, {:<4.1f}, {:<4.1f}," +
                              "{:<4.1f}, {:<4.1f}".format(
                              x[0], x[1], x[2], x[3], x[4], x[5])
                              for x in ucs.niggli_ucs],
                      leaf_font_size=8,
                      color_threshold=ucs.threshold)
    ax=fig.gca()
    if log:
      ax.set_yscale("log")
    else:
      ax.set_ylim(-ax.get_ylim()[1]/100,ax.get_ylim()[1])
    fig.show()
    fig.savefig("{}_dendogram.pdf".format(outname))

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D # Special Import
    import matplotlib.cm as mpl_cmaps
    cmap = mpl_cmaps.Paired
    fig = plt.figure('unit_cells_dimensions')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('a [A]')
    ax.set_ylabel('b [A]')
    ax.set_zlabel('c [A]')
    for i in range(len(ucs.niggli_ucs)):
      ax.scatter(np.array(ucs.niggli_ucs)[i,0],
                 np.array(ucs.niggli_ucs)[i,1],
                 np.array(ucs.niggli_ucs)[i,2],
                 c=cmap(ucs.clusters[i-1]/max(ucs.clusters)), marker='o', s=20)
    fig = plt.figure('unit_cells_angles')
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(ucs.niggli_ucs)):
      ax.scatter(np.array(ucs.niggli_ucs)[i,3],
                 np.array(ucs.niggli_ucs)[i,4],
                 np.array(ucs.niggli_ucs)[i,5],
                 c=cmap(ucs.clusters[i-1]/max(ucs.clusters)), marker='o', s=20)
    ax.set_xlabel('alpha')
    ax.set_ylabel('beta')
    ax.set_zlabel('gamma')
    plt.show()

#    centroids,_ = kmeans(reduced_data[:,:3], num_ucs.clusters)
#    idx,_ = vq(reduced_data[:,:3],centroids)
