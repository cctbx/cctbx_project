from __future__ import division
from __future__ import print_function
import os
import math
from cctbx.uctbx.determine_unit_cell import NCDist
import scipy.cluster.hierarchy as hcluster
import numpy as np
import json
from .SingleFrame import SingleFrame
import logging


class Clu_BELIEVE_THIS_WHOLE_DIRECTORY_IS_DEAD_CODE_20171120_ster:
  def __init__(self, data, cname, info, log_level='INFO'):
    """ Contains a list of SingFrame objects, as well as information about these
    as a cluster (e.g. mean unit cell)."""
    self.cname = cname
    self.members = data
    self.info = info

    # Calculate medians and stdevs
    unit_cells = np.zeros([len(self.members), 6])
    self.pg_composition = {}
    for i, member in enumerate(self.members):
      unit_cells[i, :] = member.uc
      # Calculate point group composition
      if member.pg in self.pg_composition.keys():
        self.pg_composition[member.pg] += 1
      else:
        self.pg_composition[member.pg] = 1

    self.medians = np.median(unit_cells, 0).tolist()
    self.stdevs = np.std(unit_cells, 0).tolist()
    #ToDo
    self.res = None


  @classmethod
  def from_directories(cls, path_to_integration_dir,
                       _prefix='cluster_from_file'):
    """Constructor to get a cluster from pickle files, from the recursively
    walked paths. Can take more than one argument for multiple folders.
    usage: Cluster.from_directories(..)"""
    data = []
    for arg in path_to_integration_dir:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          path = os.path.join(dirpath, filename)
          data.append(SingleFrame(path, filename))
    return cls(data, _prefix,
               'Made from files in {}'.format(path_to_integration_dir[:]))

  @classmethod
  def from_json(cls, json_file, _prefix='cluster_from_json'):
    """ Does not work!! Do not use! """
    with open(json_file, 'rb') as inFile:
      data = json.load(inFile)
    return cls(data, _prefix,
               'Made from {}'.format(json_file))

  def make_sub_cluster(self, new_members, new_prefix, new_info):
    """ Make a sub-cluster from a list of cluster indicies from the old
    SingleFrame array.
    """
    return Cluster(new_members, new_prefix,
                   '{}\n{} Next filter {}\n{}\n{} of {} images passed on to this cluster'.format(
                     self.info, '#' * 30, '#' * 30, new_info, len(new_members), len(self.members)))

  def print_ucs(self):
    outfile = "{}_niggli_ucs".format(self.cname)
    out_str = ["File name, Point group, a, b, c, alpha, beta, gamma"]
    for image in self.members:
      out_str.append("{}, {}, {}, {}, {}, {}, {}, {}".format(
        image.name, image.pg,
        image.uc[0], image.uc[1],
        image.uc[2], image.uc[3],
        image.uc[4], image.uc[5]))
    with open("{}.csv".format(outfile), 'w') as _outfile:
      _outfile.write("\n".join(out_str))

  def point_group_filer(self, point_group):
    """ Return all the SingleFrames that have a given pointgroup. """
    new_prefix = '{}_only'.format(point_group)
    new_info = 'Cluster filtered by for point group {}.'.format(point_group)
    return self.make_sub_cluster([image for image in self.members if image.pg == point_group],
                                 new_prefix,
                                 new_info)

  def total_intensity_filter(self, res='', completeness_threshold=0.95, plot=False):
    """ Creates a sub-cluster using the highest total intensity images that yield a dataset specified by:
          res -- desired resolution. Defaults to that of the dataset.
          completeness -- the desired completeness of the subset
          multiplicity -- the desired multiplicity of the subset
    """
    logging.info(("Performing intensity filtering, aiming for {}% overall completenes"
                  " at {} A resolution").format(completeness_threshold * 100, res))

    # 0. Check that the cluster has consistent point_group (completness doesn't mean much otherwise...
    assert all(i.pg == self.members[0].pg for i in self.members)

    # 1. Sort SingleFrames by total intensity
    sorted_cluster = sorted(self.members, key=lambda y: -1 * y.total_i)

    if plot:
      from matplotlib import pyplot as plt
      plt.plot([x.total_i for x in sorted_cluster])
      plt.show()

    if res == '':
      res = sorted_cluster[0].d_min()  # Use the high-res limit from the brightest image. ToDo: make this better.
      logging.warning("no resolution limit specified, using the res limit of the top-rankeed image: {} A".format(res))

    # 2. Incrementally merge frames until criterion are matched

    temp_miller_indicies = sorted_cluster[0].miller_array
    for idx, image in enumerate([x.miller_array for x in sorted_cluster[1:]]):
      temp_miller_indicies = temp_miller_indicies.concatenate(image, assert_is_similar_symmetry=False)
      current_completeness = temp_miller_indicies.merge_equivalents().array().completeness()
      logging.debug("{} images: {:.2f}% complete".format(idx, current_completeness * 100))
      if current_completeness <= completeness_threshold:
        temp_miller_indicies.concatenate(image, assert_is_similar_symmetry=False)
        if idx + 1 == len(sorted_cluster[1:]):
          logging.warning("Desired completeness could not be achieved, sorry.")
          file_threshold = idx
          break
      else:
        file_threshold = idx
        break

    return self.make_sub_cluster(sorted_cluster[:file_threshold],
                                 'I_threshold_d{}_{}comp'.format(res, completeness_threshold),
                                 ('Subset cluster made using total_intensity_filter() with'
                                  '\nRes={}\ncompleteness_threshold={}').format(res,
                                                                                completeness_threshold))

  def ab_cluster(self, threshold=10000, method='distance', linkage_method='single', log=False, plot=False):
    """ Do basic hierarchical clustering using the Andrews-Berstein distance
    on the Niggli cells """
    print("Hierarchical clustering of unit cells:")
    import scipy.spatial.distance as dist

    print("Using Andrews-Bernstein Distance from Andrews & Bernstein J Appl Cryst 47:346 (2014).")

    def make_g6(uc):
      """ Take a reduced Niggli Cell, and turn it into the G6 representation """
      a = uc[0] ** 2
      b = uc[1] ** 2
      c = uc[2] ** 2
      d = 2 * uc[1] * uc[2] * math.cos(uc[3])
      e = 2 * uc[0] * uc[2] * math.cos(uc[4])
      f = 2 * uc[0] * uc[1] * math.cos(uc[5])
      return [a, b, c, d, e, f]

    # 1. Create a numpy array of G6 cells
    g6_cells = np.array([make_g6(image.uc)
                         for image in self.members])

    # 2. Do hierarchichal clustering, using the find_distance method above.
    pair_distances = dist.pdist(g6_cells,
                                metric=lambda a, b: NCDist(a, b))
    logging.debug("Distances have been calculated")
    this_linkage = hcluster.linkage(pair_distances,
                                    method=linkage_method,
                                    metric=lambda a, b: NCDist(a, b))
    cluster_ids = hcluster.fcluster(this_linkage,
                                    threshold,
                                    criterion=method)
    logging.debug("Clusters have been calculated")
    # Create an array of sub-cluster objects from the clustering
    sub_clusters = []
    for cluster in range(max(cluster_ids)):
      info_string = ('Made using ab_cluster with t={},'
                     ' {} method, and {} linkage').format(threshold,
                                                          method,
                                                          linkage_method)
      sub_clusters.append(self.make_sub_cluster([self.members[i]
                                                 for i in
                                                 range(len(self.members))
                                                 if
                                                 cluster_ids[i] == cluster + 1],
                                                'cluster_{}'.format(
                                                  cluster + 1),
                                                info_string))

    # 3. print out some information that is useful.
    out_str = "{} clusters have been identified.".format(max(cluster_ids))
    out_str += "\n{:^5} {:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
      "C_id",
      "Num in cluster",
      "Med_a", "Med_b", "Med_c",
      "Med_alpha", "Med_beta", "Med_gamma")
    singletons = []
    for cluster in sub_clusters:
      if len(cluster.members) != 1:

        sorted_pg_comp = sorted(cluster.pg_composition.items(), key=lambda x: -1 * x[1])
        pg_strings = ["{} in {}".format(pg[1], pg[0])
                      for pg in sorted_pg_comp]
        point_group_string = ", ".join(pg_strings) + "."
        out_str += ("\n{:^5} {:^14} {:<5.1f}({:<4.1f}) {:<5.1f}({:<4.1f})"
                    " {:<5.1f}({:<4.1f}) {:<6.2f}({:<4.2f}) {:<6.2f}"
                    "({:<4.2f}) {:<6.2f}({:<4.2f})").format(
          cluster.cname,
          len(cluster.members),
          cluster.medians[0], cluster.stdevs[0],
          cluster.medians[1], cluster.stdevs[1],
          cluster.medians[2], cluster.stdevs[2],
          cluster.medians[3], cluster.stdevs[3],
          cluster.medians[4], cluster.stdevs[4],
          cluster.medians[5], cluster.stdevs[5])
        out_str += "\n" + point_group_string
      else:
        singletons.append("".join([("{:<14} {:<11.1f} {:<11.1f} {:<11.1f}"
                                    "{:<12.1f} {:<12.1f} {:<12.1f}").format(
          cluster.pg_composition.keys()[0],
          cluster.members[0].uc[0], cluster.members[0].uc[1],
          cluster.members[0].uc[2], cluster.members[0].uc[3],
          cluster.members[0].uc[4], cluster.members[0].uc[5]),
                                   '\n']))
    out_str += "\nStandard deviations are in brackets."
    out_str += "\n" + str(len(singletons)) + " singletons:"
    out_str += "\n{:^14} {:<11} {:<11} {:<11} {:<12} {:<12} {:<12}".format(
      "Point group",
      "a", "b", "c",
      "alpha", "beta", "gamma")
    out_str += "".join(singletons)
    print(out_str)

    if plot:
      import matplotlib.pyplot as plt

      fig = plt.figure("Distance Dendogram")
      hcluster.dendrogram(this_linkage,
                          labels=[image.name for image in self.members],
                          leaf_font_size=8,
                          color_threshold=threshold)
      ax = fig.gca()
      if log:
        ax.set_yscale("log")
      else:
        ax.set_ylim(-ax.get_ylim()[1] / 100, ax.get_ylim()[1])
      fig.savefig("{}_dendogram.pdf".format(self.cname))
      plt.show()

    return sub_clusters

  def dump_file_list(self, out_file_name=None):
    if out_file_name is None:
      out_file_name = self.cname

    with open("{}.members".format(out_file_name), 'wb') as outfile:
      for i in self.members:
        outfile.write(i.path + "\n")

def cluster_pca(self):
  """ Should use BLEND clustering protocols in python (Foaldi et al. Acta D.
  2013). i.e. filter for parameters that vary, do PCA, then ward linkage
  clustering on this. """
  # Will come back to this soon <-- Oli
  #columns_to_use = [True if np.std(self.niggli_ucs[:, i])
  #                  else False for i in range(6)]
