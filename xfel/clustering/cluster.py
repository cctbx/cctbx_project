""" This module is designed to provide tools to deal with groups of serial
crystallography images.

The class Cluster allows the creation, storage and manipulation of these sets
of frames. Methods exist to create sub-clusters (new cluster objects) or to act
on an existing cluster, e.g. to plot the unit cell distributions.

**Author:**   Oliver Zeldin <zeldin@stanford.edu>
"""
from __future__ import division
__author__ = 'zeldin'

from cctbx.array_family import flex
import os
import math
import logging
from xfel.clustering.singleframe import SingleFrame, ImageNode
from cctbx.uctbx.determine_unit_cell import NCDist
import numpy as np
import matplotlib.patheffects as patheffects
import matplotlib.pyplot as plt
import random


class Cluster:
  """Class for operation on groups of single XFEL images (here described by
  SingleFrame objects) as cluster objects.

  Objects can be created directily using the __init__ method, or by using a
  classmethod e.g. to create an object from a folder full of integration
  pickles. Whenever a method can plot, there is the option of
  passing it an appropriate number of matplotlib axes objects, which will then
  get returned for use in composite plots. See cluster.42 for an example.
  If no axes are passed, the methods will just plot the result to the screen.
  Clustering filters can act on these to break them up into cluster objects with
  different members. A 'filter' is just a clustering procedure that puts the
  passes and fails into different clusters. This is acheived through the
  make_sub_cluster() method. This also keeps track of a sub-clusters heritage
  through the .info string, which is appended to. The idea is to be able to
  write filter scripts for each data.
  e.g ::

      >>> test_cluster = Cluster.from_directories(["~/test_data"],
      ...                                        'test_script')
      >>> P3_only = test_cluster.point_group_filer('P3')
      >>> sub_clusters = P3_only.ab_cluster(1200)
      >>> big_cluster = max(sub_clusters, key=lambda x: len(x.members))
      >>> best_data = big_cluster.total_intensity_filter(res=6.5,
      ...                                               completeness_threshold=0.1,
      ...                                               plot=False)
      >>> print best_data.info

  Subsequent postrefinenment/merging programs can be called on an output
  cluster.lst file::
      >>> prime.postrefine params.phil $(cat cluster.lst)
  """

  def __init__(self, data, cname, info):
    """
    Builds a cluster from a list of SingleFrame objects, as well as information
    about these as a cluster (e.g. mean unit cell).

    :param data: a list of SingleFrame objects
    :param cname: the name of the cluster, as a string.
    :param info: an info-string for the cluster.
    """

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
                       _prefix='cluster_from_dir',
                       **kwargs):
    """Constructor to get a cluster from pickle files, from the recursively
    walked paths. Can take more than one argument for multiple folders.
    usage: Cluster.from_directories(..)
    :param path_to_integration_dir: list of directories containing pickle files.
    Will be searched recursively.
    :param use_b: Boolean. If True, intialise Scale and B. If false, use only
    mean intensity scalling.
    """
    data = []
    for arg in path_to_integration_dir:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          path = os.path.join(dirpath, filename)
          this_frame = SingleFrame(path, filename, **kwargs)
          if hasattr(this_frame, 'miller_array'):
            data.append(this_frame)
          else:
            logging.info('skipping file {}'.format(filename))
    return cls(data, _prefix,
               'Made from files in {}'.format(path_to_integration_dir[:]))

  @classmethod
  def from_files(cls, pickle_list,
                       _prefix='cluster_from_file',
                       use_b=True):
    """Constructor to get a cluster from a list of pickle files.
    :param pickle_list: list of pickle files
    :param use_b: Boolean. If True, intialise Scale and B. If false, use only
    mean intensity scalling.
    """
    data = []
    for filename in pickle_list:
      name_only = filename.split('/')[-1]
      this_frame = SingleFrame(filename, name_only, use_b=use_b)
      if hasattr(this_frame, 'name'):
        data.append(this_frame)
      else:
        logging.info('skipping file {}'.format(filename))
    return cls(data, _prefix, 'Made by Cluster.from_files')

  def make_sub_cluster(self, new_members, new_prefix, new_info):
    """ Make a sub-cluster from a list of SingleFrame objects from the old
    SingleFrame array.

    :param new_members: a new set of SingleFrame objects, typically a subset of
    the current cluster.
    :param new_prefix: prefix to be passed directly to __init__
    :param new_info: new information about this cluster. This is inteligently
    appended to the old info string so that the history of sub-clusters is
    tracted.
    """
    return Cluster(new_members, new_prefix,
                   ('{}\n{} Next filter {}\n{}\n{} of {} images passed'
                    'on to this cluster').format(
                     self.info, '#' * 30, '#' * 30, new_info,
                     len(new_members), len(self.members)))

  def print_ucs(self):
    """ Prints a list of all the unit cells in the cluster to CSV."""
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

  def point_group_filter(self, point_group):
    """ Return all the SingleFrames that have a given pointgroup. """
    new_prefix = '{}_only'.format(point_group)
    new_info = 'Cluster filtered by for point group {}.'.format(point_group)
    return self.make_sub_cluster([image
                                  for image
                                  in self.members
                                  if image.pg == point_group],
                                 new_prefix,
                                 new_info)

  def total_intensity_filter(self, res='',
                             completeness_threshold=0.95,
                             plot=False):
    """
    .. note::
      This is still in development, use at own risk!

    Creates an optimal sub-cluster using the fewest images that fulfil the
    criteria defined by:
    :param res: desired resolution. Defaults to that of the dataset.
    :param completeness: the desired completeness of the subset
    :param multiplicity: the desired multiplicity of the subset
    """
    logging.info(("Performing intensity filtering, aiming for {}% overall "
                  "completenes at {} A resolution").format(
      completeness_threshold * 100, res))

    # 0. Check that the cluster has consistent point_group (completness doesn't
    #  mean much otherwise...
    assert all(i.pg == self.members[0].pg for i in self.members)

    # 1. Sort SingleFrames by total intensity
    sorted_cluster = sorted(self.members, key=lambda y: -1 * y.total_i)

    if plot:
      plt.plot([x.total_i for x in sorted_cluster])
      plt.show()

    if res == '':
      res = sorted_cluster[0].d_min()  # Use the high-res limit from the
      # brightest image. ToDo: make this better
      logging.warning(("no resolution limit specified, using the res limit of"
                       "the top-rankeed image: {} A").format(res))

    # 2. Incrementally merge frames until criterion are matched

    temp_miller_indicies = sorted_cluster[0].miller_array
    for idx, image in enumerate((x.miller_array for x in sorted_cluster[1:])):
      temp_miller_indicies = temp_miller_indicies. \
        concatenate(image, assert_is_similar_symmetry=False)
      current_completeness = temp_miller_indicies.merge_equivalents() \
                                                .array() \
                                                .completeness()
      logging.debug(
        "{} images: {:.2f}% complete".format(idx, current_completeness * 100))
      if current_completeness <= completeness_threshold:
        temp_miller_indicies.concatenate(image,
                                         assert_is_similar_symmetry=False)
        if idx + 1 == len(sorted_cluster[1:]):
          logging.warning("Desired completeness could not be acheived, sorry.")
          file_threshold = idx
          break
      else:
        file_threshold = idx
        break

    return self.make_sub_cluster(sorted_cluster[:file_threshold],
                                 'I_threshold_d{}_{}comp'.format(res,
                                                        completeness_threshold),
                                 ("Subset cluster made using "
                                  "total_intensity_filter() with"
                                  "\nRes={}\ncompleteness_threshold={}").format(
                                   res,
                                   completeness_threshold))


  def ab_cluster(self, threshold=10000, method='distance',
                 linkage_method='single', log=False,
                 ax=None, write_file_lists=True, fast=False, doplot=True,
                 labels='default'):
    """
    Hierarchical clustering using the unit cell dimentions.

    :param threshold: the threshold to use for prunning the tree into clusters.
    :param fast: if True, use simple euclidian distance, otherwise, use Andrews-Berstein distance from Andrews & Bernstein J Appl Cryst 47:346 (2014) on the Niggli cells.
    :param method: which clustering method from scipy to use when creating the tree (see scipy.cluster.hierarchy)
    :param linkage_method: which linkage method from scipy to use when creating the linkages. x (see scipy.cluster.hierarchy)
    :param log: if True, use log scale on y axis.
    :param ax: if a matplotlib axes object is provided, plot to this. Otherwise, create a new axes object and display on screen.
    :param write_file_lists: if True, write out the files that make up each cluster.
    :return: A list of Clusters ordered by largest Cluster to smallest
    .. note::
      Use 'fast' option with caution, since it can cause strange behaviour
      around symmetry boundaries.
    """

    logging.info("Hierarchical clustering of unit cells")
    import scipy.spatial.distance as dist
    import scipy.cluster.hierarchy as hcluster

    # 1. Create a numpy array of G6 cells
    g6_cells = np.array([SingleFrame.make_g6(image.uc)
                         for image in self.members])

    # 2. Do hierarchichal clustering, using the find_distance method above.
    if fast:
      logging.info("Using Euclidean distance")
      pair_distances = dist.pdist(g6_cells, metric='euclidean')
      logging.info("Distances have been calculated")
      this_linkage = hcluster.linkage(pair_distances,
                                      method=linkage_method,
                                      metric='euclidean')
    else:
      logging.info("Using Andrews-Bernstein distance from Andrews & Bernstein "
                   "J Appl Cryst 47:346 (2014)")
      pair_distances = dist.pdist(g6_cells,
                                metric=lambda a, b: NCDist(a, b))
      logging.info("Distances have been calculated")
      this_linkage = hcluster.linkage(pair_distances,
                                      method=linkage_method,
                                      metric=lambda a, b: NCDist(a, b))


    cluster_ids = hcluster.fcluster(this_linkage,
                                    threshold,
                                    criterion=method)
    logging.debug("Clusters have been calculated")

    # 3. Create an array of sub-cluster objects from the clustering
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

    sub_clusters = sorted(sub_clusters, key=lambda x: len(x.members))
    # Rename to order by size
    for num, cluster in enumerate(sub_clusters):
      cluster.cname = 'cluster_{}'.format(num + 1)

    # 3.5 optionally write out the clusters to files.
    if write_file_lists:
      for cluster in sub_clusters:
        if len(cluster.members) > 1:
          cluster.dump_file_list(out_file_name="{}.lst".format(cluster.cname))

    if doplot:
      if labels == 'default':
        if len(self.members) > 100:
          labels = ['' for _ in self.members]
        else:
          labels = [image.name for image in self.members]

      # 4. Plot a dendogram to the axes if no axis is passed, otherwise just
      #    return the axes object
      if ax is None:
        fig = plt.figure("Distance Dendogram")
        ax = fig.gca()
        direct_visualisation = True
      else:
        direct_visualisation = False

      hcluster.dendrogram(this_linkage,
                          labels=labels,
                          leaf_font_size=8, leaf_rotation=90.0,
                          color_threshold=threshold, ax=ax)

      if log:
        ax.set_yscale("symlog", linthreshx=(-1,1))
      else:
        ax.set_ylim(-ax.get_ylim()[1] / 100, ax.get_ylim()[1])

      if direct_visualisation:
        fig.savefig("{}_dendogram.pdf".format(self.cname))
        plt.show()

    return sub_clusters, ax

  def dump_file_list(self, out_file_name=None):
    """ Dumps a list of paths to inegration pickle files to a file. One
    line per image. Provides easy input into post-refinement programs.

    :param out_file_name: the output file name.
    """
    if out_file_name is None:
      out_file_name = self.cname

    with open(out_file_name, 'wb') as outfile:
      for i in self.members:
        outfile.write(i.path + "\n")

  def visualise_orientational_distribution(self, axes_to_return=None,
                                           cbar=True):

    """ Creates a plot of the orientational distribution of the unit cells.

    :param axes_to_return: if None, print to screen, otherwise, requires 3 axes objects, and will return them.
    :param cbar: boolean to specify if a color bar should be used.
    """
    from mpl_toolkits.basemap import Basemap
    import scipy.ndimage as ndi

    def cart2sph(x, y, z):
      # cctbx (+z to source, y to ceiling) to
      # lab frame (+x to source, z to ceiling)
      z, x, y = x, y, z
      dxy = np.sqrt(x ** 2 + y ** 2)
      r = np.sqrt(dxy ** 2 + z ** 2)
      theta = np.arctan2(y, x)
      phi = np.arctan2(z, dxy)  # angle of the z axis relative to xy plane
      theta, phi = np.rad2deg([theta, phi])
      return theta % 360, phi, r

    def xy_lat_lon_from_orientation(orientation_array, axis_id):
      logging.debug("axis_id: {}".format(axis_id))
      dist = math.sqrt(orientation_array[axis_id][0] ** 2 +
                       orientation_array[axis_id][1] ** 2 +
                       orientation_array[axis_id][2] ** 2)
      flon, flat, bla = cart2sph(orientation_array[axis_id][0] / dist,
                                 orientation_array[axis_id][1] / dist,
                                 orientation_array[axis_id][2] / dist)
      x, y = euler_map(flon, flat)
      return x, y, flon, flat

    orientations = [flex.vec3_double(flex.double(
      image.orientation.direct_matrix()))
      for image in self.members]

    space_groups = [image.orientation.unit_cell().lattice_symmetry_group()
                    for image in self.members]

    # Now do all the plotting
    if axes_to_return is None:
      plt.figure(figsize=(10, 14))
      axes_to_return = [plt.subplot2grid((3, 1), (0, 0)),
                        plt.subplot2grid((3, 1), (1, 0)),
                        plt.subplot2grid((3, 1), (2, 0))]
      show_image = True
    else:
      assert len(axes_to_return) == 3, "If using axes option, must hand" \
                                       " 3 axes to function."
      show_image = False

    axis_ids = [0, 1, 2]
    labels = ["a",
              "b",
              "c"]

    for ax, axis_id, label in zip(axes_to_return, axis_ids, labels):

      # Lists of x,y,lat,long for the master orientation, and for all
      # symmetry mates.
      x_coords = []
      y_coords = []
      lon = []
      lat = []
      sym_x_coords = []
      sym_y_coords = []
      sym_lon = []
      sym_lat = []
      euler_map = Basemap(projection='eck4', lon_0=0, ax=ax)

      for orientation, point_group_type in zip(orientations, space_groups):

        # Get position of main spots.
        main_x, main_y, main_lon, main_lat \
          = xy_lat_lon_from_orientation(list(orientation), axis_id)
        x_coords.append(main_x)
        y_coords.append(main_y)
        lon.append(main_lon)
        lat.append(main_lat)

        # Get position of symetry mates
        symmetry_operations = list(point_group_type.smx())[1:]
        for mx in symmetry_operations:
          rotated_orientation = list(mx.r().as_double() * orientation)  # <--
          # should make sense if orientation was a vector, not clear what is
          # going on since orientation is a matrix. Or, make some test cases
          # with 'orientation' and see if the behave as desired.
          sym_x, sym_y, sym_lo, sym_la \
            = xy_lat_lon_from_orientation(rotated_orientation, axis_id)
          #assert (sym_x, sym_y) != (main_x, main_y)
          sym_x_coords.append(sym_x)
          sym_y_coords.append(sym_y)
          sym_lon.append(sym_lo)
          sym_lat.append(sym_la)

      # Plot each image as a yellow sphere
      logging.debug(len(x_coords))
      euler_map.plot(x_coords, y_coords, 'oy',
                     markersize=4,
                     markeredgewidth=0.5)

      # Plot the symetry mates as black crosses
      #euler_map.plot(sym_x_coords, sym_y_coords, 'kx')

      # Use a histogram to bin the data in lattitude/longitude space, smooth it,
      # then plot this as a contourmap. This is for all the symetry-related
      # copies
      density_hist = np.histogram2d(lat + sym_lat, lon + sym_lon,
                                    bins=[range(-90, 91), range(0, 361)])
      smoothed = ndi.gaussian_filter(density_hist[0], (15, 15), mode='wrap')
      local_intensity = []
      x_for_plot = []
      y_for_plot = []
      for _lat in range(0, 180):
        for _lon in range(0, 360):
          _x, _y = euler_map(density_hist[2][_lon], density_hist[1][_lat])
          x_for_plot.append(_x)
          y_for_plot.append(_y)
          local_intensity.append(smoothed[_lat, _lon])
      cs = euler_map.contourf(np.array(x_for_plot),
                              np.array(y_for_plot),
                              np.array(local_intensity), tri=True)

      #  Pretty up graph
      if cbar:
        _cbar = plt.colorbar(cs, ax=ax)
        _cbar.ax.set_ylabel('spot density [AU]')
      middle = euler_map(0, 0)
      path_effect = [patheffects.withStroke(linewidth=3, foreground="w")]
      euler_map.plot(middle[0], middle[1], 'o', markersize=10, mfc='none')
      euler_map.plot(middle[0], middle[1], 'x', markersize=8)
      ax.annotate("beam", xy=(0.52, 0.52), xycoords='axes fraction',
                  size='medium', path_effects=path_effect)
      euler_map.drawmeridians(np.arange(0, 360, 60),
                              labels=[0, 0, 1, 0],
                              fontsize=10)
      euler_map.drawparallels(np.arange(-90, 90, 30),
                              labels=[1, 0, 0, 0],
                              fontsize=10)
      ax.annotate(label, xy=(-0.05, 0.9), xycoords='axes fraction',
                  size='x-large', weight='demi')

    if show_image:
      plt.show()

    return axes_to_return


  def intensity_statistics(self, ax=None):
    """
    Uses the per-frame B and G fits (gradient and intercept of the ln(i) vs
    (sin(theta)/lambda)**2 plot) to create three agregate plots:
    1) histogram of standard errors on the per-frame fits
    2) histogram of B factors
    3) scatter  plot of intercept vs. gradient (G vs. B)

    :param ax: optionally hand the method three matplotlib axes objects to plot onto. If not specified, will plot the data.
    :return: the three axes, with the data plotted onto them.
    """
    if ax is None:
      plt.figure(figsize=(10, 14))
      axes_to_return = [plt.subplot2grid((3, 1), (0, 0)),
                        plt.subplot2grid((3, 1), (1, 0)),
                        plt.subplot2grid((3, 1), (2, 0))]
      show_image = True
    else:
      assert len(ax) == 3, "If using axes option, must hand" \
                                       " 3 axes to function."
      axes_to_return = ax
      show_image = False

    errors = [i.wilson_err['Standard Error'] for i in self.members]
    axes_to_return[0].hist(errors, 50, range=[0, 200])
    axes_to_return[0].set_title("Distribution of Standard Errors on the Wilson fit")

    rs = [-1 * i.minus_2B / 2 for i in self.members]
    axes_to_return[1].hist(rs, 50, range=[-50, 200])
    axes_to_return[1].set_title("Distribution of B values for the Wilson plot")

    axes_to_return[2].plot([i.G for i in self.members],
             [-1 * i.minus_2B / 2 for i in self.members], 'x')
    axes_to_return[2].set_xlabel("G")
    axes_to_return[2].set_ylabel("B")
    axes_to_return[2].set_title("G and B for all members")

    plt.tight_layout()

    if show_image:
      plt.show()

    return axes_to_return

  def all_frames_intensity_stats(self, ax=None, smoothing_width=2000):
    """
    Goes through all frames in the cluster, and plots all the partial intensites.
    Then does a linear fit and rolling average on these.

    :param smoothing_width: the width of the smoothing window.
    :param ax: Optional matplotlib axes object to plot to. Otherwise, plot to screen.
    :return: the axis, with the data plotted onto it.
    """
    from scipy.stats import linregress
    from xfel.clustering.singleframe import SingleFrame as Sf

    if ax is None:
      fig = plt.figure("All images intensity statistics")
      ax = fig.gca()
      direct_visualisation = True
    else:
      direct_visualisation = False


    all_logi = []
    all_one_over_d_squared = []

    for frame in self.members:
      all_logi.append(frame.log_i)
      all_one_over_d_squared.append(frame.sinsqtheta_over_lambda_sq)

    all_logi = np.concatenate(all_logi)
    all_one_over_d_squared = np.concatenate(all_one_over_d_squared)

    plotting_data = sorted(zip(all_logi, all_one_over_d_squared),
                           key = lambda x: x[1])

    log_i, one_over_d_square = zip(*[i for i in plotting_data
                                     if i[0] >=0])
    minus_2B, G, r_val, _, std_err = linregress(one_over_d_square, log_i)
    fit_info = "G: {}, -2B: {}, r: {}, std_err: {}".format(G, minus_2B,
                                                            r_val, std_err)

    smooth = Sf._moving_average(log_i, n=smoothing_width)
    ax.plot(one_over_d_square, log_i, 'bo', ms=1)
    ax.plot(one_over_d_square[smoothing_width - 1:], smooth,'--r', lw=2)
    plt.xlim([0, max(one_over_d_square)])
    ax.plot([0, -1 * G / minus_2B], [G, 0], 'y-', lw=2)
    plt.xlabel("(sin(theta)/lambda)^2")
    plt.ylabel("ln(I)")
    plt.title("Simple Wilson fit\n{}".format(fit_info))
    plt.tight_layout()

    if direct_visualisation:
      fig.savefig("{}_dendogram.pdf".format(self.cname))
      plt.show()

    return ax


class Edge:
  """
  .. note::
    Developmental code. Do not use without contacting zeldin@stanford.edu

  Defines an undirected edge in a graph. Contains the connecting vertices, and a
  weight.
  """

  def __init__(self, vertex_a, vertex_b, weight):
    self.vertex_a = vertex_a
    self.vertex_b = vertex_b
    self.weight = weight

  def mean_residual(self):
    """
    :return: the mean residual of the absolute value of the all the weigts on this edge.
    """
    return np.mean(np.abs(self.residuals()))

  def other_vertex(self, vertex):
    """
    Simple method to get the other vertex along and edge.

    :param vertex: a vertex that is on one end of this edge
    :return: the vertex at the other end of the edge
    """
    assert vertex == self.vertex_a or vertex == self.vertex_b
    if vertex is self.vertex_a:
      return self.vertex_b
    elif vertex is self.vertex_b:
      return self.vertex_a

  def residuals(self):
    """
    Calculates the edge residual, as defined as the sum over all common miller
    indices of:
      log(scale * partiality of a) - log(scale * partiality of b) - log(I_a/I_b)

    :return: the residual score for this edge
    """
    # 1. Create flex selection array
    # 2. Trim these so that they only contain common reflections
    # 3. Calculate residual

    partialities_a = self.vertex_a.partialities
    partialities_b = self.vertex_b.partialities
    scales_a = self.vertex_a.scales
    scales_b = self.vertex_b.scales

    mtch_indcs = self.vertex_a.miller_array. \
      match_indices(self.vertex_b.miller_array,
                    assert_is_similar_symmetry=False)

    va_selection = mtch_indcs.pair_selection(0)
    vb_selection = mtch_indcs.pair_selection(1)

    sp_a = partialities_a.select(va_selection) * scales_a.select(va_selection)
    sp_b = partialities_b.select(vb_selection) * scales_b.select(vb_selection)

    ia_over_ib = self.vertex_a.miller_array.data().select(va_selection) / \
                 self.vertex_b.miller_array.data().select(vb_selection)
    residuals = (flex.log(sp_a) - flex.log(sp_b) - flex.log(ia_over_ib))
    residuals = residuals.as_numpy_array()
    #logging.debug("Mean Residual: {}".format(np.mean(residuals)))
    return residuals[~np.isnan(residuals)]


class Graph(Cluster):
  """
  .. note::
    Developmental code. Do not use without contacting zeldin@stanford.edu

  Class for treating a graph of serial still XFEL shots as a graph.
  """
  def __init__(self, vertices, min_common_reflections=10):
    """
    Extends the constructor from cluster.Cluster to describe the cluster as a
    graph.

    :param min_common_reflections: number of reflections two images must have in
    common for an edge to be created.
    """
    Cluster.__init__(self, vertices, "Graph cluster", "made as a graph")
    self.common_miller_threshold = min_common_reflections
    self.edges = self._make_edges()


  def _make_edges(self):
    all_edges = []
    for a, vertex1 in enumerate(self.members):
      for b, vertex2 in enumerate(self.members[:a]):
        miller1 = vertex1.miller_array
        miller2 = vertex2.miller_array
        edge_weight = miller1.common_set(miller2,
                                         assert_is_similar_symmetry=False).size()
        if edge_weight >= self.common_miller_threshold:
          logging.debug("Making edge: ({}, {}) {}, {}"
                        .format(a, b, vertex1.name, vertex2.name))
          this_edge = Edge(vertex1, vertex2, edge_weight)
          logging.debug("Edge weight: {}".format(edge_weight))
          all_edges.append(this_edge)
          vertex1.edges.append(this_edge)
          vertex2.edges.append(this_edge)
    return all_edges


  @classmethod
  def from_directories(cls, folder, **kwargs):
    """
    Creates a Graph class instance from a directory of files (recursively)

    :param folder: the root directory containing integration pickles.
    :return: A Graph instance.
    """
    def add_frame(all_vertices, dirpath, filename, **kwargs):
      path = os.path.join(dirpath, filename)
      this_frame = ImageNode(path, filename, **kwargs)
      if hasattr(this_frame, 'miller_array'):
        all_vertices.append(this_frame)
      else:
        logging.info('skipping file {}'.format(filename))

    nmax = kwargs.get('nmax', None)

    all_vertices = []
    for arg in folder:
      for (dirpath, dirnames, filenames) in os.walk(arg):
        for filename in filenames:
          if not nmax:
            add_frame(all_vertices, dirpath, filename, **kwargs)
          else:
            if len(all_vertices) <= nmax:
              add_frame(all_vertices, dirpath, filename)
    return cls(all_vertices)

  def make_nx_graph(self, edge_width=False, saveas=None, show_singletons=True,
                    **kwargs):
    """
    Create a networkx Graph, and visualise it.

    :param edge_width: if edge widths should be proportional to number of common miller indices.
    :param saveas: optional filename to save the graph as a pdf as. Displays graph if not specified.
    :param pos: an optional networkx position dictionairy for plotting the graph
    :return: the position dictionairy, so that mutiple plots can be made using the same set of positions.
    """
    import networkx as nx
    nx_graph = nx.Graph()
    if show_singletons:
      nx_graph.add_nodes_from(self.members)
    else:
      nx_graph.add_nodes_from([vert for vert in self.members if vert.edges])

    logging.debug("Vertices added to graph.")
    for edge in self.edges:
      nx_graph.add_edge(edge.vertex_a, edge.vertex_b, {'n_conn': edge.weight,
                                                       'residual': np.mean(np.abs(edge.residuals()))})
    logging.debug("Edges added to graph.")

    mean_I = np.mean([node.total_i for node in self.members])
    abs_residuals = [edge[2]["residual"]
                     for edge in nx_graph.edges_iter(data=True)]
    if edge_width:
      edge_widths = [edge[2]["n_conn"]
                     for edge in nx_graph.edges_iter(data=True)]
    else:
      edge_widths = [1 for edge in nx_graph.edges_iter()]

    if not kwargs.has_key('pos') :
      kwargs['pos'] = nx.spring_layout(nx_graph, iterations=1000)


    fig = plt.figure(figsize=(12, 9))
    im_vertices = nx.draw_networkx_nodes(nx_graph, with_labels=False,
                                         node_color='0.7', **kwargs)
    im_edges = nx.draw_networkx_edges(nx_graph, vmin=0,
                                      vmax=max(abs_residuals),
                                      edge_color=abs_residuals,
                                      edge_cmap=plt.get_cmap("jet"),
                                      width=edge_widths,
                                      **kwargs)
    cb = plt.colorbar(im_edges)
    cmin, cmax = cb.get_clim()
    ticks = np.linspace(cmin, cmax, 7)
    cb.set_ticks(ticks)
    cb.set_ticklabels(['{:.3f}'.format(t) for t in ticks])
    cb.set_label("Edge Residual")
    plt.axis('off')
    plt.tight_layout()
    plt.title("Network Processing\nThickness shows number of connections")
    if saveas is not None:
      plt.savefig(saveas, bbox_inches='tight')
    else:
      plt.show()
    return kwargs['pos']

  def merge_dict(self, estimate_fullies=True):
    """ Make a dict of Miller indices with  ([list of intensities], resolution)
    value tuples for each miller index. Use the fully-corrected equivalent.

    :param estimate_fullies: Boolean for if the partialities and scales should be used.
    """
    intensity_miller = {}
    for vertex in self.members:
      if vertex.edges:
        miller_array = vertex.miller_array
        if estimate_fullies:
          miller_array = miller_array / (vertex.partialities
                                         * vertex.scales)
        #miller_array = miller_array.map_to_asu()
        d_spacings = list(miller_array.d_spacings().data())
        miller_indeces = list(miller_array.indices())
        miller_intensities = list(miller_array.data())

        for observation in zip(miller_indeces, miller_intensities, d_spacings):
          try:
            intensity_miller[observation[0]][0].append(observation[1])
          except KeyError:
            intensity_miller[observation[0]] = ([observation[1]],
                                                observation[2])
    return intensity_miller

  def cc_half(self, nbins=10, estimate_fullies=True):
    """
    Calculate the correlation coefficients for the vertices of the graph.

    :param nbins: number of bins to use for binned correlations.
    :param estimate_fullies: if True, use current orientation, spot profile and scale to estimate the 'full' intensity of each reflection.

    :return cc_half, p_value: overall CC1/2 for all reflections in the cluster.
    :return pretty_string: string with the CC by bin. Lowest resolution bin is extended if data cannot be split into n_bins evenly.
    """
    from scipy.stats import pearsonr
    import operator

    intensity_miller = self.merge_dict(estimate_fullies=estimate_fullies)

    # Remove miller indices that are only measured once.
    multiple_millers = [values for values in intensity_miller.values()
                        if len(values[0]) > 1]

    # Sort, since we will be binning later
    multiple_millers.sort(key=lambda x: x[1])

    # Avoid corner case where number of millers goes exactly into the bins
    if len(multiple_millers) % nbins == 0:
      nbins += 1

    # Figure out bin sizes:
    bin_size = int(len(multiple_millers) / nbins)
    bin_edges = [multiple_millers[i][1] for i in range(0, len(multiple_millers),
                                                       bin_size)]

    # Extend the last bin does not cover the whole resolution range
    bin_edges[-1] = multiple_millers[-1][1]

    assert bin_edges[0] == multiple_millers[0][1]
    assert bin_edges[-1] == multiple_millers[-1][1]
    assert len(bin_edges) == nbins + 1

    # For bins of size bin_size, split each miller index into two, and add half
    # the observations to first_half and half to second_half
    first_half = [[] for _ in range(nbins)]
    second_half = [[] for _ in range(nbins)]
    for indx, intensities in enumerate(multiple_millers):
      # Extending the last bin if necesarry
      bin_id = indx // bin_size if indx // bin_size < nbins else nbins - 1
      obs = intensities[0]
      random.shuffle(obs)
      middle = len(obs) // 2
      first_half[bin_id].append(np.mean(obs[:middle]))
      second_half[bin_id].append(np.mean(obs[middle:]))

    # Calculate the CC for the two halves of each resolution bin.
    logging.info("Calculating CC1/2 by resolution bin. "
                 "{} miller indices per bin".format(bin_size))
    pretty_string = ''
    for bin_id in reversed(range(nbins)):
      bin_cc, bin_p = pearsonr(first_half[bin_id], second_half[bin_id])
      bin_str = "{:6.3f} - {:6.3f}: {:4.3f}".format(bin_edges[bin_id + 1],
                                                    bin_edges[bin_id],
                                                    bin_cc)
      logging.info(bin_str)
      pretty_string += bin_str + "\n"

    # Combine all the split millers, and take the CC of everything
    first_half_all = reduce(operator.add, first_half)
    second_half_all = reduce(operator.add, second_half)
    cc_half, p_value = pearsonr(first_half_all, second_half_all)
    logging.info("CC 1/2 calculated: {:4.3f}, p-value: {}".format(cc_half,
                                                                  p_value))

    return cc_half, p_value, pretty_string

  @staticmethod
  def total_score(params, *args):
    """ return to total residual given params, which is a list of tuples, each
    containing the parameters to refine for one image. Formated for use with
    the scipy minimization routines.
    There are thus 12*N_images parameters, and sum_edges(edge_weight) variables

    :param params: a list of the parameters for the minimisation as a tuple, only including vertices that have edges. See ImageNode.get_x0() for details of parameters.
    :param *args: is the Graph to be minimised, the starting point of params for each new frame, as a list, and a list of the vertex objects that are to be used.
    :return: the total squared residual of the graph minimisation
    """
    assert len(args) == 3, "Must have 2 args: Graph object, and length of each" \
                           "set of parameters."

    graph = args[0]
    param_knots = args[1]
    x0 = args[2]

    # Scale the normalised params up by their starting values.
    params = params * x0

    # 0. break params up into a list of correct length
    param_list = []
    current_pos = 0
    for x0_length in param_knots:
      param_list.append(params[current_pos:current_pos + x0_length])
      current_pos += x0_length
    params = param_list

    # 1. Update scales and partialities for each node that has edges
    # using params.
    for v_id, vertex in enumerate([vertex for vertex in graph.members
                                   if vertex.edges]):
      vertex.partialties = vertex.calc_partiality(params[v_id])
      vertex.scales = vertex.calc_scales(params[v_id])
      vertex.G = params[v_id][0]
      vertex.B = params[v_id][1]

    # 2. Calculate all new residuals
    residuals_by_edge = []
    for edge in graph.edges:
      residuals_by_edge.append(edge.residuals())

    # 3. Return the sum squared of all residuals
    total_sum = 0
    for edge in residuals_by_edge:
      for residual in edge:
        total_sum += residual**2

    logging.debug("Total Score: {}".format(total_sum))
    return total_sum

  def global_minimise(self, nsteps=15000):
    """
    Perform a global minimisation on the total squared residuals on all the
    edges, as calculated by Edge.residuals(). Uses the L-BFGS algorithm.

    :param nsteps: max number of iterations for L-BFGS to use.
    """
    from scipy.optimize import fmin_l_bfgs_b as lbfgs

    # 1. Make a big array of all x0 values
    x0 = []
    param_knots = []
    for vertex in self.members:
      if vertex.edges:
        x0.extend(vertex.get_x0())
        param_knots.append(len(vertex.get_x0()))

    # Test for the above:
    current_pos = 0
    vertices_with_edges = [vert for vert in self.members if vert.edges]
    for im_num, x0_length in enumerate(param_knots):
      these_params = x0[current_pos:current_pos+ x0_length]
      assert all(these_params == vertices_with_edges[im_num].get_x0())
      current_pos += x0_length

    # 2. Do the magic
    final_params, min_total_res, info_dict = lbfgs(Graph.total_score,
                                                   x0,
                                                   approx_grad=True,
                                                   epsilon=0.001,
                                                   args=(self, param_knots,
                                                         np.ones(len(x0))),
                                                   factr=10**12,
                                                   iprint=0,
                                                   disp=10,
                                                   maxiter=nsteps)

    return final_params, min_total_res, info_dict
