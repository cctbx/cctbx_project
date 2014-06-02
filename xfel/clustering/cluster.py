from __future__ import division
import os
import math
import logging
import scipy.cluster.hierarchy as hcluster
import numpy as np
from cctbx.uctbx.determine_unit_cell import NCDist
from xfel.clustering.singleframe import SingleFrame

__author__ = 'zeldin'


class Cluster:
  """Groups single XFEL images (here described by SingleImage objects) as
  cluster objects. You can create a cluster object directly, by using the
  __init__ method, which takes in a list of SingleImage objects, and some
  string-info, or by using a classmethod e.g. to create an object from a folder
  full of integration pickles. SingleImage objects have most of the stuff from
  an integration pickle, but can also have methods to calculate things relating
  to one image at a time.
  Clustering filters can act on these to break them up into cluster objects with
  different members. A 'filter' is just a clustering procedure that puts the
  passes and fails into different clusters. This is acheived through the
  make_sub_cluster() method. This also keeps track of a sub-clusters heritage
  through the .info string, which is appended to. The idea is to be able to
  write filter scripts for each data. e.g:

    test_cluster = Cluster.from_directories(["~/test_data"],
                                          'test_script')
    P3_only = test_cluster.point_group_filer('P3')
    sub_clusters = P3_only.ab_cluster(1200)
    big_cluster = max(sub_clusters, key=lambda x: len(x.members))
    best_data = big_cluster.total_intensity_filter(res=6.5,
                                                   completeness_threshold=0.1,
                                                   plot=False)
    print best_data.info

  cxi.postrefine (or any other merging thing) will be able to be called on the
  ouput of a cluster object method (ToDo)
  """

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
      logging.debug("averaging unit cell {}".format(member.name))
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
          this_frame = SingleFrame(path, filename)
          if hasattr(this_frame, 'name'):
            data.append(this_frame)
          else:
            logging.info('skipping file {}'.format(filename))
    return cls(data, _prefix,
               'Made from files in {}'.format(path_to_integration_dir[:]))


  def make_sub_cluster(self, new_members, new_prefix, new_info):
    """ Make a sub-cluster from a list of SingleFrame objects from the old
    SingleFrame array.
    """
    return Cluster(new_members, new_prefix,
                   ('{}\n{} Next filter {}\n{}\n{} of {} images passed'
                    'on to this cluster').format(
                     self.info, '#' * 30, '#' * 30, new_info,
                     len(new_members), len(self.members)))


  def print_ucs(self):
    """ Prints a list of all the unit cells in the cluster."""
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
    """ Creates a sub-cluster using the highest total intensity images that
          yield a dataset specified by:
          res -- desired resolution. Defaults to that of the dataset.
          completeness -- the desired completeness of the subset
          multiplicity -- the desired multiplicity of the subset
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
      from matplotlib import pyplot as plt

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
      temp_miller_indicies = temp_miller_indicies.concatenate(image,
                                                              assert_is_similar_symmetry=False)
      current_completeness = temp_miller_indicies. \
        merge_equivalents(). \
        array(). \
        completeness()
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
                 linkage_method='single', log=False, max_only=True,
                 ax=None):
    """ Do hierarchical clustering using the Andrews-Berstein distance from
    Andrews & Bernstein J Appl Cryst 47:346 (2014) on the Niggli cells. Returns
    the largest cluster if max_only is true, otherwise a list of clusters. Also
    return a matplotlib axes object for display of a dendogram."""

    logging.info("Hierarchical clustering of unit cells using Andrews-Bernstein"
                 "Distance from Andrews & Bernstein J Appl Cryst 47:346 (2014)")
    import scipy.spatial.distance as dist
    import matplotlib.pyplot as plt

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

    # 4. Plot a dendogram to the axes if no axis is passed, otherwise just
    #    return the axes object
    if ax is None:
      fig = plt.figure("Distance Dendogram")
      ax = fig.gca()
      direct_visualisation = True
    else:
      direct_visualisation = False

    hcluster.dendrogram(this_linkage,
                        labels=[image.name for image in self.members],
                        leaf_font_size=8,
                        color_threshold=threshold, ax=ax)
    if log:
      ax.set_yscale("log")
    else:
      ax.set_ylim(-ax.get_ylim()[1] / 100, ax.get_ylim()[1])

    if direct_visualisation:
      fig.savefig("{}_dendogram.pdf".format(self.cname))
      plt.show()

    return sub_clusters, ax

  def dump_file_list(self, out_file_name=None):
    """ Simply dumps a list of paths to inegration pickle files to a file. One
    line per image
    """
    if out_file_name is None:
      out_file_name = self.cname

    with open("{}.members".format(out_file_name), 'wb') as outfile:
      for i in self.members:
        outfile.write(i.path + "\n")

  def visualise_orientational_distribution(self):
    """ Creates a plot of the orientational distribution of the unit cells.
    """
    from mpl_toolkits.basemap import Basemap
    import matplotlib.pyplot as plt
    import scipy.ndimage as ndi
    from cctbx.array_family import flex

    def cart2sph(x, y, z):
      # cctbx (+z to source, y to ceiling) to
      # lab frame (+x to source, z to ceiling)
      #z, x, y = x, y, z
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
    axis_ids = [0, 1, 2]
    labels = ["a axis in lab frame",
              "b axis in lab frame",
              "c axis in lab frame"]
    plt.figure(figsize=(10, 14))

    for plot_num, axis_id, label in zip(range(1, 4), axis_ids, labels):

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
      ax = plt.subplot(3, 1, plot_num)
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
          rotated_orientation = list(mx.r().as_double() * orientation)
          sym_x, sym_y, sym_lo, sym_la \
            = xy_lat_lon_from_orientation(rotated_orientation, axis_id)
          #assert (sym_x, sym_y) != (main_x, main_y)
          sym_x_coords.append(sym_x)
          sym_y_coords.append(sym_y)
          sym_lon.append(sym_lo)
          sym_lat.append(sym_la)

      # Plot each image as a yellow sphere
      logging.debug(len(x_coords))
      euler_map.plot(x_coords, y_coords, 'oy', markersize=4, markeredgewidth=0.5)

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
      CS = euler_map.contourf(np.array(x_for_plot),
                         np.array(y_for_plot),
                         np.array(local_intensity), tri=True)

      #  Pretty up grap
      cbar = plt.colorbar(CS)
      cbar.ax.set_ylabel('spot density')
      euler_map.drawmeridians(np.arange(0, 360, 60),
                              labels=[0, 0, 1, 0],
                              fontsize=10)
      euler_map.drawparallels(np.arange(-90, 90, 30),
                              labels=[1, 0, 0, 0],
                              fontsize=10)
      plt.title(label, y=1.08)

    plt.show()




#def cluster_pca(self)::w
#  """ Should use BLEND clustering protocols in python (Foaldi et al. Acta D.
#  2013). i.e. filter for parameters that vary, do PCA, then ward linkage
#  clustering on this. """
#  # Will come back to this soon <-- Oli
#  #columns_to_use = [True if np.std(self.niggli_ucs[:, i])
#  #                  else False for i in range(6)]
#@classmethod
#def from_json(cls, json_file, _prefix='cluster_from_json'):
#    """ Does not work!! Do not use! """
#    with open(json_file, 'rb') as inFile:
#      data = json.load(inFile)
#    return cls(data, _prefix,
#               'Made from {}'.format(json_file))
#
