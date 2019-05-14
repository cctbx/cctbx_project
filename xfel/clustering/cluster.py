""" This module is designed to provide tools to deal with groups of serial
crystallography images.

The class Cluster allows the creation, storage and manipulation of these sets
of frames. Methods exist to create sub-clusters (new cluster objects) or to act
on an existing cluster, e.g. to plot the unit cell distributions.

**Author:**   Oliver Zeldin <zeldin@stanford.edu>
"""

from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os
import math
import logging
logger = logging.getLogger(__name__)
from xfel.clustering.singleframe import SingleFrame, SingleDialsFrame, SingleDialsFrameFromFiles
from xfel.clustering.singleframe import SingleDialsFrameFromJson
from cctbx.uctbx.determine_unit_cell import NCDist
import numpy as np

class SingleMiller():
  """ Class for representing multiple measurements of a single miller index. """
  def __init__(self, index, resolution):
    self.index = index
    self.resolution = resolution
    self.intensities = []
    self.sigmas = []

  def add_obs(self, intensity, sigma):
    """ Add an observation of the miller index """
    self.intensities.append(intensity)
    self.sigmas.append(sigma)

  def weighted_mean_and_std(self):
    """ return the mean of the observations, weighted by 1/sigmas.
    :return: [weighted_mean, weighted_std_dev]
    """
    intensities = np.array(self.intensities)
    weights = 1/np.array(self.sigmas)
    w_mean = np.average(intensities, weights=weights)
    w_std = np.sqrt(np.average((intensities - w_mean)**2, weights=weights))
    return w_mean, w_std

  def nobs(self):
    """ return how many observations of index exist """
    return len(self.intensities)



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

  def __init__(self, data, cname="cluster", info=""):
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

    # Calculate medians and stdevs.  Explanation:
    # When clustering cells together, median and standard deviation parameters should always
    # be calculated on the Niggli setting, not the input cell. Admittedly, this does not cover
    # corner cases where the relationship between two Niggli cells crosses a polytope boundary,
    # but for most cases the Niggli settings are the appropriate cells to be averaged.
    # Fixes cctbx issue #97.
    # Prior to 20171121, try the member.orientation.unit_cell(), then
    # member.crystal_symmetry.unit_cell().  After 20171121, use member.uc, which is
    # supposedly always constructed to be the niggli cell.
    unit_cells = np.zeros([len(self.members), 6])
    self.pg_composition = {}
    for i, member in enumerate(self.members):
      unit_cells[i, :] = member.uc # supposed to be the Niggli setting cell parameters tuple
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
                       n_images=None,
                       dials=False,
                       **kwargs):
    """Constructor to get a cluster from pickle files, from the recursively
    walked paths. Can take more than one argument for multiple folders.
    usage: Cluster.from_directories(..)
    :param path_to_integration_dir: list of directories containing pickle files.
    Will be searched recursively.
    :param n_images: find at most this number of images.
    :param use_b: Boolean. If True, intialise Scale and B. If False, use only
    mean intensity scalling.
    """
    if dials:
      dials_refls = []
      dials_expts = []
      pickles = None
      for arg in path_to_integration_dir:
        for (dirpath, dirnames, filenames) in os.walk(arg):
          for filename in filenames:
            path = os.path.join(dirpath, filename)
            if path.endswith("integrated.pickle"):
              dials_refls.append(path)
            elif path.endswith("experiments.json"):
              dials_expts.append(path)

    else:
      pickles = []
      dials_refls = None
      dials_expts = None
      for arg in path_to_integration_dir:
        for (dirpath, dirnames, filenames) in os.walk(arg):
          for filename in filenames:
            path = os.path.join(dirpath, filename)
            if path.endswith(".pickle"):
              print(path, "ends with .pickle")
              pickles.append(path)

    return Cluster.from_files(pickle_list=pickles, dials_refls=dials_refls,
      dials_expts=dials_expts, _prefix=_prefix,
      _message='Made from files in {}'.format(path_to_integration_dir[:]),
      dials=dials, n_images=n_images, **kwargs)

  @classmethod
  def from_crystal_symmetries(cls, crystal_symmetries,
                              lattice_ids=None,
                              _prefix='cluster_from_crystal_symmetries',
                              _message='Made from list of individual cells',
                              n_images=None,
                              dials=False,
                              **kwargs):
    """Constructor to get a cluster from a list of crystal symmetries.
    """

    data = []

    from xfel.clustering.singleframe import CellOnlyFrame
    if lattice_ids is not None:
      assert len(lattice_ids) == len(crystal_symmetries)
    for j, cs in enumerate(crystal_symmetries):
      name = "lattice%07d"%j
      lattice_id = None
      if lattice_ids is not None:
        lattice_id = lattice_ids[j]
      this_frame = CellOnlyFrame(
        crystal_symmetry=cs, path=name, name=name, lattice_id=lattice_id)
      if hasattr(this_frame, 'crystal_symmetry'):
          data.append(this_frame)
      else:
          logger.info('skipping item {}'.format(item))
    logger.info("%d lattices will be analyzed"%(len(data)))

    return cls(data, _prefix, _message)

  @classmethod
  def from_list( cls,file_name,
                 raw_input=None,
                 pickle_list=[],
                 dials_refls=[],
                 dials_expts=[],
                 _prefix='cluster_from_file',
                 _message='Made from list of individual cells',
                 n_images=None,
                 dials=False,
                 **kwargs):
    """Constructor to get a cluster from a single file.  The file must list unit cell a,b,c,alpha,beta,gamma
    and space_group_type, each as a single token.
    :param file_name: pathname of the file
    """

    data = []

    from xfel.clustering.singleframe import CellOnlyFrame
    stream = open(file_name,"r").readlines()
    print("There are %d lines in the input file"%(len(stream)))
    for j,item in enumerate(stream):
      tokens = item.strip().split()
      assert len(tokens) == 7, tokens
      unit_cell_params = tuple([float(t) for t in tokens[0:5]])
      space_group_type = tokens[6]
      from cctbx.uctbx import unit_cell
      uc_init = unit_cell(unit_cell_params)
      from cctbx.sgtbx import space_group_info
      sgi = space_group_info(space_group_type)
      from cctbx import crystal
      crystal_symmetry = crystal.symmetry(unit_cell=uc_init, space_group_info=sgi)
      name = "lattice%07d"%j
      this_frame = CellOnlyFrame(crystal_symmetry, path=name, name=name)
      if hasattr(this_frame, 'crystal_symmetry'):
          data.append(this_frame)
      else:
          logger.info('skipping item {}'.format(item))
    print("%d lattices will be analyzed"%(len(data)))

    return cls(data, _prefix, _message)

  @classmethod
  def from_iterable( cls,iterable,
                     _prefix='cluster_from_iterable',
                     _message='Made from list of individual cells',
                     **kwargs):
    """Constructor to get a cluster from an iterable (a list or tuple).  The
    file must list unit cell a,b,c,alpha,beta,gamma and space_group_type,
    each as a single token.
    :param iterable: a list or a tuple
    """

    data = []
    from xfel.clustering.singleframe import CellOnlyFrame
    from cctbx.uctbx import unit_cell
    from cctbx.sgtbx import space_group_info
    from cctbx import crystal

    for j,item in enumerate(iterable):
      try:
        assert len(item) == 7
        unit_cell_params = tuple([float(t) for t in item[0:6]])
        space_group_type = item[6]
        uc_init = unit_cell(unit_cell_params)
        sgi = space_group_info(space_group_type)
        crystal_symmetry = crystal.symmetry(unit_cell=uc_init, space_group_info=sgi)
        name = "lattice%07d"%j
        this_frame = CellOnlyFrame(crystal_symmetry, path=name, name=name)
        if hasattr(this_frame, 'crystal_symmetry'):
            data.append(this_frame)
      except Exception as e:
        pass

    return cls(data, _prefix, _message)

  @classmethod
  def from_files(cls,
                 raw_input=None,
                 pickle_list=[],
                 dials_refls=[],
                 dials_expts=[],
                 _prefix='cluster_from_file',
                 _message='Made from list of individual files',
                 n_images=None,
                 dials=False,
                 json=False,
                 **kwargs):
    """Constructor to get a cluster from a list of individual files.
    :param pickle_list: list of pickle files
    :param dials_refls: list of DIALS integrated reflections
    :param dials_expts: list of DIALS experiment jsons
    :param n_images: find at most this number of images
    :param dials: use the dials_refls and dials_expts arguments to construct the clusters (default: False)
    :param use_b: Boolean. If True, intialise Scale and B. If False, use only
    mean intensity scalling.
    """

    data = []

    def sort_dials_raw_input(raw):
      expts = []
      refls = []
      for path in raw:
        if path.endswith(".pickle"):
          refls.append(path)
        elif path.endswith(".json"):
          expts.append(path)
      return (refls, expts)

    def done():
      if n_images is None:
        return False
      return len(data) >= n_images

    if dials:
      if raw_input is not None:
        r, e = sort_dials_raw_input(raw_input)
        dials_refls.extend(r)
        dials_expts.extend(e)
      for r, e in zip(dials_refls, dials_expts):
        this_frame = SingleDialsFrameFromFiles(refls_path=r, expts_path=e, **kwargs)
        if hasattr(this_frame, 'miller_array'):
          data.append(this_frame)
          if done():
            break
        else:
          logger.info('skipping reflections {} and experiments {}'.format(r, e))
    elif json:
      if raw_input is not None:
        r, e = sort_dials_raw_input(raw_input)
        dials_expts.extend(e)
      dials_expts_ids = [os.path.join(os.path.dirname(e), os.path.basename(e).split("_")[0])
                         for e in dials_expts]
      for e in dials_expts:
        name = os.path.join(os.path.dirname(e), os.path.basename(e).split("_")[0])
        this_frame = SingleDialsFrameFromJson(expts_path=e,  **kwargs)
        this_frame.name=name
        data.append(this_frame)
        if done():
            break
    else:
      if raw_input is not None:
        pickle_list.extend(raw_input)
      print("There are %d input files"%(len(pickle_list)))
      from xfel.command_line.print_pickle import generate_data_from_streams
      for data_dict in generate_data_from_streams(pickle_list):
        this_frame = SingleFrame(dicti=data_dict, **kwargs)
        if hasattr(this_frame, 'miller_array'):
          data.append(this_frame)
          if done():
            break
        else:
          logger.info('skipping file {}'.format(os.path.basename(path)))
      print("%d lattices will be analyzed"%(len(data)))

    return cls(data, _prefix, _message)

  @classmethod
  def from_expts(cls,
                 refl_table=None,
                 expts_list=None,
                 _prefix='cluster_from_file',
                 _message='Made from experiment objects',
                 n_images=None,
                 **kwargs):
    """Constructor to get a cluster from experiment and reflection list objects
    :param refl_table: DIALS integrated reflection table
    :param expts_list: DIALS experiment list
    :param n_images: find at most this number of images
    """

    data = []

    def done():
      if n_images is None:
        return False
      return len(data) >= n_images

    for i, expt in enumerate(expts_list):
      sel = refl_table['id'] == i
      refls_sel = refl_table.select(sel)
      this_frame = SingleDialsFrame(refl=refls_sel, expt=expt, id=i, **kwargs)
      if hasattr(this_frame, 'miller_array'):
        data.append(this_frame)
        if done():
          break
      else:
        logger.info('skipping invalid experiment #{}'.format(i))

    return cls(data, _prefix, _message)


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

  def best_by_CC(self, other, assert_is_similar_symmetry=False):
    """ Return the SingleFrame object with the highest CC to a reference miller array.
    :param other: miller array object to be correlated against
    :return: a SingleFrame object.
    """
    max = 0
    for sf in self.members:
      corr =  sf.miller_array.correlation(other,
                      assert_is_similar_symmetry=assert_is_similar_symmetry)
      if corr > max:
        best = sf
    return best


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
    logger.info(("Performing intensity filtering, aiming for {}% overall "
                  "completenes at {} A resolution").format(
      completeness_threshold * 100, res))

    # 0. Check that the cluster has consistent point_group (completness doesn't
    #  mean much otherwise...
    assert all(i.pg == self.members[0].pg for i in self.members)

    # 1. Sort SingleFrames by total intensity
    sorted_cluster = sorted(self.members, key=lambda y: -1 * y.total_i)

    if plot:
      import matplotlib.pyplot as plt
      plt.plot([x.total_i for x in sorted_cluster])
      plt.show()

    if res == '':
      res = sorted_cluster[0].d_min()  # Use the high-res limit from the
      # brightest image. ToDo: make this better
      logger.warning(("no resolution limit specified, using the res limit of"
                       "the top-rankeed image: {} A").format(res))

    # 2. Incrementally merge frames until criterion are matched

    temp_miller_indicies = sorted_cluster[0].miller_array
    for idx, image in enumerate((x.miller_array for x in sorted_cluster[1:])):
      temp_miller_indicies = temp_miller_indicies. \
        concatenate(image, assert_is_similar_symmetry=False)
      current_completeness = temp_miller_indicies.merge_equivalents() \
                                                .array() \
                                                .completeness()
      logger.debug(
        "{} images: {:.2f}% complete".format(idx, current_completeness * 100))
      if current_completeness <= completeness_threshold:
        temp_miller_indicies.concatenate(image,
                                         assert_is_similar_symmetry=False)
        if idx + 1 == len(sorted_cluster[1:]):
          logger.warning("Desired completeness could not be acheived, sorry.")
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
                 ax=None, write_file_lists=True, schnell=False, doplot=True,
                 labels='default'):
    """
    Hierarchical clustering using the unit cell dimentions.

    :param threshold: the threshold to use for prunning the tree into clusters.
    :param method: which clustering method from scipy to use when creating the tree (see scipy.cluster.hierarchy)
    :param linkage_method: which linkage method from scipy to use when creating the linkages. x (see scipy.cluster.hierarchy)
    :param log: if True, use log scale on y axis.
    :param ax: if a matplotlib axes object is provided, plot to this. Otherwise, create a new axes object and display on screen.
    :param write_file_lists: if True, write out the files that make up each cluster.
    :param schnell: if True, use simple euclidian distance, otherwise, use Andrews-Berstein distance from Andrews & Bernstein J Appl Cryst 47:346 (2014) on the Niggli cells.
    :param doplot: Boolean flag for if the plotting should be done at all.
    Runs faster if switched off.
    :param labels: 'default' will not display any labels for more than 100 images, but will display file names for fewer. This can be manually overidden with a boolean flag.
    :return: A list of Clusters ordered by largest Cluster to smallest

    .. note::
      Use 'schnell' option with caution, since it can cause strange behaviour
      around symmetry boundaries.
    """

    logger.info("Hierarchical clustering of unit cells")
    import scipy.spatial.distance as dist
    import scipy.cluster.hierarchy as hcluster

    # 1. Create a numpy array of G6 cells
    g6_cells = np.array([SingleFrame.make_g6(image.uc)
                         for image in self.members])

    # 2. Do hierarchichal clustering, using the find_distance method above.
    if schnell:
      logger.info("Using Euclidean distance")
      pair_distances = dist.pdist(g6_cells, metric='euclidean')
      metric = 'euclidean'
    else:
      logger.info("Using Andrews-Bernstein distance from Andrews & Bernstein "
                   "J Appl Cryst 47:346 (2014)")
      pair_distances = dist.pdist(g6_cells,
                                metric=lambda a, b: NCDist(a, b))
      metric = lambda a, b: NCDist(a, b)
    if len(pair_distances) > 0:
      logger.info("Distances have been calculated")
      this_linkage = hcluster.linkage(pair_distances,
                                      method=linkage_method,
                                      metric=metric)
      cluster_ids = hcluster.fcluster(this_linkage,
                                      threshold,
                                      criterion=method)
      logger.debug("Clusters have been calculated")
    else:
      logger.debug("No distances were calculated. Aborting clustering.")
      return [], None

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
      import matplotlib.pyplot as plt
      if labels is True:
        labels = [image.name for image in self.members]
      elif labels is False:
        labels = ['' for _ in self.members]
      elif labels == 'default':
        if len(self.members) > 100:
          labels = ['' for _ in self.members]
        else:
          labels = [image.name for image in self.members]
      else:
         labels = [getattr(v, labels, '') for v in self.members]

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
                          p=200,
                          truncate_mode='lastp', # show only the last p merged clusters
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
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as patheffects
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
      logger.debug("axis_id: {}".format(axis_id))
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
      logger.debug(len(x_coords))
      euler_map.plot(x_coords, y_coords, 'oy',
                     markersize=4,
                     markeredgewidth=0.5)

      # Plot the symetry mates as black crosses
      #euler_map.plot(sym_x_coords, sym_y_coords, 'kx')

      # Use a histogram to bin the data in lattitude/longitude space, smooth it,
      # then plot this as a contourmap. This is for all the symetry-related
      # copies
      #density_hist = np.histogram2d(lat + sym_lat, lon + sym_lon,
      #                                    bins=[range(-90, 91), range(0, 361)])
      # No symmetry mates until we can verify what the cctbx libs are doing
      density_hist = np.histogram2d(lat, lon,
                                    bins=[list(range(-90, 91)), list(range(0, 361))])
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
    import matplotlib.pyplot as plt
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
    axes_to_return[0].set_title("Standard Errors on Wilson fit")
    axes_to_return[0].set_ylabel("Count")
    axes_to_return[0].set_xlabel("Standard Error [$\AA^2$]")

    rs = [-1 * i.minus_2B / 2 for i in self.members]
    axes_to_return[1].hist(rs, 50, range=[-50, 200])
    axes_to_return[1].set_title("B values for Wilson plot")
    axes_to_return[1].set_ylabel("Count")
    axes_to_return[1].set_xlabel(r"B [$\AA^2$]")

    axes_to_return[2].plot([i.G for i in self.members],
             [-1 * i.minus_2B / 2 for i in self.members], 'x')
    axes_to_return[2].set_xlabel("G [AU]")
    axes_to_return[2].set_ylabel("B [$\AA^2$]")
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
    import matplotlib.pyplot as plt

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
    fit_info = "G: {:.2f}, -2B: {:.2f}, r: {:.2f}, std_err: {:.2f}".format(G, minus_2B,
                                                            r_val, std_err)
    smooth = Sf._moving_average(log_i, n=smoothing_width)
    ax.plot(one_over_d_square, log_i, 'bo', ms=1)
    ax.plot(one_over_d_square[smoothing_width - 1:], smooth,'--r', lw=2)
    plt.xlim([0, max(one_over_d_square)])
    ax.plot([0, -1 * G / minus_2B], [G, 0], 'y-', lw=2)
    plt.xlabel(r"$(sin(\theta)/\lambda)^2 [\AA^{-2}]$")
    plt.ylabel("ln(I)")
    plt.title("Simple Wilson fit\n{}".format(fit_info))
    plt.tight_layout()

    if direct_visualisation:
      fig.savefig("{}_dendogram.pdf".format(self.cname))
      plt.show()

    return ax

  def merge_dict(self, use_fullies=False):
    """ Make a dict of Miller indices with  ([list of intensities], resolution)
    value tuples for each miller index.
    """

    miller_dict = {}
    for m in self.members:
      # Use fullies if requested
      if use_fullies:
        if m.miller_fullies:
          miller_array = m.miller_fullies
        else:
          logger.warning("Fully recorded array has not been calculated")
      else:
        miller_array = m.miller_array

      # Create a dictionairy if SingleMiller observations
      d_spacings = list(miller_array.d_spacings().data())
      miller_indeces = list(miller_array.indices())
      miller_intensities = list(miller_array.data())
      miller_sigmas = list(miller_array.sigmas())

      for observation in zip(miller_indeces,
                             miller_intensities,
                             miller_sigmas,
                             d_spacings):
        try:
          miller_dict[observation[0]].add_obs(observation[1],
                                              observation[2])
        except KeyError:
          miller_dict[observation[0]] = SingleMiller(observation[0],
                                                     observation[3])
          miller_dict[observation[0]].add_obs(observation[1],
                                              observation[2])
    return miller_dict

  def __len__(self):
    """ Number of images in the cluster """
    return len(self.members)

  def dump_as_mtz(self, mtz_name, use_fullies=False):
    """ Merge using weighted mean and standard deviation if all miller arrays """
    from cctbx.crystal import symmetry
    from cctbx import miller

    assert all((str(m.miller_array.space_group_info())  \
                == str(self.members[0].miller_array.space_group_info())
                for m in self.members)),  \
                "All images must be in the same point group!"

    final_sym = symmetry(unit_cell=self.medians,
              space_group_info=self.members[0].miller_array.space_group_info())
    # Find mean_Iobs
    mil_dict = self.merge_dict(use_fullies=use_fullies)
    single_millers = mil_dict.values()
    indices = [md.index for md in single_millers]
    iobs, sig_iobs = zip(*[md.weighted_mean_and_std() for md in single_millers])
    all_obs = miller.array(miller_set=self.members[0] \
                                          .miller_array \
                                          .customized_copy(
                                            crystal_symmetry=final_sym,
                                            indices=flex.miller_index(indices),
                                            unit_cell=self.medians),
                                            data = flex.double(iobs),
                                            sigmas = flex.double(sig_iobs))
    all_obs = all_obs.set_observation_type_xray_intensity() \
                     .average_bijvoet_mates()

    all_obs = all_obs.select(all_obs.data() > 0)
    mtz_out = all_obs.as_mtz_dataset(column_root_label="Iobs",
                                         title=self.cname,
                                         wavelength=np.median(
                                           [m.wavelength for m in self.members]))
    mtz_out.add_miller_array(miller_array=all_obs,
                                 column_root_label="IMEAN")
    mtz_obj = mtz_out.mtz_object()
    mtz_obj.write(mtz_name)

    logger.info("MTZ file written to {}".format(mtz_name))

  def to_pandas(self):
    import pandas as pd
    return pd.DataFrame({s['name']: s
                         for s
                         in [m.to_panda() for m in self.members]}).T

  def uc_feature_vector(self):
    """ Return a len(cluster) * 6 numpy array of features for use in ML algos"""
    ucs = [c.uc for c in self.members]
    return  np.array(zip(*ucs))


  def prime_postrefine(self, inputfile):
    """ Run postrefinement from Prime on a cluster. Implements the prime API, and updates the SingleFrame objects with the attribute `miller_fullies`.
    :param cluster: a cluster object.
    :param inputfile: a Prime .inp file.
    """
    from .api import refine_many
    miller_fullies = refine_many(self.members, inputfile)
    for mil, sf in zip(miller_fullies, self.members):
      sf.miller_fullies = mil
