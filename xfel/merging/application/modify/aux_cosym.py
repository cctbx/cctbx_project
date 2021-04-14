from __future__ import absolute_import, division, print_function
import copy
import os
from scitbx.array_family import flex
from dials.util.observer import Subject
from cctbx import sgtbx, miller
from libtbx import easy_mp
from scipy import sparse
import numpy as np

# Specialization, run only a subset of cosym steps and include plot
from dials.algorithms.symmetry.cosym import CosymAnalysis as BaseClass
from dials.util.multi_dataset_handling import select_datasets_on_identifiers

from dials.algorithms.symmetry.cosym.target import Target

class CosymAnalysis(BaseClass):

  def __init__(self, *args, **kwargs):
    self.do_plot = kwargs.pop('do_plot', False)
    self.i_plot = kwargs.pop('i_plot', None)
    self.plot_fname = kwargs.pop('plot_fname', None)
    self.plot_format = kwargs.pop('plot_format', None)
    self.output_dir = kwargs.pop('output_dir', None)
    super(CosymAnalysis, self).__init__(*args, **kwargs)

  def run(self):
    super(CosymAnalysis, self).run()
    if self.do_plot: self.plot_after_cluster_analysis()

  def plot_after_optimize(self):
          print ("optimized coordinates", self.coords.shape)
          xx = []
          yy = []
          for item in range(self.coords.shape[0]):
            xx.append(self.coords[(item,0)])
            yy.append(self.coords[(item,1)])
          from matplotlib import pyplot as plt
          plt.plot(xx,yy,"r.")
          # denminator of 12 is specific to the use case of P6 (# symops in the metric superlattice)
          plt.plot(xx[::len(xx)//12],yy[::len(yy)//12],"b.")
          plt.plot(xx[:1],yy[:1],"g.")
          plt.axes().set_aspect("equal")
          circle = plt.Circle((0,0),1,fill=False,edgecolor="b")
          ax = plt.gca()
          ax.add_artist(circle)
          plt.show()

  def plot_after_cluster_analysis(self):
          xx = flex.double()
          yy = flex.double()
          for item in range(self.coords.shape[0]):
            xx.append(self.coords[(item,0)])
            yy.append(self.coords[(item,1)])
          from matplotlib import pyplot as plt
          plt.plot(xx, yy, 'r.')
          ax = plt.gca()
          ax.set_aspect("equal")
          circle = plt.Circle((0,0),1,fill=False,edgecolor="b")
          ax.add_artist(circle)
          if self.plot_fname is None:
            plt.show()
          else:
            plot_path = os.path.join(self.output_dir, self.plot_fname)
            plot_fname = "{}_{}.{}".format(
                plot_path, self.i_plot, self.plot_format)
            plt.savefig(plot_fname)


from dials.command_line.cosym import logger
from dials.command_line.cosym import cosym as dials_cl_cosym_wrapper
from dials.util.exclude_images import get_selection_for_valid_image_ranges
from dials.command_line.symmetry import (
    apply_change_of_basis_ops,
    change_of_basis_ops_to_minimum_cell,
    eliminate_sys_absent,
)
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections

class dials_cl_cosym_subclass (dials_cl_cosym_wrapper):
    def __init__(self, experiments, reflections, uuid_cache_in, params=None,
            do_plot=False, i_plot=None, output_dir=None):
        super(dials_cl_cosym_wrapper, self).__init__(
            events=["run_cosym", "performed_unit_cell_clustering"]
        )
        if params is None:
            params = phil_scope.extract()
        self.params = params

        self._reflections = []
        for refl, expt in zip(reflections, experiments):
            sel = get_selection_for_valid_image_ranges(refl, expt)
            self._reflections.append(refl.select(sel))

        self._experiments, self._reflections = self._filter_min_reflections(
            experiments, self._reflections, uuid_cache_in
        )
        self.ids_to_identifiers_map = {}
        for table in self._reflections:
            self.ids_to_identifiers_map.update(table.experiment_identifiers())
        self.identifiers_to_ids_map = {
            value: key for key, value in self.ids_to_identifiers_map.items()
        }

        if len(self._experiments) > 1:
            # perform unit cell clustering
            identifiers = self._unit_cell_clustering(self._experiments)
            if len(identifiers) < len(self._experiments):
                logger.info(
                    "Selecting subset of %i datasets for cosym analysis: %s"
                    % (len(identifiers), str(identifiers))
                )
                self._experiments, self._reflections = select_datasets_on_identifiers(
                    self._experiments, self._reflections, use_datasets=identifiers
                )

        # Map experiments and reflections to minimum cell
        cb_ops = change_of_basis_ops_to_minimum_cell(
            self._experiments,
            params.lattice_symmetry_max_delta,
            params.relative_length_tolerance,
            params.absolute_angle_tolerance,
        )
        exclude = [
            expt.identifier
            for expt, cb_op in zip(self._experiments, cb_ops)
            if not cb_op
        ]
        if len(exclude):
            logger.info(
                "Rejecting {} datasets from cosym analysis "\
                "(couldn't determine consistent cb_op to minimum cell):\n"\
                "{}".format(len(exclude), exclude)
            )
            self._experiments, self._reflections = select_datasets_on_identifiers(
                self._experiments, self._reflections, exclude_datasets=exclude
            )
            cb_ops = list(filter(None, cb_ops))

        # Eliminate reflections that are systematically absent due to centring
        # of the lattice, otherwise they would lead to non-integer miller indices
        # when reindexing to a primitive setting
        self._reflections = eliminate_sys_absent(self._experiments, self._reflections)

        self._experiments, self._reflections = apply_change_of_basis_ops(
            self._experiments, self._reflections, cb_ops
        )

        # transform models into miller arrays
        datasets = filtered_arrays_from_experiments_reflections(
            self.experiments,
            self.reflections,
            outlier_rejection_after_filter=False,
            partiality_threshold=params.partiality_threshold,
        )

        datasets = [
            ma.as_anomalous_array().merge_equivalents().array() for ma in datasets
        ]

        # opportunity here to subclass as defined above, instead of the dials-implemented version
        self.cosym_analysis = CosymAnalysis(
            datasets,
            self.params,
            do_plot=do_plot,
            i_plot=i_plot,
            plot_fname=self.params.plot.filename,
            plot_format=self.params.plot.format,
            output_dir=output_dir
            )
        #Fixed in subclass: parent class apparently erases the knowledge of input-to-minimum cb_ops.
        # without storing the op in self, we can never trace back to input setting.
        self.cb_op_to_minimum = cb_ops

        #Not sure yet, we may be assuming that all the cb_ops are the same (applicable for PSI with P63)
        assertion_set = set(cb_ops)
        assert len(assertion_set)==1 # guarantees all are the same; revisit with different use cases later

    def _filter_min_reflections(self, experiments, reflections, uuid_cache_in):
        identifiers = []
        self.uuid_cache = []
        for expt, refl, uuid in zip(experiments, reflections, uuid_cache_in):
            if len(refl) >= self.params.min_reflections:
                identifiers.append(expt.identifier)
                self.uuid_cache.append(uuid)

        return select_datasets_on_identifiers(
            experiments, reflections, use_datasets=identifiers
        )


class TargetWithFastRij(Target):
  def __init__(self, *args, **kwargs):
    test_data_path = kwargs.pop('test_data_path', None)
    if test_data_path is None:
      super(TargetWithFastRij, self).__init__(*args, **kwargs)
      return
    else:
      # This is only for unit testing
      import pickle
      import numpy as np
      self._nproc = 1
      with open(test_data_path, 'rb') as f:
        self._lattices = np.array(pickle.load(f))
        self.sym_ops = pickle.load(f)
        self._weights = pickle.load(f)
        self._data = pickle.load(f)
        self._patterson_group = pickle.load(f)
        self._min_pairs = 3 # minimum number of mutual miller indices for a match

        # truncate the input data to save time
        self._lattices = self._lattices[:10]
        i_last = self._lattices[-1]
        self._data = self._data[:i_last]

  def compute_gradients(self, x):
    grad = super(TargetWithFastRij, self).compute_gradients(x)
    return grad.A1

  def _compute_rij_wij(self, use_cache=True, use_super=False):

    if use_super:
      # for testing
      return super()._compute_rij_wij(use_cache=use_cache)

    def _compute_one_row(args):
      """
      Call the compute_one_row method of CC, which is a compute_rij_wij_detail
      instance. Args is a tuple that we unpack because easy_mp.parallel_map
      can only pass a single argument.
      """
      CC, i, NN = args
      rij_row, rij_col, rij_data, wij_row, wij_col, wij_data = [
          list(x) for x in CC.compute_one_row(self._lattices.size, i)
      ]
      rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
      wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))
      return rij, wij

    n_lattices = self._lattices.size
    n_sym_ops = len(self.sym_ops)
    NN = n_lattices * n_sym_ops
    lower_i = flex.int()
    upper_i = flex.int()
    for lidx in range(self._lattices.size):
      LL,UU = self._lattice_lower_upper_index(lidx)
      lower_i.append(int(LL))
      if UU is None:  UU = self._data.data().size()
      upper_i.append(int(UU))
    indices = {}
    space_group_type = self._data.space_group().type()
    from xfel.merging import compute_rij_wij_detail
    CC = compute_rij_wij_detail(
        lower_i,
        upper_i,
        self._data.data(),
        self._min_pairs)
    for sym_op in self.sym_ops:
        cb_op = sgtbx.change_of_basis_op(sym_op)
        indices_reindexed = cb_op.apply(self._data.indices())
        miller.map_to_asu(space_group_type, False, indices_reindexed)
        indices[cb_op.as_xyz()] = indices_reindexed
        CC.set_indices(cb_op, indices_reindexed)

    results = easy_mp.parallel_map(
        _compute_one_row,
        [(CC, i, NN) for i in range(n_lattices)],
        processes=self._nproc
    )

    rij_matrix = None
    wij_matrix = None
    for (rij, wij) in results:
        if rij_matrix is None:
            rij_matrix = rij
            wij_matrix = wij
        else:
            rij_matrix += rij
            wij_matrix += wij
    self.rij_matrix = rij_matrix.todense().astype(np.float64)
    self.wij_matrix = wij_matrix.todense().astype(np.float64)
    return self.rij_matrix, self.wij_matrix
