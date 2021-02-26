from __future__ import absolute_import, division, print_function
import copy
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

  def plot_after_optimize(self):
          print ("optimized coordinates", self.coords.focus())
          xx = []
          yy = []
          for item in range(self.coords.focus()[0]):
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
          # do another coords plot, but color by cluster
          # but the problem with this is I don't know to which cluster id==0['h,k,l'] belongs
          xx = flex.double()
          yy = flex.double()
          for item in range(self.coords.focus()[0]):
            xx.append(self.coords[(item,0)])
            yy.append(self.coords[(item,1)])
          from matplotlib import pyplot as plt
          for cl_id in [0.0,1.0]:
            xs = []
            ys = []
            for k in range(len(xx)):
              if self.cluster_labels[k]==cl_id:
                xs.append(xx[k])
                ys.append(yy[k])
            plt.plot(xs,ys,{0.0:"b.",1.0:"r."}[cl_id])
          plt.axes().set_aspect("equal")
          circle = plt.Circle((0,0),1,fill=False,edgecolor="b")
          ax = plt.gca()
          ax.add_artist(circle)
          plt.show()

  def _space_group_for_dataset(self, dataset_id, sym_ops):
      """
      Code from @rjgildea originally in DIALS commit 60986de. This was refactored
      out of DIALS so for now we are reproducing it here.
      """
      if self.input_space_group is not None:
          sg = copy.deepcopy(self.input_space_group)
      else:
          sg = sgtbx.space_group()
      ref_sym_op_id = None
      ref_cluster_id = None
      for sym_op_id in range(len(sym_ops)):
          i_cluster = self.cluster_labels[
              len(self.input_intensities) * sym_op_id + dataset_id
          ]
          if i_cluster < 0:
              continue
          if ref_sym_op_id is None:
              ref_sym_op_id = sym_op_id
              ref_cluster_id = i_cluster
              continue
          op = sym_ops[ref_sym_op_id].inverse().multiply(sym_ops[sym_op_id])
          if i_cluster == ref_cluster_id:
              sg.expand_smx(op.new_denominators(1, 12))
      return sg.make_tidy()


  def run(self): # specializes dials.algorithms.symmetry.cosym.CosymAnalysis
        #from libtbx.development.timers import Profiler
        #P = Profiler("initialize target W")
        self._intialise_target()

        #P = Profiler("optimize W")
        max_calls = self.params.minimization.max_calls
        self._optimise(
            self.params.minimization.engine,
            max_iterations=self.params.minimization.max_iterations,
            max_calls=min(20, max_calls) if max_calls else max_calls,
        )
        #self.plot_after_optimize()

        #P = Profiler("analyze symmetry W")
        self._analyse_symmetry()
        #print(str(self._symmetry_analysis))

        #P = Profiler("cluster analysis W")
        self._cluster_analysis()
        #self.plot_after_cluster_analysis()

  @Subject.notify_event(event="analysed_clusters")
  def _cluster_analysis(self):

        if self.params.cluster.n_clusters == 1:
            self.cluster_labels = flex.double(self.coords.all()[0])
        else:
            self.cluster_labels = self._do_clustering(self.params.cluster.method)

        sym_ops = [
            sgtbx.rt_mx(s).new_denominators(1, 12) for s in self.target.get_sym_ops()
        ]

        reindexing_ops = {}
        space_groups = {}
        self.cosets = {}

        for dataset_id in range(len(self.input_intensities)):
            space_groups[dataset_id] = self._space_group_for_dataset(
                dataset_id, sym_ops
            )

            cosets = sgtbx.cosets.left_decomposition(
                self.target._lattice_group, space_groups[dataset_id]
            )
            self.cosets[dataset_id] = cosets

            reindexing_ops[dataset_id] = self._reindexing_ops_for_dataset(
                dataset_id, sym_ops, cosets
            )

        self.space_groups = space_groups
        self.reindexing_ops = reindexing_ops




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
    def __init__(self, experiments, reflections, uuid_cache_in, params=None):
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
        self.cosym_analysis = CosymAnalysis(datasets, self.params)
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

    @Subject.notify_event(event="run_cosym")
    def run(self):
        self.cosym_analysis.run()

        space_groups = {}
        reindexing_ops = {}
        for dataset_id in self.cosym_analysis.reindexing_ops:
            if 0 in self.cosym_analysis.reindexing_ops[dataset_id]:
                cb_op = self.cosym_analysis.reindexing_ops[dataset_id][0]
                reindexing_ops.setdefault(cb_op, [])
                reindexing_ops[cb_op].append(dataset_id)
            if dataset_id in self.cosym_analysis.space_groups:
                space_groups.setdefault(
                    self.cosym_analysis.space_groups[dataset_id], []
                )
                space_groups[self.cosym_analysis.space_groups[dataset_id]].append(
                    dataset_id
                )

        logger.info("Space groups:")
        for sg, datasets in space_groups.items():
            logger.info(str(sg.info().reference_setting()))
            logger.info(datasets)

        logger.info("Reindexing operators:")
        for cb_op, datasets in reindexing_ops.items():
            logger.info(cb_op)
            logger.info(datasets)

        self._apply_reindexing_operators(
            reindexing_ops, subgroup=None
        # changed in xfel subclass:  set subgroup=None.
        # in dials parent class:  subgroup=self.cosym_analysis.best_subgroup
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
