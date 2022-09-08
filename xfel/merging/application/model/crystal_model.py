from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
from cctbx import miller
from cctbx.crystal import symmetry
import iotbx.pdb
import mmtbx.model
from six.moves import cStringIO as StringIO
import os

class crystal_model(worker):

  def __repr__(self):
    return 'Build crystal model'

  def __init__(self, params, purpose, mpi_helper=None, mpi_logger=None):
    super(crystal_model, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)
    self.purpose = purpose

  def run(self, experiments, reflections):
    '''Create a model for the crystal'''
    self.logger.log_step_time("CREATE_CRYSTAL_MODEL")

    # perform some run-time validation
    assert self.purpose in ["scaling", "statistics", "cosym"]
    if self.purpose == "statistics" and self.params.merging.set_average_unit_cell: # without this flag the average unit cell wouldn't have been generated
      assert 'average_unit_cell' in (self.params.statistics).__dict__ # the worker, tasked with averaging the unit cell, would put it there

    # Generate the reference, or model, intensities
    i_model = None
    if self.purpose in ["scaling","cosym"]:
      model_file_path = str(self.params.scaling.model)
      if model_file_path is not None:
        self.logger.log("Scaling model: " + model_file_path)
        if self.mpi_helper.rank == 0:
          self.logger.main_log("Scaling model: " + model_file_path)
      else:
        self.logger.log("No scaling model has been provided")
        if self.mpi_helper.rank == 0:
          self.logger.main_log("No scaling model has been provided")
    elif self.purpose == "statistics":
      model_file_path = str(self.params.statistics.cciso.mtz_file)
      if model_file_path is not None:
        self.logger.log("Reference for statistics: " + model_file_path)
        if self.mpi_helper.rank == 0:
          self.logger.main_log("Reference for statistics: " + model_file_path)
      else:
        self.logger.log("No reference for statistics has been provided")
        if self.mpi_helper.rank == 0:
          self.logger.main_log("No reference for statistics has been provided")

    if model_file_path is not None:
      if model_file_path.endswith(".mtz"):
        i_model = self.create_model_from_mtz(model_file_path)
      elif model_file_path.endswith(".pdb"):
        i_model = self.create_model_from_pdb(model_file_path)
      elif model_file_path.endswith(".cif"):
        i_model = self.create_model_from_structure_file(model_file_path)

    if self.purpose == "cosym":
      return i_model

    # Generate a full miller set. If the purpose is scaling and i_model is available, then the miller set has to be consistent with the model
    miller_set, i_model = self.consistent_set_and_model(i_model=i_model)

    # Save i_model and miller_set to the parameters
    if i_model is not None:
      if not 'i_model' in self.params.scaling.__dict__:
        self.params.scaling.__inject__('i_model', i_model)
      else:
        self.params.scaling.__setattr__('i_model', i_model)

    if not 'miller_set' in self.params.scaling.__dict__:
      self.params.scaling.__inject__('miller_set', miller_set)
    else:
      self.params.scaling.__setattr__('miller_set', miller_set)

    self.logger.log_step_time("CREATE_CRYSTAL_MODEL", True)

    # Add asymmetric unit indexes to the reflection table
    if self.purpose == "scaling":
      self.logger.log_step_time("ADD_ASU_HKL_COLUMN")
      self.add_asu_miller_indices_column(reflections)
      self.logger.log_step_time("ADD_ASU_HKL_COLUMN", True)

    return experiments, reflections

  def create_model_from_pdb(self, model_file_path):
      return self.create_model_from_structure_file(model_file_path)

  def create_model_from_structure_file(self, model_file_path):

    if self.mpi_helper.rank == 0:
      model_ext = os.path.splitext(model_file_path)[-1].lower()
      assert model_ext in [".pdb", ".cif"]
      if model_ext == ".pdb":
        from iotbx import file_reader
        pdb_in = file_reader.any_file(model_file_path, force_type="pdb")
        pdb_in.assert_file_type("pdb")
        xray_structure = pdb_in.file_object.xray_structure_simple()
      elif model_ext == ".cif":
        from libtbx.utils import Sorry
        try:
          from cctbx.xray import structure
          xs_dict = structure.from_cif(file_path=model_file_path)
          assert len(xs_dict) == 1, "CIF should contain only one xray structure"
          xray_structure = list(xs_dict.values())[0]
        except Sorry:
          inp = iotbx.pdb.input(model_file_path)
          model = mmtbx.model.manager(model_input=inp)
          xray_structure = model.get_xray_structure()
      out = StringIO()
      xray_structure.show_summary(f=out)
      self.logger.main_log(out.getvalue())
      space_group = xray_structure.crystal_symmetry().space_group().info()
      unit_cell = xray_structure.crystal_symmetry().unit_cell()
    else:
      xray_structure = space_group = unit_cell = None
    if self.purpose != "cosym":
      xray_structure, space_group, unit_cell = self.mpi_helper.comm.bcast((xray_structure, space_group, unit_cell), root=0)
    if self.purpose == "scaling":
      # save space group and unit cell as scaling targets
      self.params.scaling.space_group = space_group
      self.params.scaling.unit_cell = unit_cell

    # prepare phil parameters to generate model intensities
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()

    '''
    #
    # RB: for later
    #
    # We need to set the model resolution limits here, which will determine the set of asu HKLs,
    # for which the intensities will be calculated. So our input resolution limits need to be adjusted
    # to account for the difference between the model unit cell and the variations of the unit cells in the experiments.
    # Ideally, we should calculate the model intensity for _every_ observed HKL in every experiment, which might be a time-wasteful calculation.
    # So we use a shortcut: a ratio of the corresponding unt cell volumes as well as a fudge factor - "resolution scalar".

    model_unit_cell_volume = xray_structure.crystal_symmetry().unit_cell().volume()

    min_exp_unit_cell_volume, max_exp_unit_cell_volume = self.get_min_max_experiment_unit_cell_volume(self.experiments)

    assert self.params.scaling.resolution_scalar > 0 and self.params.scaling.resolution_scalar < 1.0 # the way this scalar is used, it should be less than unity

    params2.high_resolution = self.params.merging.d_min / math.pow(max_exp_unit_cell_volume/model_unit_cell_volume, 1./3.) * self.params.scaling.resolution_scalar
    if self.params.merging.d_max is not None:
      params2.low_resolution = self.params.merging.d_max / math.pow(min_exp_unit_cell_volume/model_unit_cell_volume, 1./3.) / self.params.scaling.resolution_scalar
    '''

    #
    # RB: for now we do it the way it's done in cxi.merge. There is a problem with that approach though,
    # because the "filter unit cell" may be different from the "model unit cell",
    # so the cell variations we are trying to account for here might be for a different unit cell.
    #
    # from cxi.merge:
    # adjust the cutoff of the generated intensities to assure that
    # statistics will be reported to the desired high-resolution limit
    # even if the observed unit cell differs slightly from the reference.
    params2.high_resolution = self.params.merging.d_min * self.params.scaling.resolution_scalar
    if self.params.merging.d_max is not None:
      params2.low_resolution = self.params.merging.d_max / self.params.scaling.resolution_scalar

    params2.output.type = "real"

    if self.params.scaling.pdb.include_bulk_solvent:
      params2.fmodel.k_sol = self.params.scaling.pdb.k_sol
      params2.fmodel.b_sol = self.params.scaling.pdb.b_sol

    # vvv These params restore the "legacy" solvent mask generation before
    # vvv cctbx commit 2243cc9a
    if self.params.scaling.pdb.solvent_algorithm == "flat":
      params2.mask.Fmask_res_high = 0
      params2.mask.grid_step_factor = 4
      params2.mask.solvent_radius = 1.11
      params2.mask.use_resolution_based_gridding = True
    # ^^^

    # Build an array of the model intensities according to the input parameters
    f_model = mmtbx.utils.fmodel_from_xray_structure(xray_structure = xray_structure,
                                                     f_obs          = None,
                                                     add_sigmas     = True,
                                                     params         = params2).f_model
    if not self.params.merging.merge_anomalous:
      f_model = f_model.generate_bijvoet_mates()

    return f_model.as_intensity_array().change_basis(self.params.scaling.model_reindex_op).map_to_asu()

  def create_model_from_mtz(self, model_file_path):

    if self.mpi_helper.rank == 0:
      from iotbx import mtz
      data_SR = mtz.object(model_file_path)
      arrays = data_SR.as_miller_arrays()
      space_group = data_SR.space_group().info()
      unit_cell   = data_SR.crystals()[0].unit_cell()
    else:
      arrays = space_group = unit_cell = None

    if self.purpose != "cosym":
      arrays, space_group, unit_cell = self.mpi_helper.comm.bcast((arrays, space_group, unit_cell), root=0)

    # save space group and unit cell as scaling targets
    if self.purpose == "scaling":
      self.params.scaling.space_group = space_group
      self.params.scaling.unit_cell   = unit_cell

    if self.purpose in ["scaling", "cosym"]:
      mtz_column_F = str(self.params.scaling.mtz.mtz_column_F.lower())
    elif self.purpose == "statistics":
      mtz_column_F = str(self.params.statistics.cciso.mtz_column_F.lower())

    for array in arrays:
      this_label = array.info().label_string().lower()
      if True not in ["sig"+tag in this_label for tag in ["iobs","imean", mtz_column_F]]:
        continue

      return array.as_intensity_array().change_basis(self.params.scaling.model_reindex_op).map_to_asu()

    raise Exception("mtz did not contain expected label Iobs or Imean")

  def consistent_set_and_model(self, i_model=None):
    assert self.params.scaling.space_group, "Space group must be specified in the input parameters or a reference file must be present"
    # which unit cell are we using?
    if self.purpose == "scaling":
      assert self.params.scaling.unit_cell is not None, "Unit cell must be specified in the input parameters or a reference file must be present"
      unit_cell = self.params.scaling.unit_cell
      self.logger.log("Using target unit cell: " + str(unit_cell))
      if self.mpi_helper.rank == 0:
        self.logger.main_log("Using target unit cell: " + str(unit_cell))
    elif self.purpose == "statistics":
      if self.params.merging.set_average_unit_cell:
        assert self.params.statistics.average_unit_cell is not None, "Average unit cell hasn't been calculated"
        unit_cell = self.params.statistics.average_unit_cell
        unit_cell_formatted = "(%.6f, %.6f, %.6f, %.3f, %.3f, %.3f)"\
                          %(unit_cell.parameters()[0], unit_cell.parameters()[1], unit_cell.parameters()[2], \
                            unit_cell.parameters()[3], unit_cell.parameters()[4], unit_cell.parameters()[5])
        self.logger.log("Using average unit cell: " + unit_cell_formatted)
        if self.mpi_helper.rank == 0:
          self.logger.main_log("Using average unit cell: " + unit_cell_formatted)
      else:
        assert self.params.scaling.unit_cell is not None, "Unit cell must be specified in the input parameters or a reference file must be present"
        unit_cell = self.params.scaling.unit_cell
        self.logger.log("Using target unit cell: " + str(unit_cell))
        if self.mpi_helper.rank == 0:
          self.logger.main_log("Using target unit cell: " + str(unit_cell))

    # create symmetry for the full miller set
    symm = symmetry(unit_cell=unit_cell, space_group_info = self.params.scaling.space_group)

    # Adjust the minimum d-spacing of the generated Miller set to assure
    # that the desired high-resolution limit is included even if the
    # observed unit cell differs slightly from the target.  Use the same
    # expansion formula as used in merging/general_fcalc.py, to assure consistency.
    # If a reference model is present, ensure that Miller indices are ordered
    # identically.

    # set up the resolution limits
    d_max = 100000 # a default like in cxi-merge
    if self.params.merging.d_max != None:
      d_max = self.params.merging.d_max
    # RB: for later
    #d_max /= self.params.scaling.resolution_scalar
    d_min = self.params.merging.d_min * self.params.scaling.resolution_scalar

    # build the full miller set
    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=d_max, d_min=d_min)
    miller_set = miller_set.change_basis(self.params.scaling.model_reindex_op).map_to_asu()

    # Handle the case where model is anomalous=False but the requested merging is anomalous=True
    if i_model is not None:
      if i_model.anomalous_flag() is False and miller_set.anomalous_flag() is True:
        i_model = i_model.generate_bijvoet_mates()

      # manage the sizes of arrays. General_fcalc assures that
      # N(i_model) >= N(miller_set) since it fills non-matches with invalid structure factors
      # However, if N(i_model) > N(miller_set), it's because this run of cxi.merge requested
      # a smaller resolution range.  Must prune off the reference model.

      #RB 10/07/2019 The old cxi.merge comment regarding N(i_model) vs. N(miller_set) - see above - refers to pdb references only.
      #Now applying the same approach to all cases when the reference is mtz.

      if self.purpose == "scaling":
        is_mtz = str(self.params.scaling.model).endswith(".mtz")
        if i_model.indices().size() > miller_set.indices().size() or is_mtz:
          matches = miller.match_indices(i_model.indices(), miller_set.indices())
          pairs = matches.pairs()
          i_model = i_model.select(pairs.column(0))
        matches = miller.match_indices(i_model.indices(), miller_set.indices())
        if is_mtz:
          assert matches.pairs().size() >= self.params.scaling.mtz.minimum_common_hkls, "Number of common HKLs between mtz reference and data (%d) is less than required (%d)."\
                 %(matches.pairs().size(), self.params.scaling.mtz.minimum_common_hkls)
        miller_set = miller_set.select(matches.permutation())

    return miller_set, i_model

  def add_asu_miller_indices_column(self, reflections):
    '''Add a "symmetry-reduced hkl" column to the reflection table'''
    if reflections.size() == 0:
      return

    import copy

    # Build target symmetry. The exact experiment unit cell values don't matter for converting HKLs to asu HKLs.
    target_unit_cell = self.params.scaling.unit_cell
    target_space_group_info = self.params.scaling.space_group
    target_symmetry = symmetry(unit_cell=target_unit_cell, space_group_info=target_space_group_info)
    target_space_group = target_symmetry.space_group()

    # generate and add an asu hkl column
    if 'miller_index_asymmetric' not in reflections:
      reflections['miller_index_asymmetric'] = copy.deepcopy(reflections['miller_index'])
      miller.map_to_asu(target_space_group.type(),
                        not self.params.merging.merge_anomalous,
                        reflections['miller_index_asymmetric'])
  '''
  def get_min_max_experiment_unit_cell_volume(self, experiments):

    vols = []
    for experiment in experiments:
      vols.append(experiment.crystal.get_crystal_symmetry().unit_cell().volume())

    return min(vols),max(vols)
  '''

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(crystal_model)
