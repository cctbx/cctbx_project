from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
from cctbx import miller
from cctbx.crystal import symmetry

class crystal_model(worker):

  def __repr__(self):
    return 'Build crystal model'

  def run(self, experiments, reflections):
    '''Create a model for the crystal'''
    self.logger.log_step_time("CREATE_CRYSTAL_MODEL")
    # Generate the reference, or model, intensities from the input parameters
    # This will save the space group and unit cell into the input parameters
    if self.params.scaling.model.endswith(".mtz"):
      i_model = self.create_model_from_mtz()
    elif self.params.scaling.model.endswith(".pdb"):
      i_model = self.create_model_from_pdb()

    # Generate a full miller set consistent with the crystal model
    miller_set, i_model = self.consistent_set_and_model(i_model=i_model)

    # Save i_model and miller_set to the parameters
    self.params.scaling.__inject__('i_model', i_model)
    self.params.scaling.__inject__('miller_set', miller_set)

    if self.mpi_helper.rank == 0:
      self.logger.main_log("Scaling model: " + self.params.scaling.model)
      self.logger.main_log("Space group: " + str(self.params.scaling.space_group))
      self.logger.main_log("Unit cell: " + str(self.params.scaling.unit_cell))

    self.logger.log_step_time("CREATE_CRYSTAL_MODEL", True)

    return experiments, reflections

  def create_model_from_pdb(self):

    from iotbx import file_reader

    pdb_in = file_reader.any_file(self.params.scaling.model, force_type="pdb")
    pdb_in.assert_file_type("pdb")

    xray_structure = pdb_in.file_object.xray_structure_simple()

    self.params.scaling.space_group = xray_structure.crystal_symmetry().space_group().info()
    self.params.scaling.unit_cell = xray_structure.crystal_symmetry().unit_cell()

    if self.mpi_helper.rank == 0:
      xray_structure.show_summary()

    # prepare phil parameters to generate model intensities
    phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
    params2 = phil2.extract()

    '''
    #
    # RB: for later
    #
    # We need to set the model resolution limits here, which will determine the set of HKLs, for which the intensities will be calculated.
    # So our input resolution limits need to be adjusted to account for the difference between the model unit cell and the variations of the unit cells in all of the experiments.
    # Ideally, we should calculate the model intensity for _every_ observed HKL in every experiment, but that may be a time-wasteful calculation. So we use a shortcut:
    # a ratio of the corresponding unt cell volumes as well as a fudge factor - "resolution scalar".

    model_unit_cell_volume = xray_structure.crystal_symmetry().unit_cell().volume()

    min_exp_unit_cell_volume, max_exp_unit_cell_volume = self.get_min_max_experiment_unit_cell_volume(self.experiments)

    assert self.params.scaling.resolution_scalar < 1.0 # the way this scalar is used, it should be less than unity

    params2.high_resolution = self.params.merging.d_min / math.pow(max_exp_unit_cell_volume/model_unit_cell_volume, 1./3.) * self.params.scaling.resolution_scalar
    if self.params.merging.d_max is not None:
      params2.low_resolution = self.params.merging.d_max / math.pow(min_exp_unit_cell_volume/model_unit_cell_volume, 1./3.) / self.params.scaling.resolution_scalar
    '''

    #
    # RB: for now do it the way it's done in cxi-merge. There is a problem with that approach though, because the "filter unit cell" may be different from the "model unit cell",
    # so the cell variations we are trying to account for here might be for a different unit cell.
    #
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

    # Build an array of the model intensities according to the input parameters
    f_model = mmtbx.utils.fmodel_from_xray_structure(xray_structure = xray_structure,
                                                     f_obs          = None,
                                                     add_sigmas     = True,
                                                     params         = params2).f_model
    if not self.params.merging.merge_anomalous:
      f_model = f_model.generate_bijvoet_mates()

    return f_model.as_intensity_array().change_basis(self.params.scaling.model_reindex_op).map_to_asu()

  def create_model_from_mtz(self):

    from iotbx import mtz
    data_SR = mtz.object(self.params.scaling.model)

    self.params.scaling.space_group = data_SR.space_group().info()
    self.params.scaling.unit_cell   = data_SR.crystals()[0].unit_cell()

    for array in data_SR.as_miller_arrays():
      this_label = array.info().label_string().lower()
      if True not in [this_label.find(tag)>=0 for tag in ["iobs","imean", self.params.scaling.mtz.mtz_column_F]]:
        continue

      return array.as_intensity_array().change_basis(self.params.scaling.model_reindex_op).map_to_asu()

    raise Exception("mtz did not contain expected label Iobs or Imean")

  def consistent_set_and_model(self, i_model=None):
    # Adjust the minimum d-spacing of the generated Miller set to assure
    # that the desired high-resolution limit is included even if the
    # observed unit cell differs slightly from the target.  Use the same
    # expansion formula as used in merging/general_fcalc.py, to assure consistency.
    # If a reference model is present, ensure that Miller indices are ordered
    # identically.

    symm = symmetry(unit_cell = self.params.scaling.unit_cell, space_group_info = self.params.scaling.space_group)

    # set up the resolution limits
    d_max = 100000 # a default like in cxi-merge
    if self.params.merging.d_max != None:
      d_max = self.params.merging.d_max
    # RB: for later
    #d_max /= self.params.scaling.resolution_scalar
    d_min = self.params.merging.d_min * self.params.scaling.resolution_scalar

    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=d_max, d_min=d_min)
    miller_set = miller_set.change_basis(self.params.scaling.model_reindex_op).map_to_asu()

    # Handle the case where model is anomalous=False but the requested merging is anomalous=True
    if i_model.anomalous_flag() is False and miller_set.anomalous_flag() is True:
      i_model = i_model.generate_bijvoet_mates()

    # manage the sizes of arrays. General_fcalc assures that
    # N(i_model) >= N(miller_set) since it fills non-matches with invalid structure factors
    # However, if N(i_model) > N(miller_set), it's because this run of cxi.merge requested
    # a smaller resolution range.  Must prune off the reference model.
    if i_model.indices().size() > miller_set.indices().size():
      matches = miller.match_indices(i_model.indices(), miller_set.indices())
      pairs = matches.pairs()
      i_model = i_model.select(pairs.column(0))

    matches = miller.match_indices(i_model.indices(), miller_set.indices())
    assert not matches.have_singles()
    miller_set = miller_set.select(matches.permutation())

    return miller_set, i_model

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
