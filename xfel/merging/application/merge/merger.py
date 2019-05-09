from __future__ import print_function, division
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from dxtbx.model.experiment_list import ExperimentList
from cctbx import miller, uctbx
from cctbx.crystal import symmetry
from xfel.command_line.cxi_merge import unit_cell_distribution
import os

try:
  import resource
  import platform
  def get_memory_usage():
    # getrusage returns kb on linux, bytes on mac
    units_per_mb = 1024
    if platform.system() == "Darwin":
      units_per_mb = 1024*1024
    return ('Memory usage: %.1f MB' % (int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb))
except ImportError:
  def debug_memory_usage():
    pass

class merger(worker):
  """
  Merges multiple measurements of symmetry-reduced hkl's.
  """

  def __repr__(self):
    return "Merge multiple measurements of symmetry-reduced hkl's"

  def merging_reflection_table(self):
    '''Create a reflection table for storing merged hkl's'''
    table = flex.reflection_table()
    table['miller_index'] = flex.miller_index()
    table['intensity']    = flex.double()
    table['esd']          = flex.double()
    table['rmsd']         = flex.double()
    table['multiplicity'] = flex.int()
    return table

  def calc_reflection_intensity_stats(self, reflections):
    '''Calculate intensity statistics for reflection table'''
    multiplicity = len(reflections)
    assert multiplicity != 0

    stats = flex.mean_and_variance(reflections['intensity.sum.value'])
    propagated_esd = (flex.sum(reflections['intensity.sum.variance']) ** 0.5)/ multiplicity

    rmsd = 0.0
    if multiplicity > 1:
      rmsd = stats.unweighted_sample_standard_deviation()

    return {'intensity'     : stats.mean(),
            'esd'           : propagated_esd,
            'rmsd'          : rmsd,
            'multiplicity'  : multiplicity}

  def output_merged_reflections(self, experiments, reflections):
    merged_reflections_file_path = os.path.join(self.params.output.output_dir + 'merge.out')
    merged_file = open(merged_reflections_file_path, 'w')

    for ref in reflections:
      merged_file.write("%s %f %f %f %d\n"%(ref.get('miller_index'), ref.get('intensity'), ref.get('esd'), ref.get('rmsd'), ref.get('multiplicity')))
    merged_file.close()

    if self.params.merging.set_average_unit_cell:
      dist = unit_cell_distribution()
      for crystal in experiments.crystals():
        dist.add_cell(crystal.get_unit_cell())
      abc = dist.get_average_cell_dimensions()
      angles = crystal.get_unit_cell().parameters()[3:]
      unit_cell = uctbx.unit_cell(list(abc) + list(angles))
    else:
      unit_cell = self.params.scaling.unit_cell
    final_symm = symmetry(
      unit_cell=unit_cell,
      space_group_info=self.params.scaling.space_group)

    wavelength = flex.mean(flex.double([b.get_wavelength() for b in experiments.beams()]))

    all_obs = miller.array(
      miller_set=miller.set(final_symm, reflections['miller_index']),
      data=reflections['intensity'],
      sigmas=reflections['esd']).resolution_filter(
        d_min=self.params.merging.d_min,
        d_max=self.params.merging.d_max).set_observation_type_xray_intensity()
    mtz_file = os.path.join(self.params.output.output_dir, "%s.mtz" % self.params.output.prefix)

    mtz_out = all_obs.as_mtz_dataset(
      column_root_label="Iobs",
      title=self.params.output.title,
      wavelength=wavelength)
    mtz_out.add_miller_array(
      miller_array=all_obs.average_bijvoet_mates(),
      column_root_label="IMEAN")
    mtz_obj = mtz_out.mtz_object()
    mtz_obj.write(mtz_file)
    self.logger.log("  Anomalous and mean data:\n    %s" % \
      os.path.abspath(mtz_file))
    self.logger.log("")
    self.logger.log("Final data:")
    #all_obs.show_summary(self.log, prefix="  ") # don't have a buffer object for this logger
    all_obs.show_summary(prefix="  ")

  def run(self, experiments, reflections):

    # merge reflection intensities: calculate the average and other statistics
    self.logger.log_step_time("AVERAGE")
    self.logger.log("Averaging intensities...")
    all_rank_merged_reflections = self.merging_reflection_table()

    if len(reflections) > 0:
      for hkl_reflection_table in reflection_table_utils.get_next_hkl_reflection_table(reflections):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=hkl_reflection_table)
        intensity_stats['miller_index'] = hkl_reflection_table[0].get('miller_index_asymmetric')
        all_rank_merged_reflections.append(intensity_stats)

    self.logger.log("Merged intensities for %d HKLs"%(all_rank_merged_reflections.size()))
    self.logger.log_step_time("AVERAGE", True)

    # gather all merged intensities at rank 0
    self.logger.log_step_time("GATHER")
    if self.mpi_helper.rank != 0:
      self.logger.log("Executing MPI gathering of all reflection tables at rank 0...")
    results = self.mpi_helper.comm.gather((experiments, all_rank_merged_reflections), root = 0)
    self.logger.log_step_time("GATHER", True)

    # rank 0: concatenate all merged intensities into the final table
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("MERGE")
      all_experiments = ExperimentList()
      final_merged_reflection_table = self.merging_reflection_table()

      self.logger.log("Performing final merging of reflection tables received from all ranks...")
      for expts, table in results:
        all_experiments.extend(expts)
        final_merged_reflection_table.extend(table)
      self.logger.main_log("Total merged HKLs: {}".format(final_merged_reflection_table.size()))
      self.logger.log_step_time("MERGE", True)

      # write the final merged reflection table out to an ASCII file
      self.logger.log_step_time("WRITE")
      self.output_merged_reflections(all_experiments, final_merged_reflection_table)
      self.logger.log_step_time("WRITE", True)

    return None, None
