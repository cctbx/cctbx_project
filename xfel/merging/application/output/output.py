from __future__ import print_function, division
from dials.array_family import flex
from xfel.merging.application.worker import worker
from cctbx import miller, uctbx
from cctbx.crystal import symmetry
from xfel.command_line.cxi_merge import unit_cell_distribution
import os

class output(worker):
  """
  Outputs merged mupltiple measurements of symmetry-reduced hkl's.
  """

  def __repr__(self):
    return "Output merged reflections"

  def output_reflections_ascii(self, reflections):
    merged_reflections_file_path = os.path.join(self.params.output.output_dir + 'merge.out')
    merged_file = open(merged_reflections_file_path, 'w')
    for ref in reflections:
      merged_file.write("%s %f %f %f %d\n"%(ref.get('miller_index'), ref.get('intensity'), ref.get('esd'), ref.get('rmsd'), ref.get('multiplicity')))
    merged_file.close()

  def output_reflections_mtz(self, experiments, reflections):
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
    if self.mpi_helper.rank == 0:
      # write the final merged reflection table out to an ASCII file -- for testing/debugging
      self.logger.log_step_time("WRITE ASCII")
      self.output_reflections_ascii(reflections)
      self.logger.log_step_time("WRITE ASCII", True)
      # write the final merged reflection table out to an MTZ file
      self.logger.log_step_time("WRITE MTZ")
      self.output_reflections_mtz(experiments, reflections)
      self.logger.log_step_time("WRITE MTZ", True)

    return None, None

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(output)
