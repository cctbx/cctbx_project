from __future__ import print_function, division
from dials.array_family import flex
from xfel.merging.application.worker import worker
from six.moves import range

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

  def get_next_hkl_reflection_table(self, reflections):
    '''Generate asu hkl slices from an asu hkl-sorted reflection table'''
    i_begin = 0
    hkl_ref = reflections[0].get('miller_index_asymmetric')

    for i in range(len(reflections)):
      hkl = reflections[i].get('miller_index_asymmetric')
      if hkl == hkl_ref:
        continue
      else:
        yield reflections[i_begin:i]
        i_begin = i
        hkl_ref = hkl

    yield reflections[i_begin:i+1]

  def output_merged_reflections(self, reflections):
    merged_reflections_file_path = self.params.output.output_dir + '/merge.out'
    merged_file = open(merged_reflections_file_path, 'w')

    for ref in reflections:
      merged_file.write("%s %f %f %f %d\n"%(ref.get('miller_index'), ref.get('intensity'), ref.get('esd'), ref.get('rmsd'), ref.get('multiplicity')))

    merged_file.close()

  def run(self, experiments, reflections):

    # merge reflection intensities: calculate the average and other statistics
    self.logger.log_step_time("AVERAGE")
    self.logger.log("Averaging intensities...")
    all_rank_merged_reflections = self.merging_reflection_table()

    if len(reflections) > 0:
      for hkl_reflection_table in self.get_next_hkl_reflection_table(reflections):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=hkl_reflection_table)
        intensity_stats['miller_index'] = hkl_reflection_table[0].get('miller_index_asymmetric')
        all_rank_merged_reflections.append(intensity_stats)

    self.logger.log("Merged intensities for %d HKLs"%(all_rank_merged_reflections.size()))
    self.logger.log_step_time("AVERAGE", True)

    # gather all merged intensities at rank 0
    self.logger.log_step_time("GATHER")
    if self.mpi_helper.rank != 0:
      self.logger.log("Executing MPI gathering of all reflection tables at rank 0...")
    all_merged_reflection_tables = self.mpi_helper.comm.gather(all_rank_merged_reflections, root = 0)
    self.logger.log_step_time("GATHER", True)

    # rank 0: concatenate all merged intensities into the final table
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("MERGE")
      final_merged_reflection_table = self.merging_reflection_table()

      self.logger.log("Performing final merging of reflection tables received from all ranks...")
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      self.logger.main_log("Total merged HKLs: {}".format(final_merged_reflection_table.size()))
      self.logger.log_step_time("MERGE", True)

      # write the final merged reflection table out to an ASCII file
      self.logger.log_step_time("WRITE")
      self.output_merged_reflections(final_merged_reflection_table)
      self.logger.log_step_time("WRITE", True)

    return None, None
