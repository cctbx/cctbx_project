from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker

class resolution_binner(worker):

  def __repr__(self):
    return 'Set up resolution bins'

  def run(self, experiments, reflections):
    '''Set up resolution bins; assign bin number to each hkl in the full miller set'''

    # Create resolution binner
    full_miller_set = self.params.scaling.miller_set
    full_miller_set.setup_binner(d_max=100000, d_min=self.params.merging.d_min, n_bins=self.params.statistics.n_bins)
    resolution_binner = full_miller_set.binner()

    # Save resolution binner to the parameters
    self.params.statistics.__inject__('resolution_binner', resolution_binner)

    # Provide resolution bin number for each asu hkl
    hkl_resolution_bins = {} # hkl vs resolution bin number
    hkls_with_assigned_bin = 0
    for i_bin in resolution_binner.range_used():
      bin_hkl_selection = resolution_binner.selection(i_bin)
      bin_hkls = full_miller_set.select(bin_hkl_selection)
      for hkl in bin_hkls.indices():
        assert not hkl in hkl_resolution_bins # each hkl should be assigned a bin number only once
        hkl_resolution_bins[hkl] = i_bin
        hkls_with_assigned_bin += 1
    self.logger.log("Provided resolution bin number for %d asu hkls"%(hkls_with_assigned_bin))

    # Save hkl bin asignments to the parameters
    self.params.statistics.__inject__('hkl_resolution_bins', hkl_resolution_bins)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(resolution_binner)
