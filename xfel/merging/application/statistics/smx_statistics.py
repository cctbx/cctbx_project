from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
import numpy as np
from dxtbx import flumpy
from dxtbx.model.experiment_list import ExperimentList
import os
from xfel.small_cell.command_line.powder_from_spots import phil_scope as pfs_phil_scope
from xfel.small_cell.powder_util import Spotfinder_radial_average
from iotbx.phil import parse

class smx_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(smx_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Prepare info for unit cell determination'

  def run(self, experiments, reflections):
    self.logger.log_step_time("SMX_STATISTICS")

    if self.params.statistics.smx.save_combined or self.params.statistics.smx.save_powder_from_spots:
      all_experiments_gathered = self.mpi_helper.comm.gather(experiments, 0)
      all_reflections_gathered = self.mpi_helper.comm.gather(reflections, 0)

      if self.mpi_helper.rank == 0:
        all_experiments = ExperimentList()
        for expts in all_experiments_gathered:
          if expts:
            all_experiments.extend(expts)
        all_reflections = flex.reflection_table.concat(all_reflections_gathered)

        all_experiments.as_file(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_combined.expt"))
        all_reflections.as_file(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_combined.refl"))

    if self.params.statistics.smx.save_powder_from_spots:
      if self.mpi_helper.rank == 0:
        overrides = "output.plot_file=%s\nplot.interactive=False"%(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_powder.png"))
        pfs_params = pfs_phil_scope.fetch(parse(overrides)).extract()
        averager = Spotfinder_radial_average(all_experiments, all_reflections, pfs_params)
        averager.calculate()
        averager.plot()

    if self.params.statistics.smx.save_triplets:
      triplet_data = TripletData(experiments, reflections, self.params.statistics.smx)
      all_triplets = self.mpi_helper.comm.gather(triplet_data.triplets, 0)
      if self.mpi_helper.rank == 0:
        triplets = np.vstack(all_triplets)
        fname = self.params.output.prefix + "_triplets.npz"
        np.savez(fname, triplets=triplets)

    return experiments, reflections


class TripletData:
  """Generate and store triplets of spots with their geometric relationships
  """
  def __init__(self, experiments, reflections, params):
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    self.triplets = None

    self.compute_triplets()

  def compute_triplets(self):
    """Compute all d1,d2,theta triplets from input data

    Fills triplets array with columns:
      (frame_id, d1, d2, theta, spot1_idx, spot2_idx, hand)
    """
    if self.triplets is not None:
      return


    triplets = []

    if 's1' not in self.reflections.keys():
      self.reflections.centroid_px_to_mm(self.experiments)
      self.reflections.map_centroids_to_reciprocal_space(self.experiments)
    for i_expt, expt in enumerate(self.experiments):

      # Get spots for this frame
      sel = flex.bool(self.reflections['id'] == i_expt)
      frame_refls = self.reflections.select(sel)

      # Get detector coordinates relative to beam center
      xyz = frame_refls['xyzobs.mm.value']
      beam_x, beam_y = expt.detector[0].get_beam_centre(expt.beam.get_s0())
      x = xyz.parts()[0] - beam_x
      y = xyz.parts()[1] - beam_y

      # Compute d-spacings and filter
      s0 = expt.beam.get_s0()
      s1_vectors = frame_refls['s1'].as_numpy_array()
      d_spacings = np.array([1/np.linalg.norm(s1 - s0) for s1 in s1_vectors])

      valid_spots = d_spacings >= self.params.triplets.d_min
      if not np.any(valid_spots):
        continue

      d_spacings = d_spacings[valid_spots]
      s1_vectors = s1_vectors[valid_spots]
      spot_indices = np.arange(len(frame_refls))[valid_spots]
      flex_valid_spots = flumpy.from_numpy(valid_spots)
      x = x.select(flex_valid_spots)
      y = y.select(flex_valid_spots)

      # Compute all pairs
      n_spots = len(d_spacings)
      for i in range(n_spots):
        for j in range(i+1, n_spots):
          # Compute angle
          vec_a = s1_vectors[i] - s0
          vec_b = s1_vectors[j] - s0
          cos_angle = np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b))
          angle = np.rad2deg(np.arccos(min(1.0, max(-1.0, cos_angle))))

          # Make sure d1 > d2 and track coordinates
          if d_spacings[i] > d_spacings[j]: # correct order
            d1 = d_spacings[i]
            d2 = d_spacings[j]
            i1 = spot_indices[i]
            i2 = spot_indices[j]
            x1, y1 = x[i], y[i]
            x2, y2 = x[j], y[j]
          else: # Switch them
            d1 = d_spacings[j]
            d2 = d_spacings[i]
            i1 = spot_indices[j]
            i2 = spot_indices[i]
            x1, y1 = x[j], y[j]
            x2, y2 = x[i], y[i]

          # Determine handedness from cross product z-component
          # For d1 > d2:
          # cross_z > 0 means counterclockwise rotation from 1 to 2 (left-handed)
          cross_z = x1*y2 - y1*x2
          hand = -1 if cross_z > 0 else 1

          triplets.append((i_expt, d1, d2, angle, i1, i2, hand))

    self.triplets = np.array(triplets)


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(smx_statistics)
