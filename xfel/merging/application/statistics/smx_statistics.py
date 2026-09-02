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
from dxtbx.model.experiment_list import DetectorComparison

class smx_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(smx_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Prepare info for unit cell determination'

  def run(self, experiments, reflections):
    self.logger.log_step_time("SMX_STATISTICS")

    if self.params.statistics.smx.group_by_identifier_prefix:
      prefix_length = self.params.statistics.smx.identifier_prefix_length
      prefixes = list(set([e.identifier[:prefix_length] for e in experiments]))
      all_prefixes_gathered = self.mpi_helper.comm.gather(prefixes, 0)
      if self.mpi_helper.rank == 0:
        all_prefixes = set(sum(all_prefixes_gathered, []))
      else:
        all_prefixes = None
      all_prefixes = self.mpi_helper.comm.bcast(all_prefixes, 0)

      root_prefix = self.params.output.prefix
      for prefix in sorted(all_prefixes):
        group_expts = ExperimentList()
        sel = flex.bool(len(reflections), False)
        for expt_id, expt in enumerate(experiments):
          if expt.identifier[:prefix_length] == prefix:
            group_expts.append(expt)
            sel |= reflections['id'] == expt_id
        self.params.output.prefix = root_prefix + "_" + prefix
        self.logger.log('Saving group %s, %d experiments'%(self.params.output.prefix, len(group_expts)))
        self.save_group(group_expts, reflections.select(sel))
      self.params.output.prefix = root_prefix
    else:
      self.save_group(experiments, reflections)

    return experiments, reflections

  def save_group(self, experiments, reflections):
    if self.params.statistics.smx.save_combined or self.params.statistics.smx.save_powder_from_spots:
      if self.mpi_helper.rank == 0:
        # Rank 0 collects from all other ranks
        all_experiments_gathered = [experiments]  # Start with rank 0's own data
        all_reflections_gathered = [reflections]
        
        for rank in range(1, self.mpi_helper.size):
          recv_experiments = self.mpi_helper.comm.recv(source=rank, tag=0)
          recv_reflections = self.mpi_helper.comm.recv(source=rank, tag=1)
          all_experiments_gathered.append(recv_experiments)
          all_reflections_gathered.append(recv_reflections)
      else:
        # All other ranks send their data to rank 0
        self.mpi_helper.comm.send(experiments, dest=0, tag=0)
        self.mpi_helper.comm.send(reflections, dest=0, tag=1)
        all_experiments_gathered = None
        all_reflections_gathered = None
      #all_experiments_gathered = self.mpi_helper.comm.gather(experiments, 0)
      #all_reflections_gathered = self.mpi_helper.comm.gather(reflections, 0)
      #all_reflections = self.mpi_helper.gather_reflection_table(reflections)

      if self.mpi_helper.rank == 0:
        all_experiments = ExperimentList()
        for expts in all_experiments_gathered:
          if expts:
            all_experiments.extend(expts)
        all_reflections = flex.reflection_table.concat(all_reflections_gathered)
        if all_experiments:
          compare_detector = DetectorComparison()
          detector = all_experiments[0].detector
          for expt in all_experiments:
            if expt.detector is detector: continue
            assert compare_detector(detector, expt.detector)
            expt.detector = detector
        all_experiments.as_file(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_combined.expt"))
        all_reflections.as_file(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_combined.refl"))

    if self.params.statistics.smx.save_powder_from_spots:
      if self.mpi_helper.rank == 0:
        overrides = "output.plot_file=%s\nplot.interactive=False\nd_min=2"%(os.path.join(self.params.output.output_dir, self.params.output.prefix + "_powder.png"))
        pfs_params = pfs_phil_scope.fetch(parse(overrides)).extract()
        averager = Spotfinder_radial_average(all_experiments, all_reflections, pfs_params)
        averager.calculate()
        averager.plot()

    if self.params.statistics.smx.save_triplets:
      #import line_profiler, io
      #lp = line_profiler.LineProfiler(TripletData.compute_triplets)
      #lp.enable()
      triplet_data = TripletData(experiments, reflections, self.params.statistics.smx)
      #lp.disable()
      #s = io.StringIO()
      #lp.print_stats(stream=s)
      #self.logger.log(s.getvalue())
      self.logger.log("Have %d triplets from %d experiments and %d reflections"%(triplet_data.triplets.size, len(experiments), len(reflections)))
      all_triplets = self.mpi_helper.comm.gather(triplet_data.triplets, 0)
      if self.mpi_helper.rank == 0:
        if any([t.size for t in all_triplets]):
          triplets = np.vstack([t for t in all_triplets if t.size])
          self.logger.main_log("Have %d total triplets"%(triplets.size))
          fname = os.path.join(self.params.output.output_dir, self.params.output.prefix + "_triplets.npz")
          np.savez(fname, triplets=triplets)

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

    if self.params.triplets.mask:
      from libtbx import easy_pickle
      import math
      mask = easy_pickle.load(self.params.triplets.mask)
    else:
      mask = None

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

      if self.params.triplets.mask:
        mask_shape = mask[0].as_numpy_array().shape
        for i_refl in range(len(frame_refls)):
          if valid_spots[i_refl]:
            refl = frame_refls[i_refl]
            panel = refl['panel']
            center_y = int(round(refl['xyzobs.px.value'][1]))
            center_x = int(round(refl['xyzobs.px.value'][0]))
            
            for dy in range(-2, 3):  # -2, -1, 0, 1, 2
              for dx in range(-2, 3):
                pixel_y = center_y + dy
                pixel_x = center_x + dx
                if 0 <= pixel_y < mask_shape[0] and 0 <= pixel_x < mask_shape[1]:
                  if mask[panel][pixel_y, pixel_x] == False:
                    valid_spots[i_refl] = False
                    break
              if not valid_spots[i_refl]:
                break
 
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
