from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from scitbx import matrix
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage

class polarization(worker):
  """
  Computes the polarization correction as defined by Kahn 1982.
  Modifies the intensity.sum.value and intensity.sum.variance columns in place.
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(polarization, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Apply polarization correction'

  def run(self, experiments, reflections):

    self.logger.log_step_time("POLARIZATION_CORRECTION")

    result = flex.reflection_table()

    for expt_id, experiment in enumerate(experiments):
      refls = reflections.select(reflections['id'] == expt_id)
      if len(refls) == 0: continue
      beam = experiment.beam
      # Remove the need for pixel size within cxi.merge.  Allows multipanel detector with dissimilar panels.
      # Relies on new frame extractor code called by dials.stills_process that writes s0, s1 and polarization normal
      # vectors all to the integration pickle.  Future path (IE THIS CODE): use dials json and reflection file.
      s0_vec = matrix.col(beam.get_s0()).normalize()
      s0_polar_norm = beam.get_polarization_normal()
      s1_vec = refls['s1']
      Ns1 = len(s1_vec)
      # project the s1_vector onto the plane normal to s0.  Get result by subtracting the
      # projection of s1 onto s0, which is (s1.dot.s0_norm)s0_norm
      s0_norm = flex.vec3_double(Ns1,s0_vec)
      s1_proj = (s1_vec.dot(s0_norm))*s0_norm
      s1_in_normal_plane = s1_vec - s1_proj
      # Now want the polar angle between the projected s1 and the polarization normal
      s0_polar_norms = flex.vec3_double(Ns1,s0_polar_norm)
      dotprod = (s1_in_normal_plane.dot(s0_polar_norms))
      costheta = dotprod/(s1_in_normal_plane.norms())
      theta = flex.acos(costheta)
      cos_two_polar_angle = flex.cos(2.0*theta)
      # gives same as old answer to ~1% but not exact.  Not sure why, should not matter.

      tt_vec = experiment.crystal.get_unit_cell().two_theta(miller_indices = refls['miller_index'],
                                                            wavelength = beam.get_wavelength())
      cos_tt_vec = flex.cos(tt_vec)
      sin_tt_vec = flex.sin(tt_vec)
      cos_sq_tt_vec = cos_tt_vec * cos_tt_vec
      sin_sq_tt_vec = sin_tt_vec * sin_tt_vec
      P_nought_vec = 0.5 * (1. + cos_sq_tt_vec)

      F_prime = -1.0 # Hard-coded value defines the incident polarization axis
      P_prime = 0.5 * F_prime * cos_two_polar_angle * sin_sq_tt_vec

      # added as a diagnostic
      #prange=P_nought_vec - P_prime
      #other_F_prime = 1.0
      #otherP_prime = 0.5 * other_F_prime * cos_two_polar_angle * sin_sq_tt_vec
      #otherprange=P_nought_vec - otherP_prime
      #diff2 = flex.abs(prange - otherprange)
      #print >> out, "mean diff is",flex.mean(diff2), "range",flex.min(diff2), flex.max(diff2)
      # done

      correction = 1 / ( P_nought_vec - P_prime )
      refls['intensity.sum.value'] = refls['intensity.sum.value'] * correction
      refls['intensity.sum.variance'] = refls['intensity.sum.variance'] * correction**2 # propagated error
      # This corrects observations for polarization assuming 100% polarization on
      # one axis (thus the F_prime = -1.0 rather than the perpendicular axis, 1.0)
      # Polarization model as described by Kahn, Fourme, Gadet, Janin, Dumas & Andre
      # (1982) J. Appl. Cryst. 15, 330-337, equations 13 - 15.

      result.extend(refls)

    if len(reflections) > 0:
      self.logger.log("Applied polarization correction. Mean intensity changed from %.2f to %.2f"%(flex.mean(reflections['intensity.sum.value']), flex.mean(result['intensity.sum.value'])))

    self.logger.log_step_time("POLARIZATION_CORRECTION", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    # Remove 's1' column from the reflection table
    from xfel.merging.application.reflection_table_utils import reflection_table_utils
    reflections = reflection_table_utils.prune_reflection_table_keys(reflections=result, keys_to_delete=['s1'],
                                                                     keys_to_ignore=self.params.input.persistent_refl_cols)
    self.logger.log("Pruned reflection table")
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(polarization)
