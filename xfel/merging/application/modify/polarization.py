from __future__ import print_function, division
from dials.array_family import flex
from scitbx import matrix

"""
Compute the polarazation correction as defined by Kahn 1982.

Modifies the intensity.sum.value and intensity.sum.variance columns
in place.
"""

class polarization(object):
  def __init__(self, params):
    self.params = params

  def __call__(self, experiments, reflections, add_correction_column=False):
    result = flex.reflection_table()

    for expt_id, experiment in enumerate(experiments):
      refls = reflections.select(reflections['id'] == expt_id)
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

      if add_correction_column:
        refls['polarization_correction'] = correction

      result.extend(refls)

    return result

from xfel.merging.application.phil.phil import Script as Script_Base
class Script(Script_Base):
  def modify(self, experiments, reflections):
    result = polarization(None)(experiments, reflections)
    print ("Mean intensity changed from %.2f to %.2f"%(flex.mean(reflections['intensity.sum.value']), flex.mean(result['intensity.sum.value'])))
    return experiments, reflections

if __name__ == '__main__':
  import sys
  from dxtbx.model.experiment_list import ExperimentListFactory
  from libtbx import easy_pickle
  experiments_filename, reflections_filename = sys.argv[1:3]
  experiments = ExperimentListFactory.from_json_file(experiments_filename, check_format=False)
  reflections = easy_pickle.load(reflections_filename)
  result = polarization(None)(experiments, reflections)
  print ("Mean intensity changed from %.2f to %.2f"%(flex.mean(reflections['intensity.sum.value']), flex.mean(result['intensity.sum.value'])))
