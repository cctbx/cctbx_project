from __future__ import absolute_import, division, print_function
from iotbx import reflection_file_reader

class mx_handler(object):
  """
  Author      : Uervirojnangkoorn, M.
  Created     : 4/15/2016
  A collection of macromolecular-crystallagphic wrapper functions
  """
  def __init__(self):
    """
    Constructor
    """

  def get_asu_contents(self, n_residues):
    asu_contents = None
    if n_residues > 0:
      asu_contents = {"H":8.0*float(n_residues),
                      "C":5.0*float(n_residues),
                      "N":1.5*float(n_residues),
                      "O":1.2*float(n_residues)
                     }
    return asu_contents

  def get_miller_array_from_reflection_file(self, hklisoin):
    #get iso if given
    flag_hklisoin_found = False
    miller_array_iso = None
    if hklisoin is not None:
      flag_hklisoin_found = True
      reflection_file_iso = reflection_file_reader.any_reflection_file(hklisoin)
      miller_arrays_iso=reflection_file_iso.as_miller_arrays()
      is_found_iso_as_intensity_array = False
      is_found_iso_as_amplitude_array = False
      for miller_array in miller_arrays_iso:
        if miller_array.is_xray_intensity_array():
          miller_array_iso = miller_array.deep_copy()
          is_found_iso_as_intensity_array = True
          break
        elif miller_array.is_xray_amplitude_array():
          is_found_iso_as_amplitude_array = True
          miller_array_converted_to_intensity = miller_array.as_intensity_array()
      if is_found_iso_as_intensity_array == False:
        if is_found_iso_as_amplitude_array:
          miller_array_iso = miller_array_converted_to_intensity.deep_copy()
        else:
          flag_hklisoin_found = False
    return flag_hklisoin_found, miller_array_iso

  def get_resolution_step_for_B(self, iparams):
    resolution_gap = 7 - iparams.scale.d_min
    resolution_step = resolution_gap/ iparams.n_bins
    return resolution_step
