from __future__ import division

def refine_many(frame_lst, input_file):
  """ Run postrefienement.
  :param frame_lst: list of dictionairies describing integration results.
  :param input_file: Prime .inp file. See XXXXXX for full description.
  Schema for elements of `frame_lst`:
        - `miller_array`: the cctbx.miller miller array of spot intensities.
        - `mapped_predictions`: the mapped_predictions locations
        - `name`: file-name, used as an identifier
        - `pg`: point group of pickle
        - `orientation`: cctbx crystal_orientation object
        - `xbeam`: x-location of beam centre
        - `ybeam`: y-location of beam centre
        - `wavelength`: the wavelength in Angstroms
  :return: a list of tuples of the form: (G, B, rotx, roty, ry, rz, re, a, b, c, alpha, beta, gamma), describing the refined parameters of the images described by the list of dictionairies.
  """
  pass
