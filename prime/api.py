from __future__ import division

def refine_many(frame_lst, input_file):
  """ Perform a full run on Prime, using the parameters specified in the `input_file`. Merging of fully-corrected intensities is not performed.

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
  :return: A list of `miller_array` objects, containing the partiality-corrected (full-intensity equivalent) intensities, and associated errors for each of the images specified in the input dictionairy.

  .. warning::
     This will modify the `orientation` object to reflect the updated parameters. Deep copy first if you do not want this!

  """
  pass
