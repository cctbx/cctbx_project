from __future__ import absolute_import, division, print_function
from abc import ABCMeta, abstractmethod
import cctbx
import scitbx_array_family_flex_ext
import cctbx_orientation_ext

class InputFrame():
  """ Wrapper class for describing single images to be refined, created from an input dictionary of schema:
        - `miller_array`: the cctbx.miller miller array of spot intensities.
        - `mapped_predictions`: the mapped_predictions locations
        - `orientation`: cctbx crystal_orientation object
        - `xbeam`: x-location of beam centre
        - `ybeam`: y-location of beam centre
        - `wavelength`: the wavelength in Angstroms

  Child classes must simply have these attributes, and will be modified in place. See xfel.cluster.singleframe.SingleFrame for example.
  """
  def __init__(self, frame_dict):
    self.__dict__.update(frame_dict)

  def check_prime_input(self):
    """ Check attribute names and types to ensure that the object is valid. """
    required = {'miller_array': cctbx.miller.array,
                'mapped_predictions': scitbx_array_family_flex_ext.vec2_double,
                'orientation': cctbx_orientation_ext.crystal_orientationa,
                'xbeam': float,
                'ybeam': float,
                'wavelength': float}
    for key in required:
      assert isinstance(getattr(self, key), required[key])


def refine_many(frame_lst, input_file):
  """ Perform a full run on Prime, using the parameters specified in the `input_file`. Merging of fully-corrected intensities is not performed.

  :param frame_lst: list of InputFrame objects describing integration results.
  :param input_file: Prime .inp file. See XXXXXX for full description.
  :return: A list of `miller_array` objects, containing the partiality-corrected (full-intensity equivalent) intensities, and associated errors for each of the images specified in the frame_lst.

  .. note::
     This will also update the InputFrame objects in place.
  """
  for frame in frame_lst:
    frame.check_prime_input()

  pass

class CustomMicrocycle():
  """ Abstract class for using the Prime microcycle algorithm to refine a single frame. Implement `func` to specify a target function."""

  __metaclass__ = ABCMeta

  def __init__(self, frame, input_file):
    """
    :param frame: an InputFrame-like object.
    :param input_file: Prime .inp file. Only uses micro-cycle parameters.

    .. note::
       This will also update the InputFrame objects in place.
    """
    frame.check_prime_input()
    pass  # Maybe create internal prime input class here?

  def __call__(self):
    """ perform post-refinement on the single frame.
    :return: a `miller_array` object containing the partiality-corrected (full-intensity equivalent) intensities, and associated errors.
    """
    pass

  @abstractmethod
  def func(self, corr_miller):
    """ Penalty function for refinement, an abstract class that must be implemented.

    :param corr_miller: miller array of correction factors (such that I_full = corr_miller * I_obs).
    :return: penalty score that should be minimized.

    """
    pass
