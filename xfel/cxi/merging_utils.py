from __future__ import division
from __future__ import print_function
from scitbx.array_family import flex
from libtbx import adopt_init_args

class intensity_data (object) :
  """
  Container for scaled intensity data.
  """
  def __init__ (self, n_refl) :
    self.n_refl = n_refl
    self.initialize()

  def initialize (self) :
    self.ISIGI        = {}
    self.completeness = flex.int(self.n_refl, 0)
    self.completeness_predictions = flex.int(self.n_refl, 0)
    self.summed_N     = flex.int(self.n_refl, 0)
    self.summed_weight= flex.double(self.n_refl, 0.)
    self.summed_wt_I  = flex.double(self.n_refl, 0.)

class frame_data (intensity_data) :
  """
  Intensity data for a single frame.
  """
  def __init__ (self, n_refl, file_name) :
    intensity_data.__init__(self, n_refl)
    self.file_name = file_name
    self.n_obs = 0
    self.n_rejected = 0
    self.corr = 0
    self.d_min = -1
    self.accept = False
    self.indexed_cell = None
    self.log_out = file_name
    self.wavelength = None

  def set_indexed_cell (self, unit_cell) :
    self.indexed_cell = unit_cell

  def set_log_out (self, out_str) :
    self.log_out = out_str

  def show_log_out (self, out) :
    print(self.log_out, file=out)

class null_data (object) :
  """
  Stand-in for a frame rejected due to conflicting symmetry.  (No flex arrays
  included, to save pickling time during multiprocessing.)
  """
  def __init__ (self, file_name, log_out,
                file_error=False,
                low_signal=False,
                wrong_bravais=False,
                wrong_cell=False,
                low_resolution=False,
                low_correlation=False,
                reason=None) :
    adopt_init_args(self, locals())

  def show_log_out (self, out) :
    print(self.log_out, file=out)
