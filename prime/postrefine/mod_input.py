from __future__ import division
import iotbx.phil
from libtbx.utils import Usage, Sorry
import sys, os

master_phil = iotbx.phil.parse("""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.  Repeat to \
    specify additional directories.
run_no = None
  .type = str
  .help = Run no. is used as folder name that stores output files.
  .optional = False
title = None
  .type = str
  .help = Title of the run.
  .multiple = False
scale
  .help = "Parameters used to generate mean-intensity scaled reference set."
{
  d_min = 0.1
    .type = float
    .help = Minimum resolution.
  d_max = 99
    .type = float
    .help = Maximum resolution.
  sigma_min = 1.5
    .type = float
    .help = Minimum I/sigI cutoff.
}
postref
  .help = Post-refinement parameters
{
  scale
    .help = Scale factors
  {
    d_min = 0.1
      .type = float
      .help = Minimum resolution.
    d_max = 99
      .type = float
      .help = Maximum resolution.
    sigma_min = 1.5
      .type = float
      .help = Minimum I/sigI cutoff.
    partiality_min = 0.1
      .type = float
      .help = Minimum partiality cutoff.
  }
  crystal_orientation
    .help = Crystal orientations
  {
    flag_on = True
      .type = bool
      .help = Set to False to turn post-refinement in this section off.
    d_min = 0.1
      .type = float
      .help = Minimum resolution.
    d_max = 99
      .type = float
      .help = Maximum resolution.
    sigma_min = 1.5
      .type = float
      .help = Minimum I/sigI cutoff.
    partiality_min = 0.1
      .type = float
      .help = Minimum partiality cutoff.
  }
  reflecting_range
    .help = Reflecting range
  {
    flag_on = True
      .type = bool
      .help = Set to False to turn post-refinement in this section off.
    d_min = 0.1
      .type = float
      .help = Minimum resolution.
    d_max = 99
      .type = float
      .help = Maximum resolution.
    sigma_min = 1.5
      .type = float
      .help = Minimum I/sigI cutoff.
    partiality_min = 0.1
      .type = float
      .help = Minimum partiality cutoff.
  }
  unit_cell
    .help = Unit-cell dimensions
  {
    flag_on = True
      .type = bool
      .help = Set to False to turn post-refinement in this section off.
    d_min = 0.1
      .type = float
      .help = Minimum resolution.
    d_max = 99
      .type = float
      .help = Maximum resolution.
    sigma_min = 1.5
      .type = float
      .help = Minimum I/sigI cutoff.
    partiality_min = 0.1
      .type = float
      .help = Minimum partiality cutoff.
    uc_tolerance = 3
      .type = float
      .help = Unit-cell tolerance in percent.
  }
  allparams
    .help = All parameters
  {
    flag_on = False
      .type = bool
      .help = Set to True to refine all parameters together.
    d_min = 0.1
      .type = float
      .help = Minimum resolution.
    d_max = 99
      .type = float
      .help = Maximum resolution.
    sigma_min = 1.5
      .type = float
      .help = Minimum I/sigI cutoff.
    partiality_min = 0.1
      .type = float
      .help = Minimum partiality cutoff.
    uc_tolerance = 3
      .type = float
      .help = Unit-cell tolerance in percent.
  }
}
merge
  .help = "Parameters used in merging"
{
  d_min = 0.1
    .type = float
    .help = Minimum resolution.
  d_max = 99
    .type = float
    .help = Maximum resolution.
  sigma_min = 1.5
    .type = float
    .help = Minimum I/sigI cutoff.
  partiality_min = 0.1
    .type = float
    .help = Minimum partiality cutoff.
  uc_tolerance = 3
    .type = float
    .help = Unit-cell tolerance in percent.
}
target_unit_cell = None
  .type = unit_cell
  .help = Target unit-cell parameters are used to discard outlier cells.
  .optional = False
target_space_group = None
  .type = str
  .help = Target space group.
  .optional = False
target_pointgroup = None
  .type = str
  .help = Target point group.
target_anomalous_flag = False
  .type = bool
  .help = Target anomalous flag (False = Not anomalous data)
  .optional = False
indexing_ambiguity
  .help = "Parameters used in resolving indexing ambiguity"
{
  flag_on = False
    .type = bool
    .help = Set to True to allow the program to read in polarity info. \
      from the pickle file specified in index_basis_in.
  index_basis_in = None
    .type = path
    .help = Pickle file storing polarity info. (output from Brehm & \
      Diederichs clustering algorithm).
}
hklisoin = None
  .type = path
  .help = Mtz file for the calculation of CCiso
hklrefin = None
  .type = path
  .help = Mtz file used as a reference in post-refinement.
flag_plot = False
  .type = bool
  .help = Normal plots.
flag_plot_expert = False
  .type = bool
  .help = Expert plots.
flag_force_no_postrefine = False
  .type = bool
  .help = Set to True to output only the mean-intensity scaled merged.
n_postref_cycle = 3
  .type = int
  .help = No. of cycles for post-refinement.
n_postref_sub_cycle = 3
  .type = int
  .help = No. of cycles for the least-squares minimization in post-refinement.
n_bins = 20
  .type = int
  .help = No. of bins used to report statistics.
pixel_size_mm = None
  .type = float
  .help = Pixel size in mm. (MAR = 0.079346)
  .optional = False
frame_accept_min_cc = 0.25
  .type = float
  .help = CC cut-off for the rejection of frames before merging.
n_residues = None
  .type = int
  .help = No. of residues in asymmetric unit.
  .optional = False
flag_apply_b_by_frame = False
  .type = bool
  .help = Set to True to apply B-factor on each frame.
wilson_b_max = 999
  .type = float
  .help = Maximum B-factor (calculated on each still) allowed.
b_refine_d_min = 99
  .type = float
  .help = Minimum resolution.
partiality_model = Lorentzian
  .type = str
  .help = Your choice of partiality model: Lorentzian (default), Disc, Kabsch, or Rossmann.
flag_LP_correction = False
  .type = bool
  .help = Do Lorentz-factor and polarization correction.
n_processors = 32
  .type = int
  .help = No. of processing units
  .optional = True
gamma_e = 0.002
  .type = float
  .help = Initial spread of the energy spectrum (1/Angstrom).
polarization_horizontal_fraction = 1.0
  .type = float
  .help = Polarization fraction in horizontal direction.
""")


def process_input(argv=None):
  if argv == None:
    argv = sys.argv[1:]

  user_phil = []
  for arg in sys.argv[1:]:
    if os.path.isfile(arg):
      user_phil.append(iotbx.phil.parse(open(arg).read()))
    elif (os.path.isdir(arg)) :
      user_phil.append(iotbx.phil.parse("""data=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(iotbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  if (len(params.data) == 0) or False in [os.path.isdir(p) for p in params.data]:
    master_phil.show()
    raise Usage("No data")

  #generate run_no folder
  if not os.path.exists(params.run_no):
    os.makedirs(params.run_no)

  #capture input read out by phil
  from cStringIO import StringIO
  class Capturing(list):
    def __enter__(self):
      self._stdout = sys.stdout
      sys.stdout = self._stringio = StringIO()
      return self
    def __exit__(self, *args):
      self.extend(self._stringio.getvalue().splitlines())
      sys.stdout = self._stdout

  with Capturing() as output:
    working_phil.show()

  txt_out = 'prime.postrefine input:\n'
  for one_output in output:
    txt_out += one_output + '\n'

  return params, txt_out
