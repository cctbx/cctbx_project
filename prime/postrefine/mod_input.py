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
icering
  .help = "Allowing exclusion of icering."
{
  flag_on = False
    .type = bool
    .help = Turn this flag on to allow exclusion of icering.
  d_upper = 3.9
    .type = float
    .help = Minimum resolution.
  d_lower = 3.85
    .type = float
    .help = Maximum resolution.
}
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
  residual_threshold = 5
    .type = float
    .help = Percent increase in residual allowed during microcycle.
  residual_threshold_xy = 5
    .type = float
    .help = Percent increase in residual (xy) allowed during microcycle.
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
flag_override_unit_cell = False
  .type = bool
  .help = Set to True to overide unit cell in observations with the target cell.
target_space_group = None
  .type = str
  .help = Target space group.
  .optional = False
target_anomalous_flag = False
  .type = bool
  .help = Target anomalous flag (False = Not anomalous data)
  .optional = False
flag_weak_anomalous = False
  .type = bool
  .help = Set to True to indicate that you have a very weak anomalous signal.
target_crystal_system = None
  .type = str
  .help = Target crystal system
  .optional = True
n_residues = None
  .type = int
  .help = No. of amino acid residues.
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
flag_outlier_rejection = True
  .type = bool
  .help = Set to True to perform advance outlier rejection.
n_postref_cycle = 3
  .type = int
  .help = No. of cycles for post-refinement.
n_postref_sub_cycle = 3
  .type = int
  .help = No. of cycles for the least-squares minimization in post-refinement.
n_rejection_cycle = 3
  .type = int
  .help = No. of cycles for the outlier rejection.
sigma_rejection = 3
  .type = float
  .help = Sigma level for outlier rejection.
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
flag_apply_b_by_frame = True
  .type = bool
  .help = Set to False to dismiss B-factor checking.
b_refine_d_min = 99
  .type = float
  .help = Minimum resolution.
partiality_model = Lorentzian
  .type = str
  .help = Your choice of partiality model: Lorentzian (default), Disc, Kabsch, or Rossmann.
flag_LP_correction = True
  .type = bool
  .help = Do polarization correction.
flag_volume_correction = True
  .type = bool
  .help = Do volume correction.
flag_beam_divergence = False
  .type = bool
  .help = Default is not to refine beam divergence. Set to True to allow gammaX and gammaY refinement.
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
flag_output_verbose = False
  .type = bool
  .help = Output full detail of the refinement results.
flag_replace_sigI = False
  .type = bool
  .help = Replace to experimental errors I with sqrt(I).
percent_cone_fraction = 5.0
  .type = float
  .help = Perent used to select reflections inside a cone centering on each crystal axis.
reindex_op = h,k,l
  .type = str
  .help = Change basis for the input observations.
reindex_apply_to = None
  .type = str
  .help = Space group which reindex will be performed on.
  .optional = False
flag_normalized = False
  .type = bool
  .help = Set to True for amplitude normalization.
sigma_global_min = 1.0
  .type = float
  .help = Minimum I/sigI cutoff (global value).
iotacc
  .help = "Parameters used in prime.iotacc selection."
{
  set_id = None
    .type = str
    .help = Crystal set id. Set to None for one single set.
  n_shots = 99
    .type = int
    .help = Maximum no. of shots per crystal.
    .optional = True
  d_min = 0.1
    .type = float
    .help = Minimum resolution.
  d_max = 99
    .type = float
    .help = Maximum resolution.
  uc_tolerance = 3
    .type = float
    .help = Unit-cell tolerance in percent.
  LEN_SHOT_SEQ = 5
    .type = int
    .help = No. of digits used to format shot no.
    .optional = True
  sigma_min = 1.5
    .type = float
    .help = Minimum I/sigI cutoff.
  percent_top_cc = 10
    .type = float
    .help = Percent of best cc.
    .optional = True
}
""")

txt_help = """**************************************************************************************************

Prime: post-refinement and merging.

For more detail and citation, see Enabling X-ray free electron laser crystallography
for challenging biological systems from a limited number of crystals
"DOI: http://dx.doi.org/10.7554/eLife.05421".

Usage: prime.postrefine parameter.phil

With this command, you can specify all parameters required by prime in your parameter.phil file.
To obtain the template of these parameters, you can perform a dry run (simply run prime.postrefine).
You can then change the values of the parameters.

For feedback, please contact monarin@stanford.edu.

**************************************************************************************************

List of available parameters:
"""

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
      print arg
      if arg == '--help' or arg == '-h':
        print txt_help
        master_phil.show(attributes_level=1)
        exit()
      try :
        user_phil.append(iotbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  if (len(params.data) == 0):
    master_phil.show()
    raise Usage("Use the above list of parameters to generate your input file (.phil). For more information, run prime.postrefine -h.")

  #generate run_no folder
  if os.path.exists(params.run_no):
    raise Sorry("The run number %s already exists."%params.run_no)

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
