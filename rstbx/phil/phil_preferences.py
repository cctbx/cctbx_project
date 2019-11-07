from __future__ import absolute_import, division, print_function
from libtbx.phil.command_line import argument_interpreter as model_argument_interpreter
from libtbx.utils import Sorry
from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs # implicit import

try:
  from dials.algorithms.integration.kapton_correction import absorption_defs
except ImportError:
  absorption_defs = ""

libtbx_misc_defs = """
predictions_file = ""
    .type = str
    .help = File has xds parameters for spot predictions in XDS XPARAM format.

parallel = 0
  .type = int

integration {
  file_template = None
    .type = str
    .help = "Full path for the files to integrate, expressed as template like lyso_1_###.img"
  file_range = None
    .type = ints (value_min=1)
    .help = First and last file number to integrate, forming a contiguous INCLUSIVE range (not like Python range).
  rocking_curve = *none gh1982a
    .type = choice
    .help = gh1982a is the Greenhough & Helliwell 1982 section I model, with epsilon from eqn (V.6)
  mosaicity_deg = 0.0
    .type = float
    .help = full width effective mosaicity, degrees, for specified rocking curve model such as gh1982a
  guard_width_sq = 11
    .type = int
    .help = Guard is the reserved area around the Bragg spot mask that cannot be used for background
    .help = plane determination because the tail of the signal distribution may leak from the spot.
    .help = Value represents max disallowed squared hypotenuse between neighboring signal & background pixels, as integer px**2.
  detector_gain = 0.0
    .type = float
    .help = Detector gain in units of Analog Digital Units per photon.
    .help = Future plan: this value overrides any given through the dxtbx format mechanism.
    .help = Present: this is a mandatory value, code throws an exception with the default value.
  background_factor = 1
    .type = int (value_min=1)
    .help = require minimum number of pixels for background fit = background_factor x # spot pixels
  model = *rossmann1979jac12-225 use_case_3_simulated_annealing use_case_3_simulated_annealing_7 use_case_3_simulated_annealing_9 user_supplied
    .type = choice
    .help = algorithm for prediction of spots
    .help = Michael Rossman (1979) J. Appl. Cryst. 12, 225-238.
  use_subpixel_translations = None
    .type = floats
    .help = list slow,fast offsets for correcting tile positions to subpixel resolution (2 numbers for each tile)
  subpixel_joint_model {
    rotations = None
      .type=floats
      .help = joint-refined tile rotations & translations [along with per-image beam,dist,rotz,wavelength, not stored]
    translations = None
      .type=floats
      .help = joint-refined tile rotations & translations [along with per-image beam,dist,rotz,wavelength, not stored]
  }
  spot_shape_verbose = False
    .type = bool
    .help = analysis of radial and azimuthal spot shapes.
  signal_penetration = 0.5
    .type = float
    .help = For computing parallax effect due to finite sensor thickness, fraction of signal to attenuate before
    .help = ignoring the remaining trailing parallax.  Small value (0.0) means do not account for parallax.
  spotfinder_subset = *inlier_spots goodspots spots_non-ice
    .type = choice
    .help = which subset to use for parameter refinement and constructing integration profiles.
    .help = subsets are nested goodspots > spots_non-ice > inlier_spots
  mask_pixel_value = None
    .type = int
    .help = pixels set to this value will be ignored during integration
  mosaic {
    refinement_target = LSQ *ML
      .type = choice
      .help = Modeling of still shots, refinement of crystal physical parameters by one of two targets given in the paper
      .help = [Sauter et al. Acta D (2014) 70:3299-3309]: least squares (LSQ) or maximum likelihood (ML).
      .help = As the paper describes, ML gives superior results and is set to be the default as of 7/2017.
    kludge1 = 1
      .type = float
      .help = bugfix 1 of 2 for protocol 6, equation 2.  When matching obs and model reflections,
      .help = must increase the provisional model mosaicity so as to model enough spots.
      .help = default 1 might be too small for a bad model, aborting indexing with too few matches.
      .help = > 2 might be too large and could lead to misindexing.
    bugfix2_enable = True
      .type = bool
      .help = bugfix 2 of 2 for protocol 6, equation 2.  Assure that the result of the arctan()
      .help = function call is placed in the correct principle value region.
      .help = choice of "False" preserves the ability to roll back program behavior
      .help = as promised to the Nature Methods editor.
    domain_size_lower_limit = 10
      .type = float
      .help = Mosaic blocks must be at least this many unit cells on edge (refers to cube-root of cell volume)
      .help = used in simple_integration.py as an initial value
      .help = too-small value predicts too many spots, leading to misindexing.
      .help = too-large value doesn't predict enough spots
      .help = choose lowest value physically reasonable, 10 unit cells.
    enable_rotational_target_highsym = True
      .type = bool
      .help = Protocol 6, use equation 2 (true) or protocol 5, use equation 1 (false) for higher Bravais-setting refinement
    enable_rotational_target_triclinic = True
      .type = bool
      .help = Protocol 5, use equation 2 (true) or protocol 4, use equation 1 (false) for triclinic Bravais-setting refinement
    enable_simplex = False
      .type = bool
      .help = simplex_optimization
    enable_AD14F7B = False
      .type = bool
      .help = enable Sauter et al. Acta D (2014) 70:3299-3309, Fig. 7(b) plot, angular excursion vs two theta
    enable_polychromatic = False
      .type = bool
      .help = followup top hat wavelength dispersion
    ewald_proximal_volume_resolution_cutoff = 2.5
      .type = float
      .help = Resolution cutoff in Angstroms.  Set roughly to the apparent cutoff of bright spots for one or a group of images.
      .help = governs the cutoff resolution for teh ewald proximal volume: that volume of reciprocal space containing
      .help =  spots predicted by the current mosaicity model.
  }
  initial_volume_factor = 1.8
      .type = float
      .help = controls the initial integration limits before filtering by sig(I) durig merging
      .help = used in rstbx/new_horizons/oscillation_shots.py -- not fully understood
      .help = a larger volume factor integrates a larger volume of reciprocal space
  enable_residual_map = False
      .type = bool
      .help = x,y model vs. spotfinder residuals plotted vs. image position
  enable_residual_map_deltapsi = False
      .type = bool
      .help = if residual map is enabled, colorcode the Bragg spots by Delta-Psi
      .help = blue, negative delta psi, outside Ewald sphere. red, positive delta psi, inside Ewald sphere.
      .help = requires also setting enable_residual_map to True and spot_prediction to dials
  enable_residual_scatter = False
      .type = bool
      .help = x,y model vs. spotfinder residuals scatter plot
  graphics_backend = *screen pdf
      .type = choice
      .help = data plots are written to the terminal screen, or to pdf files
  pdf_output_dir = "."
      .type = str
      .help = plots are written to this directory
  montecarlo_integration_limit = None
    .type = float
    .help = use None to limit the integration resolution based on Wilson plot
    .help = alternatively set integration limit to a fixed value given in Angstroms.
  greedy_integration_limit = False
    .type = bool
    .help = governs the resolution limits used for higher Bravais settings
    .help = if true, use the expanded limiting_resolution used for the last triclinic round
    .help = if false, attempt to analyze triclinic integrated spots for limits (default until Feb 2015, not so reliable)
  combine_sym_constraints_and_3D_target = False
    .type = bool
    .help = enable redesigned refinement protocol from Acta D 2014 Sauter paper
    .help = if false (default as of Feb 2015) symmetry constraints are applied after eqn (1) positional refinement
    .help =   but before indexing of high-Bravais symmetry spots, possibly leading to misindexing of high-angle spots
    .help = if true, spots are indexed once only, in triclinic setting.  After application of high-symmetry constraints
    .help =   dials is used to refine positions (eqn 1) and deltapsi angle (eqn 2).
  spot_prediction = *ucbp3 dials
      .type = choice
      .help = in high-symmetry integration protocol, predict with ucbp3 (CSPAD subpixel corrections) or dials (detector tilt, but unit translations)
  enable_one_to_one_safeguard = False
    .type = bool
    .help = Flag enables a safeguard within the stills integration algorithm, within the mapping of
    .help = predicted spots to observations prior to model refinement.  The safeguard prevents
    .help = a many-to-one predicted-to-observation mapping, assuring that each prediction is mapped
    .help = to at most one observation, one that is closest.  The legacy code did not include
    .help = this safeguard.  It is expected that True should become the default after a testing period.
  dials_refinement{
    strategy = *distance wavelength fix
      .type = choice
      .help = either refine the distance or the wavelength for XFEL stills, or fix them both in place
  }%s
}
""" % absorption_defs
indexing_defs = """
include scope spotfinder.command_line.signal_strength.master_params
spotfinder = *distl speck
  .type=choice
  .help = "Choose among spotfinder implementations [distl|speck]"

speckfinder {

  dark_stddev = ""
    .type = str
    .help = Mandatory dark standard deviation image for gain correction.
  dark_adu_scale = 100
    .type = int
    .help = "Mandatory scale at which dark was calculated; must be >1 on account of integer rounding."
}

indexing {
  data = None
    .type=str
    .multiple=True
    .help="Relative or absolute path names for raw image files to be indexed"
  indexing_pickle = None
    .type=str
    .help = "pickle file name for integration results subsequent to indexing."
  completeness_pickle = None
    .type=str
    .help = "pickle file name for HKL, I, SIGI, XY."
  open_wx_viewer = False
    .type = bool
  verbose_cv = False
    .type = bool
    .help = "screen printout of the obs vs predicted spot correction vectors,"
    .help = "for empriical repositioning of the detector tiles."
  lattice_model_scoring_cutoff = 2.0
    .type = float
    .help = Cutoff value for the <Z-score> over integrated signal from the model lattice.
    .help = Used for choosing the most accurate combination of candidate basis vectors.
  devel_algorithm = None
    .type = str
    .help = for development only, turn on whatever testing behavior

  outlier_detection {
    allow = True
      .type = bool
      .multiple=False
      .help="Algorithm (Sauter&Poon[2010] J Appl Cryst 43:611) provides superior positional fit with noisy data."
    switch=False
      .type=bool
      .multiple=False
      .help="Switch to the outlying spots to detect a second lattice. False==first lattice; True==second lattice"
    verbose=False
      .type=bool
      .multiple=False
      .help="Verbose output."
    pdf=None
      .type=str
      .multiple=False
      .help="Output file name for making graphs of |dr| vs spot number and dy vs dx."
  }
  plot_search_scope = False
    .type = bool
    .help = improvement of the model, plot target function of origin offset or S0
  mm_search_scope = 4.0
    .type = float
    .help = global radius of origin_offset search, used for plotting the search scope
  improve_local_scope = *origin_offset S0_vector
    .type = choice
    .help = improve 'beam position' according to Sauter et al (2004).  Local minimum only
    .help = specifies which parameter to optimize.
}
"""

iotbx_defs_viewer_detail = """
  powder_arcs{
    show = False
      .type=bool
      .help = "show powder arcs calculated from PDB file."
    code = None
      .type=str
      .help = "PDB code (4 characters) for file; fetch it from the Internet."
  }
  calibrate_silver = False
      .type=bool
      .help = "Open special GUI for distance/metrology from silver behenate."
  calibrate_pdb{
    code = None
      .type=str
      .help = "Open pdb code (over Internet) to get unit cell & symmetry for powder rings."
      .help = "Most useful for calibrating low-Q rings on far detector."
      .help = "Option is mutually exclusive with calibrate silver, unit cell and powder arcs options."
    d_min = 20.
      .type=float
      .help = "Limiting resolution to calculate powder rings"
  }
  calibrate_unitcell{
    unitcell = None
      .type=unit_cell
      .help = "Specify unit cell for powder rings."
      .help = "Option is mutually exclusive with calibrate silver, pdb and powder arcs options."
    d_min = 20.
      .type=float
      .help = "Limiting resolution to calculate powder rings"
    spacegroup = None
      .type=str
      .help = "Specify spacegroup for the unit cell"
  }
"""

iotbx_defs_viewer = """
viewer {
  %s
}
"""%iotbx_defs_viewer_detail

iotbx_defs_target = """
target_cell=None
  .type=unit_cell
  .multiple=False
  .help="Imperative unit cell applied at the level of DPS algorithm basis selection."
target_cell_centring_type= *P C I R F
  .type=choice
  .multiple=False
  .help="Centring symbol for the target cell"
isoforms
  .help=Constrain the unit cell to specific values during refinement
  .help=As presently implemented, applies only to dials_refinement_preceding_integration
  .help=and applies only to the higher-symmetry integration trial, not the initial triclinic
  .multiple=True
  {
    name=None
      .type=str
    cell=None
      .type=unit_cell
    lookup_symbol=None
      .type=str
      .help=The sgtbx lookup symbol of the reflections pointgroup
    rmsd_target_mm=None
      .type=float
      .help=Maximum acceptable DIALS positional rmsd, in mm
    beam_restraint=None
      .type=floats(size=2)
      .help=Known beam position in mm X,Y, rmsd_target_mm is reused here as a circle of confusion
      .help=to assure that no images are accepted where the lattice is misindexed by a unit shift.
  }
"""
libtbx_defs = indexing_defs + libtbx_misc_defs
iotbx_defs = iotbx_defs_viewer + iotbx_defs_target
indexing_api_defs = indexing_defs + iotbx_defs_target

class EffectiveParamGenerator:
  def __init__(self,libtbx_defs,iotbx_defs):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

  def master(self,package = 'iotbx'):
    libselector = {'libtbx':self.libtbx_defs,
                   'iotbx':self.iotbx_defs+self.libtbx_defs,
                  } [package]
    if (package == "libtbx"):
      from libtbx import phil
    else:
      from iotbx import phil
    return phil.parse(input_string=libselector, process_includes=True)

  def default(self,item = 'iotbx'):
    app_master = self.master(item)
    return app_master.fetch(sources=[app_master,])

  def show(self, modpython):
    modified_params = self.master().format(python_object = modpython)
    modified_params.show()

  def merge(self, args):
    #future:  this member function should be deprecated; replace with preferences.py

    effective_params = self.default()

    argument_interpreter = model_argument_interpreter(
      master_phil=self.master(),
      #home_scope =
    )
    consume = []
    for arg in args:

      try:
        command_line_params = argument_interpreter.process(
          arg=arg
        )
        effective_params = effective_params.fetch(sources=[command_line_params,])
        consume.append(arg)

      except Sorry as e:
        pass

    for item in consume:
      args.remove(item)


    # effective_params.show()

    params = effective_params.extract()

    self.validation(params)

    self.effective_params = effective_params
    return params

  def validation(self,trial_params):
    pass

effective_param_generator = EffectiveParamGenerator(libtbx_defs,iotbx_defs) #singleton
