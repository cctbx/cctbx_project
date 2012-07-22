from libtbx.phil.command_line import argument_interpreter as model_argument_interpreter
from libtbx.utils import Sorry
from spotfinder.command_line.signal_strength import additional_spotfinder_phil_defs # implicit import

libtbx_defs = """\

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

predictions_file = ""
    .type = str
    .help = File has xds parameters for spot predictions in XDS XPARAM format.

parallel = 0
  .type = int

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
}

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
  background_factor = 1
    .type = int (value_min=1)
    .help = require minimum number of pixels for background fit = background_factor x # spot pixels
  model = *rossmann1979jac12-225 use_case_3_simulated_annealing use_case_3_simulated_annealing_7 use_case_3_simulated_annealing_9
    .type = choice
    .help = algorithm for prediction of spots
    .help = Michael Rossman (1979) J. Appl. Cryst. 12, 225-238.
  use_subpixel_translations = None
    .type = floats
    .help = list slow,fast offsets for correcting tile positions to subpixel resolution (2 numbers for each tile)
  spot_shape_verbose = False
    .type = bool
    .help = analysis of radial and azimuthal spot shapes.
}

viewer {
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
}
"""

iotbx_defs = """\
"""

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

      except Sorry,e:
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
