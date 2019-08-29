from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME distl.signal_strength
import spotfinder
import libtbx.phil
from libtbx.utils import Sorry
import sys, os
additional_spotfinder_phil_defs ="""
distl {
  minimum_spot_area = None
    .type = int
    .help = "Override application default (PADs:5, others:10) set minimum spot area (in pixels) within spotfinder."
  minimum_signal_height = None
    .type = float
    .help = "Override application default (PADs:2.5, CCDs:1.5) set minimum signal height (in units of background noise sigma) within spotfinder."
  minimum_spot_height = None
    .type = float
    .help = "Expert use only; after pixels are classified as signals rather than noise, minimum height to be considered a spot maximum (in units of background noise sigma). Default=3.5"
  spot_area_maximum_factor = None
    .type = float
    .help = "Expert use only; max spot area expressed as a multiple of minimum_spot_area. Default=5.0"
  peak_intensity_maximum_factor = None
    .type = float
    .help = "Expert use only; max peak intensity filter for use by libdistl. Default=10.0"
  method2_cutoff_percentage = 20
    .type = float
    .help = As described in Zhang(2006) spotfinder paper, the resolution cutoff is made at a shell that has this percentage of
    .help = spots as compared to the two inner shells.  Needs to be set to a lower percentage for still (XFEL) data!
    .help = Low percentage (5%% as recommended for XFEL) is more permissive and accepts higher resolution Bragg spots.
  compactness_filter = False
    .type = bool
    .help = "For CXI data (and Pilatus data?), set this to True to eliminate small non-compact spots"
  detector_tiling = None
    .type = ints
    .help = "Prior analysis of active area tiling on the detector; groups of 4 integers: UL_slow,UL_fast,LR_slow,LR_fast"
    .help = "...where UL=upper left corner, LR=lower right corner; coordinates form a Python range"
    .expert_level = 2
  tile_translations = None
    .type = ints
    .help = "For each identified tile, a group of 2 integers for slow & fast translations to correct position,"
    .help = "given in pixels"
    .expert_level = 2
  tile_flags = None
    .type = ints (value_min=0,value_max=1)
    .help = "For each identified tile, an integer [0,1] acting as a bool to define whether to ignore its data."
    .expert_level = 2
  quad_translations = None
    .type = ints
    .help = "For quadrants UL,UR,LL,LR groups of 2 integers for slow & fast translations to correct position,"
    .help = "given in pixels"
    .expert_level = 2
  detector_format_version = None
    .type = str
    .help = "Whatever additional text information is necessary to get an accurate file read"
    .expert_level = 2
  %s
  pdf_output = None
    .type = str
    .help="File name for optional PDF graphical output for distl.signal_strength (*.pdf)"
  image_viewer = False
    .type = bool
    .expert_level=2
    .help="Open the image viewer to inspect the spotfinder spots."
  port = 8125
    .type = int
    .help="For the server version, port number to listen for requests"
  processors = 1
    .type = int
    .help="For the multithreaded server version, number of server processes to be used"
  nproc = 1
    .type = int
    .short_caption = Number of processes
    .help="Number of processes to be used when processing multiple images."
}
""" % spotfinder.inner_phil_str

master_params_defs = """\
distl {
  image = None
    .type = str
    .help="Image file name"
  res {
    outer=None
      .type=float
      .help="High resolution limit in angstroms"
    inner=None
      .type=float
      .help="Low resolution limit in angstroms"
  }
  verbose = False
    .type = bool
    .help="Lengthy spot printout"
  dxtbx = False
    .type = bool
    .help="Switch algorithms to dxtbx models"
  bins {
    verbose = False
      .type = bool
      .help="Additional printout binned by resolution range"
    N = 20
      .type = int
      .help="Maximum number of bins, but fewer can result if there are few spots"
    corner = True
      .type = bool
      .help="Extend the binning all the way to detector corner, otherwise to outermost spot on first image"
  }
  range = None
    .type = ints (value_min = 0, size_max = 2)
    .help = "For HDF5 data, specifies image index (0-based), or Python-style range (0,10 means first 10)"
}

%s
"""%(spotfinder.labelit_related_commands)+additional_spotfinder_phil_defs

master_params = libtbx.phil.parse(master_params_defs)

def run(args, command_name="distl.signal_strength"):
  help_str="""explanation:
Local background and background standard deviation are determined.
Pixels are classified as signal if > minimum_signal_height sigmas above background.
Signals are chosen as spot maxima if > minimum_spot_height sigmas above background.
Spots are grown around the maxima, and retained if they fit minimal area criteria.
Total number of candidates at this stage is reported as "Spot Total"
Ice rings are eliminated by at least two different algorithms (rings of high pixel
  values and rings of high spot count).
Resolution filters are applied if given on the command line.
Total number of candidates at this stage is reported as "In-Resolution Total"
Other spot-quality filters are applied to give the number of "Good Bragg Candidates".
Method 1 Resolution is a published legacy algorithm (Zhang et al, 2006) no longer used.
Method 2 Resolution reflects drop off of spot count as a function of resolution shell,
  but is overridden by command line input of distl.res.outer
Signal strength of the Good Bragg Candidates is then presented as integrated area of
  the spot above local background, expressed in pixel-analog/digital units.
Very verbose output is available by setting distl.verbose=True
Full documentation: http://cci.lbl.gov/publications/download/ccn_jul2010_page18.pdf
"""

  if (len(args) == 0 or args[0] in ["H","h","-H","-h","help","--help","-help"]):
    print("usage:   %s image_filename [parameter=value ...]" % command_name)
    print("example: %s lysozyme_001.img distl.res.outer=2.0 distl.res.inner=6.0 distl.minimum_spot_area=8"%command_name)
    master_params.show(attributes_level=1,expert_level=1)
    print(help_str)
    return

  print("%s: characterization of candidate Bragg spots"%command_name)

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")
  image_file_name = None
  moving_pdb_file_name = None
  for arg in args:
    if (os.path.isfile(arg)):
      if (image_file_name is None): image_file_name = arg
      else: raise Sorry("Too many file names.")
    else:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  def raise_missing(what):
      raise Sorry("""\
Missing file name for %(what)s structure:
  Please add
    %(what)s=file_name
  to the command line to specify the %(what)s structure.""" % vars())

  if (image_file_name is None):
    if (params.distl.image is None): raise_missing("file name")
  else:
    params.distl.image = image_file_name

  print("#Parameters used:")
  print("#phil __ON__")
  print()
  working_params = master_params.format(python_object=params)
  working_params.show(expert_level=1)
  print()
  print("#phil __OFF__")
  print()

  #Now actually run the program logic
  from spotfinder.applications import signal_strength
  return signal_strength.run_signal_strength(working_params.extract())

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
