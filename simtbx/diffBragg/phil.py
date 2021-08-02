from __future__ import absolute_import, division, print_function
from iotbx.phil import parse

#help_message = '''
#diffBragg command line utility
#'''

simulator_phil = """
simulator {
  oversample = 0
    .type = int
    .help = pixel oversample rate (0 means auto-select)
  device_id = 0
    .type = int
    .help = device id for GPU simulation
  init_scale = 1
    .type = float
    .help = initial scale factor for this crystal simulation
  total_flux = 1e12
    .type = float
    .help = total photon flux for all energies
  crystal {
    ncells_abc = (10,10,10)
      .type = floats(size=3)
      .help = number of unit cells along each crystal axis making up a mosaic domain
    ncells_def = (0,0,0)
      .type = floats(size=3)
      .help = off-diagonal components for mosaic domain model (experimental)
    has_isotropic_ncells = False
      .type = bool
      .help = if True, ncells_abc are constrained to be the same values during refinement
    mosaicity = 0
      .type = float
      .help = mosaic spread in degrees
    anisotropic_mosaicity = None
      .type = floats
      .help = mosaic spread 3-tuple or 6-tuple specifying anisotropic mosaicity
    num_mosaicity_samples = 1
      .type = int
      .help = the number of mosaic domains to use when simulating mosaic spread
    crystal_for_anisotropic_mosaicity = None
      .type = str
      .help = path to a experiment file holding a crystal (if None, defaults to experment crystal)
    mos_angles_per_axis = 10
      .type = int
      .help = if doing a uniform mosaicity sampling, use this many angles per rotation axis
    num_mos_axes = 10
      .type = int
      .help = number of sampled rot axes if doing a uniform mosaicity sampling
    mosaicity_method = 2
      .type = int
      .help = 1 or 2. 1 is random sampling, 2 is even sampling
    num_sausages = 1
      .type = int
      .help = the number of blocks that make up the mosaic texture
    rotXYZ_ucell = None
      .type = floats(size=9)
      .help = three missetting angles (about X,Y,Z axes), followed by
      .help = unit cell parameters. The crystal will be rotated according to
      .help = the matrix RotX*RotY*RotZ, and then the unit cell will be updated
  }
  structure_factors {
    mtz_name = None
      .type = str
      .help = path to an MTZ file
    mtz_column = None
      .type = str
      .help = column in an MTZ file
    dmin = 1.5
      .type = float
      .help = minimum resolution for structure factor array
    dmax = 30
      .type = float
      .help = maximum resolution for structure factor array
    default_F = 0
      .type = float
      .help = default value for structure factor amps
  }
  spectrum {
    filename = None
      .type = str
      .help = a .lam file (precognition) for inputting wavelength spectra
    stride = 1
      .type = int
      .help = stride of the spectrum (e.g. set to 10 to keep every 10th value in the spectrum file data)
    filename_list = None
      .type = str
      .help = path to a file containing 1 .lam filename per line
  }
  beam {
    size_mm = 1
      .type = float
      .help = diameter of the beam in mm
  }
  detector {
    force_zero_thickness = False
      .type = bool
      .help = if True, then set sensor thickness to 0
  }
}
"""

refiner_phil = """
refiner {
  break_signal = None
    .type = int
    .help = intended to be used to break out of a long refinement job prior to a timeout on a super computer
    .help = On summit, set this to 12 (SIGUSR2), at least thats what it was last I checked (July 2021)
  debug_pixel_panelfastslow = None
    .type = ints(size=3)
    .help = 3-tuple of panel ID, fast coord, slow coord. If set, show the state of diffBragg
    .help = for this pixel once refinement has finished
  res_ranges = None
    .type = str
    .help = resolution-defining strings, where each string is
    .help = is comma-separated substrings, formatted according to "%f-%f,%f-%f" where the first float
    .help = in each substr specifies the high-resolution for the refinement trial, and the second float
    .help = specifies the low-resolution for the refinement trial. Should be same length as max_calls
  parameter_hdf5_path = None
    .type = str
    .help = path to an output file where parameter log can be written (stage 1 only)
  dump_params = False
    .type = bool
    .help = save params to the parameter_hdf5_path folder
  n_sausages_one_scale = False
    .type = bool
    .help = force all mosaic textures blocks to have the same scale factor
  print_end = None
    .type = str
    .help = append this string to the end of each print statement during refinement
  mask = None
    .type = str
    .help = path to a dials mask flagging the trusted pixels
  force_symbol = None
    .type = str
    .help = a space group lookup symbol used to map input miller indices to ASU
  force_unit_cell = None
    .type = ints(size=6)
    .help = a unit cell tuple to use
  use_perpixel_dark_rms = False
    .type = bool
    .help = if True, try to exrtact per pixel dark data from image format
    .help = (only supported for Jungfrau16M formats)
    .help = see dxtbx FormatNexusJungfrauResRaw
  use_cuda = False
    .type = bool
    .help = Use GPU acceleration during refinement (requires a CCTBX cuda build)
  randomize_devices = True
    .type = bool
    .help = randomly select a device Id
  num_devices = 1
    .type = int
    .help = number of cuda devices on current node
  refine_blueSausages = None
    .type = ints(size_min=1)
    .help = whether to refine the mosaic texture (blue sausages)
  refine_eta = None
    .type = ints(size_min=1)
    .help = whether to refine the mosaic spread parameter, eta
  refine_Fcell = None
    .type = ints(size_min=1)
    .help = whether to refine the structure factor amplitudes
  refine_per_spot_scale = None
    .type = ints(size_min=1)
    .help = whether to refine scale factors per spot, per shot
  refine_Bmatrix = None
    .type = ints(size_min=1)
    .help = whether to refine the Bmatrix (unit cell parameters)
  refine_Umatrix = None
    .type = ints(size_min=1)
    .help = whether to refine the Umatrix (crystal orientation parameters: rotX, rotY, rotZ)
  refine_ncells = None
    .type = ints(size_min=1)
    .help = whether to refine the ncells_abc (mosaic domain parameters)
  refine_ncells_def = None
    .type = ints(size_min=1)
    .help = whether to refine the off-diagonal ncells components (ncells_def) (mosaic domain parameters)
  refine_bg = None
    .type = ints(size_min=1)
    .help = whether to refine the tilt plane coefficients
  refine_spot_scale = None
    .type = ints(size_min=1)
    .help = whether to refine the crystal scale factor
  refine_spectra = None
    .type = ints(size_min=1)
    .help = whether to refine the two spectra coefficients:  Lambda -> a*Lambda + b
  refine_panelRotO = None
    .type = ints(size_min=1)
    .help = whether to refine rotation of detector panels about their normal vectors
  refine_panelRotF = None
    .type = ints(size_min=1)
    .help = whether to refine rotation of detector panels about their fast-scan vectors
  refine_panelRotS = None
    .type = ints(size_min=1)
    .help = whether to refine rotation of detector panels about their slow-scan vectors
  refine_panelXY = None
    .type = ints(size_min=1)
    .help = whether to refine the XY component of the detector panel's origin vector
  refine_panelZ = None
    .type = ints(size_min=1)
    .help = whether to refine the Z component of the detector panel's origin vector,
    .help = (1 unique value per panel group, see panel_group_file)
  refine_detdist = None
    .type = ints(size_min=1)
    .help = whether to refine the detector distance (Z), one value refined for all panels
  max_calls = [100]
    .type = ints(size_min=1)
    .help = maximum number of calls for the refinement trial
  fix_rotZ = False
    .type = bool
    .help = whether to fix the rotation of the crystal about the Z direction (usually the beam direction)
  save_models = False
    .type = bool
    .help = whether to save a models file during refinement
  rescale_fcell_by_resolution = False
    .type = bool
    .help = whether to rescale the structure factors according to their resolution
    .help = in an attempt to refine structure factors equally
  panel_group_file = None
    .type = str
    .help = a text file with 2 columns, the first column is the panel_id and the second
    .help = column is the panel_group_id. Panels geometries in the same group are refined together
  init {
    sausages = None
      .type = floats
      .help = the initial parameters for the sausage, should be multiples of (rotX, rotY, rotZ, scale)
      .help = , one for each mosaic block (sausage) making up texture
    eta = 0
      .type = float
      .help = mosaic spread sigma in degrees
    spot_scale = 1
      .type = float
      .help = initial value for spot scale
    ncells_abc = [10, 10, 10]
      .type = floats(size=3)
      .help = initial value for ncells abc that will override params.simulator.crystal.ncells_abc
    ncells_def = (0,0,0)
      .type = floats(size=3)
      .help = initial value for ncells def
    spectra_coefficients = [0, 1]
      .type = floats(size=2)
      .help = initial value for spectrum coefficients  Lambda -> c0 + Lambda*c1
  }
  sensitivity {
    eta = 1
      .type = float
      .help = refinement sensitivity for eta, the mosaic spread parameter
    panelRotOFS = [0.01, 0.01, 0.01]
      .type = floats(size=3)
      .help = refinement sensitivity factor for panel rotations (orthogonal, fast, slow)
    panelXY = [1, 1]
      .type = floats(size=2)
      .help = refinement sensitivity factor for panel XY translations in the lab frame
    panelZ = 1
      .type = float
      .help = refinement sensitivity factor for panel Z translations in the lab frame
    rotXYZ = [0.1, 0.1, 0.1]
      .type = floats(size=3)
      .help = refinement sensitivity factor for rotation parameters
    originZ = 0.1
      .type = float
      .help = refinement sensitivity factor for origin Z parameters
    unitcell = [1, 1, 1, 0.1, 0.1, 0.1]
      .type = floats(size=6)
      .help = unit cell parameter sigma values. All 6 must be present, even
      .help = if the crystal system has higher symmetry.
    ncells_abc = [1, 1, 1]
      .type = floats(size=3)
      .help = refinement sensitivity factor for ncells abc parameters
    ncells_def = [1, 1, 1]
      .type = floats(size=3)
      .help = refinement sensitivity factor for ncells def parameters
    spot_scale = 1
      .type = float
      .help = refinement sensitivity factor for crystal scale parameters
    tilt_abc = [0.1, 0.1, 0.1]
      .type = floats(size=3)
      .help = refinement sensitivity factor for tilt plane coefficients
    fcell = 1
      .type = float
      .help = refinement sensitivity factor for fcell parameters
    spectra_coefficients = [0.1, 0.1]
      .type = floats(size=2)
      .help = refinement sensitivity factor for spectrum coefficients
  }
  ranges {
    originZ = [-0.5, 0.5]
      .type = floats(size=2)
      .help = range of values in mm for originZ shift
    spectra0 = [-0.01, 0.01]
      .type = floats(size=2)
      .help = range of values for offset to Lambda correction (Angstrom)
    spectra1 = [0.95, 1.05]
      .type = floats(size=2)
      .help = range of values for multiplicative Lambda correction
    ucell_edge_percentage = None
      .type = float
      .help = jitter percentage for unit cell edge
    ucell_angle_deviation = None
      .type = float
      .help = absolute deviation for unit cell angle
    ncells_abc = None
      .type = ints
      .help = minimum and maximum values for the parameters Ncells_abc. Should be either 2 or 6 values
      .help = e.g. [3,300] such that 3 < Na,Nb,Nc < 300, or [3,100,3,100,3,40] such that  3 < Na,Nb < 100 and 3 < Nc < 40
      .help = From James Holton, a good rule of thumb is to require that N > 3.
    ncells_def = [0,200,0,200,0,200]
      .type = ints(size=6)
      .help = minimum and maximum values for the parameters Ncells_def
      .help = (see above for ncells_abc)
    panel_X = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in meters for XY shift
    panel_Y = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in meters for XY shift
    panel_Z = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in meters for detector panel Z shift
    panel_rotO = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in degrees for rotation of panel about its normal vector
    panel_rotF = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in degrees for rotation of panel about its fast-axis
    panel_rotS = [-1e-6, 1e-6]
      .type = floats(size=2)
      .help = range of values in degrees for rotation of panel about its slow-axis
    eta = [0,1]
      .type = floats(size=2)
      .help = range of possible values for mosaicity. First valus is min, second is max.
      .help = Note the min should always be greater than or equal to 0
  }
  update_oversample_during_refinement = False
    .type = bool
    .help = whether to update the oversample parameter as ncells changes
  compute_image_model_correlation = False
    .type = bool
    .help = whether to compute model image intensity correlations
  sigma_r = 3
    .type = float
    .help = standard deviation of the dark signal fluctuation
  adu_per_photon = 1
    .type = float
    .help = how many ADUs (detector units) equal 1 photon
  plot {
    display = False
      .type = bool
      .help = whether to make plots (default is two dimensional heat map plots)
    as_residuals = False
      .type = bool
      .help = whether to plot data minus model
    iteration_stride = 1
      .type = int
      .help = skip this many iterations before updating the plots
  }
  big_dump = False
    .type = bool
    .help = whether to output parameter information
  use_curvatures_threshold = 10
    .type = int
    .help = how many consecutiv positive curvature results before switching to curvature mode
  curvatures = False
    .type = bool
    .help = whether to try using curvatures
  start_with_curvatures = False
    .type = bool
    .help = whether to try using curvatures in the first iteration
  poissononly = False
    .type = bool
    .help = whether to only use poisson statistics (if detector dark signal is negligible)
  tradeps = 1e-2
    .type = float
    .help = LBFGS termination parameter  (smaller means minimize for longer)
  io {
    restart_file = None
      .type = str
      .help = output file for re-starting a simulation
    output_dir = None
      .type = str
      .help = optional output directory
  }
  quiet = False
    .type = bool
    .help = silence the refiner
  verbose = 0
    .type = int
    .help = verbosity level (0-10) for nanoBragg
  record_xy_calc = False
    .type = bool
    .help = whether to record xy calcs for Bragg spots modeled with diffBragg
  only_predict_model = False
    .type = bool
    .help = if True, just run diffBragg once to get the model predictions
  num_macro_cycles = 1
    .type = int
    .help = keep repeating the same refinement scheme over and over, this many times
  ncells_mask = *000 110 101 011 111
    .type = choice
    .help = a mask specifying which ncells parameters should be the same
    .help = e.g. 110 specifies Na and Nb are refined together as one parameter
  reference_geom = None
    .type = str
    .help = path to expt list file containing a detector model
  stage_two {
    merge_stat_freq = 3
      .type = int
      .help = show merge stats after this many iterations
    min_multiplicity = 1
      .type = int
      .help = structure factors whose multiplicity falls below this value
      .help = will not be refined
    print_reso_bins = False
      .type = bool
      .help = whether to show merge stats in resolution bins
    Fref_mtzname = None
      .type = str
      .help = path to a reference MTZ file. if passed, this is used solely to
      .help = observe the R-factor and CC between it and the Fobs being optimized
    Fref_mtzcol = "Famp(+),Famp(-)"
      .type = str
      .help = column in the mtz file containing the data
    d_min = 2
      .type = float
      .help = high res lim for binner
    d_max = 999
      .type = float
      .help = low res lim for binner
    n_bin = 10
      .type = int
      .help = number of binner bins
  }
}
"""

roi_phil = """
roi {
  cache_dir_only = False
    .type = bool
    .help = if True, create the cache folder , populate it with the roi data, then exit
  make_cache_dir = False
    .type = bool
    .help = if True, create a directory with name supplied by parameter cachefile_dir
    .help = and store roi data in that folder (experimental)
  cachefile_dir = None
    .type = str
    .help = path to a folder where roi data can be cached for faster re-loads
    .help = (experimental)
  fit_tilt = False
    .type = bool
    .help = fit tilt plane, or else background is simply an offset
  force_negative_background_to_zero = False
    .type = bool
    .help = if True and the background model evaluates to a negative number
    .help = within an ROI, then force the background to be 0 for all pixels in that ROI
  background_threshold = 3.5
    .type = float
    .help = for determining background pixels
  pad_shoebox_for_background_estimation = None
    .type = int
    .help = shoebox_size specifies the dimenstion of the shoebox used during refinement
    .help = and this parameter is used to increase that shoebox_size only during the background
    .help = estimation stage
  shoebox_size = 10
    .type = int
    .help = roi box dimension
  deltaQ = None
    .type = float
    .help = roi dimension in inverse Angstrom, such that shoeboxes at wider angles are larger.
    .help = If this parameter is supplied, shoebox_size will be ignored.
  reject_edge_reflections = True
    .type = bool
    .help = whether to reject ROIs if they occur near the detector panel edge
  reject_roi_with_hotpix = True
    .type = bool
    .help = whether to reject an ROI if it has a bad pixel
  background_mask = None
    .type = str
    .help = path to a mask specifying background (background pixels set to True)
  hotpixel_mask = None
    .type = str
    .help = path to a hotpixel mask (hot pixels set to True)
  panels = None
    .type = str
    .help = panel list for refinement as a string, e.g. "0-8,10,32-40" . The ranges are inclusive,
    .help = e.g. 0-8 means panels 0,1,2,3,4,5,6,7,8
  fit_tilt_using_weights = True
    .type = bool
    .help = if not using robust estimation for background, and instead using tilt plane fit,
    .help = then this parameter will toggle the use of weights. Weights are the estimated
    .help = pixel variance, incuding readout and shot noises.
}
"""

philz = simulator_phil + refiner_phil + roi_phil
phil_scope = parse(philz)
