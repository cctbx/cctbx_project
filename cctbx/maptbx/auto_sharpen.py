"""Command-line method to run auto-sharpening"""
from __future__ import absolute_import, division, print_function
import sys, os
import iotbx.phil
from libtbx.utils import Sorry
import libtbx.phil

map_modification_phil_str = """

     b_iso = None
       .type = float
       .short_caption = Target b_iso
       .help = Target B-value for map (sharpening will be applied to yield \
          this value of b_iso). If sharpening method is not supplied, \
          default is to use b_iso_to_d_cut sharpening.

     b_sharpen = None
       .type = float
       .short_caption = B-sharpen to apply
       .help = Sharpen with this b-value. Contrast with b_iso that yields a \
           targeted value of b_iso

     b_blur_hires = 200
       .type = float
       .short_caption = high_resolution blurring
       .help = Blur high_resolution data (higher than d_cut) with \
             this b-value. Contrast with b_sharpen applied to data up to\
             d_cut. \
             Note on defaults: If None and b_sharpen is positive (sharpening) \
             then high-resolution data is left as is (not sharpened). \
             If None and b_sharpen is negative (blurring) high-resolution data\
             is also blurred.

     resolution_dependent_b = None
       .type = floats
       .short_caption = resolution_dependent b

       .help = If set, apply resolution_dependent_b (b0 b1 b2). \
             Log10(amplitudes) will start at 1, change to b0 at half \
             of resolution specified, changing linearly, \
             change to b1/2 at resolution specified, \
             and change to b1/2+b2 at d_min_ratio*resolution

     normalize_amplitudes_in_resdep = False
       .type = bool
       .short_caption = Normalize amplitudes in resdep
       .help = Normalize amplitudes in resolution-dependent sharpening

     d_min_ratio = 0.833
       .type = float
       .short_caption = Sharpen d_min ratio
       .help = Sharpening will be applied using d_min equal to \
             d_min_ratio times resolution. Default is 0.833

     scale_max = 100000
       .type = float
       .short_caption = Scale_max
       .help = Scale amplitudes from inverse FFT to yield maximum of this value

     input_d_cut = None
       .type = float
       .short_caption = d_cut
       .help = High-resolution limit for sharpening

     rmsd = None
       .type = float
       .short_caption = RMSD of model
       .help = RMSD of model to true model (if supplied).  Used to \
             estimate expected fall-of with resolution of correct part \
             of model-based map. If None, assumed to be resolution \
             times rmsd_resolution_factor.

     rmsd_resolution_factor = 0.25
       .type = float
       .short_caption = rmsd resolution factor
        .help = default RMSD is resolution times resolution factor


     fraction_complete = None
       .type = float
       .short_caption = Completeness model
       .help = Completness of model (if supplied).  Used to \
             estimate correct part \
             of model-based map. If None, estimated from max(FSC).

     auto_sharpen = True
       .type = bool
       .short_caption = Automatically determine sharpening
       .help = Automatically determine sharpening using kurtosis maximization\
                 or adjusted surface area. Default is True

     auto_sharpen_methods = no_sharpening b_iso *b_iso_to_d_cut \
                            resolution_dependent model_sharpening \
                            half_map_sharpening target_b_iso_to_d_cut \
                            external_map_sharpening None

       .type = choice(multi=True)
       .short_caption = Sharpening methods
       .help = Methods to use in sharpening. b_iso searches for b_iso to \
          maximize sharpening target (kurtosis or adjusted_sa). \
          b_iso_to_d_cut applies b_iso only up to resolution specified, with \
          fall-over of k_sharpen.  Resolution dependent adjusts 3 parameters \
          to sharpen variably over resolution range. Default is \
          b_iso_to_d_cut . target_b_iso_to_d_cut uses target_b_iso_ratio \
          to set b_iso.

     box_in_auto_sharpen = False
       .type = bool
       .short_caption = Use box for auto_sharpening
       .help = Use a representative box of density for initial \
                auto-sharpening instead of the entire map. Default is False.

     density_select_in_auto_sharpen = True
       .type = bool
       .short_caption = density_select to choose box
       .help = Choose representative box of density for initial \
                auto-sharpening with density_select method \
                (choose region where there is high density).  Normally use\
                False for X-ray data and True for cryo-EM.

     density_select_threshold_in_auto_sharpen = None
       .type = float
       .short_caption = density_select threshold to choose box
       .help = Threshold for density select choice of box. Default is 0.05. \
               If your map has low overall contrast you might need to make this\
               bigger such as 0.2.

     allow_box_if_b_iso_set = False
       .type = bool
       .short_caption = Allow box if b_iso set
       .help = Allow box_in_auto_sharpen (if set to True) even if \
               b_iso is set. Default is to set box_n_auto_sharpen=False \
               if b_iso is set.

    soft_mask = True
      .type = bool
      .help = Use soft mask (smooth change from inside to outside with radius\
             based on resolution of map).
      .short_caption = Soft mask

     use_weak_density = False
       .type = bool
       .short_caption = Use box with poor density
       .help = When choosing box of representative density, use poor \
               density (to get optimized map for weaker density)

     discard_if_worse = None
       .type = bool
       .short_caption = Discard sharpening if worse
       .help = Discard sharpening if worse

     local_sharpening = None
       .type = bool
       .short_caption = Local sharpening
       .help = Sharpen locally using overlapping regions. \
               NOTE: Best to turn off local_aniso_in_local_sharpening \
               if NCS is present.\
               If local_aniso_in_local_sharpening is True and NCS is \
               present this can distort the map for some NCS copies \
               because an anisotropy correction is applied\
               based on local density in one copy and is transferred without \
               rotation to other copies.

     local_aniso_in_local_sharpening = None
       .type = bool
       .short_caption = Local anisotropy
       .help = Use local anisotropy in local sharpening.  \
               Default is True unless NCS is present.

     overall_before_local = True
       .type = bool
       .short_caption = Overall before local
       .help = Apply overall scaling before local scaling

     select_sharpened_map = None
       .type = int
       .short_caption = Sharpened map to use
       .help = Select a single sharpened map to use

     read_sharpened_maps = None
       .type = bool
       .short_caption = Read sharpened maps
       .help = Read in previously-calculated sharpened maps

     write_sharpened_maps = None
       .type = bool
       .short_caption = Write sharpened maps
       .help = Write out local sharpened maps

     smoothing_radius = None
       .type = float
       .short_caption = Smoothing radius
       .help = Sharpen locally using smoothing_radius. Default is 2/3 of \
                 mean distance between centers for sharpening

     box_center = None
       .type = floats
       .short_caption = Center of box
       .help = You can specify the center of the box (A units)

     box_size = 30 30 30
       .type = ints
       .short_caption = Size of box
       .help = You can specify the size of the box (grid units)

     target_n_overlap = 10
       .type = int
       .short_caption = Target overlap of boxes
       .help = You can specify the targeted overlap of boxes in local \
           sharpening

     restrict_map_size = None
       .type = bool
       .short_caption = Restrict box map size
       .help = Restrict box map to be inside full map (required for cryo-EM data).\
               Default is True if use_sg_symmetry=False and False \
                if use_sg_symmetry=True

     remove_aniso = True
       .type = bool
       .short_caption = Remove anisotropy
       .help = You can remove anisotropy (overall and locally) during sharpening

     max_box_fraction = 0.5
       .type = float
       .short_caption = Max size of box for auto_sharpening
       .help = If box is greater than this fraction of entire map, use \
                entire map. Default is 0.5.

     density_select_max_box_fraction = 0.95
       .type = float
       .short_caption = Max size of box for density_select
       .help = If box is greater than this fraction of entire map, use \
                entire map for density_select. Default is 0.95
    cc_cut = 0.2
       .type = float
       .short_caption = Min reliable CC in half-maps
       .help = Estimate of minimum highly reliable CC in half-map FSC. Used\
               to decide at what CC value to smooth the remaining CC values.

     max_cc_for_rescale = 0.2
       .type = float
       .short_caption = Max CC for rescaleMin reliable CC in half-maps
       .help = Used along with cc_cut and scale_using_last to correct for \
               small errors in FSC estimation at high resolution.  If the \
               value of FSC near the high-resolution limit is above \
               max_cc_for_rescale, assume these values are correct and do not \
               correct them.

     scale_using_last = 3
       .type = int
       .short_caption = Last N bins in FSC assumed to be about zero
       .help = If set, assume that the last scale_using_last bins in the FSC \
          for half-map or model sharpening are about zero (corrects for  \
          errors int the half-map process).


     mask_atoms = True
       .type = bool
       .short_caption = Mask atoms
       .help = Mask atoms when using model sharpening

     mask_atoms_atom_radius = 3
       .type =float
       .short_caption = Mask radius
       .help = Mask for mask_atoms will have mask_atoms_atom_radius

     value_outside_atoms = None
       .type = str
       .short_caption = Value outside atoms
       .help = Value of map outside atoms (set to 'mean' to have mean \
                value inside and outside mask be equal)

     k_sharpen = 10
       .type = float
       .short_caption = sharpening transition
       .help = Steepness of transition between sharpening (up to resolution \
           ) and not sharpening (d < resolution).  Note: for blurring, \
           all data are blurred (regardless of resolution), while for \
           sharpening, only data with d about resolution or lower are \
           sharpened. This prevents making very high-resolution data too \
           strong.  Note 2: if k_sharpen is zero or None, then no \
           transition is applied and all data is sharpened or blurred. \

     iterate = False
       .type = bool
       .short_caption = Iterate auto-sharpening
       .help = You can iterate auto-sharpening. This is useful in cases where \
                 you do not specify the solvent content and it is not \
                 accurately estimated until sharpening is optimized.

     optimize_b_blur_hires = False
       .type = bool
       .short_caption = Optimize value of b_blur_hires
       .help = Optimize value of b_blur_hires. \
                Only applies for auto_sharpen_methods b_iso_to_d_cut and \
                b_iso. This is normally carried out and helps prevent \
                over-blurring at high resolution if the same map is \
                sharpened more than once.

     optimize_d_cut = None
       .type = bool
       .short_caption = Optimize value of d_cut
       .help = Optimize value of d_cut. \
                Only applies for auto_sharpen_methods b_iso_to_d_cut and \
                b_iso

     adjust_region_weight = True
       .type = bool
       .short_caption = Adjust region weight
       .help = Adjust region_weight to make overall change in surface area \
               equal to overall change in normalized regions over the range \
               of search_b_min to search_b_max using b_iso_to_d_cut.

     region_weight_method = initial_ratio *delta_ratio b_iso
       .type = choice
       .short_caption = Region weight method
       .help = Method for choosing region_weights. Initial_ratio uses \
               ratio of surface area to regions at low B value.  Delta \
               ratio uses change in this ratio from low to high B. B_iso \
               uses resolution-dependent b_iso (not weights) with the \
               formula b_iso=5.9*d_min**2

     region_weight_factor = 1.0
       .type = float
       .short_caption = Region weight factor
       .help = Multiplies region_weight after calculation with \
               region_weight_method above

     region_weight_buffer = 0.1
       .type = float
       .short_caption = Region weight factor buffer
       .help = Region_weight adjusted to be region_weight_buffer \
               away from minimum or maximum values

     target_b_iso_ratio = 5.9
       .type = float
       .short_caption = Target b_iso ratio
       .help = Target b_iso ratio : b_iso is estimated as \
               target_b_iso_ratio * resolution**2

     target_b_iso_model_scale = 0.
       .type = float
       .short_caption = scale on target b_iso ratio for model
       .help = For model sharpening, the target_biso is scaled \
                (normally zero).

     signal_min = 3.0
       .type = float
       .short_caption = Minimum signal
       .help = Minimum signal in estimation of optimal b_iso.  If\
                not achieved, use any other method chosen.

     search_b_min = None
       .type = float
       .short_caption = Low bound for b_iso search
       .help = Low bound for b_iso search. Default is -100.

     search_b_max = None
       .type = float
       .short_caption = High bound for b_iso search
       .help = High bound for b_iso search. Default is 300.

     search_b_n = None
       .type = int
       .short_caption = Number of b_iso values to search
       .help = Number of b_iso values to search. Default is 21.

     residual_target = None
       .type = str
       .short_caption = Residual target
       .help = Target for maximization steps in sharpening.  \
          Can be kurtosis or adjusted_sa (adjusted surface area).\
          Default is adjusted_sa.

     sharpening_target = None
       .type = str
       .short_caption = Overall sharpening target
       .help = Overall target for sharpening.  Can be kurtosis or adjusted_sa \
          (adjusted surface area).  Used to decide which sharpening approach \
          is used. Note that during optimization, residual_target is used \
          (they can be the same.) Default is adjusted_sa.

     require_improvement = None
       .type = bool
       .short_caption = Require improvement
       .help = Require improvement in score for sharpening to be applied.\
                Default is True.

     region_weight = None
       .type = float
       .short_caption = Region weighting
       .help = Region weighting in adjusted surface area calculation.\
            Score is surface area minus region_weight times number of regions.\
            Default is set automatically.  \
            A smaller value will give more sharpening.

     sa_percent = None
       .type = float
       .short_caption = Percent of target regions in adjusted_sa
       .help = Percent of target regions used in calulation of adjusted \
         surface area.  Default is 30.

     fraction_occupied = None
       .type = float
       .short_caption = Fraction of molecular volume inside contours
       .help = Fraction of molecular volume targeted to be inside contours. \
           Used to set contour level. Default is 0.20

      n_bins = None
        .type = int
        .short_caption = Resolution bins
        .help = Number of resolution bins for sharpening. Default is 20.
        .expert_level = 1

     regions_to_keep = None
       .type = int
       .short_caption = Regions to keep
       .help = You can specify a limit to the number of regions to keep\
                when generating the asymmetric unit of density.

      max_regions_to_test = None
        .type = int
        .short_caption = Max regions to test
        .help = Number of regions to test for surface area in adjusted_sa \
                scoring of sharpening. Default is 30

      eps = None
        .type = float
        .short_caption = Shift used in calculation of derivatives for \
           sharpening maximization.  Default is 0.01 for kurtosis and 0.5 for \
           adjusted_sa.

      k_sol = 0.35
        .type = float
        .help = k_sol value for model map calculation. IGNORED (Not applied)
        .short_caption = k_sol IGNORED
        .style = hidden

      b_sol = 50
        .type = float
        .help = b_sol value for model map calculation. IGNORED (Not applied)
        .short_caption = b_sol IGNORED
        .style = hidden
    """

output_files = """

    shifted_map_file = shifted_map.ccp4
      .type = str
      .help = Input map file shifted to new origin.
      .short_caption = Shifted map file
      .style = new_file

    sharpened_map_file = None
      .type = str
      .help = Sharpened input map file. In the same location as input map.
      .short_caption = Sharpened map file
      .input_size = 400
      .style = new_file
      .expert_level = 1

    shifted_sharpened_map_file = None
      .type = str
      .help = Input map file shifted to place origin at 0,0,0 and sharpened.
      .short_caption = Shifted sharpened map file
      .input_size = 400
      .style = new_file

    sharpened_map_coeffs_file = None
      .type = str
      .help = Sharpened input map \
              (shifted to new origin if original origin was not 0,0,0), \
              written out as map coefficients
      .short_caption = Sharpened map coeffs file
      .input_size = 400
      .style = new_file

    output_weight_map_pickle_file = weight_map_pickle_file.pkl
       .type = path
       .short_caption = Output weight map pickle file
       .help = Output weight map pickle file
       .style = new_file

    output_directory =  None
      .type = path
      .help = Directory where output files are to be written \
                applied.
      .short_caption = Output directory
      .style = new_file directory
"""

crystal_info = """
     is_crystal = None
       .type = bool
       .short_caption = Is a crystal
       .help = Defines whether this is a crystal (or cryo-EM).\
                Default is True if use_sg_symmetry=True and False otherwise.

     resolution = None
       .type = float
       .short_caption = Resolution
       .help = Optional nominal resolution of the map.
       .style = resolution
       .expert_level = 1


     solvent_content = None
       .type = float
       .help = Optional solvent fraction of the cell.
       .short_caption = Solvent content

     solvent_content_iterations = 3
       .type = int
       .help = Iterations of solvent fraction estimation. Used for ID of \
               solvent content in boxed maps.
       .short_caption = Solvent fraction iterations
       .style = hidden

     molecular_mass = None
       .type = float
       .help = Molecular mass of molecule in Da. Used as alternative method \
                 of specifying solvent content.
       .short_caption = Molecular mass in Da

      ncs_copies = None
        .type = int
        .help = You can specify ncs copies and seq file to define solvent \
            content
        .short_caption = NCS copies

     wang_radius = None
       .type = float
       .help = Wang radius for solvent identification. \
           Default is 1.5* resolution
       .short_caption = Wang radius

     buffer_radius = None
       .type = float
       .help = Buffer radius for mask smoothing. \
           Default is resolution
       .short_caption = Buffer radius

     pseudo_likelihood = None
       .type = bool
       .help = Use pseudo-likelihood method for half-map sharpening. \
               (In development)
       .short_caption = Pseudo-likelihood
       .style = hidden
"""

master_phil = iotbx.phil.parse("""

  input_files
    .style = menu_item auto_align
  {

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file
      .style = file_type:ccp4_map bold input_file

    half_map_file = None
      .type = path
      .multiple = True
      .short_caption = Half map
      .style = file_type:ccp4_map bold input_file
      .help = Half map (two should be supplied) for FSC calculation. Must \
               have grid identical to map_file

    external_map_file = None
      .type = path
      .short_caption = External map
      .style = file_type:ccp4_map bold input_file
      .help = External map to be used to scale map_file (power vs resolution\
              will be matched)

    map_coeffs_file = None
      .type = path
      .help = Optional file with map coefficients
      .short_caption = Map coefficients
      .style = bold file_type:hkl input_file process_hkl \
        child:map_labels:map_coeffs_labels \
        child:space_group:space_group child:unit_cell:unit_cell

    map_coeffs_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of of map coefficients \
          to use
      .short_caption = Map coeffs label
      .style = renderer:draw_map_arrays_widget

    pdb_file = None
      .type = path
      .help = If a model is supplied, the map will be adjusted to \
                maximize map-model correlation.  This can be used \
                to improve a map in regions where no model is yet \
                built.
      .style = file_type:pdb input_file
      .short_caption = Model file (optional)

    ncs_file = None
      .type = path
      .help = File with NCS information (typically point-group NCS with \
               the center specified). Typically in  PDB format. \
              Can also be a .ncs_spec file from phenix. \
              Created automatically if symmetry is specified.
      .short_caption = NCS info file

    seq_file = None
       .type = path
       .short_caption = Sequence file
       .style = file_type:seq input_file
       .help = Sequence file (unique chains only,  \
               1-letter code, chains separated by \
               blank line or greater-than sign.)  \
               Can have chains that are DNA/RNA/protein and\
               all can be present in one file.

    input_weight_map_pickle_file = None
      .type = path
      .short_caption = Input weight map pickle file
      .help = Weight map pickle file

  }

  output_files
    .style = menu_item auto_align
  {
   %s
  }

  crystal_info
    .style = menu_item auto_align
  {
   %s
  }

  map_modification {
    %s
  }


   control
     .style = menu_item auto_align
   {
     verbose = False
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output

     resolve_size = None
        .type = int
        .help = "Size of resolve to use. "
        .style = hidden

     ignore_map_limitations = None
       .type = bool
       .short_caption = Ignore map limitations
       .help = Ignore limitations such as 'map cannot be sharpened'

      multiprocessing = *multiprocessing sge lsf pbs condor pbspro slurm
        .type = choice
        .short_caption = multiprocessing type
        .help = Choices are multiprocessing (single machine) or queuing systems

      queue_run_command = None
        .type = str
        .short_caption = Queue run command
        .help = run command for queue jobs. For example qsub.

      nproc = 1
        .type = int
        .short_caption = Number of processors
        .help = Number of processors to use
        .style = renderer:draw_nproc_widget bold

   }

  include scope libtbx.phil.interface.tracking_params

  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }

""" %(output_files, crystal_info, map_modification_phil_str),
  process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):
  """Get parameters from args"""

  command_line = iotbx.phil.process_command_line_with_files(
    reflection_file_def="input_files.map_coeffs_file",
    map_file_def="input_files.map_file",
    pdb_file_def="input_files.pdb_file",
    seq_file_def="input_files.seq_file",
    ncs_file_def="input_files.ncs_file",
    args=args,
    master_phil=master_phil)


  print("\nAuto-sharpen a map\n", file=out)
  params=command_line.work.extract()
  print("Command used: %s\n" %(
   " ".join(['phenix.auto_sharpen']+args)), file=out)
  master_params.format(python_object=params).show(out=out)

  if params.output_files.output_directory is None:
    params.output_files.output_directory=os.getcwd()
  elif not os.path.isdir(params.output_files.output_directory):
    os.mkdir(params.output_files.output_directory)

  params=set_sharpen_params(params,out)
  return params


def set_sharpen_params(params,out=sys.stdout):
  """Set and check sharpening parameters"""

  if params.map_modification.resolution_dependent_b==[0,0,0]:
    params.map_modification.resolution_dependent_b=[0,0,1.e-10]
    # just so that we know it was set

  if params.map_modification.b_iso:
    if params.map_modification.k_sharpen and \
        'b_iso_to_d_cut' in params.map_modification.auto_sharpen_methods:
      params.map_modification.auto_sharpen_methods=['b_iso_to_d_cut']
    elif 'b_iso_to_d_cut' in params.map_modification.auto_sharpen_methods and\
        len(params.map_modification.auto_sharpen_methods)==1:
      params.map_modification.auto_sharpen_methods=['b_iso_to_d_cut']
    else:
      params.map_modification.auto_sharpen_methods=['b_iso']

  if (params.map_modification.b_iso and \
       params.map_modification.auto_sharpen_methods in [
       ['b_iso_to_d_cut'],['b_iso']])   or  \
     (params.map_modification.resolution_dependent_b and \
       params.map_modification.auto_sharpen_methods in [
       ['resolution_dependent']]):
    if params.map_modification.box_in_auto_sharpen and (
      not getattr(params.map_modification,'allow_box_if_b_iso_set',None)):
      params.map_modification.box_in_auto_sharpen=False
      print("Set box_in_auto_sharpen=False as parameters are set "+\
        "\nand sharpening method is %s" %(
        params.map_modification.auto_sharpen_methods[0]), file=out)
    if params.map_modification.density_select_in_auto_sharpen and (
      not getattr(params.map_modification,'allow_box_if_b_iso_set',None)):
      params.map_modification.density_select_in_auto_sharpen=False
      print("Set density_select_in_auto_sharpen=False as parameters are set "+\
        "\nand sharpening method is %s" %(
        params.map_modification.auto_sharpen_methods[0]), file=out)

  if params.map_modification.optimize_b_blur_hires and \
    not 'b_iso_to_d_cut' in params.map_modification.auto_sharpen_methods and \
       not 'b_iso' in params.map_modification.auto_sharpen_methods:
     print("Set optimize_b_blur_hires=False as neither b_iso_to_d_cut nor"+\
         " b_iso are used", file=out)
     params.map_modification.optimize_b_blur_hires=False
  if params.map_modification.iterate and \
    not 'b_iso_to_d_cut' in params.map_modification.auto_sharpen_methods and \
       not 'b_iso' in params.map_modification.auto_sharpen_methods:
     print("Set iterate=False as neither b_iso_to_d_cut nor"+\
         " b_iso are used", file=out)
     params.map_modification.iterate=False

  if 'half_map_sharpening' in params.map_modification.auto_sharpen_methods and \
   ((not params.input_files.half_map_file) or
   (not len(params.input_files.half_map_file))==2):
    print("Skipping half-map sharpening as there are no half-maps supplied...",
      file = out)
    params.map_modification.auto_sharpen_methods.remove('half_map_sharpening')

  return params


def get_map_coeffs_from_file(
      map_coeffs_file=None,
      map_coeffs_labels=None):
    """Get map coeffs from a file with reflection_file_reader. Use
     data_manager instead."""
    from iotbx import reflection_file_reader
    reflection_file=reflection_file_reader.any_reflection_file(
        map_coeffs_file)
    mtz_content=reflection_file.file_content()
    for ma in reflection_file.as_miller_arrays(merge_equivalents=True):
      if not ma.is_complex_array(): continue
      labels=",".join(ma.info().labels)
      if not map_coeffs_labels or labels==map_coeffs_labels:  # take it
         return ma

def map_inside_cell(pdb_hierarchy,crystal_symmetry=None):
  """Map (move) coordinates in hierarchy to center them in the cell"""
  pa=pdb_hierarchy.atoms()
  sites_cart=pa.extract_xyz()
  from cctbx.maptbx.segment_and_split_map import move_xyz_inside_cell
  new_sites_cart=move_xyz_inside_cell(xyz_cart=sites_cart,
     crystal_symmetry=crystal_symmetry)
  pa.set_xyz(new_sites_cart)
  return pdb_hierarchy


def get_map_and_model(params=None,
    map_data=None,
    crystal_symmetry=None,
    pdb_hierarchy=None,
    ncs_obj=None,
    half_map_data_list=None,
    map_coords_inside_cell=True,
    get_map_labels=None,
    out=sys.stdout):

  """Get maps and model for auto_sharpen.
  Returns   pdb_hierarchy,
       map_data,half_map_data_list,ncs_obj,crystal_symmetry,acc,
       original_crystal_symmetry,original_unit_cell_grid.
  if get_map_labels is set, also return map_labels
  """


  acc=None # accessor used to shift map back to original location if desired
  origin_frac=(0,0,0)
  map_labels=None
  read_map_coeffs = False
  if map_data and crystal_symmetry:
    original_crystal_symmetry=crystal_symmetry
    original_unit_cell_grid=None
    acc=map_data.accessor()
    shift_needed = not \
       (map_data.focus_size_1d() > 0 and map_data.nd() == 3 and
        map_data.is_0_based())
    if(shift_needed):
      origin_shift=(
        map_data.origin()[0]/map_data.all()[0],
        map_data.origin()[1]/map_data.all()[1],
        map_data.origin()[2]/map_data.all()[2])
      origin_frac=origin_shift  # NOTE: fraction of NEW cell
      map_data = map_data.shift_origin()
    else:
      origin_frac=(0.,0.,0.)

  elif params.input_files.map_file:
    print("\nReading map from %s\n" %( params.input_files.map_file), file=out)
    from cctbx.maptbx.segment_and_split_map import get_map_object
    map_data,space_group,unit_cell,crystal_symmetry,origin_frac,acc,\
        original_crystal_symmetry,original_unit_cell_grid,map_labels=\
      get_map_object(file_name=params.input_files.map_file,
      must_allow_sharpening=(not params.control.ignore_map_limitations),
      get_map_labels=True,out=out)
    map_data=map_data.as_double()
    if origin_frac != (0,0,0) and acc is None:
      print("\nWARNING: Unable to place output map at position of "+\
        "input map though input map has non-zero origin at %s\n" %(
        str(origin_frac)), file=out)

  elif params.input_files.map_coeffs_file:
    map_coeffs=get_map_coeffs_from_file(
      map_coeffs_file=params.input_files.map_coeffs_file,
      map_coeffs_labels=params.input_files.map_coeffs_labels)

    if not map_coeffs:
      raise Sorry("Could not get map coeffs from %s with labels %s" %(
        params.input_files.map_coeffs_file,params.input_files.map_coeffs_labels))
    read_map_coeffs = True
    print("Map coefficients read from %s with labels %s" %(
         params.input_files.map_coeffs_file,
         str(params.input_files.map_coeffs_labels)), file=out)
    crystal_symmetry=map_coeffs.crystal_symmetry()
    from cctbx.development.create_models_or_maps import get_map_from_map_coeffs
    map_data=get_map_from_map_coeffs(
      map_coeffs=map_coeffs,crystal_symmetry=crystal_symmetry)
    acc=map_data.accessor()
    original_crystal_symmetry=crystal_symmetry
    original_unit_cell_grid=None
    if not params.crystal_info.resolution:
      params.crystal_info.resolution=map_coeffs.d_min()
      print("Resolution from map_coeffs is %7.2f A" %(
          params.crystal_info.resolution), file=out)

  else:
    raise Sorry("Need ccp4 map or map_coeffs")
  if params.crystal_info.is_crystal is None:
    if read_map_coeffs:
      print("Assuming this is crystallographic data as map coefficients read",
         file = out)
      params.crystal_info.is_crystal = True
    elif (not (crystal_symmetry.space_group().type().number() in [0, 1])):
      print("Assuming this is crystallographic data as map is not P1",
         file = out)
      params.crystal_info.is_crystal = True
    else:
      print("Assuming this is cryo-EM data as map read in space group P1",
         file = out)
      params.crystal_info.is_crystal = False
  elif params.crystal_info.is_crystal:
      print("Assuming this is crystallographic data as is_crystal=True",
         file = out)
  else:
      print("Assuming this is cryo-EM data as is_crystal=False",
         file = out)


  if params.input_files.external_map_file and not \
      params.map_modification.auto_sharpen_methods==['external_map_sharpening']:
    raise Sorry("For external map sharpening, please do not set any other "+
        "options in 'Sharpening methods' (options set: %s)" %(
      ", ".join(params.map_modification.auto_sharpen_methods)))


  if params.map_modification.auto_sharpen_methods==['external_map_sharpening']:
    if not params.input_files.external_map_file:
      raise Sorry("Need external_map_file for external_map")
    if params.input_files.half_map_file:
      raise Sorry("Cannot use half_map_file with external_map")
    half_map_data_list=[]
    file_name=params.input_files.external_map_file
    print("\nReading external map from %s\n" %(file_name), file=out)
    half_map_data,half_map_space_group,half_map_unit_cell,\
        half_map_crystal_symmetry,half_map_origin_frac,half_map_acc,\
         half_map_original_crystal_symmetry,half_map_original_unit_cell_grid=\
        get_map_object(file_name=file_name,out=out)
    half_map_data=half_map_data.as_double()
    assert half_map_crystal_symmetry.is_similar_symmetry(crystal_symmetry)
    half_map_data_list.append(half_map_data) # save it in half_map_data_list
    # and add another just so we have two..
    half_map_data_list.append(half_map_data.deep_copy())
    if map_data and half_map_data.size()!=map_data.size():
      raise Sorry("Map data and external_map_file data must be the same size")


  if params.input_files.half_map_file:
    if len(params.input_files.half_map_file) != 2:
      raise Sorry("Please supply zero or two half_map files")
    half_map_data_list=[]
    from cctbx.maptbx.segment_and_split_map import get_map_object
    for file_name in params.input_files.half_map_file:
      print("\nReading half-map from %s\n" %(file_name), file=out)
      half_map_data,half_map_space_group,half_map_unit_cell,\
        half_map_crystal_symmetry,half_map_origin_frac,half_map_acc,\
         half_map_original_crystal_symmetry,half_map_original_unit_cell_grid=\
        get_map_object(file_name=file_name,out=out)
      half_map_data=half_map_data.as_double()
      assert half_map_crystal_symmetry.is_similar_symmetry(crystal_symmetry)

      half_map_data_list.append(half_map_data)

  if params.crystal_info.resolution is None:
    raise Sorry("Need resolution for b-sharpening if "+
       "map is supplied and no model is supplied")

  if params.crystal_info.resolution >= 10:
    print("\n** WARNING: auto_sharpen is designed for maps at a "+\
      "resolution of about 4.5 A\nor better.  Sharpening may be"+\
      "poor at %7.0f A" %(params.crystal_info.resolution), file=out)


  if params.input_files.pdb_file and not pdb_hierarchy: # get model
    model_file=params.input_files.pdb_file
    if not os.path.isfile(model_file):
      raise Sorry("Missing the model file: %s" %(model_file))
    from iotbx.pdb.utils import get_pdb_hierarchy
    pdb_hierarchy= get_pdb_hierarchy(file_name = model_file)
  if pdb_hierarchy: # XXX added 2019-05-05
    if origin_frac != (0,0,0):
      print("Shifting model by %s" %(str(origin_frac)), file=out)
      from cctbx.maptbx.segment_and_split_map import \
         apply_shift_to_pdb_hierarchy
      origin_shift=crystal_symmetry.unit_cell().orthogonalize(
         (-origin_frac[0],-origin_frac[1],-origin_frac[2]))
      pdb_hierarchy=apply_shift_to_pdb_hierarchy(
       origin_shift=origin_shift,
       crystal_symmetry=crystal_symmetry,
       pdb_hierarchy=pdb_hierarchy,
       out=out)
    if map_coords_inside_cell:
      # put inside (0,1)
      pdb_hierarchy=map_inside_cell(pdb_hierarchy,crystal_symmetry=crystal_symmetry)

  if params.input_files.ncs_file and not ncs_obj: # NCS
    from cctbx.maptbx.segment_and_split_map import get_ncs
    ncs_obj,dummy_tracking_data=get_ncs(params,out=out)
    if origin_frac != (0,0,0):
      origin_shift=crystal_symmetry.unit_cell().orthogonalize(
         (-origin_frac[0],-origin_frac[1],-origin_frac[2]))
      print("Shifting NCS by (%7.2f,%7.2f,%7.2f) " %((origin_shift)), file=out)
      from scitbx.math import  matrix
      ncs_obj=ncs_obj.coordinate_offset(
       coordinate_offset=matrix.col(origin_shift))

  if get_map_labels:
    return pdb_hierarchy,\
       map_data,half_map_data_list,ncs_obj,crystal_symmetry,acc,\
       original_crystal_symmetry,original_unit_cell_grid,map_labels
  else:
    return pdb_hierarchy,\
      map_data,half_map_data_list,ncs_obj,crystal_symmetry,acc,\
       original_crystal_symmetry,original_unit_cell_grid


def run(args=None,params=None,
    map_data=None,crystal_symmetry=None,
    wrapping = None,
    write_output_files=True,
    pdb_hierarchy=None,
    ncs_obj=None,
    return_map_data_only=False,
    return_unshifted_map=False,
    half_map_data_list=None,
    ncs_copies=None,
    n_residues=None,
    out=sys.stdout):
  """Run auto-sharpening"""

  # Get the parameters
  if not params:
    params=get_params(args,out=out)

  if not ncs_copies:
    ncs_copies=params.crystal_info.ncs_copies

  # get map_data and crystal_symmetry




  pdb_hierarchy,map_data,half_map_data_list,ncs_obj,crystal_symmetry,acc,\
       original_crystal_symmetry,original_unit_cell_grid,map_labels=\
        get_map_and_model(
     map_data=map_data,
     half_map_data_list=half_map_data_list,
     pdb_hierarchy=pdb_hierarchy,
     ncs_obj=ncs_obj,
     map_coords_inside_cell=False,
     crystal_symmetry=crystal_symmetry,
     get_map_labels=True,
     params=params,out=out)
  # NOTE: map_data is now relative to origin at (0,0,0).
  # Use map_data.reshape(acc) to put it back where it was if acc is not None

  # auto-sharpen the map
  from cctbx.maptbx.segment_and_split_map import auto_sharpen_map_or_map_coeffs
  si=auto_sharpen_map_or_map_coeffs(
        resolution=params.crystal_info.resolution, # required
        crystal_symmetry=crystal_symmetry,
        is_crystal=params.crystal_info.is_crystal,
        verbose=params.control.verbose,
        resolve_size=params.control.resolve_size,
        multiprocessing=params.control.multiprocessing,
        nproc=params.control.nproc,
        queue_run_command=params.control.queue_run_command,
        map=map_data,
        wrapping=wrapping,
        half_map_data_list=half_map_data_list,
        solvent_content=params.crystal_info.solvent_content,
        molecular_mass=params.crystal_info.molecular_mass,
        input_weight_map_pickle_file=\
            params.input_files.input_weight_map_pickle_file,
        output_weight_map_pickle_file=\
            params.output_files.output_weight_map_pickle_file,
        read_sharpened_maps=params.map_modification.read_sharpened_maps,
        write_sharpened_maps=params.map_modification.write_sharpened_maps,
        select_sharpened_map=params.map_modification.select_sharpened_map,
        auto_sharpen=params.map_modification.auto_sharpen,
        local_sharpening=params.map_modification.local_sharpening,
        output_directory=params.output_files.output_directory,
        smoothing_radius=params.map_modification.smoothing_radius,
        local_aniso_in_local_sharpening=\
           params.map_modification.local_aniso_in_local_sharpening,
        overall_before_local=\
           params.map_modification.overall_before_local,
        box_in_auto_sharpen=params.map_modification.box_in_auto_sharpen,
        density_select_in_auto_sharpen=params.map_modification.density_select_in_auto_sharpen,
        density_select_threshold_in_auto_sharpen=params.map_modification.density_select_threshold_in_auto_sharpen,
        use_weak_density=params.map_modification.use_weak_density,
        discard_if_worse=params.map_modification.discard_if_worse,
        box_center=params.map_modification.box_center,
        box_size=params.map_modification.box_size,
        target_n_overlap=params.map_modification.target_n_overlap,
        restrict_map_size=params.map_modification.restrict_map_size,
        remove_aniso=params.map_modification.remove_aniso,
        auto_sharpen_methods=params.map_modification.auto_sharpen_methods,
        residual_target=params.map_modification.residual_target,
        region_weight=params.map_modification.region_weight,
        sa_percent=params.map_modification.sa_percent,
        eps=params.map_modification.eps,
        n_bins=params.map_modification.n_bins,
        max_regions_to_test=params.map_modification.max_regions_to_test,
        regions_to_keep=params.map_modification.regions_to_keep,
        fraction_occupied=params.map_modification.fraction_occupied,
        sharpening_target=params.map_modification.sharpening_target,
        d_min_ratio=params.map_modification.d_min_ratio,
        scale_max=params.map_modification.scale_max,
        input_d_cut=params.map_modification.input_d_cut,
        b_blur_hires=params.map_modification.b_blur_hires,
        max_box_fraction=params.map_modification.max_box_fraction,
        cc_cut=params.map_modification.cc_cut,
        max_cc_for_rescale=params.map_modification.max_cc_for_rescale,
        scale_using_last=params.map_modification.scale_using_last,
        density_select_max_box_fraction=params.map_modification.density_select_max_box_fraction,
        mask_atoms=params.map_modification.mask_atoms,
        mask_atoms_atom_radius=params.map_modification.mask_atoms_atom_radius,
        value_outside_atoms=params.map_modification.value_outside_atoms,
        k_sharpen=params.map_modification.k_sharpen,
        optimize_b_blur_hires=params.map_modification.optimize_b_blur_hires,
        iterate=params.map_modification.iterate,
        optimize_d_cut=params.map_modification.optimize_d_cut,
        soft_mask=params.map_modification.soft_mask,
        allow_box_if_b_iso_set=params.map_modification.allow_box_if_b_iso_set,
        search_b_min=params.map_modification.search_b_min,
        search_b_max=params.map_modification.search_b_max,
        search_b_n=params.map_modification.search_b_n,
        adjust_region_weight=params.map_modification.adjust_region_weight,
        region_weight_method=params.map_modification.region_weight_method,
        region_weight_factor=params.map_modification.region_weight_factor,
        region_weight_buffer=\
            params.map_modification.region_weight_buffer,
        target_b_iso_ratio=params.map_modification.target_b_iso_ratio,
        signal_min=params.map_modification.signal_min,
        buffer_radius=params.crystal_info.buffer_radius,
        wang_radius=params.crystal_info.wang_radius,
        pseudo_likelihood=params.crystal_info.pseudo_likelihood,
        target_b_iso_model_scale=params.map_modification.target_b_iso_model_scale,
        b_iso=params.map_modification.b_iso,
        b_sharpen=params.map_modification.b_sharpen,
        resolution_dependent_b=\
           params.map_modification.resolution_dependent_b,
        normalize_amplitudes_in_resdep=\
           params.map_modification.normalize_amplitudes_in_resdep,
        pdb_hierarchy=pdb_hierarchy,
        ncs_obj=ncs_obj,
        rmsd=params.map_modification.rmsd,
        rmsd_resolution_factor=params.map_modification.rmsd_resolution_factor,
        b_sol=params.map_modification.b_sol,
        k_sol=params.map_modification.k_sol,
        fraction_complete=params.map_modification.fraction_complete,
        seq_file=params.input_files.seq_file,
        ncs_copies=ncs_copies,
        n_residues=n_residues,
        out=out)

  # get map_data and map_coeffs of final map

  new_map_data=si.as_map_data()
  new_map_coeffs=si.as_map_coeffs()

  from cctbx.maptbx.segment_and_split_map import get_b_iso,map_coeffs_as_fp_phi
  f,phi=map_coeffs_as_fp_phi(new_map_coeffs)
  temp_b_iso=get_b_iso(f,d_min=params.crystal_info.resolution)

  if not si.is_model_sharpening():
    print(file=out)
    print(80*"=","\n",80*"=", file=out)
    print("\n           Summary of sharpening information\n ", file=out)
    si.show_summary(verbose=params.control.verbose,out=out)
    print(80*"=","\n",80*"=", file=out)

  # write out the new map_coeffs and map if requested:

  offset_map_data=new_map_data.deep_copy()
  if acc is not None:  # offset the map to match original if possible
    offset_map_data.reshape(acc)

  if write_output_files and (not params.output_files.sharpened_map_file):
    params.output_files.sharpened_map_file = 'sharpened_map.ccp4'
  if write_output_files and params.output_files.sharpened_map_file and \
      offset_map_data:
    output_map_file=os.path.join(params.output_files.output_directory,
        params.output_files.sharpened_map_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    if acc is not None:  # we offset the map to match original
      print("\nWrote sharpened map in original location with origin at %s\nto %s" %(
         str(offset_map_data.origin()),output_map_file), file=out)
      write_ccp4_map(original_crystal_symmetry, output_map_file,
        offset_map_data,
        output_unit_cell_grid=original_unit_cell_grid)
    else:
      print("\nWrote sharpened map with origin at 0,0,0 "+\
        "(NOTE: may be boxed map and may not be "+\
        "\nsame as original location) to %s\n" %(
         output_map_file), file=out)

    from iotbx.mrcfile import create_output_labels
    program_name='auto_sharpen'
    limitations=["map_is_sharpened"]
    labels=create_output_labels(program_name=program_name,
       input_file_name=params.input_files.map_file,
       input_labels=map_labels,
       limitations=limitations,
       output_labels=None)

    write_ccp4_map(crystal_symmetry, output_map_file, offset_map_data,
        labels=labels)

  if write_output_files and params.output_files.shifted_sharpened_map_file:
    output_map_file=os.path.join(params.output_files.output_directory,
        params.output_files.shifted_sharpened_map_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    write_ccp4_map(crystal_symmetry, output_map_file, new_map_data)
    print("\nWrote sharpened map (origin at %s)\nto %s" %(
     str(new_map_data.origin()),output_map_file), file=out)

  if write_output_files and (not params.output_files.sharpened_map_coeffs_file):
     params.output_files.sharpened_map_coeffs_file = 'sharpened_map_coeffs.mtz'
  if write_output_files and params.output_files.sharpened_map_coeffs_file and \
      new_map_coeffs:
    output_map_coeffs_file=os.path.join(params.output_files.output_directory,
        params.output_files.sharpened_map_coeffs_file)
    new_map_coeffs.as_mtz_dataset(column_root_label='FWT').mtz_object().write(
       file_name=output_map_coeffs_file)
    print("\nWrote sharpened map_coeffs (origin at 0,0,0)\n to %s\n" %(
       output_map_coeffs_file), file=out)

  if return_unshifted_map:
    map_to_return=offset_map_data
  else:
    map_to_return=new_map_data

  if return_map_data_only:
    return map_to_return
  else:  #usual
    return map_to_return,new_map_coeffs,crystal_symmetry,si

# =============================================================================
# GUI-specific bits for running command
from libtbx import runtime_utils
class launcher(runtime_utils.target_with_save_result):
  """Launch this method from GUI"""
  def run(self):
    import os
    from wxGUI2 import utils
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    map_to_return, new_map_coeffs, crystal_symmetry, si = \
      run(args=self.args, out=sys.stdout)
    si.map_data = None
    result = (None, None, crystal_symmetry, si)
    return result

def validate_params(params):
  """Check parameters for GUI"""
  if ( (params.input_files.map_coeffs_file is None) and
       (params.input_files.map_file is None) ):
    raise Sorry('Please provide a map file.')
  if ( (params.input_files.map_coeffs_file is not None) and
       (params.input_files.map_coeffs_labels is None) ):
    raise Sorry('Please select the label for the map coefficients.')
  if ( (params.input_files.map_file is not None) and
       (params.crystal_info.resolution is None) ):
    raise Sorry('Please provide a resolution limit.')
  return True

# =============================================================================

if __name__=="__main__":
  run(sys.argv[1:])
