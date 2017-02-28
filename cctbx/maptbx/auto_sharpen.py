from __future__ import division
import sys, os
import iotbx.phil
from libtbx.utils import Sorry

master_phil = iotbx.phil.parse("""

  input_files
    .style = menu_item auto_align
  {

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

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
      .short_caption = Model file
  }

  output_files
    .style = menu_item auto_align
  {

    shifted_map_file = shifted_map.ccp4
      .type = str
      .help = Input map file shifted to new origin.
      .short_caption = Shifted map file

    sharpened_map_file = sharpened_map.ccp4
      .type = str
      .help = Sharpened input map file. In the same location as input map.
      .short_caption = Sharpened map file
      .input_size = 400

    shifted_sharpened_map_file = None
      .type = str
      .help = Input map file shifted to place origin at 0,0,0 and sharpened.
      .short_caption = Shifted sharpened map file
      .input_size = 400

    sharpened_map_coeffs_file = sharpened_map_coeffs.mtz
      .type = str
      .help = Sharpened input map \
              (shifted to new origin if original origin was not 0,0,0), \
              written out as map coefficients
      .short_caption = Sharpened map coeffs file
      .input_size = 400

    output_directory =  None
      .type = path
      .help = Directory where output files are to be written \
                applied.
      .short_caption = Output directory
      .style = directory

  }

  crystal_info
    .style = menu_item auto_align
  {

     use_sg_symmetry = False
       .type = bool
       .short_caption = Use space-group symmetry
       .help = If you set use_sg_symmetry=True then the symmetry of the space\
               group will be used. For example in P1 a point at one end of \
               the \
               unit cell is next to a point on the other end.  Normally for \
               cryo-EM data this should be set to False and for crystal data \
               it should be set to True.

     resolution = None
       .type = float
       .short_caption = Resolution
       .help = Optional nominal resolution of the map.

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
  }

  map_modification
    .style = menu_item auto_align
  {

     b_iso = None
       .type = float
       .short_caption = Target b_iso
       .help = Target B-value for map (sharpening will be applied to yield \
          this value of b_iso)

     b_sharpen = None
       .type = float
       .short_caption = Sharpening
       .help = Sharpen with this b-value. Contrast with b_iso that yield a \
           targeted value of b_iso

     resolution_dependent_b = None
       .type = floats
       .short_caption = resolution_dependent b
       .help = If set, apply resolution_dependent_b (b0 b1 b2). \
             Log10(amplitudes) will start at 1, change to b0 at half \
             of resolution specified, changing linearly, \
             change to b1 at resolution specified, \
             and change to b2 at high-resolution limit of map

     d_min_ratio = 0.833
       .type = float
       .short_caption = Sharpen d_min ratio
       .help = Sharpening will be applied using d_min equal to \
             d_min_ratio times resolution. Default is 0.833

     rmsd = None
       .type = float
       .short_caption = RMSD of model
       .help = RMSD of model to true model (if supplied).  Used to \
             estimate expected fall-of with resolution of correct part \
             of model-based map. If None, assumed to be resolution/3.

     fraction_complete = None
       .type = float
       .short_caption = Completeness model
       .help = Completness of model (if supplied).  Used to \
             estimate correct part \
             of model-based map. If None, estimated from max(FSC).

     auto_sharpen = None
       .type = bool
       .short_caption = Automatically determine sharpening
       .help = Automatically determine sharpening using kurtosis maximization\
                 or adjusted surface area. Default is True

     auto_sharpen_methods = no_sharpening b_iso b_iso_to_d_cut \
                            resolution_dependent *None
       .type = choice(multi=True)
       .short_caption = Sharpening methods
       .help = Methods to use in sharpening. b_iso searches for b_iso to \
          maximize sharpening target (kurtosis or adjusted_sa). \
          b_iso_to_d_cut applies b_iso only up to resolution specified, with \
          fall-over of k_sharpen.  Resolution dependent adjusts 3 parameters \
          to sharpen variably over resolution range. Default is all.

     box_in_auto_sharpen = None
       .type = bool
       .short_caption = Use box for auto_sharpening
       .help = Use a representative box of density for initial \
                auto-sharpening instead of the entire map. Default is True.

     max_box_fraction = None
       .type = float
       .short_caption = Max size of box for auto_sharpening
       .help = If box is greater than this fraction of entire map, use \
                entire map. Default is 0.5.

     k_sharpen = None
       .type = float
       .short_caption = sharpening transition
       .help = Steepness of transition between sharpening (up to resolution \
           ) and not sharpening (d < resolution).  Note: for blurring, \
           all data are blurred (regardless of resolution), while for \
           sharpening, only data with d about resolution or lower are \
           sharpened. This prevents making very high-resolution data too \
           strong.  Note 2: if k_sharpen is zero or None, then no \
           transition is applied and all data is sharpened or blurred. \
           Default is 10.


     maximum_low_b_adjusted_sa = 0.
       .type = float
       .short_caption = Max low-B adjusted_sa
       .help = Require adjusted surface area to be this value or less \
               when map is highly sharpened (at value of search_b_min).

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
            Default is 40. A smaller value will give more sharpening.

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
   }

  include scope libtbx.phil.interface.tracking_params

  gui
    .help = "GUI-specific parameter required for output directory"
  {
    output_dir = None
    .type = path
    .style = output_dir
  }

""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):

  command_line = iotbx.phil.process_command_line_with_files(
    reflection_file_def="input_files.map_coeffs_file",
    map_file_def="input_files.map_file",
    pdb_file_def="input_files.pdb_file",
    args=args,
    master_phil=master_phil)


  print >>out,"\nAuto-sharpen a map\n"
  params=command_line.work.extract()
  print >>out,"Command used: %s\n" %(
   " ".join(['phenix.auto_sharpen']+args))
  master_params.format(python_object=params).show(out=out)

  if params.output_files.output_directory is None:
    params.output_files.output_directory=os.getcwd()
  elif not os.path.isdir(params.output_files.output_directory):
    os.mkdir(params.output_files.output_directory)


  return params

def get_map_coeffs_from_file(
      map_coeffs_file=None,
      map_coeffs_labels=None):
    from iotbx import reflection_file_reader
    reflection_file=reflection_file_reader.any_reflection_file(
        map_coeffs_file)
    mtz_content=reflection_file.file_content()
    for ma in reflection_file.as_miller_arrays(merge_equivalents=True):
      if not ma.is_complex_array(): continue
      labels=",".join(ma.info().labels)
      if not map_coeffs_labels or labels==map_coeffs_labels:  # take it
         return ma

def get_map_and_model(params=None,out=sys.stdout):

  acc=None # accessor used to shift map back to original location if desired
  if params.input_files.map_file:
    from cctbx.maptbx.segment_and_split_map import get_map_object
    map_data,space_group,unit_cell,crystal_symmetry,origin_frac,acc=\
      get_map_object(file_name=params.input_files.map_file,out=out)
    map_data=map_data.as_double()
    if origin_frac != (0,0,0) and acc is None:
      print >>out,"\nWARNING: Unable to place output map at position of "+\
        "input map though input map has non-zero origin at %s\n" %(
        str(origin_frac))

  elif params.input_files.map_coeffs_file:
    map_coeffs=get_map_coeffs_from_file(
      map_coeffs_file=params.input_files.map_coeffs_file,
      map_coeffs_labels=params.input_files.map_coeffs_labels)

    if not map_coeffs:
      raise Sorry("Could not get map coeffs from %s with labels %s" %(
        params.input_files.map_coeffs_file,params.input_files.map_coeffs_labels))
    print >>out,"Map coefficients read from %s with labels %s" %(
         params.input_files.map_coeffs_file,
         str(params.input_files.map_coeffs_labels))
    crystal_symmetry=map_coeffs.crystal_symmetry()
    from cctbx.maptbx.segment_and_split_map import get_map_from_map_coeffs
    map_data=get_map_from_map_coeffs(
      map_coeffs=map_coeffs,crystal_symmetry=crystal_symmetry)
    acc=map_data.accessor()
    if not params.crystal_info.resolution:
      params.crystal_info.resolution=map_coeffs.d_min()
      print >>out,"Resolution from map_coeffs is %7.2f A" %(
          params.crystal_info.resolution)
  else:
    raise Sorry("Need ccp4 map or map_coeffs")

  if params.crystal_info.resolution is None:
    raise Sorry("Need resolution if map is supplied")

  if params.input_files.pdb_file: # get model
    model_file=params.input_files.pdb_file
    if not os.path.isfile(model_file):
      raise Sorry("Missing the model file: %s" %(model_file))
    pdb_inp=iotbx.pdb.input(file_name=model_file)
    if origin_frac != (0,0,0):
      print >>out,"Shifting model by %s" %(str(origin_frac))
      from cctbx.maptbx.segment_and_split_map import \
         apply_shift_to_pdb_hierarchy
      origin_shift=crystal_symmetry.unit_cell().orthogonalize(
         (-origin_frac[0],-origin_frac[1],-origin_frac[2]))
      pdb_inp=apply_shift_to_pdb_hierarchy(
       origin_shift=origin_shift,
       crystal_symmetry=crystal_symmetry,
       pdb_hierarchy=pdb_inp.construct_hierarchy(),
       out=out).as_pdb_input()

  else:
    pdb_inp=None

  return pdb_inp,map_data,crystal_symmetry,acc


def run(args,out=sys.stdout):
  # Get the parameters
  params=get_params(args,out=out)

  # get map_data and crystal_symmetry

  pdb_inp,map_data,crystal_symmetry,acc=get_map_and_model(
     params=params,out=out)

  # NOTE: map_data is now relative to origin at (0,0,0).
  # Use map_data.reshape(acc) to put it back where it was if acc is not None


  # auto-sharpen the map
  from cctbx.maptbx.segment_and_split_map import auto_sharpen_map_or_map_coeffs
  si=auto_sharpen_map_or_map_coeffs(
        resolution=params.crystal_info.resolution, # required
        crystal_symmetry=crystal_symmetry,
        map=map_data,
        solvent_content=params.crystal_info.solvent_content,
        box_in_auto_sharpen=params.map_modification.box_in_auto_sharpen,
        auto_sharpen_methods=params.map_modification.auto_sharpen_methods,
        residual_target=params.map_modification.residual_target,
        region_weight=params.map_modification.region_weight,
        sa_percent=params.map_modification.sa_percent,
        eps=params.map_modification.eps,
        n_bins=params.map_modification.n_bins,
        max_regions_to_test=params.map_modification.max_regions_to_test,
        fraction_occupied=params.map_modification.fraction_occupied,
        sharpening_target=params.map_modification.sharpening_target,
        d_min_ratio=params.map_modification.d_min_ratio,
        max_box_fraction=params.map_modification.max_box_fraction,
        k_sharpen=params.map_modification.max_box_fraction,
        search_b_min=params.map_modification.search_b_min,
        search_b_max=params.map_modification.search_b_max,
        search_b_n=params.map_modification.search_b_n,
        maximum_low_b_adjusted_sa=\
           params.map_modification.maximum_low_b_adjusted_sa,
        b_iso=params.map_modification.b_iso,
        b_sharpen=params.map_modification.b_sharpen,
        resolution_dependent_b=\
           params.map_modification.resolution_dependent_b,
        pdb_inp=pdb_inp,
        rmsd=params.map_modification.rmsd,
        fraction_complete=params.map_modification.fraction_complete,
        out=out)

  # get map_data and map_coeffs of final map

  new_map_data=si.as_map_data()
  new_map_coeffs=si.as_map_coeffs()

  if not si.is_model_sharpening():
    print >>out
    print >>out,80*"=","\n",80*"="
    print >>out,"\n           Final sharpening information\n "
    si.show_summary(verbose=params.control.verbose,out=out)
    print >>out,80*"=","\n",80*"="

  # write out the new map_coeffs and map if requested:

  if params.output_files.sharpened_map_file:
    output_map_file=os.path.join(params.output_files.output_directory,
        params.output_files.sharpened_map_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    offset_map_data=new_map_data.deep_copy()
    if acc is not None:  # offset the map to match original if possible
      offset_map_data.reshape(acc)
      print >>out,\
       "\nWrote sharpened map in original location with origin at %s\nto %s" %(
         str(offset_map_data.origin()),output_map_file)
    else:
      print >>out,"\nWrote sharpened map with origin at 0,0,0 "+\
        "(NOTE: may not be \nsame as original location) to %s\n" %(
         output_map_file)
    write_ccp4_map(crystal_symmetry, output_map_file, offset_map_data)

  if params.output_files.shifted_sharpened_map_file:
    output_map_file=os.path.join(params.output_files.output_directory,
        params.output_files.shifted_sharpened_map_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    write_ccp4_map(crystal_symmetry, output_map_file, new_map_data)
    print >>out,"\nWrote sharpened map (origin at %s)\nto %s" %(
     str(new_map_data.origin()),output_map_file)

  if params.output_files.sharpened_map_coeffs_file:
    output_map_coeffs_file=os.path.join(params.output_files.output_directory,
        params.output_files.sharpened_map_coeffs_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    new_map_coeffs.as_mtz_dataset(column_root_label='FWT').mtz_object().write(
       file_name=output_map_coeffs_file)
    print >>out,"\nWrote sharpened map_coeffs (origin at 0,0,0)\n to %s\n" %(
       output_map_coeffs_file)

  return new_map_data,new_map_coeffs,crystal_symmetry,si

# =============================================================================
# GUI-specific bits for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    import os
    from wxGUI2 import utils
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args=self.args, out=sys.stdout)
    return result

def validate_params(params):
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
