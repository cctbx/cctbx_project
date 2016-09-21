from __future__ import division
import sys, os
import iotbx.phil
from libtbx.utils import Sorry

master_phil = iotbx.phil.parse("""

  input_files {

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    map_coeffs_file = None
      .type = path
      .help = Optional file with map coefficients
      .short_caption = Map coefficients
      .style = bold file_type:hkl input_file process_hkl child:fobs:data_labels\
        child:space_group:space_group child:unit_cell:unit_cell

    map_coeffs_labels = None
      .type = str
      .input_size = 160
      .help = Optional label specifying which columns of of map coefficients \
          to use
      .short_caption = Map coeffs label
      .style = bold renderer:draw_fobs_label_widget


  }

  output_files {

    shifted_map_file = shifted_map.ccp4
      .type = path
      .help = Input map file shifted to new origin.
      .short_caption = Shifted map file

    shifted_sharpened_map_file = shifted_sharpened_map.ccp4
      .type = path
      .help = Input map file shifted to new origin and sharpened.
      .short_caption = Shifted sharpened map file

    shifted_sharpened_map_coeffs_file = shifted_sharpened_map_coeffs.mtz
      .type = path
      .help = Input map coeffs shifted to new origin and sharpened.
      .short_caption = Shifted sharpened map coeffs file

    output_directory =  None
      .type = path
      .help = Directory where output files are to be written \
                applied.
      .short_caption = Output directory

  }

  crystal_info {

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
       .short_caption = resolution
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

  map_modification {

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

   control {
     verbose = False
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output

     resolve_size = None
        .type = int
        .help = "Size of resolve to use. "
        .style = hidden
   }
""", process_includes=True)
master_params = master_phil

def get_params(args,out=sys.stdout):

  command_line = iotbx.phil.process_command_line_with_files(
    reflection_file_def="input_files.map_coeffs_file",
    map_file_def="input_files.map_file",
    args=args,
    master_phil=master_phil)


  print >>out,"\nAuto-sharpen a map\n"
  params=command_line.work.extract()
  print >>out,"Command used: %s\n" %(
   " ".join(['phenix.auto_sharpen']+args))
  master_params.format(python_object=params).show(out=out)

  if params.output_files.output_directory is None:
    params.output_files.output_directory=os.getcwd()

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

def get_map(params=None,out=sys.stdout):
  if params.input_files.map_file:
    from cctbx.maptbx.segment_and_split_map import get_map_object
    map_data,space_group,unit_cell,crystal_symmetry,origin_frac=get_map_object(
       file_name=params.input_files.map_file,out=out)
    map_data=map_data.as_double()

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
    if not params.crystal_info.resolution:
      params.crystal_info.resolution=map_coeffs.d_min()
      print >>out,"Resolution from map_coeffs is %7.2f A" %(
          params.crystal_info.resolution)
  else:
    raise Sorry("Need ccp4 map or map_coeffs")

  if params.crystal_info.resolution is None:
    raise Sorry("Need resolution if map is supplied")

  return map_data,crystal_symmetry


def run(args,out=sys.stdout):
  # Get the parameters
  params=get_params(args,out=out)

  # get map_data and crystal_symmetry
  map_data,crystal_symmetry=get_map(params=params,out=out)


  # auto-sharpen the map
  from cctbx.maptbx.segment_and_split_map import auto_sharpen_map_or_map_coeffs
  new_map_data=auto_sharpen_map_or_map_coeffs(
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
        b_iso=params.map_modification.b_iso,
        b_sharpen=params.map_modification.b_sharpen,
        resolution_dependent_b=\
           params.map_modification.resolution_dependent_b,
        out=out)

  # convert to map_coeffs also

  from cctbx.maptbx.segment_and_split_map import get_f_phases_from_map
  new_map_coeffs=get_f_phases_from_map(map_data=new_map_data,
       crystal_symmetry=crystal_symmetry,
       d_min=params.crystal_info.resolution,
       d_min_ratio=params.map_modification.d_min_ratio,
       return_as_map_coeffs=True,
       out=out)

  # write out the new map_coeffs and map if requested:

  if params.output_files.shifted_sharpened_map_file:
    output_map_file=os.path.join(params.output_files.output_directory,
        params.output_files.shifted_sharpened_map_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    write_ccp4_map(crystal_symmetry, output_map_file, new_map_data)
    print >>out,"\nWrote sharpened map to %s" %(output_map_file)

  if params.output_files.shifted_sharpened_map_coeffs_file:
    output_map_coeffs_file=os.path.join(params.output_files.output_directory,
        params.output_files.shifted_sharpened_map_coeffs_file)
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    new_map_coeffs.as_mtz_dataset(column_root_label='FWT').mtz_object().write(
       file_name=output_map_coeffs_file)
    print >>out,"\nWrote sharpened map_coeffs to %s" %(output_map_coeffs_file)

  return new_map_data,new_map_coeffs,crystal_symmetry

if __name__=="__main__":
  run(sys.argv[1:])
