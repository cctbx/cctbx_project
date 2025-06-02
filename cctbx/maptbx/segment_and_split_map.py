"""
Methods for map segmentation and auto-sharpening.
"""
from __future__ import absolute_import, division, print_function
from operator import itemgetter
import iotbx.phil
import iotbx.mrcfile
from cctbx import crystal
from cctbx import maptbx
from libtbx.utils import Sorry
import sys, os
from cctbx.array_family import flex
from scitbx.math import matrix
from copy import deepcopy
from libtbx.utils import null_out
import libtbx.callbacks # import dependency
from libtbx import group_args
from six.moves import range
from six.moves import zip
from cctbx.development.create_models_or_maps import get_map_from_map_coeffs

master_phil = iotbx.phil.parse("""

  input_files {

    seq_file = None
       .type = path
       .short_caption = Sequence file
       .help = Sequence file (unique chains only,  \
               1-letter code, chains separated by \
               blank line or greater-than sign.)  \
               Can have chains that are DNA/RNA/protein and\
               all can be present in one file. \
               If not supplied, must supply molecular mass or \
               solvent content.

    map_file = None
      .type = path
      .help = File with CCP4-style map
      .short_caption = Map file

    half_map_file = None
      .type = path
      .multiple = True
      .short_caption = Half map
      .help = Half map (two should be supplied) for FSC calculation. Must \
               have grid identical to map_file

    external_map_file = None
      .type = path
      .short_caption = External map
      .style = file_type:ccp4_map bold input_file
      .help = External map to be used to scale map_file (power vs resolution\
              will be matched). Not used in segment_and_split_map

    ncs_file = None
      .type = path
      .help = File with symmetry information (typically point-group NCS with \
               the center specified). Typically in  PDB format. \
              Can also be a .ncs_spec file from phenix. \
              Created automatically if symmetry is specified.
      .short_caption = symmetry file

    pdb_file = None
      .type = path
      .help = Optional PDB file matching map_file to be offset

    pdb_to_restore = None
      .type = path
      .help = Optional PDB file to restore to position matching original \
              map_file.  Used in combination with info_file = xxx.pkl \
              and restored_pdb = yyyy.pdb
      .short_caption = PDB to restore

    info_file = None
      .type = path
      .help = Optional pickle file with information from a previous run.\
              Can be used with pdb_to_restore to restore a PDB file to \
              to position matching original \
              map_file.
      .short_caption = Info file

    target_ncs_au_file = None
      .type = path
      .help = Optional PDB file to partially define the ncs asymmetric \
               unit of the map. The coordinates in this file will be used \
               to mark part of the ncs au and all points nearby that are \
               not part of another ncs au will be added.

     input_weight_map_pickle_file = None
       .type = path
       .short_caption = Input weight map pickle file
       .help = Weight map pickle file
  }

  output_files {

    magnification_map_file = magnification_map.ccp4
      .type = path
      .help = Input map file with magnification applied.  Only written if\
                magnification is applied.
      .short_caption = Magnification map file

    magnification_ncs_file = magnification_ncs.ncs_spec
      .type = path
      .help = Input NCS with magnification applied.  Only written if\
                magnification is applied.
      .short_caption = Magnification NCS file

    shifted_map_file = shifted_map.ccp4
      .type = path
      .help = Input map file shifted to new origin.
      .short_caption = Shifted map file

    shifted_sharpened_map_file = shifted_sharpened_map.ccp4
      .type = path
      .help = Input map file shifted to new origin and sharpened.
      .short_caption = Shifted sharpened map file

    sharpened_map_file = sharpened_map.ccp4
      .type = str
      .short_caption = Sharpened map file
      .help = Output sharpened map file, superimposed on the original map.
      .input_size = 400

    shifted_pdb_file = shifted_pdb.pdb
      .type = path
      .help = Input pdb file shifted to new origin.
      .short_caption = Shifted pdb file

    shifted_ncs_file = shifted_ncs.ncs_spec
      .type = path
      .help = NCS information shifted to new origin.
      .short_caption = Output NCS info file

    shifted_used_ncs_file = shifted_used_ncs.ncs_spec
      .type = path
      .help = NCS information (just the part that is used) shifted \
               to new origin.
      .short_caption = Output used NCS info file

    output_directory = segmented_maps
      .type = path
      .help = Directory where output files are to be written \
                applied.
      .short_caption = Output directory

    box_map_file = box_map_au.ccp4
      .type = path
      .help = Output map file with one NCS asymmetric unit, cut out box
      .short_caption = Box NCS map file

    box_mask_file = box_mask_au.ccp4
      .type = path
      .help = Output mask file with one NCS asymmetric unit, cut out box
      .short_caption = Box NCS mask file

    box_buffer = 5
      .type = int
      .help = Buffer (grid units) around NCS asymmetric unit in box_mask and map
      .short_caption = Box buffer size

    au_output_file_stem = shifted_au
      .type = str
      .help = File stem for output map files with one NCS asymmetric unit
      .short_caption = Output au file stem

    write_intermediate_maps = False
      .type = bool
      .help = Write out intermediate maps and masks for visualization
      .short_caption = Write intermediate maps

    write_output_maps = True
      .type = bool
      .help = Write out maps
      .short_caption = Write maps

    remainder_map_file = remainder_map.ccp4
      .type = path
      .help = output map file with remainder after initial regions identified
      .short_caption = Output remainder map file

    output_info_file = segment_and_split_map_info.pkl
      .type = path
      .help = Output pickle file with information about map and masks
      .short_caption = Output pickle file

    restored_pdb = None
      .type = path
      .help = Output name of PDB restored to position matching original \
              map_file.  Used in combination with info_file = xxx.pkl \
              and pdb_to_restore = xxxx.pdb
      .short_caption = Restored PDB file

    output_weight_map_pickle_file = weight_map_pickle_file.pkl
       .type = path
       .short_caption = Output weight map pickle file
       .help = Output weight map pickle file
  }

  crystal_info {

     chain_type = *None PROTEIN RNA DNA
       .type = choice
       .short_caption = Chain type
       .help = Chain type. Determined automatically from sequence file if \
               not given. Mixed chain types are fine (leave blank if so).

     sequence = None
       .type = str
       .short_caption = Sequence
       .help = Sequence as string

     is_crystal = None
       .type = bool
       .short_caption = Is a crystal
       .help = Defines whether this is a crystal (or cryo-EM).\
                Default is True if use_sg_symmetry = True and False otherwise.

     use_sg_symmetry = False
       .type = bool
       .short_caption = Use space-group symmetry
       .help = If you set use_sg_symmetry = True then the symmetry of the space\
               group will be used. For example in P1 a point at one end of \
               the \
               unit cell is next to a point on the other end.  Normally for \
               cryo-EM data this should be set to False and for crystal data \
               it should be set to True. This will normally also set the \
               value of is_crystal (same value as use_sg_symmetry) and \
               restrict_map_size (False if use_sg_symmetry = True).

     resolution = None
       .type = float
       .short_caption = resolution
       .help = Nominal resolution of the map. This is used later to decide on\
               resolution cutoffs for Fourier inversion of the map. Note: \
               the resolution is not cut at this value, it is cut at \
               resolution*d_min_ratio if at all.

     space_group = None
       .type = space_group
       .help = Space group (used for boxed maps)
       .style = hidden

     unit_cell = None
       .type = unit_cell
       .help = Unit Cell (used for boxed maps)
       .style = hidden

     original_unit_cell = None
       .type = unit_cell
       .help = Original unit cell (of input map). Used internally
       .style = hidden

     original_unit_cell_grid = None
       .type = ints
       .help = Original unit cell grid (of input map). Used internally
       .style = hidden

     molecular_mass = None
       .type = float
       .help = Molecular mass of molecule in Da. Used as alternative method \
                 of specifying solvent content.
       .short_caption = Molecular mass in Da

     solvent_content = None
       .type = float
       .help = Solvent fraction of the cell. Used for ID of \
               solvent content in boxed maps.
       .short_caption = Solvent content

     solvent_content_iterations = 3
       .type = int
       .help = Iterations of solvent fraction estimation. Used for ID of \
               solvent content in boxed maps.
       .short_caption = Solvent fraction iterations
       .style = hidden

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

  }

  reconstruction_symmetry {

     symmetry = None
       .type = str
       .short_caption = Symmetry type
       .help = Symmetry used in reconstruction. For example D7, C3, C2\
          I (icosahedral), T (tetrahedral), or ANY (try everything and \
          use the highest symmetry found). Not needed if ncs_file is supplied. \

     include_helical_symmetry = True
       .type = bool
       .short_caption = Include helical symmetry
       .help = You can include or exclude searches for helical symmetry

     must_be_consistent_with_space_group_number = None
       .type = int
       .short_caption = Space group to match
       .help = Searches for symmetry must be compatible with this space group\
               number.

     symmetry_center = None
       .type = floats
       .short_caption = symmetry center
       .help = Center (in A) for symmetry operators (if symmetry is found \
          automatically). \
          If set to None, first guess is the center of the cell and then \
          if that fails, found automatically as the center of the \
          density in the map.

     optimize_center = False
       .type = bool
       .short_caption = Optimize symmetry center
       .help = Optimize position of symmetry center. Also checks for center \
                at (0, 0, 0) vs center of map

     helical_rot_deg = None
       .type = float
       .short_caption = helical rotation
       .help = helical rotation about z in degrees

     helical_trans_z_angstrom = None
       .type = float
       .short_caption = helical translation
       .help = helical translation along z in Angstrom units

     max_helical_optimizations = 2
       .type = int
       .short_caption = Max helical optimizations
       .help = Number of optimizations of helical parameters\
               when finding symmetry

     max_helical_ops_to_check = 5
       .type = int
       .short_caption = Max helical ops to check
       .help = Number of helical operations in each direction to check \
               when finding symmetry

     max_helical_rotations_to_check = None
       .type = int
       .short_caption = Max helical rotations
       .help = Number of helical rotations to check \
               when finding symmetry

     two_fold_along_x = None
       .type = bool
       .short_caption = D two-fold along x
       .help = Specifies if D or I two-fold is along x (True) or y (False). \
               If None, both are tried.

     smallest_object = None
       .type = float
       .short_caption = Smallest object to consider
       .help = Dimension of smallest object to consider\
               when finding symmetry. Default is 5 * resolution

     score_basis = ncs_score cc *None
       .type = choice
       .short_caption = Symmetry score basis
       .help = Symmetry score basis. Normally ncs_score (sqrt(n)* cc) is \
               used except for identification of helical symmetry

     scale_weight_fractional_translation =  1.05
       .type = float
       .short_caption = Scale on fractional translation
       .help =  Give slight increase in weighting in helical symmetry \
               search to translations that are a fraction (1/2, 1/3) of \
               the d-spacing of the peak of intensity in the fourier \
               transform of the density.


     random_points = 100
       .type = int
       .short_caption = Random points
       .help = Number of random points in map to examine in finding symmetry

     identify_ncs_id = True
       .type = bool
       .short_caption = Identify NCS ID
       .help = If symmetry is not point-group symmetry, try each possible \
               operator when evaluating symmetry and choose the one that  \
               results in the most uniform density at symmetry-related points.

     min_ncs_cc = 0.75
       .type = float
       .short_caption = Minimum symmetry CC to keep it
       .help =  Minimum symmetry CC to keep operators when identifying \
                 automatically

     n_rescore = 5
       .type = int
       .short_caption = symmetry operators to rescore
       .help = Number of symmetry operators to rescore

     op_max = 14
       .type = int
       .short_caption = Max operators to try
       .help = If symmetry is ANY, try up to op_max-fold symmetries


    tol_r = 0.02
      .type = float
      .help = tolerance in rotations for point group or helical symmetry
      .short_caption = Rotation tolerance

    abs_tol_t = 2
      .type = float
      .help = tolerance in translations (A) for point group or helical symmetry
      .short_caption = Translation tolerance absolute

    max_helical_operators =  None
      .type = int
      .help = Maximum helical operators (if extending existing helical\
             operators)
      .short_caption = Maximum helical operators

    rel_tol_t = .05
      .type = float
      .help = tolerance in translations (fractional) for point group or \
             helical symmetry
      .short_caption = Translation tolerance fractional

    require_helical_or_point_group_symmetry = False
      .type = bool
      .help = normally helical or point-group symmetry (or none) is expected. \
             However in some cases (helical + rotational symmetry for \
             example) this is not needed and is not the case.
      .short_caption = Require helical or point-group or no symmetry


     }

  map_modification {

     magnification = None
       .type = float
       .short_caption = Magnification
       .help = Magnification to apply to input map.  Input map grid will be \
                scaled by magnification factor before anything else is done.

     b_iso = None
       .type = float
       .short_caption = Target b_iso
       .help = Target B-value for map (sharpening will be applied to yield \
          this value of b_iso). If sharpening method is not supplied, \
          default is to use b_iso_to_d_cut sharpening.

     b_sharpen = None
       .type = float
       .short_caption = Sharpening
       .help = Sharpen with this b-value. Contrast with b_iso that yield a \
           targeted value of b_iso. B_sharpen greater than zero is sharpening.\
           Less than zero is blurring.

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

     regions_to_keep = None
       .type = int
       .short_caption = Regions to keep
       .help = You can specify a limit to the number of regions to keep\
                when generating the asymmetric unit of density.

     auto_sharpen = True
       .type = bool
       .short_caption = Automatically determine sharpening
       .help = Automatically determine sharpening using kurtosis maximization\
                 or adjusted surface area

     auto_sharpen_methods = no_sharpening b_iso *b_iso_to_d_cut \
                            resolution_dependent model_sharpening \
                            half_map_sharpening target_b_iso_to_d_cut None

       .type = choice(multi = True)
       .short_caption = Sharpening methods
       .help = Methods to use in sharpening. b_iso searches for b_iso to \
          maximize sharpening target (kurtosis or adjusted_sa). \
          b_iso_to_d_cut applies b_iso only up to resolution specified, with \
          fall-over of k_sharpen.  Resolution dependent adjusts 3 parameters \
          to sharpen variably over resolution range. Default is \
          b_iso_to_d_cut .  target_b_iso_to_d_cut uses target_b_iso_ratio \
          to set b_iso.

     box_in_auto_sharpen = False
       .type = bool
       .short_caption = Use box for auto_sharpening
       .help = Use a representative box of density for initial \
                auto-sharpening instead of the entire map.

     density_select_in_auto_sharpen = True
       .type = bool
       .short_caption = density_select to choose box
       .help = Choose representative box of density for initial \
                auto-sharpening with density_select method \
                (choose region where there is high density). \
               Normally use this as well as density_select = True which \
               carries out density_select at start of segmentation.

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
               b_iso is set. Default is to set box_n_auto_sharpen = False \
               if b_iso is set.

     soft_mask = True
       .type = bool
       .help = Use soft mask (smooth change from inside to outside with radius\
             based on resolution of map). Required if you use half-map \
             sharpening without a model, otherwise optional.
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
               if symmetry is present.\
               If local_aniso_in_local_sharpening is True and symmetry is \
               present this can distort the map for some symmetry copies \
               because an anisotropy correction is applied\
               based on local density in one copy and is transferred without \
               rotation to other copies.

     local_aniso_in_local_sharpening = None
       .type = bool
       .short_caption = Local anisotropy
       .help = Use local anisotropy in local sharpening.  \
               Default is True unless symmetry is present.

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

     box_size = 40 40 40
       .type = ints
       .short_caption = Size of box
       .help = You can specify the size of the boxes to use (grid units)

     target_n_overlap = 10
       .type = int
       .short_caption = Target overlap of boxes
       .help = You can specify the targeted overlap of boxes in local \
           sharpening

     restrict_map_size = None
       .type = bool
       .short_caption = Restrict box map size
       .help = Restrict box map to be inside full map (required for cryo-EM data).\
               Default is True if use_sg_symmetry = False.

     restrict_z_turns_for_helical_symmetry = 1
       .type = float
       .short_caption = Restrict Z turns for helical symmetry
       .help = Restrict Z turns for helical symmetry.  Number of \
               turns of helix going each direction in Z is specified.

     restrict_z_distance_for_helical_symmetry = None
       .type = float
       .short_caption = Restrict Z distance for helical symmetry
       .help = Restrict Z distance (+/- this distance from center) \
              for helical symmetry.

     remove_aniso = True
       .type = bool
       .short_caption = Remove aniso
       .help = You can remove anisotropy (overall and locally) during sharpening

     cc_cut = 0.2
       .type = float
       .short_caption = Min reliable CC in half-maps
       .help = Estimate of minimum highly reliable CC in half-map FSC. Used\
               to decide at what CC value to smooth the remaining CC values.

     max_cc_for_rescale = 0.2
       .type = float
       .short_caption = Max CC for rescale
       .help =  Min reliable CC in half-maps. \
               Used along with cc_cut and scale_using_last to correct for \
               small errors in FSC estimation at high resolution.  If the \
               value of FSC near the high-resolution limit is above \
               max_cc_for_rescale, assume these values are correct and do not \
               correct them. To keep all original values use\
               max_cc_for_rescale = 1

     scale_using_last = 3
       .type = int
       .short_caption = Last N bins in FSC assumed to be about zero
       .help = If set, assume that the last scale_using_last bins in the FSC \
          for half-map or model sharpening are about zero (corrects for  \
          errors in the half-map process).

     max_box_fraction = 0.5
       .type = float
       .short_caption = Max size of box for auto_sharpening
       .help = If box is greater than this fraction of entire map, use \
                entire map.

     density_select_max_box_fraction = 0.95
       .type = float
       .short_caption = Max size of box for density_select
       .help = If box is greater than this fraction of entire map, use \
                entire map for density_select. Default is 0.95

     mask_atoms = True
       .type = bool
       .short_caption = Mask atoms
       .help = Mask atoms when using model sharpening

     mask_atoms_atom_radius = 3
       .type  = float
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
           strong.  Note 2: if k_sharpen is zero, then no \
           transition is applied and all data is sharpened or blurred. \
           Note 3: only used if b_iso is set.

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
                b_iso. Not normally carried out.

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
               formula b_iso = 5.9*d_min**2

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

     region_weight_default = 30.
       .type = float
       .short_caption = Region weight default
       .help = Region_weight adjusted to be region_weight_default\
               if no information available

     target_b_iso_ratio = 5.9
       .type = float
       .short_caption = Target b_iso ratio
       .help = Target b_iso ratio : b_iso is estimated as \
               target_b_iso_ratio * resolution**2

     signal_min = 3.0
       .type = float
       .short_caption = Minimum signal
       .help = Minimum signal in estimation of optimal b_iso.  If\
                not achieved, use any other method chosen.

     target_b_iso_model_scale = 0.
       .type = float
       .short_caption = scale on target b_iso ratio for model
       .help = For model sharpening, the target_biso is scaled \
                (normally zero).

     search_b_min = -100
       .type = float
       .short_caption = Low bound for b_iso search
       .help = Low bound for b_iso search.

     search_b_max = 300
       .type = float
       .short_caption = High bound for b_iso search
       .help = High bound for b_iso search.

     search_b_n = 21
       .type = int
       .short_caption = Number of b_iso values to search
       .help = Number of b_iso values to search.

     residual_target = 'adjusted_sa'
       .type = str
       .short_caption = Residual target
       .help = Target for maximization steps in sharpening.  \
          Can be kurtosis or adjusted_sa (adjusted surface area)

     sharpening_target = 'adjusted_sa'
       .type = str
       .short_caption = Overall sharpening target
       .help = Overall target for sharpening.  Can be kurtosis or adjusted_sa \
          (adjusted surface area) or adjusted_path_length.  \
            Used to decide which sharpening approach \
          is used. Note that during optimization, residual_target is used \
          (they can be the same.)

     region_weight = 40
       .type = float
       .short_caption = Region weighting
       .help = Region weighting in adjusted surface area calculation.\
            Score is surface area minus region_weight times number of regions.\
            Default is 40. A smaller value will give more sharpening.

     sa_percent = 30.
       .type = float
       .short_caption = Percent of target regions in adjusted_sa
       .help = Percent of target regions used in calulation of adjusted \
         surface area.  Default is 30.

     fraction_occupied = 0.20
       .type = float
       .short_caption = Fraction of molecular volume inside contours
       .help = Fraction of molecular volume targeted to be inside contours. \
           Used to set contour level. Default is 0.20

      n_bins = 20
        .type = int
        .short_caption = Resolution bins
        .help = Number of resolution bins for sharpening. Default is 20.

      max_regions_to_test = 30
        .type = int
        .short_caption = Max regions to test
        .help = Number of regions to test for surface area in adjusted_sa \
                scoring of sharpening

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
  }

  segmentation {

    select_au_box = None
      .type = bool
      .help = Select box containing at least one representative region of \
            the map. Also select just symmetry operators relevant to that box. \
              Default is true if number of operators is at least \
              n_ops_to_use_au_box
      .short_caption = select au box

    n_ops_to_use_au_box = 25
      .type = int
      .help = If number of operators is this big or more and \
              select_au_box is None, set it to True.
      .short_caption = N ops to use au_box

    n_au_box = 5
      .type = int
      .help = Number of symmetry copies to try and get inside au_box
      .short_caption = N au box


    lower_bounds = None
      .type = ints
      .help = You can select a part of your map for analysis with \
              lower_bounds and upper_bounds.
      .short_caption = Lower bounds

    upper_bounds = None
      .type = ints
      .help = You can select a part of your map for analysis with \
              lower_bounds and upper_bounds.
      .short_caption = Upper bounds


    density_select = True
      .type = bool
      .help = Run map_box with density_select = True to cut out the region \
              in the input map that contains density. Useful if the input map \
              is much larger than the structure. Done before segmentation is\
              carried out.
      .short_caption = Trim map to density

    density_select_threshold = 0.05
      .type = float
      .help = Choose region where density is this fraction of maximum or greater
      .short_caption = threshold for density_select

    get_half_height_width = None
      .type = bool
      .help = Use 4 times half-width at half-height as estimate of max size
      .short_caption = Half-height width estimation

    box_ncs_au = True
      .type = bool
      .help = Box the map containing just the au of the map
      .short_caption = Box NCS au

    cell_cutoff_for_solvent_from_mask = 150
      .type = float
      .help = For cells with average edge over this cutoff, use the\
               low resolution mask (backup) method for solvent estimation
      .short_caption = Cell cutoff for solvent_from_mask

    mask_padding_fraction = 0.025
      .type = float
      .help = Adjust threshold of standard deviation map in low resolution \
            mask identification of solvent content to give this much more \
            inside mask than would be obtained with the value of\
             fraction_of_max_mask_threshold.
      .short_caption = Mask padding fraction

    fraction_of_max_mask_threshold = .05
      .type = float
      .help = threshold of standard deviation map in low resolution mask \
             identification of solvent content.
      .short_caption = Fraction of max mask_threshold

    mask_threshold = None
      .type = float
      .help = threshold in identification of overall mask. If None, guess \
               volume of molecule from sequence and symmetry copies.
      .short_caption = Density select threshold

    grid_spacing_for_au = 3
      .type = int
      .help = Grid spacing for asymmetric unit when constructing asymmetric unit.
      .short_caption = Grid spacing for constructing asymmetric unit

    radius = None
      .type = float
      .help = Radius for constructing asymmetric unit.
      .short_caption = Radius for constructing asymmetric unit


    value_outside_mask = 0.0
      .type = float
      .help = Value to assign to density outside masks
      .short_caption = Value outside mask

    density_threshold = None
      .type = float
      .short_caption = Density threshold
      .help = Threshold density for identifying regions of density. \
             Applied after normalizing the density in the region of \
             the molecule to an rms of 1 and mean of zero.

    starting_density_threshold = None
      .type = float
      .short_caption = Starting density threshold
      .help = Optional guess of threshold density

    iteration_fraction = 0.2
      .type = float
      .short_caption = Iteration fraction
      .help = On iteration of finding regions, assume target volume is \
              this fraction of the value on previous iteration

    max_overlap_fraction = 0.05
      .type = float
      .short_caption = Max overlap
      .help = Maximum fractional overlap allowed to density in another \
              asymmetric unit. Definition of a bad region.

    remove_bad_regions_percent = 1
      .type = float
      .short_caption = Remove worst overlapping regions
      .help = Remove the worst regions that are part of more than one NCS \
           asymmetric unit, up to remove_bad_regions_percent of the total

    require_complete = True
      .type = bool
      .short_caption = Require all symmetry copies to be represented for a region
      .help =  Require all symmetry copies to be represented for a region

    split_if_possible = True
      .type = bool
      .short_caption = Split regions if mixed
      .help = Split regions that are split in some symmetry copies.\
              If None, split if most copies are split.

    write_all_regions = False
      .type = bool
      .short_caption = Write all regions
      .help = Write all regions to ccp4 map files.

    max_per_au = None
      .type = int
      .short_caption = Max regions in au
      .help = Maximum number of regions to be kept in the NCS asymmetric unit

    max_per_au_ratio = 5.
      .type = int
      .short_caption = Max ratio of regions to expected
      .help = Maximum ratio of number of regions to be kept in the \
         NCS asymmetric unit to those expected

    min_ratio_of_ncs_copy_to_first = 0.5
      .type = float
      .short_caption = Minimum ratio of ncs copy to first
      .help = Minimum ratio of the last ncs_copy region size to maximum

    min_ratio = 0.1
      .type = float
      .short_caption = Minimum ratio to keep
      .help = Minimum ratio of region size to maximum to keep it

    max_ratio_to_target = 3
      .type = float
      .help = Maximum ratio of grid points in top region to target
      .short_caption = Max ratio to target

    min_ratio_to_target = 0.3
      .type = float
      .help = Minimum ratio of grid points in top region to target
      .short_caption = Min ratio to target

    min_volume = 10
      .type = int
      .help = Minimum region size to consider (in grid points)
      .short_caption = Minimum region size

    residues_per_region = 50
      .type = float
      .help = Target number of residues per region
      .short_caption = Residues per region

    seeds_to_try = 10
      .type = int
      .help = Number of regions to try as centers
      .short_caption = Seeds to try

    iterate_with_remainder = True
      .type = bool
      .short_caption = Iterate
      .help = Iterate looking for regions based on remainder from first analysis

    weight_rad_gyr = 0.1
      .type = float
      .short_caption = Weight on radius of gyration
      .help = Weight on radius of gyration of group of regions in NCS AU \
               relative to weight on closeness to neighbors.  Normalized to\
               largest cell dimension with weight = weight_rad_gyr*300/cell_max

    expand_size = None
      .type = int
      .help = Grid points to expand size of regions when excluding for next \
               round. If None, set to approx number of grid points to get \
               expand_target below
      .short_caption = Expand size

    expand_target = 1.5
      .type = float
      .help = Target expansion of regions (A)
      .short_caption = Expand target

    mask_additional_expand_size = 1
      .type = int
      .help = Mask expansion in addition to expand_size for final map
      .short_caption = Mask additional expansion

    mask_expand_ratio = 1
      .type = int
      .help = Mask expansion relative to resolution for save_box_map_ncs_au
      .short_caption = Mask expand ratio

    exclude_points_in_ncs_copies = True
      .type = bool
      .help = Exclude points that are in symmetry copies when creating NCS au. \
               Does not apply if add_neighbors = True
      .short_caption = Exclude points in symmetry copies

    add_neighbors = True
      .type = bool
      .help = Add neighboring regions around the au. Turns off \
           exclude_points_in_ncs_copies also.
      .short_caption = Add neighbors

    add_neighbors_dist = 1.
      .type = float
      .help = Max increase in radius of gyration by adding region to keep it.
      .short_caption = Add neighbors dist
  }

   control {
     verbose = False
        .type = bool
        .help = '''Verbose output'''
        .short_caption = Verbose output

     shift_only = None
        .type = bool
        .short_caption = Shift only
        .help = Shift map and half_maps and stop

     sharpen_only = None
        .type = bool
        .short_caption = Sharpen only
        .help = Sharpen map and stop

     check_ncs = None
       .type = bool
       .short_caption = Check NCS
       .help = Check the NCS symmetry by estimating NCS correlation and stop

     resolve_size = None
        .type = int
        .help = "Size of resolve to use. "
        .style = hidden

      quick = True
        .type = bool
        .help = Run quickly if possible
        .short_caption = Quick run

     memory_check = True
        .type = bool
        .help = Map-to-model checks to make sure you have enough memory on \
                  your machine to run.  You can disable this by setting this \
                  keyword to False. The estimates are approximate so it is \
                  possible your job could run even if the check fails.  Note \
                  the check does not take any other uses of the memory on \
                  your machine into account.
        .short_caption = Memory check

     save_box_map_ncs_au = False
       .type = bool
       .help = Controls whether the map_box ncs_au is saved. Internal use only
       .style = hidden

     write_files = True
       .type = bool
       .help = Controls whether files are written
       .short_caption = Write files

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
""", process_includes = True)
master_params = master_phil

class map_and_b_object:
  """Holder for map_data and starting and final b_iso values"""
  def __init__(self,
       map_data = None,
      starting_b_iso = None,
      final_b_iso = None):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())


class pdb_info_object:
  """Holder for PDB file name and number of residues"""
  def __init__(self,
    file_name = None,
    n_residues = None,
    ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    import time
    self.init_asctime = time.asctime()

  def show_summary(self, out = sys.stdout):
    """Summarize pdb_info_object"""
    print("PDB file:%s" %(self.file_name), end = ' ', file = out)
    if self.n_residues:
      print("   Residues: %d" %(self.n_residues), file = out)
    else:
      print(file = out)

class seq_info_object:
  """Holder for sequence file name, sequence and number of residues"""
  def __init__(self,
    file_name = None,
    sequence = None,
    n_residues = None,
    ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    import time
    self.init_asctime = time.asctime()

  def show_summary(self, out = sys.stdout):
    """Summarize seq_info_object"""
    if self.file_name:
      print("Sequence file:%s" %(self.file_name), end = ' ', file = out)
    if self.n_residues:
      print("   Residues: %d" %(self.n_residues), file = out)
    else:
      print(file = out)


class ncs_info_object:
  """Holder for NCS file name, number of operators, whether it is helical
    symmetry, original number of operators"""
  def __init__(self,
    file_name = None,
    number_of_operators = None,
    is_helical_symmetry = None,
    original_number_of_operators = None,
    ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    import time
    self.init_asctime = time.asctime()
    if original_number_of_operators is None:
       self.original_number_of_operators = number_of_operators

    self._has_updated_operators = False

  def show_summary(self, out = sys.stdout):
    """Summarize ncs_info_object"""
    print("NCS file:%s   Operators: %d" %(self.file_name,
      self.number_of_operators), file = out)
    if self.is_helical_symmetry:
      print("Helical symmetry is present", file = out)

  def has_updated_operators(self):
    """Return True if operators have been updated"""
    return self._has_updated_operators

  def update_number_of_operators(self, number_of_operators = None):
    """Update number of operators in ncs_info_object"""
    self.number_of_operators = number_of_operators
    self._has_updated_operators = True

  def update_is_helical_symmetry(self, is_helical_symmetry = None):
    """Update is_helical_symmetry in ncs_info_object"""
    self.is_helical_symmetry = is_helical_symmetry
    self._has_updated_operators = True


class map_info_object:
  """Holder for map file name, origin, all(), crystal_symmetry, b_sharpen,
      whether it is a map (alternative is mask), map_id, and id"""
  def __init__(self,
    file_name = None,
    origin = None,
    all = None,
    crystal_symmetry = None,
    is_map = None,
    map_id = None,
    b_sharpen = None,
    id = None,
    ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    import time
    self.init_asctime = time.asctime()

  def show_summary(self, out = sys.stdout):
    """Summarize map_info_object"""
    if self.is_map:
      print("Map file:%s" %(self.file_name), end = ' ', file = out)
    else:
      print("Mask file:%s" %(self.file_name), end = ' ', file = out)
    if self.id is not None:
      print("ID: %d" %(self.id), end = ' ', file = out)
    if self.b_sharpen is not None:
      print("B-sharpen: %7.2f" %(self.b_sharpen), end = ' ', file = out)
    if self.map_id is not None:
      print("Map ID: %s" %(self.map_id), file = out)
    else:
      print(file = out)
    if self.origin and self.all:
      print("   Origin: %d  %d  %d   Extent: %d  %d  %d" %(
       tuple(self.origin)+tuple(self.all)), file = out)
    if self.crystal_symmetry:
      print("   Map unit cell: %.1f  %.1f  %.1f    %.1f  %.1f  %.1f " %(
        self.crystal_symmetry.unit_cell().parameters()), file = out)

  def lower_upper_bounds(self):
    """Return lower and upper bounds of this map_info_object"""
    lower_bounds = self.origin
    upper_bounds = []
    for a, b in zip(self.origin, self.all):
      upper_bounds.append(a+b-1) # 2019-11-05 upper bound is na-1
    return list(self.origin), list(upper_bounds)

class info_object:
  """Holder for information about a map and its crystal symmetry, NCS,
    B-values, and segmentation"""
  def __init__(self,
      acc = None,
      ncs_obj = None,
      min_b = None,
      max_b = None,
      b_sharpen = None,  # b_sharpen applied to map
      ncs_group_list = None,
      origin_shift = None,
      crystal_symmetry = None, # after density_select
      original_crystal_symmetry = None, # before density_select
      full_crystal_symmetry = None, # from real_map object
      full_unit_cell_grid = None,  # from real_map object
      edited_volume_list = None,
      region_range_dict = None,
      selected_regions = None,
      ncs_related_regions = None,
      self_and_ncs_related_regions = None,
      map_files_written = None,
      bad_region_list = None,
      region_centroid_dict = None,
      original_id_from_id = None,
      remainder_id_dict = None,  # dict relating regions in a remainder object to
      params = None, # input params
      input_pdb_info = None,
      input_map_info = None,
      input_ncs_info = None,
      input_seq_info = None,
      shifted_pdb_info = None,
      shifted_map_info = None,
      shifted_ncs_info = None,
      shifted_used_ncs_info = None,
      n_residues = None,
      solvent_fraction = None,
      output_ncs_au_map_info = None,
      output_ncs_au_mask_info = None,
      output_ncs_au_pdb_info = None,
      output_box_map_info = None,
      output_box_mask_info = None,
      output_region_map_info_list = None,
      output_region_pdb_info_list = None,
      sharpening_info_obj = None,
      box_map_bounds_first = None,
      box_map_bounds_last = None,
      final_output_sharpened_map_file = None,
      box_map_ncs_au = None,
      box_map_ncs_au_crystal_symmetry = None,
    ):
    if not selected_regions: selected_regions = []
    if not ncs_related_regions: ncs_related_regions = []
    if not self_and_ncs_related_regions: self_and_ncs_related_regions = []
    if not map_files_written: map_files_written = []
    if not output_region_map_info_list: output_region_map_info_list = []
    if not output_region_pdb_info_list: output_region_pdb_info_list = []
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

    self.object_type = "segmentation_info"
    import time
    self.init_asctime = time.asctime()

  def set_box_map_ncs_au_map_data(self,
       box_map_ncs_au_map_data = None,
       box_mask_ncs_au_map_data = None,
       box_map_ncs_au_half_map_data_list = None,
       box_map_ncs_au_crystal_symmetry = None):
    """Set box_map_ncs information"""
    self.box_map_ncs_au_map_data = box_map_ncs_au_map_data.deep_copy()
    self.box_mask_ncs_au_map_data = box_mask_ncs_au_map_data.deep_copy()
    self.box_map_ncs_au_half_map_data_list = []
    for hm in box_map_ncs_au_half_map_data_list:
       self.box_map_ncs_au_half_map_data_list.append(hm.deep_copy())
    self.box_map_ncs_au_crystal_symmetry = box_map_ncs_au_crystal_symmetry
    if self.origin_shift and self.origin_shift !=  (0, 0, 0):
      self.box_map_ncs_au_map_data = self.shift_map_back(
        map_data = self.box_map_ncs_au_map_data,
        crystal_symmetry = self.box_map_ncs_au_crystal_symmetry,
        shift_cart = self.origin_shift)
      self.box_mask_ncs_au_map_data = self.shift_mask_back(
        mask_data = self.box_mask_ncs_au_map_data,
        crystal_symmetry = self.box_mask_ncs_au_crystal_symmetry,
        shift_cart = self.origin_shift)

      new_hm_list = []
      for hm in self.box_map_ncs_au_half_map_data_list:
        hm = self.shift_map_back(
          map_data = hm,
          crystal_symmetry = self.box_map_ncs_au_crystal_symmetry,
          shift_cart = self.origin_shift)
        new_hm_list.append(hm)
      self.box_map_ncs_au_half_map_data_list = new_hm_list

  def shift_map_back(self, map_data = None,
      crystal_symmetry = None, shift_cart = None):
    """Shift map back to original location"""
    from scitbx.matrix import col
    new_origin = self.origin_shift_grid_units(crystal_symmetry = crystal_symmetry,
      map_data = map_data, shift_cart = shift_cart, reverse = True)
    new_all = list(col(map_data.all())+col(new_origin))
    shifted_map_data = map_data.deep_copy()
    shifted_map_data.resize(flex.grid(new_origin, new_all))
    return shifted_map_data

  def origin_shift_grid_units(self, crystal_symmetry = None, map_data = None,
       shift_cart = None, reverse = False):
    """Get origin shift in grid units from shift_cart"""
    from scitbx.matrix import col
    cell = crystal_symmetry.unit_cell().parameters()[:3]
    origin_shift_grid = []
    for s, c, a in zip(shift_cart, cell, map_data.all()):
      if s<0:
        delta = -0.5
      else:
        delta = 0.5
      origin_shift_grid.append( int(delta+ a*s/c))
    if reverse:
      return list(-col(origin_shift_grid))
    else:
      return origin_shift_grid


  def is_segmentation_info_object(self):
    """Return True, this is a segmentation_info_object"""
    return True

  def set_params(self, params):
    """Set self.params from params"""
    self.params = deepcopy(params)

  def set_input_seq_info(self, file_name = None, sequence = None, n_residues = None):
    """Set input sequence information"""
    self.input_seq_info = seq_info_object(file_name = file_name,
       sequence = sequence,
       n_residues = n_residues)

  def set_input_pdb_info(self, file_name = None, n_residues = None):
    """Set input PDB information"""
    self.input_pdb_info = pdb_info_object(file_name = file_name,
     n_residues = n_residues)

  def set_input_ncs_info(self, file_name = None, number_of_operators = None):
    """Set input NCS information"""
    self.input_ncs_info = ncs_info_object(file_name = file_name,
      number_of_operators = number_of_operators)

  def update_ncs_info(self, number_of_operators = None, is_helical_symmetry = None,
      shifted = False):
    """Update NCS information"""
    if shifted:
      ncs_info = self.shifted_ncs_info
    else:
      ncs_info = self.input_ncs_info
    assert ncs_info
    if number_of_operators is not None:
      ncs_info.update_number_of_operators(
        number_of_operators = number_of_operators)
    if is_helical_symmetry is not None:
      ncs_info.update_is_helical_symmetry(
        is_helical_symmetry = is_helical_symmetry)

  def set_sharpening_info(self, sharpening_info_obj = None):
     """Set sharpening information """
     self.sharpening_info_obj = sharpening_info_obj

  def set_input_map_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None):
    """Create map_info object and use it to set input_map info"""
    self.input_map_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      is_map = True)

  def set_ncs_obj(self, ncs_obj = None):
    """Set NCS object"""
    self.ncs_obj = ncs_obj

  def set_origin_shift(self, origin_shift = None):
    """Set origin shift"""
    if not origin_shift: origin_shift = (0, 0, 0)
    self.origin_shift = tuple(origin_shift)

  def set_crystal_symmetry(self, crystal_symmetry):
    """Set crystal_symmetry"""
    self.crystal_symmetry = deepcopy(crystal_symmetry)

  def set_original_crystal_symmetry(self, crystal_symmetry):
    """Set original_crystal_symmetry"""
    self.original_crystal_symmetry = deepcopy(crystal_symmetry)

  def set_full_crystal_symmetry(self, crystal_symmetry):
    """Set full_crystal_symmetry"""
    self.full_crystal_symmetry = deepcopy(crystal_symmetry)

  def set_full_unit_cell_grid(self, unit_cell_grid):
    """Set full unit cell grid"""
    self.full_unit_cell_grid = deepcopy(unit_cell_grid)

  def set_box_map_bounds_first_last(self, box_map_bounds_first,
      box_map_bounds_last):
    """Set box map_bounds first and last"""
    self.box_map_bounds_first = box_map_bounds_first
    self.box_map_bounds_last = []
    for l in box_map_bounds_last:
      self.box_map_bounds_last.append(l+1)  # it is one bigger...

  def set_accessor(self, acc):
    """Set accessor"""
    self.acc = acc

  def set_shifted_map_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None, b_sharpen = None):
    """Create map_info object and use it to set shifted map_info"""
    self.shifted_map_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      b_sharpen = b_sharpen,
      is_map = True)

  def set_shifted_pdb_info(self, file_name = None, n_residues = None):
    """Set shifted_pdb info"""
    self.shifted_pdb_info = pdb_info_object(file_name = file_name,
     n_residues = n_residues)

  def set_shifted_ncs_info(self, file_name = None, number_of_operators = None,
       is_helical_symmetry = None):
    """Set shifted_ncs info"""
    self.shifted_ncs_info = ncs_info_object(file_name = file_name,
      number_of_operators = number_of_operators,
      is_helical_symmetry = is_helical_symmetry)

  def set_shifted_used_ncs_info(self, file_name = None, number_of_operators = None,
       is_helical_symmetry = None):
    """Create ncs_info object and use it to set shifted_used_ncs_info"""
    self.shifted_used_ncs_info = ncs_info_object(file_name = file_name,
      number_of_operators = number_of_operators,
      is_helical_symmetry = is_helical_symmetry)

  def set_solvent_fraction(self, solvent_fraction):
    """Set the solvent fraction"""
    self.solvent_fraction = solvent_fraction

  def set_n_residues(self, n_residues): # may not be the same as seq file
    """Set the number of residues"""
    self.n_residues = n_residues

  def set_output_ncs_au_map_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None):
    """Create map_info object and set output_ncs_au_map_info"""
    self.output_ncs_au_map_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      is_map = True)

  def set_output_ncs_au_mask_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None):
    """Create map_info object and set output_ncs_au_mask_info"""
    self.output_ncs_au_mask_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      is_map = False)

  def set_output_ncs_au_pdb_info(self, file_name = None, n_residues = None):
    """Create pdb_info object and use it to set output_ncs_au_pdb_info"""
    self.output_ncs_au_pdb_info = pdb_info_object(file_name = file_name,
     n_residues = n_residues)

  def set_output_box_map_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None):
    """Create map_info object and set output_ncs_box_map_info"""
    self.output_box_map_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      is_map = True)

  def set_output_box_mask_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None):
    """Create map_info object and set output_ncs_box_mask_info"""
    self.output_box_mask_info = map_info_object(file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      is_map = False)

  def add_output_region_map_info(self, file_name = None, crystal_symmetry = None,
    origin = None, all = None, map_id = None):
    """Create map_info object and set output_region_map_info"""
    self.output_region_map_info_list.append(map_info_object(
      file_name = file_name,
      crystal_symmetry = crystal_symmetry,
      origin = origin,
      all = all,
      id = len(self.output_region_map_info_list)+1,
      map_id = map_id,
      is_map = True)
     )

  def add_output_region_pdb_info(self, file_name = None, n_residues = None):
    """Create pdb_info_object and append to output_region_pdb_info_list"""
    self.output_region_pdb_info_list.append(pdb_info_object(
      file_name = file_name,
      n_residues = n_residues)
     )


  def show_summary(self, out = sys.stdout):
    """Summarize this info_object"""
    print("\n ==========  Summary of %s:  ======== \n" %(self.object_type), file = out)
    print("Created: %s" %(self.init_asctime), file = out)
    print("\nInput files used:\n", file = out)
    if self.input_map_info:
      self.input_map_info.show_summary(out = out)
    if self.input_pdb_info:
      self.input_pdb_info.show_summary(out = out)
    if self.input_ncs_info:
      self.input_ncs_info.show_summary(out = out)
    if self.input_seq_info:
      self.input_seq_info.show_summary(out = out)

    print(file = out)

    if self.crystal_symmetry:
      print("Working unit cell: %.1f  %.1f  %.1f    %.1f  %.1f  %.1f " %(
        self.crystal_symmetry.unit_cell().parameters()), file = out)

    if self.n_residues:
      print("Estimated total number of residues: %d" %(self.n_residues), file = out)

    if self.solvent_fraction:
      print("Estimated solvent fraction: %5.3f" %(self.solvent_fraction), file = out)

    if self.origin_shift and self.origin_shift !=  (0, 0, 0):
      print("\nOrigin offset applied: %.1f  %.1f  %.1f" %(self.origin_shift), file = out)
    else:
      print("\nNo origin offset applied", file = out)

    if self.shifted_map_info:
      print("\nShifted/sharpened map, pdb and ncs files created "+\
         "(after origin offset):\n", file = out)
      if self.shifted_map_info:
        self.shifted_map_info.show_summary(out = out)
      if self.shifted_pdb_info:
        self.shifted_pdb_info.show_summary(out = out)
      if self.shifted_ncs_info:
        self.shifted_ncs_info.show_summary(out = out)

    if self.output_ncs_au_pdb_info:
      print("\nOutput PDB file with dummy atoms representing the NCS AU:", file = out)
      self.output_ncs_au_pdb_info.show_summary(out = out)

    if self.output_ncs_au_mask_info or self.output_ncs_au_map_info:
      print("\nOutput map files showing just the NCS AU (same size", end = ' ', file = out)
      if self.origin_shift and self.origin_shift !=  (0, 0, 0):
        print("\nand location as shifted map files:\n", file = out)
      else:
        print("\nand location as input map:\n", file = out)

      if self.output_ncs_au_mask_info:
        self.output_ncs_au_mask_info.show_summary(out = out)
      if self.output_ncs_au_map_info:
        self.output_ncs_au_map_info.show_summary(out = out)

    if self.output_box_mask_info or self.output_box_map_info:
      print("\nOutput cut-out map files trimmed to contain just "+\
        "the \nNCS AU (superimposed on", end = ' ', file = out)
      if self.origin_shift and self.origin_shift !=  (0, 0, 0):
        print("shifted map files, note origin offset):\n", file = out)
      else:
        print("input map, note origin offset):\n", file = out)

      if self.output_box_mask_info:
        self.output_box_mask_info.show_summary(out = out)
      if self.output_box_map_info:
        self.output_box_map_info.show_summary(out = out)

    if self.output_region_pdb_info_list:
      print("\nOutput PDB files representing one region of connected"+\
        " density.\nThese are useful for marking where to look in cut-out map"+\
        " files.", file = out)
      for output_region_pdb_info in self.output_region_pdb_info_list:
        output_region_pdb_info.show_summary(out = out)

    if self.output_region_map_info_list:
      print("\nOutput cut-out map files trimmed to contain just "+\
        "one region of \nconnected density (superimposed on", end = ' ', file = out)
      if self.origin_shift and self.origin_shift !=  (0, 0, 0):
        print("shifted map files, note origin offset):\n", file = out)
      else:
        print(" input map, note origin offset):\n", file = out)
      for output_region_map_info in self.output_region_map_info_list:
        output_region_map_info.show_summary(out = out)

    print("\n"+50*"="+"\n", file = out)

class make_ccp4_map: # just a holder so map_to_structure_factors will run
  """Holder to run map_to_structure_factors.  Deprecated,
  replace with map_manager."""
  def __init__(self, map = None, unit_cell = None):
    self.data = map
    self.unit_cell_parameters = unit_cell.parameters()
    self.space_group_number = 1
    self.unit_cell_grid = map.all()

  def unit_cell(self):
    """Return the unit_cell of this map"""
    return self.crystal_symmetry().unit_cell()

  def unit_cell_crystal_symmetry(self):
    """Return the unit_cell_crystal_symmetry of this map.
       Same as crystal_symmetry"""
    return self.crystal_symmetry()

  def map_data(self):
    """Return the map_data for this map"""
    return self.data

  def crystal_symmetry(self):
    """Return the crystal_symmetry for this map"""
    return crystal.symmetry(self.unit_cell_parameters,
        self.space_group_number)


class b_vs_region_info:
  """Holder for b_iso values in regions"""
  def __init__(self):
    self.b_iso = 0.
    self.b_vs_region_dict = {}
    self.sa_sum_v_vs_region_dict = {}
    self.sa_nn_vs_region_dict = {}
    self.sa_ratio_b_vs_region_dict = {}

class box_sharpening_info:
  """Holder for information about map sharpening for a box"""
  def __init__(self, tracking_data = None,
      crystal_symmetry = None,
      solvent_fraction = None,
      b_iso = None,
      resolution = None,
      d_min_ratio = None,
      scale_max = None,
      lower_bounds = None,
      upper_bounds = None,
      wrapping = None,
      n_real = None,
      n_buffer = None,
      map_data = None,
      smoothing_radius = None,
      smoothed_box_mask_data = None,
      original_box_map_data = None,
       ):
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    del self.tracking_data # do not save it
    if tracking_data:
      self.crystal_symmetry = tracking_data.crystal_symmetry
      self.solvent_fraction = tracking_data.solvent_fraction
      self.wrapping = tracking_data.params.crystal_info.use_sg_symmetry

  def get_gaussian_weighting(self, out = sys.stdout):
    """Return a gaussian function centered on center of the map, fall-off
     based on smoothing_radius"""

    # Calculate weight map, max near location of centers_ncs_cart
    # U = rmsd**2
    # (b_eff = 8*3.14159**2*U)
    #  rmsd is at least distance between centers, not too much bigger than
    #  unit cell size, typically 10-20 A,
    print("\nFall-off of local weight is 1/%6.1f A\n" %(
       self.smoothing_radius), file = out)
    u = self.smoothing_radius**2

    from cctbx import xray
    xrs, scatterers = set_up_xrs(crystal_symmetry = self.crystal_symmetry)
    unit_cell = self.crystal_symmetry.unit_cell()
    for xyz_fract in [(0.5, 0.5, 0.5, )]:
      scatterers.append( xray.scatterer(scattering_type = "H", label = "H",
        site = xyz_fract, u = u, occupancy = 1.0))
    xrs = xray.structure(xrs, scatterers = scatterers)
    f_array, phases = get_f_phases_from_map(map_data = self.map_data,
       crystal_symmetry = self.crystal_symmetry,
       d_min = self.resolution,
       scale_max = self.scale_max,
       d_min_ratio = self.d_min_ratio,
       get_remove_aniso_object = False, # don't need it
       out = out)

    weight_f_array = f_array.structure_factors_from_scatterers(
      algorithm = 'direct',
      xray_structure = xrs).f_calc()

    weight_map = get_map_from_map_coeffs(map_coeffs = weight_f_array,
      crystal_symmetry = self.crystal_symmetry, n_real = self.map_data.all())
    min_value = weight_map.as_1d().min_max_mean().min
    weight_map = weight_map-min_value # all positive or zero
    max_value = weight_map.as_1d().min_max_mean().max
    weight_map = weight_map/max(1.e-10, max_value)  # normalize; max = 1 now
    min_value = 1.e-10  # just a small value for all distances far from center
    s = (weight_map <min_value )  # make extra sure every point is above this
    weight_map = weight_map.set_selected(s, min_value)
    return weight_map


  def remove_buffer_from_bounds(self, minimum = 1):
    """Back off by n_buffer in each direction, leave at
      least minimum grid on either side of center"""

    adjusted_lower_bounds, adjusted_upper_bounds = [], []
    delta_lower_bounds, delta_upper_bounds = [], []
    for lb, ub in zip(self.lower_bounds, self.upper_bounds):
      sum = lb+ub
      if sum >= 0:
        mid = (1+sum)//2
      else:
        mid = (-1+sum)//2
      alb = min(mid-minimum, lb+self.n_buffer)
      aub = max(mid+minimum, ub-self.n_buffer)
      adjusted_lower_bounds.append(alb)
      adjusted_upper_bounds.append(aub)
      delta_lower_bounds.append(alb-lb)
      delta_upper_bounds.append(aub-ub)
    return adjusted_lower_bounds, adjusted_upper_bounds, \
      delta_lower_bounds, delta_upper_bounds


  def merge_into_overall_map(self, overall_map = None):
    """Smoothly fill out edges of the small map with overall_map"""

    assert self.smoothed_box_mask_data is not None
    assert self.original_box_map_data is not None

    self.map_data =  (self.map_data * self.smoothed_box_mask_data) + \
       (self.original_box_map_data * (1-self.smoothed_box_mask_data))

  def remove_buffer(self, out = sys.stdout):
    """Remove the buffer from this box"""

    new_lower_bounds, new_upper_bounds, delta_lower, delta_upper = \
       self.remove_buffer_from_bounds()

    cut_out_lower_bounds = []
    cut_out_upper_bounds = []
    for o, a, dlb, dub in zip(self.map_data.origin(), self.map_data.all(),
      delta_lower, delta_upper):
     cut_out_lower_bounds.append(o+dlb)
     cut_out_upper_bounds.append(a+dub-1)
    self.map_data, self.crystal_symmetry, \
       self.smoothed_box_mask_data, self.original_box_map_data = cut_out_map(
          map_data = self.map_data,
          crystal_symmetry = self.crystal_symmetry,
          soft_mask = False,
          resolution = self.resolution,
          shift_origin = True,
          min_point = cut_out_lower_bounds,
          max_point = cut_out_upper_bounds, out = out)
    self.lower_bounds = new_lower_bounds
    self.upper_bounds = new_upper_bounds

class sharpening_info:
  """Holder for sharpening information about a map"""
  def __init__(self,
      tracking_data = None,
      crystal_symmetry = None,
      is_crystal = None,
      sharpening_method = None,
      solvent_fraction = None,
      n_residues = None,
      ncs_copies = None,
      ncs_file = None,
      seq_file = None,
      sequence = None,
      n_real = None,
      region_weight = None,
      n_bins = None,
      eps = None,
      d_min = None,
      d_min_ratio = None,
      scale_max = None,
      input_d_cut = None,
      b_blur_hires = None,
      rmsd = None,
      rmsd_resolution_factor = None,
      k_sol = None,
      b_sol = None,
      fraction_complete = None,
      wrapping = None,
      sharpening_target = None,
      residual_target = None,
      fraction_occupied = None,
      nproc = None,
      multiprocessing = None,
      queue_run_command = None,
      resolution = None, # changed from d_cut
      resolution_dependent_b = None,  # linear sharpening
      normalize_amplitudes_in_resdep = None,  # linear sharpening
      b_sharpen = None,
      b_iso = None,  # expected B_iso after applying b_sharpen
      k_sharpen = None,
      optimize_b_blur_hires = None,
      iterate = None,
      optimize_d_cut = None,
      kurtosis = None,
      adjusted_sa = None,
      sa_ratio = None,
      normalized_regions = None,
      score = None,
      input_weight_map_pickle_file = None,
      output_weight_map_pickle_file = None,
      read_sharpened_maps = None,
      write_sharpened_maps = None,
      select_sharpened_map = None,
      output_directory = None,
      smoothing_radius = None,
      local_sharpening = None,
      local_aniso_in_local_sharpening = None,
      overall_before_local = None,
      use_local_aniso = None,
      original_aniso_obj = None,
      auto_sharpen = None,
      box_in_auto_sharpen = None,
      density_select_in_auto_sharpen = None,
      density_select_threshold_in_auto_sharpen = None,
      use_weak_density = None,
      discard_if_worse = None,
      max_box_fraction = None,
      cc_cut = None,
      max_cc_for_rescale = None,
      scale_using_last = None,
      density_select_max_box_fraction = None,
      mask_atoms = None,
      mask_atoms_atom_radius = None,
      value_outside_atoms = None,
      soft_mask = None,
      allow_box_if_b_iso_set = None,
      search_b_min = None,
      search_b_max = None,
      search_b_n = None,
      adjust_region_weight = None,
      region_weight_method = None,
      region_weight_factor = None,
      region_weight_buffer = None,
      region_weight_default = None,
      target_b_iso_ratio = None,
      signal_min = None,
      target_b_iso_model_scale = None,
      box_sharpening_info_obj = None,
      chain_type = None,
      target_scale_factors = None,
      remove_aniso = None,
      d_min_list = None,
      verbose = None,
      resolve_size = None,
      pdb_hierarchy = None,  # XXX just used to set params
      local_solvent_fraction = None,
      wang_radius = None,
      buffer_radius = None,
      pseudo_likelihood = None,
      preliminary_sharpening_done = False,
      adjusted_path_length = None,
        ):

    from libtbx import adopt_init_args
    adopt_init_args(self, locals())
    del self.tracking_data  # don't need it as part of the object
    del self.box_sharpening_info_obj# don't need it as part of the object

    if tracking_data:  # use tracking data information
      self.update_with_tracking_data(tracking_data = tracking_data)

    if box_sharpening_info_obj: # update information
      self.update_with_box_sharpening_info(
         box_sharpening_info_obj = box_sharpening_info_obj)

    if self.resolution_dependent_b is None:
      self.resolution_dependent_b = [0, 0, 0]

    if self.target_scale_factors and \
        self.sharpening_method!= 'model_sharpening' \
        and self.sharpening_method!= 'half_map_sharpening':
      assert self.sharpening_method is None # XXX may want to print out error
      self.sharpening_method = 'model_sharpening'

    if self.sharpening_method == 'b_iso' and self.k_sharpen is not None:
      self.k_sharpen = None

    if pdb_hierarchy:
        self.sharpening_method = 'model_sharpening'
        self.box_in_auto_sharpen = True
        self.density_select_in_auto_sharpen = False
        self.sharpening_target = 'model'

  def get_d_cut(self):
    """Get a reasonable valuen of d_cut from input_d_cut and resolution"""
    if self.input_d_cut is not None:
       return self.input_d_cut
    else:
       return self.resolution

  def get_target_b_iso(self):
    """Get a reasonable target_b_iso
       from self.target_b_iso_ratio*self.resolution**2"""
    if self.target_b_iso_ratio is None:
      return None
    if self.resolution is None:
      return None
    return self.target_b_iso_ratio*self.resolution**2

  def set_resolution_dependent_b(self,
      resolution_dependent_b = None,
      sharpening_method = 'resolution_dependent'):
    """Set the value of resolution_dependent_b and sharpening_method"""
    if resolution_dependent_b:
      self.resolution_dependent_b = resolution_dependent_b
    if sharpening_method:
      self.sharpening_method = sharpening_method

  def sharpening_is_defined(self):
    """Return True if sharpening is defined"""
    if self.sharpening_method is None:
      return False

    if self.target_scale_factors:
      return True

    if self.sharpening_method == 'target_b_iso_to_d_cut':
      return True

    if self.b_iso is not None or \
       self.b_sharpen is not None or \
       (self.resolution_dependent_b is not None and
        self.resolution_dependent_b!= [0, 0, 0]):
      return True

    return False

  def update_with_box_sharpening_info(self, box_sharpening_info_obj = None):
      """Update stored information using info in box_sharpening_info_obj"""
      if not box_sharpening_info_obj:
        return self
      self.crystal_symmetry = box_sharpening_info_obj.crystal_symmetry
      self.solvent_fraction = box_sharpening_info_obj.solvent_fraction
      self.wrapping = box_sharpening_info_obj.wrapping
      self.n_real = box_sharpening_info_obj.n_real
      return self

  def update_with_tracking_data(self, tracking_data = None):
      """Update stored information using info in tracking_data"""
      self.update_with_params(params = tracking_data.params,
         crystal_symmetry = tracking_data.crystal_symmetry,
         solvent_fraction = tracking_data.solvent_fraction,
         n_residues = tracking_data.n_residues,
         ncs_copies = tracking_data.input_ncs_info.number_of_operators)
      return self

  def update_with_params(self, params = None,
     crystal_symmetry = None,
     is_crystal = None,
     solvent_fraction = None,
     auto_sharpen = None,
     sharpening_method = None,
     pdb_hierarchy = None,
     half_map_data_list = None,
     n_residues = None, ncs_copies = None):
      """Update stored information using info in params"""
      self.crystal_symmetry = crystal_symmetry
      self.is_crystal = is_crystal
      self.solvent_fraction = solvent_fraction
      self.auto_sharpen = auto_sharpen
      self.n_residues = n_residues
      self.ncs_copies = ncs_copies
      self.seq_file = params.input_files.seq_file
      self.chain_type = params.crystal_info.chain_type
      self.verbose = params.control.verbose
      self.resolve_size = params.control.resolve_size
      self.multiprocessing = params.control.multiprocessing
      self.nproc = params.control.nproc
      self.queue_run_command = params.control.queue_run_command

      self.wrapping = params.crystal_info.use_sg_symmetry
      self.fraction_occupied = params.map_modification.fraction_occupied
      self.sa_percent = params.map_modification.sa_percent
      self.region_weight = params.map_modification.region_weight
      self.max_regions_to_test = params.map_modification.max_regions_to_test
      self.regions_to_keep = params.map_modification.regions_to_keep
      self.d_min_ratio = params.map_modification.d_min_ratio
      self.scale_max = params.map_modification.scale_max

      self.input_d_cut = params.map_modification.input_d_cut
      self.b_blur_hires = params.map_modification.b_blur_hires
      self.rmsd = params.map_modification.rmsd
      self.rmsd_resolution_factor = params.map_modification.rmsd_resolution_factor
      self.k_sol = params.map_modification.k_sol
      self.b_sol = params.map_modification.b_sol
      self.fraction_complete = params.map_modification.fraction_complete
      self.resolution = params.crystal_info.resolution  # changed from d_cut
      #  NOTE:
      #  resolution = X-ray resolution or nominal resolution of cryoEM map
      #  high-res cutoff of reflections is d_min*d_min_ratio
      self.buffer_radius = params.crystal_info.buffer_radius
      self.wang_radius = params.crystal_info.wang_radius
      self.pseudo_likelihood = params.crystal_info.pseudo_likelihood

      self.max_box_fraction = params.map_modification.max_box_fraction
      self.cc_cut = params.map_modification.cc_cut
      self.max_cc_for_rescale = params.map_modification.max_cc_for_rescale
      self.scale_using_last = params.map_modification.scale_using_last
      self.density_select_max_box_fraction = params.map_modification.density_select_max_box_fraction
      self.mask_atoms = params.map_modification.mask_atoms
      self.mask_atoms_atom_radius = params.map_modification.mask_atoms_atom_radius
      self.value_outside_atoms = params.map_modification.value_outside_atoms
      self.soft_mask = params.map_modification.soft_mask
      self.allow_box_if_b_iso_set = params.map_modification.allow_box_if_b_iso_set
      self.k_sharpen = params.map_modification.k_sharpen
      self.optimize_b_blur_hires = params.map_modification.optimize_b_blur_hires
      self.iterate = params.map_modification.iterate
      self.optimize_d_cut = params.map_modification.optimize_d_cut
      self.sharpening_target = params.map_modification.sharpening_target
      self.residual_target = params.map_modification.residual_target
      self.eps = params.map_modification.eps
      self.n_bins = params.map_modification.n_bins
      self.input_weight_map_pickle_file = params.input_files.input_weight_map_pickle_file
      self.output_weight_map_pickle_file = params.output_files.output_weight_map_pickle_file
      self.read_sharpened_maps = params.map_modification.read_sharpened_maps
      self.write_sharpened_maps = params.map_modification.write_sharpened_maps
      self.select_sharpened_map = params.map_modification.select_sharpened_map
      self.output_directory = params.output_files.output_directory
      self.smoothing_radius = params.map_modification.smoothing_radius
      self.local_sharpening = params.map_modification.local_sharpening
      self.local_aniso_in_local_sharpening = \
         params.map_modification.local_aniso_in_local_sharpening
      self.overall_before_local = \
         params.map_modification.overall_before_local
      self.box_in_auto_sharpen = params.map_modification.box_in_auto_sharpen
      self.density_select_in_auto_sharpen = params.map_modification.density_select_in_auto_sharpen
      self.density_select_threshold_in_auto_sharpen = params.map_modification.density_select_threshold_in_auto_sharpen
      self.use_weak_density = params.map_modification.use_weak_density
      self.discard_if_worse = params.map_modification.discard_if_worse
      self.box_center = params.map_modification.box_center
      self.box_size = params.map_modification.box_size
      self.target_n_overlap = params.map_modification.target_n_overlap
      self.restrict_map_size = params.map_modification.restrict_map_size
      self.remove_aniso = params.map_modification.remove_aniso
      self.min_ratio_of_ncs_copy_to_first = \
         params.segmentation.min_ratio_of_ncs_copy_to_first
      self.max_ratio_to_target = params.segmentation.max_ratio_to_target
      self.min_ratio_to_target = params.segmentation.min_ratio_to_target
      self.residues_per_region = params.segmentation.residues_per_region
      self.mask_padding_fraction = \
         params.segmentation.mask_padding_fraction
      self.fraction_of_max_mask_threshold = \
         params.segmentation.fraction_of_max_mask_threshold
      self.cell_cutoff_for_solvent_from_mask = \
         params.segmentation.cell_cutoff_for_solvent_from_mask
      self.starting_density_threshold = \
         params.segmentation.starting_density_threshold
      self.density_threshold = params.segmentation.density_threshold
      self.min_ratio = params.segmentation.min_ratio
      self.min_volume = params.segmentation.min_volume
      self.search_b_min = params.map_modification.search_b_min
      self.search_b_max = params.map_modification.search_b_max
      self.search_b_n = params.map_modification.search_b_n
      self.adjust_region_weight = params.map_modification.adjust_region_weight
      self.region_weight_method = params.map_modification.region_weight_method
      self.region_weight_factor = params.map_modification.region_weight_factor
      self.region_weight_buffer = params.map_modification.region_weight_buffer
      self.region_weight_default = params.map_modification.region_weight_default
      self.target_b_iso_ratio = params.map_modification.target_b_iso_ratio
      self.signal_min = params.map_modification.signal_min
      self.target_b_iso_model_scale = params.map_modification.target_b_iso_model_scale

      if sharpening_method is not None:
        self.sharpening_method = sharpening_method

      if not self.sharpening_method and \
         len(params.map_modification.auto_sharpen_methods) == 1:
        self.sharpening_method = params.map_modification.auto_sharpen_methods[0]

      if half_map_data_list or self.sharpening_method == 'half_map_sharpening':
        self.sharpening_method = 'half_map_sharpening'
        self.sharpening_target = 'half_map'

      elif pdb_hierarchy or self.sharpening_method == 'model_sharpening':
        self.sharpening_method = 'model_sharpening'
        self.box_in_auto_sharpen = True
        self.density_select_in_auto_sharpen = False
        self.sharpening_target = 'model'

      elif params.map_modification.b_iso is not None or \
          params.map_modification.b_sharpen is not None:
        if self.sharpening_method is None:
          raise Sorry("b_iso is not set")
        # if sharpening values are specified, set them
        if params.map_modification.b_iso is not None:
          self.b_iso = params.map_modification.b_iso # but we need b_sharpen
        elif params.map_modification.b_sharpen is not None:
          self.b_sharpen = params.map_modification.b_sharpen
      elif (params.map_modification.resolution_dependent_b is not None
        and params.map_modification.resolution_dependent_b!= [0, 0, 0]):
        self.sharpening_method = 'resolution_dependent'
        self.resolution_dependent_b = \
            params.map_modification.resolution_dependent_b

      if self.sharpening_method == 'b_iso' and self.k_sharpen is not None:
        self.k_sharpen = None
      return self

  def show_summary(self, verbose = False, list_scale_factors = True,
      out = sys.stdout):
    """Summarize this sharpening_info object"""
    method_summary_dict = {
       'b_iso':"Overall b_iso sharpening",
       'b_iso_to_d_cut':"b_iso sharpening to high_resolution cutoff",
       'resolution_dependent':"Resolution-dependent sharpening",
       'model_sharpening':"Model sharpening",
       'half_map_sharpening':"Half-map sharpening",
       'no_sharpening':"No sharpening",
       None:"No sharpening",
        }

    target_summary_dict = {
       'adjusted_sa':"Adjusted surface area",
       'adjusted_path_length':"Adjusted path length",
       'kurtosis':"Map kurtosis",
       'model':"Map-model CC",
      }
    print("\nSummary of sharpening:\n", file = out)

    if not hasattr(self, 'sharpening_method'):
      return # backwards compatibility

    print("Sharpening method used:         %s\n" %(
       method_summary_dict.get(self.sharpening_method)), file = out)

    if self.sharpening_method == "b_iso":
      if self.b_sharpen is not None:
        print("Overall b_sharpen applied:      %7.2f A**2" %(
          self.b_sharpen), file = out)
      if self.b_iso is not None:
        print("Final b_iso obtained:           %7.2f A**2" %(self.b_iso), file = out)
    elif self.sharpening_method == "b_iso_to_d_cut":
      if self.b_sharpen is not None:
        print("Overall b_sharpen applied:      %7.2f A**2" %(
          self.b_sharpen), file = out)
      if self.b_iso is not None:
        print("Final b_iso obtained:           %7.2f A**2" %(self.b_iso), file = out)
      if self.input_d_cut:
        print("High-resolution cutoff:         %7.2f A" %(self.input_d_cut), file = out)
      else:
        print("High-resolution cutoff:         %7.2f A" %(self.resolution), file = out)
    elif self.sharpening_method == "resolution_dependent":
      print("Resolution-dependent b values (%7.2f, %7.2f, %7.2f)\n" %(
        tuple(self.resolution_dependent_b)), file = out)

      print("Effective b_iso vs resolution obtained:", file = out)
      from cctbx.maptbx.refine_sharpening import get_effective_b_values
      d_min_values, b_values = get_effective_b_values(
        d_min_ratio = self.d_min_ratio,
         resolution_dependent_b = self.resolution_dependent_b,
         resolution = self.resolution)
      print("                                Resolution  Effective B-iso", file = out)
      print("                                    (A)         (A**2)", file = out)
      for dd, b in zip(d_min_values, b_values):
        print("                                 %7.1f       %7.2f " %(
         dd, b), file = out)

    elif self.sharpening_method == "model_sharpening":
      print("Resolution-dependent model sharpening", file = out)
      if self.d_min_list and self.target_scale_factors and list_scale_factors:
        print("Scale vs resolution:", file = out)
        for d_min, sc in zip(
          self.d_min_list,
          self.target_scale_factors):
          print("Dmin: %7.2f  Scale: %9.6f" %(d_min, sc), file = out)

    elif self.sharpening_method == "half_map_sharpening":
      print("Resolution-dependent half-map sharpening", file = out)
      if self.d_min_list and self.target_scale_factors and list_scale_factors:
        print("Scale vs resolution:", file = out)
        for d_min, sc in zip(
          self.d_min_list,
          self.target_scale_factors):
          print("Dmin: %7.2f  Scale: %9.6f" %(d_min, sc), file = out)

    if self.sharpening_method in ["b_iso_to_d_cut"] and \
      self.k_sharpen and self.resolution:
        print("Transition from sharpening"+\
        " to not sharpening (k_sharpen):%7.2f " %(self.k_sharpen), file = out)

    print("\nSharpening target used:         %s" %(
       target_summary_dict.get(self.sharpening_target)), file = out)
    if self.adjusted_sa is not None:
      print("Final adjusted map surface area:  %7.2f" %(self.adjusted_sa), file = out)
    if self.kurtosis is not None:
      print("Final map kurtosis:               %7.2f" %(self.kurtosis), file = out)
    if hasattr(self, 'adjusted_path_length') and \
         self.adjusted_path_length is not None:
      print("Final adjusted path length:        %7.2f A" %(
         self.adjusted_path_length), file = out)

    print(file = out)

    if verbose:
      for x in dir(self):
        if x.startswith("__"): continue
        if type(getattr(self, x)) in [type('a'), type(1), type(1.), type([]),
          type((1, 2, ))]:
          print("%s : %s" %(x, getattr(self, x)), file = out)

  def get_effective_b_iso(self, map_data = None, out = sys.stdout):
    map_coeffs_ra, map_coeffs, f_array, phases = effective_b_iso(
      map_data = map_data,
      resolution = self.resolution,
      d_min_ratio = self.d_min_ratio,
      scale_max = self.scale_max,
      crystal_symmetry = self.crystal_symmetry,
       out = out)
    return map_coeffs_ra.b_iso

  def sharpen_and_score_map(self, map_data = None, set_b_iso = False, out = sys.stdout):
    if self.n_real is None: # need to get it
      self.n_real = map_data.all()
    map_and_b = sharpen_map_with_si(
        sharpening_info_obj = self,
        map_data = map_data,
        resolution = self.resolution, out = out)
    self.map_data = map_and_b.map_data
    if set_b_iso:
      self.b_iso = map_and_b.final_b_iso
    score_map(map_data = self.map_data,
        sharpening_info_obj = self,
        out = null_out())
    return self

  def show_score(self, out = sys.stdout):
    print("Adjusted surface area: %7.3f  Kurtosis: %7.3f  Score: %7.3f\n" %(
       self.adjusted_sa, self.kurtosis, self.score), file = out)

  def is_target_b_iso_to_d_cut(self):
    if self.sharpening_method == 'target_b_iso_to_d_cut':
       return True
    else:
       return False

  def is_b_iso_sharpening(self):
    if self.is_resolution_dependent_sharpening():
       return False
    if self.is_model_sharpening():
       return False
    if self.is_half_map_sharpening():
       return False
    if self.is_external_map_sharpening():
       return False
    return True

  def is_resolution_dependent_sharpening(self):
    if self.sharpening_method == 'resolution_dependent':
       return True
    else:
       return False

  def is_model_sharpening(self):
    if self.sharpening_method == 'model_sharpening':
       return True
    else:
       return False

  def is_external_map_sharpening(self):
    if self.sharpening_method == 'external_map_sharpening':
       return True
    else:
       return False

  def is_half_map_sharpening(self):
    if self.sharpening_method == 'half_map_sharpening':
       return True
    else:
       return False

  def as_map_coeffs(self, out = sys.stdout):
    map_data = getattr(self, 'map_data', None)
    if map_data:
      map_coeffs, dummy = get_f_phases_from_map(map_data = self.map_data,
       crystal_symmetry = self.crystal_symmetry,
       d_min = self.resolution,
       d_min_ratio = self.d_min_ratio,
       scale_max = self.scale_max,
       return_as_map_coeffs = True,
       out = out)
      return map_coeffs
    else:
      return None

  def as_map_data(self):
    return getattr(self, 'map_data', None)


class ncs_group_object:
  """Holder for NCS information for one group of NCS-related regions in map"""
  def __init__(self,
      ncs_obj = None,
      ncs_ops_used = None,
      ncs_group_list = None,
      edited_mask = None,
      crystal_symmetry = None,
      max_cell_dim = None,
      origin_shift = None,
      edited_volume_list = None,
      region_range_dict = None,
      selected_regions = None,
      ncs_related_regions = None,
      self_and_ncs_related_regions = None,
      equiv_dict = None,
      map_files_written = None,
      bad_region_list = None,
      region_centroid_dict = None,
      region_scattered_points_dict = None,
      shared_group_dict = None,
      co = None,
      min_b = None,
      max_b = None,
      original_id_from_id = None,
      remainder_id_dict = None,  # dict relating regions in a remainder object to
                               #  those in the original map
         ):
    if not selected_regions: selected_regions = []
    if not ncs_related_regions: ncs_related_regions = []
    if not self_and_ncs_related_regions: self_and_ncs_related_regions = []
    if not map_files_written: map_files_written = []
    if not max_cell_dim: max_cell_dim = 0.
    from libtbx import adopt_init_args
    adopt_init_args(self, locals())

    if self.crystal_symmetry and not self.max_cell_dim:
      self.max_cell_dim = 0.
      for x in self.crystal_symmetry.unit_cell().parameters()[:3]:
        self.max_cell_dim = max(max_cell_dim, x)

  def as_info_object(self):
    return info_object(
      ncs_obj = self.ncs_obj,
      max_b = self.max_b,
      min_b = self.min_b,
      ncs_group_list = self.ncs_group_list,
      origin_shift = self.origin_shift,
      edited_volume_list = self.edited_volume_list,
      region_range_dict = self.region_range_dict,
      selected_regions = self.selected_regions,
      ncs_related_regions = self.ncs_related_regions,
      self_and_ncs_related_regions = self.self_and_ncs_related_regions,
      bad_region_list = self.bad_region_list,
      region_centroid_dict = self.region_centroid_dict,
      original_id_from_id = self.original_id_from_id,
      map_files_written = self.map_files_written,
     )

  def set_ncs_ops_used(self, ncs_ops_used):
    self.ncs_ops_used = deepcopy(ncs_ops_used)

  def set_selected_regions(self, selected_regions):
    self.selected_regions = deepcopy(selected_regions)

  def set_ncs_related_regions(self, ncs_related_regions):
    self.ncs_related_regions = deepcopy(ncs_related_regions)

  def set_self_and_ncs_related_regions(self, self_and_ncs_related_regions):
    self.self_and_ncs_related_regions = deepcopy(self_and_ncs_related_regions)

  def set_map_files_written(self, map_files_written):
    self.map_files_written = deepcopy(map_files_written)

def zero_if_none(x):
  """Return 0 if x is None"""
  if not x:
    return 0
  else:
    return x
def scale_map(map, scale_rms = 1.0, out = sys.stdout):
    """Scale map to rms value of scale_rms. Does not subtract mean"""
    sd = map.as_double().as_1d().sample_standard_deviation()
    if (sd > 1.e-10):
      scale = scale_rms/sd
      if 0: print("Scaling map by %7.3f to set SD = 1" %(scale), file = out)
      map = map*scale
    else:
      print("Cannot scale map...all zeros", file = out)
    return map

def scale_map_coeffs(map_coeffs, scale_max = None, out = sys.stdout):
  """Scale map_coeffs to yield maximum value of scale_max"""
  f_array, phases = map_coeffs_as_fp_phi(map_coeffs)
  max_value = f_array.data().min_max_mean().max
  if scale_max and max_value is not None:
    scale = scale_max/max(1.e-10, max_value)
  else:
    scale = 1.0
  if 0:
    print("Scaling map_coeffs by %9.3f to yield maximum of %7.0f" %(
     scale, scale_max), file = out)
  return f_array.array(data = f_array.data()*scale
       ).phase_transfer(phase_source = phases, deg = True)


def get_map_object(file_name = None, must_allow_sharpening = None,
      get_map_labels = None, out = sys.stdout):
  """Read a ccp4 map file and return sg, cell and map objects"""
  # 2012-01-16
  if not os.path.isfile(file_name):
    raise Sorry("The map file %s is missing..." %(file_name))
  map_labels = None
  if file_name.endswith(".xplor"):
    raise Sorry("Unable to read xplor maps with segment_and_split_map")
  from iotbx.data_manager import DataManager
  dm = DataManager()
  m = dm.get_real_map(file_name)
  m.shift_origin()

  print("MIN MAX MEAN RMS of map: %7.2f %7.2f  %7.2f  %7.2f " %(
      m.header_min, m.header_max, m.header_mean, m.header_rms), file = out)
  print("grid: ", m.unit_cell_grid, file = out)
  print("cell:  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  " %tuple(
       m.unit_cell().parameters()), file = out)
  print("SG: ", m.unit_cell_crystal_symmetry().space_group_number(), file = out)
  if must_allow_sharpening and m.cannot_be_sharpened():
      raise Sorry("Input map is already modified and should not be sharpened")
  if get_map_labels:
      map_labels = m.labels
  print("ORIGIN: ", m.map_data().origin(), file = out)
  print("EXTENT: ", m.map_data().all(), file = out)
  print("IS PADDED: ", m.map_data().is_padded(), file = out)

  map_data = m.data
  acc = map_data.accessor()
  origin_frac = (
     m.origin_shift_grid_units[0]/m.map_data().all()[0],
     m.origin_shift_grid_units[1]/m.map_data().all()[1],
     m.origin_shift_grid_units[2]/m.map_data().all()[2])
  # determine if we need to trim off the outer part of the map duplicating inner
  offsets = []
  need_offset = False
  for g, e in zip(m.unit_cell_grid, map_data.all() ):
    offset = e-g
    offsets.append(offset)
  if  offsets  ==  [1, 1, 1]:
    if origin_frac!= (0., 0., 0.):  # this was a shifted map...we can't do this
      raise Sorry("Sorry if a CCP4 map has an origin other than (0, 0, 0) "+
        "the extent \nof the map must be the same as the grid or 1 "+
        "\ngreater for "+
        "segment_and_split_map routines."+
       "The file %s has a grid of %s and extent of %s" %(
       file_name, str(m.unit_cell_grid), str(map_data.all())))
    map = map_data[:-1, :-1, :-1]
    acc = map.accessor()
  else:
    map = map_data

  # now get space group and cell
  from cctbx import crystal
  from cctbx import sgtbx
  if m.unit_cell_crystal_symmetry().space_group_number() == 0:
    n = 1 # fix mrc formatting
  else:
    n = m.unit_cell_crystal_symmetry().space_group_number()
  if hasattr(m, 'crystal_symmetry'):
    space_group_info = sgtbx.space_group_info(number = n)
    unit_cell = m.unit_cell_crystal_symmetry().unit_cell()
    original_unit_cell_grid = m.unit_cell_grid
    original_crystal_symmetry = crystal.symmetry(
      unit_cell = unit_cell, space_group_info = space_group_info)
    if original_crystal_symmetry and map.all() == m.unit_cell_grid:
      crystal_symmetry = original_crystal_symmetry
      print("\nUnit cell crystal symmetry used: ", file = out)
    else:
      crystal_symmetry = m.crystal_symmetry()
      print("\nBox crystal symmetry used: ", file = out)
    crystal_symmetry.show_summary(f = out)
    space_group = crystal_symmetry.space_group()
    unit_cell = crystal_symmetry.unit_cell()
  else:
    space_group = None
    unit_cell = None
    crystal_symmetry = None
    original_crystal_symmetry, original_unit_cell_grid = None, None

  map = scale_map(map, out = out)
  if get_map_labels:
    return map, space_group, unit_cell, crystal_symmetry, origin_frac, acc, \
      original_crystal_symmetry, original_unit_cell_grid, map_labels
  else:
    return map, space_group, unit_cell, crystal_symmetry, origin_frac, acc, \
      original_crystal_symmetry, original_unit_cell_grid

def write_ccp4_map(crystal_symmetry, file_name, map_data,
   output_unit_cell_grid = None, labels = None):
  """Write a ccp4-style map with map data to file_name.  Deprecated. Use
  data_manager instead"""
  if output_unit_cell_grid is None:
    output_unit_cell_grid = map_data.all()
  if labels is None:
    labels = flex.std_string([""])

  iotbx.mrcfile.write_ccp4_map(
      file_name = file_name,
      unit_cell = crystal_symmetry.unit_cell(),
      space_group = crystal_symmetry.space_group(),
      unit_cell_grid = output_unit_cell_grid,
      map_data = map_data.as_double(),
      labels = labels)

def set_up_xrs(crystal_symmetry = None):  # dummy xrs to write out atoms
  """Create a dummy xray_structure to write out atoms. Use with write_atoms"""

  lines = ["ATOM     92  SG  CYS A  10       8.470  28.863  18.423  1.00 22.05           S"] # just a random line to set up x-ray structure
  from iotbx.pdb.utils import get_pdb_input
  pdb_inp = get_pdb_input(lines = lines)
  xrs = pdb_inp.xray_structure_simple(crystal_symmetry = crystal_symmetry)
  scatterers = flex.xray_scatterer()
  return xrs, scatterers

def write_atoms(tracking_data = None, sites = None, file_name = None,
      crystal_symmetry = None,
      atom_name = None, resname = None, atom_type = None, occ = None,
      out = sys.stdout):
    """Write out a file with coordinates in sites and supplied atom_name,
    resname, atom_type, and occupancy."""
    if crystal_symmetry is None:
       crystal_symmetry = tracking_data.crystal_symmetry
    xrs, scatterers = set_up_xrs(crystal_symmetry = crystal_symmetry)
    from cctbx import xray
    unit_cell = crystal_symmetry.unit_cell()
    for xyz_cart in sites:
      scatterers.append( xray.scatterer(scattering_type = "O",
         label = "O",
        site = unit_cell.fractionalize(xyz_cart), u = 0.38, occupancy = 1.0))
    text = write_xrs(xrs = xrs, scatterers = scatterers, file_name = file_name, out = out)
    if atom_name and resname and atom_type:
      text = text.replace("O      O  ", " %2s  %3s A" %(atom_name, resname) )
      text = text.replace("           O", "           %1s" %(atom_type))
    if occ:
      text = text.replace(" 1.00 ", " %.2f " %(occ))
    return text


def write_xrs(
    xrs = None, scatterers = None, file_name = "atoms.pdb", out = sys.stdout):
  """Write an xray_structure to a PDB file."""
  from cctbx import xray
  xrs = xray.structure(xrs, scatterers = scatterers)
  text = xrs.as_pdb_file()  # PDB OK just writing out some atoms
  if file_name:
    f = open(file_name, 'w')
    print(text, file = f)
    f.close()
    print("Atoms written to %s" %file_name, file = out)
  return text

def get_b_iso(miller_array, d_min = None, return_aniso_scale_and_b = False,
    d_max = 100000.):
  """Estimate b_iso value from a miller array with amplitudes"""

  if d_min:
    res_cut_array = miller_array.resolution_filter(d_max = d_max,
       d_min = d_min)
  else:
    res_cut_array = miller_array

  from mmtbx.scaling import absolute_scaling
  try:
    aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(
      miller_array = res_cut_array, n_residues = 200, n_bases = 0, ignore_errors = True)
    b_cart = aniso_scale_and_b.b_cart
  except Exception as e:
    b_cart = [0, 0, 0]
    aniso_scale_and_b = None
  b_aniso_mean = 0.
  if b_cart:
    for k in [0, 1, 2]:
      b_aniso_mean+= b_cart[k]
  if return_aniso_scale_and_b:
    return b_aniso_mean/3.0, aniso_scale_and_b
  else: # usual
    return b_aniso_mean/3.0

def map_coeffs_as_fp_phi(map_coeffs):
  """Return amplitudes and phases arrays based on map_coeffs array"""
  amplitudes = map_coeffs.amplitudes()
  amplitudes.set_observation_type_xray_amplitude()
  assert amplitudes.is_real_array()
  phases = map_coeffs.phases(deg = True)
  return amplitudes, phases

def map_coeffs_to_fp(map_coeffs):
  """Return amplitudes array based on map_coeffs array"""
  amplitudes = map_coeffs.amplitudes()
  amplitudes.set_observation_type_xray_amplitude()
  assert amplitudes.is_real_array()
  return amplitudes

def get_f_phases_from_model(f_array = None, pdb_hierarchy= None,
      overall_b = None,
     k_sol = None, b_sol = None, out = sys.stdout):
  """Calculate amplitudes and phases from a model.
  XXX ignores values of k_sol and b_sol.
  Returns complex array"""
  xray_structure = pdb_hierarchy.extract_xray_structure(
     crystal_symmetry = f_array.crystal_symmetry())
  print("Getting map coeffs from model with %s atoms.." %(
    xray_structure.sites_frac().size()), file = out)


  if overall_b is not None:
    print("Setting overall b_iso to %7.1f for model " %(
      overall_b), file = out)
    xray_structure.set_b_iso(value = overall_b)
  model_f_array = f_array.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()

  return model_f_array

def get_f_phases_from_map(
      map_data = None, crystal_symmetry = None, d_min = None,
      d_max = 100000.,
      d_min_ratio = None,
      return_as_map_coeffs = False,
      remove_aniso = None,
      get_remove_aniso_object = True,
      scale_max = None,
      origin_frac = None,
        out = sys.stdout):

    """ Calculate structure factors from a map.  Returns map coefficients
        in space group P1 always (regardless of crystal_symmetry).
    """

    if d_min is not None:
      d_min_use = d_min
      if d_min_ratio is not None:
        d_min_use = d_min*d_min_ratio
    else:
      d_min_use = None
    if map_data.origin() != (0,0,0):
      map_data = map_data.shift_origin()
    from iotbx.map_manager import map_manager
    mm = map_manager(map_data = map_data,
       unit_cell_grid = map_data.all(),
       unit_cell_crystal_symmetry = crystal_symmetry,
       wrapping = False)
    map_coeffs = mm.map_as_fourier_coefficients(d_min = d_min_use,
       d_max = d_max if d_min_use is not None else None)

    if map_coeffs.data().size() < 1:
      raise Sorry("No map coefficients found in calculation of map at "+
         "resolution of %.2f A." %( d_min_use))
    if origin_frac and tuple(origin_frac) !=  (0., 0., 0.):  # shift origin
      map_coeffs = map_coeffs.translational_shift(origin_frac, deg = False)

    map_coeffs = scale_map_coeffs(map_coeffs, scale_max = scale_max, out = out)

    if remove_aniso:
      print("\nRemoving aniso in data before analysis\n", file = out)
      get_remove_aniso_object = True

    from cctbx.maptbx.refine_sharpening import analyze_aniso
    map_coeffs, map_coeffs_ra = analyze_aniso(
         remove_aniso = remove_aniso,
         get_remove_aniso_object = get_remove_aniso_object,
         map_coeffs = map_coeffs, resolution = d_min, out = out)

    if return_as_map_coeffs:
      return map_coeffs, map_coeffs_ra
    else:
      return map_coeffs_as_fp_phi(map_coeffs)

def apply_sharpening(map_coeffs = None,
    sharpening_info_obj = None,
    n_real = None, b_sharpen = None, crystal_symmetry = None,
    target_scale_factors = None,
    f_array = None, phases = None, d_min = None, k_sharpen = None,
    b_blur_hires = None,
    include_sharpened_map_coeffs = False,
    out = sys.stdout):
    """Apply sharpening information in sharpening_info_obj to map_coeffs
    and return map_and_b object with map_data, starting and final b_iso"""
    if map_coeffs and f_array is None and phases is None:
      f_array, phases = map_coeffs_as_fp_phi(map_coeffs)

    if sharpening_info_obj is not None:
      b_sharpen = sharpening_info_obj.b_sharpen
      b_blur_hires = sharpening_info_obj.b_blur_hires
      k_sharpen = sharpening_info_obj.k_sharpen
      if sharpening_info_obj.input_d_cut:
        d_min = sharpening_info_obj.input_d_cut
      else:
        d_min = sharpening_info_obj.resolution# changed from d_cut
      n_real = sharpening_info_obj.n_real
      target_scale_factors = sharpening_info_obj.target_scale_factors
      n_bins = sharpening_info_obj.n_bins
      remove_aniso = sharpening_info_obj.remove_aniso
      resolution = sharpening_info_obj.resolution

    if target_scale_factors:
      assert sharpening_info_obj is not None
      print("\nApplying target scale factors vs resolution", file = out)
      if not map_coeffs:
        map_coeffs = f_array.phase_transfer(phase_source = phases, deg = True)
      f_array, phases = map_coeffs_as_fp_phi(map_coeffs)
      f_array_b_iso = get_b_iso(f_array, d_min = d_min)
      if not f_array.binner():
        (local_d_max, local_d_min) = f_array.d_max_min(
          d_max_is_highest_defined_if_infinite=True)
        f_array.setup_binner(n_bins = n_bins, d_max = local_d_max,
        d_min = local_d_min)

      from cctbx.maptbx.refine_sharpening import apply_target_scale_factors
      map_and_b = apply_target_scale_factors(f_array = f_array, phases = phases,
        resolution = d_min,
        target_scale_factors = target_scale_factors,
        n_real = n_real,
        out = out)
      return map_and_b

    elif b_sharpen is None or (
        b_sharpen in [0, None] and k_sharpen in [0, None]):
      if not map_coeffs:
        map_coeffs = f_array.phase_transfer(phase_source = phases, deg = True)
      map_data = get_map_from_map_coeffs(map_coeffs = map_coeffs,
        crystal_symmetry = crystal_symmetry, n_real = n_real)
      return map_and_b_object(map_data = map_data)

    elif k_sharpen is None or d_min is None or k_sharpen<= 0 or \
        ( b_blur_hires is None and b_sharpen < 0):
      # 2016-08-10 original method: apply b_sharpen to all data
      # Use this if blurring (b_sharpen<0) or if k_sharpen is not set
      from cctbx import adptbx # next lines from xtriage (basic_analysis.py)
      b_cart_aniso_removed = [ b_sharpen, b_sharpen, b_sharpen, 0, 0, 0]
      from mmtbx.scaling import absolute_scaling
      u_star_aniso_removed = adptbx.u_cart_as_u_star(
        f_array.unit_cell(), adptbx.b_as_u( b_cart_aniso_removed  ) )
      f_array_sharpened = absolute_scaling.anisotropic_correction(
        f_array, 0.0, u_star_aniso_removed, must_be_greater_than = -0.0001)
    else:
      # Apply sharpening only to data from infinity to d_min, with transition
      # steepness of k_sharpen.
      # 2017-08-21 if b_blur_hires is set, sharpen with
      # b_sharpen-b_blur_hires data beyond d_min (with same
      # transition, so transition goes from b_sharpen TO b_sharpen-b_blur_hires
      data_array = f_array.data()
      sthol_array = f_array.sin_theta_over_lambda_sq()
      d_spacings = f_array.d_spacings()
      scale_array = flex.double()
      import math

      if b_blur_hires is not None:
        b_sharpen_hires_use = b_sharpen-b_blur_hires
      else:
        b_sharpen_hires_use = 0.
      for x, (ind, sthol), (ind1, d) in zip(data_array, sthol_array, d_spacings):
        # for small value b = b_sharpen
        # for large value b = -b_sharpen_hires_use
        # transition is determined by k_sharpen
        value = min(20., max(-20., k_sharpen*(d_min-d)))
        lowres_weight = 1./(1.+math.exp(value))
        hires_weight = max(0., 1-lowres_weight)
        b_sharpen_use = b_sharpen*lowres_weight+b_sharpen_hires_use*hires_weight
        log_scale = sthol*b_sharpen_use
        scale_array.append(math.exp(log_scale))
      data_array = data_array*scale_array
      f_array_sharpened = f_array.customized_copy(data = data_array)

    actual_b_iso = get_b_iso(f_array_sharpened, d_min = d_min)
    print("B-iso after sharpening by b_sharpen = %6.1f is %7.2f\n" %(
      b_sharpen, actual_b_iso), file = out)
    sharpened_map_coeffs = f_array_sharpened.phase_transfer(
      phase_source = phases, deg = True)

    # And get new map
    map_data = get_map_from_map_coeffs(map_coeffs = sharpened_map_coeffs,
      crystal_symmetry = crystal_symmetry,
       n_real = n_real)
    mb = map_and_b_object(map_data = map_data, final_b_iso = actual_b_iso)
    if include_sharpened_map_coeffs:
      mb.sharpened_map_coeffs = sharpened_map_coeffs
    return mb

def find_symmetry_center(map_data, crystal_symmetry = None, out = sys.stdout):
  """Find center of symmetry for map_data"""
  # find center if necessary:
  origin = list(map_data.origin())
  all = list(map_data.all())
  centroid_wx = {}
  centroid_w = {}
  from cctbx import maptbx
  for ai in [0, 1, 2]:
    centroid_wx[ai] = 0.
    centroid_w[ai] = 0.
    for i in range(0, all[ai]):
      if ai == 0:
        start_tuple = tuple((i, 0, 0))
        end_tuple = tuple((i, all[1]-1, all[2]-1))  #2019-11-05 not beyond na-1
      elif ai == 1:
         start_tuple = tuple((0, i, 0))
         end_tuple = tuple((all[0]-1, i, all[2]-1))
      elif ai == 2:
         start_tuple = tuple((0, 0, i))
         end_tuple = tuple((all[0]-1, all[1]-1, i))
      new_map_data = maptbx.copy(map_data,
         start_tuple, end_tuple)
      mean_value = max(0., new_map_data.as_1d().as_double().min_max_mean().mean)
      centroid_wx[ai]+= mean_value*(i-origin[ai])
      centroid_w[ai]+= mean_value
    if centroid_w[ai]>0:
      centroid_wx[ai] = centroid_wx[ai]/centroid_w[ai]
  print("CENTROID OF DENSITY: (%7.2f, %7.2f, %7.2f) (grid units) " %(
    tuple((centroid_wx[0], centroid_wx[1], centroid_wx[2], ))), file = out)
  xyz_fract = matrix.col((centroid_wx[0]/all[0], centroid_wx[1]/all[1], centroid_wx[2]/all[2], ))
  xyz_cart = crystal_symmetry.unit_cell().orthogonalize(xyz_fract)
  print("CENTROID (A): (%7.3f, %7.3f, %7.3f) " %(
    tuple(xyz_cart)), file = out)
  return xyz_cart

def get_center_of_map(map_data, crystal_symmetry,
    round_down = True,
    place_on_grid_point = True):
  """Find center of map. Deprecated, use map_manager.absolute_center_cart()"""
  all = list(map_data.all())
  origin = list(map_data.origin())
  if place_on_grid_point:
    if round_down:
      sx, sy, sz = [int(all[0]/2)+origin[0], int(all[1]/2)+origin[1],
         int(all[2]/2)+origin[2]]
    else:
      sx, sy, sz = [int(1+all[0]/2)+origin[0], int(1+all[1]/2)+origin[1],
         int(1+all[2]/2)+origin[2]]
  else:
    sx, sy, sz = [all[0]/2+origin[0], all[1]/2+origin[1],
  all[2]/2+origin[2]]
  site_fract = matrix.col((sx/all[0], sy/all[1], sz/all[2], ))
  return crystal_symmetry.unit_cell().orthogonalize(site_fract)


def select_remaining_ncs_ops( map_data = None,
    crystal_symmetry = None,
    random_points = None,
    closest_sites = None,
    ncs_object = None,
    out = sys.stdout):

  """Identify which NCS ops still apply.  Choose the ones that maximize
  scoring with score_ncs_in_map."""

  if ncs_object.max_operators()<1:
    return ncs_object

  used_ncs_id_list = [ncs_object.ncs_groups()[0].identity_op_id()]
  ncs_copies = ncs_object.max_operators()

  # find ncs_id that maximizes score (if any)
  improving = True
  from copy import deepcopy
  best_ops_to_keep = deepcopy(used_ncs_id_list)
  working_best_ops_to_keep = None
  best_score = None

  while improving:
    improving = False
    working_best_ops_to_keep = deepcopy(best_ops_to_keep)
    working_score = None
    for ncs_id in range(ncs_copies):
      if ncs_id in best_ops_to_keep:continue
      ops_to_keep = deepcopy(best_ops_to_keep)
      ops_to_keep.append(ncs_id)
      ncs_used_obj = ncs_object.deep_copy(ops_to_keep = ops_to_keep)
      score, ncs_cc = score_ncs_in_map(map_data = map_data, ncs_object = ncs_used_obj,
        ncs_in_cell_only = True,
        allow_score_with_pg = False,
        sites_orth = closest_sites,
        crystal_symmetry = crystal_symmetry, out = null_out())
      if score is None: continue
      if working_score is None or score >working_score:
        working_score = score
        working_best_ops_to_keep = deepcopy(ops_to_keep)
    if working_score is not None and (
        best_score is None or working_score>best_score):
      improving = True
      best_score = working_score
      best_ops_to_keep = deepcopy(working_best_ops_to_keep)

  ncs_used_obj = ncs_object.deep_copy(ops_to_keep = best_ops_to_keep)
  return ncs_used_obj

def run_get_ncs_from_map(params = None,
      map_data = None,
      crystal_symmetry = None,
      map_symmetry_center = None,
      ncs_obj = None,
      fourier_filter = False,
      out = sys.stdout,
      ):

  """Run the get_ncs_from_map method, to get or check NCS operators
  Try various possibilities for center of NCS."""

  # First Fourier filter map if resolution is set
  if fourier_filter and params.crystal_info.resolution:
    print("Fourier filtering at resolution of %.2f A" %(
      params.crystal_info.resolution), file = out)
    from iotbx.map_manager import map_manager
    mm = map_manager(map_data= map_data,
      unit_cell_crystal_symmetry = crystal_symmetry,
      unit_cell_grid = map_data.all(),
      wrapping=False)
    mm.resolution_filter(d_min=params.crystal_info.resolution)
    map_data = mm.map_data()
  sym_cen = params.reconstruction_symmetry.symmetry_center

  ncs_obj_to_check = None

  if params.reconstruction_symmetry.symmetry and (
     not ncs_obj or ncs_obj.max_operators()<2):
    if params.reconstruction_symmetry.optimize_center:
                       # [use_center_of_map, round_down, symmetry_center]
      center_try_list = [[True,True,sym_cen],
                         [False,True,sym_cen],
                         [True,True,None],
                         [True,False,None],
                        ]
    else:
      center_try_list = [[True,True,sym_cen],
                         [True,True,None],
                         [True,False,None],
                        ]
  elif ncs_obj:
    center_try_list = [[True,True,sym_cen],
                         [True,True,None],
                         [True,False,None],
                        ]
    ncs_obj_to_check = ncs_obj
  elif params.reconstruction_symmetry.optimize_center:
    center_try_list = [[True,True,sym_cen],
                         [True,True,None],
                         [True,False,None],
                        ]
  else:
    return None, None, None # did not even try
  # check separately for helical symmetry
  if params.reconstruction_symmetry.symmetry and \
       params.reconstruction_symmetry.symmetry.lower() == 'helical':
    helical_list = [True]
  elif params.reconstruction_symmetry.symmetry and \
       params.reconstruction_symmetry.symmetry.lower() in ['all', 'any'] and\
      params.reconstruction_symmetry.include_helical_symmetry:
    helical_list = [False, True]
  else:
    helical_list = [False]

  new_ncs_obj, ncs_cc, ncs_score = None, None, None
  for [use_center_of_map, round_down, symmetry_center] in center_try_list:
   for include_helical in helical_list:
    local_params = deepcopy(params)
    local_params.reconstruction_symmetry.symmetry_center=None
    local_params.reconstruction_symmetry.include_helical_symmetry = \
       include_helical
    new_ncs_obj, ncs_cc, ncs_score = get_ncs_from_map(params = local_params,
      map_data = map_data,
      map_symmetry_center = map_symmetry_center,
      use_center_of_map_as_center = use_center_of_map,
      crystal_symmetry = crystal_symmetry,
      ncs_obj_to_check = ncs_obj_to_check,
      symmetry_center = symmetry_center,
      round_down = round_down,
      out = out
      )
    if new_ncs_obj:
      return new_ncs_obj, ncs_cc, ncs_score
  return new_ncs_obj, ncs_cc, ncs_score

def get_ncs_from_map(params = None,
      map_data = None,
      map_symmetry_center = None,
      symmetry = None,
      symmetry_center = None,
      helical_rot_deg = None,
      helical_trans_z_angstrom = None,
      two_fold_along_x = None,
      op_max = None,
      crystal_symmetry = None,
      optimize_center = None,
      sites_orth = None,
      random_points = None,
      n_rescore = None,
      use_center_of_map_as_center = None,
      min_ncs_cc = None,
      identify_ncs_id = None,
      ncs_obj_to_check = None,
      ncs_in_cell_only = False,
      round_down = True,
      out = sys.stdout):


  """Check through standard point groups and helical symmetry to see
   if map has symmetry. If symmetry == ANY then take highest symmetry that fits
   Otherwise limit to the one specified with symmetry.
    Use a library of symmetry matrices.  For helical symmetry generate it
    along the z axis.
   Center of symmetry is as supplied, or center of map or center of density
    If center is not supplied and use_center_of_map_as_center, try that
    and return None if it fails to achieve a map cc of min_ncs_cc"""

  # round_down defines where to guess center if gridding has even number of pts

  if ncs_in_cell_only is None:
    ncs_in_cell_only = (not params.crystal_info.use_sg_symmetry)
  if symmetry is None:
    symmetry = params.reconstruction_symmetry.symmetry
  if symmetry_center is None:
    symmetry_center = params.reconstruction_symmetry.symmetry_center
  if optimize_center is None:
    optimize_center = params.reconstruction_symmetry.optimize_center
  if helical_rot_deg is None:
    helical_rot_deg = params.reconstruction_symmetry.helical_rot_deg
  if helical_trans_z_angstrom is None:
    helical_trans_z_angstrom = \
      params.reconstruction_symmetry.helical_trans_z_angstrom
  if n_rescore is None:
    n_rescore = params.reconstruction_symmetry.n_rescore
  if random_points is None:
    random_points = params.reconstruction_symmetry.random_points
  if op_max is None:
    op_max = params.reconstruction_symmetry.op_max
  if two_fold_along_x is None:
    two_fold_along_x = params.reconstruction_symmetry.two_fold_along_x
  if identify_ncs_id is None:
    identify_ncs_id = params.reconstruction_symmetry.identify_ncs_id
  if min_ncs_cc is None:
    min_ncs_cc = params.reconstruction_symmetry.min_ncs_cc

  # if ncs_obj_to_check is supplied...just use that ncs
  if ncs_obj_to_check and ncs_obj_to_check.max_operators()>1:
    symmetry = "SUPPLIED NCS"

  if map_symmetry_center is None:
    map_symmetry_center = get_center_of_map(map_data, crystal_symmetry,
      round_down = round_down)
  if optimize_center is None:
    if symmetry_center is None and (not use_center_of_map_as_center):
      optimize_center = True
      print("Setting optimize_center = True as no symmetry_center is supplied", file = out)
    else:
      optimize_center = False

  if symmetry_center is not None:
    symmetry_center = matrix.col(symmetry_center)
  elif use_center_of_map_as_center:
    print("Using center of map as NCS center", file = out)
    symmetry_center = map_symmetry_center
  else: # Find it
    if not ncs_obj_to_check:
      print("Finding NCS center as it is not supplied", file = out)
    symmetry_center = find_symmetry_center(
    map_data, crystal_symmetry = crystal_symmetry,
       out = out)
  print("Center of NCS (A): (%7.3f, %7.3f, %7.3f) " %(
    tuple(symmetry_center)), file = out)

  ncs_list = get_ncs_list(params = params,
    symmetry = symmetry,
   symmetry_center = symmetry_center,
   helical_rot_deg = helical_rot_deg,
   two_fold_along_x = two_fold_along_x,
   op_max = op_max,
   helical_trans_z_angstrom = helical_trans_z_angstrom,
   ncs_obj_to_check = ncs_obj_to_check,
   map_data = map_data,
   crystal_symmetry = crystal_symmetry,
   out = out,
   )

  print("Total of %d NCS types to examine..." %(len(ncs_list)), file = out)
  if not sites_orth:
    sites_orth = get_points_in_map(
     map_data, n = random_points, crystal_symmetry = crystal_symmetry)
  # some random points in the map

  # Now make sure symmetry applied to points in points_list gives similar values

  results_list = []
  for ncs_obj in ncs_list:
    symmetry = ncs_obj.get_ncs_name()
    score, cc_avg = score_ncs_in_map(map_data = map_data, ncs_object = ncs_obj,
       identify_ncs_id = identify_ncs_id,
       ncs_in_cell_only = ncs_in_cell_only,
      sites_orth = sites_orth, crystal_symmetry = crystal_symmetry, out = out)
    if cc_avg is None or cc_avg < min_ncs_cc:
      score = 0. # Do not allow low CC values to be used
    if score is None:
      print("symmetry:", symmetry, " no score", ncs_obj.max_operators(), file = out)
    else:
      results_list.append([score, cc_avg, ncs_obj, symmetry])
  if not results_list:
    return None, None, None

  results_list.sort(key=itemgetter(0))
  results_list.reverse()

  # Rescore top n_rescore
  if n_rescore and not ncs_obj_to_check:
    print("Rescoring top %d results" %(min(n_rescore, len(results_list))), file = out)
    rescore_list = results_list[n_rescore:]
    all_rescore = []
    new_sites_orth = get_points_in_map(
      map_data, n = 10*random_points, crystal_symmetry = crystal_symmetry)
    new_sites_orth.extend(sites_orth)
    found_ok = False
    for orig_score, orig_cc_avg, ncs_obj, symmetry in results_list[:n_rescore]:
      score, cc_avg = score_ncs_in_map(map_data = map_data, ncs_object = ncs_obj,
        identify_ncs_id = identify_ncs_id,
        ncs_in_cell_only = ncs_in_cell_only,
        sites_orth = new_sites_orth, crystal_symmetry = crystal_symmetry, out = out)
      all_rescore.append([score, cc_avg, ncs_obj, symmetry])
      if cc_avg is None or cc_avg < min_ncs_cc:
        score = 0. # Do not allow low CC values to be used
      else:
        found_ok = True
      if score is None:
        print("symmetry:", symmetry, " no score", ncs_obj.max_operators(), file = out)
      else:
        rescore_list.append([score, cc_avg, ncs_obj, symmetry])
    if not found_ok:
      rescore_list = all_rescore
    rescore_list.sort(key=itemgetter(0))
    rescore_list.reverse()
    results_list = rescore_list
  if len(results_list) == 1:
    # check for C1
    score, cc_avg, ncs_obj, symmetry = results_list[0]
    if symmetry and symmetry.strip() == 'C1':
      score = 1.
      cc_avg = 1.
      results_list = [[score, cc_avg, ncs_obj, symmetry], ]


  print("Ranking of NCS types:", file = out)
  if min_ncs_cc is not None:
    print("NOTE: any NCS type with CC < %.2f (min_ncs_cc) is unscored " %(
      min_ncs_cc), file = out)
  print("\n  SCORE    CC   OPERATORS     SYMMETRY", file = out)
  for score, cc_avg, ncs_obj, symmetry in results_list:
    if not symmetry: symmetry = ""
    if not cc_avg: cc_avg = 0.0


    print(" %6.2f  %5.2f    %2d          %s" %(
       score, cc_avg, ncs_obj.max_operators(), symmetry.strip(), ), file = out)

  score, cc_avg, ncs_obj, ncs_info = results_list[0]
  # check for offset by gridding
  if  hasattr(params.reconstruction_symmetry,'check_grid_offset') and \
     params.reconstruction_symmetry.check_grid_offset:
    symmetry_center, cc_avg, score, ncs_obj = optimize_center_position(
       map_data, sites_orth,
       crystal_symmetry,
       ncs_info, symmetry_center, ncs_obj, score, cc_avg,
       params = params,
       helical_rot_deg = helical_rot_deg,
       two_fold_along_x = two_fold_along_x,
       op_max = op_max,
       min_ncs_cc = min_ncs_cc,
       identify_ncs_id = identify_ncs_id,
       ncs_obj_to_check = ncs_obj_to_check,
       ncs_in_cell_only = ncs_in_cell_only,
       helical_trans_z_angstrom = helical_trans_z_angstrom,
       check_grid_offset = True,
       out = out)

  # Optimize center if necessary
  if optimize_center:
    symmetry_center, cc_avg, score, ncs_obj = optimize_center_position(
       map_data, sites_orth,
       crystal_symmetry,
       ncs_info, symmetry_center, ncs_obj, score, cc_avg,
       params = params,
       helical_rot_deg = helical_rot_deg,
       two_fold_along_x = two_fold_along_x,
       op_max = op_max,
       min_ncs_cc = min_ncs_cc,
       identify_ncs_id = identify_ncs_id,
       ncs_obj_to_check = ncs_obj_to_check,
       ncs_in_cell_only = ncs_in_cell_only,
       helical_trans_z_angstrom = helical_trans_z_angstrom, out = out)
    print("New center: (%7.3f, %7.3f, %7.3f)" %(tuple(symmetry_center)), file = out)

  if (not cc_avg) or (cc_avg < min_ncs_cc):
    print("No suitable symmetry found", file = out)
    return None, None, None

  print("\nBest NCS type is: ", end = ' ', file = out)
  print("\n  SCORE    CC   OPERATORS     SYMMETRY", file = out)
  if not ncs_info: ncs_info = ""
  print(" %6.2f  %5.2f    %2d          %s  Best NCS type" %(
       score, cc_avg, ncs_obj.max_operators(), ncs_info.strip(), ), file = out)
  return ncs_obj, cc_avg, score


def optimize_center_position(map_data, sites_orth, crystal_symmetry,
     ncs_info, symmetry_center, ncs_obj, score, cc_avg,
     params = None,
     helical_rot_deg = None,
     two_fold_along_x = None,
     op_max = None,
     identify_ncs_id = None,
     ncs_obj_to_check = None,
     ncs_in_cell_only = None,
     min_ncs_cc = None,
     helical_trans_z_angstrom = None,
     check_grid_offset = False,
     out = sys.stdout):

  if ncs_info is None:
    ncs_info = "None"

  symmetry = ncs_info.split()[0]
  if check_grid_offset:
    print(
      "Checking for improvement with offset by +/-1 grid unit...type is %s" %(
       ncs_info), file = out)
  else:
    print("Optimizing center position...type is %s" %(ncs_info), file = out)

  if len(ncs_info.split())>1 and ncs_info.split()[1] == '(a)':
    two_fold_along_x = True
  elif len(ncs_info.split())>1 and ncs_info.split()[1] == '(b)':
    two_fold_along_x = False
  else:
    two_fold_along_x = None

  best_center = matrix.col(symmetry_center)
  best_ncs_obj = ncs_obj
  best_score = score
  best_cc_avg = cc_avg
  if (ncs_info.find("Helical") > -1) and check_grid_offset:
    return best_center, best_cc_avg, best_score, best_ncs_obj
  print("Starting center: (%7.3f, %7.3f, %7.3f)" %(tuple(best_center)), file = out)
  if check_grid_offset:
    n_range = 1
    k_range = list(range(-1,2))
    i_range = list(range(-1,2))
    j_range = list(range(-1,2))
    abc = crystal_symmetry.unit_cell().parameters()[:3]
    N_ = map_data.all()

    scale_x = abc[0]/N_[0]
    scale_y = abc[1]/N_[1]
    scale_z = abc[2]/N_[2]
  else:
    n_range = 6
    k_range = range(0,1)
    i_range = range(-4,5)
    j_range = range(-4,5)
    scale_x = 1
    scale_y = 1
    scale_z = 1
  from libtbx.utils import null_out
  scale = 5.
  for itry in range(n_range):
   scale = scale/5.
   for k in k_range:
    for i in i_range:
     for j in j_range:
      local_center = matrix.col(symmetry_center)+matrix.col((scale_x*i, scale_y*j, scale_z*k, ))
      ncs_list = get_ncs_list(params = params, symmetry = symmetry,
       symmetry_center = local_center,
       helical_rot_deg = helical_rot_deg,
       two_fold_along_x = two_fold_along_x,
       op_max = op_max,
       helical_trans_z_angstrom = helical_trans_z_angstrom,
       ncs_obj_to_check = ncs_obj_to_check,
       map_data = map_data,
       crystal_symmetry = crystal_symmetry,
       out = null_out(),
       )
      if ncs_list:
        ncs_obj = ncs_list[0]
        score, cc_avg = score_ncs_in_map(map_data = map_data, ncs_object = ncs_obj,
          identify_ncs_id = identify_ncs_id,
          ncs_in_cell_only = ncs_in_cell_only,
          sites_orth = sites_orth, crystal_symmetry = crystal_symmetry, out = out)
        if cc_avg < min_ncs_cc:
          score = 0. # Do not allow low CC values to be used
      else:
        ncs_obj = None
        score, cc_avg = None, None
      if best_score is None or score>best_score:
        best_cc_avg = cc_avg
        best_score = score
        best_center = local_center
        best_ncs_obj = ncs_obj

  symmetry_center = best_center
  cc_avg = best_cc_avg
  score = best_score
  ncs_obj = best_ncs_obj
  print("Optimized center: (%7.3f, %7.3f, %7.3f)" %(tuple(best_center)), file = out)
  return best_center, best_cc_avg, best_score, best_ncs_obj


def score_ncs_in_map_point_group_symmetry(
    map_data = None, ncs_object = None, sites_orth = None,
     crystal_symmetry = None, out = sys.stdout):
  ncs_group = ncs_object.ncs_groups()[0]
  all_value_lists = []
  for c, t, r in zip(ncs_group.centers(),
                       ncs_group.translations_orth(),
                       ncs_group.rota_matrices()):
    new_sites_cart = flex.vec3_double()
    r_inv = r.inverse()
    for site in sites_orth:
      new_sites_cart.append(r_inv * (matrix.col(site) - t))
    # get value at new_sites cart and make sure they are all the same...
    new_sites_fract = crystal_symmetry.unit_cell().fractionalize(new_sites_cart)
    values = flex.double()
    for site_fract in new_sites_fract:
      values.append(map_data.value_at_closest_grid_point(site_fract))
    all_value_lists.append(values)
  return get_cc_among_value_lists(all_value_lists)

def get_cc_among_value_lists(all_value_lists):
  if not all_value_lists:
    return None, None
  a = all_value_lists[0]
  cc_avg = 0.
  cc_low = None
  cc_n = 0.
  for j in range(1, len(all_value_lists)):
      b = all_value_lists[j]
      cc = flex.linear_correlation(a, b).coefficient()
      cc_avg+= cc
      cc_n+= 1.
      if cc_low is None or cc<cc_low:
        cc_low = cc
  cc_avg = cc_avg/max(1., cc_n)
  if cc_n>0:
    import math
    return cc_low*math.sqrt(len(all_value_lists)), cc_avg
  else:
    return None, None

def score_ncs_in_map(map_data = None, ncs_object = None, sites_orth = None,
     identify_ncs_id = None,
     ncs_in_cell_only = None,
     allow_score_with_pg = True,
     crystal_symmetry = None, out = sys.stdout):
  if not ncs_object or ncs_object.max_operators()<2:
    return None, None

  ncs_group = ncs_object.ncs_groups()[0]
      # don't use point-group symmetry if we have only some of the ops
  if allow_score_with_pg and (
     (not identify_ncs_id) or ncs_group.is_point_group_symmetry()):
    return score_ncs_in_map_point_group_symmetry(
     map_data = map_data, ncs_object = ncs_object,
     sites_orth = sites_orth, crystal_symmetry = crystal_symmetry, out = out)

        # This version does not assume point-group symmetry: find the NCS
        #  operator that maps each point on to all others the best, then save
  #  that list of values

  if (not ncs_group) or not ncs_group.n_ncs_oper():
    identify_ncs_id_list = [None]
  else:
    identify_ncs_id_list = list(range(ncs_group.n_ncs_oper()))+[None]

  all_value_lists = []

  if not sites_orth:
    sites_orth = get_points_in_map(map_data, n = 100,
    minimum_fraction_of_max = 0.05,
    crystal_symmetry = crystal_symmetry)

  for site in sites_orth:
    best_id = 0
    best_score = None
    best_values = None
    for site_ncs_id in identify_ncs_id_list:  #last is real one
      if site_ncs_id is None:
        site_ncs_id = best_id
        real_thing = True
      else:
        real_thing = False

      if identify_ncs_id and site_ncs_id:
        local_site = ncs_group.rota_matrices()[site_ncs_id] * matrix.col(site) + \
           ncs_group.translations_orth()[site_ncs_id]
      else:
        local_site = site
      new_sites_cart = flex.vec3_double()
      for c, t, r in zip(ncs_group.centers(),
                       ncs_group.translations_orth(),
                       ncs_group.rota_matrices()):
        r_inv = r.inverse()
        new_sites_cart.append(r_inv * (matrix.col(local_site) - t))
      new_sites_fract = crystal_symmetry.unit_cell().fractionalize(
         new_sites_cart)
      values = flex.double()
      for site_frac in new_sites_fract:
        if (not ncs_in_cell_only) or (
              site_frac[0]>= 0 and site_frac[0]<= 1 and  \
              site_frac[1]>= 0 and site_frac[1]<= 1 and  \
              site_frac[2]>= 0 and site_frac[2]<= 1):
          values.append(map_data.value_at_closest_grid_point(site_frac))
        else:
          values.append(0.)

      score = values.standard_deviation_of_the_sample()
      if real_thing or (best_score is None or score < best_score):
          best_score = score
          best_id = site_ncs_id
          best_values = values
    all_value_lists.append(best_values)

  if not all_value_lists:
    return None, None

  values_by_site_dict = {} # all_value_lists[j][i] -> values_by_site_dict[i][j]
  # there are sites_orth.size() values of j
  # there are len(ncs_group_centers) == len(all_value_lists[0]) values of i
  for i in range(len(all_value_lists[0])):
    values_by_site_dict[i] = flex.double() # value_list[0][1]
    for j in range(sites_orth.size()):
      values_by_site_dict[i].append(all_value_lists[j][i])
  new_all_values_lists = []
  for i in range(len(all_value_lists[0])):
    new_all_values_lists.append(values_by_site_dict[i])
  score, cc = get_cc_among_value_lists(new_all_values_lists)
  return score, cc

def get_points_in_map(map_data, n = None,
      minimum_fraction_of_max = 0.,
      random_xyz = None,
      max_tries_ratio = 100, crystal_symmetry = None):
  map_1d = map_data.as_1d()
  map_mean = map_1d.min_max_mean().mean
  map_max = map_1d.min_max_mean().max
  minimum_value = map_mean+minimum_fraction_of_max*(map_max-map_mean)
  points_list = flex.vec3_double()
  import random
  random.seed(1)

  nu, nv, nw = map_data.all()
  xyz_fract = crystal_symmetry.unit_cell().fractionalize(
       tuple((17.4, 27.40128571, 27.32985714, )))
  for i in range(int(max_tries_ratio*n)): # max tries
    ix = random.randint(0, nu-1)
    iy = random.randint(0, nv-1)
    iz = random.randint(0, nw-1)
    xyz_fract = matrix.col((ix/nu, iy/nv, iz/nw, ))
    value = map_data.value_at_closest_grid_point(xyz_fract)
    if value > minimum_value and value <map_max:
      if random_xyz:
         offset = []
         for i in range(3):
           offset.append((random.random()-0.5)*2.*random_xyz)
         offset = crystal_symmetry.unit_cell().fractionalize(matrix.col(offset))
         new_xyz_fract = []
         for x, o in zip(xyz_fract, offset):
           new_xyz_fract.append(max(0, min(1, x+o)))
         xyz_fract = matrix.col(new_xyz_fract)
      points_list.append(xyz_fract)
      if points_list.size()>= n: break
  sites_orth = crystal_symmetry.unit_cell().orthogonalize(points_list)
  return sites_orth


def get_ncs_list(params = None, symmetry = None,
   symmetry_center = None,
   helical_rot_deg = None,
   helical_trans_z_angstrom = None,
   op_max = None,
   two_fold_along_x = None,
   ncs_obj_to_check = None,
   crystal_symmetry = None,
   map_data = None,
   include_helical_symmetry = None,
   max_helical_ops_to_check = None,
   require_helical_or_point_group_symmetry = None,
    out = sys.stdout):
  """Generate NCS operations specified by symmetry, if requested, offset
  center of symmetry and generate helical symmetry"""
  # params.reconstruction_symmetry.require_helical_or_point_group_symmetry
  # params.reconstruction_symmetry.include_helical_symmetry):
  # params.reconstruction_symmetry.max_helical_ops_to_check))

  if ncs_obj_to_check and ncs_obj_to_check.max_operators()>1:
    return [ncs_obj_to_check] # , ["SUPPLIED NCS"]

  from mmtbx.ncs.ncs import generate_ncs_ops
  ncs_list = generate_ncs_ops(
   must_be_consistent_with_space_group_number = \
      params.reconstruction_symmetry.must_be_consistent_with_space_group_number,
   symmetry = symmetry,
   helical_rot_deg = helical_rot_deg,
   helical_trans_z_angstrom = helical_trans_z_angstrom,
   op_max = op_max,
   two_fold_along_x = two_fold_along_x,
   include_helical_symmetry = \
       params.reconstruction_symmetry.include_helical_symmetry,
   max_helical_ops_to_check = \
       params.reconstruction_symmetry.max_helical_ops_to_check,
   require_helical_or_point_group_symmetry = \
      params.reconstruction_symmetry.require_helical_or_point_group_symmetry,
   out = out)

  # Generate helical symmetry from map if necessary
  if symmetry.lower() == 'helical' or (
      symmetry.lower() in ['all', 'any'] and
      params.reconstruction_symmetry.include_helical_symmetry):
     if helical_rot_deg is None or helical_trans_z_angstrom is None:

     # returns ncs for symmetry_center at symmetry_center
      ncs_object, helical_rot_deg, helical_trans_z_angstrom = \
         find_helical_symmetry(params = params,
        symmetry_center = symmetry_center,
        map_data = map_data,
        crystal_symmetry = crystal_symmetry, out = out)

      if ncs_object:
        ncs_name = "Type: Helical %5.2f deg  %6.2f Z-trans " %(
          helical_rot_deg, helical_trans_z_angstrom)
        ncs_object.set_ncs_name(ncs_name)
        ncs_list.append(ncs_object)

  if symmetry_center and tuple(symmetry_center) !=  (0, 0, 0, ):
    print("Offsetting NCS center by (%.2f, %.2f, %.2f) A " %(tuple(symmetry_center)), file = out)
    new_list = []
    for ncs_obj in ncs_list:
      new_list.append(ncs_obj.coordinate_offset(coordinate_offset = symmetry_center))
    ncs_list = new_list

  if (require_helical_or_point_group_symmetry):
    for ncs_obj in ncs_list:
      assert ncs_obj.is_helical_along_z() or ncs_obj.is_point_group_symmetry()

  return ncs_list

def find_helical_symmetry(params = None,
        symmetry_center = None,
        map_data = None,
        crystal_symmetry = None,
        max_z_to_test = 2, max_peaks_to_score = 5, out = sys.stdout):

  """Find helical symmetry in a map"""

  params = deepcopy(params) # so changing them does not go back
  if not params.crystal_info.resolution:
    from cctbx.maptbx import d_min_from_map
    params.crystal_info.resolution = d_min_from_map(
      map_data, crystal_symmetry.unit_cell(), resolution_factor = 1./4.)


  if str(params.reconstruction_symmetry.score_basis) == 'None':
    params.reconstruction_symmetry.score_basis = 'cc'
  if params.reconstruction_symmetry.smallest_object is None:
    params.reconstruction_symmetry.smallest_object = \
      5*params.crystal_info.resolution

  print("\nLooking for helical symmetry with center at (%.2f, %.2f, %.2f) A\n" %(
      tuple(symmetry_center)), file = out)
  print("\nFinding likely translations along z...", file = out)

  map_coeffs, dummy = get_f_phases_from_map(map_data = map_data,
       crystal_symmetry = crystal_symmetry,
       d_min = params.crystal_info.resolution,
       return_as_map_coeffs = True,
       out = out)
  f_array, phases = map_coeffs_as_fp_phi(map_coeffs)
  from cctbx.maptbx.refine_sharpening import quasi_normalize_structure_factors
  (d_max, d_min) = f_array.d_max_min()
  f_array.setup_binner(d_max = d_max, d_min = d_min, n_bins = 20)
  normalized = quasi_normalize_structure_factors(
          f_array, set_to_minimum = 0.01)

  # Now look along c* and get peaks
  zero_index_a = 0  # look along c*
  zero_index_b = 1
  c_star_list = get_c_star_list(f_array = normalized,
     zero_index_a = zero_index_a,
     zero_index_b = zero_index_b)

  likely_z_translations = get_helical_trans_z_angstrom(params = params,
     c_star_list = c_star_list,
    crystal_symmetry = crystal_symmetry,
    out = out)

  # Now for each z...get the rotation. Try +/- z and try rotations
  abc = crystal_symmetry.unit_cell().parameters()[:3]
  max_z = max(abc[0], abc[1])
  if params.reconstruction_symmetry.max_helical_rotations_to_check:
    min_delta_rot = 360/max(1,
    params.reconstruction_symmetry.max_helical_rotations_to_check)
  else:
    min_delta_rot = 0.01

  delta_rot = max(min_delta_rot,
     (180./3.141659)*params.crystal_info.resolution/max_z)

  delta_z = params.crystal_info.resolution/4.
  print("\nOptimizing helical paramers:", file = out)
  print("\n Rot   Trans  Score    CC", file = out)

  n_rotations = int(0.5+360/delta_rot)
  rotations = []
  for k in range(1, n_rotations):
    helical_rot_deg = k*delta_rot
    if helical_rot_deg > 180: helical_rot_deg = helical_rot_deg-360
    rotations.append(helical_rot_deg)

  done = False
  n_try = 0

  score_list = []
  quick = params.control.quick
  best_ncs_cc = None
  best_ncs_obj = None
  best_score = None
  best_helical_trans_z_angstrom = None
  best_helical_rot_deg = None
  import math
  for helical_trans_z_angstrom  in likely_z_translations[:max_z_to_test]:
    n_try+= 1
    if done: break
    for helical_rot_deg in rotations:
      if done: break
      if abs(helical_trans_z_angstrom)+abs(math.sin(
        helical_rot_deg*(3.14159/180.))*max_z) < \
         params.reconstruction_symmetry.smallest_object:
         continue  # this is identity
      new_ncs_obj, ncs_cc, ncs_score, \
          new_helical_trans_z_angstrom, new_helical_rot_deg = try_helical_params(
        optimize = 0,
        helical_rot_deg = helical_rot_deg,
        helical_trans_z_angstrom = helical_trans_z_angstrom,
        params = params,
        map_data = map_data,
        map_symmetry_center = symmetry_center, # should not be needed XXX
        symmetry_center = symmetry_center,
        crystal_symmetry = crystal_symmetry,
        out = null_out())

      if not ncs_score: continue
      if ncs_cc> 0.95 or \
        ncs_cc > 2*params.reconstruction_symmetry.min_ncs_cc or \
        (quick and (ncs_cc> 0.90 or
         ncs_cc > 1.5*params.reconstruction_symmetry.min_ncs_cc)):
        print(" %.2f   %.2f   %.2f   %.2f (ok to go on)" %(
           new_helical_rot_deg, new_helical_trans_z_angstrom,
           ncs_score, ncs_cc), file = out)
        done = True
      if params.control.verbose:
        print(" %.2f   %.2f   %.2f   %.2f" %(
          new_helical_rot_deg, new_helical_trans_z_angstrom, ncs_score, ncs_cc), file = out)
      score_list.append(
       [ncs_score, ncs_cc, new_helical_rot_deg, new_helical_trans_z_angstrom])
    score_list.sort(key=itemgetter(0))
    score_list.reverse()
    done = False
    for ncs_score, ncs_cc, helical_rot_deg, helical_trans_z_angstrom in \
       score_list[:max_peaks_to_score]:
      if done: continue
      # rescore and optimize:

      new_ncs_obj, ncs_cc, ncs_score, \
          new_helical_trans_z_angstrom, new_helical_rot_deg = try_helical_params(
        params = params,
        best_score = best_score,
        helical_rot_deg = helical_rot_deg,
        helical_trans_z_angstrom = helical_trans_z_angstrom,
        delta_z = delta_z/2.,
        delta_rot = delta_rot/2.,
        map_data = map_data,
        map_symmetry_center = symmetry_center,
        symmetry_center = symmetry_center,
        crystal_symmetry = crystal_symmetry,
        out = null_out())
      if not ncs_score or ncs_score <0:
        continue
      if best_score is None or ncs_score>best_score:
        best_ncs_cc = ncs_cc
        best_ncs_obj = new_ncs_obj
        best_score = ncs_score
        best_helical_trans_z_angstrom = new_helical_trans_z_angstrom
        best_helical_rot_deg = new_helical_rot_deg


      # after trying out a range quit if good enough
      if best_ncs_cc> 0.90 or \
          best_ncs_cc > 1.5*params.reconstruction_symmetry.min_ncs_cc or \
          ( (quick or n_try>1) and \
             ncs_cc>= params.reconstruction_symmetry.min_ncs_cc):
        print(" %.2f   %.2f   %.2f   %.2f (high enough to go on)" %(
           best_helical_rot_deg, best_helical_trans_z_angstrom,
           best_score, best_ncs_cc), file = out)


        done = True

  # Optimize one more time trying fractional values, but only if
  #  that makes the delta less than the resolution

  print("\nTrying fraction of rot/trans", file = out)
  for iter in [0, 1]:
    if not best_helical_rot_deg: continue
    for ifract in range(2, 11):
      if iter == 0:  # try fractional
        test_helical_rot_deg = best_helical_rot_deg/ifract
        test_helical_trans_z_angstrom = best_helical_trans_z_angstrom/ifract
        if test_helical_trans_z_angstrom > params.crystal_info.resolution*1.1:
          continue # skip it...would have been a peak if ok
      else: # iter >0
        if ifract > 0:
           continue # skip these
        else: # optimize current best
          test_helical_rot_deg = best_helical_rot_deg
          test_helical_trans_z_angstrom = best_helical_trans_z_angstrom
      new_ncs_obj, new_ncs_cc, new_ncs_score, \
            new_helical_trans_z_angstrom, new_helical_rot_deg = \
        try_helical_params(
          params = params,
          helical_rot_deg = test_helical_rot_deg,
          helical_trans_z_angstrom = test_helical_trans_z_angstrom,
          delta_z = delta_z/2.,
          delta_rot = delta_rot/2.,
          map_data = map_data,
          map_symmetry_center = symmetry_center,
          symmetry_center = symmetry_center,
          crystal_symmetry = crystal_symmetry,
          out = out)
      if new_ncs_score is not None:
        new_ncs_score = new_ncs_score*\
           params.reconstruction_symmetry.scale_weight_fractional_translation
           # give slight weight to smaller

        if best_score is None or new_ncs_score > best_score:
          best_ncs_cc = new_ncs_cc
          best_ncs_obj = new_ncs_obj
          best_score = new_ncs_score
          best_helical_trans_z_angstrom = new_helical_trans_z_angstrom
          best_helical_rot_deg = new_helical_rot_deg
          print(" %.2f   %.2f   %.2f   %.2f (improved fractions)" %(
             best_helical_rot_deg, best_helical_trans_z_angstrom,
             best_score, best_ncs_cc), file = out)
        else:
          print(" %.2f   %.2f   %.2f   %.2f (worse with fractions)" %(
             new_helical_rot_deg, new_helical_trans_z_angstrom,
             new_ncs_score, new_ncs_cc), file = out)


  # Optimize one more time trying multiples of values to get better param
  imult = int(0.5+
     0.33*max_z/params.reconstruction_symmetry.max_helical_ops_to_check)


  working_ncs_cc = best_ncs_cc
  working_ncs_obj = best_ncs_obj
  working_score = None
  if best_helical_rot_deg:
    working_helical_rot_deg = best_helical_rot_deg*imult
    working_helical_trans_z_angstrom = best_helical_trans_z_angstrom*imult
  else:
    working_helical_rot_deg = None
    working_helical_trans_z_angstrom = None
    imult = 0

  if imult > 1:
    print("\nTrying %sx multiples of rot/trans" %(imult), file = out)
    improved = False
    for iter in [1, 2, 3]:
      if iter > 1 and not improved: break
      improved = False

      new_ncs_obj, new_ncs_cc, new_ncs_score, \
            new_helical_trans_z_angstrom, new_helical_rot_deg = \
        try_helical_params(
          params = params,
          helical_rot_deg = working_helical_rot_deg,
          helical_trans_z_angstrom = working_helical_trans_z_angstrom,
          delta_z = delta_z/(2.*iter),
          delta_rot = delta_rot/(2.*iter),
          map_data = map_data,
          map_symmetry_center = symmetry_center,
          symmetry_center = symmetry_center,
          crystal_symmetry = crystal_symmetry,
          out = out)
      if new_ncs_score is not None:
        if working_score is None or new_ncs_score > working_score:
          working_ncs_cc = new_ncs_cc
          working_ncs_obj = new_ncs_obj
          working_score = new_ncs_score
          working_helical_trans_z_angstrom = new_helical_trans_z_angstrom
          working_helical_rot_deg = new_helical_rot_deg
          print(" %.2f   %.2f   %.2f   %.2f (Scoring for multiple)" %(
               working_helical_rot_deg, working_helical_trans_z_angstrom,
               working_score, working_ncs_cc), file = out)

      #  and rescore with this:
      working_helical_rot_deg = working_helical_rot_deg/imult
      working_helical_trans_z_angstrom = working_helical_trans_z_angstrom/imult
      for iter in [1, 2, 3]:
        new_ncs_obj, new_ncs_cc, new_ncs_score, \
            new_helical_trans_z_angstrom, new_helical_rot_deg = \
        try_helical_params(
          params = params,
          helical_rot_deg = working_helical_rot_deg,
          helical_trans_z_angstrom = working_helical_trans_z_angstrom,
          delta_z = delta_z/(2.*iter),
          delta_rot = delta_rot/(2.*iter),
          map_data = map_data,
          map_symmetry_center = symmetry_center,
          symmetry_center = symmetry_center,
          crystal_symmetry = crystal_symmetry,
          out = null_out())
        if new_ncs_score is not None:
          working_ncs_obj, working_ncs_cc, working_ncs_score, \
            working_helical_trans_z_angstrom, working_helical_rot_deg = \
            new_ncs_obj, new_ncs_cc, new_ncs_score, \
            new_helical_trans_z_angstrom, new_helical_rot_deg
          if best_score is None or working_ncs_score > best_score:
            if params.control.verbose:
              print("\nTaking parameters from multiples", file = out)
            best_ncs_cc = working_ncs_cc
            best_ncs_obj = working_ncs_obj
            best_score = working_ncs_score
            best_helical_trans_z_angstrom = working_helical_trans_z_angstrom
            best_helical_rot_deg = working_helical_rot_deg
            print(" %.2f   %.2f   %.2f   %.2f (improved by multiples)" %(
               best_helical_rot_deg, best_helical_trans_z_angstrom,
               best_score, best_ncs_cc), file = out)
            improved = True

            working_ncs_cc = best_ncs_cc
            working_ncs_obj = best_ncs_obj
            working_score = None
            working_helical_trans_z_angstrom = best_helical_trans_z_angstrom*imult
            working_helical_rot_deg = best_helical_rot_deg*imult


  if best_helical_rot_deg and best_helical_trans_z_angstrom and best_score and best_ncs_cc:
    # Check to make sure there is no overlap

    print(" %.2f   %.2f   %.2f   %.2f (Final)" %(
     best_helical_rot_deg, best_helical_trans_z_angstrom, \
         best_score, best_ncs_cc), file = out)

    from mmtbx.ncs.ncs import get_helical_symmetry

    ncs_object = get_helical_symmetry(
       helical_rot_deg = best_helical_rot_deg,
       helical_trans_z_angstrom = best_helical_trans_z_angstrom,
       max_ops = params.reconstruction_symmetry.max_helical_ops_to_check)
  else:
    ncs_object = None

  return ncs_object, best_helical_rot_deg, best_helical_trans_z_angstrom

def try_helical_params(
    optimize = None,
    best_score = None,
    delta_z = None,
    delta_rot = None,
    helical_rot_deg = None,
    helical_trans_z_angstrom = None,
    params = None,
    map_data = None,
    map_symmetry_center = None,
    symmetry_center = None,
    crystal_symmetry = None,
    out = sys.stdout):

  """Try to find helical symmetry with specified parameters"""
  if delta_z is None or delta_rot is None:
    assert not optimize
  local_params = deepcopy(params)
  local_params.reconstruction_symmetry.min_ncs_cc = -100
  local_params.reconstruction_symmetry.n_rescore = 0
  local_params.reconstruction_symmetry.symmetry = 'helical'
  local_params.reconstruction_symmetry.helical_rot_deg = helical_rot_deg
  local_params.reconstruction_symmetry.helical_trans_z_angstrom = \
     helical_trans_z_angstrom
  abc = crystal_symmetry.unit_cell().parameters()[:3]
  max_z = max(abc[0], abc[1])
  import math
  if abs(helical_trans_z_angstrom)+abs(math.sin(
        helical_rot_deg*(3.14159/180.))*max_z) < \
         params.reconstruction_symmetry.smallest_object:
         return None, None, None, \
                None, None

  best_helical_trans_z_angstrom, best_helical_rot_deg = \
     helical_trans_z_angstrom, helical_rot_deg
  best_ncs_score = best_score
  best_ncs_obj = None
  best_ncs_cc = None

  test_ncs_obj, test_ncs_cc, test_ncs_score = get_ncs_from_map(params = local_params,
    map_data = map_data,
    map_symmetry_center = symmetry_center,
    symmetry_center = symmetry_center,
    crystal_symmetry = crystal_symmetry,
    out = null_out())
  if params.reconstruction_symmetry.score_basis == 'cc':
     test_ncs_score = test_ncs_cc

  if best_ncs_score is None or test_ncs_score>best_ncs_score:
    best_ncs_score = test_ncs_score
    best_ncs_cc = test_ncs_cc
    best_ncs_obj = test_ncs_obj

  if optimize is None:
    optimize = params.reconstruction_symmetry.max_helical_optimizations


  # save in case we need to go back
  working_helical_trans_z_angstrom, working_helical_rot_deg = \
     helical_trans_z_angstrom, helical_rot_deg
  working_ncs_score = best_ncs_score
  working_ncs_cc = best_ncs_cc
  working_ncs_obj = best_ncs_obj


  if optimize and (best_score is None or best_ncs_score > best_score):
    # try with few to many operators..
    if params.control.verbose:
      print("\nOptimizing by varying number of operators", file = out)
      print("Starting score: %.2f" %(working_ncs_score), file = out)
    for k in range(optimize):
     local_params.reconstruction_symmetry.max_helical_ops_to_check = min(k+1,
      params.reconstruction_symmetry.max_helical_ops_to_check)
     best_score = None # start over for each number of operators
     for i in [0, -1, 1]:
       for j in [0, -1, 1]:
        new_ncs_obj, new_ncs_cc, new_ncs_score, \
            new_helical_trans_z_angstrom, new_helical_rot_deg = try_helical_params(
          optimize = False,
          helical_rot_deg = max(0.1, best_helical_rot_deg+i*delta_rot),
          helical_trans_z_angstrom = max(0.01, best_helical_trans_z_angstrom+j*delta_z),
          delta_z = delta_z,
          delta_rot = delta_rot,
          params = params,
          map_data = map_data,
          map_symmetry_center = symmetry_center,
          symmetry_center = symmetry_center,
          crystal_symmetry = crystal_symmetry,
          out = out)
        if new_ncs_score and new_ncs_score>0 and (
            best_score is None or new_ncs_score>best_score):
          if params.control.verbose:
            print("Working score for %s ops: %.2f" %(
              local_params.reconstruction_symmetry.max_helical_ops_to_check,
              new_ncs_score), file = out)
          best_score = new_ncs_score
          best_ncs_score = new_ncs_score
          best_helical_trans_z_angstrom = new_helical_trans_z_angstrom
          best_helical_rot_deg = new_helical_rot_deg
          best_ncs_obj = new_ncs_obj
          best_ncs_cc = new_ncs_cc
     delta_rot = delta_rot/2
     delta_z = delta_z/2


    #rescore with what we now have (best values) and compare with working
    local_params.reconstruction_symmetry.max_helical_ops_to_check = \
       params.reconstruction_symmetry.max_helical_ops_to_check
    if params.control.verbose:
      print("Rescoring with original number of operators (%s)" %(
        local_params.reconstruction_symmetry.max_helical_ops_to_check), file = out)
    local_params.reconstruction_symmetry.helical_rot_deg = best_helical_rot_deg
    local_params.reconstruction_symmetry.helical_trans_z_angstrom = \
      best_helical_trans_z_angstrom
    best_ncs_obj, best_ncs_cc, best_ncs_score = get_ncs_from_map(
         params = local_params,
      map_data = map_data,
      map_symmetry_center = symmetry_center,
      symmetry_center = symmetry_center,
      crystal_symmetry = crystal_symmetry,
      out = null_out())
    if params.reconstruction_symmetry.score_basis == 'cc':
      best_ncs_score = best_ncs_cc

    # now take it if better
    if best_ncs_cc and best_ncs_cc>working_ncs_cc:
      if params.control.verbose:
        print("Using optimized values (%.2f > %.2f)" %(
         best_ncs_cc, working_ncs_cc), file = out)
      # keep these (best)
    else:
      if params.control.verbose:
        print("Rejecting optimized values (%.2f <=  %.2f)" %(
         best_ncs_cc, working_ncs_cc), file = out)

      # resture working values
      best_helical_trans_z_angstrom, best_helical_rot_deg = \
         working_helical_trans_z_angstrom, working_helical_rot_deg
      best_ncs_score = working_ncs_score
      best_ncs_obj = working_ncs_obj
      best_ncs_cc = working_ncs_cc

  if params.reconstruction_symmetry.score_basis == 'cc':
     best_ncs_score = best_ncs_cc

  return best_ncs_obj, best_ncs_cc, best_ncs_score, \
            best_helical_trans_z_angstrom, best_helical_rot_deg


def get_d_and_value_list(c_star_list):
  """Convert c_star_list tuples to separate lists of d and value.
  Reject (skip) values < maximum value/1000."""
  d_list = []
  from scitbx.array_family import flex
  value_list = flex.double()
  for hkl, value, d in c_star_list:
    d_list.append(d)
    value_list.append(value)
  max_value = value_list.min_max_mean().max
  if value_list.size()>3:
    max_value = value_list[2]
  new_d_list = []
  new_value_list = []
  for d, value in zip(d_list, value_list):
    if value < max_value/1000: # reject those that are really zero
       continue
    new_d_list.append(d)
    new_value_list.append(value)
  return new_d_list, new_value_list

def get_helical_trans_z_angstrom(params = None,
  c_star_list = None, crystal_symmetry = None,
  minimum_ratio = 2.,
  max_first_peak = 1, out = sys.stdout):

  """Find helical translation along z based on
  highest-resolution peak along c*...guess it is n*z_translation
  where n is small"""

  max_z = flex.double(crystal_symmetry.unit_cell().parameters()[:3]).min_max_mean().max

  if params.control.verbose:
    print("Values along c*: ", file = out)
  d_list, value_list = get_d_and_value_list(c_star_list)

  for d, value in zip(d_list, value_list):
    if params.control.verbose:
      print(" %.2f A : %.2f " %(d, value), file = out)

  d_list, value_list = get_max_min_list(
     d_list = d_list, value_list = value_list, minimum_ratio = 1.0)

  d_list, value_list = get_max_min_list(
     d_list = d_list, value_list = value_list, minimum_ratio = 2.0,
    maximum_only = True)
  sort_list = []
  for d, value in zip(d_list, value_list):
    sort_list.append([value, d])
  sort_list.sort(key=itemgetter(0))
  sort_list.reverse()
  likely_z_translations = []
  dis_min = params.crystal_info.resolution/4
  for value, d in sort_list:
    likely_z_translations_local = []
    for i in range(1, max_first_peak+1):
      z = d/i
      if z > max_z: continue
      if z < dis_min: continue # no way
      if is_close_to_any(z = z, others = likely_z_translations,
        dis_min = dis_min): continue

      likely_z_translations_local.append(z)
      likely_z_translations.append(z)

  print("\nMaximal values along c* and likely z-translations: ", file = out)
  for z in likely_z_translations:
    print(" %.2f A " %(z), end = ' ', file = out)
  print(file = out)
  return likely_z_translations

def small_deltas(z_translations, dis_min = None):
  """Find unique differences between sequential values of z_translations"""
  delta_list = []
  for z, z1 in zip(z_translations, z_translations[1:]):
     delta = abs(z-z1)
     if not is_close_to_any(z = delta, others = z_translations+delta_list,
        dis_min = dis_min):
       delta_list.append(abs(delta))
  return delta_list

def is_close_to_any(z = None, others = None,
    dis_min = None):
  """Return True if z is within dis_min of any others"""
  for x in others:
    if abs(x-z)<dis_min:
       return True
  return False

def get_max_min_list(d_list = None, value_list = None,
   minimum_ratio = None, maximum_only = False):
  """Find local maxima and minima in value_list"""
  max_min_list = []
  max_min_d_list = []
  for value_prev, d_prev, \
       value, d, \
       value_next, d_next in zip(
         value_list+[0, 0],
         d_list+[0, 0],
         [0]+value_list+[0],
         [0]+d_list+[0],
         [0, 0]+value_list,
         [0, 0]+d_list):

    if d and ( value >=   value_prev *minimum_ratio) and (
       value >=  value_next*minimum_ratio):  # local max
      max_min_list.append(value)
      max_min_d_list.append(d)
    if (not maximum_only) and d and ( value <=   value_prev ) and (
       value <=  value_next):  # local min
      max_min_list.append(value)
      max_min_d_list.append(d)
  return max_min_d_list, max_min_list

def get_c_star_list(f_array = None,
   zero_index_a = 0, zero_index_b = 1, zero_index_c = 2):
  """Create c_star_list of tuples (indices, value, d)"""
  c_star_list = []
  for value, (indices, d) in zip(f_array.data(),
     f_array.d_spacings()):
    if indices[zero_index_a] == 0 and indices[zero_index_b] == 0 and \
      indices[zero_index_c] >= 4:
      c_star_list.append([tuple(indices), value, d])
  c_star_list.sort()
  return c_star_list

def get_params_from_args(args):
  """Get a params object from a list of args for map_file, seq_file, pdb_file
  and ncs_file and any other args"""
  command_line = iotbx.phil.process_command_line_with_files(
    map_file_def = "input_files.map_file",
    seq_file_def = "input_files.seq_file",
    pdb_file_def = "input_files.pdb_in",
    ncs_file_def = "input_files.ncs_file",
    args = args,
    master_phil = master_phil)

  return command_line.work.extract()


def get_mask_around_molecule(map_data = None,
        wang_radius = None,
        buffer_radius = None,
        force_buffer_radius = None,
        return_masked_fraction = False,
        minimum_fraction_of_max = 0.01,
        solvent_content = None,
        solvent_content_iterations = None,
        crystal_symmetry = None, out = sys.stdout):
  """Use iterated solvent fraction tool to identify mask around molecule"""
  try:
    from phenix.autosol.map_to_model import iterated_solvent_fraction
    solvent_fraction, mask = iterated_solvent_fraction(
      crystal_symmetry = crystal_symmetry,
      wang_radius = wang_radius,
      solvent_content = solvent_content,
      solvent_content_iterations = solvent_content_iterations,
      map_as_double = map_data,
      out = out)
  except Exception as e:
    print("No mask obtained...", file = out)
    return None, None

  if not mask:
    print("No mask obtained...", file = out)
    return None, None

  # Now expand the mask to increase molecular region
  expand_size = estimate_expand_size(
       crystal_symmetry = crystal_symmetry,
       map_data = map_data,
       expand_target = buffer_radius,
       minimum_expand_size = 0,
       out = out)

  print("Target mask expand size is %d based on buffer_radius of %7.1f A" %(
     expand_size, buffer_radius), file = out)

  co, sorted_by_volume, min_b, max_b = get_co(map_data = mask,
     threshold = 0.5, wrapping = False)
  if len(sorted_by_volume)<2:
    print ("\nSkipping expansion as no space is available\n", file = out)
    return None, None

  masked_fraction = sorted_by_volume[1][0]/mask.size()

  if expand_size <=  0:
    expanded_fraction = masked_fraction
  else: # Try to get expanded fraction < 0.5*(masked_fraction+1)
    upper_limit = 0.5*(masked_fraction+1)

    bool_region_mask = co.expand_mask(id_to_expand = sorted_by_volume[1][1],
         expand_size = expand_size)
    s = (bool_region_mask == True)
    expanded_fraction = s.count(True)/s.size()
    if expanded_fraction>upper_limit:
      amount_too_big = max(1.e-10, expanded_fraction-upper_limit)/max(1.e-10,
         expanded_fraction-masked_fraction)**0.667
      # cut back
      original_expand_size = expand_size
      expand_size = max(1, int(0.5+expand_size * min(1, max(0, (1-amount_too_big)))))
      if expand_size !=  original_expand_size:

        print ("\nCutting back expand size to try and get "+
           "fraction < about %.2f . New expand_size: %s" %(
          upper_limit, expand_size), file = out)
        bool_region_mask = co.expand_mask(id_to_expand = sorted_by_volume[1][1],
           expand_size = expand_size)
        s = (bool_region_mask == True)
        expanded_fraction = s.count(True)/s.size()
    if expanded_fraction > 0.9999:
          print ("\nSkipping expansion as no space is available\n", file = out)
          return None, None
  print("\nLargest masked region before buffering: %7.2f" %(masked_fraction),
      file = out)
  print("\nLargest masked region after buffering: %7.2f" %(expanded_fraction),
     file = out)

  if solvent_content and (not force_buffer_radius):
    delta_as_is = abs(solvent_content- (1-masked_fraction))
    delta_expanded = abs(solvent_content- (1-expanded_fraction))
    if delta_expanded > delta_as_is:
      # already there
      expand_size = 0
      print ("Setting expand size to zero as masked fraction already ",
         "close to solvent_content", file = out)

  s = None
  minimum_size = sorted_by_volume[1][0] * minimum_fraction_of_max
  if expand_size == 0:
    result = co.result()
  else:
    result = None

  for v1, i1 in sorted_by_volume[1:]:
    if v1 < minimum_size: break
    if expand_size > 0:
      bool_region_mask = co.expand_mask(
        id_to_expand = i1, expand_size = expand_size)
    else:
      bool_region_mask = (result == i1)
    if s is None:
      s = (bool_region_mask == True)
    else:
      s |=  (bool_region_mask == True)
  mask.set_selected(s, 1)
  mask.set_selected(~s, 0)
  masked_fraction = mask.count(1)/mask.size()
  print("Masked fraction after buffering:  %7.2f" %(masked_fraction), file = out)
  if return_masked_fraction:
    return mask.as_double(), 1-masked_fraction
  else: # usual return solvent fraction estimate
    return mask.as_double(), solvent_fraction  # This is solvent fraction est.

def get_mean_in_and_out(sel = None,
    map_data = None,
    verbose = False,
    out = sys.stdout):
  """Get mean value of map inside and outside of selection sel"""
  mean_value_in, fraction_in = get_mean_in_or_out(sel = sel,
    map_data = map_data,
    out = out)

  mean_value_out, fraction_out = get_mean_in_or_out(sel =  ~sel,
    map_data = map_data,
    out = out)

  if mean_value_out is None:
    mean_value_out = mean_value_in
  if mean_value_in is None:
    mean_value_in = mean_value_out
  if verbose:
    print("\nMean inside mask: %7.2f  Outside mask: %7.2f  Fraction in: %7.2f" %(
     mean_value_in, mean_value_out, fraction_in), file = out)
  return mean_value_in, mean_value_out, fraction_in

def get_mean_in_or_out(sel = None,
    map_data = None,
    out = sys.stdout):
  """Get mean value of map and fraction_in_or_out,
     either inside or outside of selection sel"""
  masked_map = map_data.deep_copy()
  if sel.count(False) != 0:
    masked_map.as_1d().set_selected(~sel.as_1d(), 0)
  mean_after_zeroing_in_or_out = masked_map.as_1d().min_max_mean().mean
  if sel.count(True) !=  0:
    masked_map.as_1d().set_selected(sel.as_1d(), 1)
  fraction_in_or_out = masked_map.as_1d().min_max_mean().mean
  if fraction_in_or_out >1.e-10:
    mean_value = mean_after_zeroing_in_or_out/fraction_in_or_out
  else:
    mean_value = None

  return mean_value, fraction_in_or_out

def apply_soft_mask(map_data = None,
          mask_data = None,
          rad_smooth = None,
          crystal_symmetry = None,
          set_outside_to_mean_inside = False,
          set_mean_to_zero = False,
          threshold = 0.5,
          verbose = False,
          out = sys.stdout):

  """Apply a soft mask based on mask_data to map_data.
  set value outside mask == mean value inside mask or mean value outside mask"""

  original_mask_data = mask_data.deep_copy()

  smoothed_mask_data = smooth_mask_data(mask_data = mask_data,
    crystal_symmetry = crystal_symmetry,
    threshold = threshold,
    rad_smooth = rad_smooth)

  masked_map = apply_mask_to_map(mask_data = original_mask_data,
     smoothed_mask_data = smoothed_mask_data,
     map_data = map_data,
     set_outside_to_mean_inside = set_outside_to_mean_inside,
     set_mean_to_zero = set_mean_to_zero,
     threshold = threshold,
     verbose = verbose,
     out = out)
  return masked_map, smoothed_mask_data

def smooth_one_map(map_data, crystal_symmetry = None, smoothing_radius = None,
     non_negative = False, method = 'exp'):
  """Smooth a map with radius smoothing radius. Returns new map object"""
  from cctbx.maptbx import smooth_map, unpad_in_place
  unpad_in_place(map = map_data.deep_copy())
  smoothed_map = smooth_map(
      map              = map_data,
      crystal_symmetry = crystal_symmetry,
      rad_smooth       = smoothing_radius,
      non_negative     = non_negative,
      method = method)
  return smoothed_map

def get_smoothed_cc_map(map_data_1, map_data_2,
   crystal_symmetry = None, weighting_radius = None,
   method = 'top_hat'):
  """Calculate smoothed version of map correlation (local map correlation"""
  avg1 = map_data_1.as_1d().min_max_mean().mean
  avg2 = map_data_2.as_1d().min_max_mean().mean
  avg_product_map = smooth_one_map((map_data_1 - avg1) * (map_data_2 - avg2),
     crystal_symmetry = crystal_symmetry,
      smoothing_radius = weighting_radius,
      non_negative = True, method = method)

  sm1 = smooth_one_map (flex.pow2(map_data_1 - avg1),
      crystal_symmetry = crystal_symmetry,
      smoothing_radius = weighting_radius,
      non_negative = True, method = method)
  sm2 = smooth_one_map(flex.pow2(map_data_2 - avg2),
      crystal_symmetry = crystal_symmetry,
      smoothing_radius = weighting_radius,
      non_negative = True, method = method)
  sm = flex.sqrt(sm1) * flex.sqrt(sm2)
  sm.as_1d().set_selected(sm.as_1d() < 1.e-20, 1.e-20)
  cc_map = avg_product_map/sm
  cc_map.as_1d().set_selected(cc_map.as_1d() < 0, 0)
  cc_map.as_1d().set_selected(cc_map.as_1d() > 1, 1)

  return cc_map


def smooth_mask_data(mask_data = None,
    crystal_symmetry = None,
    threshold = None,
    rad_smooth = None):

  """Smooth a mask in place"""
  if threshold is None:
    threshold = mask_data.as_1d().min_max_mean().max * 0.5

  # Smooth a mask in place. First make it a binary mask
  s = mask_data > threshold  # s marks inside mask
  mask_data = mask_data.set_selected(~s, 0)  # outside mask == 0
  mask_data = mask_data.set_selected( s, 1)
  if mask_data.count(1)  and mask_data.count(0): # something to do
    maptbx.unpad_in_place(map = mask_data)
    mask_data = maptbx.smooth_map(
      map              = mask_data,
      crystal_symmetry = crystal_symmetry,
      rad_smooth       = rad_smooth)

    # Make sure that mask_data max value is now 1, scale if not
    max_mask_data_value = mask_data.as_1d().min_max_mean().max
    if max_mask_data_value > 1.e-30 and max_mask_data_value!= 1.0:
      mask_data = mask_data*(1./max_mask_data_value)
  else:
    pass
  return mask_data


def apply_mask_to_map(mask_data = None,
    smoothed_mask_data = None,
    set_outside_to_mean_inside = None,
    set_mean_to_zero = None,
    map_data = None,
    threshold = 0.5,
    verbose = None,
    out = sys.stdout):

  """Appy a mask to map. Deprecated. Use map_manager instead"""
  if smoothed_mask_data and not mask_data:
    mask_data = smoothed_mask_data
  elif mask_data and not smoothed_mask_data:
    smoothed_mask_data = mask_data

  s = mask_data > threshold  # s marks inside mask
  # get mean inside or outside mask
  if verbose:
    print("\nStarting map values inside and outside mask:", file = out)

  mean_value_in, mean_value_out, fraction_in = get_mean_in_and_out(sel = s,
    verbose = verbose, map_data = map_data, out = out)

  if verbose:
    print("\nMask inside and outside values", file = out)

  mask_mean_value_in, mask_mean_value_out, mask_fraction_in = get_mean_in_and_out(
      sel = s, map_data = mask_data, verbose = verbose, out = out)

  if verbose:
    print("\nSmoothed mask inside and outside values", file = out)

  smoothed_mean_value_in, smoothed_mean_value_out, smoothed_fraction_in = \
     get_mean_in_and_out(sel = s, map_data = smoothed_mask_data,
       verbose = verbose, out = out)

  # Now replace value outside mask with mean_value, value inside with current,
  #   smoothly going from one to the other based on mask_data

  # set_to_mean will be a constant map with value equal to inside or outside

  if (set_outside_to_mean_inside is False):
    target_value_for_outside = 0
  elif set_outside_to_mean_inside or mean_value_out is None:
    target_value_for_outside = mean_value_in
    if verbose:
      print("Setting value outside mask to mean inside (%.2f)" %(
      target_value_for_outside), file = out)
  else:
    target_value_for_outside = mean_value_out
    if verbose:
      print("Setting value outside mask to mean outside (%.2f)" %(
        target_value_for_outside), file = out)
  set_to_mean = mask_data.deep_copy()
  ss = set_to_mean > -1.e+30 # select everything
  set_to_mean.set_selected(ss, target_value_for_outside)

  masked_map =  (map_data * smoothed_mask_data ) +  (set_to_mean * (1-smoothed_mask_data))

  if set_mean_to_zero:  # remove average
    masked_map = masked_map - masked_map.as_1d().min_max_mean().mean

  if verbose:
    print("\nFinal mean value inside and outside mask:", file = out)
  mean_value_in, mean_value_out, fraction_in = get_mean_in_and_out(sel = s,
    map_data = masked_map, verbose = verbose, out = out)

  return masked_map

def estimate_expand_size(
       crystal_symmetry = None,
       map_data = None,
       expand_target = None,
       minimum_expand_size = 1,
       out = sys.stdout):
    """Estimate value of expand_size (dimension in grid units to expand
    masked region around a molecule to have a good chance of keeping most of
    the molecule inside the mask.)
    Uses minimum_expand_size of 1, typical value of  1"""

    if not expand_target:
     return minimum_expand_size
    abc = crystal_symmetry.unit_cell().parameters()[:3]
    N_ = map_data.all()
    nn = 0.
    for i in range(3):
      delta = abc[i]/N_[i]
      nn+= expand_target/delta
    nn = max(1, int(0.5+nn/3.))
    print("Expand size (grid units): %d (about %4.1f A) " %(
      nn, nn*abc[0]/N_[0]), file = out)
    return max(minimum_expand_size, nn)

def get_max_z_range_for_helical_symmetry(params, out = sys.stdout):
  """Identify maximum z range to consider for helical symmetry based on
  params.map_modification.restrict_z_turns_for_helical_symmetry"""

  if not params.input_files.ncs_file: return
  ncs_obj, dummy_tracking_data = get_ncs(params = params, out = out)
  if not ncs_obj.is_helical_along_z(): return
  if params.map_modification.restrict_z_distance_for_helical_symmetry:  #take it
     return params.map_modification.restrict_z_distance_for_helical_symmetry

  if not params.map_modification.restrict_z_turns_for_helical_symmetry: return

  print("Identifying maximum z-range for helical symmetry", file = out)
  print("Maximum of %7.1f turns up and down in Z allowed..." %(
     params.map_modification.restrict_z_turns_for_helical_symmetry), file = out)
  r, t = ncs_obj.ncs_groups()[0].helix_rt_forwards()
  cost = r[0]
  sint = r[1]
  import math
  theta = abs(180.*math.atan2(sint, cost)/3.14159)
  trans = abs(t)
  pitch = trans*360./max(0.1, theta)
  max_z = params.map_modification.restrict_z_turns_for_helical_symmetry*pitch
  print("Z-values restricted to +/- %7.1f A" %(max_z), file = out)
  print("\nRunning map-box once to get position of molecule, again to"+\
      " apply\n Z restriction\n", file = out)
  return max_z

def dist(x, y):
  """Return Euclidian distance between x and y."""
  dd = 0.
  for a, b in zip(x, y):
    dd+= (a-b)**2
  return dd**0.5

def get_ncs_closest_sites(
    closest_sites = None,
    sites_cart = None,
    used_ncs_id_list = None,
    box_ncs_object = None,
    box_crystal_symmetry = None,
    out = sys.stdout):

  """Try to find NCS ops mapping sites_cart close to closest_sites"""

  best_id = None
  best_rms = None
  best_sites = closest_sites.deep_copy()
  for ncs_id in range(box_ncs_object.max_operators()):
    if ncs_id in used_ncs_id_list: continue

    test_sites = closest_sites.deep_copy()
    ncs_sites_cart = get_ncs_sites_cart(sites_cart = sites_cart,
       ncs_obj = box_ncs_object, unit_cell = box_crystal_symmetry.unit_cell(),
       ncs_id = ncs_id,
       ncs_in_cell_only = False)
    test_sites.extend(ncs_sites_cart)
    rms = radius_of_gyration_of_vector(test_sites)
    if best_rms is None or rms < best_rms:
      best_rms = rms
      best_ncs_id = ncs_id
      best_sites = test_sites.deep_copy()
  used_ncs_id_list.append(best_ncs_id)
  return best_sites, used_ncs_id_list

def get_closest_sites(
    high_points = None,
    sites_cart = None,
    box_ncs_object = None,
    box_crystal_symmetry = None,
    out = sys.stdout):
  """Get closest site to high_points of
   any site ncs-related to each member of sites_cart."""
  if not box_ncs_object.is_point_group_symmetry() and not \
      box_ncs_object.is_helical_along_z():
    # extract point_group symmetry if present and box_ncs_object doesn't have it
    print("Trying to extract point-group symmetry from box_ncs_object "+\
       "with %d ops" %( box_ncs_object.max_operators()), file = out)

    ncs_object = box_ncs_object.deep_copy(extract_point_group_symmetry = True)
    if ncs_object:
      print("New number of operators satisfying point-group symmetry: %d" %(
        ncs_object.max_operators()), file = out)
      box_ncs_object = ncs_object

    else:
      print("No point-group symmetry found", file = out)


  ncs_copies = box_ncs_object.max_operators()
  closest_sites = high_points
  from scitbx.matrix import col
  for id in range(sites_cart.size()):
    local_sites_cart = sites_cart[id:id+1]
    local_sites_cart.extend(get_ncs_sites_cart(sites_cart = local_sites_cart,
       ncs_obj = box_ncs_object, unit_cell = box_crystal_symmetry.unit_cell(),
       ncs_in_cell_only = True))
    if local_sites_cart.size() <ncs_copies: continue  # some were out of range

    xx = col((0., 0., 0., ))
    for site in closest_sites:
      xx+= col(site)
    xx = xx/max(1, closest_sites.size())
    target = flex.vec3_double()
    target.append(xx)

    dd, id1, id2 = target.min_distance_between_any_pair_with_id(
       local_sites_cart)
    best_points = local_sites_cart[id2:id2+1]
    closest_sites.extend(best_points)
  return closest_sites[1:]

def get_range(sites = None, unit_cell = None, map_data = None,
   boundary_tolerance = None, out = sys.stdout):
  """Create a box that can just hold all sites. Move boundaries just inside
   unit cell in all directions (boundary_tolerance)"""
  x_values = flex.double()
  y_values = flex.double()
  z_values = flex.double()
  for site_cart in sites:
    (x, y, z) = tuple(site_cart)
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)
  x_min_max_mean = x_values.min_max_mean()
  x_min = x_min_max_mean.min
  x_max = x_min_max_mean.max
  y_min_max_mean = y_values.min_max_mean()
  y_min = y_min_max_mean.min
  y_max = y_min_max_mean.max
  z_min_max_mean = z_values.min_max_mean()
  z_min = z_min_max_mean.min
  z_max = z_min_max_mean.max
  print("\nRange for box:", file = out)
  print("            X        Y        Z", file = out)
  print(" LOW:  %7.1f    %7.1f     %7.1f " %(tuple([x_min, y_min, z_min])), file = out)
  print(" HIGH: %7.1f    %7.1f     %7.1f \n" %(tuple([x_max, y_max, z_max])), file = out)

  # move to 0, 1 if near ends
  if x_min<= boundary_tolerance: x_min = 0.
  if y_min<= boundary_tolerance: y_min = 0.
  if z_min<= boundary_tolerance: z_min = 0.
  a, b, c, al, bet, gam = unit_cell.parameters()
  if x_min>= a-boundary_tolerance: x_min = a
  if y_min>= b-boundary_tolerance: y_min = b
  if z_min>= c-boundary_tolerance: z_min = c
  print("\nAdjusted range for box:", file = out)
  print("            X        Y        Z", file = out)
  print(" LOW:  %7.1f    %7.1f     %7.1f " %(tuple([x_min, y_min, z_min])), file = out)
  print(" HIGH: %7.1f    %7.1f     %7.1f \n" %(tuple([x_max, y_max, z_max])), file = out)

  nx, ny, nz = map_data.all()
  # convert to grid units
  i_min = max(0, min(nx, int(0.5+nx*x_min/a)))
  j_min = max(0, min(ny, int(0.5+ny*y_min/b)))
  k_min = max(0, min(nz, int(0.5+nz*z_min/c)))
  i_max = max(0, min(nx, int(0.5+nx*x_max/a)))
  j_max = max(0, min(ny, int(0.5+ny*y_max/b)))
  k_max = max(0, min(nz, int(0.5+nz*z_max/c)))
  lower_bounds = [i_min, j_min, k_min]
  upper_bounds = [i_max, j_max, k_max]

  print("\nGrid bounds for box:", file = out)
  print("            X        Y        Z", file = out)
  print(" LOW:  %7d    %7d     %7d " %(tuple([i_min, j_min, k_min])), file = out)
  print(" HIGH: %7d    %7d     %7d \n" %(tuple([i_max, j_max, k_max])), file = out)

  return lower_bounds, upper_bounds

def get_bounds_for_au_box(params,
     box = None, out = sys.stdout):

  """Try to get bounds for a box that include one au"""

  if not box.ncs_object or box.ncs_object.max_operators()<2:
    return None, None, None

  box_ncs_object = box.ncs_object
  box_map_data = box.map_box
  box_crystal_symmetry = box.box_crystal_symmetry
  random_points = 10*params.reconstruction_symmetry.random_points

  sites_cart = get_points_in_map(box_map_data, n = random_points,
    minimum_fraction_of_max = params.segmentation.density_select_threshold,
    random_xyz = params.crystal_info.resolution*2.,
    crystal_symmetry = box_crystal_symmetry)
  assert sites_cart.size() >0

  # apply symmetry to sites_orth and see where the molecule is
  ncs_sites_cart = get_ncs_sites_cart(sites_cart = sites_cart,
       ncs_obj = box_ncs_object, unit_cell = box_crystal_symmetry.unit_cell(),
       ncs_in_cell_only = True)

  # generate this in a lower-memory way XXX
  low_res_map_data = get_low_res_map_data(sites_cart = ncs_sites_cart,
    d_min = params.crystal_info.resolution*7.,
    crystal_symmetry = box_crystal_symmetry,
    out = out)

  high_points = get_high_points_from_map(  # actually returns just one.
        map_data = low_res_map_data,
        unit_cell = box_crystal_symmetry.unit_cell(), out = out)
  from scitbx.matrix import col
  high_points = high_points[0:1]
  cutout_center = col(high_points[0])
  print("Center of box will be near (%7.1f, %7.1f, %7.1f)" %(
     tuple(cutout_center)))

  # now figure out box that contains at least one copy of each ncs-related
  #  point.

  # Find closest ncs-related points for each unique random point to this
  # center. Starting box is the box that contains all of these.

  closest_sites = get_closest_sites(
    high_points = high_points,
    sites_cart = sites_cart,
    box_ncs_object = box_ncs_object,
    box_crystal_symmetry = box_crystal_symmetry,
    out = out)
  if not closest_sites or closest_sites.size()<1:
    print("\nNo sites representing au of map found...skipping au box\n", file = out)
    return None, None, None

  print("\nTotal of %d sites representing 1 au found" %(
      closest_sites.size()), file = out)

  # write out closest_sites to match original position
  coordinate_offset = -1*matrix.col(box.shift_cart)
  write_atoms(file_name = 'one_au.pdb', # PDB OK just writing out atoms
      crystal_symmetry = box_crystal_symmetry, sites = closest_sites+coordinate_offset)

  unique_closest_sites = closest_sites.deep_copy()
  # Now if desired, find NCS-related groups of sites
  if params.segmentation.n_au_box >1:
    print("\nFinding up to %d related au" %(params.segmentation.n_au_box), file = out)
    print("Starting RMSD of sites: %7.1f A " %(
       radius_of_gyration_of_vector(closest_sites)), file = out)

    sites_orig = closest_sites.deep_copy()
    used_ncs_id_list = [box_ncs_object.ncs_groups()[0].identity_op_id()]
    for i in range(params.segmentation.n_au_box-1):
      closest_sites, used_ncs_id_list = get_ncs_closest_sites(
        used_ncs_id_list = used_ncs_id_list,
        closest_sites = closest_sites,
        sites_cart = sites_orig,
        box_ncs_object = box_ncs_object,
        box_crystal_symmetry = box_crystal_symmetry,
        out = out)
    print("\nNew total of %d sites representing %d au found" %(
      closest_sites.size(), params.segmentation.n_au_box), file = out)
    print("New rmsd: %7.1f A " %(
       radius_of_gyration_of_vector(closest_sites)), file = out)

  lower_bounds, upper_bounds = get_range(
     sites = closest_sites, map_data = box_map_data,
     boundary_tolerance = params.crystal_info.resolution,
     unit_cell = box_crystal_symmetry.unit_cell(),
     out = out)
  return lower_bounds, upper_bounds, unique_closest_sites+coordinate_offset

def get_low_res_map_data(sites_cart = None,
    crystal_symmetry = None,
    d_min = None,
    out = sys.stdout):
    """Get low resolution map data based on sites_cart"""

    from cctbx import xray
    xrs, scatterers = set_up_xrs(crystal_symmetry = crystal_symmetry)
    unit_cell = crystal_symmetry.unit_cell()
    sites_fract = unit_cell.fractionalize(sites_cart)
    for xyz_fract in sites_fract:
      scatterers.append( xray.scatterer(scattering_type = "H", label = "H",
        site = xyz_fract, u = 0, occupancy = 1.0))
    xrs = xray.structure(xrs, scatterers = scatterers)
    # generate f_array to d_min with xrs
    f_array =  xrs.structure_factors(d_min = d_min, anomalous_flag = False).f_calc()
    weight_f_array = f_array.structure_factors_from_scatterers(
      algorithm = 'direct',
      xray_structure = xrs).f_calc()

    return get_map_from_map_coeffs(map_coeffs = weight_f_array,
      crystal_symmetry = crystal_symmetry)

def get_bounds_for_helical_symmetry(params,
     box = None, crystal_symmetry = None):
  """Get bounds for searching for helical symmetry"""
  original_cell = box.map_data.all()
  new_cell = box.map_box.all()
  z_first = box.gridding_first[2]
  z_last = box.gridding_last[2]
  assert z_last>= z_first
  z_middle = (z_first+z_last)//2
  delta_z = crystal_symmetry.unit_cell().parameters()[5]/box.map_data.all()[2]
  n_z_max =  int(0.5+
   params.map_modification.restrict_z_distance_for_helical_symmetry/delta_z)
  new_z_first = max(z_first, z_middle-n_z_max)
  new_z_last = min(z_last, z_middle+n_z_max)
  lower_bounds = list(deepcopy(box.gridding_first))
  upper_bounds = list(deepcopy(box.gridding_last))
  lower_bounds[2] = new_z_first
  upper_bounds[2] = new_z_last

  return lower_bounds, upper_bounds

def check_memory(map_data, ratio_needed, maximum_fraction_to_use = 0.90,
    maximum_map_size = 1,
    out = sys.stdout):
  """Check memory used and estimate how much more is needed"""
  map_size = map_data.size()/(1024*1024*1024)
  if maximum_map_size and map_size>maximum_map_size:
     raise Sorry("Maximum map size for this tool is %s GB" %(maximum_map_size))
  needed_memory = ratio_needed*map_size
  from libtbx.utils import guess_total_memory # returns total memory

  bytes_total_memory = guess_total_memory()
  if bytes_total_memory:
    total_memory = bytes_total_memory/(1024*1024*1024)
  else:
    total_memory = None
  print("\nMap size is " +\
      "%.2f GB.  This will require about %.1f GB of memory" %(
      map_size, needed_memory) +"\nfor this stage of analysis\n", file = out)
  if total_memory:
    print("Total memory on this computer is about %.1f GB." %(
      total_memory), file = out)
    if (needed_memory>=  0.5* total_memory):
      print("\n ** WARNING:  It is possible that this computer may not"+\
       " have **\n *** sufficient memory to complete this job. ***\n", file = out)
    if (needed_memory >=  maximum_fraction_to_use*total_memory):
      raise Sorry("This computer does not have sufficient "+
        "memory (%.0f GB needed) \nto run this job" %(needed_memory))

def get_params(args, map_data = None, crystal_symmetry = None,
    half_map_data_list = None,
    sharpening_target_pdb_hierarchy = None,
    ncs_object = None,
    write_files = None,
    auto_sharpen = None,
    density_select = None,
    add_neighbors = None,
    save_box_map_ncs_au = None,
    sequence = None,
    wrapping = None,
    target_ncs_au_file = None,
    regions_to_keep = None,
    solvent_content = None,
    resolution = None,
    molecular_mass = None,
    symmetry = None,
    chain_type = None,
    keep_low_density = None,
    box_buffer = None,
    soft_mask_extract_unique = None,
    mask_expand_ratio = None,
    include_helical_symmetry = None,
    min_ncs_cc = None,
    symmetry_center = None,
    return_params_only = None,
    out = sys.stdout):

  """Get parameters from args for auto_sharpening"""
  params = get_params_from_args(args)

  # Set params specifically if coming in from call
  if sequence is not None:
    params.crystal_info.sequence = sequence
  if (wrapping is not None) and (params.crystal_info.use_sg_symmetry is None):
    params.crystal_info.use_sg_symmetry =  wrapping
  if target_ncs_au_file is not None:
    params.input_files.target_ncs_au_file = target_ncs_au_file
  if regions_to_keep is not None:
    params.map_modification.regions_to_keep = regions_to_keep
  if solvent_content is not None:
    params.crystal_info.solvent_content = solvent_content
  if resolution is not None:
    params.crystal_info.resolution = resolution
  if molecular_mass is not None:
    params.crystal_info.molecular_mass = molecular_mass
  if symmetry is not None:
    params.reconstruction_symmetry.symmetry = symmetry
  if min_ncs_cc is not None:
    params.reconstruction_symmetry.min_ncs_cc = min_ncs_cc
  if symmetry_center is not None:
    params.reconstruction_symmetry.symmetry_center = symmetry_center
  if include_helical_symmetry is not None:
    params.reconstruction_symmetry.include_helical_symmetry = \
      include_helical_symmetry
  if chain_type is not None:
    params.crystal_info.chain_type = chain_type

  if regions_to_keep is not None:
    params.map_modification.regions_to_keep = regions_to_keep
    params.segmentation.iterate_with_remainder = False
  elif keep_low_density:
    params.segmentation.iterate_with_remainder = True
  elif keep_low_density is False:
    params.segmentation.iterate_with_remainder = False
  else:
    pass # just so it is clear this was considered

  if box_buffer is not None:
    params.output_files.box_buffer = box_buffer
  if soft_mask_extract_unique is not None:
    params.map_modification.soft_mask = soft_mask_extract_unique

  if mask_expand_ratio is not None:
    params.segmentation.mask_expand_ratio = mask_expand_ratio
  if write_files is not None:
    params.control.write_files = write_files
  if auto_sharpen is not None:
    params.map_modification.auto_sharpen = auto_sharpen
  if density_select is not None:
    params.segmentation.density_select = density_select
  if add_neighbors is not None:
    params.segmentation.add_neighbors = add_neighbors
  if save_box_map_ncs_au is not None:
    params.control.save_box_map_ncs_au = save_box_map_ncs_au


  if return_params_only:
    return params

  print("\nSegment_and_split_map\n", file = out)
  print("Command used: %s\n" %(
   " ".join(['segment_and_split_map']+args)), file = out)
  master_params.format(python_object = params).show(out = out)

  # Set space-group defaults
  if params.crystal_info.use_sg_symmetry:
    if params.map_modification.restrict_map_size is None:
      params.map_modification.restrict_map_size = False
    if params.crystal_info.is_crystal is None:
      params.crystal_info.is_crystal = True
  else:
    if params.map_modification.restrict_map_size is None:
      params.map_modification.restrict_map_size = True
    if params.crystal_info.is_crystal is None:
      params.crystal_info.is_crystal = False

  # Turn off files if desired
  if params.control.write_files is False:
    params.output_files.magnification_map_file = None

    params.output_files.magnification_map_file = None
    params.output_files.magnification_ncs_file = None
    params.output_files.shifted_map_file = None
    params.output_files.shifted_sharpened_map_file = None
    params.output_files.sharpened_map_file = None
    params.output_files.shifted_pdb_file = None
    params.output_files.shifted_ncs_file = None
    params.output_files.shifted_used_ncs_file = None
    params.output_files.box_map_file = None
    params.output_files.box_mask_file = None
    params.output_files.write_output_maps = False
    params.output_files.remainder_map_file = None
    params.output_files.output_info_file = None
    params.output_files.restored_pdb = None
    params.output_files.output_weight_map_pickle_file = None


  from cctbx.maptbx.auto_sharpen import set_sharpen_params
  params = set_sharpen_params(params, out)


  if params.input_files.seq_file and not params.crystal_info.sequence and \
      not sequence:
    if not params.crystal_info.sequence:
      if sequence:
        params.crystal_info.sequence = sequence
      else:
        params.crystal_info.sequence = open(params.input_files.seq_file).read()
    print("Read sequence from %s" %(params.input_files.seq_file), file = out)


  if not params.crystal_info.resolution and (
     params.map_modification.b_iso is not None or \
      params.map_modification.auto_sharpen
      or params.map_modification.resolution_dependent_b or
      params.map_modification.b_sharpen):
    raise Sorry("Need resolution for segment_and_split_map with sharpening")

  if params.map_modification.auto_sharpen and (
      params.map_modification.b_iso is not None or
      params.map_modification.b_sharpen is not None or
      params.map_modification.resolution_dependent_b is not None):
    print("Turning off auto_sharpen as it is not compatible with "+\
        "b_iso, \nb_sharpen, or resolution_dependent_b", file = out)
    params.map_modification.auto_sharpen = False

  if params.control.write_files and \
       params.output_files.output_directory and  \
     (not os.path.isdir(params.output_files.output_directory)):
      os.mkdir(params.output_files.output_directory)
  if not params.output_files.output_directory:
    params.output_files.output_directory = ""

  # Test to see if we can use adjusted_sa as target and use box_map with it
  if (params.map_modification.residual_target == 'adjusted_sa' or
     params.map_modification.sharpening_target == 'adjusted_sa') and \
     (params.map_modification.box_in_auto_sharpen or
       params.map_modification.density_select_in_auto_sharpen) and (
      params.map_modification.auto_sharpen):

    print("Checking to make sure we can use adjusted_sa as target...",
        end = ' ', file = out)
    try:
      from phenix.autosol.map_to_model import iterated_solvent_fraction
      dummy = iterated_solvent_fraction # just to import it
    except Exception as e:
      raise Sorry("Please either set box_in_auto_sharpen = False and "+
       "\ndensity_select_in_auto_sharpen = False or \n"+\
      "set residual_target = kurtosis and sharpening_target = kurtosis")
    print("OK", file = out)

  if not half_map_data_list: half_map_data_list = []

  if params.input_files.info_file:
    map_data = None
    pdb_hierarchy = None
    from libtbx import easy_pickle
    print("Loading tracking data from %s" %(
      params.input_files.info_file), file = out)
    tracking_data = easy_pickle.load(params.input_files.info_file)
    return params, map_data, half_map_data_list, pdb_hierarchy, tracking_data, None
  else:
    tracking_data = info_object()
    tracking_data.set_params(params)

  # PDB file
  if params.input_files.pdb_file:
    print("\nInput PDB file to be used to identify region to work with: %s\n" %(
       params.input_files.pdb_file), file = out)
    from iotbx.pdb.utils import get_pdb_hierarchy
    pdb_hierarchy = get_pdb_hierarchy(file_name = params.input_files.pdb_file)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    tracking_data.set_input_pdb_info(file_name = params.input_files.pdb_file,
      n_residues = pdb_hierarchy.overall_counts().n_residues)
  else:
    pdb_hierarchy = None

  if map_data:
    pass # ok
  elif params.input_files.map_file:
    ccp4_map = iotbx.mrcfile.map_reader(
    file_name = params.input_files.map_file)
    if not crystal_symmetry:
      crystal_symmetry = ccp4_map.crystal_symmetry() # 2018-07-18
      tracking_data.set_full_crystal_symmetry(
         ccp4_map.unit_cell_crystal_symmetry())
      tracking_data.set_full_unit_cell_grid(ccp4_map.unit_cell_grid)
    map_data = ccp4_map.data.as_double()
  else:
    raise Sorry("Need ccp4 map")
  if not crystal_symmetry:
    raise Sorry("Need crystal_symmetry")

  if params.input_files.half_map_file:
    if len(params.input_files.half_map_file) !=  2:
      raise Sorry("Please supply none or two half_map_file values")

    half_map_data_list = []
    half_map_data_list.append(iotbx.mrcfile.map_reader(
       file_name = params.input_files.half_map_file[0]).data.as_double())
    half_map_data_list.append(iotbx.mrcfile.map_reader(
       file_name = params.input_files.half_map_file[1]).data.as_double())

  # Get the NCS object
  ncs_obj, dummy_tracking_data = get_ncs(params = params,
    ncs_object = ncs_object, out = out)

  if (not params.map_modification.auto_sharpen or
       params.map_modification.b_iso is not None) and (
      not params.crystal_info.molecular_mass and
      not params.crystal_info.solvent_content and
      not params.input_files.seq_file and not params.crystal_info.sequence and
      not sequence):
    params.crystal_info.solvent_content = get_iterated_solvent_fraction(
        crystal_symmetry = crystal_symmetry,
        verbose = params.control.verbose,
        resolve_size = params.control.resolve_size,
        mask_padding_fraction = \
           params.segmentation.mask_padding_fraction,
        fraction_of_max_mask_threshold = \
           params.segmentation.fraction_of_max_mask_threshold,
        cell_cutoff_for_solvent_from_mask = \
           params.segmentation.cell_cutoff_for_solvent_from_mask,
        mask_resolution = params.crystal_info.resolution,
        map = map_data,
        out = out)
    if params.crystal_info.solvent_content:
      print("Estimated solvent content: %.2f" %(
        params.crystal_info.solvent_content), file = out)
    else:
      raise Sorry("Unable to estimate solvent content...please supply "+
        "solvent_content \nor molecular_mass")

  if params.map_modification.auto_sharpen or \
        params.map_modification.b_iso is not None or \
        params.map_modification.b_sharpen is not None or \
        params.map_modification.resolution_dependent_b is not None:
      # Sharpen the map
      print("Auto-sharpening map before using it", file = out)
      local_params = deepcopy(params)
      if tracking_data.solvent_fraction: # XXX was previously always done but may not have been set
        local_params.crystal_info.solvent_content = tracking_data.solvent_fraction
      from cctbx.maptbx.auto_sharpen import run as auto_sharpen
      acc = map_data.accessor()
      map_data, new_map_coeffs, new_crystal_symmetry, new_si = auto_sharpen(
         args = [], params = local_params,
        map_data = map_data,
        wrapping = wrapping,
        crystal_symmetry = crystal_symmetry,
        write_output_files = False,
        pdb_hierarchy = sharpening_target_pdb_hierarchy,
        ncs_obj = ncs_obj,
        return_map_data_only = False,
        return_unshifted_map = True,
        half_map_data_list = half_map_data_list,
        n_residues = tracking_data.n_residues,
        ncs_copies = ncs_obj.max_operators(),
        out = out)
      tracking_data.b_sharpen = new_si.b_sharpen
      if not tracking_data.solvent_fraction:
        tracking_data.solvent_fraction = new_si.solvent_fraction

      if tracking_data.params.output_files.sharpened_map_file:
        sharpened_map_file = os.path.join(
            tracking_data.params.output_files.output_directory,
            tracking_data.params.output_files.sharpened_map_file)
        sharpened_map_data = map_data.deep_copy()
        if acc is not None:  # offset the map to match original if possible
          sharpened_map_data.reshape(acc)

        print("Gridding of sharpened map:", file = out)
        print("Origin: ", sharpened_map_data.origin(), file = out)
        print("All: ", sharpened_map_data.all(), file = out)
        print("\nWrote sharpened map in original location with "+\
             "origin at %s\nto %s" %(
           str(sharpened_map_data.origin()), sharpened_map_file), file = out)

        # NOTE: original unit cell and grid
        write_ccp4_map(tracking_data.full_crystal_symmetry,
          sharpened_map_file, sharpened_map_data,
          output_unit_cell_grid = tracking_data.full_unit_cell_grid, )
        params.input_files.map_file = sharpened_map_file # overwrite map_file name here

      # done with any sharpening
      params.map_modification.auto_sharpen = False# so we don't do it again later
      params.map_modification.b_iso = None
      params.map_modification.b_sharpen = None
      params.map_modification.resolution_dependent_b = None
      if params.control.sharpen_only:
        print("Stopping after sharpening", file = out)
        return

  # check on size right away
  if params.control.memory_check:
    # map_box and mask generation use about 50GB of memory for
    #    map with 1 billion elements
    check_memory(map_data = map_data, maximum_map_size = None,
      ratio_needed = 50, out = out)

  if params.map_modification.magnification and \
       params.map_modification.magnification!= 1.0:
    print("\nAdjusting magnification by %7.3f\n" %(
       params.map_modification.magnification), file = out)

    if ncs_obj:
      # Magnify ncs
      print("NCS before applying magnification...", file = out)
      ncs_obj.format_all_for_group_specification(out = out)
      ncs_obj = ncs_obj.adjust_magnification(
        magnification = params.map_modification.magnification)
      if params.output_files.magnification_ncs_file:
        file_name = os.path.join(params.output_files.output_directory,
          params.output_files.magnification_ncs_file)
        print("Writing NCS after magnification of %7.3f to %s" %(
          params.map_modification.magnification, file_name), file = out)
        ncs_obj.format_all_for_group_specification(out = out)
        ncs_obj.format_all_for_group_specification(file_name = file_name)
        params.input_files.ncs_file = file_name
      else:
        raise Sorry("Need magnification_ncs_file defined if magnification is"+
          " applied \nto input NCS file")

    # Magnify map
    shrunk_uc = []
    for i in range(3):
      shrunk_uc.append(
       crystal_symmetry.unit_cell().parameters()[i] *
          params.map_modification.magnification )
    uc_params = crystal_symmetry.unit_cell().parameters()
    from cctbx import uctbx
    new_unit_cell = uctbx.unit_cell(
      parameters = (shrunk_uc[0], shrunk_uc[1], shrunk_uc[2],
          uc_params[3], uc_params[4], uc_params[5]))
    print("Original unit cell: (%7.4f, %7.4f, %7.4f, %7.4f, %7.4f, %7.4f)" %(
      crystal_symmetry.unit_cell().parameters()), file = out)
    crystal_symmetry = crystal.symmetry(
      unit_cell = new_unit_cell,
      space_group = crystal_symmetry.space_group())
    print("New unit cell:      (%7.4f, %7.4f, %7.4f, %7.4f, %7.4f, %7.4f)" %(
      crystal_symmetry.unit_cell().parameters()), file = out)

    # magnify original unit cell too..
    cell = list(tracking_data.full_crystal_symmetry.unit_cell().parameters())
    for i in range(3):
      cell[i] = cell[i]*params.map_modification.magnification
    tracking_data.set_full_crystal_symmetry(
        crystal.symmetry(tuple(cell), ccp4_map.space_group_number))
    print("New original (full unit cell): "+\
        "  (%7.4f, %7.4f, %7.4f, %7.4f, %7.4f, %7.4f)" %(
      tracking_data.full_crystal_symmetry.unit_cell.parameters()), file = out)

    if params.output_files.magnification_map_file:
      file_name = os.path.join(params.output_files.output_directory,
        params.output_files.magnification_map_file)
      # write out magnified map (our working map) (before shifting it)
      print("\nWriting magnification map (input map with "+\
        "magnification of %7.3f \n" %(params.map_modification.magnification) +\
        "applied) to %s \n" %(file_name), file = out)
      #write_ccp4_map(crystal_symmetry, file_name, map_data)
      #  NOTE: original unit cell and grid
      write_ccp4_map(tracking_data.full_crystal_symmetry,
        file_name, map_data,
        output_unit_cell_grid = tracking_data.original_unit_cell_grid, )
      params.input_files.map_file = file_name
    else:
      raise Sorry("Need a file name to write out magnification_map_file")
    params.map_modification.magnification = None  # no longer need it.

  tracking_data.set_input_map_info(file_name = params.input_files.map_file,
    crystal_symmetry = crystal_symmetry,
    origin = map_data.origin(),
    all = map_data.all())
  tracking_data.set_crystal_symmetry(crystal_symmetry = crystal_symmetry)
  tracking_data.set_original_crystal_symmetry(crystal_symmetry = crystal_symmetry)
  tracking_data.set_accessor(acc = map_data.accessor())

  # Save center of map
  map_symmetry_center = get_center_of_map(map_data, crystal_symmetry)

  # Check for helical ncs...if present we may try to cut map at +/- 1 turn
  params.map_modification.restrict_z_distance_for_helical_symmetry = \
     get_max_z_range_for_helical_symmetry(params, out = out)

  # either use map_box with density_select = True or just shift the map
  if  params.segmentation.density_select:
    print("\nTrimming map to density...", file = out)
    args =  ["output_format = ccp4"]
    if params.segmentation.density_select_threshold is not None:
      print("Threshold for density selection will be: %6.2f \n"%(
       params.segmentation.density_select_threshold), file = out)
      args.append("density_select_threshold = %s" %(
         params.segmentation.density_select_threshold))
    if params.segmentation.get_half_height_width is not None:
      args.append("get_half_height_width = %s" %(
        params.segmentation.get_half_height_width))
    if params.input_files.ncs_file:
      args.append("symmetry_file = %s" %(params.input_files.ncs_file))
    if params.input_files.pdb_file:
      args.append("pdb_file = %s" %(params.input_files.pdb_file))
    args.append("ccp4_map_file = %s" %(params.input_files.map_file))
    file_name_prefix = os.path.join(params.output_files.output_directory,
       "density_select")
    args.append("output_file_name_prefix = %s" %(file_name_prefix))
    from mmtbx.command_line.map_box import run as run_map_box
    args.append("keep_input_unit_cell_and_grid = False") # for new defaults

    assert params.crystal_info.use_sg_symmetry is not None
    wrapping = params.crystal_info.use_sg_symmetry
    if params.segmentation.lower_bounds and params.segmentation.upper_bounds:
      bounds_supplied = True
      print("\nRunning map_box with supplied bounds", file = out)
      box = run_map_box(args,
          map_data = map_data,
          ncs_object = ncs_obj,
          crystal_symmetry = crystal_symmetry,
          lower_bounds = params.segmentation.lower_bounds,
          upper_bounds = params.segmentation.upper_bounds,
          write_output_files = params.output_files.write_output_maps,
          wrapping = wrapping,
          log = out)
    else:
      bounds_supplied = False
      box = run_map_box(["density_select = True"]+args,
       map_data = map_data,
       crystal_symmetry = crystal_symmetry,
       ncs_object = ncs_obj,
       write_output_files = params.output_files.write_output_maps,
       wrapping = wrapping,
       log = out)

    # Run again to select au box
    shifted_unique_closest_sites = None
    selected_au_box = None
    if params.segmentation.select_au_box is None and  box.ncs_object and \
      box.ncs_object.max_operators() >=  params.segmentation.n_ops_to_use_au_box:
      params.segmentation.select_au_box = True
      print("Setting select_au_box to True as there are %d operators" %(
        box.ncs_object.max_operators()), file = out)
    if params.segmentation.select_au_box and not bounds_supplied:
      lower_bounds, upper_bounds, unique_closest_sites = get_bounds_for_au_box(
         params, box = box, out = out) #unique_closest_sites relative to original map
      if lower_bounds and upper_bounds:
        bounds_supplied = True
        selected_au_box = True
        score, ncs_cc = score_ncs_in_map(map_data = box.map_box,
          allow_score_with_pg = False,
          sites_orth = unique_closest_sites+box.shift_cart,
          ncs_object = box.ncs_object, ncs_in_cell_only = True,
          crystal_symmetry = box.box_crystal_symmetry, out = null_out())
        print("NCS CC before rerunning box: %7.2f   SCORE: %7.1f OPS: %d " %(
         ncs_cc, score, box.ncs_object.max_operators()), file = out)

        print("\nRunning map-box again with boxed range ...", file = out)
        del box
        box = run_map_box(args, lower_bounds = lower_bounds,
          map_data = map_data,
          crystal_symmetry = crystal_symmetry,
          ncs_object = ncs_obj,
          wrapping = wrapping,
          upper_bounds = upper_bounds, log = out)
        box.map_box = box.map_box.as_double()  # Do we need double?
        shifted_unique_closest_sites = unique_closest_sites+box.shift_cart

    # Or run again for helical symmetry
    elif params.map_modification.restrict_z_distance_for_helical_symmetry and \
       not bounds_supplied:
      bounds_supplied = True
      lower_bounds, upper_bounds = get_bounds_for_helical_symmetry(params,
         crystal_symmetry = crystal_symmetry, box = box)
      print("\nRunning map-box again with restricted Z range ...", file = out)
      box = run_map_box(args,
        map_data = map_data,
        crystal_symmetry = crystal_symmetry,
        ncs_object = ncs_obj,
        lower_bounds = lower_bounds, upper_bounds = upper_bounds,
        wrapping = wrapping,
        write_output_files = params.output_files.write_output_maps,
        log = out)

    #-----------------------------

    if bounds_supplied and box.ncs_object:
      print("Selecting remaining NCS operators", file = out)
      box.ncs_object = select_remaining_ncs_ops(
        map_data = box.map_box,
        crystal_symmetry = box.box_crystal_symmetry,
        closest_sites = shifted_unique_closest_sites,
        random_points = params.reconstruction_symmetry.random_points,
        ncs_object = box.ncs_object,
        out = out)

      score, ncs_cc = score_ncs_in_map(map_data = box.map_box,
          allow_score_with_pg = False,
          ncs_object = box.ncs_object, ncs_in_cell_only = True,
          sites_orth = shifted_unique_closest_sites,
          crystal_symmetry = box.box_crystal_symmetry, out = null_out())

      if score is not None:
        print("NCS CC after selections: %7.2f   SCORE: %7.1f  OPS: %d" %(
         ncs_cc, score, box.ncs_object.max_operators()), file = out)
    #-----------------------------

    origin_shift = box.shift_cart
    # Note: moving cell with (0, 0, 0) in middle to (0, 0, 0) at corner means
    #   total_shift_cart and origin_shift both positive
    map_data = box.map_box
    map_data = scale_map(map_data, out = out)
    crystal_symmetry = box.box_crystal_symmetry
    print("New unit cell: %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f " %(
      crystal_symmetry.unit_cell().parameters()), file = out)
    tracking_data.set_crystal_symmetry(crystal_symmetry = crystal_symmetry)
    print("Moving origin to (0, 0, 0)", file = out)
    print("Adding (%8.2f, %8.2f, %8.2f) to all coordinates\n"%(
        tuple(origin_shift)), file = out)
    # NOTE: size and cell params are now different!
    tracking_data.set_box_map_bounds_first_last(
       box.gridding_first, box.gridding_last)

    new_half_map_data_list = []
    ii = 0
    for hm in half_map_data_list:
      ii+= 1
      hm = hm.shift_origin() # shift if necessary
      hm = box.cut_and_copy_map(map_data = hm).as_double()
      hm.reshape(flex.grid(hm.all()))
      new_half_map_data_list.append(hm)
      cutout_half_map_file = os.path.join(params.output_files.output_directory,
       "cutout_half_map_%s.ccp4" %(ii))
      print("Writing cutout half_map data to %s" %(cutout_half_map_file), file = out)
      write_ccp4_map(crystal_symmetry, cutout_half_map_file, new_half_map_data_list[-1])

    half_map_data_list = new_half_map_data_list
    if params.map_modification.soft_mask:
      mask_data, map_data, half_map_data_list, \
      soft_mask_solvent_fraction, smoothed_mask_data, \
      original_box_map_data = \
       get_and_apply_soft_mask_to_maps(
        resolution = params.crystal_info.resolution,
        wang_radius = params.crystal_info.wang_radius,
        buffer_radius = params.crystal_info.buffer_radius,
        map_data = map_data, crystal_symmetry = crystal_symmetry,
        half_map_data_list = half_map_data_list,
        out = out)
      print("\nSolvent fraction from soft mask procedure: %7.2f (not used)\n" %(
        soft_mask_solvent_fraction), file = out)
    shifted_ncs_object = box.ncs_object
    if not shifted_ncs_object or shifted_ncs_object.max_operators()<2:
      from mmtbx.ncs.ncs import ncs
      shifted_ncs_object = ncs()
      shifted_ncs_object.set_unit_ncs()
  else:  # shift if necessary...
    shift_needed = not \
        (map_data.focus_size_1d() > 0 and map_data.nd()  ==  3 and
         map_data.is_0_based())

    a, b, c = crystal_symmetry.unit_cell().parameters()[:3]
    N_ = map_data.all()
    O_  = map_data.origin()
    sx, sy, sz =  O_[0]/N_[0], O_[1]/N_[1], O_[2]/N_[2]
    # Note: If (0, 0, 0) is in the middle of the box, origin at sx, sy, sz
    #  is negative, shift of coordinates will be positive
    sx_cart, sy_cart, sz_cart = crystal_symmetry.unit_cell().orthogonalize(
       [sx, sy, sz])
    print("Origin for input map is at (%8.2f, %8.2f, %8.2f)" % (
      sx_cart, sy_cart, sz_cart), file = out)
    print("Cell dimensions of this map are: (%8.2f, %8.2f, %8.2f)" % (a, b, c), file = out)
    if shift_needed:
      if(not crystal_symmetry.space_group().type().number() in [0, 1]):
          raise RuntimeError("Not implemented")
      origin_shift = [-sx_cart, -sy_cart, -sz_cart] # positive if (0, 0, 0) in middle
      print("Adding (%8.2f, %8.2f, %8.2f) to all coordinates"%(
        tuple(origin_shift))+" to put origin at (0, 0, 0)\n", file = out)

      map_data = map_data.shift_origin()
      new_half_map_data_list = []
      for hm in half_map_data_list:
        new_half_map_data_list.append(hm.shift_origin())
      half_map_data_list = new_half_map_data_list
    else:
      origin_shift = (0., 0., 0.)

    # Get NCS object if any
    if params.input_files.ncs_file and not ncs_obj:
      ncs_obj, dummy_obj = get_ncs(file_name = params.input_files.ncs_file)
    if ncs_obj:
      shifted_ncs_object = ncs_obj.coordinate_offset(
        coordinate_offset = matrix.col(origin_shift)) # shift to match shifted map
    else:
      from mmtbx.ncs.ncs import ncs
      shifted_ncs_object = ncs()
      shifted_ncs_object.set_unit_ncs()


  update_tracking_data_with_sharpening(
             map_data = map_data,
             tracking_data = tracking_data, out = out)

  # Set origin shift now
  tracking_data.set_origin_shift(origin_shift)

  map_symmetry_center = matrix.col(map_symmetry_center)+matrix.col(origin_shift) # New ctr

  if shifted_ncs_object and params.control.check_ncs:
    ncs_obj_to_check = shifted_ncs_object
  else:
    ncs_obj_to_check = None

  found_ncs = False
  if params.reconstruction_symmetry.symmetry or ncs_obj_to_check or \
     params.reconstruction_symmetry.optimize_center:
    looking_for_ncs = True
    new_ncs_obj, ncs_cc, ncs_score = run_get_ncs_from_map(params = params,
      map_data = map_data,
      map_symmetry_center = map_symmetry_center,
      crystal_symmetry = crystal_symmetry,
      ncs_obj = ncs_obj_to_check,
      out = out,
      )

    if new_ncs_obj:
      found_ncs = True
      shifted_ncs_object = new_ncs_obj.deep_copy()
      # offset this back to where it would have been before the origin offset..
      new_ncs_obj = new_ncs_obj.coordinate_offset(
       coordinate_offset = -1*matrix.col(origin_shift))
      # XXX save it in tracking_data

      if params.output_files.output_directory:
        if not os.path.isdir(params.output_files.output_directory):
          os.mkdir(params.output_files.output_directory)
        file_name = os.path.join(params.output_files.output_directory,
          'ncs_from_map.ncs_spec')
        f = open(file_name, 'w')
        new_ncs_obj.format_all_for_group_specification(out = f)
        f.close()
        print("Wrote NCS operators (for original map) to %s" %(file_name), file = out)
        if not params.control.check_ncs:
          params.input_files.ncs_file = file_name # set it

  else:
    looking_for_ncs = False

  if params.control.check_ncs:
    print("Done checking NCS", file = out)
    return params, map_data, half_map_data_list, pdb_hierarchy, tracking_data, None

  if looking_for_ncs and (not found_ncs) and \
         params.reconstruction_symmetry.symmetry.upper() not in ['ANY', 'ALL']:
      raise Sorry(
        "Unable to identify %s symmetry automatically in this map." %(
        params.reconstruction_symmetry.symmetry)+
        "\nPlease supply a symmetry file with symmetry matrices.")

  if params.segmentation.expand_size is None:
    params.segmentation.expand_size = estimate_expand_size(
       crystal_symmetry = crystal_symmetry,
       map_data = map_data,
       expand_target = params.segmentation.expand_target,
       out = out)

  if params.output_files.output_info_file and params.control.shift_only:
    write_info_file(params = params, tracking_data = tracking_data, out = out)

  return params, map_data, half_map_data_list, pdb_hierarchy, \
     tracking_data, shifted_ncs_object

def write_info_file(params = None, tracking_data = None, out = sys.stdout):
    """Write out summary info file based on tracking_data"""
    from libtbx import easy_pickle
    tracking_data.show_summary(out = out)
    print("\nWriting summary information to: %s" %(
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.output_info_file)), file = out)
    print("\nTo restore original position of a PDB file built into these maps, use:", file = out)
    print("phenix.segment_and_split_map info_file = %s" %(
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.output_info_file))+" pdb_to_restore = mypdb.pdb\n", file = out)
    easy_pickle.dump(os.path.join(tracking_data.params.output_files.output_directory, params.output_files.output_info_file),
       tracking_data)


def get_and_apply_soft_mask_to_maps(
    resolution = None,  #params.crystal_info.resolution
    wang_radius = None, #params.crystal_info.wang_radius
    buffer_radius = None, #params.crystal_info.buffer_radius
    force_buffer_radius = None, # apply buffer radius always
    map_data = None, crystal_symmetry = None,
    solvent_content = None,
    solvent_content_iterations = None,
    return_masked_fraction = True,
    rad_smooth = None,
    half_map_data_list = None,
    out = sys.stdout):
  """Create a soft mask (Gaussian fall-off at edges) based on region
  inside map containing macromolecule as determined by values of map_data.
  Then apply the soft mask around the macromolecular region.
  Returns: mask_data, map_data, half_map_data_list,
    solvent_fraction, smoothed_mask_data, original_map_data.
  """

  smoothed_mask_data = None
  if not resolution:
    from cctbx.maptbx import d_min_from_map
    resolution = d_min_from_map(
      map_data, crystal_symmetry.unit_cell(), resolution_factor = 1./4.)

  if not rad_smooth:
    rad_smooth = resolution

  if rad_smooth:
    print("\nApplying soft mask with smoothing radius of %.2f A\n" %(
      rad_smooth), file = out)
  if wang_radius:
    wang_radius = wang_radius
  else:
    wang_radius = 1.5*resolution

  if buffer_radius is not None:
    buffer_radius = buffer_radius
  else:
    buffer_radius = 2.*resolution
  original_map_data = map_data.deep_copy()
  # Check to make sure this is possible
  cell_dims = crystal_symmetry.unit_cell().parameters()[:3]
  min_cell_dim = min(cell_dims)
  if wang_radius > 0.25 * min_cell_dim or buffer_radius > 0.25 * min_cell_dim:
    new_wang_radius = min(wang_radius, 0.25 * min_cell_dim)
    new_buffer_radius = min(buffer_radius, 0.25 * min_cell_dim)
    print ("Cell is too small to get solvent fraction ...resetting "+
       "values of wang_radius \n"+
      "(was %.3f A now %.3f A) and buffer_radius (was %.3f A now %.3f A)" %(
       wang_radius, new_wang_radius, buffer_radius, new_buffer_radius), file = out)
    wang_radius = new_wang_radius
    buffer_radius = new_buffer_radius

  mask_data, solvent_fraction = get_mask_around_molecule(map_data = map_data,
    crystal_symmetry = crystal_symmetry,
    wang_radius = wang_radius,
    solvent_content = solvent_content,
    solvent_content_iterations = solvent_content_iterations,
    buffer_radius = buffer_radius,
    force_buffer_radius = force_buffer_radius,
    return_masked_fraction = return_masked_fraction,
    out = out)
  if mask_data:
    map_data, smoothed_mask_data = apply_soft_mask(map_data = map_data,
      mask_data = mask_data.as_double(),
      rad_smooth = rad_smooth,
      crystal_symmetry = crystal_symmetry,
      out = out)

    new_half_map_data_list = []
    if not half_map_data_list: half_map_data_list = []
    for half_map in half_map_data_list:
      assert half_map.size() == mask_data.size()
      half_map, smoothed_mask_data = apply_soft_mask(map_data = half_map,
        mask_data = mask_data.as_double(),
        rad_smooth = rad_smooth,
        crystal_symmetry = crystal_symmetry,
        out = out)
      new_half_map_data_list.append(half_map)
    half_map_data_list = new_half_map_data_list
  else:
    print("Unable to get mask...skipping", file = out)
  return mask_data, map_data, half_map_data_list, \
    solvent_fraction, smoothed_mask_data, original_map_data

def get_ncs(params = None, tracking_data = None, file_name = None,
     ncs_object = None, out = sys.stdout):
  """Return ncs object based on parameters, tracking_data, file_name,
   or supplied ncs_object"""
  if not file_name:
    file_name = params.input_files.ncs_file
  if (not ncs_object or ncs_object.max_operators()<2) and file_name: print("Reading ncs from %s" %(file_name), file = out)
  is_helical_symmetry = None
  if (not ncs_object or ncs_object.max_operators()<2) and not file_name: # No ncs supplied...use just 1 ncs copy..
    from mmtbx.ncs.ncs import ncs
    ncs_object = ncs()
    ncs_object.set_unit_ncs()
    #ncs_object.display_all(log = out)
  elif (not ncs_object or ncs_object.max_operators()<2) and \
      not os.path.isfile(file_name):
    raise Sorry("The ncs file %s is missing" %(file_name))
  else: # get the ncs
    if not ncs_object:
      from mmtbx.ncs.ncs import ncs
      ncs_object = ncs()
      try: # see if we can read biomtr records
        from iotbx.pdb.utils import get_pdb_input
        pdb_inp = get_pdb_input(file_name = file_name)
        ncs_object.ncs_from_pdb_input_BIOMT(pdb_inp = pdb_inp, log = out)
      except Exception as e: # try as regular ncs object
        ncs_object.read_ncs(file_name = file_name, log = out)
      #ncs_object.display_all(log = out)
    ncs_object.select_first_ncs_group()
    if ncs_object.max_operators()<1:
      from mmtbx.ncs.ncs import ncs
      ncs_object = ncs()
      ncs_object.set_unit_ncs()
    print("\nTotal of %d NCS operators read\n" %(
      ncs_object.max_operators()), file = out)
    if not tracking_data or not params:
      return ncs_object, None
    if ncs_object.max_operators()<2:
       print("No NCS present", file = out)
    elif ncs_object.is_helical_along_z(
       abs_tol_t = tracking_data.params.reconstruction_symmetry.abs_tol_t,
       rel_tol_t = tracking_data.params.reconstruction_symmetry.rel_tol_t,
       tol_r = tracking_data.params.reconstruction_symmetry.tol_r):
      print("This NCS is helical symmetry", file = out)
      is_helical_symmetry = True
    elif ncs_object.is_point_group_symmetry(
       abs_tol_t = tracking_data.params.reconstruction_symmetry.abs_tol_t,
       rel_tol_t = tracking_data.params.reconstruction_symmetry.rel_tol_t,
       tol_r = tracking_data.params.reconstruction_symmetry.tol_r):
      print("This NCS is point-group symmetry", file = out)
    elif params.crystal_info.is_crystal:
      print("This NCS is crystal symmetry", file = out)
    elif not (
      params.reconstruction_symmetry.require_helical_or_point_group_symmetry):
      print("WARNING: NCS is not crystal symmetry nor point-group "+\
         "symmetry nor helical symmetry", file = out)
    else:
      raise Sorry("Need point-group or helical symmetry.")
  if not ncs_object or ncs_object.max_operators()<1:
    raise Sorry("Need ncs information from an ncs_info file")
  if tracking_data:
    tracking_data.set_input_ncs_info(file_name = file_name,  # XXX may be updated ops
      number_of_operators = ncs_object.max_operators())

  if tracking_data and is_helical_symmetry: # update shifted_ncs_info
    if tracking_data.shifted_ncs_info: # XXX may not be needed
       shifted = True
    else:
       shifted = False
    print("Updating NCS info (shifted = %s)" %(shifted), file = out)
    tracking_data.update_ncs_info(is_helical_symmetry = True, shifted = shifted)

    if tracking_data.input_map_info and tracking_data.input_map_info.all:
      z_range = tracking_data.crystal_symmetry.unit_cell(). \
         parameters()[2]
      print("Extending NCS operators to entire cell (z_range = %.1f)" %(
         z_range), file = out)
      max_operators =  \
          tracking_data.params.reconstruction_symmetry.max_helical_operators
      if max_operators:
        print("Maximum new number of NCS operators will be %s" %(
         max_operators), file = out)
      ncs_object.extend_helix_operators(z_range = z_range,
        max_operators = max_operators)
      #ncs_object.display_all()
      print("New number of NCS operators is: %s " %(
        ncs_object.max_operators()), file = out)
      tracking_data.update_ncs_info(
        number_of_operators = ncs_object.max_operators(), is_helical_symmetry = True,
        shifted = shifted)
  return ncs_object, tracking_data

def score_threshold(b_vs_region = None, threshold = None,
     sorted_by_volume = None, n_residues = None,
     ncs_copies = None,
     fraction_occupied = None,
     solvent_fraction = None,
     map_data = None,
     residues_per_region = 50,
     min_volume = None,
     min_ratio = None,
     max_ratio_to_target = None,
     min_ratio_to_target = None,
     weight_score_grid_points = 1.,
     weight_score_ratio = 1.0,
     weight_near_one = 0.1,
     min_ratio_of_ncs_copy_to_first = None,
     target_in_all_regions = None,
     crystal_symmetry = None,
     chain_type = None,
     out = sys.stdout):
   """Score a value of threshold for identifying regions containing
    macromolecule. The threshold is used in the connectivity method to
    segment a map into connected regions.
    We want about 1 region per 50-100 residues for the biggest region.
    One possibility is to try to maximize the median size of the N top
    regions, where N = number of expected regions =  n_residues/residues_per_region

    Also note we have an idea how big a region should be (how many
    grid points) if we make an assumption about the fractional volume that
    should be inside a region compared to the total volume of protein/nucleic
    acid in the region...this gives us target_in_top_regions points.
    So using this, make the median size as close to target_in_top_regions as
    we can.
   """

   # If we have solvent fraction but not ncs_copies or n_residues, guess the
   #  number of residues and ncs copies from the volume
   if ncs_copies is not None and n_residues is not None:
     expected_regions = max(ncs_copies,
      max(1, int(0.5+n_residues/residues_per_region)))
   else:
      if chain_type in [None, 'None']: chain_type = "PROTEIN"
      assert crystal_symmetry is not None
      assert solvent_fraction is not None
      volume_per_residue, nres, chain_type = get_volume_of_seq(
          "A", chain_type = chain_type, out = out)

      expected_regions = max(1, int(0.5+(1-solvent_fraction)*\
        crystal_symmetry.unit_cell().volume()/volume_per_residue ))
      # NOTE: This is expected residues. expected_regions should be this
      #  divided by residues_per_region
      expected_regions = max(1, int(0.5+expected_regions/residues_per_region))
      ncs_copies = 1

   target_in_top_regions = target_in_all_regions/expected_regions

   nn = len(sorted_by_volume)-1 # first one is total
   ok = True

   too_low = None  # marker for way too low
   too_high = None

   if nn < ncs_copies:
     ok = False #return  # not enough

   v1, i1 = sorted_by_volume[1]
   if v1 < min_volume:
     ok = False #return

   if v1 > max_ratio_to_target*target_in_top_regions:
     ok = False #return
     too_low = True

   if v1 < min_volume or v1 < 0.1*min_ratio_to_target*target_in_top_regions:
     # way too high
     too_high = True


   # there should be about ncs_copies copies of each size region if ncs_copies>1
   if ncs_copies>1:
     v2, i2 = sorted_by_volume[max(1, min(ncs_copies, nn))]
     score_ratio = v2/v1  # want it to be about 1
     if score_ratio < min_ratio_of_ncs_copy_to_first:
       ok = False #return  # not allowed
   else:
     score_ratio = 1.0 # for ncs_copies = 1

   nn2 = min(nn, max(1, (expected_regions+1)//2))
   median_number, iavg = sorted_by_volume[nn2]

   # number in each region should be about target_in_top_regions

   if median_number > target_in_top_regions:
     score_grid_points = target_in_top_regions/max(1., median_number)
   else:
     score_grid_points = median_number/target_in_top_regions

   if v1> target_in_top_regions:
     score_grid_points_b = target_in_top_regions/max(1., v1)
   else:
     score_grid_points_b = v1/target_in_top_regions

   score_grid_points = 0.5*(score_grid_points+score_grid_points_b)

   score_grid_points = score_grid_points**2  # maybe even **3

   if threshold>1.:
     score_near_one = 1./threshold
   else:
     score_near_one = threshold

   # Normalize weight_score_ratio by target_in_top_regions:
   sc = min(1., 0.5*median_number/max(1, target_in_top_regions))
   overall_score = (
     (sc*weight_score_ratio*score_ratio+
     weight_score_grid_points*score_grid_points+
     weight_near_one*score_near_one
       ) /
     (weight_score_ratio+weight_score_grid_points+weight_near_one))

   half_expected_regions = max(1, (1+expected_regions)//2)
   ratio = sorted_by_volume[min(len(sorted_by_volume)-1, half_expected_regions)][0]/v1

   if ok and v1 >=  target_in_top_regions/2 and \
        len(sorted_by_volume)>half_expected_regions:
     last_volume = sorted_by_volume[half_expected_regions][0]
     if ratio >= min_ratio and \
         last_volume>= min_volume:
       has_sufficient_regions = True
     else:
       has_sufficient_regions = False
   else:
       has_sufficient_regions = False


   print("%7.2f  %5.2f   %5d     %4d    %5d     %5d     %6.3f   %5s    %5.3f  %s  %s" %(
       b_vs_region.b_iso, threshold, target_in_top_regions, expected_regions,
       v1, median_number, ratio, has_sufficient_regions, overall_score, ok, nn), file = out)

   if not b_vs_region.b_iso in b_vs_region.b_vs_region_dict.keys():
     b_vs_region.b_vs_region_dict[b_vs_region.b_iso] = {}
     b_vs_region.sa_sum_v_vs_region_dict[b_vs_region.b_iso] = {}
     b_vs_region.sa_nn_vs_region_dict[b_vs_region.b_iso] = {}
     b_vs_region.sa_ratio_b_vs_region_dict[b_vs_region.b_iso] = {}
   b_vs_region.b_vs_region_dict[b_vs_region.b_iso][threshold] = nn
   b_vs_region.sa_nn_vs_region_dict[b_vs_region.b_iso][threshold] = None
   b_vs_region.sa_ratio_b_vs_region_dict[b_vs_region.b_iso][threshold] = None

   return overall_score, has_sufficient_regions, \
      too_low, too_high, expected_regions, ok


def choose_threshold(b_vs_region = None, map_data = None,
     fraction_occupied = None,
     solvent_fraction = None,
     n_residues = None,
     ncs_copies = None,
     scale = 0.95,
     calculate_sa = None, # calculate surface area of top sa_percent of target
     sa_percent = None, # calculate surface area of top sa_fraction of target
     density_threshold = None,
     starting_density_threshold = None,
     wrapping = None,
     residues_per_region = None,
     min_volume = None,
     min_ratio = None,
     max_ratio_to_target = None,
     min_ratio_to_target = None,
     min_ratio_of_ncs_copy_to_first = None,
     verbose = None,
     crystal_symmetry = None,
     chain_type = None,
     out = sys.stdout):

  """Choose a threshold for calculating connectivity"""
  best_threshold = None
  best_threshold_has_sufficient_regions = None
  best_score = None
  best_ok = None

  if not ncs_copies: ncs_copies = 1

  print("\nChecking possible cutoffs for region identification", file = out)
  print("Scale: %7.3f" %(scale), file = out)
  used_ranges = []

  # Assume any threshold that is lower than a threshold that gave a non-zero value
  #  and is zero is an upper bound on the best value.  Same the other way around
  upper_bound = 1000
  lower_bound = 0.0001
  best_nn = None

  if density_threshold is not None: # use it
     print("\nUsing input threshold of %5.2f " %(
      density_threshold), file = out)
     n_range_low_high_list = [[0, 0]] # use as is
  else:
    n_range_low_high_list = [[-16, 4], [-32, 16], [-64, 80]]
    if starting_density_threshold is not None:
      starting_density_threshold = starting_density_threshold
      print("Starting density threshold is: %7.3f" %(
         starting_density_threshold), file = out)
    else:
      starting_density_threshold = 1.0
  if verbose:
    local_out = out
  else:
    from libtbx.utils import null_out
    local_out = null_out()

  target_in_all_regions = map_data.size()*fraction_occupied*(1-solvent_fraction)
  print("\nTarget number of points in all regions: %.0f" %(
    target_in_all_regions), file = local_out)


  local_threshold = find_threshold_in_map(target_points = int(
       target_in_all_regions), map_data = map_data)
  print("Cutoff will be threshold of %7.2f marking %7.1f%% of cell" %(
            local_threshold, 100.*(1.-solvent_fraction)), file = out)

  print("B-iso  Threshold  Target    N     Biggest   Median     Ratio   Enough  Score   OK  Regions", file = local_out)
  unique_expected_regions = None
  for n_range_low, n_range_high in n_range_low_high_list:
    last_score = None
    for nn in range(n_range_low, n_range_high+1):
      if nn in used_ranges: continue
      used_ranges.append(nn)
      if density_threshold is not None:
        threshold = density_threshold
      else:
        threshold = starting_density_threshold*(scale**nn)
      if threshold < lower_bound or threshold > upper_bound:
        continue

      co, sorted_by_volume, min_b, max_b = get_co(
        map_data = map_data.deep_copy(),
         threshold = threshold, wrapping = wrapping)

      if len(sorted_by_volume)<2:
        score, has_sufficient_regions, too_low, too_high, expected_regions, ok = \
          None, None, None, None, None, None
        continue # don't go on
      else:
        score, has_sufficient_regions, too_low, too_high, expected_regions, ok = \
           score_threshold(b_vs_region = b_vs_region,
         threshold = threshold,
         sorted_by_volume = sorted_by_volume,
         fraction_occupied = fraction_occupied,
         solvent_fraction = solvent_fraction,
         residues_per_region = residues_per_region,
         min_volume = min_volume,
         min_ratio = min_ratio,
         max_ratio_to_target = max_ratio_to_target,
         min_ratio_to_target = min_ratio_to_target,
         min_ratio_of_ncs_copy_to_first = min_ratio_of_ncs_copy_to_first,
         ncs_copies = ncs_copies,
         n_residues = n_residues,
         map_data = map_data,
         target_in_all_regions = target_in_all_regions,
         crystal_symmetry = crystal_symmetry,
         chain_type = chain_type,
         out = local_out)
      if expected_regions:
        unique_expected_regions = max(1,
         (ncs_copies-1+expected_regions)//ncs_copies)
      if too_high and threshold<upper_bound:
        upper_bound = threshold
      if too_low and threshold>lower_bound:
        lower_bound = threshold
      if score is None:
        if best_threshold and best_threshold_has_sufficient_regions:
          if threshold >best_threshold: # new upper bound
           upper_bound = threshold
          elif threshold <best_threshold: # new lower bound
           lower_bound = threshold
      elif  (ok or not best_ok) and  \
            (best_score is None or score > best_score):
        best_threshold = threshold
        best_threshold_has_sufficient_regions = has_sufficient_regions
        best_score = score
        best_ok = ok

  if best_threshold is not None:
    print("\nBest threshold: %5.2f\n" %(best_threshold), file = out)
    return best_threshold, unique_expected_regions, best_score, best_ok
  elif density_threshold is not None: # use it anyhow
    return density_threshold, unique_expected_regions, None, None
  else:
    return None, unique_expected_regions, None, None

def get_co(map_data = None, threshold = None, wrapping = None):
  """Run the maptbx.connectivity class to get a connectivity object
  containing regions (lists of grid indices) that are connected at this
  threshold value"""
  co = maptbx.connectivity(map_data = map_data, threshold = threshold,
         wrapping = wrapping)
  regions = co.regions()
  rr = list(range(0, co.regions().size()))

  regions_0 = regions[0]
  rr_0 = rr[0]
  regions = regions[1:]
  rr = rr[1:]
  if rr:
    z = zip(regions, rr)
    sorted_by_volume = sorted(z, key=itemgetter(0), reverse = True)
  else:
    sorted_by_volume = []
  sorted_by_volume = [(regions_0, rr_0)]+sorted_by_volume

  min_b, max_b = co.get_blobs_boundaries_tuples() # As grid points, not A
  return co, sorted_by_volume, min_b, max_b

def get_connectivity(b_vs_region = None,
     map_data = None,
     solvent_fraction = None,
     n_residues = None,
     ncs_copies = None,
     fraction_occupied = None,
     iterate_with_remainder = None,
     min_volume = None,
     min_ratio = None,
     wrapping = None,
     residues_per_region = None,
     max_ratio_to_target = None,
     min_ratio_to_target = None,
     min_ratio_of_ncs_copy_to_first = None,
     starting_density_threshold = None,
     density_threshold = None,
     crystal_symmetry = None,
     chain_type = None,
     verbose = None,
     out = sys.stdout):
  """Get the connectivity in a map.
  Optimizes threshold for calculation of connectivity.
  Returns a connectivity object co, sorted_by_volume, min_b, max_b,
      best_unique_expected_regions,
      best_score, threshold, starting_density_threshold"""
  print("\nGetting connectivity", file = out)
  libtbx.call_back(message = 'segment', data = None)


  # Normalize map data now to SD of the part that is not solvent
  map_data = renormalize_map_data(
     map_data = map_data, solvent_fraction = solvent_fraction)

  # Try connectivity at various thresholds
  # Choose one that has about the right number of grid points in top regions
  scale = 0.95
  best_threshold = None
  best_scale = scale
  best_score = None
  best_ok = None
  best_unique_expected_regions = None
  for ii in range(3):
    threshold, unique_expected_regions, score, ok = choose_threshold(
     density_threshold = density_threshold,
     starting_density_threshold = starting_density_threshold,
     b_vs_region = b_vs_region,
     map_data = map_data,
     n_residues = n_residues,
     ncs_copies = ncs_copies,
     fraction_occupied = fraction_occupied,
     solvent_fraction = solvent_fraction,
     scale = scale,
     wrapping = wrapping,
     residues_per_region = residues_per_region,
     min_volume = min_volume,
     min_ratio = min_ratio,
     max_ratio_to_target = max_ratio_to_target,
     min_ratio_to_target = min_ratio_to_target,
     min_ratio_of_ncs_copy_to_first = min_ratio_of_ncs_copy_to_first,
     crystal_symmetry = crystal_symmetry,
     chain_type = chain_type,
     verbose = verbose,
     out = out)
    # Take it if it improves (score, ok)
    if threshold is not None:
     if best_score is None or  \
      ((ok or not best_ok) and (score > best_score)):
      best_score = score
      best_unique_expected_regions = unique_expected_regions
      best_ok = ok
      best_threshold = threshold
      best_scale = scale
    if best_ok or density_threshold is not None:
      break
    else:
      scale = scale**0.333 # keep trying

  if best_threshold is None or (
      density_threshold is not None and best_score is None):
    if iterate_with_remainder: # on first try failed
      raise Sorry("No threshold found...try with density_threshold = xxx")
    else: # on iteration...ok
      print("Note: No threshold found", file = out)
      return None, None, None, None, None, None, None, None
  else:
    starting_density_threshold = best_threshold
    # try it next time

  co, sorted_by_volume, min_b, max_b = get_co(
    map_data = map_data, threshold = best_threshold, wrapping = wrapping)

  return co, sorted_by_volume, min_b, max_b, best_unique_expected_regions, \
      best_score, threshold, starting_density_threshold

def get_volume_of_seq(text, chain_type = None, out = sys.stdout):
  """Estimate partial molar volume occupied by a molecule with this
  sequence and chain_type"""
  from iotbx.bioinformatics import chain_type_and_residues
  # get chain type and residues (or use given chain type and count residues)
  chain_type, n_residues = chain_type_and_residues(text = text, chain_type = chain_type)
  if chain_type is None and n_residues is None:
    return None, None, None
  if chain_type == 'PROTEIN':
    mw_residue = 110.0  # from $CDOC/matthews.doc
    density_factor = 1.23   # 1.66/DENSITY-OF-PROTEIN = 1.66/1.35
  else:
    mw_residue = 330.0  # guess for DNA/RNA
    density_factor = 1.15   # 1.66/DENSITY-OF-DNA = 1.66/1.45
  return len(text)*density_factor*mw_residue, len(text), chain_type

def get_solvent_content_from_seq_file(params,
    sequence = None,
    seq_file = None,
    overall_chain_type = None,
    ncs_copies = None,
    map_volume = None,
    out = sys.stdout):

  """Estimate solvent content from map volume, sequence or sequence file,
    and ncs_copies"""
  if params and not overall_chain_type:
     overall_chain_type = params.crystal_info.chain_type

  if not sequence and not os.path.isfile(seq_file):
    raise Sorry(
     "The sequence file '%s' is missing." %(seq_file))
  if not sequence:
    print("\nReading sequence from %s " %(seq_file), file = out)
    sequence = open(seq_file).read()
  from iotbx.bioinformatics import get_sequences
  sequences = get_sequences(text = sequence)
  # get unique part of these sequences

  from mmtbx.validation.chain_comparison import \
       extract_unique_part_of_sequences as eups
  print("Unique part of sequences:", file = out)
  copies_in_unique, base_copies, unique_sequence_dict = eups(sequences,
     out = out)
  all_unique_sequence = []
  for seq in copies_in_unique.keys():
    print("Copies: %s  base copies: %s  Sequence: %s" %(
       copies_in_unique[seq], base_copies, seq), file = out)
    all_unique_sequence.append(seq)
  if base_copies !=  ncs_copies:
    print("NOTE: %s copies of unique portion but ncs_copies = %s" %(
       base_copies, ncs_copies), file = out)
    if ncs_copies == 1:
      ncs_copies = base_copies
      print("Using ncs_copies = %s instead" %(ncs_copies), file = out)
    else:
      print("Still using ncs_copies = %s" %(ncs_copies), file = out)

  volume_of_chains = 0.
  n_residues = 0
  chain_types_considered = []
  for seq in all_unique_sequence:
    volume, nres, chain_type = get_volume_of_seq(seq,
      chain_type = overall_chain_type, out = out)
    if volume is None: continue
    volume_of_chains+= volume
    n_residues+= nres
    if not chain_type in chain_types_considered:
      chain_types_considered.append(chain_type)
  chain_types_considered.sort()
  print("\nChain types considered: %s\n" %(
      " ".join(chain_types_considered)), file = out)
  volume_of_molecules = volume_of_chains*ncs_copies
  n_residues_times_ncs = n_residues*ncs_copies
  solvent_fraction = 1.-(volume_of_molecules/map_volume)
  solvent_fraction = max(0.001, min(0.999, solvent_fraction))
  if solvent_fraction == 0.001 or solvent_fraction == 0.999:
    print("NOTE: solvent fraction of %7.2f very unlikely..." %(
        solvent_fraction) + "please check ncs_copies and sequence ", file = out)
  print("Solvent content from composition: %7.2f" %(solvent_fraction), file = out)
  print("Cell volume: %.1f  NCS copies: %d   Volume of unique chains: %.1f" %(
     map_volume, ncs_copies, volume_of_chains), file = out)
  print("Total residues: %d  Volume of all chains: %.1f  Solvent fraction: %.3f "%(
       n_residues_times_ncs, volume_of_molecules, solvent_fraction), file = out)
  return solvent_fraction, n_residues, n_residues_times_ncs

def get_solvent_fraction(params,
     ncs_object = None, ncs_copies = None,
     do_not_adjust_dalton_scale = None,
     sequence = None,
     seq_file = None,
     molecular_mass = None,
     solvent_content = None,
     crystal_symmetry = None, tracking_data = None, out = sys.stdout):
  """Estimate solvent fraction from params, tracking_data,
     sequence, seq_file, or molecular mass.
  Returns tracking_data object if one is supplied, otherwise returns
  solvent_content"""

  if tracking_data and not crystal_symmetry:
    #crystal_symmetry = tracking_data.original_crystal_symmetry not used
    crystal_symmetry = tracking_data.crystal_symmetry
  map_volume = crystal_symmetry.unit_cell().volume()
  if tracking_data and not ncs_copies:
    #ncs_copies = tracking_data.input_ncs_info.original_number_of_operators
    ncs_copies = tracking_data.input_ncs_info.number_of_operators # 2018-01-29 put back
  if not ncs_copies: ncs_copies = 1


  if params and not solvent_content:
    solvent_content = params.crystal_info.solvent_content
  if params and not molecular_mass:
    molecular_mass = params.crystal_info.molecular_mass
  if params and not seq_file:
    seq_file = params.input_files.seq_file
  if params and not sequence:
    sequence = params.crystal_info.sequence

  if seq_file or sequence:
    solvent_content, n_residues, n_residues_times_ncs = \
         get_solvent_content_from_seq_file(
     params,
     sequence = sequence,
     seq_file = seq_file,
     ncs_copies = ncs_copies,
     map_volume = map_volume,
     out = out)
    if params and not params.crystal_info.solvent_content:
      params.crystal_info.solvent_content = solvent_content
      print("Solvent fraction from composition: %7.2f "%(
       params.crystal_info.solvent_content), file = out)
    elif params:
      print("Solvent content from parameters: %7.2f" %(
        params.crystal_info.solvent_content), file = out)

  else:
    if params and params.crystal_info.solvent_content:
      print("Solvent content from parameters: %7.2f" %(
        params.crystal_info.solvent_content), file = out)
    elif molecular_mass:
       solvent_content = \
         get_solvent_fraction_from_molecular_mass(
        crystal_symmetry = crystal_symmetry,
        do_not_adjust_dalton_scale = do_not_adjust_dalton_scale,
        molecular_mass = molecular_mass,
        out = out)
       if params:
         params.crystal_info.solvent_content = solvent_content

    else:
      print("Getting solvent content automatically.", file = out)

  if tracking_data:
    if params.input_files.seq_file or params.crystal_info.sequence:
      tracking_data.set_input_seq_info(file_name = params.input_files.seq_file,
       sequence = params.crystal_info.sequence,
        n_residues = n_residues)
      tracking_data.set_n_residues(
        n_residues = n_residues_times_ncs)
    if params.crystal_info.solvent_content:
      tracking_data.set_solvent_fraction(params.crystal_info.solvent_content)

    return tracking_data
  else:
    return solvent_content

def top_key(dd):
  """Return key from dict dd that has the largest value of dd[key]"""
  if not dd:
    return None, None
  elif len(dd) == 1:
    return list(dd.items())[0]
  else:
    best_key = None
    best_n = None
    for key in dd.keys():
      if not best_n or dd[key] > best_n:
        best_n = dd[key]
        best_key = key
    return best_key, best_n

def choose_max_regions_to_consider(params,
    sorted_by_volume = None,
    ncs_copies = None):

  """Sort and eliminate regions with few points and those at end of list"""
  max_per_au = params.segmentation.max_per_au
  min_ratio = params.segmentation.min_ratio
  min_volume = params.segmentation.min_volume
  if len(sorted_by_volume)<2:
    return 0
  max_grid_points = sorted_by_volume[1][0]
  cntr = 0
  for p in sorted_by_volume[1:]:
    cntr+= 1
    if max_per_au and (cntr>max_per_au*ncs_copies):
      cntr-= 1
      break
    v, i = p  # v = volume in grid points, i = id
    if v/max_grid_points<min_ratio or v < min_volume:
      cntr-= 1
      break
  return cntr

def get_edited_mask(sorted_by_volume = None,
    max_regions_to_consider = None,
    co = None,
    out = sys.stdout):
  """Analyze connectivity object co and create an edited mask with values at
  each grid point that specify which region it is associated with"""
  conn_obj = co.result()
  origin = list(conn_obj.accessor().origin())
  all = list(conn_obj.accessor().all())
  conn_obj.accessor().show_summary(out)
  edited_mask = conn_obj.deep_copy()
  first = True
  edited_volume_list = []
  original_id_from_id = {}
  for i in range(1, max_regions_to_consider+1):
    v, id = sorted_by_volume[i]
    original_id_from_id[i] = id
    edited_volume_list.append(v)
    s = (conn_obj == id)
    if first:
      edited_mask = edited_mask.set_selected(~s, 0)
      first = False
    edited_mask = edited_mask.set_selected(s, i)   # edited mask has ID of
         # regions, labeled in decreasing size from 1 to max_regions_to_consider
  return edited_mask, edited_volume_list, original_id_from_id

def choose_subset(a, target_number = 1):
  """Choose target_number of evenly-spaced elements of a"""
  new_array = flex.vec3_double()
  assert type(new_array) == type(a)
  n = a.size()
  nskip = max(1, n//target_number)
  i = 0
  for x in a:
    if i%nskip == 0 or i == n-1:
     new_array.append(x)
    i+= 1
  return new_array

def run_get_duplicates_and_ncs(
   ncs_obj = None,
   min_b = None,
   max_b = None,
   edited_mask = None,
   original_id_from_id = None,
   edited_volume_list = None,
   max_regions_to_consider = None,
   regions_left = None,
   tracking_data = None,
   out = sys.stdout,
   ):
  """Run the get_duplicates_and_ncs method and make sure
   we have region_centroid for all values by varying max_regions_to_consider.
  Returns duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict,
        region_range_dict, region_centroid_dict,
        region_scattered_points_dict """

  duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict, region_range_dict, \
     region_centroid_dict, region_scattered_points_dict = \
      get_duplicates_and_ncs(
        ncs_obj = ncs_obj,
        min_b = min_b,
        max_b = max_b,
        edited_mask = edited_mask,
        edited_volume_list = edited_volume_list,
        original_id_from_id = original_id_from_id,
        max_regions_to_consider = max_regions_to_consider,
        tracking_data = tracking_data,
        out = out)

  # check that we have region_centroid for all values
  complete = True
  missing = []
  for i in range(1, max_regions_to_consider+1):
    if not i in region_centroid_dict.keys():
      if (regions_left is None) or (i in regions_left):
         complete = False
         missing.append(i)
  if complete:
       return duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict, \
        region_range_dict, region_centroid_dict, \
        region_scattered_points_dict
  else:
    raise Sorry("Cannot find region-centroid for all regions? Missing: %s" %(
      missing))

def copy_dict_info(from_dict, to_dict):
  """Copy from_dict to to_dict.  Use instead to_dict = from_dict.copy()"""
  for key in from_dict.keys():
    to_dict[key] = from_dict[key]

def get_centroid_from_blobs(min_b = None, max_b = None,
    id = None, original_id_from_id = None):
  """Get the centroid of a blob with bounds min_b, max_b"""
  orig_id = original_id_from_id[id]
  upper = max_b[orig_id]
  lower = min_b[orig_id]
  avg = []
  for u, l in zip(upper, lower):
    avg.append(0.5*(u+l))
  return avg

def get_duplicates_and_ncs(
   ncs_obj = None,
   min_b = None,
   max_b = None,
   edited_mask = None,
   original_id_from_id = None,
   edited_volume_list = None,
   max_regions_to_consider = None,
   target_points_per_region = 30,
   minimum_points_per_region = 10,
   maximum_points_per_region = 100,
   tracking_data = None,
   out = sys.stdout,
   ):

  """
  Analyze a list of regions in edited_volume_list to identify NCS and
  duplicates. Use the NCS relationships in ncs_obj.

  Return: duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict, \
      region_range_dict, region_centroid_dict, region_scattered_points_dict
  """
  unit_cell = tracking_data.crystal_symmetry.unit_cell()

  # Get scattered points in each region
  region_scattered_points_dict = get_region_scattered_points_dict(
     edited_volume_list = edited_volume_list,
     edited_mask = edited_mask,
     unit_cell = unit_cell,
     target_points_per_region = target_points_per_region,
     minimum_points_per_region = minimum_points_per_region,
     maximum_points_per_region = maximum_points_per_region)


  # Now just use the scattered points to get everything else:
  region_n_dict = {}  # count of points used by region (differs from volume due
     # to the sampling)
  region_range_dict = {} # keyed by region in edited_mask; range for x, y, z
  region_centroid_dict = {} # keyed by region in edited_mask; range for x, y, z
  for id in region_scattered_points_dict.keys():
    sites = region_scattered_points_dict[id]
    region_n_dict[id] = sites.size()
    if region_n_dict[id]:
      region_centroid_dict[id] = list(sites.mean())
    else: # No points...use bounds from object
      region_centroid_dict[id] = get_centroid_from_blobs(min_b = min_b,
        max_b = max_b,
        id = id, original_id_from_id = original_id_from_id)

  # Now get NCS relationships

  ncs_group = ncs_obj.ncs_groups()[0]
  duplicate_dict = {}  # keyed by id, number of duplicates for that region
  equiv_dict = {}  # equiv_dict[id][other_id] = number_of points other_id matches
                 #  id through an ncs relationship
  equiv_dict_ncs_copy_dict = {}
  for id in region_scattered_points_dict.keys():
    duplicate_dict[id] = 0
    equiv_dict[id] = {}
    equiv_dict_ncs_copy_dict[id] = {}

  # Figure out which ncs operator is the identity
  identity_op = ncs_group.identity_op_id()
  print("Identity operator is %s" %(identity_op), file = out)

  # 2017-12-16 Score poorly if it involves a cell translation unless it
  #  is a crystal

  if len(ncs_group.translations_orth())>1:
    # Skip if no ncs...
    for id in region_scattered_points_dict.keys():
      for xyz_cart in region_scattered_points_dict[id]:
        n = 0
        for i0 in range(len(ncs_group.translations_orth())):
          if i0 == identity_op: continue
          r = ncs_group.rota_matrices_inv()[i0] # inverse maps pos 0 on to pos i
          t = ncs_group.translations_orth_inv()[i0]

          n+= 1
          new_xyz_cart = r * matrix.col(xyz_cart) + t
          new_xyz_frac = unit_cell.fractionalize(new_xyz_cart)
          if tracking_data.params.crystal_info.use_sg_symmetry or \
            (new_xyz_frac[0]>= 0 and new_xyz_frac[0]<= 1 and \
             new_xyz_frac[1]>= 0 and new_xyz_frac[1]<= 1 and \
             new_xyz_frac[2]>= 0 and new_xyz_frac[2]<= 1):
            value = edited_mask.value_at_closest_grid_point(new_xyz_frac)
          else:
            value = 0  # value for nothing there 2017-12-16
          if value == id:
            duplicate_dict[id]+= 1
            break # only count once
          elif value>0:  # notice which one is matched
            if not value in equiv_dict[id]:
              equiv_dict[id][value] = 0
              equiv_dict_ncs_copy_dict[id][value] = {}
            equiv_dict[id][value]+= 1
            if not n in equiv_dict_ncs_copy_dict[id][value]:
              equiv_dict_ncs_copy_dict[id][value][n] = 0
            equiv_dict_ncs_copy_dict[id][value][n]+= 1  # how many are ncs copy n
  return duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict, \
      region_range_dict, region_centroid_dict, region_scattered_points_dict

def get_region_scattered_points_dict(
  edited_volume_list = None,
  edited_mask = None,
  unit_cell = None,
  sampling_rate = None,
  target_points_per_region = None,
  minimum_points_per_region = None,
  maximum_points_per_region = None):

  """Get a sample of points in each region"""
  sample_dict = {}
  region_scattered_points_dict = {} # some points in each region
  if not sampling_rate:
    sampling_rate = edited_volume_list[0]//target_points_per_region
    sampling_rate_set = False
  else:
    sampling_rate_set = True
  volumes = flex.int()
  sampling_rates = flex.int()
  id_list = []
  # have to set up dummy first set:
  volumes.append(0)
  sampling_rates.append(0)
  id_list.append(0)

  for i in range(len(edited_volume_list)):
    id = i+1
    v = edited_volume_list[i]

    region_scattered_points_dict[id] = flex.vec3_double()

    volumes.append(v)
    if sampling_rate_set:
      sample_dict[id] = sampling_rate
      sampling_rates.append(sampling_rate)
    else:
      sample_dict[id] = max(1,
        max(v//maximum_points_per_region,
          min(v//minimum_points_per_region,
              sampling_rate)  ))
      sampling_rates.append(max(1,
      max(v//maximum_points_per_region,
          min(v//minimum_points_per_region,
              sampling_rate)  )))
    id_list.append(id)

  sample_regs_obj = maptbx.sample_all_mask_regions(
      mask = edited_mask,
      volumes = volumes,
      sampling_rates = sampling_rates,
      unit_cell = unit_cell)

  for id in id_list[1:]:  # skip the dummy first set
    region_scattered_points_dict[id] = sample_regs_obj.get_array(id)

  return region_scattered_points_dict

def remove_bad_regions(params = None,
  duplicate_dict = None,
  edited_volume_list = None,
  out = sys.stdout):
  """Remove regions that are in duplicate_dict and have more than
    max_overlap_fraction overlap.
  Returns region_list, region_volume_dict, new_sorted_by_volume, bad_region_list
  """
  worst_list = []
  for id in list(duplicate_dict.keys()):
    fract = duplicate_dict[id]/edited_volume_list[id-1]
    if duplicate_dict[id] and fract >= params.segmentation.max_overlap_fraction:
      worst_list.append([fract, id])
    else:
      del duplicate_dict[id]
  worst_list.sort()
  worst_list.reverse()

  bad_region_list = []
  max_number_to_remove = int(0.5+
    0.01*params.segmentation.remove_bad_regions_percent*len(edited_volume_list))
  if worst_list:
    print("\nRegions that span multiple NCS au:", file = out)
    for fract, id in worst_list:
      print("ID: %d  Duplicate points: %d (%.1f %%)" %(
        id, duplicate_dict[id], 100.*fract), end = ' ', file = out)
      if  len(bad_region_list)<max_number_to_remove:
         bad_region_list.append(id)
         print(" (removed)", file = out)
      else:
         print(file = out)

  new_sorted_by_volume = []
  region_list = []
  region_volume_dict = {}
  for i in range(len(edited_volume_list)):
    id = i+1
    v = edited_volume_list[i]
    new_sorted_by_volume.append([v, id])
    region_list.append(id)
    region_volume_dict[id] = v
  if bad_region_list:
    print("Bad regions (excluded)", bad_region_list, file = out)
  return region_list, region_volume_dict, new_sorted_by_volume, bad_region_list

def sort_by_ncs_overlap(matches, equiv_dict_ncs_copy_dict_id):
    """Sort matches (pairs of region IDs) by NCS overlap"""
    sort_list = []
    for id1 in matches:
      key, n = top_key(equiv_dict_ncs_copy_dict_id[id1]) # Take top ncs_copy
      sort_list.append([n, id1])
    sort_list.sort(key=itemgetter(0))
    sort_list.reverse()
    key_list = []
    for n, id1 in sort_list:
      key_list.append(id1)
    return key_list


def get_ncs_equivalents(
    bad_region_list = None,
    region_list = None,
    region_scattered_points_dict = None,
    equiv_dict = None,
    ncs_copies = None,
    equiv_dict_ncs_copy_dict = None,
    min_coverage = .10,
    out = sys.stdout):

  """Identify regions that are NCS-related.
  Returns  equiv_dict_ncs_copy"""
  equiv_dict_ncs_copy = {}
  for id in region_list:
    if id in bad_region_list: continue
    match_dict = equiv_dict.get(id, {}) # which are matches
    matches = list(match_dict.keys())
    if not matches: continue
    key_list = sort_by_ncs_overlap(matches, equiv_dict_ncs_copy_dict[id])
    n_found = 0
    for id1 in key_list:
      #     id matches id1 N = match_dict[id1]
      #  2017-12-16 Do not include if there is a cell translation

      key, n = top_key(equiv_dict_ncs_copy_dict[id][id1]) # ncs_copy, n-overlap
      if n<min_coverage*region_scattered_points_dict[id].size():
        break
      else:
        if not id in equiv_dict_ncs_copy:equiv_dict_ncs_copy[id] = {}
        equiv_dict_ncs_copy[id][id1] = key
        n_found+= 1
        if n_found>= ncs_copies-1:
          break

  return equiv_dict_ncs_copy

def get_overlap(l1, l2):
  """Return list of unique members of l2 that are present in l1"""
  overlap_list = []
  l1a = single_list(l1)
  l2a = single_list(l2)
  for i in l1a:
    if i in l2a and not i in overlap_list: overlap_list.append(i)
  return overlap_list

def group_ncs_equivalents(params,
    region_list = None,
    region_volume_dict = None,
    equiv_dict_ncs_copy = None,
    tracking_data = None,
    split_if_possible = None,
    out = sys.stdout):

  """Group NCS equivalent regions together
  Group together all the regions that are related to region 1...etc
   if split_if_possible then skip all groups with multiple entries
   equiv_dict_ncs_copy[id][id1] = ncs_copy
  Returns: ncs_group_list, shared_group_dict"""

  ncs_equiv_groups_as_list = []
  ncs_equiv_groups_as_dict = {}
  for id in region_list:
    equiv_group = {}  #equiv_group[ncs_copy] = [id1, id2, id3...]
    equiv_group[0] = [id] # always
    for id1 in equiv_dict_ncs_copy.get(id, {}).keys():
      ncs_copy = equiv_dict_ncs_copy[id][id1]
      if not ncs_copy in equiv_group: equiv_group[ncs_copy] = []
      equiv_group[ncs_copy].append(id1) # id1 is ncs_copy of id
    all_single = True
    equiv_group_as_list = []
    total_grid_points = 0
    missing_ncs_copies = []
    present_ncs_copies = []
    for ncs_copy in range(tracking_data.input_ncs_info.number_of_operators):
        # goes 0 to ncs_copies-1 (including extra ones if present)
      local_equiv_group = equiv_group.get(ncs_copy, [])
      if local_equiv_group:
        equiv_group_as_list.append(local_equiv_group)
        present_ncs_copies.append(ncs_copy)
        if ncs_copy > 0 and \
          len(local_equiv_group)>1 and len(equiv_group.get(0, [])) == 1:
          all_single = False
        for id in equiv_group.get(ncs_copy, []):
          total_grid_points+= region_volume_dict[id]
      else:
        missing_ncs_copies.append(ncs_copy)
    equiv_group_as_list.sort()
    if tracking_data.input_ncs_info.is_helical_symmetry:
      # complete if we have original_number_of_operators worth
      if (not params.segmentation.require_complete) or \
         len(present_ncs_copies)>=  \
         tracking_data.input_ncs_info.original_number_of_operators:
        complete = True
      else:
        complete = False
    else:
      if len(missing_ncs_copies) == 0:
        complete = True
      else:
        complete = False
    if complete and \
        (not str(equiv_group_as_list) in ncs_equiv_groups_as_dict or
         total_grid_points>ncs_equiv_groups_as_dict[str(equiv_group_as_list)]) \
        and (all_single or (not split_if_possible)):
      ncs_equiv_groups_as_dict[str(equiv_group_as_list)] = total_grid_points
      ncs_equiv_groups_as_list.append([total_grid_points, equiv_group_as_list])

  ncs_equiv_groups_as_list.sort()
  ncs_equiv_groups_as_list.reverse()

  # Now remove any group that duplicates a previous group
  # 2015-11-07 allow a member to be in multiple groups though (for example
  #   one that spans several groups because it contains 2 region in other ncs
  #   copies)
  #  Make sure that if there are duplicates they are all in the leading
  #    positions of the list (these must be very big ones as they match 2
  #    regions in other ncs copies)

  max_duplicates = tracking_data.input_ncs_info.number_of_operators-1 # not all duplicates
  ncs_group_list = []
  used_list = []
  print("All equiv groups:", file = out)
  used_regions = []
  for total_grid_points, equiv_group_as_list in ncs_equiv_groups_as_list:
    duplicate = False
    n_dup = 0
    for equiv_group in equiv_group_as_list:
      for x in equiv_group:
        if x in used_list:
          n_dup+= 1
    if n_dup>max_duplicates or n_dup >len(equiv_group_as_list)-1:
      duplicate = True
    if not duplicate and n_dup>0:  # check carefully to make sure that all
      # are leading entries
      for ncs_group in ncs_group_list:
        overlaps = get_overlap(ncs_group, equiv_group_as_list)
        if not overlaps: continue
        overlaps.sort()
        expected_match = single_list(equiv_group_as_list)[:len(overlaps)]
        expected_match.sort()
        if overlaps!= expected_match: # not leading entries
          duplicate = True
          break

    if not duplicate:
      #print >>out, "NCS GROUP:", equiv_group_as_list, ":", total_grid_points

      ncs_group_list.append(equiv_group_as_list)
      for equiv_group in equiv_group_as_list:
        for x in equiv_group:
          if not x in used_list: used_list.append(x)
  print("Total NCS groups: %d" %len(ncs_group_list), file = out)

  # Make a dict that lists all ids that are in the same group as region x
  shared_group_dict = {}
  for ncs_group in ncs_group_list:
    for group_list in ncs_group:
      for id1 in group_list:
        if not id1 in shared_group_dict: shared_group_dict[id1] = []
        for other_group_list in ncs_group:
          if other_group_list is group_list:continue
          for other_id1 in other_group_list:
            if not other_id1 in shared_group_dict [id1]:
              shared_group_dict[id1].append(other_id1)

  return ncs_group_list, shared_group_dict

def identify_ncs_regions(params,
     sorted_by_volume = None,
     co = None,
     min_b = None,
     max_b = None,
     ncs_obj = None,
     tracking_data = None,
     out = sys.stdout):

  """Identify NCS regions in map
  1.choose top regions to work with
  2.remove regions that are in more than one au of the NCS
  3.identify groups of regions that are related by NCS
  Also note the centers and bounds of each region.
  """

  # Choose number of top regions to consider

  max_regions_to_consider = choose_max_regions_to_consider(params,
    sorted_by_volume = sorted_by_volume,
    ncs_copies = tracking_data.input_ncs_info.original_number_of_operators)

  print("\nIdentifying NCS-related regions.Total regions to consider: %d" %(
    max_regions_to_consider), file = out)
  if max_regions_to_consider<1:
    print("\nUnable to identify any NCS regions", file = out)
    return None, tracking_data, None

  # Go through all grid points; discard if not in top regions
  #  Renumber regions in order of decreasing size

  load_saved_files = False  # set to True to load results from previous run
  dump_files = False        # set to True to dump results and speed up next run
  if not load_saved_files:
    edited_mask, edited_volume_list, original_id_from_id = get_edited_mask(
     sorted_by_volume = sorted_by_volume,
     co = co,
     max_regions_to_consider = max_regions_to_consider, out = out)
    if dump_files:
      from libtbx import easy_pickle
      easy_pickle.dump("edited_mask.pkl",
        [edited_mask, edited_volume_list, original_id_from_id])
  else:
    from libtbx import easy_pickle
    [edited_mask, edited_volume_list, original_id_from_id
        ] = easy_pickle.load("edited_mask.pkl")
    print("Loading edited_mask.pkl", file = out)

  # edited_mask contains re-numbered region id's

  # Identify duplicate and ncs relationships between regions
  # duplicate_dict[id] =  number of duplicates for that region
  # equiv_dict[id][other_id] = number_of points other_id matches
                   #  id through an ncs relationship
  if not load_saved_files:
    duplicate_dict, equiv_dict, equiv_dict_ncs_copy_dict, \
      region_range_dict, region_centroid_dict, \
      region_scattered_points_dict = \
    run_get_duplicates_and_ncs(
      ncs_obj = ncs_obj,
      min_b = min_b,
      max_b = max_b,
      edited_mask = edited_mask,
      original_id_from_id = original_id_from_id,
      edited_volume_list = edited_volume_list,
      max_regions_to_consider = max_regions_to_consider,
      tracking_data = tracking_data,
      out = out)
    # Remove any bad regions
    region_list, region_volume_dict, new_sorted_by_volume, \
      bad_region_list = remove_bad_regions(
    params = params,
    duplicate_dict = duplicate_dict,
    edited_volume_list = edited_volume_list,
    out = out)
    # Identify groups of regions that are ncs-related
    # equiv_dict_ncs_copy[id][id1] = ncs_copy of id that corresponds to id1
    equiv_dict_ncs_copy = get_ncs_equivalents(
    region_list = region_list,
    bad_region_list = bad_region_list,
    region_scattered_points_dict = region_scattered_points_dict,
    equiv_dict = equiv_dict,
    ncs_copies = tracking_data.input_ncs_info.number_of_operators,
    equiv_dict_ncs_copy_dict = equiv_dict_ncs_copy_dict,
    out = out)
    if dump_files:
      from libtbx import easy_pickle
      easy_pickle.dump("save.pkl", [duplicate_dict, equiv_dict, region_range_dict, region_centroid_dict, region_scattered_points_dict, region_list, region_volume_dict, new_sorted_by_volume, bad_region_list, equiv_dict_ncs_copy, tracking_data])
      print("Dumped save.pkl", file = out)
  else:
    from libtbx import easy_pickle
    [duplicate_dict, equiv_dict, region_range_dict, region_centroid_dict, region_scattered_points_dict, region_list, region_volume_dict, new_sorted_by_volume, bad_region_list, equiv_dict_ncs_copy, tracking_data] = easy_pickle.load("save.pkl")
    print("Loaded save.pkl", file = out)

  # Group together regions that are ncs-related. Also if one ncs
  #   copy has 2 or more regions linked together, group the other ones.

  # each entry in ncs_group_list is a list of regions for each ncs_copy:
  #  e.g.,  [[8], [9, 23], [10, 25], [11, 27], [12, 24], [13, 22], [14, 26]]
  #  May contain elements that are in bad_region_list (to exclude later)
  if not load_saved_files:
    ncs_group_list, shared_group_dict = group_ncs_equivalents(params,
    split_if_possible = params.segmentation.split_if_possible,
    tracking_data = tracking_data,
    region_volume_dict = region_volume_dict,
    region_list = region_list,
    equiv_dict_ncs_copy = equiv_dict_ncs_copy,
    out = out)
    if dump_files:
      from libtbx import easy_pickle
      easy_pickle.dump("group_list.pkl", [ncs_group_list, shared_group_dict])
      print("Dumped to group_list.pkl", file = out)
  else:
    from libtbx import easy_pickle
    [ncs_group_list, shared_group_dict] = easy_pickle.load("group_list.pkl")
    print("Loaded group_list.pkl", file = out)

  ncs_group_obj = ncs_group_object(
     ncs_group_list = ncs_group_list,
     shared_group_dict = shared_group_dict,
     ncs_obj = ncs_obj,
     crystal_symmetry = tracking_data.crystal_symmetry,
     edited_mask = edited_mask,
     origin_shift = tracking_data.origin_shift,
     co = co,
     min_b = min_b,
     max_b = max_b,
     equiv_dict = equiv_dict,
     bad_region_list = bad_region_list,
     original_id_from_id = original_id_from_id,
     edited_volume_list = edited_volume_list,
     region_range_dict = region_range_dict,
     region_scattered_points_dict = region_scattered_points_dict,
     region_centroid_dict = region_centroid_dict)

  return ncs_group_obj, tracking_data, equiv_dict_ncs_copy

def get_center_list(regions,
    region_centroid_dict = None):
  """Get list of centers of each region"""
  center_list = []
  for region in regions:
    center_list.append(region_centroid_dict[region])
  return center_list

def get_average_center(regions,
    region_centroid_dict = None):
  """Get average center of all regions"""
  center_list = get_center_list(regions, region_centroid_dict = region_centroid_dict)
  for region in regions:
    center_list.append(region_centroid_dict[region])
  average_center = deepcopy(center_list[0])
  if len(center_list)>1:
    for r in center_list[1:]:
      for i in range(3):
        average_center[i]+= r[i]
    for i in range(3):
      average_center[i]/= len(center_list)
  return average_center

def get_dist(r, s):
  """Get distance between r and s"""
  dd = 0.
  for i in range(3):
    dd+= (r[i]-s[i])**2
  return dd**0.5

def has_intersection(set1, set2):
  """Return True if intersection of set1 and set2 is not empty"""
  set1a = single_list(set1)
  set2a = single_list(set2)
  for x in set1a:
    if x in set2a:
      return True
  return False

def get_scattered_points_list(other_regions,
       region_scattered_points_dict = None):
  """Get scattered_points_list from other_regions
   using region_scattered_points_dict"""
  scattered_points_list = flex.vec3_double()
  for x in other_regions:
    scattered_points_list.extend(region_scattered_points_dict[x])
  return scattered_points_list

def get_inter_region_dist_dict(ncs_group_obj = None,
    selected_regions = None, target_scattered_points = None):
  """Get dictionary with inter-region distances"""
  dd = {}
  for i in range(len(selected_regions)):
    id = selected_regions[i]
    if not id in dd: dd[id] = {}
    test_centers = ncs_group_obj.region_scattered_points_dict[id]
    for j in range(i+1, len(selected_regions)):
      id1 = selected_regions[j]
      test_centers1 = ncs_group_obj.region_scattered_points_dict[id1]
      dist = get_closest_dist(test_centers, test_centers1)
      dd[id][id1] = dist
      if not id1 in dd: dd[id1] = {}
      dd[id1][id] = dist
  return dd

def get_dist_to_first_dict(ncs_group_obj = None,
     selected_regions = None,
     inter_region_dist_dict = None,
     target_scattered_points = None):

  """Get distance to region 0 ( or to target_scattered_points if supplied)
  for each region in selected_regions"""
  dist_to_first_dict = {}
  if target_scattered_points:
    start_region = 0
    for x in selected_regions:
      dist_to_first_dict[x] = get_closest_dist(
        ncs_group_obj.region_scattered_points_dict[x],
        target_scattered_points)
  else:
    start_region = 1
    x0 = selected_regions[0]
    dist_to_first_dict[x0] = 0
    for x in selected_regions[1:]:
      dist_to_first_dict[x] = inter_region_dist_dict[x0][x]
  changing = True
  while changing:
    changing = False
    for x in selected_regions[start_region:]:
      for y in selected_regions[start_region:]:
        if x == y: continue
        if dist_to_first_dict[y]<dist_to_first_dict[x] and \
            inter_region_dist_dict[x][y]<dist_to_first_dict[x]:
          dist_to_first_dict[x] = max(
            dist_to_first_dict[y], inter_region_dist_dict[x][y])
          changing = True
  return dist_to_first_dict

def radius_of_gyration_of_vector(xyz):
  """Return radius of gyration of flex.vec3_double object xyz"""
  return (xyz-xyz.mean()).rms_length()

def get_radius_of_gyration(ncs_group_obj = None,
    selected_regions = None):
  """Return radius of gyration of all points contained in selected regions"""
  centers = flex.vec3_double()
  for s in selected_regions:
    centers.append(ncs_group_obj.region_centroid_dict[s])
  centers = centers-centers.mean()
  return centers.rms_length()


def get_closest_neighbor_rms(ncs_group_obj = None, selected_regions = None,
    target_scattered_points = None, verbose = False, out = sys.stdout):

  """Return rms closest distance of each region center to
     lowest_numbered region, allowing sequential tracking taking max
     of inter-region distances
  """

  # XXX can't we save some of this for next time?

  inter_region_dist_dict = get_inter_region_dist_dict(ncs_group_obj = ncs_group_obj,
     selected_regions = selected_regions)
  if verbose:
    print("Inter-region distance dict:", file = out)
    keys = list(inter_region_dist_dict.keys())
    keys.sort()
    for key in keys:
      for key2 in inter_region_dist_dict[key].keys():
        print("%s  %s  : %.1f " %(key, key2, inter_region_dist_dict[key][key2]), file = out)

  dist_to_first_dict = get_dist_to_first_dict(ncs_group_obj = ncs_group_obj,
     selected_regions = selected_regions,
     inter_region_dist_dict = inter_region_dist_dict,
     target_scattered_points = target_scattered_points)

  if verbose:
    print("Distance-to-first dict:", file = out)
    keys = list(dist_to_first_dict.keys())
    keys.sort()
    for key in keys: print("\n %s:  %.1f " %(key, dist_to_first_dict[key]), file = out)

  if target_scattered_points:
    start_region = 0 # we are getting dist to target_scattered_points
  else:
    start_region = 1 # we are getting dist to region 0

  rms = 0.
  rms_n = 0.
  for x in selected_regions[start_region:]:
    dist = dist_to_first_dict[x]
    rms+= dist**2
    rms_n+= 1.
  if rms_n>1:
    rms/= rms_n
  rms = rms**0.5
  return rms


def get_rms(selected_regions = None,
    region_centroid_dict = None):
  """Return rms distance of each region center from average of all others"""
  rms = 0.
  rms_n = 0.
  for x in selected_regions:
    other_regions = remove_one_item(selected_regions, item_to_remove = x)
    current_center = get_average_center(other_regions,
       region_centroid_dict = region_centroid_dict)
    test_center = region_centroid_dict[x]
    dist = get_dist(current_center, test_center)
    rms+= dist**2
    rms_n+= 1.
  if rms_n>1:
    rms/= rms_n
  return rms**0.5

def single_list(list_of_lists):
  """Convert list_of_lists to a single list"""
  single = []
  for x in list_of_lists:
    if type(x) == type([1, 2, 3]):
      single+= single_list(x)
    else:
      single.append(x)
  return single

def get_closest_dist(test_center, target_centers):
  """Get the closest distance between test_center and target_centers.
  Converts to flex.vece_double() and uses min_distance_between_any_pair() """

  # make sure we have target_centers = vec3_double and not a list,
  #  and vec3_double or tuple for test_center

  if type(test_center) == type([1, 2, 3]):
    test_center = flex.vec3_double(test_center)
  if type(target_centers) == type([1, 2, 3]):
    target_centers = flex.vec3_double(target_centers)
  if test_center.size()<1 or target_centers.size()<1: return None
  closest_dist = test_center.min_distance_between_any_pair(target_centers)
  return closest_dist

def region_lists_have_ncs_overlap(set1, set2, ncs_group_obj = None, cutoff = 0):
  """Return True if region_lists (set1, set2) share members through
    shared_group_dict"""
  for id1 in set1:
    for id2 in set2:
      if id2 in ncs_group_obj.shared_group_dict.get(id1, []):
        return True
  return False

def get_effective_radius(ncs_group_obj = None,
    target_scattered_points = None,
    weight_rad_gyr = None,
    selected_regions = None):
  """Get the effective radius of gyration of all points contained
   in selected_regions"""
  sr = deepcopy(selected_regions)
  sr.sort()
  rad_gyr = get_radius_of_gyration(ncs_group_obj = ncs_group_obj,
     selected_regions = sr)
  rms = get_closest_neighbor_rms(ncs_group_obj = ncs_group_obj,
    target_scattered_points = target_scattered_points,
    selected_regions = sr)
  max_cell_dim = 0.
  if ncs_group_obj.max_cell_dim and ncs_group_obj.max_cell_dim > 1.0:
    wrg = weight_rad_gyr*(300/ncs_group_obj.max_cell_dim)  # have a consistent scale
  else:
    wrg = weight_rad_gyr
  effective_radius = (rms+wrg*rad_gyr)/(1.+wrg)
  return effective_radius

def add_neighbors(params,
      selected_regions = None,
      max_length_of_group = None,
      target_scattered_points = None,
      tracking_data = None,
      equiv_dict_ncs_copy = None,
      ncs_group_obj = None, out = sys.stdout):

  """Add neighboring regions on to selected_regions.
  Same rules as select_from_seed"""

  selected_regions = single_list(deepcopy(selected_regions))

  added_regions = []
  start_dist = get_effective_radius(ncs_group_obj = ncs_group_obj,
        target_scattered_points = target_scattered_points,
        weight_rad_gyr = params.segmentation.weight_rad_gyr,
        selected_regions = selected_regions)
  delta_dist = params.segmentation.add_neighbors_dist
  max_dist = start_dist+delta_dist

  starting_selected_regions = deepcopy(selected_regions)

  for x in selected_regions:  # delete, add in alternatives one at a time and
    #  keep all the ok ones
    ncs_groups_to_use = get_ncs_related_regions(
      ncs_group_obj = ncs_group_obj,
      selected_regions = [x],
      include_self = False)

    for x in ncs_groups_to_use: # try adding from each group
      if x in selected_regions+added_regions:
        continue
      ncs_group = [[x]]
      current_scattered_points_list = get_scattered_points_list(selected_regions,
        region_scattered_points_dict = ncs_group_obj.region_scattered_points_dict)

      for ncs_set in ncs_group: # pick the best ncs_set from this group
        if has_intersection(ncs_group_obj.bad_region_list, ncs_set):
          continue

        dist = get_effective_radius(ncs_group_obj = ncs_group_obj,
          target_scattered_points = target_scattered_points,
          weight_rad_gyr = params.segmentation.weight_rad_gyr,
          selected_regions = selected_regions+ncs_set)

        if dist <=  max_dist:
          added_regions.append(x)

  selected_regions = selected_regions+added_regions
  dist = get_effective_radius(ncs_group_obj = ncs_group_obj,
          target_scattered_points = target_scattered_points,
          weight_rad_gyr = params.segmentation.weight_rad_gyr,
          selected_regions = selected_regions)

  # Identify all the NCS operators required to map final to starting
  # equiv_dict_ncs_copy[id][id1] = ncs_copy of id that corresponds to id1
  ncs_group = ncs_group_obj.ncs_obj.ncs_groups()[0]
  identity_op = ncs_group.identity_op_id()
  ncs_ops_used = [identity_op]

  for id in selected_regions:
    related_regions = get_ncs_related_regions(
      ncs_group_obj = ncs_group_obj,
      selected_regions = [id],
      include_self = False)
    for id1 in selected_regions:
      if not id1 in related_regions: continue
      ncs_copy1 = equiv_dict_ncs_copy.get(id, {}).get(id1, None)
      ncs_copy2 = equiv_dict_ncs_copy.get(id1, {}).get(id, None)
      for a in [ncs_copy1, ncs_copy2]:
        if a is not None and not a in ncs_ops_used:
            ncs_ops_used.append(a)
  selected_regions.sort()
  ncs_ops_used.sort()
  for x in selected_regions:
    print("GROUP ", x, ":", ncs_group_obj.shared_group_dict.get(x, []), file = out)

  return selected_regions, dist, ncs_ops_used

def select_from_seed(params,
      starting_regions,
      target_scattered_points = None,
      max_length_of_group = None,
      ncs_groups_to_use = None,
      tracking_data = None,
      ncs_group_obj = None):
  """
  Select additional regions to group with starting_regions
  Do not allow any region in ncs_group_obj.bad_region_list
  also do not allow any region that is in an ncs-related group to any region
  already used.  Use ncs_group_obj.equiv_dict to identify these.
  Return selected_regions, dist"""

  selected_regions = single_list(deepcopy(starting_regions))
  if not ncs_groups_to_use:
    ncs_groups_to_use = ncs_group_obj.ncs_group_list

  for ncs_group in ncs_groups_to_use: # try adding from each group
    if max_length_of_group is not None and \
       len(selected_regions)>= max_length_of_group:
      break
    best_ncs_set = None
    best_dist = None
    if has_intersection(ncs_group, selected_regions):
      continue
    current_scattered_points_list = get_scattered_points_list(selected_regions,
       region_scattered_points_dict = ncs_group_obj.region_scattered_points_dict)
    if target_scattered_points:
      current_scattered_points_list.extend(target_scattered_points)

    for ncs_set in ncs_group: # pick the best ncs_set from this group
      if has_intersection(ncs_group_obj.bad_region_list, ncs_set): continue

      # does any ncs copy of anything in selected_regions actually overlap
      #  with any member of ncs_set... might be efficient to delete the entire
      #   ncs_group if any ncs_set overlaps, but could lose some.
      if region_lists_have_ncs_overlap(ncs_set, selected_regions,
          ncs_group_obj = ncs_group_obj):
        continue

      dist = get_effective_radius(ncs_group_obj = ncs_group_obj,
        target_scattered_points = target_scattered_points,
        weight_rad_gyr = params.segmentation.weight_rad_gyr,
        selected_regions = selected_regions+ncs_set)

      if best_dist is None or dist<best_dist:
        best_dist = dist
        best_ncs_set = ncs_set
    if best_ncs_set is not None:
      selected_regions+= best_ncs_set

  dist = get_effective_radius(ncs_group_obj = ncs_group_obj,
    target_scattered_points = target_scattered_points,
    weight_rad_gyr = params.segmentation.weight_rad_gyr,
    selected_regions = selected_regions)

  return selected_regions, dist

def remove_one_item(input_list, item_to_remove = None):
  """Remove one item from list. Deprecated.
  Use instead
    if input_list.find(item_to_remove) > -1:
      input_list.remove(item_to_remove)
  """
  new_list = []
  for item in input_list:
    if item !=  item_to_remove:
      new_list.append(item)
  return new_list


def get_ncs_related_regions_specific_list(
    ncs_group_obj = None,
    target_regions = None,
    include_self = False):
  """Return all regions ncs-related to target_regions"""
  all_regions = []
  for target_region in target_regions:
    all_regions+= get_ncs_related_regions_specific_target(
      ncs_group_obj = ncs_group_obj,
      target_region = target_region,
      other_regions = remove_one_item(
         target_regions, item_to_remove = target_region),
      include_self = include_self)
  return all_regions

def get_ncs_related_regions_specific_target(
          ncs_group_obj = None,
          target_region = None,
          other_regions = None,
          include_self = False):
  """Similar to get_ncs_related_regions, but find just one  ncs group that
  contains x but does not contain any member of other_regions
  """
  for ncs_group in ncs_group_obj.ncs_group_list: # might this be the group
    ids_in_group = single_list(ncs_group)
    if not target_region in ids_in_group: continue # does not contain target
    contains_others = False
    for other_id in other_regions:
      if other_id in ids_in_group:
        contains_other = True
        break# contains other members
    if not contains_others:
      # this is the group
      if include_self:
        return ids_in_group
      else:
        return remove_one_item(ids_in_group, item_to_remove = target_region)
  return []


def get_ncs_related_regions(
    ncs_group_obj = None,
    selected_regions = None,
    include_self = False):
  """Returns a simple list of region ids NCS-related to selected_regions.
  if include_self then include selected regions and all ncs-related
  otherwise do not include selected regions or anything that might
  overlap with them"""

  ncs_related_regions = []
  if include_self:
    for id in selected_regions:
      if not id in ncs_related_regions:
        ncs_related_regions.append(id)
      for ncs_group in ncs_group_obj.ncs_group_list:
        ids_in_group = single_list(ncs_group)
        if id in ids_in_group:  # this group contains this selected id
          for i in ids_in_group:
            if not i in ncs_related_regions:
              ncs_related_regions.append(i)

  else:
    for id in selected_regions:
      found = False
      for ncs_group in ncs_group_obj.ncs_group_list:
        ids_in_group = single_list(ncs_group)
        if id in ids_in_group:  # this group contains this selected id
          found = True
          for i in ids_in_group:
            if (not i == id) and (not i in selected_regions) and \
               (not i in ncs_related_regions):
              ncs_related_regions.append(i)
          break # don't look at any more ncs groups

  return ncs_related_regions

def all_elements_are_length_one(list_of_elements):
  """Return True if all elements have a length of one"""
  for x in list_of_elements:
    if type(x) == type([1, 2, 3]):
      if len(x)!= 1: return False
  return True

def as_list_of_lists(ll):
  """Return list ll as list-of-lists"""
  new_list = []
  for x in ll:
    new_list.append([x])
  return new_list

def select_regions_in_au(params,
     ncs_group_obj = None,
     target_scattered_points = None,
     unique_expected_regions = None,
     equiv_dict_ncs_copy = None,
     tracking_data = None,
     out = sys.stdout):
  """Choose one region or set of regions from each ncs_group
  up to about unique_expected_regions.
  Optimize closeness of centers...
  If target scattered_points is supplied, include them as allowed target"""

  if not ncs_group_obj.ncs_group_list:
    return ncs_group_obj, []

  max_length_of_group = max(1, unique_expected_regions*
     params.segmentation.max_per_au_ratio)
  print("Maximum length of group: %d" %(max_length_of_group), file = out)

  if all_elements_are_length_one(ncs_group_obj.ncs_group_list):
    # This is where there is no ncs. Basically skipping everything
    best_selected_regions = single_list(ncs_group_obj.ncs_group_list)
    best_rms = None
  else:
    #--------------  Find initial set of regions --------------------
    # Seed with members of the first NCS group or with the target points
    #  and find the member of each NCS group that is closest

    if target_scattered_points:
      starting_regions = [None]
    else:
      starting_regions = ncs_group_obj.ncs_group_list[0]

    best_selected_regions = None
    best_rms = None
    ok_seeds_examined = 0
    for starting_region in starting_regions: # NOTE starting_region is a list
      if not starting_region and not target_scattered_points:continue
      if ok_seeds_examined >=  params.segmentation.seeds_to_try:
        break # don't bother to keep trying
      if starting_region and starting_region in ncs_group_obj.bad_region_list:
        continue # do not use
      if starting_region: # NOTE: starting_region is a list itself
        starting_region_list = [starting_region]
      else:
        starting_region_list = []
      selected_regions, rms = select_from_seed(params,
        starting_region_list,
        target_scattered_points = target_scattered_points,
        max_length_of_group = max_length_of_group,
        tracking_data = tracking_data,
        ncs_group_obj = ncs_group_obj)
      if not selected_regions:
        continue
      ok_seeds_examined+= 1
      if best_rms is None or rms<best_rms:
        best_rms = rms
        best_selected_regions = selected_regions
        print("New best selected: rms: %7.1f: %s " %(
           rms, str(selected_regions)), file = out)

    if best_rms is not None:
      print("Best selected so far: rms: %7.1f: %s " %(
            best_rms, str(best_selected_regions)), file = out)

    if not best_selected_regions:
      print("\nNo NCS regions found ...", file = out)
      return ncs_group_obj, []

    # Now we have a first version of best_rms, best_selected_regions

    #--------------  END Find initial set of regions --------------------


    #--------------  Optimize choice of regions -------------------------
    max_tries = 10
    improving = True
    itry = 0
    while improving and itry<max_tries:
      itry+= 1
      improving = False
      previous_selected_regions = deepcopy(best_selected_regions)
      previous_selected_regions.sort()
      print("\nTry %d for optimizing regions" %(itry), file = out)
      # Now see if replacing any regions with alternatives would improve it
      for x in previous_selected_regions:
        starting_regions = remove_one_item(previous_selected_regions,
          item_to_remove = x)
        # identify ncs_related regions to x, but not to other members of
        #  selected_regions
        ncs_related_regions = get_ncs_related_regions_specific_list(
          ncs_group_obj = ncs_group_obj,
          include_self = True,
          target_regions = [x])
        if not ncs_related_regions: continue
        ncs_groups_to_use = [as_list_of_lists(ncs_related_regions)]
        new_selected_regions, rms = select_from_seed(params, starting_regions,
          target_scattered_points = target_scattered_points,
          max_length_of_group = max_length_of_group,
          tracking_data = tracking_data,
          ncs_groups_to_use = ncs_groups_to_use,
          ncs_group_obj = ncs_group_obj)

        if not new_selected_regions: continue
        if best_rms is None or rms<best_rms:
          best_selected_regions = new_selected_regions
          best_selected_regions.sort()
          best_rms = rms
          improving = True
      print("Optimized best selected: rms: %7.1f: %s " %(
          best_rms, str(best_selected_regions)), file = out)

      # Done with this try

  selected_regions = best_selected_regions
  selected_regions.sort()
  available_selected_regions = len(selected_regions)
  print("\nAvailable selected regions: %s ..." %(available_selected_regions), file = out)
  if tracking_data:
    tracking_data.available_selected_regions = available_selected_regions

  if params.map_modification.regions_to_keep is not None:
    if params.map_modification.regions_to_keep <= 0:
       # keep just region abs(regions_to_keep)
       ii = min(len(selected_regions)-1, abs(params.map_modification.regions_to_keep))
       selected_regions = selected_regions[ii:ii+1]
    else: # usual
      selected_regions = selected_regions[:params.map_modification.regions_to_keep]

  rms = get_closest_neighbor_rms(ncs_group_obj = ncs_group_obj,
    selected_regions = selected_regions, verbose = False, out = out)

  if params.segmentation.add_neighbors and \
       ncs_group_obj.ncs_obj.max_operators()>1:
    print("\nAdding neighbor groups...", file = out)
    selected_regions, rms, ncs_ops_used = add_neighbors(params,
          selected_regions = selected_regions,
          max_length_of_group = max_length_of_group,
          target_scattered_points = target_scattered_points,
          equiv_dict_ncs_copy = equiv_dict_ncs_copy,
          tracking_data = tracking_data,
          ncs_group_obj = ncs_group_obj, out = out)
  else:
    ncs_ops_used = None

  print("\nFinal selected regions with rms of %6.2f: " %(rms), end = ' ', file = out)
  for x in selected_regions:
    print(x, end = ' ', file = out)
  if ncs_ops_used:
    print("\nNCS operators used: ", end = ' ', file = out)
    for op in ncs_ops_used:  print(op, end = ' ', file = out)
    print(file = out)
  # Save an ncs object containing just the ncs_ops_used
  ncs_group_obj.set_ncs_ops_used(ncs_ops_used)

  # Identify scattered points for all selected regions:

  scattered_points = get_scattered_points_list(selected_regions,
     region_scattered_points_dict = ncs_group_obj.region_scattered_points_dict)

  # Identify ncs-related regions for all the selected regions
  self_and_ncs_related_regions = get_ncs_related_regions(
    ncs_group_obj = ncs_group_obj,
    selected_regions = selected_regions,
    include_self = True)

  ncs_related_regions = get_ncs_related_regions(
    ncs_group_obj = ncs_group_obj,
    selected_regions = selected_regions,
    include_self = False)

  print("NCS-related regions (not used): %d " %(len(ncs_related_regions)), file = out)
  ncs_group_obj.set_selected_regions(selected_regions)
  ncs_group_obj.set_self_and_ncs_related_regions(self_and_ncs_related_regions)
  ncs_group_obj.set_ncs_related_regions(ncs_related_regions)

  return ncs_group_obj, scattered_points

def get_bool_mask_as_int(ncs_group_obj = None, mask_as_int = None, mask_as_bool = None):
  """Convert bool mask into integer mask"""
  if mask_as_int:
    mask_as_int = mask_as_int.deep_copy()
  else:
    mask_as_int = ncs_group_obj.edited_mask.deep_copy()
  s = (mask_as_bool == True)
  mask_as_int = mask_as_int.set_selected(s, 1)
  mask_as_int = mask_as_int.set_selected(~s, 0)
  return mask_as_int

def get_bool_mask_of_regions(ncs_group_obj = None, region_list = None,
    expand_size = None):
  """Return a bool mask marking location of regions in region_list"""
  s = (ncs_group_obj.edited_mask  ==  -1)
  if region_list is None: region_list = []
  for id in region_list:

    if not expand_size:
      s |=  (ncs_group_obj.edited_mask == id)  # just take this region

    else:  # expand the size of the regions...use expand_mask which operates
         # on the original id numbers and uses the co
      bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand = ncs_group_obj.original_id_from_id[id],
        expand_size = expand_size)
      s |=  (bool_region_mask ==  True)

  bool_mask = ncs_group_obj.co.expand_mask(id_to_expand = 1, expand_size = 1) # just to get bool mask
  bool_mask = bool_mask.set_selected(s, True)
  bool_mask = bool_mask.set_selected(~s, False)

  return bool_mask

def create_remaining_mask_and_map(params,
    ncs_group_obj = None,
    map_data = None,
    crystal_symmetry = None,
    out = sys.stdout):

  """Return remainder map after
  removing ncs_group_obj.selected_regions"""

  if not ncs_group_obj.selected_regions:
    print("No regions selected", file = out)
    return map_data

  # create new remaining_map containing everything except the part that
  # has been interpreted (and all points in interpreted NCS-related copies)

  bool_all_used = get_bool_mask_of_regions(ncs_group_obj = ncs_group_obj,
   region_list = ncs_group_obj.selected_regions+
       ncs_group_obj.self_and_ncs_related_regions,
   expand_size = params.segmentation.expand_size)
  map_data_remaining = map_data.deep_copy()
  s = (bool_all_used == True)

  map_data_remaining = map_data_remaining.set_selected(s,
    params.segmentation.value_outside_mask)
  return map_data_remaining

def get_lower(lower_bounds, lower):
  """Update lower_bounds by min(lower_bounds, lower)"""
  new_lower = []
  for i in range(3):
    if lower_bounds[i] is None:
      new_lower.append(lower[i])
    elif lower[i] is None:
      new_lower.append(lower_bounds[i])
    else:
      new_lower.append(min(lower_bounds[i], lower[i]))
  return new_lower

def get_upper(upper_bounds, upper):
  """Update upper_bounds by max(upper_bounds, upper)"""
  new_upper = []
  for i in range(3):
    if upper_bounds[i] is None:
      new_upper.append(upper[i])
    elif upper[i] is None:
      new_upper.append(upper_bounds[i])
    else:
      new_upper.append(max(upper_bounds[i], upper[i]))
  return new_upper

def get_bounds(ncs_group_obj = None, id = None):
  """Get bounds (grid units) for ncs_group id"""
  orig_id = ncs_group_obj.original_id_from_id[id]
  lower = ncs_group_obj.min_b[orig_id]
  upper = ncs_group_obj.max_b[orig_id]
  return lower, upper

def get_selected_and_related_regions(params,
    ncs_group_obj = None):
  """Identify all points in the targeted regions"""
  bool_selected_regions = get_bool_mask_of_regions(
       ncs_group_obj = ncs_group_obj,
     region_list = ncs_group_obj.selected_regions,
     expand_size = params.segmentation.expand_size+\
      params.segmentation.mask_additional_expand_size)
  # and all points in NCS-related copies (to be excluded)
  if params.segmentation.exclude_points_in_ncs_copies and (
     not params.segmentation.add_neighbors):
    bool_ncs_related_mask = get_bool_mask_of_regions(ncs_group_obj = ncs_group_obj,
       region_list = ncs_group_obj.ncs_related_regions)
     # NOTE: using ncs_related_regions here NOT self_and_ncs_related_regions
  else:
    bool_ncs_related_mask = None

  lower_bounds = [None, None, None]
  upper_bounds = [None, None, None]
  if ncs_group_obj.selected_regions:
    for id in ncs_group_obj.selected_regions:
      lower, upper = get_bounds(
        ncs_group_obj = ncs_group_obj, id = id)
      lower_bounds = get_lower(lower_bounds, lower)
      upper_bounds = get_upper(upper_bounds, upper)

  return bool_selected_regions, bool_ncs_related_mask, lower_bounds, upper_bounds

def adjust_bounds(params,
   lower_bounds, upper_bounds, map_data = None, out = sys.stdout):
  """Adjust lower and upper bounds to add params.output_files.box_buffer"""
  # range is lower_bounds to upper_bounds
  lower_bounds = list(lower_bounds)
  upper_bounds = list(upper_bounds)
  if params is None or params.output_files.box_buffer is None:
     box_buffer = 0
  else:
     box_buffer = int(0.5+params.output_files.box_buffer)
  for i in range(3):
    if lower_bounds[i] is None: lower_bounds[i] = 0
    if upper_bounds[i] is None: upper_bounds[i] = 0
    lower_bounds[i]-= box_buffer
    lower_bounds[i] = max(0, lower_bounds[i])
    upper_bounds[i]+= box_buffer
    upper_bounds[i] = min(map_data.all()[i]-1, upper_bounds[i])


  """
  print >>out, "\nRange:  X:(%6d, %6d)    Y:(%6d, %6d)    Z:(%6d, %6d)" %(
     lower_bounds[0], upper_bounds[0],
     lower_bounds[1], upper_bounds[1],
     lower_bounds[2], upper_bounds[2])
  """

  return lower_bounds, upper_bounds

def write_region_maps(params,
    ncs_group_obj = None,
    map_data = None,
    tracking_data = None,
    remainder_ncs_group_obj = None,
    regions_to_skip = None,
    out = sys.stdout):
  """Write maps marking each region (segmented region)"""
  remainder_regions_written = []
  map_files_written = []
  if not ncs_group_obj:
    return map_files_written, remainder_regions_written

  if not ncs_group_obj.selected_regions:
    return map_files_written, remainder_regions_written

  for id in ncs_group_obj.selected_regions:


    if regions_to_skip and id in regions_to_skip:
      print("Skipping remainder region %d (already written out)" %(id), file = out)
      continue
    print("Writing region %d" %(id), end = ' ', file = out)

    # dummy atoms representing this region
    sites = ncs_group_obj.region_scattered_points_dict[id]

    bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand = ncs_group_obj.original_id_from_id[id],
        expand_size = params.segmentation.expand_size)

    s = (bool_region_mask == True)

    lower_bounds, upper_bounds = get_bounds(ncs_group_obj = ncs_group_obj, id = id)

    if remainder_ncs_group_obj:
      for remainder_id in remainder_ncs_group_obj.remainder_id_dict.keys():
        if remainder_ncs_group_obj.remainder_id_dict[remainder_id] == id:
          remainder_regions_written.append(remainder_id)

          sites.extend(
            remainder_ncs_group_obj.region_scattered_points_dict[remainder_id])

          print("(including remainder region %d)" %(remainder_id), end = ' ', file = out)
          remainder_bool_region_mask = remainder_ncs_group_obj.co.expand_mask(
           id_to_expand = remainder_ncs_group_obj.original_id_from_id[remainder_id],
           expand_size = params.segmentation.expand_size)
          s|=  (remainder_bool_region_mask == True)
          lower, upper = get_bounds(
            ncs_group_obj = remainder_ncs_group_obj, id = remainder_id)
          lower_bounds = get_lower(lower_bounds, lower)
          upper_bounds = get_upper(upper_bounds, upper)


    region_mask = map_data.deep_copy()
    region_mask = region_mask.set_selected(s, 1)
    region_mask = region_mask.set_selected(~s, 0)
    local_map_data = map_data.deep_copy()
    local_map_data = local_map_data * region_mask.as_double()

    # Now cut down the map to the size we want
    lower_bounds, upper_bounds = adjust_bounds(params, lower_bounds, upper_bounds,
      map_data = map_data, out = out)
    box_map, box_crystal_symmetry, \
      dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
      map_data = local_map_data, \
       crystal_symmetry = tracking_data.crystal_symmetry,
       min_point = lower_bounds, max_point = upper_bounds, out = out)

    if remainder_ncs_group_obj:
      text = ""
    else:
      text = "_r"
    base_file = 'map%s_%d.ccp4' %(text, id)
    base_pdb_file = 'atoms%s_%d.pdb' %(text, id) # PDB OK just atoms
    if tracking_data.params.output_files.output_directory:
      if not os.path.isdir(tracking_data.params.output_files.output_directory):
        os.mkdir(tracking_data.params.output_files.output_directory)
      file_name = os.path.join(
        tracking_data.params.output_files.output_directory, base_file)
      pdb_file_name = os.path.join(
        tracking_data.params.output_files.output_directory, base_pdb_file)
    else:
      file_name = base_file
      pdb_file_name = base_pdb_file
    write_ccp4_map(box_crystal_symmetry, file_name, box_map)
    print("to %s" %(file_name), file = out)
    map_files_written.append(file_name)

    tracking_data.add_output_region_map_info(
      file_name = file_name,
      crystal_symmetry = box_crystal_symmetry,
      origin = box_map.origin(),
      all = box_map.all(),
      map_id = base_file)

    print("Atoms representation written to %s" %(pdb_file_name), file = out)
    write_atoms(tracking_data = tracking_data, sites = sites, file_name = pdb_file_name,
       out = out)
    tracking_data.add_output_region_pdb_info(
      file_name = pdb_file_name)

  return map_files_written, remainder_regions_written

def get_bounds_from_sites(sites_cart = None, map_data = None,
    unit_cell = None):
  """Calculate bounds (grid units) containing sites_cart"""
  lower_bounds = [None, None, None]
  upper_bounds = [None, None, None]
  sites_frac = unit_cell.fractionalize(sites_cart)
  nx, ny, nz = map_data.all()
  for x_frac in sites_frac:
    x = [
      int(0.5+nx*x_frac[0]),
      int(0.5+ny*x_frac[1]),
      int(0.5+nz*x_frac[2])]

    if lower_bounds[0] is None or x[0]<lower_bounds[0]: lower_bounds[0] = x[0]
    if lower_bounds[1] is None or x[1]<lower_bounds[1]: lower_bounds[1] = x[1]
    if lower_bounds[2] is None or x[2]<lower_bounds[2]: lower_bounds[2] = x[2]

    if upper_bounds[0] is None or x[0]>upper_bounds[0]: upper_bounds[0] = x[0]
    if upper_bounds[1] is None or x[1]>upper_bounds[1]: upper_bounds[1] = x[1]
    if upper_bounds[2] is None or x[2]>upper_bounds[2]: upper_bounds[2] = x[2]
  return lower_bounds, upper_bounds

def write_output_files(params,
    tracking_data = None,
    map_data = None,
    half_map_data_list = None,
    ncs_group_obj = None,
    remainder_ncs_group_obj = None,
    pdb_hierarchy = None,
    removed_ncs = None,
    out = sys.stdout):

  """Write output files based on tracking_data."""
  half_map_data_list_au = []
  if not half_map_data_list: half_map_data_list = []

  if params.output_files.au_output_file_stem:
    au_mask_output_file = os.path.join(tracking_data.params.output_files.output_directory, params.output_files.au_output_file_stem+"_mask.ccp4")
    au_map_output_file = os.path.join(tracking_data.params.output_files.output_directory, params.output_files.au_output_file_stem+"_map.ccp4")
    au_atom_output_file = os.path.join(tracking_data.params.output_files.output_directory, params.output_files.au_output_file_stem+"_atoms.pdb")
  else:
    au_mask_output_file = None
    au_map_output_file = None
    au_atom_output_file = None

  # Write out pdb file with dummy atoms for the AU to au_atom_output_file
  if au_atom_output_file and params.output_files.write_output_maps:
    sites = flex.vec3_double()
    for id in ncs_group_obj.selected_regions:
      sites.extend(ncs_group_obj.region_scattered_points_dict[id])
    if remainder_ncs_group_obj:
      for id in remainder_ncs_group_obj.selected_regions:
        sites.extend(remainder_ncs_group_obj.region_scattered_points_dict[id])
    write_atoms(tracking_data = tracking_data, sites = sites,
      file_name = au_atom_output_file, out = out)
    tracking_data.set_output_ncs_au_pdb_info(file_name = au_atom_output_file)


  # Write out mask and map representing one NCS copy and none of
  #   other NCS copies.  Expand the mask to include neighboring points (but
  #   not those explicitly in other NCS copies

  if params.map_modification.soft_mask and params.control.save_box_map_ncs_au:
    mask_expand_size = estimate_expand_size(
       crystal_symmetry = tracking_data.crystal_symmetry,
       map_data = map_data,
       expand_target = tracking_data.params.segmentation.mask_expand_ratio*\
          tracking_data.params.crystal_info.resolution,
          out = out)
    params.segmentation.mask_additional_expand_size = max(mask_expand_size,
      params.segmentation.mask_additional_expand_size, )

  bool_selected_regions, bool_ncs_related_mask, lower_bounds, upper_bounds = \
     get_selected_and_related_regions(
      params, ncs_group_obj = ncs_group_obj)
  if bool_ncs_related_mask is not None:
    s_ncs_related =  (bool_ncs_related_mask == True)
  else:
    s_ncs_related =  None

  # Add in remainder regions if present
  if remainder_ncs_group_obj:
    bool_remainder_selected_regions, bool_remainder_ncs_related_mask, \
      remainder_lower_bounds, remainder_upper_bounds = \
       get_selected_and_related_regions(
       params, ncs_group_obj = remainder_ncs_group_obj)

    lower_bounds = get_lower(lower_bounds, remainder_lower_bounds)
    upper_bounds = get_upper(upper_bounds, remainder_upper_bounds)

    s_remainder_au =  (bool_remainder_selected_regions == True)
    bool_selected_regions = bool_selected_regions.set_selected(
       s_remainder_au, True)
    if s_ncs_related is not None and \
         bool_remainder_ncs_related_mask is not None:
      s_ncs_related |=   (bool_remainder_ncs_related_mask == True)

  # Now create NCS mask by eliminating all points in target (expanded) in
  #   NCS-related copies
  if s_ncs_related is not None:
    bool_selected_regions = bool_selected_regions.set_selected(
       s_ncs_related, False)

  if tracking_data.params.map_modification.regions_to_keep is None:
    # Identify full (possibly expanded) ncs au starting with what we have
    au_mask = get_one_au(tracking_data = tracking_data,
       starting_mask = bool_selected_regions,
      removed_ncs = removed_ncs,
      ncs_obj = ncs_group_obj.ncs_obj, map_data = map_data, out = out)
    print("\nExpanding NCS AU if necessary...", file = out)
    print("Size of AU mask: %s  Current size of AU: %s" %(
      au_mask.count(True), bool_selected_regions.count(True)), file = out)
    bool_selected_regions = (bool_selected_regions | au_mask)
    print("New size of AU mask: %s" %(bool_selected_regions.count(True)), file = out)

  sites_cart = get_marked_points_cart(mask_data = bool_selected_regions,
     unit_cell = ncs_group_obj.crystal_symmetry.unit_cell(),
     every_nth_point = tracking_data.params.segmentation.grid_spacing_for_au,
     boundary_radius = tracking_data.params.segmentation.radius)
  sites_lower_bounds, sites_upper_bounds = get_bounds_from_sites(
      unit_cell = ncs_group_obj.crystal_symmetry.unit_cell(),
      sites_cart = sites_cart, map_data = map_data)
  print("Original bounds: %5s  %5s  %5s  to %5s  %5s  %5s" %(
    tuple(lower_bounds+upper_bounds)), file = out)
  lower_bounds = get_lower(lower_bounds, sites_lower_bounds)
  upper_bounds = get_upper(upper_bounds, sites_upper_bounds)
  print("Updated bounds:  %5s  %5s  %5s  to %5s  %5s  %5s" %(
    tuple(lower_bounds+upper_bounds)), file = out)

  lower_bounds, upper_bounds = adjust_bounds(params, lower_bounds, upper_bounds,
    map_data = map_data, out = out)
  box_ncs_au = params.segmentation.box_ncs_au
  if (not box_ncs_au):
    print("Using entire input map (box_ncs_au = False)", file = out)
    lower_bounds = map_data.origin()
    upper_bounds = tuple(matrix.col(map_data.all())+
        matrix.col(map_data.origin())-matrix.col((1, 1, 1)))


  print("\nMaking two types of maps for AU of NCS mask and map with "+\
      "buffer of %d grid units \nin each direction around AU" %(
      params.output_files.box_buffer), file = out)
  if params.output_files.write_output_maps:
   print("Both types of maps have the same origin and overlay on %s" %(
   os.path.join(tracking_data.params.output_files.output_directory,
     params.output_files.shifted_map_file)), file = out)


   print("\nThe standard maps (%s, %s) have the \noriginal cell dimensions." %(
   os.path.join(tracking_data.params.output_files.output_directory, au_mask_output_file),
   os.path.join(tracking_data.params.output_files.output_directory, au_map_output_file))+\
   "\nThese maps show only the unique (NCS AU) part of the map.", file = out)

   print("\nThe cut out box_maps (%s, %s) have \nsmaller cell dimensions." %(
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_mask_file),
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_map_file), ) +\
   "\nThese maps also show only the unique part of the map and have this"+\
   "\nunique part cut out.\n", file = out)


  # Write out NCS AU with shifted origin but initial crystal_symmetry
  # Mask
  mask_data_ncs_au = get_bool_mask_as_int(
     ncs_group_obj = ncs_group_obj, mask_as_bool = bool_selected_regions)

  if au_mask_output_file and params.output_files.write_output_maps:
    # Write out the mask (as int)
    write_ccp4_map(tracking_data.crystal_symmetry,
      au_mask_output_file, mask_data_ncs_au)
    print("Output NCS AU mask:  %s" %(au_mask_output_file), file = out)
    tracking_data.set_output_ncs_au_mask_info(
      file_name = au_mask_output_file,
      crystal_symmetry = tracking_data.crystal_symmetry,
      origin = mask_data_ncs_au.origin(),
      all = mask_data_ncs_au.all())

  # Map
  map_data_ncs_au = map_data.deep_copy()
  s = (bool_selected_regions == True)
  mask = map_data.deep_copy()
  mask = mask.set_selected(s, 1)
  mask = mask.set_selected(~s, 0)
  if params.map_modification.soft_mask:
    # buffer and smooth the mask
    map_data_ncs_au, smoothed_mask_data = apply_soft_mask(map_data = map_data_ncs_au,
      mask_data = mask.as_double(),
      rad_smooth = tracking_data.params.crystal_info.resolution,
      crystal_symmetry = tracking_data.crystal_symmetry,
      out = out)
    half_map_data_list_au = []
    for hm in half_map_data_list:  # apply mask to half maps
      hm_data_ncs_au, hm_smoothed_mask_data = apply_soft_mask(
        map_data = hm.deep_copy().as_double(),
        mask_data = mask.as_double(),
        rad_smooth = tracking_data.params.crystal_info.resolution,
        crystal_symmetry = tracking_data.crystal_symmetry,
        out = out)
      half_map_data_list_au.append(hm_data_ncs_au)


  elif (box_ncs_au): # usual.  If box_ncs_au is False, do not mask

    map_data_ncs_au = map_data_ncs_au*mask

    one_d = map_data_ncs_au.as_1d()
    n_zero = mask.count(0)
    n_tot = mask.size()
    mean_in_box = one_d.min_max_mean().mean*n_tot/(n_tot-n_zero)
    map_data_ncs_au = map_data_ncs_au+(1-mask)*mean_in_box
    half_map_data_list_au = []
    for hm in half_map_data_list:  # apply mask to half maps
      one_d = hm.as_1d()
      mean_in_box = one_d.min_max_mean().mean*n_tot/(n_tot-n_zero)
      hm_data_ncs_au = hm+(1-mask)*mean_in_box
      half_map_data_list_au.append(hm_data_ncs_au)

    del one_d, mask

  if au_map_output_file and params.output_files.write_output_maps:
    # Write out the NCS au of density
    write_ccp4_map(tracking_data.crystal_symmetry, au_map_output_file,
      map_data_ncs_au)
    print("Output NCS AU map:  %s" %(au_map_output_file), file = out)
    tracking_data.set_output_ncs_au_map_info(
      file_name = au_map_output_file,
      crystal_symmetry = tracking_data.crystal_symmetry,
      origin = map_data_ncs_au.origin(),
      all = map_data_ncs_au.all())

  # Now box_map of cut out AU

  box_mask_ncs_au, box_crystal_symmetry, \
        dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
       map_data = mask_data_ncs_au.as_double(),
       crystal_symmetry = tracking_data.crystal_symmetry,
       min_point = lower_bounds, max_point = upper_bounds, out = out)

  # Mask
  if params.output_files.box_mask_file and params.output_files.write_output_maps:
    # write out box_map NCS mask representing one AU of the NCS
    write_ccp4_map(
     box_crystal_symmetry,
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_mask_file),
      box_mask_ncs_au)
    print("Output NCS au as box (cut out) mask:  %s " %(
      os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_mask_file)), file = out)
    tracking_data.set_output_box_mask_info(
      file_name = os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_mask_file),
      crystal_symmetry = box_crystal_symmetry,
      origin = box_mask_ncs_au.origin(),
      all = box_mask_ncs_au.all())

  # Map
  box_map_ncs_au, box_crystal_symmetry, \
       dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
       soft_mask = tracking_data.params.map_modification.soft_mask,
       resolution = tracking_data.params.crystal_info.resolution,
       map_data = map_data_ncs_au.as_double(),
       crystal_symmetry = tracking_data.crystal_symmetry,
       min_point = lower_bounds, max_point = upper_bounds, out = out)

  half_map_data_list_au_box = []
  for hmdlu in half_map_data_list_au:
    hm_box_map_ncs_au, dummy_box_crystal_symmetry, \
       dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
       soft_mask = tracking_data.params.map_modification.soft_mask,
       resolution = tracking_data.params.crystal_info.resolution,
       map_data = hmdlu.as_double(),
       crystal_symmetry = tracking_data.crystal_symmetry,
       min_point = lower_bounds, max_point = upper_bounds, out = out)
    half_map_data_list_au_box.append(hm_box_map_ncs_au)

  if params.control.save_box_map_ncs_au:
       tracking_data.set_box_map_ncs_au_map_data(
       box_map_ncs_au_crystal_symmetry = box_crystal_symmetry,
       box_map_ncs_au_map_data = box_map_ncs_au,
       box_mask_ncs_au_map_data = box_mask_ncs_au,
       box_map_ncs_au_half_map_data_list = half_map_data_list_au_box,
       )

  if params.output_files.box_map_file:
    # write out NCS map as box_map (cut out region of map enclosed in box_mask)
    if params.output_files.write_output_maps:
      write_ccp4_map(box_crystal_symmetry,
        os.path.join(tracking_data.params.output_files.output_directory,
          params.output_files.box_map_file), box_map_ncs_au)
      print("Output NCS au as box (cut out) map:  %s " %(
      os.path.join(tracking_data.params.output_files.output_directory,
          params.output_files.box_map_file)), file = out)
      tracking_data.set_output_box_map_info(
        file_name = os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_map_file),
        crystal_symmetry = box_crystal_symmetry,
        origin = box_map_ncs_au.origin(),
        all = box_map_ncs_au.all())


  # Write out all the selected regions
  if params.output_files.write_output_maps:
    print("\nWriting out region maps. "+\
      "These superimpose on the NCS AU map \nand "+\
      "mask %s, %s\n" %(
        os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_map_file),
        os.path.join(tracking_data.params.output_files.output_directory, params.output_files.box_mask_file), ), file = out)

    map_files_written, remainder_regions_written = write_region_maps(params,
      map_data = map_data,
      tracking_data = tracking_data,
      ncs_group_obj = ncs_group_obj,
      remainder_ncs_group_obj = remainder_ncs_group_obj,
      out = out)

    # and pick up the remainder regions not already written
    remainder_map_files_written, dummy_remainder = write_region_maps(params,
      map_data = map_data,
      tracking_data = tracking_data,
      ncs_group_obj = remainder_ncs_group_obj,
      regions_to_skip = remainder_regions_written,
      out = out)
    map_files_written+= remainder_map_files_written
  else:
    map_files_written = []
  return map_files_written

def write_intermediate_maps(params,
    map_data = None,
    map_data_remaining = None,
    ncs_group_obj = None,
    tracking_data = None,
    out = sys.stdout):

  """Write intermediate map files (remainder maps)"""
  if map_data_remaining and params.output_files.remainder_map_file:
    write_ccp4_map(
       tracking_data.crystal_symmetry, params.output_files.remainder_map_file,
      map_data_remaining)
    print("Wrote output remainder map to %s" %(
       params.output_files.remainder_map_file), file = out)

  if params.segmentation.write_all_regions:
    for id in ncs_group_obj.selected_regions:
      region_mask = ncs_group_obj.edited_mask.deep_copy()
      s = (ncs_group_obj.edited_mask  ==  -1)
      s |=  (ncs_group_obj.edited_mask == id)
      region_mask = region_mask.set_selected(s, 1)
      region_mask = region_mask.set_selected(~s, 0)

      write_ccp4_map(tracking_data.crystal_symmetry,
          'mask_%d.ccp4' %id, region_mask)
      print("Wrote output mask for region %d to %s" %(id,
        "mask_%d.ccp4" %(id)), file = out)


def iterate_search(params,
      map_data_remaining = None,
      map_data = None,
      ncs_obj = None,
      ncs_group_obj = None,
      scattered_points = None,
      tracking_data = None,
      out = sys.stdout):

  """Iterate the process of segmenting map. First run standard run,
  then run again after removing all the parts that are already marked.
  """

  # Write out intermediate maps if desired
  if params.output_files.write_intermediate_maps:
    write_intermediate_maps(params,
      map_data = map_data,
      map_data_remaining = map_data_remaining,
      ncs_group_obj = ncs_group_obj,
      tracking_data = tracking_data,
      out = out)
  new_params = deepcopy(params)
  new_params.segmentation.iterate_with_remainder = False
  new_params.segmentation.density_threshold = None
  new_params.output_files.write_output_maps = False
  new_params.output_files.output_info_file = None
  if params.output_files.write_intermediate_maps:
    new_params.output_files.au_output_file_stem = \
      params.output_files.au_output_file_stem+"_cycle_2"
  else:
    new_params.output_files.au_output_file_stem = None

  fraction = params.segmentation.iteration_fraction
  if tracking_data.n_residues:
    new_n_residues = int(tracking_data.n_residues*fraction)
  new_solvent_fraction = max(0.001, min(0.999,
      1- (1-tracking_data.solvent_fraction)*fraction))

  new_tracking_data = deepcopy(tracking_data)
  if new_tracking_data.n_residues:
    new_tracking_data.set_n_residues(new_n_residues)
  new_tracking_data.set_solvent_fraction(new_solvent_fraction)
  new_tracking_data.set_origin_shift() # sets it to zero
  new_tracking_data.params.segmentation.starting_density_threshold = new_params.segmentation.starting_density_threshold # this is new
  print("\nIterating with remainder density", file = out)
  # NOTE: do not include pdb_hierarchy here unless you deep_copy it
  remainder_ncs_group_obj, dummy_remainder, remainder_tracking_data = run(
    None, params = new_params,
    map_data = map_data_remaining,
    ncs_obj = ncs_obj,
    target_scattered_points = scattered_points,
    tracking_data = new_tracking_data,
    is_iteration = True,
    out = out)
  if not remainder_ncs_group_obj: # Nothing to do
    return None

  # Combine the results to get remainder_id_dict
  #   remainder_id_dict[id_remainder] = id_nearby

  remainder_ncs_group_obj = combine_with_iteration(params,
     map_data = map_data,
     crystal_symmetry = tracking_data.crystal_symmetry,
     ncs_group_obj = ncs_group_obj,
     remainder_ncs_group_obj = remainder_ncs_group_obj,
     out = out)

  return remainder_ncs_group_obj

def bounds_overlap(lower = None, upper = None,
          other_lower = None, other_upper = None, tol = 1):
   """Return True if bounds (lower, uppper) vs (other_lower, other_upper)
   overlap"""
   for i in range(3):
     if       upper[i]+tol<other_lower[i]: return False
     if other_upper[i]+tol<lower[i]:       return False
   return True

def combine_with_iteration(params,
    map_data = None,
    crystal_symmetry = None,
    ncs_group_obj = None,
    remainder_ncs_group_obj = None,
    out = sys.stdout):

  """Combine original segmented regions with those found in iteration.
  """
  if not ncs_group_obj.selected_regions or not remainder_ncs_group_obj \
      or not remainder_ncs_group_obj.selected_regions:
    return None

  # see if any regions in ncs_obj overlap with remainder_ncs_group_obj...
  #   If so, combine


  remainder_id_dict = {}
  for id_remainder in remainder_ncs_group_obj.selected_regions:
    best_id = None
    best_overlaps = None
    remainder_centers = \
        remainder_ncs_group_obj.region_scattered_points_dict[id_remainder]
    # figure out typical distance between scattered_points...
    touching_dist = get_touching_dist(remainder_centers)

    # Notice bounds of remainder region:
    r_lower, r_upper = get_bounds(
       ncs_group_obj = remainder_ncs_group_obj, id = id_remainder)

    for id in ncs_group_obj.selected_regions:
      # Skip if not likely to be very close...
      lower, upper = get_bounds(ncs_group_obj = ncs_group_obj, id = id)
      if not bounds_overlap(lower = lower, upper = upper,
          other_lower = r_lower, other_upper = r_upper):
        continue

      test_centers = ncs_group_obj.region_scattered_points_dict[id]
      dist = get_closest_dist(test_centers, remainder_centers)
      if touching_dist is not None and dist>touching_dist:
        continue

      bool_region_mask = ncs_group_obj.co.expand_mask(
        id_to_expand = ncs_group_obj.original_id_from_id[id],
        expand_size = params.segmentation.expand_size+1) # just touching
      s = (bool_region_mask ==  True)
      s &=   (remainder_ncs_group_obj.edited_mask == id_remainder)
      overlaps = s.count(True)
      if best_overlaps is None or overlaps>best_overlaps:
        best_overlaps = overlaps
        best_id = id
    if best_overlaps:
      print("\nCombining remainder id %d with original id %d (overlaps = %d)" %(
        id_remainder, best_id, best_overlaps), file = out)
      remainder_id_dict[id_remainder] = best_id
  remainder_ncs_group_obj.remainder_id_dict = remainder_id_dict
  return remainder_ncs_group_obj

def get_touching_dist(centers, default = 100., min_dist = 8.):
  """Figure out typical distance between centers"""
  mean_dist = 0.
  mean_dist_n = 0.
  nskip = max(1, len(centers)//10) # try to get 10
  for i in range(0, len(centers), nskip):
     if i == 0:
       target = centers[1:]
     elif i == len(centers)-1:
       target = centers[:-1]
     else:
       target = centers[:i]
       target.extend(centers[i+1:])
     other = centers[i:i+1]
     if not target or not other: continue
     dist = get_closest_dist(target, other)
     if dist is not None:
       mean_dist+= dist
       mean_dist_n+= 1.
  if mean_dist_n>0:
    return max(min_dist, 2.0*mean_dist/mean_dist_n)
  else:
    return default

def get_grid_units(map_data = None, crystal_symmetry = None, radius = None,
     out = sys.stdout):
    """Get number of grid units representing the distance radius"""
    N_ = map_data.all()
    sx, sy, sz =  1/N_[0], 1/N_[1], 1/N_[2]
    sx_cart, sy_cart, sz_cart = crystal_symmetry.unit_cell().orthogonalize(
       [sx, sy, sz])
    grid_spacing = (sx_cart+sy_cart+sz_cart)/3.
    grid_units = int(radius/grid_spacing)
    min_cell_grid_units = min(N_[0], N_[1], N_[2])
    grid_units = max(1, min(grid_units, int(min_cell_grid_units/3)))
    print("Grid units representing %7.1f A will be %d" %(
       radius, grid_units), file = out)
    return grid_units

def cut_out_map(map_data = None, crystal_symmetry = None,
    soft_mask = None, soft_mask_radius = None, resolution = None,
    shift_origin = None,
    min_point = None, max_point = None, out = sys.stdout):
  """Cut out map from min_point to max_point, optionally smooth
  edge of new map.
  Returns new_map_data, new_crystal_symmetry,
    smoothed_mask_data, original_map_data

  NOTE: end point of map is max_point, so size of map (new all()) is
  (max_point-min_point+ (1, 1, 1))
  shrink unit cell, angles are the same
  NOTE 2: the origin of output map will be min_point (not 0, 0, 0).
  """
  from cctbx import uctbx
  from cctbx import maptbx
  na = map_data.all() # tuple with dimensions
  for i in range(3):
    assert min_point[i] >=  0
    assert max_point[i] < na[i]  # 2019-11-05 just na-1
  new_map_data = maptbx.copy(map_data, tuple(min_point), tuple(max_point))

  shrunk_uc = []
  for i in range(3):
    shrunk_uc.append(
     crystal_symmetry.unit_cell().parameters()[i]*new_map_data.all()[i]/na[i] )
  uc_params = crystal_symmetry.unit_cell().parameters()
  new_unit_cell_box = uctbx.unit_cell(
    parameters = (shrunk_uc[0], shrunk_uc[1], shrunk_uc[2],
        uc_params[3], uc_params[4], uc_params[5]))
  new_crystal_symmetry = crystal.symmetry(
    unit_cell = new_unit_cell_box, space_group = 'p1')

  if soft_mask:
    if soft_mask_radius is None:
       soft_mask_radius = resolution
       assert soft_mask_radius is not None

    original_map_data = new_map_data.deep_copy()
    new_map_data, smoothed_mask_data = set_up_and_apply_soft_mask(
       map_data = new_map_data,
       shift_origin = shift_origin,
       crystal_symmetry = new_crystal_symmetry,
       resolution = resolution,
       radius = soft_mask_radius, out = out)
  else:
    original_map_data = None
    smoothed_mask_data = None

  return new_map_data, new_crystal_symmetry, \
    smoothed_mask_data, original_map_data

def get_zero_boundary_map(
   map_data = None,
   grid_units_for_boundary = None,
   crystal_symmetry = None,
   radius = None):

    """Get a map with a zero around the boundary using
    maptbx.zero_boundary_box_map"""

    assert grid_units_for_boundary or (crystal_symmetry and radius)

    # grid_units is how many grid units are about equal to soft_mask_radius
    if grid_units_for_boundary is None:
      grid_units = get_grid_units(map_data = map_data,
        crystal_symmetry = crystal_symmetry, radius = radius, out = null_out())
      grid_units = int(0.5+0.5*grid_units)
    else:
      grid_units = grid_units_for_boundary

    from cctbx import maptbx
    zero_boundary_map = maptbx.zero_boundary_box_map(
       map_data, grid_units).result()
    return zero_boundary_map

def set_up_and_apply_soft_mask(map_data = None, shift_origin = None,
  crystal_symmetry = None, resolution = None,
  grid_units_for_boundary = None,
  radius = None, out = None):
    """Set up and apply a soft boundary mask to map_data"""
    if out is None:
      from libtbx.utils import null_out
      out = null_out()

    acc = map_data.accessor()
    map_data = map_data.shift_origin()
    new_acc = map_data.accessor()

    # Add soft boundary to mean around outside of mask

    zero_boundary_map = get_zero_boundary_map(
     map_data = map_data,
     grid_units_for_boundary = grid_units_for_boundary,
     crystal_symmetry = crystal_symmetry,
     radius = radius)

    # this map is zero's around the edge and 1 in the middle
    # multiply zero_boundary_map--smoothed & new_map_data and return
    print("Applying soft mask to boundary of cut out map", file = out)
    new_map_data, smoothed_mask_data = apply_soft_mask(map_data = map_data,
          mask_data = zero_boundary_map,
          rad_smooth = resolution,
          crystal_symmetry = crystal_symmetry,
          out = out)
    if new_acc !=  acc:
      new_map_data.reshape(acc)
      smoothed_mask_data.reshape(acc)
    return new_map_data, smoothed_mask_data

def apply_shift_to_pdb_hierarchy(
    origin_shift = None,
    crystal_symmetry = None,
    pdb_hierarchy = None, out = sys.stdout):
  """Apply a shift to a pdb_hierarchy.  Deprecated. Crystal_symmetry
  is ignored"""

  if origin_shift is not None:
    sites_cart = pdb_hierarchy.atoms().extract_xyz()
    sites_cart_shifted = sites_cart+\
      flex.vec3_double(sites_cart.size(), origin_shift)
    pdb_hierarchy.atoms().set_xyz(sites_cart_shifted)

  return pdb_hierarchy

def apply_origin_shift(origin_shift = None,
    ncs_object = None,
    shifted_ncs_object = None,
    pdb_hierarchy = None,
    target_hierarchy = None,
    map_data = None,
    shifted_map_file = None,
    shifted_pdb_file = None,
    shifted_ncs_file = None,
    tracking_data = None,
    out = sys.stdout):

  """Update tracking_data and write shifted_map_file, shifted_pdb_file
  after origin shift is applied (it is applied before coming in).

  Returns: shifted_pdb_file, ncs_obj, pdb_hierarchy, target_hierarchy,
      tracking_data
  """

  if shifted_map_file:
      write_ccp4_map(tracking_data.crystal_symmetry,
      shifted_map_file,
      map_data)
      print("Wrote shifted map to %s" %(
        shifted_map_file), file = out)
      tracking_data.set_shifted_map_info(file_name =
        shifted_map_file,
        crystal_symmetry = tracking_data.crystal_symmetry,
        origin = map_data.origin(),
        all = map_data.all())
  if origin_shift: # Note origin shift does not change crystal_symmetry

    if pdb_hierarchy:
      pdb_hierarchy = apply_shift_to_pdb_hierarchy(
       origin_shift = origin_shift,
       crystal_symmetry = tracking_data.crystal_symmetry,
       pdb_hierarchy = pdb_hierarchy,
       out = out)

    if target_hierarchy:
      target_hierarchy = apply_shift_to_pdb_hierarchy(
       origin_shift = origin_shift,
       crystal_symmetry = tracking_data.crystal_symmetry,
       pdb_hierarchy = target_hierarchy,
       out = out)

    from scitbx.math import  matrix
    if ncs_object and not shifted_ncs_object:
      shfted_ncs_object = ncs_object.coordinate_offset(
         coordinate_offset = matrix.col(origin_shift))


  if shifted_pdb_file and pdb_hierarchy:
      shifted_pdb_file = pdb_hierarchy.write_pdb_or_mmcif_file(
        target_filename = shifted_pdb_file,
        crystal_symmetry = tracking_data.crystal_symmetry)

      print("Wrote shifted pdb file to %s" %(
        shifted_pdb_file), file = out)
      tracking_data.set_shifted_pdb_info(file_name = shifted_pdb_file,
      n_residues = pdb_hierarchy.overall_counts().n_residues)


  if shifted_ncs_file and shifted_ncs_object:
      shifted_ncs_object.format_all_for_group_specification(
         file_name = shifted_ncs_file)
      print("Wrote %s NCS operators for shifted map to %s" %(
         shifted_ncs_object.max_operators(),
         shifted_ncs_file), file = out)
      if tracking_data.input_ncs_info.has_updated_operators():
        print("NOTE: these may include additional operators added to fill the cell"+\
        " or\nhave fewer operators if not all applied.", file = out)
      tracking_data.set_shifted_ncs_info(file_name = shifted_ncs_file,
        number_of_operators = shifted_ncs_object.max_operators(),
        is_helical_symmetry = tracking_data.input_ncs_info.is_helical_symmetry)
      tracking_data.shifted_ncs_info.show_summary(out = out)

  return shifted_pdb_file, shifted_ncs_object, pdb_hierarchy, \
     target_hierarchy, tracking_data

def restore_pdb(params, tracking_data = None, out = sys.stdout):
  """Take params.input_files.pdb_to_restore and restore it to
   its original location (before origin shifts) and write to
   params.output_files.restored_pdb"""

  if not params.output_files.restored_pdb:
    params.output_files.restored_pdb = \
       params.input_files.pdb_to_restore[:-4]+"_restored.pdb"
  print("Shifting origin of %s and writing to %s" %(
    params.input_files.pdb_to_restore,
    params.output_files.restored_pdb), file = out)
  os = tracking_data.origin_shift
  origin_shift = (-os[0], -os[1], -os[2])
  print("Origin shift will be: %.1f  %.1f  %.1f "%(origin_shift), file = out)

  from iotbx.pdb.utils import get_pdb_hierarchy
  pdb_hierarchy = get_pdb_hierarchy(
     file_name = params.input_files.pdb_to_restore)

  pdb_hierarchy = apply_shift_to_pdb_hierarchy(
    origin_shift = origin_shift,
    crystal_symmetry = tracking_data.crystal_symmetry,
    pdb_hierarchy = pdb_hierarchy,
    out = out)

  params.output_files.restored_pdb = pdb_hierarchy.write_pdb_or_mmcif_file(
      crystal_symmetry = tracking_data.crystal_symmetry,
      target_filename = params.output_files.restored_pdb)

  print("Wrote restored pdb file to %s" %(
     params.output_files.restored_pdb), file = out)

def find_threshold_in_map(target_points = None,
      map_data = None,
      require_at_least_target_points = None,
      iter_max = 10):

  """Find threshold value such that target_points in map_data are greater
  than threshold."""

  map_1d = map_data.as_1d()
  map_mean = map_1d.min_max_mean().mean
  map_max = map_1d.min_max_mean().max
  map_min = map_1d.min_max_mean().min

  cutoff = map_mean
  low = map_min
  high = map_max

  best_cutoff = None
  best_score = None
  for iter in range(iter_max):
    s = (map_1d >cutoff)
    n_cutoff = s.count(True)
    if (not require_at_least_target_points) or (n_cutoff >=  target_points):
      score = abs(n_cutoff-target_points)
    else:
      score = 1.e+10 # allow it but anything above cutoff will be better

    if best_score is None or score < best_score:
      best_cutoff = cutoff
      best_score = score

    if n_cutoff  ==  target_points:
      return best_cutoff
    elif n_cutoff < target_points: # lower it
      high = cutoff
      cutoff = 0.5*(cutoff+low)
    else:  # raise it
      low = cutoff
      cutoff = 0.5*(cutoff+high)
  return best_cutoff


def remove_points(mask, remove_points = None):
  """Remove all points in remove_points from mask. Points are
  indices obtained with selections"""
  keep_points = (remove_points == False)
  new_mask = (mask & keep_points)
  return new_mask

def get_ncs_sites_cart(sites_cart = None, ncs_id = None,
     ncs_obj = None, unit_cell = None, ncs_in_cell_only = True):
  """Return points ncs-related to sites_cart"""
  ncs_sites_cart = flex.vec3_double()
  if not ncs_obj or not ncs_obj.ncs_groups() or not ncs_obj.ncs_groups()[0] or \
     not ncs_obj.ncs_groups()[0].translations_orth():
    return ncs_sites_cart

  # identify ncs-related points
  ncs_group = ncs_obj.ncs_groups()[0]
  identity_op = ncs_group.identity_op_id()
  ncs_sites_cart = flex.vec3_double()
  for xyz_cart in sites_cart:
    for i0 in range(len(ncs_group.translations_orth())):
      if i0 == identity_op: continue
      if ncs_id is not None and i0!= ncs_id: continue
      r = ncs_group.rota_matrices_inv()[i0] # inverse maps pos 0 on to pos i
      t = ncs_group.translations_orth_inv()[i0]
      new_xyz_cart = r * matrix.col(xyz_cart) + t
      ncs_sites_cart.append(new_xyz_cart)
  if ncs_in_cell_only:
    new_sites_cart = flex.vec3_double()
    ncs_sites_frac = unit_cell.fractionalize(ncs_sites_cart)
    for site_frac, site_cart in zip(ncs_sites_frac, ncs_sites_cart):
      if site_frac[0]>= 0 and site_frac[0]<= 1 and  \
         site_frac[1]>= 0 and site_frac[1]<= 1 and  \
         site_frac[2]>= 0 and site_frac[2]<= 1:
        new_sites_cart.append(site_cart)
    ncs_sites_cart = new_sites_cart

  return ncs_sites_cart

def get_ncs_mask(map_data = None, unit_cell = None, ncs_object = None,
   starting_mask = None, radius = None, expand_radius = None,
    overall_mask = None,
   every_nth_point = None):

  """Get a mask marking the asymmetric unit of the NCS of this map (smallest
  region which when NCS is applied, generates the entire map)

  Returns working_au_mask, working_ncs_mask
  """
  assert every_nth_point is not None
  if not expand_radius: expand_radius = 2.*radius

  working_au_mask = starting_mask.deep_copy()

  working_ncs_mask = mask_from_sites_and_map(  # empty ncs mask
    map_data = map_data, unit_cell = unit_cell,
    sites_cart = flex.vec3_double(), radius = radius, overall_mask = overall_mask)

  au_points_last = working_au_mask.count(True)
  ncs_points_last = working_ncs_mask.count(True)

  max_tries = 10000
  for ii in range(max_tries): # just a big number; should take just a few

    # Find all points in au (sample every_nth_point in grid)

    au_sites_cart = get_marked_points_cart(mask_data = working_au_mask,
     unit_cell = unit_cell, every_nth_point = every_nth_point,
     boundary_radius = radius)

    # Find all points ncs-related to marked point in mask
    ncs_sites_cart = get_ncs_sites_cart(sites_cart = au_sites_cart,
       ncs_obj = ncs_object, unit_cell = unit_cell, ncs_in_cell_only = True)

    # Expand au slightly with all points near to au_sites_cart
    new_au_mask = mask_from_sites_and_map(
      map_data = map_data, unit_cell = unit_cell,
      sites_cart = au_sites_cart, radius = radius, overall_mask = overall_mask)
    working_au_mask = (working_au_mask | new_au_mask) # add on to existing
    keep_points = (working_ncs_mask == False)  # cross off those in  ncs
    working_au_mask = (working_au_mask & keep_points)

    # mark ncs au with all points not in au that are close to ncs_sites_cart
    new_ncs_mask = mask_from_sites_and_map(
      map_data = map_data, unit_cell = unit_cell,
      sites_cart = ncs_sites_cart, radius = radius, overall_mask = overall_mask)
    keep_points = (working_au_mask == False)  # cross off those in au
    new_ncs_mask = (new_ncs_mask & keep_points)
    working_ncs_mask = (new_ncs_mask | working_ncs_mask) # add on to existing

    au_points = working_au_mask.count(True)
    ncs_points = working_ncs_mask.count(True)
    if au_points == au_points_last and ncs_points == ncs_points_last:
      break
    au_points_last = au_points
    ncs_points_last = ncs_points
    # Now expand the au and repeat

    working_au_mask = mask_from_sites_and_map(
      map_data = map_data, unit_cell = unit_cell,
      sites_cart = au_sites_cart, radius = expand_radius, overall_mask = overall_mask)
    keep_points = (working_ncs_mask == False)  # cross off those in  ncs
    working_au_mask = (working_au_mask & keep_points)

  return working_au_mask, working_ncs_mask

def renormalize_map_data(
  map_data = None, solvent_fraction = None):

  """Return map_data object normalized so the RMS of the non-solvent
  region is about 1. Guessed without actually masking the solvent region."""

  sd = max(0.0001, map_data.sample_standard_deviation())
  if solvent_fraction >=  10.: solvent_fraction = solvent_fraction/100.
  solvent_fraction = min(0.999, max(0.001, solvent_fraction))
  scaled_sd = sd/(1-solvent_fraction)**0.5
  map_data = (map_data-map_data.as_1d().min_max_mean().mean)/scaled_sd
  return map_data


def mask_from_sites_and_map(
    map_data = None, unit_cell = None,
    sites_cart = None, radius = None, overall_mask = None):
  """Generate a mask around sites"""
  assert radius is not None
  from cctbx import maptbx

  sel = maptbx.grid_indices_around_sites(
      unit_cell  = unit_cell,
      fft_n_real = map_data.focus(),
      fft_m_real = map_data.all(),
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), radius))
  map_data_1d = map_data.as_1d()
  mask = (map_data_1d == 0 and map_data_1d == 1)  # 1D bool array all False
  mask.set_selected(sel, True)  # mark points around sites
  mask.reshape(map_data.accessor())
  if overall_mask:
    assert overall_mask.all() == mask.all()
    mask = (mask & overall_mask)
  return mask

def set_radius(unit_cell = None, map_data = None, every_nth_point = None):
  """Set radius so that radius will capture all points on grid if sampled
  on every_nth_point """
  a, b, c = unit_cell.parameters()[:3]
  nx, ny, nz = map_data.all()
  # furthest possible minimum distance between grid points
  max_diagonal_between_sampled = every_nth_point*(
      (a/nx)**2+(b/ny)**2+(c/nz)**2)**0.5
  radius = max_diagonal_between_sampled*0.55  # big enough to cover everything
  return radius

def get_marked_points_cart(mask_data = None, unit_cell = None,
   every_nth_point = 3, boundary_radius = None):
  """Return list of cartesian coordinates of grid points that are marked.
  Only sample every every_nth_point in each direction..."""
  assert mask_data.origin()  ==  (0, 0, 0)
  nx, ny, nz = mask_data.all()
  if boundary_radius:
    # How far from edges shall we stay:
    grid_frac = (1./nx, 1./ny, 1./nz)
    grid_orth = unit_cell.orthogonalize(grid_frac)
    boundary_grid_points = 0
    for go in grid_orth:
      bgp = int(0.99+boundary_radius/go)
      boundary_grid_points = max(boundary_grid_points, bgp)
  else:
    boundary_grid_points = 0

  marked_points = maptbx.marked_grid_points(
    map_data = mask_data,
    every_nth_point = every_nth_point).result()
  sites_frac = flex.vec3_double()
  boundary_points_skipped = 0
  for grid_point in marked_points:
    if boundary_grid_points:
      if \
         grid_point[0]<boundary_grid_points or \
         grid_point[0]>nx-boundary_grid_points or \
         grid_point[1]<boundary_grid_points or \
         grid_point[1]>ny-boundary_grid_points or \
         grid_point[2]<boundary_grid_points or \
         grid_point[2]>nz-boundary_grid_points:  # XXX was typo previously
        boundary_points_skipped+= 1
        continue
    sites_frac.append(
        (grid_point[0]/nx,
         grid_point[1]/ny,
         grid_point[2]/nz))

  sites_cart = unit_cell.orthogonalize(sites_frac)
  return sites_cart

def get_overall_mask(
    map_data = None,
    mask_threshold = None,
    fraction_of_max_mask_threshold = None,
    use_solvent_content_for_threshold = None, # use instead of fraction_of
    mask_padding_fraction = None,
    solvent_fraction = None,
    crystal_symmetry = None,
    radius = None,
    resolution = None,
    d_max = 100000.,
    out = sys.stdout):


  """Get a mask around region containing macromolecule. Based on
  finding region of sd_map > threshold, setting threshold to get
  about solvent_fraction outside of marked region.
  Returns overall_mask, max_in_sd_map, sd_map.
  """

  # This routine cannot use mask_data with origin != (0,0,0)
  if map_data.origin() != (0,0,0):
    print("Map origin must be at (0,0,0) for get_overall_mask")
    assert map_data.origin() == (0,0,0)  # Map origin must be at (0,0,0)

  # Make a local SD map from our map-data
  from cctbx.maptbx import crystal_gridding
  from cctbx import sgtbx
  cg = crystal_gridding(
        unit_cell = crystal_symmetry.unit_cell(),
        space_group_info = sgtbx.space_group_info(number = 1), # Always
        pre_determined_n_real = map_data.all())

  if not resolution:
    from cctbx.maptbx import d_min_from_map
    resolution = d_min_from_map(
      map_data, crystal_symmetry.unit_cell(), resolution_factor = 1./4.)
    print("\nEstimated resolution of map: %6.1f A\n" %(
     resolution), file = out)

  if radius:
    smoothing_radius = 2.*radius
  else:
    smoothing_radius = 2.*resolution

  from iotbx.map_manager import map_manager
  mm = map_manager(map_data = map_data,
       unit_cell_grid = map_data.all(),
       unit_cell_crystal_symmetry = crystal_symmetry,
       wrapping = False)
  map_coeffs = mm.map_as_fourier_coefficients(
       d_min = resolution, d_max = d_max)

  if not map_coeffs:
    raise Sorry("No map coeffs obtained")

  complete_set = map_coeffs.complete_set()
  stol = flex.sqrt(complete_set.sin_theta_over_lambda_sq().data())
  import math
  w = 4 * stol * math.pi * smoothing_radius
  sphere_reciprocal = 3 * (flex.sin(w) - w * flex.cos(w))/flex.pow(w, 3)
  try:
    temp = complete_set.structure_factors_from_map(
      flex.pow2(map_data-map_data.as_1d().min_max_mean().mean))
  except Exception as e:
    print(e, file = out)
    print ("The sampling of the map appears to be too low for a "+
      "\nresolution of %s. Using a larger value for resolution" %(
       resolution), file = out)
    from cctbx.maptbx import d_min_from_map
    resolution = d_min_from_map(
      map_data, crystal_symmetry.unit_cell(), resolution_factor = 1./4.)
    print("\nEstimated resolution of map: %6.1f A\n" %(
     resolution), file = out)
    map_coeffs = map_coeffs.resolution_filter(d_min = resolution, d_max = d_max)
    complete_set = map_coeffs.complete_set()
    stol = flex.sqrt(complete_set.sin_theta_over_lambda_sq().data())
    import math
    w = 4 * stol * math.pi * smoothing_radius
    sphere_reciprocal = 3 * (flex.sin(w) - w * flex.cos(w))/flex.pow(w, 3)
    temp = complete_set.structure_factors_from_map(
      flex.pow2(map_data-map_data.as_1d().min_max_mean().mean))


  fourier_coeff = complete_set.array(data = temp.data()*sphere_reciprocal)
  sd_map = fourier_coeff.fft_map(
      crystal_gridding = cg,
      ).apply_volume_scaling().real_map_unpadded()
  assert sd_map.all() == map_data.all()
  # now use sd_map

  # First mask out the map based on threshold
  mm = sd_map.as_1d().min_max_mean()
  max_in_sd_map = mm.max
  mean_in_map = mm.mean
  min_in_map = mm.min
  print("Highest value in SD map is %7.2f. Mean is %7.2f .  Lowest is %7.2f " %(
    max_in_sd_map,
    mean_in_map,
    min_in_map), file = out)
  if fraction_of_max_mask_threshold and (
        (not solvent_fraction) or (not use_solvent_content_for_threshold)):
    mask_threshold = fraction_of_max_mask_threshold*max_in_sd_map
    print("Using fraction of max as threshold: %.3f " %(
        fraction_of_max_mask_threshold), \
        "which is threshold of %.3f" %(mask_threshold), file = out)
    if mask_padding_fraction:
      # Adjust threshold to increase by mask_padding_fraction, proportional
      #  to fraction available
      overall_mask = (sd_map>=  mask_threshold)
      current_above_threshold = overall_mask.count(True)/overall_mask.size()
      # current+(1-current)*pad
      additional_padding = (1-current_above_threshold)*mask_padding_fraction
      target_above_threshold = min(
        0.99, current_above_threshold+additional_padding)
      print("Target with padding of %.2f will be %.2f" %(
         mask_padding_fraction, target_above_threshold), file = out)
      solvent_fraction = (1-target_above_threshold)
      mask_threshold = None

  if mask_threshold:
    print("Cutoff for mask will be input threshold", file = out)
    threshold = mask_threshold
  else:  # guess based on solvent_fraction
    if solvent_fraction is None:
       print("Guessing solvent fraction of 0.9", file = out)
       solvent_fraction = 0.9 # just guess
    threshold = find_threshold_in_map(target_points = int(
      (1.-solvent_fraction)*sd_map.size()),
      map_data = sd_map)
    print("Cutoff will be threshold marking about %7.1f%% of cell" %(
      100.*(1.-solvent_fraction)), file = out)

  overall_mask = (sd_map>=  threshold)
  print("Model region of map "+\
    "(density above %7.3f )" %( threshold) +" includes %7.1f%% of map" %(
      100.*overall_mask.count(True)/overall_mask.size()), file = out)
  return overall_mask, max_in_sd_map, sd_map

def get_skew(data = None):
  """Get the skew of flex.double() array data"""
  mean = data.min_max_mean().mean
  sd = data.standard_deviation_of_the_sample()
  x = data-mean
  return (x**3).min_max_mean().mean/sd**3

def get_kurtosis(data = None):
  """Get the kurtosis of flex.double() array data"""
  mean = data.min_max_mean().mean
  sd = data.standard_deviation_of_the_sample()
  x = data-mean
  return (x**4).min_max_mean().mean/sd**4

def score_map(map_data = None,
        sharpening_info_obj = None,
        solvent_fraction = None,
        fraction_occupied = None,
        wrapping = None,
        sa_percent = None,
        region_weight = None,
        max_regions_to_test = None,
        scale_region_weight = False,
        out = sys.stdout):
  """Score a map based on adjusted surface area, kurtosis,
   or adjusted_path_length.
   Return updated sharpening_info_obj
  """
  if sharpening_info_obj:
    solvent_fraction = sharpening_info_obj.solvent_fraction
    wrapping = sharpening_info_obj.wrapping
    fraction_occupied = sharpening_info_obj.fraction_occupied
    sa_percent = sharpening_info_obj.sa_percent
    region_weight = sharpening_info_obj.region_weight
    max_regions_to_test = sharpening_info_obj.max_regions_to_test
  else:
    sharpening_info_obj = sharpening_info()
  if solvent_fraction is None:  # skip SA score
    sharpening_info_obj.adjusted_sa = 0.
    assert sharpening_info_obj.sharpening_target == 'kurtosis'
  else:  # usual
    map_data = renormalize_map_data(
       map_data = map_data, solvent_fraction = solvent_fraction)

    target_in_all_regions = map_data.size()*fraction_occupied*(1-solvent_fraction)
    print("\nTarget number of points in all regions: %.0f" %(
      target_in_all_regions), file = out)

    threshold = find_threshold_in_map(target_points = int(
         target_in_all_regions), map_data = map_data)
    print("Cutoff will be threshold of %7.2f marking %7.1f%% of cell" %(
              threshold, 100.*(1.-solvent_fraction)), file = out)

    co, sorted_by_volume, min_b, max_b = get_co(
      map_data = map_data.deep_copy(),
       threshold = threshold, wrapping = wrapping)

    if len(sorted_by_volume)<2:
      return sharpening_info_obj# skip it, nothing to do

    target_sum =  sa_percent* target_in_all_regions*0.01
    print("Points for %.1f percent of target in all regions: %.1f" %(
        sa_percent, target_sum), file = out)

    cntr = 0
    sum_v = 0.
    sum_new_v = 0.
    for p in sorted_by_volume[1:max_regions_to_test+2]:
      cntr+= 1
      v, i = p
      sum_v+= v
      bool_expanded = co.expand_mask(id_to_expand = i, expand_size = 1)
      new_v = bool_expanded.count(True)-v
      sum_new_v+= new_v
      sa_ratio = new_v/v
      if sum_v>= target_sum: break
    sa_ratio = sum_new_v/max(1., sum_v) # ratio of SA to volume
    regions = len(sorted_by_volume[1:])
    normalized_regions = regions/max(1, target_in_all_regions)
    skew = get_skew(map_data.as_1d())

    if scale_region_weight:
      solvent_fraction_std = 0.85 # typical value, ends up as scale on weight
      region_weight_scale = (1.-solvent_fraction)/(1.-solvent_fraction_std)
      region_weight_use = region_weight*region_weight_scale
    else:
      region_weight_use = region_weight
    sharpening_info_obj.adjusted_sa = \
       sa_ratio - region_weight_use*normalized_regions
    sharpening_info_obj.sa_ratio = sa_ratio
    sharpening_info_obj.normalized_regions = normalized_regions

  sharpening_info_obj.kurtosis = get_kurtosis(map_data.as_1d())

  if sharpening_info_obj.sharpening_target == 'adjusted_path_length':
    sharpening_info_obj.adjusted_path_length = get_adjusted_path_length(
      map_data = map_data,
      resolution = sharpening_info_obj.resolution,
      crystal_symmetry = sharpening_info_obj.crystal_symmetry,
      out = out)
  else:
    sharpening_info_obj.adjusted_path_length = None

  if sharpening_info_obj.sharpening_target == 'kurtosis':
    sharpening_info_obj.score = sharpening_info_obj.kurtosis
  if sharpening_info_obj.sharpening_target == 'adjusted_path_length':
    sharpening_info_obj.score = sharpening_info_obj.adjusted_path_length
  else:
    sharpening_info_obj.score = sharpening_info_obj.adjusted_sa

  return sharpening_info_obj

def get_adjusted_path_length(
    map_data = None,
    crystal_symmetry = None,
    resolution = None,
    out = sys.stdout):
  """Calculate adjusted path length for tubes of density obtained after
  segmenting map. Uses trace_and_build"""

  try:
    from phenix.autosol.trace_and_build import trace_and_build
  except Exception as e: # Not available
    return 0

  from phenix.programs.trace_and_build import master_phil_str
  import iotbx.phil
  tnb_params = iotbx.phil.parse(master_phil_str).extract()
  tnb_params.crystal_info.resolution = resolution
  tnb_params.strategy.retry_long_branches = False
  tnb_params.strategy.correct_segments = False
  tnb_params.strategy.split_and_join = False
  tnb_params.strategy.vary_sharpening = []
  tnb_params.strategy.get_path_length_only = True
  tnb_params.trace_and_build.find_helices_strands = False
  tnb = trace_and_build(
      params = tnb_params,
      map_data = map_data,
      crystal_symmetry = crystal_symmetry,
      origin_cart = (0, 0, 0),
      origin = (0, 0, 0),
      log = out)
  tnb.run()
  return tnb.adjusted_path_length

def sharpen_map_with_si(sharpening_info_obj = None,
     f_array_normalized = None,
     f_array = None, phases = None,
     map_data = None,
     overall_b = None,
     resolution = None,
     out = sys.stdout):
  """Sharpen a map using information in si object.
   Return map_and_b object with map_data and b_iso"""
  si = sharpening_info_obj


  if si.sharpening_method == 'no_sharpening':
     return map_and_b_object(map_data = map_data)

  if map_data and (not f_array or not phases):
    map_coeffs, dummy = get_f_phases_from_map(map_data = map_data,
       crystal_symmetry = si.crystal_symmetry,
       d_min = si.resolution,
       d_min_ratio = si.d_min_ratio,
       return_as_map_coeffs = True,
       scale_max = si.scale_max,
       out = out)
    f_array, phases = map_coeffs_as_fp_phi(map_coeffs)

  if si.remove_aniso:
    if si.use_local_aniso and \
      (si.local_aniso_in_local_sharpening or
       (si.local_aniso_in_local_sharpening is None and si.ncs_copies == 1)) and \
         si.original_aniso_obj: # use original
      aniso_obj = si.original_aniso_obj
      print("\nRemoving aniso from map using saved aniso object before sharpening\n", file = out)
    else:
      print("\nRemoving aniso from map before sharpening\n", file = out)
      aniso_obj = None
    from cctbx.maptbx.refine_sharpening import analyze_aniso
    f_array, f_array_ra = analyze_aniso(
        aniso_obj = aniso_obj,
        remove_aniso = si.remove_aniso,
        f_array = f_array, resolution = si.resolution, out = out)

  if si.is_model_sharpening() or si.is_half_map_sharpening():
    from cctbx.maptbx.refine_sharpening import scale_amplitudes
    ff = f_array.phase_transfer(phase_source = phases, deg = True)
    map_and_b = scale_amplitudes(
      map_coeffs = f_array.phase_transfer(phase_source = phases, deg = True),
      si = si, overall_b = overall_b, out = out)
    return map_and_b

  elif si.is_resolution_dependent_sharpening():
    if f_array_normalized is None:

      from cctbx.maptbx.refine_sharpening import get_sharpened_map, \
       quasi_normalize_structure_factors
      (d_max, d_min) = f_array.d_max_min()
      if not f_array.binner():
        f_array.setup_binner(n_bins = si.n_bins, d_max = d_max, d_min = d_min)
      f_array_normalized = quasi_normalize_structure_factors(
          f_array, set_to_minimum = 0.01)
    map_data = get_sharpened_map(ma = f_array_normalized, phases = phases,
       b = si.resolution_dependent_b, resolution = si.resolution, n_real = si.n_real,
       d_min_ratio = si.d_min_ratio)
    return map_and_b_object(map_data = map_data)

  else:
    map_and_b = apply_sharpening(n_real = si.n_real,
          f_array = f_array, phases = phases,
          sharpening_info_obj = si,
          crystal_symmetry = si.crystal_symmetry,
          out = null_out())
    return map_and_b

def put_bounds_in_range(
     lower_bounds = None, upper_bounds = None,
     box_size = None,
     n_buffer = None,
     n_real = None, out = sys.stdout):
  """Put lower and upper inside (0, n_real-1) and try to
    make size at least minimum"""

  new_lb = []
  new_ub = []
  print("Putting bounds in range...(%s, %s, %s) to (%s, %s, %s)" %(
       tuple(list(lower_bounds)+list(upper_bounds))), file = out)
  if n_buffer:
     print("Buffer of %s added" %(n_buffer), file = out)
  for lb, ub, ms, nr in zip(lower_bounds, upper_bounds, box_size, n_real):
    if ub<lb: ub = lb
    if lb >ub: lb = ub
    extra = (ms-(ub-lb))//2
    lb = lb-extra
    ub = ub+extra

    if n_buffer:
       lb = lb-n_buffer
       ub = ub+n_buffer

    if lb<0:
      shift = -lb
      lb+= shift
      ub+= shift
    boundary = int(ms-(ub-lb+1))//2
    if boundary>0:
       lb = lb-boundary
       ub = ub+boundary
    if lb<0: lb = 0
    if ub>= nr: ub = nr-1  # 2019-11-05 cannot go beyond na-1
    new_lb.append(lb)
    new_ub.append(ub)
  print("New bounds ...(%s, %s, %s) to (%s, %s, %s)" %(
       tuple(list(new_lb)+list(new_ub))), file = out)
  return tuple(new_lb), tuple(new_ub)

def get_iterated_solvent_fraction(map = None,
    verbose = None,
    resolve_size = None,
    crystal_symmetry = None,
    mask_padding_fraction = None,
    fraction_of_max_mask_threshold = None,
    solvent_content = None,
    use_solvent_content_for_threshold = None, # use instead of fraction_of_..
    cell_cutoff_for_solvent_from_mask = None,
    mask_resolution = None,
    return_mask_and_solvent_fraction = None,
    out = sys.stdout):
  """Use the iterated_solvent_fraction function to estimate solvent
  fraction in a map"""

  if cell_cutoff_for_solvent_from_mask and \
   crystal_symmetry.unit_cell().volume() > cell_cutoff_for_solvent_from_mask**3:
     #go directly to low_res_mask
     return get_solvent_fraction_from_low_res_mask(
      crystal_symmetry = crystal_symmetry,
      map_data = map.deep_copy(),
      mask_padding_fraction = mask_padding_fraction,
      fraction_of_max_mask_threshold = fraction_of_max_mask_threshold,
      solvent_content = solvent_content,
      use_solvent_content_for_threshold = use_solvent_content_for_threshold,
       return_mask_and_solvent_fraction = return_mask_and_solvent_fraction,
      mask_resolution = mask_resolution, out = out)

  try:
    from phenix.autosol.map_to_model import iterated_solvent_fraction
    solvent_fraction, overall_mask = iterated_solvent_fraction(
      crystal_symmetry = crystal_symmetry,
      map_as_double = map,
      verbose = verbose,
      resolve_size = resolve_size,
      out = out)
    if solvent_fraction<= 0.989:  # means that it was 0.99 which is hard limit
      if return_mask_and_solvent_fraction:
        return overall_mask, solvent_fraction
      else:
        return solvent_fraction
    else:  # use backup method
      return get_solvent_fraction_from_low_res_mask(
        crystal_symmetry = crystal_symmetry,
        map_data = map.deep_copy(),
        mask_padding_fraction = mask_padding_fraction,
        fraction_of_max_mask_threshold = fraction_of_max_mask_threshold,
      solvent_content = solvent_content,
       use_solvent_content_for_threshold = use_solvent_content_for_threshold,
        return_mask_and_solvent_fraction = return_mask_and_solvent_fraction,
        mask_resolution = mask_resolution, out = out)
  except Exception as e:
    # catch case where map was not on proper grid
    if str(e).find("sym equiv of a grid point must be a grid point")>-1:
      print("\nSpace group:%s \n Unit cell: %s \n Gridding: %s \nError message: %s" %(
        crystal_symmetry.space_group().info(),
        str(crystal_symmetry.unit_cell().parameters()),
        str(map.all()), str(e)), file = out)
      raise Sorry(
      "The input map seems to be on a grid incompatible with crystal symmetry"+
         "\n(symmetry equivalents of a grid point must be on "+
          "an integer grid point)")
    elif str(e).find("maximum size for resolve is")>-1:
      raise Sorry(str(e)+
       "\nIt may be possible to go on by supplying solvent content"+
      "or molecular_mass")
    # Try to get solvent fraction with low_res mask
    return get_solvent_fraction_from_low_res_mask(
      crystal_symmetry = crystal_symmetry,
      map_data = map.deep_copy(),
      mask_padding_fraction = mask_padding_fraction,
      fraction_of_max_mask_threshold = fraction_of_max_mask_threshold,
      solvent_content = solvent_content,
      use_solvent_content_for_threshold = use_solvent_content_for_threshold,
      return_mask_and_solvent_fraction = return_mask_and_solvent_fraction,
      mask_resolution = mask_resolution, out = out)

def get_solvent_fraction_from_low_res_mask(
      crystal_symmetry = None, map_data = None,
      fraction_of_max_mask_threshold = None,
      solvent_content = None,
      use_solvent_content_for_threshold = None, # use instead of fraction_of
      mask_padding_fraction = None,
      return_mask_and_solvent_fraction = None,
      mask_resolution = None,
      out = sys.stdout):
  """Estimate solvent fraction using volume of a low-resolution mask"""
  overall_mask, max_in_sd_map, sd_map = get_overall_mask(map_data = map_data,
    fraction_of_max_mask_threshold = fraction_of_max_mask_threshold,
    use_solvent_content_for_threshold = use_solvent_content_for_threshold,
    mask_padding_fraction = mask_padding_fraction,
    solvent_fraction = solvent_content,  # note name change XXX
    crystal_symmetry = crystal_symmetry,
    resolution = mask_resolution,
    out = out)
  if overall_mask is None:
    if return_mask_and_solvent_fraction:
      return None, None
    else:
      return None

  solvent_fraction = overall_mask.count(False)/overall_mask.size()
  print("Solvent fraction from overall mask: %.3f " %(solvent_fraction), file = out)
  if return_mask_and_solvent_fraction:
    mask_data = map_data.deep_copy()
    mask_data.as_1d().set_selected(overall_mask.as_1d(), 1)
    mask_data.as_1d().set_selected(~overall_mask.as_1d(), 0)
    return mask_data, solvent_fraction
  else:
    return solvent_fraction


def get_solvent_fraction_from_molecular_mass(
        do_not_adjust_dalton_scale = None,
        crystal_symmetry = None, molecular_mass = None, out = sys.stdout):
     """Estimate solvent fraction from molecular mass"""
     map_volume = crystal_symmetry.unit_cell().volume()
     density_factor = 1000*1.23 # just protein density, close enough...
     mm = molecular_mass
     molecule_fraction =  mm*density_factor/map_volume
     if do_not_adjust_dalton_scale or molecule_fraction > 1 and mm > 1000:
        mm = mm/1000  # was in Da

     solvent_fraction = max(0.01, min(1., 1 - (
         mm*density_factor/map_volume)))
     print("Solvent content of %7.2f from molecular mass of %7.1f kDa" %(
     solvent_fraction, mm), file = out)
     return solvent_fraction


def set_up_si(var_dict = None, crystal_symmetry = None,
      is_crystal = None,
      ncs_copies = None, n_residues = None,
      solvent_fraction = None, molecular_mass = None,
      pdb_hierarchy = None, map = None,
      auto_sharpen = True, half_map_data_list = None, verbose = None,
      out = sys.stdout):
    """Set up a sharpening info object (si)"""
    si = sharpening_info(n_real = map.all())
    args = []
    auto_sharpen_methods = var_dict.get('auto_sharpen_methods')
    if auto_sharpen_methods and auto_sharpen_methods !=  ['None'] and \
        len(auto_sharpen_methods) == 1:
      sharpening_method = auto_sharpen_methods[0]
    else:
      sharpening_method = None

    for param in [
       'verbose', 'resolve_size', 'seq_file', 'sequence',
       'box_size',
       'target_n_overlap',
       'restrict_map_size',
       'box_center', 'remove_aniso',
       'input_weight_map_pickle_file', 'output_weight_map_pickle_file',
       'read_sharpened_maps', 'write_sharpened_maps', 'select_sharpened_map',
       'output_directory',
       'smoothing_radius', 'use_local_aniso',
       'local_aniso_in_local_sharpening',
       'overall_before_local',
       'local_sharpening',
       'box_in_auto_sharpen',
       'density_select_in_auto_sharpen',
       'density_select_threshold_in_auto_sharpen',
       'use_weak_density',
       'resolution',
       'd_min_ratio',
       'scale_max',
       'input_d_cut',
       'b_blur_hires',
       'discard_if_worse',
       'mask_atoms', 'mask_atoms_atom_radius', 'value_outside_atoms',
       'soft_mask',
       'tol_r', 'abs_tol_t',
       'rel_tol_t',
       'require_helical_or_point_group_symmetry',
       'max_helical_operators',
       'allow_box_if_b_iso_set',
       'max_box_fraction',
       'cc_cut',
       'max_cc_for_rescale',
       'scale_using_last',
       'density_select_max_box_fraction',
       'k_sharpen',
       'optimize_b_blur_hires',
       'iterate',
       'optimize_d_cut',
        'residual_target', 'sharpening_target',
       'search_b_min', 'search_b_max', 'search_b_n', 'adjust_region_weight',
       'region_weight_method',
        'region_weight_factor',
        'region_weight_buffer',
        'region_weight_default',
        'target_b_iso_ratio',
        'signal_min',
        'buffer_radius',
        'wang_radius',
        'pseudo_likelihood',
        'target_b_iso_model_scale',
       'b_iso', 'b_sharpen',
       'resolution_dependent_b',
       'normalize_amplitudes_in_resdep',
       'region_weight',
       'sa_percent',
       'n_bins',
       'eps',
       'max_regions_to_test',
       'regions_to_keep',
       'fraction_occupied',
       'rmsd',
       'rmsd_resolution_factor',
       'k_sol',
       'b_sol',
       'fraction_complete',
       'nproc',
       'multiprocessing',
       'queue_run_command',
       'verbose',
         ]:
     x = var_dict.get(param)
     if x is not None:
       if type(x) == type([1, 2, 3]):
         xx = []
         for k in x:
           xx.append(str(k))
         args.append("%s = %s" %(param, " ".join(xx)))
       else:
         args.append("%s = %s" %(param, x))
    local_params = get_params_from_args(args)

    # Set solvent content from molecular_mass if present
    if molecular_mass and not solvent_fraction:
       solvent_fraction = get_solvent_fraction_from_molecular_mass(
        crystal_symmetry = crystal_symmetry, molecular_mass = molecular_mass,
        out = out)

    # Test to see if we can use adjusted_path_length as target
    if local_params.map_modification.sharpening_target == \
            'adjusted_path_length':
      print(
       "Checking to make sure we can use 'adjusted_path_length'as target...",
          end = ' ', file = out)
      try:
        from phenix.autosol.trace_and_build import trace_and_build
      except Exception as e:
        raise Sorry("Please set sharpening target to something other than "+
          "adjusted_path_length (not available)")
      print("OK", file = out)



    if (local_params.input_files.seq_file or
       local_params.crystal_info.sequence) and \
        not local_params.crystal_info.solvent_content and \
        not solvent_fraction: # 2017-12-19
        solvent_fraction = get_solvent_fraction(local_params,
          crystal_symmetry = crystal_symmetry,
          ncs_copies = ncs_copies, out = out)
    si.update_with_params(params = local_params,
      crystal_symmetry = crystal_symmetry,
      is_crystal = is_crystal,
      solvent_fraction = solvent_fraction,
      ncs_copies = ncs_copies,
      n_residues = n_residues,
      auto_sharpen = auto_sharpen,
      sharpening_method = sharpening_method,
      pdb_hierarchy = pdb_hierarchy,
      half_map_data_list = half_map_data_list,
      )
    return si

def bounds_to_frac(b, map_data):
  """Convert bounds (grid units) to fractional bounds"""
  a = map_data.all()
  return b[0]/a[0], b[1]/a[1], b[2]/a[2]

def bounds_to_cart(b, map_data, crystal_symmetry = None):
  """Convert bounds (grid units) to cartesian bounds"""
  bb = bounds_to_frac(b, map_data)
  a, b, c = crystal_symmetry.unit_cell().parameters()[:3]

  return a*bb[0], b*bb[1], c*bb[2]

def select_inside_box(lower_bounds = None, upper_bounds = None, xrs = None,
     hierarchy = None):
  """Select atoms inside lower_bounds to upper_bounds from xray_structure
    object xrs"""
  if not hierarchy or not xrs:
    return None
  selection = flex.bool(xrs.scatterers().size())
  for atom_group in hierarchy.atom_groups():
    for atom in atom_group.atoms():
      if atom.xyz[0]>= lower_bounds[0] and \
         atom.xyz[0]<= upper_bounds[0] and \
         atom.xyz[1]>= lower_bounds[1] and \
         atom.xyz[1]<= upper_bounds[1] and \
         atom.xyz[2]>= lower_bounds[2] and \
         atom.xyz[2]<= upper_bounds[2]:
        selection[atom.i_seq] = True

  asc1 = hierarchy.atom_selection_cache()
  return hierarchy.select(selection)

def make_empty_map(template_map = None, value = 0.):
    """Create empty map original_map_in_box"""
    empty_map = flex.double(template_map.as_1d().as_double().size(), value)
    empty_map.reshape(flex.grid(template_map.all()))
    return empty_map

def sum_box_data(starting_map = None, box_map = None,
      lower_bounds = None, upper_bounds = None):
    """sum box data into starting_map"""


    #Pull out current starting_map data
    starting_box_data =  maptbx.copy(starting_map, tuple(lower_bounds), tuple(upper_bounds))
    assert starting_box_data.all() == box_map.all()

    # Add to box_map
    starting_box_data = starting_box_data.as_1d()
    box_map_data = box_map.as_1d()
    starting_box_data+= box_map_data

    # put back into shape
    starting_box_data.reshape(flex.grid(box_map.all()))

    maptbx.set_box(
        map_data_from = starting_box_data,
        map_data_to   = starting_map,
        start         = lower_bounds,
        end           = upper_bounds)
    return starting_map


def copy_box_data(starting_map = None, box_map = None,
      lower_bounds = None, upper_bounds = None):
    """Copy box data into original_map_in_box"""
    maptbx.set_box(
        map_data_from = box_map,
        map_data_to   = starting_map,
        start         = lower_bounds,
        end           = upper_bounds)
    return starting_map

def select_box_map_data(si = None,
           map_data = None,
           first_half_map_data = None,
           second_half_map_data = None,
           pdb_hierarchy = None,
           get_solvent_fraction = True, # XXX test not doing this...
           n_min = 30, # at least 30 atoms to run model sharpening
           restrict_map_size = None,
           out = sys.stdout, local_out = sys.stdout):

  """Select data corresponding to a boxed map"""
  box_solvent_fraction = None
  solvent_fraction = si.solvent_fraction
  crystal_symmetry = si.crystal_symmetry
  box_size = si.box_size
  lower_bounds = None
  upper_bounds = None
  smoothed_box_mask_data = None
  original_box_map_data = None #
  n_buffer = None

  if  (not pdb_hierarchy and not si.box_in_auto_sharpen) and (
      first_half_map_data and second_half_map_data):
      print("Creating density-based soft mask and applying to half-map data", file = out)

      if not si.soft_mask:
        raise Sorry(
         "Need to set soft_mask = True for half-map sharpening without model")
      # NOTE: could precede this by density_select on map_data, save bounds and
      #   cut out half-maps with those bounds. That case could
      #   cover si.soft_mask = False

      half_map_data_list = [first_half_map_data, second_half_map_data]
      box_mask_data, box_map, half_map_data_list, \
       box_solvent_fraction, smoothed_box_mask_data, original_box_map_data = \
       get_and_apply_soft_mask_to_maps(
        resolution = si.resolution,
        wang_radius = si.wang_radius,
        buffer_radius = si.buffer_radius,
        map_data = map_data, crystal_symmetry = crystal_symmetry,
        half_map_data_list = half_map_data_list,
        out = out)

      box_first_half_map, box_second_half_map = half_map_data_list
      box_crystal_symmetry = crystal_symmetry
      box_pdb_hierarchy = None

  elif pdb_hierarchy or (
     si.density_select_in_auto_sharpen and not si.box_in_auto_sharpen):

  #   use map_box for pdb_hierarchy (mask with model)
  #   also use map_box for density_select_in_auto_sharpen sharpening because
  #     need to use the same density select for all 3 maps.


    assert not si.local_sharpening

    if pdb_hierarchy:
      print("Using map_box based on input model", file = out)
      max_box_fraction = si.max_box_fraction
      si.density_select_in_auto_sharpen = False
    else:
      #print >>out, "Using density_select in map_box"
      assert si.density_select_in_auto_sharpen
      max_box_fraction = si.density_select_max_box_fraction

    #----------------------trimming model-------------------------------
    if si.box_center:  # Have model but center at box_center and trim hierarchy
      lower_bounds, upper_bounds = box_from_center(si = si,
        map_data = map_data, out = out)
      if si.soft_mask:
        n_buffer = get_grid_units(map_data = map_data,
          crystal_symmetry = crystal_symmetry,
          radius = si.resolution, out = out)
        n_buffer = int(0.5+n_buffer*1.5)
      else:
        n_buffer = 0
      lower_bounds, upper_bounds = put_bounds_in_range(
       lower_bounds = lower_bounds, upper_bounds = upper_bounds,
       box_size = box_size, n_buffer = n_buffer,
       n_real = map_data.all(), out = out)
      lower_frac = bounds_to_frac(lower_bounds, map_data)
      upper_frac = bounds_to_frac(upper_bounds, map_data)
      lower_cart = bounds_to_cart(
         lower_bounds, map_data, crystal_symmetry = crystal_symmetry)
      upper_cart = bounds_to_cart(
        upper_bounds, map_data, crystal_symmetry = crystal_symmetry)

      if hierarchy: # trimming hierarchy to box and then using trimmed
                    # hierarchy in map_box to create actual box
        xrs = hierarchy.extract_xray_structure(
          crystal_symmetry = si.crystal_symmetry)
        # find everything in box
        sel_hierarchy = select_inside_box(lower_bounds = lower_cart,
         upper_bounds = upper_cart, xrs = xrs, hierarchy = hierarchy)
        n = sel_hierarchy.overall_counts().n_atoms
        print("Selected atoms inside box: %d" %(n), file = out)
        if n<n_min:
          print("Skipping...using entire structure", file = out)
        else:
          hierarchy = sel_hierarchy
    #----------------------end trimming model-------------------------------


    from mmtbx.command_line.map_box import run as run_map_box
    args = ["keep_input_unit_cell_and_grid = False"]
    if si.density_select_in_auto_sharpen:
      args.append('density_select = True')
      #print >>out, "Using density_select in map_box"
      if si.density_select_threshold_in_auto_sharpen is not None:
        args.append('density_select_threshold = %s' %(
          si.density_select_threshold_in_auto_sharpen))
    elif si.box_in_auto_sharpen and not si.mask_atoms:
      print("Using map_box with model", file = out)
    elif si.mask_atoms:
      print("Using map_box with model and mask_atoms", file = out)
      args.append('mask_atoms = True')
      if si.mask_atoms_atom_radius:
        args.append('mask_atoms_atom_radius = %s' %(si.mask_atoms_atom_radius))
      if si.value_outside_atoms:
        args.append('value_outside_atoms = %s' %(si.value_outside_atoms))
      if si.soft_mask:
        print("Using soft mask", file = out)
        args.append('soft_mask = %s' %(si.soft_mask))
        args.append('soft_mask_radius = %s' %(si.resolution))
    else:
      raise Sorry("Unknown choice in select_box_data")

    if restrict_map_size:
      args.append('restrict_map_size = True')
    print("Getting map as box now", file = out)
    local_hierarchy = None
    if pdb_hierarchy:
       local_hierarchy = pdb_hierarchy.deep_copy()
       # run_map_box modifies its argument
    assert isinstance(si.wrapping, bool) # wrapping must be defined
    box = run_map_box(args,
        map_data = map_data, pdb_hierarchy = local_hierarchy,
       write_output_files = False,
       wrapping = si.wrapping,
       crystal_symmetry = crystal_symmetry, log = out)

    lower_bounds = box.gridding_first
    upper_bounds = box.gridding_last

    box_map = box.map_box
    box_map = scale_map(box_map, out = out)
    box_crystal_symmetry = box.box_crystal_symmetry
    if box.hierarchy:
      box_pdb_hierarchy = box.hierarchy
    else:
      box_pdb_hierarchy = None
    if first_half_map_data:
      print("Getting first map as box", file = out)
      if pdb_hierarchy:
        local_hierarchy = pdb_hierarchy.deep_copy() # required
      box_first = run_map_box(args,
        map_data = first_half_map_data, pdb_hierarchy = local_hierarchy,
       write_output_files = False,
       lower_bounds = lower_bounds,
       upper_bounds = upper_bounds,
       crystal_symmetry = crystal_symmetry,
       wrapping = si.wrapping,
       log = out)
      box_first_half_map = box_first.map_box.as_double()
    else:
      box_first_half_map = None

    if second_half_map_data:
      print("Getting second map as box", file = out)
      if pdb_hierarchy:
        local_hierarchy = pdb_hierarchy.deep_copy() # required
      box_second = run_map_box(args,
        map_data = second_half_map_data, pdb_hierarchy = local_hierarchy,
       write_output_files = False,
       lower_bounds = lower_bounds,
       upper_bounds = upper_bounds,
       crystal_symmetry = crystal_symmetry,
       wrapping = si.wrapping,
       log = out)
      box_second_half_map = box_second.map_box.as_double()
    else:
      box_second_half_map = None

  else: # cut out box based on box_center or regions

    if si.box_center:  # center at box_center
      print("Cutting out box based on center at (%.2f, %.2f, %.2f) " %( si.box_center), file = out)
      lower_bounds, upper_bounds = box_from_center(si = si,
        map_data = map_data, out = out)
    elif si.use_weak_density:
      print("Cutting out box based on centering on weak density", file = out)
      lower_bounds, upper_bounds = box_of_smallest_region(si = si,
           map_data = map_data,
           out = out)
    else:
      print("Cutting out box based on centering on strong density", file = out)
      lower_bounds, upper_bounds = box_of_biggest_region(si = si,
           map_data = map_data,
           out = out)
    if si.soft_mask:
      n_buffer = get_grid_units(map_data = map_data,
        crystal_symmetry = crystal_symmetry,
        radius = si.resolution, out = out)
      n_buffer = int(0.5+n_buffer*1.5)
    else:
      n_buffer = 0
    lower_bounds, upper_bounds = put_bounds_in_range(
       lower_bounds = lower_bounds, upper_bounds = upper_bounds,
       box_size = box_size, n_buffer = n_buffer,
       n_real = map_data.all(), out = out)

    # select map data inside this box
    print("\nSelecting map data inside box", file = out)

    box_map, box_crystal_symmetry, \
         smoothed_box_mask_data, original_box_map_data = cut_out_map(
       map_data = map_data.as_double(),
       crystal_symmetry = crystal_symmetry,
       soft_mask = si.soft_mask,
       soft_mask_radius = si.resolution,
       resolution = si.resolution,
       shift_origin = True,
       min_point = lower_bounds, max_point = upper_bounds, out = out)
    box_pdb_hierarchy = None

    if first_half_map_data:
      box_first_half_map, box_first_crystal_symmetry, \
        dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
       map_data = first_half_map_data.as_double(),
       crystal_symmetry = crystal_symmetry,
       soft_mask = si.soft_mask,
       soft_mask_radius = si.resolution,
       resolution = si.resolution,
       shift_origin = True,
       min_point = lower_bounds, max_point = upper_bounds, out = local_out)
    else:
      box_first_half_map = None

    if second_half_map_data:
      box_second_half_map, box_second_crystal_symmetry, \
         dummy_smoothed_box_mask_data, dummy_original_box_map_data = cut_out_map(
       map_data = second_half_map_data.as_double(),
       crystal_symmetry = crystal_symmetry,
       soft_mask = si.soft_mask,
       soft_mask_radius = si.resolution,
       resolution = si.resolution,
       shift_origin = True,
       min_point = lower_bounds, max_point = upper_bounds, out = local_out)
    else:
      box_second_half_map = None

  if not box_map or (
       (not pdb_hierarchy and not second_half_map_data) and \
      box_map.size() > si.max_box_fraction* map_data.size()):

    return None, map_data, first_half_map_data, \
        second_half_map_data, crystal_symmetry, None, \
        smoothed_box_mask_data, None, None # no point

  # figure out solvent fraction in this box...

  if get_solvent_fraction: #
    if box_solvent_fraction is None:
      box_solvent_fraction = get_iterated_solvent_fraction(
        crystal_symmetry = box_crystal_symmetry,
        map = box_map,
        fraction_of_max_mask_threshold = si.fraction_of_max_mask_threshold,
        mask_resolution = si.resolution,
        out = out)
    if box_solvent_fraction is None:
      box_solvent_fraction = si.solvent_fraction
      print("Using overall solvent fraction for box", file = out)
    print("Local solvent fraction: %7.2f" %(box_solvent_fraction), file = out)
  else:
    box_solvent_fraction = None

  box_sharpening_info_obj = box_sharpening_info(
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    n_real = box_map.all(),
    scale_max = si.scale_max,
    wrapping = False,
    crystal_symmetry = box_crystal_symmetry,
    solvent_fraction = box_solvent_fraction)
  return box_pdb_hierarchy, box_map, box_first_half_map, box_second_half_map, \
      box_crystal_symmetry, box_sharpening_info_obj, \
      smoothed_box_mask_data, original_box_map_data, n_buffer

def inside_zero_one(xyz):
  """Move fractional xyz to center it at fractional coordinates (1/2,1/2,1/2)"""
  from scitbx.matrix import col
  offset = xyz-col((0.5, 0.5, 0.5))
  lower_int = offset.iround().as_vec3_double()
  return xyz-lower_int

def move_xyz_inside_cell(xyz_cart = None, crystal_symmetry = None):
  """Try to move coordinates xyz_cart inside a box represented by the
  unit_cell of crystal_symmetry"""
  xyz_local = flex.vec3_double()
  if type(xyz_cart) == type(xyz_local):
    xyz_local = xyz_cart
    is_single = False
  else:
    is_single = True
    xyz_local.append(xyz_cart)

  xyz_frac = crystal_symmetry.unit_cell().fractionalize(xyz_local)
  new_xyz_frac = inside_zero_one(xyz_frac)
  new_xyz_cart = crystal_symmetry.unit_cell().orthogonalize(new_xyz_frac)
  if is_single:
    return new_xyz_cart[0]
  else:
    return new_xyz_cart

def box_from_center( si = None,
           map_data = None,
           out = sys.stdout):
    """Calculate lower, upper bounds of box at si.box_center.
    NOTE: lower and upper are identical always and mark the center"""
    cx, cy, cz = si.crystal_symmetry.unit_cell().fractionalize(si.box_center)
    if cx<0 or cx>1 or cy<0 or cy>1 or cz<0 or cz>1:
       print("Moving box center inside (0, 1)", file = out)
       si.box_center = move_xyz_inside_cell(
           xyz_cart = si.box_center, crystal_symmetry = si.crystal_symmetry)
    cx, cy, cz = si.crystal_symmetry.unit_cell().fractionalize(si.box_center)
    print("\nBox centered at (%7.2f, %7.2f, %7.2f) A" %(
      tuple(si.box_center)), file = out)

    ax, ay, az = map_data.all()
    cgx, cgy, cgz = int(0.5+ax*cx), int(0.5+ay*cy), int(0.5+az*cz),
    print("Box grid centered at (%d, %d, %d)\n" %(cgx, cgy, cgz), file = out)
    return (cgx, cgy, cgz), (cgx, cgy, cgz)

def box_of_smallest_region(si = None,
           map_data = None,
           return_as_list = None,
           out = sys.stdout):
  """Return lower and upper bounds representing the smallest region
  in the si object.
  If return_as_list, return instead a list of fractional centers of all the
  regions"""

  return box_of_biggest_region(si = si,
           map_data = map_data,
           return_as_list = return_as_list,
           use_smallest = True,
           out = out)

def box_of_biggest_region(si = None,
           map_data = None,
           return_as_list = None,
           use_smallest = False,
           out = sys.stdout):
    n_residues = si.n_residues
    ncs_copies = si.ncs_copies
    solvent_fraction = si.solvent_fraction
    """Return lower and upper bounds representing the biggest region
    in the si object.
    If return_as_list, return instead a list of fractional centers of all the
     regions.
    If use_smallest, return instead values for smallest region"""

    b_vs_region = b_vs_region_info()
    co, sorted_by_volume, min_b, max_b, unique_expected_regions, best_score, \
       new_threshold, starting_density_threshold = \
        get_connectivity(
           b_vs_region = b_vs_region,
           map_data = map_data,
           n_residues = n_residues,
           ncs_copies = ncs_copies,
           solvent_fraction = solvent_fraction,
           min_volume = si.min_volume,
           min_ratio = si.min_ratio,
           fraction_occupied = si.fraction_occupied,
           wrapping = si.wrapping,
           residues_per_region = si.residues_per_region,
           max_ratio_to_target = si.max_ratio_to_target,
           min_ratio_to_target = si.min_ratio_to_target,
           min_ratio_of_ncs_copy_to_first = si.min_ratio_of_ncs_copy_to_first,
           starting_density_threshold = si.starting_density_threshold,
           density_threshold = si.density_threshold,
           crystal_symmetry = si.crystal_symmetry,
           chain_type = si.chain_type,
           verbose = si.verbose,
           out = out)

    if len(sorted_by_volume)<2:
      return # nothing to do

    if use_smallest:
      small_ratio = 0.25
      maximum_position_ratio = 0.75
      v1, i1 = sorted_by_volume[1]
      v_small = small_ratio*v1
      maximum_position_small = maximum_position_ratio*(len(sorted_by_volume)-1)+1

      best_pos = 1
      ii = 0
      for v, i in sorted_by_volume[1:]:
        ii+= 1
        if v < v_small: continue
        if ii > maximum_position_small: continue
        best_pos = ii

      v, i = sorted_by_volume[best_pos]
      print("\nVolume of target region %d: %d grid points: "%(best_pos, v), file = out)
    else: # usual
      v, i = sorted_by_volume[1]
      print("\nVolume of largest region: %d grid points: "%(v), file = out)

    print("Region %3d (%3d)  volume:%5d  X:%6d - %6d   Y:%6d - %6d  Z:%6d - %6d "%(
     1, i, v,
     min_b[i][0], max_b[i][0],
     min_b[i][1], max_b[i][1],
     min_b[i][2], max_b[i][2]), file = out)

    if (not return_as_list):
      return min_b[i], max_b[i]

    else: # return a list of centers
      centers_frac = flex.vec3_double()
      a1, a2, a3 = map_data.all()

      for v, i in sorted_by_volume[1:]:
        centers_frac.append(
          tuple((
          (min_b[i][0]+max_b[i][0])/(2.*a1),
          (min_b[i][1]+max_b[i][1])/(2.*a2),
          (min_b[i][2]+max_b[i][2])/(2.*a3),
               ))
                      )
      return centers_frac

def get_fft_map(n_real = None, map_coeffs = None):
    """Calculate an fft_map from map_coeffs with gridding n_real"""
    from cctbx import maptbx
    from cctbx.maptbx import crystal_gridding
    if n_real:
      cg = crystal_gridding(
        unit_cell = map_coeffs.crystal_symmetry().unit_cell(),
        space_group_info = map_coeffs.crystal_symmetry().space_group_info(),
        pre_determined_n_real = n_real)
    else:
      cg = None
    ccs = map_coeffs.crystal_symmetry()
    fft_map = map_coeffs.fft_map( resolution_factor = 0.25,
       crystal_gridding = cg,
       symmetry_flags = maptbx.use_space_group_symmetry)
    fft_map.apply_sigma_scaling()
    return fft_map

def average_from_bounds(lower, upper, grid_all = None):
  """Return average of lower and upper, as fraction of grid_all"""
  avg = []
  for u, l in zip(upper, lower):
    avg.append(0.5*(u+l))
  if grid_all:
     avg_fract = []
     for a, g in zip(avg, grid_all):
       avg_fract.append(a/g)
     avg = avg_fract
  return avg

def get_ncs_copies(site_cart, ncs_object = None,
   only_inside_box = None, unit_cell = None):

  """ Get NCS copies of site_cart, optionally only including those inside
   unit cell box"""
  ncs_group = ncs_object.ncs_groups()[0]
  from scitbx.array_family import flex
  from scitbx.matrix import col
  sites_cart_ncs = flex.vec3_double()

  for t, r in zip(ncs_group.translations_orth_inv(),
                 ncs_group.rota_matrices_inv()):

    sites_cart_ncs.append(r * col(site_cart)  + t)

  if only_inside_box:
    assert unit_cell is not None
    sites_frac_ncs = unit_cell.fractionalize(sites_cart_ncs)
    new_sites_frac = flex.vec3_double()
    for x in sites_frac_ncs:
      if  x[0]>= 0 and x[0]<= 1  and \
          x[1]>= 0 and x[1]<= 1  and \
          x[2]>= 0 and x[2]<= 1:
        new_sites_frac.append(x)
    sites_cart_ncs = unit_cell.orthogonalize(new_sites_frac)
  return sites_cart_ncs


def fit_bounds_inside_box(lower, upper, box_size = None, all = None):
  """Adjust bounds so upper>lower and box size is at least box_size"""
  new_lower = []
  new_upper = []
  if box_size:
   for u, l, s, a in zip(upper, lower, box_size, all):
     ss = u-l+1
     delta = int((1+s-ss)/2) # desired increase in size, to subtract from l
     l = max(0, l-delta)
     ss = u-l+1
     delta = (s-ss) # desired increase in size, to add to u
     u = min(a-1, u+delta)
     l = max(0, l)
     new_lower.append(l)
     new_upper.append(u)
  else:
   for u, l, a in zip(upper, lower, all):
     u = min(a-1, u)
     l = max(0, l)
     new_lower.append(l)
     new_upper.append(u)

  return new_lower, new_upper

def split_boxes(lower = None, upper = None, target_size = None, target_n_overlap = None,
     fix_target_size_and_overlap = None):
  """Split the region from lower to upper into overlapping
  boxes of size about target_size
  NOTE: does not actually use target_n_overlap unless
  fix_target_size_and_overlap is set"""

  all_box_list = []
  for l, u, ts in zip (lower, upper, target_size):
    n = u+1-l
    # offset defined by: ts-offset + ts-offset+...+ts  = n
    #                 ts*n_box-offset*(n_box-1) = n approx
    #             n_box (ts-offset) +offset = n
    #             n_box =  (n-offset)/(ts-offset)
    assert ts>target_n_overlap
    if fix_target_size_and_overlap:
      n_box = (n-target_n_overlap-1)/(ts-target_n_overlap)
      if n_box>int(n_box):
        n_box = int(n_box)+1
      else:
        n_box = int(n_box)

      new_target_size = ts
      offset = ts-target_n_overlap
    else: # original version
      n_box = max(1, (n+(3*ts//4))//ts)
      new_target_size = int(0.9+n/n_box)
      offset = int(0.9+(n-new_target_size)/max(1, n_box-1))
    box_list = []
    for i in range(n_box):
      start_pos = max(l, l+i*offset)
      end_pos = min(u, start_pos+new_target_size)
      if fix_target_size_and_overlap:
        start_pos = max(0, min(start_pos, end_pos-new_target_size))
      box_list.append([start_pos, end_pos])
    all_box_list.append(box_list)
  new_lower_upper_list = []
  for xs, xe in all_box_list[0]:
    for ys, ye in all_box_list[1]:
      for zs, ze in all_box_list[2]:
        new_lower_upper_list.append([(xs, ys, zs, ), (xe, ye, ze, )])
  return new_lower_upper_list

def get_target_boxes(si = None, ncs_obj = None, map = None,
    pdb_hierarchy = None, out = sys.stdout):

  """Identify suitable parts of map to sharpen by segmentation of map.
  Represent parts to sharpen as boxes with lower, uppper bounds.
  Return: upper_bounds_list, lower_bounds_list,
     centers_cart_ncs_list, centers_cart, all_cart
  """
  print(80*"-", file = out)
  print("Getting segmented map to ID locations for sharpening", file = out)
  print(80*"-", file = out)

  if si.input_weight_map_pickle_file:
    from libtbx import easy_pickle
    file_name = si.input_weight_map_pickle_file
    print("Loading segmentation data from %s" %(file_name), file = out)
    tracking_data = easy_pickle.load(file_name)

  else:
    args = [
        'resolution = %s' %(si.resolution),
        'seq_file = %s' %(si.seq_file),
        'sequence = %s' %(si.sequence),
        'solvent_content = %s' %(si.solvent_fraction),
        'auto_sharpen = False', # XXX could sharpen overall
        'write_output_maps = True',
        'add_neighbors = False',
        'density_select = False', ]
    if si.is_crystal:
      args.append("is_crystal = True")
    ncs_group_obj, remainder_ncs_group_obj, tracking_data = run(
     args,
     map_data = map.deep_copy(),
     ncs_obj = ncs_obj,
     crystal_symmetry = si.crystal_symmetry)

  if si.output_weight_map_pickle_file:
    from libtbx import easy_pickle
    file_name = os.path.join(si.output_directory, si.output_weight_map_pickle_file)
    print("Dumping segmentation data to %s" %(file_name), file = out)
    easy_pickle.dump(file_name, tracking_data)

  if not ncs_obj or ncs_obj.max_operators() == 0:
    from mmtbx.ncs.ncs import ncs
    ncs_obj = ncs()
    ncs_obj.set_unit_ncs()

  print("Regions in this map:", file = out)
  centers_frac = flex.vec3_double()
  upper_bounds_list = []
  lower_bounds_list = []
  if pdb_hierarchy and pdb_hierarchy.atoms().extract_xyz().size()>1:
    xyz_list = pdb_hierarchy.atoms().extract_xyz()
    i_end = xyz_list.size()
    n_centers = min(i_end, max(1, len(tracking_data.output_region_map_info_list)))
    n_steps = min(n_centers, xyz_list.size())
    i_step = int(0.5+min(i_end/2, i_end/n_steps)) # about n_centers but up to n_atoms
    i_start = max(1, int(0.5+i_step/2))
    from scitbx.matrix import col
    ma = map.all()
    for i in range(i_start, i_end, i_step):
      lower_cart = col(xyz_list[i])
      lower_frac = si.crystal_symmetry.unit_cell().fractionalize(lower_cart)
      lower = [
        int(0.5+ma[0]*lower_frac[0]),
        int(0.5+ma[1]*lower_frac[1]),
        int(0.5+ma[2]*lower_frac[2])]
      lower, upper = fit_bounds_inside_box(
        lower, lower, box_size = si.box_size, all = map.all())
      upper_bounds_list.append(upper)
      lower_bounds_list.append(lower)
      average_fract = average_from_bounds(lower, upper, grid_all = map.all())
      centers_frac.append(average_fract)
  else:
    for map_info_obj in tracking_data.output_region_map_info_list:
      lower, upper = map_info_obj.lower_upper_bounds()
      lower, upper = fit_bounds_inside_box(
        lower, upper,
        box_size = None, # take the whole region, not just center
        all = map.all())
      for lower, upper in split_boxes(lower = lower, upper = upper,
         target_size = si.box_size,
         target_n_overlap = si.target_n_overlap):
         upper_bounds_list.append(upper)
         lower_bounds_list.append(lower)
         average_fract = average_from_bounds(lower, upper, grid_all = map.all())
         centers_frac.append(average_fract)


  centers_cart = si.crystal_symmetry.unit_cell().orthogonalize(centers_frac)


  #  Make ncs-related centers
  print("NCS ops:", ncs_obj.max_operators(), file = out)
  centers_cart_ncs_list = []
  for i in range(centers_cart.size()):
    centers_cart_ncs_list.append(get_ncs_copies(
       centers_cart[i], ncs_object = ncs_obj, only_inside_box = True,
       unit_cell = si.crystal_symmetry.unit_cell()) )

  all_cart = flex.vec3_double()
  for center_list in centers_cart_ncs_list:
    all_cart.extend(center_list)

  sharpening_centers_file = os.path.join(
      si.output_directory, "sharpening_centers.pdb")
  write_atoms(file_name = sharpening_centers_file, # PDB OK
    crystal_symmetry = si.crystal_symmetry, sites = centers_cart)
  ncs_sharpening_centers_file = os.path.join(
      si.output_directory, "ncs_sharpening_centers.pdb")
  write_atoms(file_name = ncs_sharpening_centers_file, # PDB OK
    crystal_symmetry = si.crystal_symmetry, sites = all_cart)

  print("\nSharpening centers (matching shifted_map_file).\n\n "+\
      "Written to: \n%s \n%s\n"%(
      sharpening_centers_file, ncs_sharpening_centers_file), file = out)

  for i in range(centers_cart.size()):
    print("Center: %s (%7.2f, %7.2f, %7.2f)  Bounds: %s :: %s " %(
        i, centers_cart[i][0], centers_cart[i][1], centers_cart[i][2],
          str(lower_bounds_list[i]), str(upper_bounds_list[i])), file = out)

  print(80*"-", file = out)
  print("Done getting segmented map to ID locations for sharpening", file = out)
  print(80*"-", file = out)


  return upper_bounds_list, lower_bounds_list, \
     centers_cart_ncs_list, centers_cart, all_cart

def get_box_size(lower_bound = None, upper_bound = None):
  """Get box size from lower and upper bounds"""
  box_size = []
  for lb, ub in zip(lower_bound, upper_bound):
    box_size.append(ub-lb+1)
  return box_size

def mean_dist_to_nearest_neighbor(all_cart):
  """Get mean distance to nearest-neighbors for a set of points in all_cart"""
  if all_cart.size()<2:  # nothing to check
    return None
  sum_dist = 0.
  sum_n = 0.
  for i in range(all_cart.size()):
    xyz = all_cart[i:i+1]
    others = all_cart[:i]
    others.extend(all_cart[i+1:])
    sum_dist+= get_closest_dist(xyz, others)
    sum_n+= 1.
  return sum_dist/max(1., sum_n)


def run_local_sharpening(si = None,
    auto_sharpen_methods = None,
    map = None,
    ncs_obj = None,
    half_map_data_list = None,
    pdb_hierarchy = None,
    out = sys.stdout):
  """Run local sharpening.
  Run auto_sharpen_map_or_map_coeffs with box_in_auto_sharpen = True and
  centered at different places.  Identify the places as centers of regions.
  Run on au of NCS and apply NCS to get remaining positions
  """

  print(80*"-", file = out)
  print("Running local sharpening", file = out)
  print(80*"-", file = out)

  if si.overall_before_local:
    # first do overall sharpening of the map to get it about right
    print(80*"*", file = out)
    print("\nSharpening map overall before carrying out local sharpening\n", file = out)
    overall_si = deepcopy(si)
    overall_si.local_sharpening = False  # don't local sharpen
    overall_si = auto_sharpen_map_or_map_coeffs(si = overall_si,
          auto_sharpen_methods = auto_sharpen_methods,
          map = map,
          half_map_data_list = half_map_data_list,
          pdb_hierarchy = pdb_hierarchy,
          out = out)
    sharpened_map = overall_si.map_data
    print("\nDone sharpening map overall before carrying out local sharpening\n", file = out)
    print(80*"*", file = out)
  else:
    sharpened_map = map


  # Accumulate sums
  starting_weight = 0.01 # put starting map everywhere with low weight
  sum_weight_map = make_empty_map(template_map = map, value = starting_weight)

  # in case a pixel is not covered...
  sum_weight_value_map = starting_weight*sharpened_map.deep_copy()

  print("\nUsing overall map for any regions where "+\
     "no local information is present", file = out)

  id_list = []
  b_iso_list = flex.double()
  starting_b_iso_list = flex.double()

  # use sharpened map here
  upper_bounds_list, lower_bounds_list, \
     centers_cart_ncs_list, centers_cart, all_cart = \
     get_target_boxes(si = si, map = sharpened_map, ncs_obj = ncs_obj,
       pdb_hierarchy = pdb_hierarchy, out = out)

  dist = mean_dist_to_nearest_neighbor(all_cart)
  if not dist:
    dist = 10.
    if not si.smoothing_radius:
      print("No nearest neighbors...best to set smoothing radius", file = out)
  print("\nMean distance to nearest center is %7.2f A " %(
    dist), file = out)
  if not si.smoothing_radius:
    si.smoothing_radius = float("%.0f" %(dist*2/3)) # 10A from nearest neighbor
    print("Using %s A for smoothing radius" %(si.smoothing_radius), file = out)

  i = -1
  for ub, lb, centers_ncs_cart, center_cart in zip(
    upper_bounds_list, lower_bounds_list, centers_cart_ncs_list, centers_cart):
    i+= 1
    if si.select_sharpened_map is not None and i !=  si.select_sharpened_map:
      continue
    map_file_name = 'sharpened_map_%s.ccp4' %(i)
    if 0 and si.read_sharpened_maps: # cannot do this as no bounds
      print("\nReading sharpened map directly from %s" %(map_file_name), file = out)
      result = get_map_object(file_name = map_file_name,
        out = out)
      local_map_data = result[0]
    else:

      local_si = deepcopy(si)
      local_si.local_sharpening = False  # don't do it again
      local_si.box_size = get_box_size(lower_bound = lb, upper_bound = ub)
      local_si.box_center = center_cart
      local_si.box_in_auto_sharpen = True
      local_si.density_select_in_auto_sharpen = False
      local_si.use_local_aniso = si.local_aniso_in_local_sharpening
      local_si.remove_aniso = si.local_aniso_in_local_sharpening
      local_si.max_box_fraction = 999 # just bigger than 1
      local_si.density_select_max_box_fraction = 999
      local_si.nproc = 1
      print(80*"+", file = out)
      print("Getting local sharpening for box %s" %(i), file = out)
      print(80*"+", file = out)
      bsi = auto_sharpen_map_or_map_coeffs(si = local_si,
        auto_sharpen_methods = auto_sharpen_methods,
        map = sharpened_map,
        half_map_data_list = half_map_data_list,
        pdb_hierarchy = pdb_hierarchy,
        return_bsi = True, # just return the bsi of sharpened data
        out = out)

      if not bsi or not bsi.map_data:
        print("\nNo result for local map %s ...skipping" %(i), file = out)
        continue

      # merge with background using bsi.smoothed_box_mask_data
      if bsi.smoothed_box_mask_data:
        print("Merging small map into overall map in soft-mask region", file = out)
        bsi.merge_into_overall_map(overall_map = map) # XXX overall_map not used

      # Now remove buffer region
      if bsi.n_buffer: # extract just the good part
         print("Removing buffer from small map", file = out)
         bsi.remove_buffer(out = out)


    weight_data = bsi.get_gaussian_weighting(out = out)
    weighted_data = bsi.map_data*weight_data
    sum_weight_value_map = sum_box_data(starting_map = sum_weight_value_map,
       box_map = weighted_data,
       lower_bounds = bsi.lower_bounds,
       upper_bounds = bsi.upper_bounds)

    sum_weight_map = sum_box_data(starting_map = sum_weight_map,
       box_map = weight_data,
       lower_bounds = bsi.lower_bounds,
       upper_bounds = bsi.upper_bounds)

    id_list.append(i)
    starting_b_iso_list.append(bsi.starting_b_iso)
    b_iso_list.append(bsi.b_iso)

    print(80*"+", file = out)
    print("End of getting local sharpening for small box %s" %(i), file = out)
    print(80*"+", file = out)

  print("\nOverall map created from total of %s local maps" %(i), file = out)
  if si.overall_before_local:
    print("Note: overall map already sharpened with global sharpening", file = out)

  if starting_b_iso_list.size()<1:
    print("No results for local sharpening...", file = out)
  else:
    print("Summary of b_iso values by local map:", file = out)
    print(" ID   Starting B-iso    Sharpened B-iso", file = out)
    for i, starting_b_iso, b_iso in zip(id_list, starting_b_iso_list, b_iso_list):
      print(" %4s    %7.2f        %7.2f" %(i, starting_b_iso, b_iso), file = out)
    print("\nMean    %7.2f        %7.2f" %(
     starting_b_iso_list.min_max_mean().mean,
     b_iso_list.min_max_mean().mean), file = out)

  si.map_data = sum_weight_value_map/sum_weight_map

  # Get overall b_iso...
  print("\nGetting overall b_iso of composite map...", file = out)
  map_coeffs_aa, map_coeffs, f_array, phases = effective_b_iso(
     map_data = si.map_data,
      resolution = si.resolution,
      d_min_ratio = si.d_min_ratio,
      scale_max = si.scale_max,
      crystal_symmetry = si.crystal_symmetry,
      out = out)

  print(80*"+", file = out)
  print("End of getting local sharpening ", file = out)
  print(80*"+", file = out)

  return si

def auto_sharpen_map_or_map_coeffs(
        si = None,
        resolution = None,        # resolution is required
        crystal_symmetry = None,  # supply crystal_symmetry and map or
        map = None,               #  map and n_real
        wrapping = None,
        half_map_data_list = None,     #  two half-maps matching map
        is_crystal = None,
        map_coeffs = None,
        pdb_hierarchy = None,
        ncs_obj = None,
        seq_file = None,
        sequence = None,
        rmsd = None,
        rmsd_resolution_factor = None,
        k_sol = None,
        b_sol = None,
        fraction_complete = None,
        n_real = None,
        solvent_content = None,
        molecular_mass = None,
        region_weight = None,
        sa_percent = None,
        n_bins = None,
        eps = None,
        max_regions_to_test = None,
        regions_to_keep = None,
        fraction_occupied = None,
        input_weight_map_pickle_file = None,
        output_weight_map_pickle_file = None,
        read_sharpened_maps = None,
        write_sharpened_maps = None,
        select_sharpened_map = None,
        output_directory = None,
        smoothing_radius = None,
        local_sharpening = None,
        local_aniso_in_local_sharpening = None,
        overall_before_local = None,
        use_local_aniso = None,
        auto_sharpen = None,
        box_in_auto_sharpen = None, # n_residues, ncs_copies required if not False
        density_select_in_auto_sharpen = None,
        density_select_threshold_in_auto_sharpen = None,
        allow_box_if_b_iso_set = None,
        use_weak_density = None,
        discard_if_worse = None,
        n_residues = None,
        ncs_copies = None,
        box_center = None,
        remove_aniso = None,
        box_size = None,
        target_n_overlap = None,
        lower_bounds = None,
        upper_bounds = None,
        restrict_map_size = None,
        auto_sharpen_methods = None,
        residual_target = None,
        sharpening_target = None,
        d_min_ratio = None,
        scale_max = None,
        input_d_cut = None,
        b_blur_hires = None,
        max_box_fraction = None,
        cc_cut = None,
        max_cc_for_rescale = None,
        scale_using_last = None,
        density_select_max_box_fraction = None,
        mask_atoms = None,
        mask_atoms_atom_radius = None,
        value_outside_atoms = None,
        soft_mask = None,
        tol_r = None,
        abs_tol_t = None,
        rel_tol_t = None,
        require_helical_or_point_group_symmetry = None,
        max_helical_operators = None,
        k_sharpen = None,
        optimize_d_cut = None,
        optimize_b_blur_hires = None,
        iterate = None,
        search_b_min = None,
        search_b_max = None,
        search_b_n = None,
        adjust_region_weight = None,
        region_weight_method = None,
        region_weight_factor = None,
        region_weight_buffer = None,
        region_weight_default = None,
        target_b_iso_ratio = None,
        signal_min = None,
        buffer_radius = None,
        wang_radius = None,
        pseudo_likelihood = None,
        target_b_iso_model_scale = None,
        b_iso = None, # if set, use it
        b_sharpen = None, # if set, use it
        resolution_dependent_b = None, # if set, use it
        normalize_amplitudes_in_resdep = None, # if set, use it
        return_bsi = False,
        verbose = None,
        resolve_size = None,
        nproc = None,
        multiprocessing = None,
        queue_run_command = None,
        out = sys.stdout):

    """Auto-sharpen a map or map coefficients"""

    if si:  #
      resolution = si.resolution
      crystal_symmetry = si.crystal_symmetry
      if not auto_sharpen:
        auto_sharpen = si.auto_sharpen
      if verbose is None:
        verbose = si.verbose
      if resolve_size is None:
        resolve_size = si.resolve_size

    if auto_sharpen is None:
      auto_sharpen = True

    if map_coeffs and not resolution:
       resolution = map_coeffs.d_min()
    if map_coeffs and not crystal_symmetry:
       crystal_symmetry = map_coeffs.crystal_symmetry()

    assert resolution is not None

    if map:
      return_as_map = True
    else:  # convert from structure factors to create map if necessary
      map = get_fft_map(n_real = n_real, map_coeffs = map_coeffs).real_map_unpadded()
      return_as_map = False

    # Set ncs_copies if possible
    if ncs_copies is None and ncs_obj and ncs_obj.max_operators():
      ncs_copies = ncs_obj.max_operators()
      print("Set ncs copies based on ncs_obj to %s" %(ncs_copies), file = out)

    # Determine if we are running model_sharpening
    if half_map_data_list and len(half_map_data_list) == 2:
      if auto_sharpen_methods !=  ['external_map_sharpening']:
        auto_sharpen_methods = ['half_map_sharpening']
    elif pdb_hierarchy:
      auto_sharpen_methods = ['model_sharpening']
    if not si:
      # Copy parameters to si (sharpening_info_object)
      si = set_up_si(var_dict = locals(),
        crystal_symmetry = crystal_symmetry,
        is_crystal = is_crystal,
        solvent_fraction = solvent_content,
        molecular_mass = molecular_mass,
        auto_sharpen = auto_sharpen,
        map = map,
        verbose = verbose,
        half_map_data_list = half_map_data_list,
        pdb_hierarchy = pdb_hierarchy,
        ncs_copies = ncs_copies,
        n_residues = n_residues, out = out)
    if wrapping is not None:
      si.wrapping = wrapping
    # Figure out solvent fraction
    if si.solvent_fraction is None:
      si.solvent_fraction = get_iterated_solvent_fraction(
        crystal_symmetry = crystal_symmetry,
        verbose = si.verbose,
        resolve_size = si.resolve_size,
        fraction_of_max_mask_threshold = si.fraction_of_max_mask_threshold,
        mask_resolution = si.resolution,
        map = map,
        out = out)
    if si.solvent_fraction:
      print("Estimated solvent content: %.2f" %(si.solvent_fraction), file = out)
    else:
      raise Sorry("Unable to estimate solvent content...please supply "+
        "solvent_content \nor molecular_mass")
    # Determine if we are running half-map or model_sharpening
    if half_map_data_list and len(half_map_data_list) == 2:
      first_half_map_data = half_map_data_list[0]
      second_half_map_data = half_map_data_list[1]
    else:
      first_half_map_data = None
      second_half_map_data = None

    # Decide if we are running local sharpening (overlapping set of sharpening
    #   runs at various locations)
    libtbx.call_back(message = 'sharpen', data = None)
    if si.local_sharpening:
      return run_local_sharpening(si = si,
         auto_sharpen_methods = auto_sharpen_methods,
         map = map,
         ncs_obj = ncs_obj,
         half_map_data_list = half_map_data_list,
         pdb_hierarchy = pdb_hierarchy,
         out = out)

    # Get preliminary values of sharpening
    working_map = map  # use another name for map XXX
    if si.iterate and not si.preliminary_sharpening_done:
      si.preliminary_sharpening_done = True
      si.iterate = False
      # first do overall sharpening of the map to get it about right
      print(80*"*", file = out)
      print("\nSharpening map overall before carrying out final sharpening\n", file = out)
      overall_si = deepcopy(si)
      overall_si.local_sharpening = False  # don't local sharpen
      overall_si = auto_sharpen_map_or_map_coeffs(si = overall_si,
            auto_sharpen_methods = auto_sharpen_methods,
            map = map,
            half_map_data_list = half_map_data_list,
            pdb_hierarchy = pdb_hierarchy,
            out = out)
      working_map = overall_si.map_data
      # Get solvent content again
      overall_si.solvent_content = None
      overall_si.solvent_fraction = get_iterated_solvent_fraction(
        crystal_symmetry = crystal_symmetry,
        verbose = overall_si.verbose,
        resolve_size = overall_si.resolve_size,
        fraction_of_max_mask_threshold = si.fraction_of_max_mask_threshold,
        mask_resolution = si.resolution,
        map = working_map,
        out = out)
      print("Resetting solvent fraction to %.2f " %(
         overall_si.solvent_fraction), file = out)
      si.solvent_fraction = overall_si.solvent_fraction

      print("\nDone sharpening map overall before carrying out final sharpening\n", file = out)
      print(80*"*", file = out)
      si.b_blur_hires = 0.  # from now on, don't apply extra blurring
    else:
      working_map = map

    # Now identify optimal sharpening params
    print(80*"=", file = out)
    print("\nRunning auto_sharpen to get sharpening parameters\n", file = out)
    print(80*"=", file = out)
    result = run_auto_sharpen( # get sharpening parameters standard run
      si = si,
      map_data = working_map,
      first_half_map_data = first_half_map_data,
      second_half_map_data = second_half_map_data,
      pdb_hierarchy = pdb_hierarchy,
      auto_sharpen_methods = auto_sharpen_methods,
      print_result = False,
      return_bsi = return_bsi,
      out = out)
    if return_bsi:
      return result  # it is a box_sharpening_info object
    else:
      si = result
    print(80*"=", file = out)
    print("\nDone running auto_sharpen to get sharpening parameters\n", file = out)
    print(80*"=", file = out)

    # Apply the optimal sharpening values and save map in si.map_data
    # First test without sharpening if sharpening_method is b_iso, b and
    # b_iso is not set
    if si.sharpening_method in [
       'b_iso', 'b_iso_to_d_cut', 'resolution_dependent'] and b_iso is None:
      local_si = deepcopy(si)
      local_si.sharpening_method = 'no_sharpening'
      local_si.sharpen_and_score_map(map_data = working_map, out = null_out())
      print("\nScore for no sharpening: %7.2f " %(local_si.score), file = out)
    else:
      local_si = None

    print(80*"=", file = out)
    print("\nApplying final sharpening to entire map", file = out)
    print(80*"=", file = out)
    si.sharpen_and_score_map(map_data = working_map, set_b_iso = True, out = out)
    if si.discard_if_worse and local_si and local_si.score > si.score:
       print("Sharpening did not improve map "+\
        "(%7.2f sharpened, %7.2f unsharpened). Discarding sharpened map" %(
        si.score, local_si.score), file = out)
       print("\nUse discard_if_worse = False to keep the sharpening", file = out)
       local_si.sharpen_and_score_map(map_data = working_map, out = out)
       si = local_si
    if not si.is_model_sharpening() and not si.is_half_map_sharpening():
      si.show_score(out = out)
      si.show_summary(out = out)

    return si  # si.map_data is the sharpened map

def estimate_signal_to_noise(value_list = None, minimum_value_to_include = 0):
  """Estimate signal and noise in a set of values.
  Get "noise" from rms value of value_list(n) compared with
   average of n-2, n-1, n+1, n+2.
  assumes middle is the high part of the very smooth curve.
  Don't include values < minimum_value_to_include.
  """
  mean_square_diff = 0.
  mean_square_diff_n = 0.
  for b2, b1, value, p1, p2 in zip(value_list,
    value_list[1:],
    value_list[2:],
    value_list[3:],
    value_list[4:]):
    too_low = False
    for xx in [b2, b1, value, p1, p2]:
      if xx <minimum_value_to_include:
        too_low = True
    if not too_low:
      mean_square_diff_n+= 1
      mean_square_diff+= ( (b2+b1+p1+p2)*0.25 - value)**2
  rmsd = (mean_square_diff/max(1, mean_square_diff_n))**0.5
  if value_list.size()>0:
    min_adj_sa = max(value_list[0], value_list[-1])
    max_adj_sa = value_list.min_max_mean().max
    signal_to_noise = (max_adj_sa-min_adj_sa)/max(1.e-10, rmsd)
  else:
    signal_to_noise = 0.
  return signal_to_noise

def optimize_b_blur_or_d_cut_or_b_iso(
           optimization_target = 'b_blur_hires',
           local_best_si = None,
           local_best_map_and_b = None,
           si_id_list = None,
           si_score_list = None,
           delta_b = None,
           original_b_iso = None,
           f_array = None,
           phases = None,
           delta_b_blur_hires = 100,
           delta_d_cut = 0.25,
           n_cycle_optimize = 5,
           min_cycles = 2,
           n_range = 5,
           out = sys.stdout):

  """Optimize b_blur or d_cut and b_iso in a map. Score
   with adjusted surface area, kurtosis, or adjusted path length"""

  assert optimization_target in ['b_blur_hires', 'd_cut', 'b_iso']
  if optimization_target == 'b_blur_hires':
    print("\nOptimizing b_blur_hires. ", file = out)
  elif optimization_target == 'd_cut':
    print("\nOptimizing d_cut. ", file = out)
    local_best_si.input_d_cut = local_best_si.get_d_cut()
  elif optimization_target == 'b_iso':
    print("\nOptimizing b_iso. ", file = out)

  local_best_si.show_summary(out = out)

  print("Current best score = %7.3f b_iso = %5.1f  b_blur_hires = %5.1f d_cut = %5.1f" %(
       local_best_si.score, local_best_si.b_iso,
       local_best_si.b_blur_hires,
       local_best_si.get_d_cut()), file = out)


  # existing values:
  value_dict = {}
  for id, score in zip(si_id_list, si_score_list):
    value_dict[id] = score
  best_score = local_best_si.score
  delta_b_iso = delta_b

  local_best_score = best_score
  improving = True
  working_best_si = deepcopy(local_best_si)
  for cycle in range(n_cycle_optimize):
    if not improving: break
    print("Optimization cycle %s" %(cycle), file = out)
    print("Current best score = %7.3f b_iso = %5.1f  b_blur_hires = %5.1f d_cut = %5.1f" %(
     working_best_si.score, working_best_si.b_iso,
        working_best_si.b_blur_hires,
       working_best_si.get_d_cut()), file = out)
    if working_best_si.verbose:
      print(" B-sharpen B-iso B-blur   Adj-SA    "+\
           "Kurtosis  SA-ratio   Regions   d_cut   b_blur_hires", file = out)
    local_best_working_si = deepcopy(working_best_si)
    improving = False
    for jj in range(-n_range, n_range+1):
        if optimization_target == 'b_blur_hires': # try optimizing b_blur_hires
          test_b_blur_hires = max(0., working_best_si.b_blur_hires+jj*delta_b_blur_hires)
          test_d_cut = working_best_si.get_d_cut()
          test_b_iso = working_best_si.b_iso
        elif optimization_target == 'd_cut':
          test_b_blur_hires = working_best_si.b_blur_hires
          test_d_cut = working_best_si.get_d_cut()+jj*delta_d_cut
          test_b_iso = working_best_si.b_iso
        elif optimization_target == 'b_iso':
          test_b_blur_hires = working_best_si.b_blur_hires
          test_d_cut = working_best_si.get_d_cut()
          test_b_iso = working_best_si.b_iso+jj*delta_b_iso

        id = "%.3f_%.3f_%.3f" %(test_b_iso, test_b_blur_hires, test_d_cut)
        if id in value_dict:
          score = value_dict[id]
        else:
          local_si = deepcopy(local_best_si)
          local_f_array = f_array
          local_phases = phases
          local_si.b_blur_hires = test_b_blur_hires
          local_si.input_d_cut = test_d_cut
          local_si.b_iso = test_b_iso
          local_si.b_sharpen = original_b_iso-local_si.b_iso

          local_map_and_b = apply_sharpening(
            f_array = local_f_array, phases = local_phases,
            sharpening_info_obj = local_si,
            crystal_symmetry = local_si.crystal_symmetry,
            out = null_out())
          local_si = score_map(
             map_data = local_map_and_b.map_data, sharpening_info_obj = local_si,
            out = null_out())
          value_dict[id] = local_si.score
          if local_si.verbose:
            print(" %6.1f     %6.1f  %5s   %7.3f  %7.3f" %(
              local_si.b_sharpen, local_si.b_iso,
               local_si.b_blur_hires,
               local_si.adjusted_sa, local_si.kurtosis) + \
              "  %7.1f         %7.3f   %7.3f %7.3f " %(
               zero_if_none(local_si.adjusted_path_length), #local_si.sa_ratio,
               local_si.normalized_regions,
               test_d_cut,
               test_b_blur_hires
                ), file = out)

          if local_si.score > local_best_score:
            local_best_score = local_si.score
            local_best_working_si = deepcopy(local_si)
    if local_best_score > best_score:
      best_score = local_best_score
      working_best_si = deepcopy(local_best_working_si)
      delta_b_iso = delta_b_iso/2
      delta_b_blur_hires = delta_b_blur_hires/2
      delta_d_cut = delta_d_cut/2
      print("Current working best "+\
          "score = %7.3f b_iso = %5.1f  b_blur_hires = %5.1f d_cut = %5.1f" %(
            working_best_si.score, working_best_si.b_iso,
        working_best_si.b_blur_hires,
        working_best_si.get_d_cut()), file = out)
      improving = True

  if working_best_si and working_best_si.score > local_best_si.score:
    print("Using new values of b_iso and b_blur_hires and d_cut", file = out)
    local_best_si = working_best_si

  local_best_si.show_summary(out = out)

  return local_best_si, local_best_map_and_b

def set_mean_sd_of_map(map_data = None, target_mean = None, target_sd = None):
    """Set mean of map to 0, normalize SD to 1"""
    if not map_data: return None
    new_mean = map_data.as_1d().min_max_mean().mean
    new_sd = max(1.e-10, map_data.sample_standard_deviation())
    map_data = (map_data-new_mean)/new_sd # normalized
    return map_data*target_sd + target_mean # restore original

def run_auto_sharpen(
      si = None,
      map_data = None,
      first_half_map_data = None,
      second_half_map_data = None,
      pdb_hierarchy = None,
      auto_sharpen_methods = None,
      print_result = True,
      return_bsi = False,
      out = sys.stdout):

  """Run auto sharpening.  Returns an si sharpening_info object.

  Identifies parameters for optimal map sharpening using analysis of density,
  model-correlation, or half-map correlation (first_half_map_data vs
  vs second_half_map_data).

  NOTE: We can apply this to any map_data (a part or whole of the map)
  BUT: need to update n_real if we change the part of the map!
  change with map data: crystal_symmetry, solvent_fraction, n_real, wrapping,
  """

  if si.verbose:
     local_out = out
  else:
     local_out = null_out()

  smoothed_box_mask_data = None
  original_box_map_data = None
  if si.auto_sharpen and (
    si.box_in_auto_sharpen or
    si.density_select_in_auto_sharpen or pdb_hierarchy):

    original_box_sharpening_info_obj = deepcopy(si)  # should really not be box
    box_pdb_hierarchy, box_map_data, box_first_half_map_data, \
         box_second_half_map_data, \
         box_crystal_symmetry, box_sharpening_info_obj, \
         smoothed_box_mask_data, original_box_map_data, n_buffer = \
       select_box_map_data(si = si,
           map_data = map_data,
           first_half_map_data = first_half_map_data,
           second_half_map_data = second_half_map_data,
           pdb_hierarchy = pdb_hierarchy,
           restrict_map_size = si.restrict_map_size,
           out = out, local_out = local_out)

    if box_sharpening_info_obj is None: # did not do it
      print("Box map is similar in size to entire map..."+\
         "skipping representative box of density", file = out)
      original_box_sharpening_info_obj = None
      crystal_symmetry = si.crystal_symmetry
    else:
      print("Using small map to identify optimal sharpening", file = out)
      print("Box map grid: %d  %d  %d" %(
         box_map_data.all()), file = out)
      print("Box map cell: %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f  "%(
        box_crystal_symmetry.unit_cell().parameters()), file = out)
      original_map_data = map_data
      original_crystal_symmetry = si.crystal_symmetry

      map_data = box_map_data
      pdb_hierarchy = box_pdb_hierarchy
      if si.density_select_in_auto_sharpen and ( # catch empty pdb_hierarchy
         (not pdb_hierarchy) or (not
           pdb_hierarchy.overall_counts().n_residues)):
        pdb_hierarchy = None

      crystal_symmetry = box_crystal_symmetry
      if box_first_half_map_data:
        first_half_map_data = box_first_half_map_data
      if box_second_half_map_data:
        second_half_map_data = box_second_half_map_data
      # SET si for box now...
      si = deepcopy(si).update_with_box_sharpening_info(
        box_sharpening_info_obj = box_sharpening_info_obj)
  else:
    original_box_sharpening_info_obj = None
    box_sharpening_info_obj = None
    crystal_symmetry = si.crystal_symmetry

  starting_mean = map_data.as_1d().min_max_mean().mean
  starting_sd = map_data.sample_standard_deviation()

  print("\nGetting original b_iso...", file = out)
  map_coeffs_aa, map_coeffs, f_array, phases = effective_b_iso(
     map_data = map_data,
      resolution = si.resolution,
      d_min_ratio = si.d_min_ratio,
      scale_max = si.scale_max,
      remove_aniso = si.remove_aniso,
      crystal_symmetry = si.crystal_symmetry,
      out = out)
  original_b_iso = map_coeffs_aa.b_iso
  if original_b_iso is None:
    print("Could not determine original b_iso...setting to 200", file = out)
    original_b_iso = 200.

  si.original_aniso_obj = map_coeffs_aa # set it so we can apply it later if desired

  if first_half_map_data:
    first_half_map_coeffs, dummy = get_f_phases_from_map(
          map_data = first_half_map_data,
       crystal_symmetry = si.crystal_symmetry,
       d_min = si.resolution,
       d_min_ratio = si.d_min_ratio,
       remove_aniso = si.remove_aniso,
       scale_max = si.scale_max,
       return_as_map_coeffs = True,
       out = local_out)
  else:
    first_half_map_coeffs = None

  if second_half_map_data:
    second_half_map_coeffs, dummy = get_f_phases_from_map(
      map_data = second_half_map_data,
       crystal_symmetry = si.crystal_symmetry,
       d_min = si.resolution,
       d_min_ratio = si.d_min_ratio,
       scale_max = si.scale_max,
       remove_aniso = si.remove_aniso,
       return_as_map_coeffs = True,
       out = local_out)
  else:
    second_half_map_coeffs = None

  if pdb_hierarchy:
    # Getting model information if pdb_hierarchy present -------------------
    from cctbx.maptbx.refine_sharpening import get_model_map_coeffs_normalized
    model_map_coeffs = get_model_map_coeffs_normalized(
     pdb_hierarchy = pdb_hierarchy,
       si = si,
       f_array = f_array,
       resolution = si.resolution,
       out = out)
    if not model_map_coeffs:  # give up
      pdb_hierarchy = None
      if si.is_model_sharpening():
        raise Sorry("Cannot carry out model sharpening without a model."+
            " It could be that the model was outside the map")

  else:
    model_map_coeffs = None

  # Try various methods for sharpening. # XXX fix this up

  local_si = deepcopy(si).update_with_box_sharpening_info(
      box_sharpening_info_obj = box_sharpening_info_obj)
  if si.adjust_region_weight and \
      (not si.sharpening_is_defined()) and (not si.is_model_sharpening()) \
     and (not si.is_half_map_sharpening()) and (
      not si.is_target_b_iso_to_d_cut()) and (
      si.sharpening_target == 'adjusted_sa'):
   for iii in range(1):  # just so we can break

    local_si = deepcopy(si).update_with_box_sharpening_info(
      box_sharpening_info_obj = box_sharpening_info_obj)
    local_si.sharpening_target = 'adjusted_sa'
    local_si.sharpening_method = 'b_iso_to_d_cut'


    sa_ratio_list = []
    normalized_regions_list = []

    if 0: #si.resolution:
      # 2017-07-26 reset b_low, b_mid, b_high, using 5.9*resolution**2 for b_mid
      delta_search = si.search_b_max-si.search_b_min
      b_mid = si.get_target_b_iso()
      b_low = b_mid-150*delta_search/400
      b_high = b_mid+250*delta_search/400
      print("Centering search on b_iso = %7.2f" %(b_mid), file = out)
    else:
      b_low = min(original_b_iso, si.search_b_min)
      b_high = max(original_b_iso, si.search_b_max)
      b_mid = b_low+0.375*(b_high-b_low)

    ok_region_weight = True
    results_list = []
    kw_list = []
    first = True

    id = 0
    for b_iso in [b_low, b_high, b_mid]:
      id+= 1

      if first and local_si.multiprocessing == 'multiprocessing' or \
          local_si.nproc == 1:  # can do anything
        local_log = out
      else:  # skip log entirely
        local_log = None  # will set this later and return as r.log_as_text
      first = False
      lsi = deepcopy(local_si)

      lsi.b_sharpen = original_b_iso-b_iso
      lsi.b_iso = b_iso

      # ------ SET UP RUN HERE ----------
      kw_list.append(
      {
      'f_array':f_array,
      'phases':phases,
        'crystal_symmetry':lsi.crystal_symmetry,
        'local_si':lsi,
        'id':id,
        'out':local_log,
       })
      # We are going to call autosharpening with this
      # ------ END OF SET UP FOR RUN ----------

      """

      local_map_and_b = apply_sharpening(
            f_array = f_array, phases = phases,
            sharpening_info_obj = lsi,
            crystal_symmetry = lsi.crystal_symmetry,
            out = null_out())
      local_si = score_map(map_data = local_map_and_b.map_data,
          sharpening_info_obj = local_si,
          out = null_out())
      """
      # This is the actual run here  =============

    from libtbx.easy_mp import run_parallel
    results_list = run_parallel(
       method = si.multiprocessing,
       qsub_command = si.queue_run_command,
       nproc = si.nproc,
       target_function = run_sharpen_and_score, kw_list = kw_list)
    # results looks like: [result, result2]

    sort_list = []
    for result in results_list:
      sort_list.append([result.id, result])
    sort_list.sort(key=itemgetter(0))
    for id, result in sort_list:
      local_si = result.local_si

      if local_si.sa_ratio is None or local_si.normalized_regions is None:
        ok_region_weight = False
      sa_ratio_list.append(local_si.sa_ratio)
      normalized_regions_list.append(local_si.normalized_regions)
    if not ok_region_weight:
      break # skip it

    # Set region weight so that either:
    #  (1) delta_sa_ratio == region_weight*delta_normalized_regions
    #  (2) sa_ratio = region_weight*normalized_regions (at low B)

    # region weight from change over entire region
    d_sa_ratio = sa_ratio_list[0]-sa_ratio_list[1]
    d_normalized_regions = normalized_regions_list[0]-normalized_regions_list[1]
    delta_region_weight = si.region_weight_factor*d_sa_ratio/max(
          1.e-10, d_normalized_regions)
    if d_sa_ratio < 0 or d_normalized_regions < 0:
      print("Not using delta_region_weight with unusable values", file = out)
      ok_region_weight = False


    # region weight from initial values
    init_region_weight = si.region_weight_factor* \
          sa_ratio_list[0]/max(1.e-10, normalized_regions_list[0])

    # Ensure that adjusted_sa at b_mid is > than either end
    #  adjusted_sa = sa_ratio - region_weight*normalized_regions

    # sa[2] = sa_ratio_list[2]-region_weight*normalized_regions[2]
    # sa[1] = sa_ratio_list[1]-region_weight*normalized_regions[1]
    # sa[0] = sa_ratio_list[0]-region_weight*normalized_regions[0]
    #  sa[2] >=  sa[1] and sa[2] >=  sa[0]
    # sa_ratio_list[2]-region_weight*normalized_regions[2] >=
    #    sa_ratio_list[1]-region_weight*normalized_regions[1]
    # NOTE: sa_ratio_list and normalized_regions both  decrease in order:
    #     low med high or [0] [2] [1]
    max_region_weight = (sa_ratio_list[2]- sa_ratio_list[1])/max(0.001,
         normalized_regions_list[2]-normalized_regions_list[1])
    min_region_weight = (sa_ratio_list[0]- sa_ratio_list[2])/max(0.001,
       normalized_regions_list[0]-normalized_regions_list[2])

    min_region_weight = max(1.e-10, min_region_weight) # positive
    max_region_weight = max(1.e-10, max_region_weight) # positive

    delta_weight = max(0., max_region_weight-min_region_weight)
    min_buffer = delta_weight*si.region_weight_buffer
    min_region_weight+= min_buffer
    max_region_weight-= min_buffer
    min_max_region_weight = True
    if min_region_weight >=  max_region_weight:
      print("Warning: min_region_weight >=  max_region_weight...", file = out)
      min_max_region_weight = False
      #ok_region_weight = False


    print("Region weight bounds: Min: %7.1f  Max: %7.1f " %(
      min_region_weight, max_region_weight), file = out)
    print("Region weight estimates:", file = out)
    print("From ratio of low-B surface area to regions: %7.1f" %(
     init_region_weight), file = out)
    print("Ratio of change in surface area to change in regions: %7.1f" %(
     delta_region_weight), file = out)

    # put them in bounds but note if we did it

    out_of_range = False
    if ok_region_weight and si.region_weight_method == 'initial_ratio':
      if min_max_region_weight and (
          init_region_weight > max_region_weight or \
          init_region_weight<min_region_weight):
        init_region_weight = max(
          min_region_weight, min(max_region_weight, init_region_weight))
        out_of_range = True
      print("\nRegion weight adjusted to %7.1f using initial ratio" %(
          init_region_weight), file = out)
      si.region_weight = init_region_weight

    elif ok_region_weight and si.region_weight_method == 'delta_ratio':
      if min_max_region_weight and (
          delta_region_weight > max_region_weight or \
          delta_region_weight<min_region_weight):
        delta_region_weight = max(
          min_region_weight, min(max_region_weight, delta_region_weight))
        out_of_range = True
      si.region_weight = delta_region_weight
      print("\nRegion weight set to %7.1f using overall ratio and " %(
          si.region_weight) +\
          "\nfactor of %5.1f" %(si.region_weight_factor), file = out)

    else: # just use default target for b_iso
      si.region_weight = si.region_weight_default

      print("Skipping region_weight analysis as signal-to-noise is zero ("+\
           "adjusted sa\nvs b_iso does not have low values at extremes and "+\
           "clear maximum in the middle.)", file = out)

      print("\nUnable to set region_weight ... using value of %7.2f" % (
          si.region_weight), file = out)
      if si.discard_if_worse:
        print("Setting discard_if_worse = False as region_weight failed ", file = out)
        si.discard_if_worse = False

    if out_of_range and auto_sharpen_methods and \
        'resolution_dependent' in auto_sharpen_methods:
      new_list = []
      have_something_left = False
      for x in auto_sharpen_methods:
        if x !=  'resolution_dependent':
          if str(x) !=  'None':
            have_something_left = True
          new_list.append(x)
      if have_something_left:
        auto_sharpen_methods = new_list
        print("Removed resolution_dependent sharpening ( "+\
          "weights were out of range)", file = out)

    if box_sharpening_info_obj:
      si.local_solvent_fraction = box_sharpening_info_obj.solvent_fraction
    else:
      si.local_solvent_fraction = si.solvent_fraction


  null_si = None
  best_si = deepcopy(si).update_with_box_sharpening_info(
      box_sharpening_info_obj = box_sharpening_info_obj)
  best_map_and_b = map_and_b_object()

  if si.sharpening_is_defined():  # Use this if come in with method
    print("\nUsing specified sharpening", file = out)
    best_si = set_up_sharpening(si = si, map_data = map_data, out = out)
    best_si.sharpen_and_score_map(map_data = map_data,
          out = out).show_score(out = out)
    best_si.show_summary(out = out)

  else:
    if best_si.is_model_sharpening():
      print("\nSetting up model sharpening", file = out)
    elif best_si.is_half_map_sharpening():
      print("\nSetting up half-map sharpening", file = out)
    elif best_si.is_external_map_sharpening():
      print("\nSetting up external map sharpening", file = out)
    else:
      print("\nTesting sharpening methods with target of %s" %(
        best_si.sharpening_target), file = out)
    if not auto_sharpen_methods or auto_sharpen_methods == ['None']:
      auto_sharpen_methods = ['no_sharpening']
    for m in auto_sharpen_methods:
      # ------------------------
      if m in ['no_sharpening', 'resolution_dependent', 'model_sharpening',
          'half_map_sharpening', 'target_b_iso_to_d_cut',
           'external_map_sharpening']:
        if m == 'target_b_iso_to_d_cut':
          b_min = si.get_target_b_iso()
          b_max = si.get_target_b_iso()
        else:
          b_min = original_b_iso
          b_max = original_b_iso
        b_n = 1
        k_sharpen = 0.
        delta_b = 0
        if m in ['resolution_dependent', 'model_sharpening',
           'half_map_sharpening', 'external_map_sharpening']:
          pass # print out later
        else:
          print("\nB-sharpen   B-iso   k_sharpen   SA   "+\
             "Kurtosis  Path len  Normalized regions", file = out)
      # ------------------------
      # ------------------------
      else:  #  ['b_iso', 'b_iso_to_d_cut']:
        if si.search_b_n>1:
          b_min = min(original_b_iso, si.search_b_min)
          b_max = max(original_b_iso, si.search_b_max)
        else: # for just one, take it
          b_min = si.search_b_min
          b_max = si.search_b_max
        b_n = si.search_b_n
        delta_b = (b_max-b_min)/max(1, b_n-1)
        print("\nTesting %s with b_iso from %7.1f to %7.1f in %d steps of %7.1f" %(
          m, b_min, b_max, b_n, delta_b), file = out)
        print("(b_sharpen from %7.1f to %7.1f ) " %(
           original_b_iso-b_min, original_b_iso-b_max), file = out)
        if m == 'b_iso':
          k_sharpen = 0.
        else:
          k_sharpen = si.k_sharpen

        print("\nB-sharpen   B-iso   k_sharpen   SA   "+\
             "Kurtosis  Path len  Normalized regions", file = out)
      # ------------------------
      local_best_map_and_b = map_and_b_object()
      local_best_si = deepcopy(si).update_with_box_sharpening_info(
        box_sharpening_info_obj = box_sharpening_info_obj)

      si_b_iso_list = flex.double()
      si_score_list = flex.double()
      si_id_list = []

      kw_list = []
      first = True
      if return_bsi: assert local_si.nproc == 1

      results_list = []
      for i in range(b_n):
      #============================================
        local_si = deepcopy(si).update_with_box_sharpening_info(
          box_sharpening_info_obj = box_sharpening_info_obj)
        local_si.sharpening_method = m
        local_si.n_real = map_data.all()
        local_si.k_sharpen = k_sharpen

        if first and local_si.multiprocessing == 'multiprocessing' or \
            local_si.nproc == 1:  # can do anything
          local_log = out
        else:  # skip log entirely
          local_log = None  # will set this later and return as r.log_as_text
        first = False


        if m == 'resolution_dependent':
          print("\nRefining resolution-dependent sharpening based on %s" %(
            local_si.residual_target), file = out)
          local_si.b_sharpen = 0
          local_si.b_iso = original_b_iso
          from cctbx.maptbx.refine_sharpening import run as refine_sharpening
          local_f_array, local_phases = refine_sharpening(
             map_coeffs = map_coeffs,
             sharpening_info_obj = local_si,
             out = out)
        elif m == 'model_sharpening':
          print("\nUsing model-based sharpening", file = out)
          local_si.b_sharpen = 0
          local_si.b_iso = original_b_iso
          from cctbx.maptbx.refine_sharpening import scale_amplitudes
          scale_amplitudes(
            model_map_coeffs = model_map_coeffs, map_coeffs = map_coeffs,
            si = local_si, out = out)

          # local_si contains target_scale_factors now
          local_f_array = f_array
          local_phases = phases
        elif m == 'half_map_sharpening':
          print("\nUsing half-map-based sharpening", file = out)
          local_si.b_sharpen = 0
          local_si.b_iso = original_b_iso
          from cctbx.maptbx.refine_sharpening import scale_amplitudes
          scale_amplitudes(
            model_map_coeffs = model_map_coeffs,
            map_coeffs = map_coeffs,
            first_half_map_coeffs = first_half_map_coeffs,
            second_half_map_coeffs = second_half_map_coeffs,
            si = local_si, out = out)
          # local_si contains target_scale_factors now
          local_f_array = f_array
          local_phases = phases
        elif m == 'external_map_sharpening':
          print("\nUsing external-map-based sharpening", file = out)
          local_si.b_sharpen = 0
          local_si.b_iso = original_b_iso
          from cctbx.maptbx.refine_sharpening import scale_amplitudes
          scale_amplitudes(
            model_map_coeffs = model_map_coeffs,
            map_coeffs = map_coeffs,
            external_map_coeffs = first_half_map_coeffs,
            si = local_si, out = out)
          # local_si contains target_scale_factors now
          local_f_array = f_array
          local_phases = phases

        else:
          local_f_array = f_array
          local_phases = phases
          b_iso = b_min+i*delta_b
          local_si.b_sharpen = original_b_iso-b_iso
          local_si.b_iso = b_iso

        # ------ SET UP RUN HERE ----------
        kw_list.append(
        {
        'f_array':local_f_array,
        'phases':local_phases,
          'crystal_symmetry':local_si.crystal_symmetry,
          'original_b_iso':original_b_iso,
          'local_si':local_si,
          'm':m,
          'return_bsi':return_bsi,
          'out':local_log,
          'id':i+1,
         })
        # We are going to call autosharpening with this
        # ------ END OF SET UP FOR RUN ----------


        # This is the actual run here  =============

      from libtbx.easy_mp import run_parallel
      results_list = run_parallel(
         method = si.multiprocessing,
         qsub_command = si.queue_run_command,
         nproc = si.nproc,
         target_function = run_sharpen_and_score, kw_list = kw_list)
      # results looks like: [result, result2]

      sort_list = []
      for result in results_list:
        sort_list.append([result.id, result])
      sort_list.sort(key=itemgetter(0))
      for id, result in sort_list:
        local_si = result.local_si
        local_map_and_b = result.local_map_and_b
        if result.text:
          print(result.text, file = out)
        # Run through all result to get these
        if local_si.b_sharpen is not None and local_si.b_iso is not None and\
           local_si.k_sharpen is not None and local_si.kurtosis is not None \
           and local_si.adjusted_sa is not None and local_si.score is not None:
            si_b_iso_list.append(local_si.b_iso)
            si_score_list.append(local_si.score)
            if local_si.k_sharpen is not None:
             si_id_list.append("%.3f_%.3f_%.3f" %(
               local_si.b_iso, local_si.k_sharpen,
               local_si.get_d_cut()))

        if m == 'no_sharpening':
          null_si = local_si
        if local_best_si.score is None or local_si.score>local_best_si.score:
          local_best_si = local_si
          local_best_map_and_b = local_map_and_b

        #  ============================================
        # DONE WITH ALL RUNS

      if not local_best_si.is_model_sharpening() and \
          not local_best_si.is_half_map_sharpening():
        if local_best_si.sharpening_method == 'resolution_dependent':
          print("\nBest scores for sharpening with "+\
            "b[0] = %6.2f b[1] = %6.2f b[2] = %6.2f: " %(
            local_best_si.resolution_dependent_b[0],
            local_best_si.resolution_dependent_b[1],
            local_best_si.resolution_dependent_b[2]), file = out)
        else:
          print("\nBest scores for sharpening with "+\
            "b_iso = %6.1f b_sharpen = %6.1f k_sharpen = %s: " %(
            local_best_si.b_iso, local_best_si.b_sharpen,
             local_best_si.k_sharpen), file = out)
        if local_best_si.score is not None:
          local_best_si.show_summary(out = out)
          print("Adjusted surface area: %7.3f  Kurtosis: %7.3f  Score: %7.3f\n" %(
           local_best_si.adjusted_sa, local_best_si.kurtosis, local_best_si.score), file = out)

      if  si_score_list.size()>1: # test for signal
        signal_to_noise = estimate_signal_to_noise(value_list = si_score_list)
        print("Estimated signal-to-noise in ID of optimal sharpening: %5.1f" %(
           signal_to_noise), file = out)
        if signal_to_noise<local_best_si.signal_min and \
            'target_b_iso_to_d_cut' in auto_sharpen_methods:
          print("Skipping this analysis as signal-to-noise is less than %5.1f " %(
           local_best_si.signal_min), file = out)
          local_best_si.score = None

      optimize_b_blur_hires = False
      optimize_d_cut = False
      n_cycles = 0
      if local_best_si.score is not None and local_best_si.optimize_d_cut and \
        local_best_si.sharpening_method in ['b_iso_to_d_cut', 'b_iso']:
        optimize_d_cut = True
        n_cycles+= 1
      if local_best_si.score is not None and \
        local_best_si.optimize_b_blur_hires and \
        local_best_si.k_sharpen is not None and \
        local_best_si.sharpening_method in ['b_iso_to_d_cut', 'b_iso']:
        optimize_b_blur_hires = True
        n_cycles+= 1

      ##########################################
      optimize_b_iso = True
      for cycle in range(n_cycles):
        if optimize_b_blur_hires:
          local_best_si, local_best_map_and_b = optimize_b_blur_or_d_cut_or_b_iso(
           #optimization_target = 'k_sharpen',
           optimization_target = 'b_blur_hires',
           local_best_si = local_best_si,
           local_best_map_and_b = local_best_map_and_b,
           si_id_list = si_id_list,
           si_score_list = si_score_list,
           delta_b = delta_b,
           original_b_iso = original_b_iso,
           f_array = f_array,
           phases = phases,
           out = out)

        if optimize_d_cut:
          local_best_si, local_best_map_and_b = optimize_b_blur_or_d_cut_or_b_iso(
           optimization_target = 'd_cut',
           local_best_si = local_best_si,
           local_best_map_and_b = local_best_map_and_b,
           si_id_list = si_id_list,
           si_score_list = si_score_list,
           delta_b = delta_b,
           original_b_iso = original_b_iso,
           f_array = f_array,
           phases = phases,
           out = out)

        if optimize_b_iso:
          local_best_si, local_best_map_and_b = optimize_b_blur_or_d_cut_or_b_iso(
           optimization_target = 'b_iso',
           local_best_si = local_best_si,
           local_best_map_and_b = local_best_map_and_b,
           si_id_list = si_id_list,
           si_score_list = si_score_list,
           delta_b = delta_b,
           original_b_iso = original_b_iso,
           f_array = f_array,
           phases = phases,
           out = out)

      ##########################################
      if (local_best_si.score is not None or
         local_best_si.is_model_sharpening()) and (
          best_si.score is None or local_best_si.score > best_si.score):
        best_si = local_best_si
        best_map_and_b = local_best_map_and_b
        if not best_si.is_model_sharpening() and \
            not best_si.is_half_map_sharpening():
          print("This is the current best score\n", file = out)

  if (best_si.score is not None )  and (
     not best_si.is_model_sharpening() )  and (
     not best_si.is_half_map_sharpening()):
    print("\nOverall best sharpening method: %s Score: %7.3f\n" %(
       best_si.sharpening_method, best_si.score), file = out)
    best_si.show_summary(out = out)
  if (not best_si.is_model_sharpening()) and \
       (not best_si.is_half_map_sharpening()) and null_si:
    if best_si.score>null_si.score:  # we improved them..
      print("Improved score with sharpening...", file = out)
    else:
      print("Did not improve score with sharpening...", file = out)
  if return_bsi:
    map_data = best_map_and_b.map_data
    if not map_data: # no result
       return None
    map_data = set_mean_sd_of_map(map_data = map_data,
      target_mean = starting_mean, target_sd = starting_sd)
    box_sharpening_info_obj.map_data = map_data
    box_size= map_data.all()
    calculated_box_size=tuple([i-j+1 for i,j in zip(
       box_sharpening_info_obj.upper_bounds,
       box_sharpening_info_obj.lower_bounds)])
    calculated_box_size_minus_one=tuple([i-j for i,j in zip(
       box_sharpening_info_obj.upper_bounds,
       box_sharpening_info_obj.lower_bounds)])
    if calculated_box_size_minus_one == box_size: # one too big
      box_sharpening_info_obj.upper_bounds =tuple(
       [i-1 for i in box_sharpening_info_obj.upper_bounds]
      ) #  work-around for upper bounds off by one in model sharpening
    box_sharpening_info_obj.smoothed_box_mask_data = smoothed_box_mask_data
    box_sharpening_info_obj.original_box_map_data = original_box_map_data
    box_sharpening_info_obj.n_buffer = n_buffer
    box_sharpening_info_obj.crystal_symmetry = best_si.crystal_symmetry
    box_sharpening_info_obj.resolution = best_si.resolution
    box_sharpening_info_obj.d_min_ratio = best_si.d_min_ratio
    box_sharpening_info_obj.scale_max = best_si.scale_max
    box_sharpening_info_obj.smoothing_radius = best_si.smoothing_radius
    box_sharpening_info_obj.b_iso = best_map_and_b.final_b_iso
    box_sharpening_info_obj.starting_b_iso = best_map_and_b.starting_b_iso
    return box_sharpening_info_obj

  if original_box_sharpening_info_obj:
      # Put back original crystal_symmetry with original_box_sharpening_info_obj
      print("\nRestoring original symmetry to best sharpening info", file = out)
      best_si.update_with_box_sharpening_info(
        box_sharpening_info_obj = original_box_sharpening_info_obj)
      print("(%7.3f, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f) "%(tuple(
        best_si.crystal_symmetry.unit_cell().parameters())), file = out)
      # and set tracking data with result
  return best_si

def run_sharpen_and_score(f_array = None,
  phases = None,
  local_si = None,
  crystal_symmetry = None,
  original_b_iso = None,
  m = None,
  return_bsi = None,
  id = None,
  out = sys.stdout):
        """Apply sharpening to a map and score it"""

        local_map_and_b = apply_sharpening(
            f_array = f_array, phases = phases,
            sharpening_info_obj = local_si,
            crystal_symmetry = crystal_symmetry,
            out = null_out())
        local_si = score_map(map_data = local_map_and_b.map_data,
          sharpening_info_obj = local_si,
          out = null_out())
        # Record b_iso values
        if not local_map_and_b.starting_b_iso:
          local_map_and_b.starting_b_iso = original_b_iso
        if not local_map_and_b.final_b_iso:
          local_map_and_b.final_b_iso = local_si.b_iso

        # This is printout below here  ===============

        if m == 'resolution_dependent':
          text = \
           "\nb[0]   b[1]   b[2]   SA   Kurtosis   sa_ratio  Normalized regions"
          text+= "\n"+\
            "\nB-sharpen   B-iso   k_sharpen   SA   "+\
             "Kurtosis  Path len  Normalized regions"
          text+= "\n"+" %6.2f  %6.2f  %6.2f  " %(
              local_si.resolution_dependent_b[0],
              local_si.resolution_dependent_b[1],
              local_si.resolution_dependent_b[2]) +\
            "  %7.3f  %7.3f  " %(
                local_si.adjusted_sa, local_si.kurtosis)+\
            " %7.1f  %7.3f" %(
             zero_if_none(local_si.adjusted_path_length), #local_si.sa_ratio,
             local_si.normalized_regions)
        elif local_si.b_sharpen is not None and local_si.b_iso is not None and\
           local_si.k_sharpen is not None and local_si.kurtosis is not None \
           and local_si.adjusted_sa is not None:
          text = \
           " %6.1f     %6.1f  %5s   %7.3f  %7.3f" %(
            local_si.b_sharpen, local_si.b_iso,
             local_si.k_sharpen, local_si.adjusted_sa, local_si.kurtosis) + \
            "  %7.1f         %7.3f" %(
             zero_if_none(local_si.adjusted_path_length), #local_si.sa_ratio,
             local_si.normalized_regions)
        else:
          text = ""

        if return_bsi:
          r = group_args(
            local_si = local_si,
            local_map_and_b = local_map_and_b,
            text = text,
            id = id)
        else:
          r = group_args(
            local_si = local_si,
            local_map_and_b = None,
            text = text,
            id = id)
        return r

def effective_b_iso(map_data = None, tracking_data = None,
      box_sharpening_info_obj = None,
      crystal_symmetry = None,
      resolution = None,
      remove_aniso = None,
      d_min_ratio = None,
      scale_max = None,
      out = sys.stdout):
    """Estimate effective b_iso from a map"""
    if not crystal_symmetry:
      if box_sharpening_info_obj:
        crystal_symmetry = box_sharpening_info_obj.crystal_symmetry
      else:
        crystal_symmetry = tracking_data.crystal_symmetry
    if resolution:
       d_min = resolution
    else:
       d_min = tracking_data.params.crystal_info.resolution

    if not d_min_ratio:
       d_min_ratio = tracking_data.params.map_modification.d_min_ratio

    map_coeffs, map_coeffs_ra = get_f_phases_from_map(map_data = map_data,
       crystal_symmetry = crystal_symmetry,
       d_min = d_min,
       d_min_ratio = d_min_ratio,
       scale_max = scale_max,
       remove_aniso = remove_aniso,
       return_as_map_coeffs = True,
       out = out)

    f_array, phases = map_coeffs_as_fp_phi(map_coeffs)
    if map_coeffs_ra:
      b_iso = map_coeffs_ra.b_iso
    else:
      b_iso = None
    if b_iso is not None:
      print("Effective B-iso = %7.2f\n" %(b_iso), file = out)
    else:
      print("Effective B-iso not determined\n", file = out)
    return map_coeffs_ra, map_coeffs, f_array, phases

def update_tracking_data_with_sharpening(map_data = None, tracking_data = None,
       si = None, out = sys.stdout):
    """Update tracking data with sharpening data from si"""

    # Set shifted_map_info if map_data is new
    if tracking_data.params.output_files.shifted_sharpened_map_file:
      shifted_sharpened_map_file = os.path.join(
          tracking_data.params.output_files.output_directory,
          tracking_data.params.output_files.shifted_sharpened_map_file)
    else:
      shifted_sharpened_map_file = None
    from cctbx.maptbx.segment_and_split_map import write_ccp4_map
    if shifted_sharpened_map_file:
      write_ccp4_map(tracking_data.crystal_symmetry,
          shifted_sharpened_map_file, map_data)
      print("Wrote shifted, sharpened map to %s" %(
          shifted_sharpened_map_file), file = out)
    tracking_data.set_shifted_map_info(file_name =
          shifted_sharpened_map_file,
          crystal_symmetry = tracking_data.crystal_symmetry,
          origin = map_data.origin(),
          all = map_data.all(),
          b_sharpen = None)


def get_high_points_from_map(
     map_data = None,
     boundary_radius = 5.,
     unit_cell = None,
     out = sys.stdout):
    """Find the highest point in a map not within boundary_radius of edge
     of map"""

    max_in_map_data = map_data.as_1d().min_max_mean().max
    for cutoff in [0.99, 0.98, 0.95, 0.90, 0.50]:
      high_points_mask = (map_data>=  cutoff*max_in_map_data)
      sda = map_data.as_1d().min_max_mean().max
      for nth_point in [4, 2, 1]:
        sites_cart = get_marked_points_cart(mask_data = high_points_mask,
          unit_cell = unit_cell, every_nth_point = nth_point,
          boundary_radius = boundary_radius)
        if sites_cart.size()>0: break
      if sites_cart.size()>0: break
    assert sites_cart.size()>0
    del high_points_mask
    sites_cart = sites_cart[:1]
    xyz_frac = unit_cell.fractionalize(sites_cart[0])
    value = map_data.value_at_closest_grid_point(xyz_frac)
    print("High point in map at (%7.2f %7.2f %7.2f) with value of %7.2f " %(
        sites_cart[0][0], sites_cart[0][1], sites_cart[0][2], value), file = out)
    return sites_cart

def get_one_au(tracking_data = None,
    sites_cart = None,
    ncs_obj = None,
    map_data = None,
    starting_mask = None,
    radius = None,
    every_nth_point = None,
    removed_ncs = None,
    out = sys.stdout):
  """Return mask marking one asymmetric unit (of NCS)"""
  unit_cell = tracking_data.crystal_symmetry.unit_cell()

  if removed_ncs: # take everything left
    mm = map_data.as_1d().min_max_mean()
    mask_threshold = mm.min+max(0.00001, 0.0001*(mm.mean-mm.min)) # just above min
  else:
    mask_threshold = tracking_data.params.segmentation.mask_threshold

  every_nth_point = tracking_data.params.segmentation.grid_spacing_for_au
  radius = tracking_data.params.segmentation.radius
  if not radius:
    radius = set_radius(unit_cell = unit_cell, map_data = map_data,
     every_nth_point = every_nth_point)
    tracking_data.params.segmentation.radius = radius
  print("\nRadius for AU identification: %7.2f A" %(radius), file = out)

  overall_mask, max_in_sd_map, sd_map = get_overall_mask(map_data = map_data,
    mask_threshold = mask_threshold,
    crystal_symmetry = tracking_data.crystal_symmetry,
    resolution = tracking_data.params.crystal_info.resolution,
    solvent_fraction = tracking_data.solvent_fraction,
    radius = radius,
    out = out)

  if starting_mask:
    print("Points in starting mask:", starting_mask.count(True), file = out)
    print("Points in overall mask:", overall_mask.count(True), file = out)
    print("Points in both:", (starting_mask & overall_mask).count(True), file = out)
    if tracking_data.params.crystal_info.is_crystal:
      # take starting mask as overall...
      overall_mask =  starting_mask
    else: # usual
      # make sure overall mask is at least as big..
      overall_mask = (overall_mask | starting_mask)
    print("New size of overall mask: ", overall_mask.count(True), file = out)
  else:
    if not sites_cart: # pick top of map
      sites_cart = get_high_points_from_map(
        boundary_radius = radius,
        map_data = sd_map,
        unit_cell = unit_cell, out = out)

    starting_mask = mask_from_sites_and_map( # starting au mask
      map_data = sd_map, unit_cell = unit_cell,
      sites_cart = sites_cart, radius = radius,
      overall_mask = overall_mask)

  del sd_map

  au_mask, ncs_mask = get_ncs_mask(
    map_data = map_data, unit_cell = unit_cell, ncs_object = ncs_obj,
    starting_mask = starting_mask,
    radius = radius,
    overall_mask = overall_mask,
    every_nth_point = every_nth_point)

  print("Points in au: %d  in ncs: %d  (total %7.1f%%)   both: %d Not marked: %d" %(
     au_mask.count(True), ncs_mask.count(True),
     100.*float(au_mask.count(True)+ncs_mask.count(True))/au_mask.size(),
     (au_mask & ncs_mask).count(True),
     au_mask.size()-au_mask.count(True)-ncs_mask.count(True), ), file = out)

  return au_mask

def set_up_sharpening(si = None, map_data = None, out = sys.stdout):
         """Set up sharpening specified by si object"""
         print("\nCarrying out specified sharpening/blurring of map", file = out)
         check_si = si  # just use input information
         check_si.show_summary(out = out)
         if check_si.is_target_b_iso_to_d_cut():
           check_si.b_iso = check_si.get_target_b_iso()
           check_si.b_sharpen = None
           print("Setting target b_iso of %7.1f " %(check_si.b_iso), file = out)
         if check_si.b_sharpen is None and check_si.b_iso is not None:
           # need to figure out b_sharpen
           print("\nGetting b_iso of map", file = out)
           b_iso = check_si.get_effective_b_iso(map_data = map_data, out = out)
           check_si.b_sharpen = b_iso-check_si.b_iso # sharpen is what to
           print("Value of b_sharpen to obtain b_iso of %s is %5.2f" %(
             check_si.b_iso, check_si.b_sharpen), file = out)
         elif check_si.b_sharpen is not None:
           print("Sharpening b_sharpen will be %s" %(check_si.b_sharpen), file = out)
         elif check_si.resolution_dependent_b:
           print("Resolution-dependent b_sharpening values:" +\
              "b0: %7.2f  b1: %7.2f  b2: %7.2f " %(
             tuple(check_si.resolution_dependent_b)), file = out)
         elif check_si.target_scale_factors:
           print("Model sharpening scale values:", file = out)
           for x in check_si.target_scale_factors: print(x, end = ' ', file = out)
           print(file = out)
         return check_si

def run(args,
     params = None,
     map_data = None,
     crystal_symmetry = None,
     write_files = None,
     auto_sharpen = None,
     density_select = None,
     add_neighbors = None,
     save_box_map_ncs_au = None,
     resolution = None,
     sequence = None,
     half_map_data_list = None,
     ncs_obj = None,
     tracking_data = None,
     target_scattered_points = None,
     is_iteration = False,
     pdb_hierarchy = None,
     target_xyz = None,
     target_hierarchy = None,
     target_model = None,
     sharpening_target_pdb_hierarchy = None,
     wrapping = None,
     target_ncs_au_file = None,
     regions_to_keep = None,
     solvent_content = None,
     molecular_mass = None,
     symmetry = None,
     chain_type = None,
     keep_low_density = None,
     box_buffer = None,
     soft_mask_extract_unique = None,
     mask_expand_ratio = None,
     out = sys.stdout):

  """Segment and split a map into regions"""

  if is_iteration:
    print("\nIteration tracking data:", file = out)
    tracking_data.show_summary(out = out)
  else:
    # get the parameters and map_data (sharpened, magnified, shifted...)
    params, map_data, half_map_data_list, pdb_hierarchy, tracking_data, \
        shifted_ncs_object = get_params( #
       args, map_data = map_data, crystal_symmetry = crystal_symmetry,
       half_map_data_list = half_map_data_list,
       ncs_object = ncs_obj,
       write_files = write_files,
       auto_sharpen = auto_sharpen,
       density_select = density_select,
       add_neighbors = add_neighbors,
       save_box_map_ncs_au = save_box_map_ncs_au,
       sequence = sequence,
       wrapping = wrapping,
       target_ncs_au_file = target_ncs_au_file,
       regions_to_keep = regions_to_keep,
       solvent_content = solvent_content,
       resolution = resolution,
       molecular_mass = molecular_mass,
       symmetry = symmetry,
       chain_type = chain_type,
       keep_low_density = keep_low_density,
       box_buffer = box_buffer,
       soft_mask_extract_unique = soft_mask_extract_unique,
       mask_expand_ratio = mask_expand_ratio,
       sharpening_target_pdb_hierarchy = sharpening_target_pdb_hierarchy,
        out = out)
    if params.control.shift_only:
      return map_data, ncs_obj, tracking_data
    elif params.control.check_ncs or \
        params.control.sharpen_only:
      return None, None, tracking_data

    if params.input_files.pdb_to_restore:
      restore_pdb(params, tracking_data = tracking_data, out = out)
      return None, None, tracking_data
    # read and write the ncs (Normally point-group NCS)
    ncs_obj, tracking_data = get_ncs(params = params, tracking_data = tracking_data,
       ncs_object = shifted_ncs_object,
       out = out)

    if target_model:
      target_hierarchy = target_model.get_hierarchy()
    elif params.input_files.target_ncs_au_file: # read in target
      from iotbx.pdb.utils import get_pdb_hierarchy
      pdb_hierarchy = get_pdb_hierarchy(
         file_name = params.input_files.target_ncs_au_file)

    print("\nShifting model based on origin shift (if any)", file = out)
    print("Coordinate shift is (%7.2f, %7.2f, %7.2f)" %(
        tuple(tracking_data.origin_shift)), file = out)
    if not map_data:
       raise Sorry("Need map data for segment_and_split_map")

    if params.output_files.shifted_map_file:
        shifted_map_file = os.path.join(
          tracking_data.params.output_files.output_directory,
          params.output_files.shifted_map_file)
    else:
      shifted_map_file = None
    if params.output_files.shifted_ncs_file:
        shifted_ncs_file = os.path.join(
          tracking_data.params.output_files.output_directory,
          params.output_files.shifted_ncs_file)
    else:
      shifted_ncs_file = None
    if params.output_files.shifted_ncs_file:
        shifted_pdb_file = os.path.join(
          tracking_data.params.output_files.output_directory,
          params.output_files.shifted_pdb_file)
    else:
      shifted_pdb_file = None
    shifted_pdb_file, ncs_obj, pdb_hierarchy, target_hierarchy, \
      tracking_data= apply_origin_shift(
        shifted_map_file = shifted_map_file,
        shifted_pdb_file = shifted_pdb_file,
        shifted_ncs_file = shifted_ncs_file,
        origin_shift = tracking_data.origin_shift,
        shifted_ncs_object = shifted_ncs_object,
        pdb_hierarchy = pdb_hierarchy,
        target_hierarchy = target_hierarchy,
        map_data = map_data,
        tracking_data = tracking_data,
        out = out)
    if shifted_pdb_file:
      params.output_files.shifted_pdb_file = os.path.split(shifted_pdb_file)[-1]

    if target_hierarchy:
      target_xyz = target_hierarchy.atoms().extract_xyz()
      del target_hierarchy

    # We can use params.input_files.target_ncs_au_file here to define ncs au
    if target_xyz and not target_scattered_points:
       target_scattered_points = flex.vec3_double()
       target_scattered_points.append(target_xyz.mean())

    # get the chain types and therefore (using ncs_copies) volume fraction
    tracking_data = get_solvent_fraction(params,
      ncs_object = ncs_obj, tracking_data = tracking_data, out = out)

    # Done with getting params and maps
    # Summarize after any sharpening
    tracking_data.show_summary(out = out)

  original_ncs_obj = ncs_obj # in case we need it later...
  original_input_ncs_info = tracking_data.input_ncs_info
  removed_ncs = False

  n_residues = tracking_data.n_residues
  ncs_copies = tracking_data.input_ncs_info.number_of_operators
  if (not tracking_data.solvent_fraction) and \
      params.crystal_info.molecular_mass:
    tracking_data.solvent_fraction = get_solvent_fraction_from_molecular_mass(
        crystal_symmetry = tracking_data.crystal_symmetry,
        molecular_mass = params.crystal_info.molecular_mass,
        out = out)
  if tracking_data.solvent_fraction:
    solvent_fraction = tracking_data.solvent_fraction
  else:
    raise Sorry("Need solvent fraction or molecular mass or sequence file")

  # Now usual method, using our new map...should duplicate best result above
  for itry in range(2):
    # get connectivity  (conn = connectivity_object.result)
    b_vs_region = b_vs_region_info()
    si = sharpening_info(tracking_data = tracking_data)
    co, sorted_by_volume, min_b, max_b, unique_expected_regions, best_score, \
       new_threshold, starting_density_threshold = \
         get_connectivity(
           b_vs_region = b_vs_region,
           map_data = map_data,
           iterate_with_remainder = params.segmentation.iterate_with_remainder,
           n_residues = n_residues,
           ncs_copies = ncs_copies,
           solvent_fraction = solvent_fraction,
           fraction_occupied = si.fraction_occupied,
           min_volume = si.min_volume,
           min_ratio = si.min_ratio,
           wrapping = si.wrapping,
           residues_per_region = si.residues_per_region,
           max_ratio_to_target = si.max_ratio_to_target,
           min_ratio_to_target = si.min_ratio_to_target,
           min_ratio_of_ncs_copy_to_first = si.min_ratio_of_ncs_copy_to_first,
           starting_density_threshold = si.starting_density_threshold,
           density_threshold = si.density_threshold,
           crystal_symmetry = si.crystal_symmetry,
           chain_type = si.chain_type,
           verbose = si.verbose,
           out = out)
    params.segmentation.starting_density_threshold = starting_density_threshold # have to set tracking data as we are passing that above
    tracking_data.params.segmentation.starting_density_threshold = starting_density_threshold # have to set tracking data as we are passing that above
    if new_threshold:
      print("\nNew threshold is %7.2f" %(new_threshold), file = out)
    if co is None: # no luck
      return None, None, tracking_data

    # Check to see which regions are in more than one au of the NCS
    #   and set them aside.  Group ncs-related regions together

    ncs_group_obj, tracking_data, equiv_dict_ncs_copy = identify_ncs_regions(
       params, sorted_by_volume = sorted_by_volume,
       co = co,
       min_b = min_b,
       max_b = max_b,
       ncs_obj = ncs_obj,
       tracking_data = tracking_data,
       out = out)
    if ncs_group_obj and ncs_group_obj.ncs_group_list: # ok
      break
    elif ncs_obj and itry == 0 and not is_iteration:# try again
      print("No NCS groups identified on first try...taking entire NCS AU.", file = out)
      # Identify ncs au
      au_mask = get_one_au(tracking_data = tracking_data,
        ncs_obj = ncs_obj,
        map_data = map_data, out = out)
      s = (au_mask == False)
      min_in_map = map_data.as_1d().min_max_mean().min
      map_data.set_selected(s, min_in_map)  # mask out all but au
      from mmtbx.ncs.ncs import ncs
      ncs_obj = ncs()
      ncs_obj.set_unit_ncs()
      tracking_data.set_ncs_obj(ncs_obj = None)
      tracking_data.update_ncs_info(number_of_operators = 1)
      if n_residues:
        n_residues = n_residues/ncs_copies
      solvent_fraction = max(0.001, min(0.999,
       1-((1-solvent_fraction)/ncs_copies)))
      ncs_copies = 1
      params.segmentation.require_complete = False
      params.segmentation.iterate_with_remainder = False # so we do not iterate
      removed_ncs = True
      # Run again
    else: # tried twice, give up
      return None, None, tracking_data

  # Choose one region or group of regions from each ncs_group in the list
  #  Optimize the closeness of centers

  # Select group of regions that are close together and represent one au

  ncs_group_obj, scattered_points = \
     select_regions_in_au(
     params,
     ncs_group_obj = ncs_group_obj,
     equiv_dict_ncs_copy = equiv_dict_ncs_copy,
     tracking_data = tracking_data,
     target_scattered_points = target_scattered_points,
     unique_expected_regions = unique_expected_regions,
     out = out)

  # write out mask and map for all the selected regions...

  # Iterate if desired
  if params.segmentation.iterate_with_remainder and \
      ncs_group_obj.selected_regions:

    print("\nCreating remaining mask and map", file = out)
    map_data_remaining = create_remaining_mask_and_map(params,
      ncs_group_obj = ncs_group_obj,
      map_data = map_data,
      crystal_symmetry = tracking_data.crystal_symmetry,
      out = out)

    remainder_ncs_group_obj = iterate_search(params,
      map_data = map_data,
      map_data_remaining = map_data_remaining,
      ncs_obj = ncs_obj,
      ncs_group_obj = ncs_group_obj,
      scattered_points = scattered_points,
      tracking_data = tracking_data,
      out = out)
  else:
    remainder_ncs_group_obj = None

  # collect all NCS ops that are needed to relate all the regions
  #  that are used
  ncs_ops_used = ncs_group_obj.ncs_ops_used
  if remainder_ncs_group_obj and remainder_ncs_group_obj.ncs_ops_used:
    for x in remainder_ncs_group_obj.ncs_ops_used:
      if not x in ncs_ops_used: ncs_ops_used.append(x)
  if ncs_ops_used:
    ncs_ops_used.sort()
    print("Final NCS ops used: ", ncs_ops_used, file = out)

  # Save the used NCS ops
  ncs_used_obj = ncs_group_obj.ncs_obj.deep_copy(ops_to_keep = ncs_ops_used)
  if params.output_files.shifted_used_ncs_file:
    shifted_used_ncs_file = os.path.join(
      tracking_data.params.output_files.output_directory,
      params.output_files.shifted_used_ncs_file)
    ncs_used_obj.format_all_for_group_specification(
         file_name = shifted_used_ncs_file)
    tracking_data.set_shifted_used_ncs_info(file_name = shifted_used_ncs_file,
      number_of_operators = ncs_used_obj.max_operators(),
      is_helical_symmetry = tracking_data.input_ncs_info.is_helical_symmetry)
    tracking_data.shifted_used_ncs_info.show_summary(out = out)

  # Write out final maps and dummy atom files
  if params.output_files.write_output_maps:
    print("\nWriting output maps", file = out)
  else:
    print("\nSetting up but not writing output maps", file = out)
  map_files_written = write_output_files(params,
      tracking_data = tracking_data,
      map_data = map_data,
      half_map_data_list = half_map_data_list,
      ncs_group_obj = ncs_group_obj,
      remainder_ncs_group_obj = remainder_ncs_group_obj,
      pdb_hierarchy = pdb_hierarchy,
      removed_ncs = removed_ncs,
      out = out)
  ncs_group_obj.set_map_files_written(map_files_written)

  # Restore ncs info if we removed it
  if removed_ncs:
    print("\nRestoring original NCS info to tracking_data", file = out)
    tracking_data.input_ncs_info = original_input_ncs_info


  if params.output_files.output_info_file and ncs_group_obj:
    write_info_file(params = params, tracking_data = tracking_data, out = out)

  return ncs_group_obj, remainder_ncs_group_obj, tracking_data


if __name__ == "__main__":
  run(args = sys.argv[1:])
