# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.merge
# LIBTBX_SET_DISPATCHER_NAME xfel.merge
#
# $Id$

from __future__ import division
from __future__ import print_function
from six.moves import range

from rstbx.dials_core.integration_core import show_observations
import iotbx.phil
from iotbx import data_plots
from cctbx.array_family import flex
from cctbx import miller
from cctbx.crystal import symmetry
from cctbx.sgtbx.bravais_types import bravais_lattice
from cctbx import uctbx
from libtbx.str_utils import format_value
from libtbx.utils import Usage, Sorry, multi_out
from libtbx import easy_pickle
from libtbx import adopt_init_args, group_args, Auto
from six.moves import cStringIO as StringIO
import os
import math
import time
import sys
import glob
from scitbx import matrix
op = os.path

from xfel.merging.database.merging_database import mysql_master_phil
master_phil="""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.  Repeat to \
    specify additional directories.
targlob = None
  .type = str
  .multiple = True
  .help = new feature, instead of data records giving directories containing integration pickles
  .help = give a single blob giving the paths of tar files where the pickles are packaged up.
  .help = This reduces the number of files to be read in.  But as currently implemented
  .help = it does not reduce the number of file opens.
hash_filenames = False
  .type = bool
  .help = For CC1/2, instead of using odd/even filenames to split images into two sets, \
          hash the filename using md5 and split the images using odd/even hashes.
predictions_to_edge {
  apply = False
    .type = bool
    .help = If True and key 'indices_to_edge' not found in integration pickles, predictions
    .help = will be made to the edge of the detector based on current unit cell, orientation,
    .help = and mosaicity.
  image = None
    .type = path
    .help = Path to an example image from which to extract active areas and pixel size.
  detector_phil = None
    .type = path
    .help = Path to the detector version phil file used to generate the selected data.
}
short_circuit = False
  .type = bool
  .help = Assess per-image resolution limits and exit early.
a_list = None
  .type = path
  .multiple = False # for now XXX possibly make it multiple later
  .help = Text file containing a list of acceptable integration pickles, that is, not ones
  .help = that are misindexed, wrong type, or otherwise rejected as determined separately by user
filename_extension = "pickle"
  .type = str
  .help = Filename extension for integration pickle files. Usually pickle but can be otherwise.
data_subset = 0
  .type = int
  .help = 0: use all data / 1: use odd-numbered frames / 2: use even-numbered frames
validation {
  exclude_CSPAD_sensor = None
    .type = int
    .help = Index in range(32) of a sensor on the CSPAD to exclude from merging, for the purposes
    .help = of testing whether an individual sensor is poorly calibrated.
}
model = None
  .type = str
  .help = PDB filename containing atomic coordinates & isomorphous cryst1 record
  .help = or MTZ filename from a previous cycle of cxi.merge (not yet tested with prime MTZ)
model_reindex_op = h,k,l
  .type = str
  .help = Kludge for cases with an indexing ambiguity, need to be able to adjust scaling model
data_reindex_op = h,k,l
  .type = str
  .help = Reindex, e.g. to change C-axis of an orthorhombic cell to align Bravais lattice from indexing with actual space group
target_unit_cell = None
  .type = unit_cell
  .help = leave as None, program uses the model PDB file cryst1 record
target_space_group = None
  .type = space_group
  .help = leave  as None, program uses the model PDB file cryst1 record
set_average_unit_cell = False
  .type = bool
rescale_with_average_cell = False
  .type = bool
d_min = None
  .type = float
  .help = limiting resolution for scaling and merging
d_max = None
  .type = float
  .help = limiting resolution for scaling and merging.  Implementation currently affects only the CCiso cal
k_sol = 0.35
  .type = float
  .help = bulk solvent scale factor - approximate mean value in PDB \
    (according to Pavel)
b_sol = 46.00
  .type = float
  .help = bulk solvent B-factor - approximate mean value in PDB \
    (according to Pavel)
merge_anomalous = False
  .type = bool
  .help = Merge anomalous contributors
include_bulk_solvent = True
  .type = bool
wavelength = None
  .type = float
pixel_size = None
  .type = float
  .help = Detector-specific parameter for pixel size in mm
elements = None
  .type = str
  .multiple = True
significance_filter {
  apply = True
    .type = bool
    .help = Apply a sigma cutoff (on unmerged data) to limit resolution from each diffraction pattern
  n_bins = 12
    .type = int (value_min=2)
    .help = Initial target number of resolution bins for sigma cutoff calculation
  min_ct = 10
    .type = int
    .help = Decrease number of resolution bins to require mean bin population >= min_ct
  max_ct = 50
    .type = int
    .help = Increase number of resolution bins to require mean bin population <= max_ct
  sigma = 0.5
    .type = float
    .help = Remove highest resolution bins such that all accepted bins have <I/sigma> >= sigma
}
min_corr = 0.1
  .type = float
  .help = Correlation cutoff for rejecting individual frames.
  .help = This filter is not applied if model==None.
unit_cell_length_tolerance = 0.1
  .type = float
  .help = Fractional change in unit cell dimensions allowed (versus target \
    cell).
unit_cell_angle_tolerance = 2.
  .type = float
nproc = None
  .type = int
raw_data {
  sdfac_auto = False
    .type = bool
    .help = apply sdfac to each-image data assuming negative intensities are normally distributed noise
  sdfac_refine = False
    .type = bool
    .help = Correct merged sigmas by refining sdfac, sdb and sdadd according to Evans 2011. Not \
            compatible with sdfac_auto.
  errors_from_sample_residuals = False
    .type = bool
    .help = Use sample residuals as error estimates. Not compatible with sdfac_auto or sdfac_refine.
  propagate_errors = False
    .type = bool
    .help = Propagate errors from estimated parameters
  error_models {
    sdfac_refine {
      random_seed = None
        .help = Random seed. May be int or None. Only used for the simplex minimizer
        .type = int
        .expert_level = 1
      minimizer = *simplex lbfgs LevMar
        .type = choice
        .help = Which minimizer to use while refining the Sdfac terms
      refine_propagated_errors = False
        .type = bool
        .help = If True and if propagate_errors is True, then during sdfac refinement, also \
                refine the estimated error used for error propagation.
      show_finite_differences = False
        .type = bool
        .help = If True and minimizer is lbfgs, show the finite vs. analytical differences
      plot_refinement_steps = False
        .type = bool
        .help = If True, plot refinement steps during refinement.
      apply_to_I_only = False
        .type = bool
        .expert_level = 4
        .help = If True, after sdfac refinement, apply the new errors only to the merged I \
                during the weighted sum, keeping the old errors for the sigI value.        \
                Flag for research only, not recommended for general use.
    }
  }
}
output {
  n_bins = 10
    .type = int
    .help = Number of resolution bins in statistics table
  prefix = iobs
    .type = str
    .help = Prefix for all output file names
  title = None
    .type = str
    .help = Title for run - will appear in MTZ file header
  unit_cell = None
    .type = unit_cell
    .help = override the average or reference unit cell
  space_group = None
    .type = space_group
    .help = override the identified space group
}
merging {
  refine_G_Imodel = False
    .type = bool
    .help = "Refine per-frame scaling factors and model intensities
             after initial merging"
  reverse_lookup = None
    .type = str
    .help = filename, pickle format, generated by the cxi.brehm_diederichs program.  Contains a
    .help = (key,value) dictionary where key is the filename of the integrated data pickle file (supplied
    .help = with the data phil parameter and value is the h,k,l reindexing operator that resolves the
    .help = indexing ambiguity.
  minimum_multiplicity = None
    .type = int
    .help = If defined, merged structure factors not produced for the Miller indices below this threshold.
}
scaling {
  mtz_file = None
    .type = str
    .help = for Riso/ CCiso, the reference structure factors, must have data type F
    .help = a fake file is written out to this file name if model is None
  mtz_column_F = fobs
    .type = str
    .help = for Riso/ CCiso, the column name containing reference structure factors
  log_cutoff = None
    .type = float
    .help = for CC calculation, log(intensity) cutoff, ignore values less than this
  show_plots = False
    .type = bool
  enable = True
    .type = bool
    .help = enable the mark0 algorithm, otherwise individual-image scale factors are set to 1.0
    .expert_level = 3
  algorithm = *mark0 mark1 levmar
    .type = choice
    .help = "mark0: original per-image scaling by reference to
             isomorphous PDB model"
    .help = "mark1: no scaling, just averaging (i.e. Monte Carlo
             algorithm).  Individual image scale factors are set to 1."
    .help = "Sparse-matrix Levenberg-Marquardt scaling and merging.  Under development, not available in cxi.merge"
  simulation = None
    .type = str
    .help = To test scaling, use simulated data from the model instead of the actual observations
    .help = String value governs how the sim data are calculated
  simulation_data = None
    .type = floats
    .help = Extra parameters for the simulation, exact meaning depends on calculation method
  report_ML = False
    .type = bool
    .help = Report statistics on per-frame attributes modeled by max-likelihood fit (expert only)
}
postrefinement {
  enable = False
    .type = bool
    .help = enable the preliminary postrefinement algorithm (monochromatic)
    .expert_level = 3
  algorithm = *rs rs2 rs_hybrid eta_deff
    .type = choice
    .help = rs only, eta_deff protocol 7
    .expert_level = 3
  rs2
    .help = Reimplement postrefinement with the following (Oct 2016):
    .help = Refinement engine now work on analytical derivatives instead of finite differences
    .help = Better convergence using "traditional convergence test"
    .help = Use a streamlined frame_db schema, currently only supported for FS (filesystem) backend
    {}
  rs_hybrid
    .help = More aggressive postrefinement with the following (Oct 2016):
    .help = One round of 'rs2' using LBFGS minimizer as above to refine G,B,rotx,roty
    .help = Gentle weighting rather than unit weighting for the postrefinement target
    .help = Second round of LevMar adding an Rs refinement parameter
    .help = Option of weighting the merged terms by partiality
    {
    partiality_threshold = 0.2
      .type = float ( value_min = 0.01 )
      .help = throw out observations below this value. Hard coded as 0.2 for rs2, allow value for hybrid
      .help = must enforce minimum positive value because partiality appears in the denominator
    }
  target_weighting = *unit variance gentle extreme
    .type = choice
    .help = weights for the residuals in the postrefinement target (for rs2 or rs_hybrid)
    .help = Unit: each residual weighted by 1.0
    .help = Variance: weighted by 1/sigma**2.  Doesn't seem right, constructive feedback invited
    .help = Gentle: weighted by |I|/sigma**2.  Seems like best option
    .help = Extreme: weighted by (I/sigma)**2.  Also seems right, but severely downweights weak refl
  merge_weighting = *variance
    .type = choice
    .help = assumed that individual reflections are weighted by the counting variance
  merge_partiality_exponent = 0
    .type = float
    .help = additionally weight each measurement by partiality**exp when merging
    .help = 0 is no weighting, 1 is partiality weighting, 2 is weighting by partiality-squared
  lineshape = *lorentzian gaussian
    .type = choice
    .help = Soft sphere RLP modeled with Lorentzian radial profile as in prime
    .help = or Gaussian radial profile. (for rs2 or rs_hybrid)
  show_trumpet_plot = False
    .type = bool
    .help = each-image trumpet plot showing before-after plot. Spot color warmth indicates I/sigma
    .help = Spot radius for lower plot reflects partiality. Only implemented for rs_hybrid
}
include_negatives = False
  .type = bool
  .help = Whether to include negative intensities during scaling and merging
include_negatives_fix_27May2018 = True
  .type = bool
  .help = Bugfix for include negatives. Affects cxi.xmerge.
plot_single_index_histograms = False
  .type = bool
data_subsubsets {
  subsubset = None
    .type = int
  subsubset_total = None
    .type = int
}
isoform_name = None
  .type = str
  .help = Only accept this isoform
memory {
  shared_array_allocation = None
    .type = int
    .help = Yes, it's true. Python shared arrays don't seem to be dynamic in size,
    .help = so we have to allocate them just like in Fortran.  Take your best guess
    .help = as to how many total measurements you are joining.  Only for the levmar
    .help = algorithm for now. Insufficient allocation generates
    .help = ValueError: Can only assign sequence of same size
}
levmar {
  compute_cc_half = True
    .type = bool
    .help = Double the work, but only way to get reportable results. False for debug.
  sdfac_value = 1.0
    .type = float
    .help = Multiply all input sigmas by a constant factor so as to adjust the final
    .help = values of chisq/dof to about 1.0.  Still trial and error after 30 years.
  termination
    .help = Adjust the termination criteria for LevMar non-linear least squares refinement.
    {
    step_threshold = 0.0001
      .type = float
      .help = threshold for ||step||/||parameters||.  Very sensitive, 0.0001 is comprehensive
      .help = refinement, 0.001 is short refinement, 0.01 is abortive.
    objective_decrease_threshold = None
      .type = float
      .help = threshold for the fractional decrease of the objective function. Convenient
      .help = parameter for fine tuning very large parameter optimizations.  1.E-7 is a
      .help = reasonably comprehensive refinement, 1.E-4 is a good (short) choice for >10E5 parameters.
    }
  parameter_flags = PartialityDeff PartialityEtaDeff Bfactor Deff Eta Rxy
    .help = choices for the refinement
    .help = PartialityDeff and PartialityEtaDeff are mutually exclusize mosaic models (or don't use)
    .help = Others are refineable parameters to choose from, default is only refine G and I.
    .type = choice
    .multiple = True
}
lattice_rejection {
  unit_cell = Auto
    .type = unit_cell
    .help = unit_cell for filtering crystals with the given unit cell params. If Auto will automatically choose PDB model unit cell.
  space_group = Auto
    .type = space_group
    .help = space_group for filtering crystals with the given space group params. If Auto will automatically choose PDB space group.
  d_min = None
    .type = float
    .help = minimum resolution for lattices to be merged
}
""" + mysql_master_phil

def get_observations (work_params):
  try:
    data_dirs = work_params.data
    data_subset = work_params.data_subset
    subsubset = work_params.data_subsubsets.subsubset
    subsubset_total = work_params.data_subsubsets.subsubset_total
    extension = work_params.filename_extension
  except Exception as e:
    exit("Changed the interface for get_observations, please contact authors "+str(e))

  print("Step 1.  Get a list of all files")
  if work_params.a_list is not None:
    permissible_file_names = [a.strip() for a in open(work_params.a_list,"r").readlines()]
    permissible_file_hash = dict( zip(permissible_file_names, [None]*len(permissible_file_names)) )
    n_sorry = 0
  file_names = []
  if work_params.targlob:
    tar_list = [tar for tg in work_params.targlob for tar in glob.glob(tg)]
    for tarname in tar_list:
      import tarfile
      T = tarfile.open(name=tarname, mode='r')
      K = T.getmembers()
      NT = len(K)
      for nt in range(NT):
        k = os.path.basename(K[nt].path)
        file_names.append("%s;member%05d;timestamp%s"%(tarname,nt,k))
      print(tarname,NT)
  else:
    for dir_name in data_dirs :
      if not os.path.isdir(dir_name):
        if os.path.isfile(dir_name):
          #check if list-of-pickles text file is given
          pickle_list_file = open(dir_name,'r')
          pickle_list = pickle_list_file.read().split("\n")
        else:
          pickle_list = glob.glob(dir_name)
        for pickle_filename in pickle_list:
          if work_params.a_list is not None and pickle_filename not in permissible_file_hash:
            # use A_list mechanism to reject files not on the "acceptable" list
            #print "SORRY--%s FILE NOT ON THE A-List"%(pickle_filename)
            n_sorry+=1
            continue
          if os.path.isfile(pickle_filename) and pickle_filename.endswith("."+extension):
            if data_subset==0 or \
              (data_subset==1 and (int(os.path.basename(pickle_filename).split("."+extension)[0][-1])%2==1)) or \
              (data_subset==2 and (int(os.path.basename(pickle_filename).split("."+extension)[0][-1])%2==0)):
              file_names.append(pickle_filename)
        continue
      for file_name in os.listdir(dir_name):
        if work_params.a_list is not None and os.path.join(dir_name, file_name) not in permissible_file_hash:
          # use A_list mechanism to reject files not on the "acceptable" list
          print("SORRY--%s FILE NOT ON THE A-List"%(os.path.join(dir_name, file_name)))
          n_sorry+=1
          continue
        if (file_name.endswith("_00000."+extension)):
          if data_subset==0 or \
            (data_subset==1 and (int(os.path.basename(file_name).split("_00000."+extension)[0][-1])%2==1)) or \
            (data_subset==2 and (int(os.path.basename(file_name).split("_00000."+extension)[0][-1])%2==0)):
            file_names.append(os.path.join(dir_name, file_name))
        elif (file_name.endswith("."+extension)):
          if data_subset==0 or \
            (data_subset==1 and (int(os.path.basename(file_name).split("."+extension)[0][-1])%2==1)) or \
            (data_subset==2 and (int(os.path.basename(file_name).split("."+extension)[0][-1])%2==0)):
            file_names.append(os.path.join(dir_name, file_name))
    if work_params.a_list is not None:
      print ("A_LIST: %d names rejected for not being on the a_list, leaving %d accepted"%(n_sorry, len(file_names)))
    if subsubset is not None and subsubset_total is not None:
      file_names = [file_names[i] for i in range(len(file_names)) if (i+subsubset)%subsubset_total == 0]
  print("Number of pickle files found:", len(file_names))
  print()
  return file_names

class WrongBravaisError (Exception) :
  pass

class OutlierCellError (Exception) :
  pass

def get_boundaries_from_sensor_ID (sensor_ID,
                                   detector_version_phil=None,
                                   image_with_header=None) :
  if detector_version_phil is not None and image_with_header is not None:
    from xfel.command_line.cctbx_integration_pickle_viewer import get_CSPAD_active_areas
    active_areas = get_CSPAD_active_areas(image_with_header, detector_version_phil)
  else:
    from xfel.command_line.cctbx_integration_pickle_viewer import LG36_active_areas
    active_areas = LG36_active_areas
  return active_areas[8*sensor_ID:8*sensor_ID+8]

def load_result (file_name,
                 ref_bravais_type,
                 reference_cell,
                 params,
                 reindex_op,
                 out,
                 get_predictions_to_edge=False,
                 image_info=None,
                 exclude_CSPAD_sensor=None) :
  # If @p file_name cannot be read, the load_result() function returns
  # @c None.

  print("-" * 80, file=out)
  print("Step 2.  Load pickle file into dictionary obj and filter on lattice & cell with",reindex_op, file=out)
  print(file_name, file=out)
  """
  Take a pickle file, confirm that it contains the appropriate data, and
  check the lattice type and unit cell against the reference settings - if
  rejected, raises an exception (for tracking statistics).
  """
  # Ignore corrupted pickle files.
  if params.targlob:
    file_name,imember = file_name.split(";member")
    imember,timestamp = imember.split(";timestamp")
    import sys,tarfile
    T = tarfile.open(name=file_name, mode='r')
    K = T.getmembers()
    this_member = K[int(imember)]
    fileIO = T.extractfile(member=this_member)
    from six.moves import cPickle as pickle
    try:
      obj = pickle.load(fileIO)
    except Exception:
      return None
  else:
    try:
      obj = easy_pickle.load(file_name=file_name)
    except Exception:
      return None
  if ("observations" not in obj) :
    return None
  if params.isoform_name is not None:
    if not "identified_isoform" in obj:
      return None
    if obj["identified_isoform"] != params.isoform_name:
      return None
  if reindex_op == "h,k,l":
    pass
    #raise WrongBravaisError("Skipping file with h,k,l")
  else:
    for idx in range(len(obj["observations"])):
      obj["observations"][idx] = obj["observations"][idx].change_basis(reindex_op)
    from cctbx.sgtbx import change_of_basis_op
    cb_op = change_of_basis_op(reindex_op)
    ORI = obj["current_orientation"][0]
    Astar = matrix.sqr(ORI.reciprocal_matrix())
    CB_OP_sqr = matrix.sqr(cb_op.c().r().as_double()).transpose()
    from cctbx.crystal_orientation import crystal_orientation
    ORI_new = crystal_orientation((Astar*CB_OP_sqr).elems, True)
    obj["current_orientation"][0] = ORI_new
    # unexpected, first since cxi_merge had a previous postrefinement implementation that
    #  apparently failed to apply the reindex_op (XXX still must be fixed), and secondly
    #  that the crystal_orientation.change basis() implementation fails to work because
    #  it does not use the transpose:  obj["current_orientation"][0].change_basis(cb_op)
    pass
    #raise WrongBravaisError("Skipping file with alternate indexing %s"%reindex_op)

  result_array = obj["observations"][0]
  unit_cell = result_array.unit_cell()
  sg_info = result_array.space_group_info()
  print("", file=out)

  print(sg_info, file=out)
  print(unit_cell, file=out)

  if obj.get('beam_s0',None) is not None:
    # Remove the need for pixel size within cxi.merge.  Allows multipanel detector with dissimilar panels.
    # Relies on new frame extractor code called by dials.stills_process that writes s0, s1 and polarization normal
    # vectors all to the integration pickle.  Future path: use dials json and reflection file.
    s0_vec = matrix.col(obj["beam_s0"]).normalize()
    s0_polar_norm = obj["beam_polarization_normal"]
    s1_vec = obj["s1_vec"][0]
    Ns1 = len(s1_vec)
    # project the s1_vector onto the plane normal to s0.  Get result by subtracting the
    # projection of s1 onto s0, which is (s1.dot.s0_norm)s0_norm
    s0_norm = flex.vec3_double(Ns1,s0_vec)
    s1_proj = (s1_vec.dot(s0_norm))*s0_norm
    s1_in_normal_plane = s1_vec - s1_proj
    # Now want the polar angle between the projected s1 and the polarization normal
    s0_polar_norms = flex.vec3_double(Ns1,s0_polar_norm)
    dotprod = (s1_in_normal_plane.dot(s0_polar_norms))
    costheta = dotprod/(s1_in_normal_plane.norms())
    theta = flex.acos(costheta)
    prospective = flex.cos(2.0*theta)
    obj["cos_two_polar_angle"] = prospective
    # gives same as old answer to ~1% but not exact.  Not sure why, should not matter.

  else:
    #Check for pixel size (at this point we are assuming we have square pixels, all experiments described in one
    #refined_experiments.json file use the same detector, and all panels on the detector have the same pixel size)
    if params.pixel_size is not None:
      pixel_size = params.pixel_size
    elif "pixel_size" in obj:
      pixel_size = obj["pixel_size"]
    else:
      raise Sorry("Cannot find pixel size. Specify appropriate pixel size in mm for your detector in phil file.")

    #Calculate displacements based on pixel size
    assert obj['mapped_predictions'][0].size() == obj["observations"][0].size()
    mm_predictions = pixel_size*(obj['mapped_predictions'][0])
    mm_displacements = flex.vec3_double()
    cos_two_polar_angle = flex.double()
    for pred in mm_predictions:
      mm_displacements.append((pred[0]-obj["xbeam"],pred[1]-obj["ybeam"],0.0))
      cos_two_polar_angle.append( math.cos( 2. * math.atan2(pred[1]-obj["ybeam"],pred[0]-obj["xbeam"]) ) )
    obj["cos_two_polar_angle"] = cos_two_polar_angle
    #then convert to polar angle and compute polarization correction

  if (not bravais_lattice(sg_info.type().number()) == ref_bravais_type) :
    raise WrongBravaisError("Skipping cell in different Bravais type (%s)" %
      str(sg_info))
  if (not unit_cell.is_similar_to(
      other=reference_cell,
      relative_length_tolerance=params.unit_cell_length_tolerance,
      absolute_angle_tolerance=params.unit_cell_angle_tolerance)) :
    raise OutlierCellError(
      "Skipping cell with outlier dimensions (%g %g %g %g %g %g" %
      unit_cell.parameters())
  # Illustrate how a unit cell filter would be implemented.
  #ucparams = unit_cell.parameters()
  #if not (130.21 < ucparams[2] < 130.61) or not (92.84 < ucparams[0] < 93.24):
  #  print >> out, "DOES NOT PASS ERSATZ UNIT CELL FILTER"
  #  return None
  print("Integrated data:", file=out)
  result_array.show_summary(f=out, prefix="  ")
  # XXX don't force reference setting here, it will be done later, after the
  # original unit cell is recorded

  #Remove observations on a selected sensor if requested
  if exclude_CSPAD_sensor is not None:
    print("excluding CSPAD sensor %d" % exclude_CSPAD_sensor, file=out)
    fast_min1, slow_min1, fast_max1, slow_max1, \
    fast_min2, slow_min2, fast_max2, slow_max2 = \
      get_boundaries_from_sensor_ID(exclude_CSPAD_sensor)
    fast_min = min(fast_min1, fast_min2)
    slow_min = min(slow_min1, slow_min2)
    fast_max = max(fast_max1, fast_max2)
    slow_max = max(slow_max1, slow_max2)
    accepted = flex.bool()
    px_preds = obj['mapped_predictions'][0]
    for idx in range(px_preds.size()):
      pred_fast, pred_slow = px_preds[idx]
      if (pred_fast < fast_min) or (pred_fast > fast_max) or (pred_slow < slow_min) or (pred_slow > slow_max):
        accepted.append(True)
      else:
        accepted.append(False)
    obj['mapped_predictions'][0] = obj['mapped_predictions'][0].select(accepted)
    obj['observations'][0] = obj['observations'][0].select(accepted)
    obj['cos_two_polar_angle'] = obj["cos_two_polar_angle"].select(accepted)
  if not 'indices_to_edge' in obj.keys():
    if get_predictions_to_edge:
      from xfel.merging.predictions_to_edges import extend_predictions
      extend_predictions(obj, file_name, image_info, dmin=1.5, dump=False)
    else:
      obj['indices_to_edge'] = None
  return obj

from xfel.cxi.merging_utils import intensity_data, frame_data, null_data

class unit_cell_distribution (object) :
  """
  Container for collecting unit cell edge length statistics - both for frames
  included in the final dataset, and those rejected due to poor correlation.
  (Frames with incompatible indexing solutions will not be included.)
  """
  # TODO make this more general - currently assumes that angles are fixed,
  # which is true for the systems studied so far
  def __init__ (self) :
    self.uc_a_values = flex.double()
    self.uc_b_values = flex.double()
    self.uc_c_values = flex.double()
    self.all_uc_a_values = flex.double()
    self.all_uc_b_values = flex.double()
    self.all_uc_c_values = flex.double()

  def add_cell (self, unit_cell, rejected=False) :
    if (unit_cell is None) :
      return
    (a,b,c,alpha,beta,gamma) = unit_cell.parameters()
    if (not rejected) :
      self.uc_a_values.append(a)
      self.uc_b_values.append(b)
      self.uc_c_values.append(c)
    self.all_uc_a_values.append(a)
    self.all_uc_b_values.append(b)
    self.all_uc_c_values.append(c)

  def add_cells(self, uc) :
    """Addition operation for unit cell statistics."""
    self.uc_a_values.extend(uc.uc_a_values)
    self.uc_b_values.extend(uc.uc_b_values)
    self.uc_c_values.extend(uc.uc_c_values)
    self.all_uc_a_values.extend(uc.all_uc_a_values)
    self.all_uc_b_values.extend(uc.all_uc_b_values)
    self.all_uc_c_values.extend(uc.all_uc_c_values)

  def show_histograms (self, reference, out, n_slots=20) :
    [a0,b0,c0,alpha0,beta0,gamma0] = reference.parameters()
    print("", file=out)
    labels = ["a","b","c"]
    ref_edges = [a0,b0,c0]
    def _show_each (edges) :
      for edge, ref_edge, label in zip(edges, ref_edges, labels) :
        h = flex.histogram(edge, n_slots=n_slots)
        smin, smax = flex.min(edge), flex.max(edge)
        stats = flex.mean_and_variance(edge)
        print("  %s edge" % label, file=out)
        print("     range:     %6.2f - %.2f" % (smin, smax), file=out)
        print("     mean:      %6.2f +/- %6.2f on N = %d" % (
          stats.mean(), stats.unweighted_sample_standard_deviation(), edge.size()), file=out)
        print("     reference: %6.2f" % ref_edge, file=out)
        h.show(f=out, prefix="    ", format_cutoffs="%6.2f")
        print("", file=out)
    edges = [self.all_uc_a_values, self.all_uc_b_values, self.all_uc_c_values]
    print("Unit cell length distribution (all frames with compatible indexing):", file=out)
    _show_each(edges)
    edges = [self.uc_a_values, self.uc_b_values, self.uc_c_values]
    print("Unit cell length distribution (frames with acceptable correlation):", file=out)
    _show_each(edges)

  def get_average_cell_dimensions (self) :
    a = flex.mean(self.uc_a_values)
    b = flex.mean(self.uc_b_values)
    c = flex.mean(self.uc_c_values)
    return a,b,c

#-----------------------------------------------------------------------
class scaling_manager (intensity_data) :
  def __init__ (self, miller_set, i_model, params, log=None) :
    if (log is None) :
      log = sys.stdout
    self.log = log
    self.params = params
    self.miller_set = miller_set
    self.i_model = i_model
    self.ref_bravais_type = bravais_lattice(
      miller_set.space_group_info().type().number())
    intensity_data.__init__(self, miller_set.size())
    self.reverse_lookup = None
    if params.merging.reverse_lookup is not None:
      self.reverse_lookup = easy_pickle.load(params.merging.reverse_lookup)
    self.reset()

  def reset (self) :
    self.n_processed = 0
    self.n_accepted = 0
    self.n_file_error = 0
    self.n_low_signal = 0
    self.n_wrong_bravais = 0
    self.n_wrong_cell = 0
    self.n_low_resolution = 0
    self.n_low_corr = 0
    self.failure_modes = {}
    self.observations = flex.int()
    self.predictions_to_edge = flex.int()
    self.corr_values = flex.double()
    self.rejected_fractions = flex.double()
    self.uc_values = unit_cell_distribution()
    self.d_min_values = flex.double()
    self.wavelength = flex.double()
    self.initialize()

  @staticmethod
  def single_reflection_histograms(obs, ISIGI):
    # Per-bin sum of I and I/sig(I) for each observation.
    for i in obs.binner().array_indices(i_bin) :
      import numpy as np
      index = obs.indices()[i]
      if (index in ISIGI) :
        # Compute m, the "merged" intensity, as the average intensity
        # of all observations of the reflection with the given index.
        N = 0
        m = 0
        for t in ISIGI[index] :
          N += 1
          m += t[0]
          print("Miller %20s n-obs=%4d  sum-I=%10.0f"%(index, N, m))
          plot_n_bins = N//10
          hist,bins = np.histogram([t[0] for t in ISIGI[index]],bins=25)
          width = 0.7*(bins[1]-bins[0])
          center = (bins[:-1]+bins[1:])/2
          import matplotlib.pyplot as plt
          plt.bar(center, hist, align="center", width=width)
          plt.show()

  def scale_all (self, file_names) :
    t1 = time.time()
    if self.params.backend == 'MySQL':
      from xfel.merging.database.merging_database import manager
    elif self.params.backend == 'SQLite':
      from xfel.merging.database.merging_database_sqlite3 import manager
    else:
      from xfel.merging.database.merging_database_fs import manager

    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())

    # Unless the number of requested processes is greater than one,
    # try parallel multiprocessing on a parallel host.  Block until
    # all database commands have been processed.
    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto):
      import libtbx.introspection
      nproc = libtbx.introspection.number_of_processors()
    if nproc > 1:
      try :
        import multiprocessing
        self._scale_all_parallel(file_names, db_mgr)
      except ImportError as e :
        print("multiprocessing module not available (requires Python >= 2.6)\n" \
          "will scale frames serially", file=self.log)
        self._scale_all_serial(file_names, db_mgr)
    else:
      self._scale_all_serial(file_names, db_mgr)
    db_mgr.join()

    t2 = time.time()
    print("", file=self.log)
    print("#" * 80, file=self.log)
    print("FINISHED MERGING", file=self.log)
    print("  Elapsed time: %.1fs" % (t2 - t1), file=self.log)
    print("  %d of %d integration files were accepted" % (
      self.n_accepted, len(file_names)), file=self.log)
    print("  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais, file=self.log)
    print("  %d rejected for unit cell outliers" % \
      self.n_wrong_cell, file=self.log)
    print("  %d rejected for low resolution" % \
      self.n_low_resolution, file=self.log)
    print("  %d rejected for low signal" % \
      self.n_low_signal, file=self.log)
    print("  %d rejected due to up-front poor correlation under min_corr parameter" % \
      self.n_low_corr, file=self.log)
    print("  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error, file=self.log)
    for key in self.failure_modes.keys():
      print("  %d rejected due to %s"%(self.failure_modes[key], key), file=self.log)

    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_corr + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell \
               + self.n_low_resolution \
               + sum([val for val in self.failure_modes.itervalues()])
    assert checksum == len(file_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print("Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min), file=self.log)

    if self.params.raw_data.propagate_errors and not self.params.raw_data.error_models.sdfac_refine.refine_propagated_errors:
      assert self.params.postrefinement.enable
      from xfel.merging.algorithms.error_model.sdfac_propagate import sdfac_propagate
      error_modeler = sdfac_propagate(self)
      error_modeler.adjust_errors()

    if self.params.raw_data.sdfac_refine or self.params.raw_data.errors_from_sample_residuals or \
        self.params.raw_data.error_models.sdfac_refine.refine_propagated_errors:
      if self.params.raw_data.sdfac_refine:
        if self.params.raw_data.error_models.sdfac_refine.minimizer == 'simplex':
          from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine as error_modeler
        elif self.params.raw_data.error_models.sdfac_refine.minimizer == 'lbfgs':
          if self.params.raw_data.error_models.sdfac_refine.refine_propagated_errors:
            from xfel.merging.algorithms.error_model.sdfac_propagate_and_refine import sdfac_propagate_and_refine as error_modeler
          else:
            from xfel.merging.algorithms.error_model.sdfac_refine_lbfgs import sdfac_refine_refltable_lbfgs as error_modeler
        elif self.params.raw_data.error_models.sdfac_refine.minimizer == 'LevMar':
          from xfel.merging.algorithms.error_model.sdfac_refine_levmar import sdfac_refine_refltable_levmar as error_modeler

      if self.params.raw_data.errors_from_sample_residuals:
        from xfel.merging.algorithms.error_model.errors_from_residuals import errors_from_residuals as error_modeler

      error_modeler(self).adjust_errors()

  def _scale_all_parallel (self, file_names, db_mgr) :
    import multiprocessing
    import libtbx.introspection

    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto) :
      nproc = libtbx.introspection.number_of_processors()

    # Input files are supplied to the scaling processes on demand by
    # means of a queue.
    #
    # XXX The input queue may need to either allow non-blocking
    # put():s or run in a separate process to prevent the procedure
    # from blocking here if the list of file paths does not fit into
    # the queue's buffer.
    input_queue = multiprocessing.Manager().JoinableQueue()
    for file_name in file_names:
      input_queue.put(file_name)

    pool = multiprocessing.Pool(processes=nproc)
    # Each process accumulates its own statistics in serial, and the
    # grand total is eventually collected by the main process'
    # _add_all_frames() function.
    for i in range(nproc) :
      sm = scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[input_queue, db_mgr],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

    # Block until the input queue has been emptied.
    input_queue.join()


  def _scale_all_serial (self, file_names, db_mgr) :
    """
    Scale frames sequentially (single-process).  The return value is
    picked up by the callback.
    """
    for file_name in file_names :
      scaled = self.scale_frame(file_name, db_mgr)
      if (scaled is not None) :
        self.add_frame(scaled)
    return (self)

  def add_frame (self, data) :
    """
    Combine the scaled data from a frame with the current overall dataset.
    Also accepts None or null_data objects, when data are unusable but we
    want to record the file as processed.
    """
    self.n_processed += 1
    if (data is None) :
      return None
    #data.show_log_out(self.log)
    #self.log.flush()
    if (isinstance(data, null_data)) :
      if (data.file_error) :
        self.n_file_error += 1
      elif (data.low_signal) :
        self.n_low_signal += 1
      elif (data.wrong_bravais) :
        self.n_wrong_bravais += 1
      elif (data.wrong_cell) :
        self.n_wrong_cell += 1
      elif (data.low_resolution) :
        self.n_low_resolution += 1
      elif (data.low_correlation) :
        self.n_low_corr += 1
      elif (getattr(data,"reason",None) is not None):
        if str(data.reason)!="":
          self.failure_modes[str(data.reason)] = self.failure_modes.get(str(data.reason),0) + 1
        elif repr(type(data.reason))!="":
          self.failure_modes[repr(type(data.reason))] = self.failure_modes.get(repr(type(data.reason)),0) + 1
        else:
          self.failure_modes["other reasons"] = self.failure_modes.get("other reasons",0) + 1
      return
    if (data.accept) :
      self.n_accepted    += 1
      self.completeness  += data.completeness
      #print "observations count increased by %d" % flex.sum(data.completeness)
      self.completeness_predictions += data.completeness_predictions
      #print "predictions count increased by %d" % flex.sum(data.completeness_predictions)
      self.summed_N      += data.summed_N
      self.summed_weight += data.summed_weight
      self.summed_wt_I   += data.summed_wt_I
      for index, isigi in data.ISIGI.iteritems() :
        if (index in self.ISIGI):
          self.ISIGI[index] += isigi
        else:
          self.ISIGI[index] = isigi
    else :
      self.n_low_corr += 1 # FIXME this is no longer the right default
    self.uc_values.add_cell(data.indexed_cell,
      rejected=(not data.accept))
    if not self.params.short_circuit:
      self.observations.append(data.n_obs)
    if (data.n_obs > 0) :
      frac_rejected = data.n_rejected / data.n_obs
      self.rejected_fractions.append(frac_rejected)
      self.d_min_values.append(data.d_min)
    self.corr_values.append(data.corr)
    self.wavelength.append(data.wavelength)

  def _add_all_frames (self, data) :
    """The _add_all_frames() function collects the statistics accumulated
    in @p data by the individual scaling processes in the process
    pool.  This callback function is run in serial, so it does not
    need a lock.
    """
    self.n_accepted += data.n_accepted
    self.n_file_error += data.n_file_error
    self.n_low_corr += data.n_low_corr
    self.n_low_signal += data.n_low_signal
    self.n_processed += data.n_processed
    self.n_wrong_bravais += data.n_wrong_bravais
    self.n_wrong_cell += data.n_wrong_cell
    self.n_low_resolution += data.n_low_resolution
    for key in data.failure_modes.keys():
      self.failure_modes[key] = self.failure_modes.get(key,0) + data.failure_modes[key]

    for index, isigi in data.ISIGI.iteritems() :
      if (index in self.ISIGI):
        self.ISIGI[index] += isigi
      else:
        self.ISIGI[index] = isigi

    self.completeness += data.completeness
    #print "observations count increased by %d" % flex.sum(data.completeness)
    self.completeness_predictions += data.completeness_predictions
    #print "predictions count increased by %d" % flex.sum(data.completeness_predictions)
    self.summed_N += data.summed_N
    self.summed_weight += data.summed_weight
    self.summed_wt_I += data.summed_wt_I

    self.corr_values.extend(data.corr_values)
    self.d_min_values.extend(data.d_min_values)
    if not self.params.short_circuit:
      self.observations.extend(data.observations)
    self.rejected_fractions.extend(data.rejected_fractions)
    self.wavelength.extend(data.wavelength)

    self.uc_values.add_cells(data.uc_values)

  def show_unit_cell_histograms (self) :
    self.uc_values.show_histograms(
      reference=self.miller_set.unit_cell(),
      out=self.log)

  def get_plot_statistics (self) :
    return plot_statistics(
      prefix=self.params.output.prefix,
      unit_cell_statistics=self.uc_values,
      reference_cell=self.miller_set.unit_cell(),
      correlations=self.corr_values,
      min_corr=self.params.min_corr,
      rejected_fractions=self.rejected_fractions,
      frame_d_min=self.d_min_values)

  def get_overall_correlation (self, sum_I) :
    """
    Correlate the averaged intensities to the intensities from the
    reference data set.  XXX The sum_I argument is really a kludge!
    """
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x  = 0
    sum_y  = 0
    N      = 0
    for i in range(len(self.summed_N)):
      if (self.summed_N[i] <= 0):
        continue
      # skip structure factor if i_model.sigma is invalid (intentionally < 0)
      if self.i_model.sigmas()[i] < 0: continue
      I_r       = self.i_model.data()[i]
      I_o       = sum_I[i]/self.summed_N[i]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
             math.sqrt(N * sum_yy - sum_y**2))
    print("SUMMARY: For %d reflections, got slope %f, correlation %f" \
        % (N, slope, corr), file=self.log)
    return N, corr

  def finalize_and_save_data (self) :
    """
    Assemble a Miller array with the summed data, setting the unit cell to
    the consensus average if desired, and write to an MTZ file (including
    merged/non-anomalous data too).
    """
    print("", file=self.log)
    print("#" * 80, file=self.log)
    print("OUTPUT FILES", file=self.log)
    if self.params.merging.minimum_multiplicity is None:
      multiplicity_flag = flex.bool(len(self.summed_N),True)
    else:
      multiplicity_flag = (self.summed_N > self.params.merging.minimum_multiplicity)
    Iobs_all = flex.double(self.miller_set.size())
    SigI_all = flex.double(self.miller_set.size())
    for i in range(len(Iobs_all)):
      if (self.summed_weight[i] > 0.):
       if (multiplicity_flag[i]):
        Iobs_all[i] = self.summed_wt_I[i] / self.summed_weight[i]
        if hasattr(self, 'summed_weight_uncorrected'):
          SigI_all[i] = math.sqrt(1. / self.summed_weight_uncorrected[i])
        else:
          SigI_all[i] = math.sqrt(1. / self.summed_weight[i])
    if (self.params.set_average_unit_cell) :
      # XXX since XFEL crystallography runs at room temperature, it may not
      # be appropriate to use the cell dimensions from a cryo structure.
      # also, some runs seem to have huge variance in the indexed cell
      # dimensions, so downstream programs (MR, refinement) may run better
      # with the cell set to the mean edge lengths.
      abc = self.uc_values.get_average_cell_dimensions()
      print("  (will use final unit cell edges %g %g %g)" % abc, file=self.log)
      angles = self.miller_set.unit_cell().parameters()[3:]
      unit_cell = uctbx.unit_cell(list(abc) + list(angles))
      final_symm = symmetry(
        unit_cell=unit_cell,
        space_group_info=self.miller_set.space_group_info())
    else :
      final_symm = self.miller_set
    if (self.params.output.unit_cell is not None or self.params.output.space_group is not None) :
      output_uc = self.params.output.unit_cell or final_symm.unit_cell()
      output_sg = self.params.output.space_group or final_symm.space_group_info()
      final_symm = symmetry(unit_cell=output_uc,
                            space_group_info=output_sg)
    all_obs = miller.array(
      miller_set=self.miller_set.customized_copy(
        crystal_symmetry=final_symm),
      data=Iobs_all,
      sigmas=SigI_all).resolution_filter(
        d_min=self.params.d_min).set_observation_type_xray_intensity()
    mtz_file = "%s.mtz" % self.params.output.prefix
    if not self.params.include_negatives:
      all_obs = all_obs.select(all_obs.data() > 0)

    mtz_out = all_obs.as_mtz_dataset(
      column_root_label="Iobs",
      title=self.params.output.title,
      wavelength=flex.mean(self.wavelength))
    mtz_out.add_miller_array(
      miller_array=all_obs.average_bijvoet_mates(),
      column_root_label="IMEAN")
    mtz_obj = mtz_out.mtz_object()
    mtz_obj.write(mtz_file)
    print("  Anomalous and mean data:\n    %s" % \
      os.path.abspath(mtz_file), file=self.log)
    print("", file=self.log)
    print("Final data:", file=self.log)
    all_obs.show_summary(self.log, prefix="  ")
    return mtz_file, all_obs

  def __call__ (self, input_queue, db_mgr) :
    # Scale frames sequentially within the current process.  The
    # return value is picked up by the callback.  See also
    # self.scale_all_serial()
    from Queue import Empty

    try :
      while True:
        try:
          file_name = input_queue.get_nowait()
        except Empty:
          return self

        scaled = self.scale_frame(file_name, db_mgr)
        if scaled is not None:
          self.add_frame(scaled)
        input_queue.task_done()

    except Exception as e :
      print(str(e), file=self.log)
      return None

  def scale_frame (self, file_name, db_mgr) :
    """The scale_frame() function populates a back end database with
    appropriately scaled intensities derived from a single frame.  The
    mark0 scaling algorithm determines the scale factor by correlating
    the frame's corrected intensities to those of a reference
    structure, while the mark1 algorithm applies no scaling at all.
    The scale_frame() function can be called either serially or via a
    multiprocessing map() function.

    @note This function must not modify any internal data or the
          parallelization will not yield usable results!

    @param file_name Path to integration pickle file
    @param db_mgr    Back end database manager
    @return          An intensity_data object
    """

    out = StringIO()
    wrong_cell = wrong_bravais = False
    reindex_op = self.params.data_reindex_op
    if self.reverse_lookup is not None:
      reindex_op = self.reverse_lookup.get(file_name, None)
      if reindex_op is None:
        return null_data(file_name=file_name, log_out=out.getvalue(), file_error=True)
    if (self.params.predictions_to_edge.apply) and (self.params.predictions_to_edge.image is not None):
      from xfel.merging.predictions_to_edges import ImageInfo
      image_info = ImageInfo(self.params.predictions_to_edge.image, detector_phil=self.params.predictions_to_edge.detector_phil)
    else:
      image_info = None

    if self.params.lattice_rejection.unit_cell == Auto:
      self.params.lattice_rejection.unit_cell = self.params.target_unit_cell

    try :
      result = load_result(
        file_name=file_name,
        reference_cell=self.params.lattice_rejection.unit_cell,
        ref_bravais_type=self.ref_bravais_type,
        params=self.params,
        reindex_op = reindex_op,
        out=out,
        get_predictions_to_edge=self.params.predictions_to_edge.apply,
        image_info=image_info,
        exclude_CSPAD_sensor=self.params.validation.exclude_CSPAD_sensor)
      if result is None:
        return null_data(
          file_name=file_name, log_out=out.getvalue(), file_error=True)
    except OutlierCellError as e :
      print(str(e), file=out)
      return null_data(
        file_name=file_name, log_out=out.getvalue(), wrong_cell=True)
    except WrongBravaisError as e :
      print(str(e), file=out)
      return null_data(
        file_name=file_name, log_out=out.getvalue(), wrong_bravais=True)
    return self.scale_frame_detail(result, file_name, db_mgr, out)

  def scale_frame_detail(self, result, file_name, db_mgr, out):
    # If the pickled integration file does not contain a wavelength,
    # fall back on the value given on the command line.  XXX The
    # wavelength parameter should probably be removed from master_phil
    # once all pickled integration files contain it.
    if ("wavelength" in result):
      wavelength = result["wavelength"]
    elif (self.params.wavelength is not None):
      wavelength = self.params.wavelength
    else:
      # XXX Give error, or raise exception?
      return None
    assert (wavelength > 0)

    observations = result["observations"][0]

    if len(observations.data()) == 0:
      print("skipping image: no observations", file=out)
      return null_data(
        file_name=file_name, log_out=out.getvalue(), low_signal=True)

    cos_two_polar_angle = result["cos_two_polar_angle"]
    indices_to_edge = result["indices_to_edge"]

    assert observations.size() == cos_two_polar_angle.size()
    tt_vec = observations.two_theta(wavelength)
    #print >> out, "mean tt degrees",180.*flex.mean(tt_vec.data())/math.pi
    cos_tt_vec = flex.cos( tt_vec.data() )
    sin_tt_vec = flex.sin( tt_vec.data() )
    cos_sq_tt_vec = cos_tt_vec * cos_tt_vec
    sin_sq_tt_vec = sin_tt_vec * sin_tt_vec
    P_nought_vec = 0.5 * (1. + cos_sq_tt_vec)

    F_prime = -1.0 # Hard-coded value defines the incident polarization axis
    P_prime = 0.5 * F_prime * cos_two_polar_angle * sin_sq_tt_vec
    # XXX added as a diagnostic
    prange=P_nought_vec - P_prime

    other_F_prime = 1.0
    otherP_prime = 0.5 * other_F_prime * cos_two_polar_angle * sin_sq_tt_vec
    otherprange=P_nought_vec - otherP_prime
    diff2 = flex.abs(prange - otherprange)
    print("mean diff is",flex.mean(diff2), "range",flex.min(diff2), flex.max(diff2), file=out)
    # XXX done
    observations = observations / ( P_nought_vec - P_prime )
    # This corrects observations for polarization assuming 100% polarization on
    # one axis (thus the F_prime = -1.0 rather than the perpendicular axis, 1.0)
    # Polarization model as described by Kahn, Fourme, Gadet, Janin, Dumas & Andre
    # (1982) J. Appl. Cryst. 15, 330-337, equations 13 - 15.

    print("Step 3. Correct for polarization.", file=out)
    indexed_cell = observations.unit_cell()

    observations_original_index = observations.deep_copy()
    if result.get("model_partialities",None) is not None and result["model_partialities"][0] is not None:
      # some recordkeeping useful for simulations
      partialities_original_index = observations.customized_copy(
        crystal_symmetry=self.miller_set.crystal_symmetry(),
        data = result["model_partialities"][0]["data"],
        sigmas = flex.double(result["model_partialities"][0]["data"].size()), #dummy value for sigmas
        indices = result["model_partialities"][0]["indices"],
        ).resolution_filter(d_min=self.params.d_min)

    assert len(observations_original_index.indices()) == len(observations.indices())

    # Now manipulate the data to conform to unit cell, asu, and space group
    # of reference.  The resolution will be cut later.
    # Only works if there is NOT an indexing ambiguity!
    if indices_to_edge is not None:
      predictions = self.miller_set.customized_copy(
        anomalous_flag=not self.params.merge_anomalous,
        crystal_symmetry=self.miller_set.crystal_symmetry(),
        indices=indices_to_edge
        ).map_to_asu()
    else:
      predictions = None

    observations = observations.customized_copy(
      anomalous_flag=not self.params.merge_anomalous,
      crystal_symmetry=self.miller_set.crystal_symmetry()
      ).map_to_asu()

    observations_original_index = observations_original_index.customized_copy(
      anomalous_flag=not self.params.merge_anomalous,
      crystal_symmetry=self.miller_set.crystal_symmetry()
      )
    print("Step 4. Filter on global resolution and map to asu", file=out)
    print("Data in reference setting:", file=out)
    #observations.show_summary(f=out, prefix="  ")
    show_observations(observations, out=out)

    if self.params.significance_filter.apply is True: #------------------------------------
      print("Step 5. Frame by frame resolution filter", file=out)
      # Apply an I/sigma filter ... accept resolution bins only if they
      #   have significant signal; tends to screen out higher resolution observations
      #   if the integration model doesn't quite fit
      N_obs_pre_filter = observations.size()
      N_bins_small_set = N_obs_pre_filter // self.params.significance_filter.min_ct
      N_bins_large_set = N_obs_pre_filter // self.params.significance_filter.max_ct

      # Ensure there is at least one bin.
      N_bins = max(
        [min([self.params.significance_filter.n_bins,N_bins_small_set]),
         N_bins_large_set, 1]
      )
      print("Total obs %d Choose n bins = %d"%(N_obs_pre_filter,N_bins), file=out)
      if indices_to_edge is not None:
        print("Total preds %d to edge of detector"%indices_to_edge.size(), file=out)
      bin_results = show_observations(observations, out=out, n_bins=N_bins)

      if result.get("fuller_kapton_absorption_correction", None) is not None:
        fkac = result["fuller_kapton_absorption_correction"][0]
        #assert that we have absorption corrections for all observations
        if len(fkac)!=N_obs_pre_filter: raise Sorry(
           "Fuller corrections %d don't match obs %d"%(len(fkac),N_obs_pre_filter))
        unobstructed = (fkac==1.0)

        xypred = result["mapped_predictions"][0]
        if len(xypred)!=N_obs_pre_filter: raise Sorry(
           "Mapped predictions %d don't match obs %d"%(len(xypred),N_obs_pre_filter))

        print("unobstructed",unobstructed.count(True), "obstructed", unobstructed.count(False), file=out)

        if False:
          from matplotlib import pyplot as plt
          xy_unobstructed = xypred.select(unobstructed)
          xy_obstructed = xypred.select(~unobstructed)
          plt.plot([a[0] for a in xy_unobstructed], [a[1] for a in xy_unobstructed], "r.")
          plt.plot([a[0] for a in xy_obstructed], [a[1] for a in xy_obstructed], "b.")
          #plt.plot([a[0] for a in xypred], [a[1] for a in xypred], "r.")
          plt.axes().set_aspect("equal")
          plt.show()

        from xfel.merging.absorption import show_observations as aso
        try:
          ASO = aso(observations, unobstructed, self.params, out=out, n_bins=N_bins)
        except Exception as e:
          # in development encountered:
          # RuntimeError, flex.mean() of empty array
          # ValueError, max() arg is empty sequence
          print("skipping image: could not process obstructed/unobstructed bins", file=out)
          return null_data(
            file_name=file_name, log_out=out.getvalue(), low_signal=True)
        observations = observations.select(ASO.master_selection)
        observations_original_index = observations_original_index.select(ASO.master_selection)

        if False:
          from matplotlib import pyplot as plt
          unobstructed_subset = (fkac.select(ASO.master_selection)==1.0)
          xy_unobstructed = xypred.select(ASO.master_selection).select(unobstructed_subset)
          xy_obstructed = xypred.select(ASO.master_selection).select(~unobstructed_subset)
          plt.plot([a[0] for a in xy_unobstructed], [a[1] for a in xy_unobstructed], "r.")
          plt.plot([a[0] for a in xy_obstructed], [a[1] for a in xy_obstructed], "b.")
          #plt.plot([a[0] for a in xypred], [a[1] for a in xypred], "r.")
          plt.axes().set_aspect("equal")
          plt.show()

      else:
        acceptable_resolution_bins = [
          bin.mean_I_sigI > self.params.significance_filter.sigma for bin in bin_results]
        acceptable_nested_bin_sequences = [i for i in range(len(acceptable_resolution_bins))
                                           if False not in acceptable_resolution_bins[:i+1]]
        if len(acceptable_nested_bin_sequences)==0:
          return null_data(
            file_name=file_name, log_out=out.getvalue(), low_signal=True)
        else:
          N_acceptable_bins = max(acceptable_nested_bin_sequences) + 1
          imposed_res_filter = float(bin_results[N_acceptable_bins-1].d_range.split()[2])
          imposed_res_sel = observations.resolution_filter_selection(
            d_min=imposed_res_filter)
          observations = observations.select(
            imposed_res_sel)
          observations_original_index = observations_original_index.select(
            imposed_res_sel)
          print("New resolution filter at %7.2f"%imposed_res_filter,file_name, file=out)
        print("N acceptable bins",N_acceptable_bins, file=out)
      print("Old n_obs: %d, new n_obs: %d"%(N_obs_pre_filter,observations.size()), file=out)
      if indices_to_edge is not None:
        print("Total preds %d to edge of detector"%indices_to_edge.size(), file=out)
      # Finished applying the binwise I/sigma filter---------------------------------------
    if self.params.raw_data.sdfac_auto is True:
      I_over_sig = observations.data()/observations.sigmas()
      #assert that at least a few I/sigmas are less than zero
      Nlt0 = I_over_sig.select(I_over_sig<0.).size()
      if Nlt0 > 2:
        # get a rough estimate for the SDFAC, assuming that negative measurements
        # represent false predictions and therefore normally distributed noise.
        no_signal = I_over_sig.select(I_over_sig<0.)
        for xns in range(len(no_signal)):
          no_signal.append(-no_signal[xns])
        Stats = flex.mean_and_variance(no_signal)
        SDFAC = Stats.unweighted_sample_standard_deviation()
      else: SDFAC=1.
      print("The applied SDFAC is %7.4f"%SDFAC, file=out)
      corrected_sigmas = observations.sigmas() * SDFAC
      observations = observations.customized_copy(sigmas = corrected_sigmas)
      observations_original_index = observations_original_index.customized_copy(
        sigmas = observations_original_index.sigmas() * SDFAC)

    print("Step 6.  Match to reference intensities, filter by correlation, filter out negative intensities.", file=out)
    assert len(observations_original_index.indices()) \
      ==   len(observations.indices())

    data = frame_data(self.n_refl, file_name)
    data.set_indexed_cell(indexed_cell)
    data.current_orientation = result['current_orientation'][0]
    data.d_min = observations.d_min()

    # Ensure that match_multi_indices() will return identical results
    # when a frame's observations are matched against the
    # pre-generated Miller set, self.miller_set, and the reference
    # data set, self.i_model.  The implication is that the same match
    # can be used to map Miller indices to array indices for intensity
    # accumulation, and for determination of the correlation
    # coefficient in the presence of a scaling reference.
    if self.i_model is not None:
      assert len(self.i_model.indices()) == len(self.miller_set.indices()) \
        and  (self.i_model.indices() ==
              self.miller_set.indices()).count(False) == 0

    matches = miller.match_multi_indices(
      miller_indices_unique=self.miller_set.indices(),
      miller_indices=observations.indices())

    if predictions is not None:
      matches_predictions = miller.match_multi_indices(
        miller_indices_unique=self.miller_set.indices(),
        miller_indices=predictions.indices())
    else:
      matches_predictions = None

    # matches_preds_obs = miller.match_multi_indices(
    #   miller_indices_unique=indices_to_edge,
    #   miller_indices=observations_original_index.indices())

    use_weights = False # New facility for getting variance-weighted correlation
    if self.params.scaling.algorithm in ['mark1','levmar']:
      # Because no correlation is computed, the correlation
      # coefficient is fixed at zero.  Setting slope = 1 means
      # intensities are added without applying a scale factor.
      sum_x = 0
      sum_y = 0
      for pair in matches.pairs():
        data.n_obs += 1
        if not self.params.include_negatives and observations.data()[pair[1]] <= 0:
          data.n_rejected += 1
        else:
          sum_y += observations.data()[pair[1]]
      N = data.n_obs - data.n_rejected

      slope = 1
      offset = 0
      corr = 0

    else:
      sum_xx = 0
      sum_xy = 0
      sum_yy = 0
      sum_x = 0
      sum_y = 0
      sum_w = 0.

      for pair in matches.pairs():
        if self.params.scaling.simulation is not None:
          observations.data()[pair[1]] = self.i_model.data()[pair[0]]     # SIM
          observations.sigmas()[pair[1]] = self.i_model.sigmas()[pair[0]] # SIM

        data.n_obs += 1
        if not self.params.include_negatives and observations.data()[pair[1]] <= 0:
          data.n_rejected += 1
          continue
        # Update statistics using reference intensities (I_r), and
        # observed intensities (I_o).
        if use_weights: I_w = 1./(observations.sigmas()[pair[1]])**2 #variance weighting
        else: I_w = 1.

        # skip the structure factor if i_model.sigma is invalid (intentionally < 0)
        if self.i_model.sigmas()[pair[0]] < 0:
          I_w = 0.

        I_r = self.i_model.data()[pair[0]]
        I_o = observations.data()[pair[1]]
        sum_xx += I_w * I_r**2
        sum_yy += I_w * I_o**2
        sum_xy += I_w * I_r * I_o
        sum_x += I_w * I_r
        sum_y += I_w * I_o
        sum_w += I_w
      # Linearly fit I_r to I_o, i.e. find slope and offset such that
      # I_o = slope * I_r + offset, optimal in a least-squares sense.
      # XXX This is backwards, really.
      N = data.n_obs - data.n_rejected
      DELTA = sum_w * sum_xx - sum_x**2 # see p. 105 in Bevington & Robinson
      if (DELTA) == 0:
        print("Skipping frame with",sum_w,sum_xx,sum_x**2, file=out)
        return null_data(file_name=file_name,
                         log_out=out.getvalue(),
                         low_signal=True)
      slope = (sum_w * sum_xy - sum_x * sum_y) / DELTA
      offset = (sum_xx * sum_y - sum_x * sum_xy) / DELTA
      corr = (sum_w * sum_xy - sum_x * sum_y) / (math.sqrt(sum_w * sum_xx - sum_x**2) *
                                                 math.sqrt(sum_w * sum_yy - sum_y**2))

    # Early return if there are no positive reflections on the frame.
    if data.n_obs <= data.n_rejected:
      return null_data(
        file_name=file_name, log_out=out.getvalue(), low_signal=True)

    # Update the count for each matched reflection.  This counts
    # reflections with non-positive intensities, too.
    data.completeness += matches.number_of_matches(0).as_int()
    if matches_predictions is not None:
      data.completeness_predictions += matches_predictions.number_of_matches(0).as_int()
    # print "updated observations count by %d and preds count by %d" % \
    #   (flex.sum(matches.number_of_matches(0).as_int()), flex.sum(matches_predictions.number_of_matches(0).as_int()))
    # print "preds matching indices:", len(matches_predictions.pairs())
    # print "preds total:", len(indices_to_edge)
    data.corr = corr
    data.wavelength = wavelength

    # Apply the correlation coefficient threshold, if appropriate.
    if self.params.scaling.algorithm == 'mark0' and \
       corr <= self.params.min_corr:
      print("Skipping these data - correlation too low.", file=out)
      data.set_log_out(out.getvalue())
      data.show_log_out(sys.stdout)
      return null_data(file_name=file_name, log_out=out.getvalue(), low_correlation=True)
    # Apply a resolution filter, if appropriate.
    if self.params.lattice_rejection.d_min and \
      observations.d_min() >= self.params.lattice_rejection.d_min:
      print("Skipping these data - diffraction worse than %.2f Angstrom" % self.params.lattice_rejection.d_min, file=out)
      data.set_log_out(out.getvalue())
      data.show_log_out(sys.stdout)
      return null_data(file_name=file_name, log_out=out.getvalue(), low_resolution=True)

    from xfel.cxi.postrefinement_factory import factory
    PF = factory(self.params)
    postrefinement_algorithm = PF.postrefinement_algorithm()

    if self.params.postrefinement.enable:
      # Refactorization of the Sauter(2015) code; result should be same to 5 significant figures.
      # Lack of binary identity is due to the use of Python for old-code weighted correlation,
      #   contrasted with flex.double arithmetic for new-code.
      postx=postrefinement_algorithm(observations_original_index, self.params,
           self.i_model, self.miller_set, result, out)
      try:
        postx.run_plain()
        observations_original_index,observations,matches = postx.result_for_cxi_merge(file_name)
      except (AssertionError,ValueError,RuntimeError) as e:
        return null_data(file_name=file_name, log_out=out.getvalue(), reason=e)

      self.postrefinement_params = postx.parameterization_class(postx.MINI.x)

      if self.params.postrefinement.show_trumpet_plot is True:
        from xfel.cxi.trumpet_plot import trumpet_wrapper
        trumpet_wrapper(result, postx, file_name, self.params, out)

    if not self.params.scaling.enable or self.params.postrefinement.enable: # Do not scale anything
      print("Scale factor to an isomorphous reference PDB will NOT be applied.", file=out)
      slope = 1.0
      offset = 0.0

    if db_mgr is None: return unpack(MINI.x) # special exit for two-color indexing

    frame_id_0_base = PF.insert_frame_call(locals())

    observations_original_index_indices = observations_original_index.indices()
    xypred = result["mapped_predictions"][0]
    indices = flex.size_t([pair[1] for pair in matches.pairs()])

    sel_observations = flex.intersection(
      size=observations.data().size(),
      iselections=[indices])

    if self.params.include_negatives_fix_27May2018:
      # Super-rare exception. If saved sigmas instead of I/sigmas in the ISIGI dict, this wouldn't be needed.
      sel_observations &= observations.data() != 0

    set_original_hkl = observations_original_index_indices.select(
      flex.intersection(
        size=observations_original_index_indices.size(),
        iselections=[indices]))
    set_xypred = xypred.select(
      flex.intersection(
        size=xypred.size(),
        iselections=[indices]))

    kwargs = {'hkl_id_0_base': [pair[0] for pair in matches.pairs()],
              'i': observations.data().select(sel_observations),
              'sigi': observations.sigmas().select(sel_observations),
              'detector_x': [xy[0] for xy in set_xypred],
              'detector_y': [xy[1] for xy in set_xypred],
              'frame_id_0_base': [frame_id_0_base] * len(matches.pairs()),
              'overload_flag': [0] * len(matches.pairs()),
              'original_h': [hkl[0] for hkl in set_original_hkl],
              'original_k': [hkl[1] for hkl in set_original_hkl],
              'original_l': [hkl[2] for hkl in set_original_hkl]}

    db_mgr.insert_observation(**kwargs)

    print("For %d reflections, got slope %f, correlation %f" % \
        (data.n_obs - data.n_rejected, slope, corr), file=out)
    print("average obs", sum_y / (data.n_obs - data.n_rejected), \
      "average calc", sum_x / (data.n_obs - data.n_rejected), file=out)
    print("Rejected %d reflections with negative intensities" % \
        data.n_rejected, file=out)

    data.extra_stuff = {}


    data.accept = True
    if self.params.postrefinement.enable and self.params.postrefinement.algorithm in ["rs_hybrid"]:
      assert slope == 1.0
      assert self.params.include_negatives
    for pair in matches.pairs():
      if not self.params.include_negatives and (observations.data()[pair[1]] <= 0) :
        continue
      Intensity = observations.data()[pair[1]] / slope
      # Super-rare exception. If saved sigmas instead of I/sigmas in the ISIGI dict, this wouldn't be needed.
      if Intensity == 0:
        continue

      # Add the reflection as a two-tuple of intensity and I/sig(I)
      # to the dictionary of observations.
      index = self.miller_set.indices()[pair[0]]
      isigi = (Intensity,
               observations.data()[pair[1]] / observations.sigmas()[pair[1]],
               slope)
      if index in data.ISIGI:
        data.ISIGI[index].append(isigi)
      else:
        data.ISIGI[index] = [isigi]
        data.extra_stuff[index] = flex.double(), flex.miller_index()

      data.extra_stuff[index][0].append(observations_original_index.data()[pair[1]])
      data.extra_stuff[index][1].append(observations_original_index.indices()[pair[1]])

      sigma = observations.sigmas()[pair[1]] / slope
      variance = sigma * sigma
      data.summed_N[pair[0]] += 1
      data.summed_wt_I[pair[0]] += Intensity / variance
      data.summed_weight[pair[0]] += 1 / variance
    print("Selected file %s to %5.2f Angstrom resolution limit" % (file_name, observations.d_min()), file=out)
    data.set_log_out(out.getvalue())
    data.show_log_out(sys.stdout)
    return data

  def sum_intensities(self):
    sum_I = flex.double(self.miller_set.size(), 0.)
    sum_I_SIGI = flex.double(self.miller_set.size(), 0.)
    for i in range(self.miller_set.size()) :
      index = self.miller_set.indices()[i]
      if index in self.ISIGI :
        for t in self.ISIGI[index]:
          sum_I[i] += t[0]
          sum_I_SIGI[i] += t[1]
    return sum_I, sum_I_SIGI

def consistent_set_and_model(work_params,i_model=None):
  # Adjust the minimum d-spacing of the generated Miller set to assure
  # that the desired high-resolution limit is included even if the
  # observed unit cell differs slightly from the target.  Use the same
  # expansion formula as used in merging/general_fcalc.py, to assure consistency.
  # If a reference model is present, ensure that Miller indices are ordered
  # identically.
  miller_set = symmetry(
      unit_cell=work_params.target_unit_cell,
      space_group_info=work_params.target_space_group
    ).build_miller_set(
      anomalous_flag=not work_params.merge_anomalous,
      d_max=work_params.d_max,
      d_min=work_params.d_min / math.pow(
        1 + work_params.unit_cell_length_tolerance, 1 / 3))
  miller_set = miller_set.change_basis(
    work_params.model_reindex_op).map_to_asu()

  if i_model is not None:
    # Handle the case where model is anomalous=False but the requested merging is anomalous=True
    if i_model.anomalous_flag() is False and miller_set.anomalous_flag() is True:
      i_model = i_model.generate_bijvoet_mates()
    # manage the sizes of arrays.  General_fcalc assures that
    # N(i_model) >= N(miller_set) since it fills non-matches with invalid structure factors
    # However, if N(i_model) > N(miller_set) it's because this run of cxi.merge requested
    # a smaller resolution range.  Must prune off the reference model.

    if i_model.indices().size() > miller_set.indices().size():
      matches = miller.match_indices(i_model.indices(), miller_set.indices())
      pairs = matches.pairs()
      i_model = i_model.select(pairs.column(0))

    matches = miller.match_indices(i_model.indices(), miller_set.indices())
    assert not matches.have_singles()
    miller_set = miller_set.select(matches.permutation())

  return miller_set, i_model
#-----------------------------------------------------------------------
def run(args):
  if ("--help" in args) :
    iotbx.phil.parse(master_phil).show(attributes_level=2)
    return
  processor = iotbx.phil.process_command_line(args=args, master_string=master_phil)
  if len(processor.remaining_args) > 0:
    print("The following arguments were not recognized:")
    print("\n".join(processor.remaining_args))
    iotbx.phil.parse(master_phil).show(attributes_level=2)
    return
  phil = processor.show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application
  application(work_params)

  if ((work_params.d_min is None) or
      (work_params.data is None) or
      ( (work_params.model is None) and work_params.scaling.algorithm != "mark1") ) :
    command_name = os.environ["LIBTBX_DISPATCHER_NAME"]
    raise Usage(command_name + " "
                "d_min=4.0 "
                "data=~/scratch/r0220/006/strong/ "
                "model=3bz1_3bz2_core.pdb")
  if ((work_params.rescale_with_average_cell) and
      (not work_params.set_average_unit_cell)) :
    raise Usage("If rescale_with_average_cell=True, you must also specify "+
      "set_average_unit_cell=True.")
  if [work_params.raw_data.sdfac_auto, work_params.raw_data.sdfac_refine, work_params.raw_data.errors_from_sample_residuals].count(True) > 1:
    raise Usage("Specify only one of sdfac_auto, sdfac_refine or errors_from_sample_residuals.")

  # Read Nat's reference model from an MTZ file.  XXX The observation
  # type is given as F, not I--should they be squared?  Check with Nat!
  log = open("%s.log" % work_params.output.prefix, "w")
  out = multi_out()
  out.register("log", log, atexit_send_to=None)
  out.register("stdout", sys.stdout)
  print("I model", file=out)
  if work_params.model is not None:
    from xfel.merging.general_fcalc import run
    i_model = run(work_params)
    work_params.target_unit_cell = i_model.unit_cell()
    work_params.target_space_group = i_model.space_group_info()
    i_model.show_summary()
  else:
    i_model = None

  print("Target unit cell and space group:", file=out)
  print("  ", work_params.target_unit_cell, file=out)
  print("  ", work_params.target_space_group, file=out)

  miller_set, i_model = consistent_set_and_model(work_params,i_model)

  frame_files = get_observations(work_params)
  scaler = scaling_manager(
    miller_set=miller_set,
    i_model=i_model,
    params=work_params,
    log=out)
  scaler.scale_all(frame_files)
  if scaler.n_accepted == 0:
    return None
  scaler.show_unit_cell_histograms()
  if (work_params.rescale_with_average_cell) :
    average_cell_abc = scaler.uc_values.get_average_cell_dimensions()
    average_cell = uctbx.unit_cell(list(average_cell_abc) +
      list(work_params.target_unit_cell.parameters()[3:]))
    work_params.target_unit_cell = average_cell
    print("", file=out)
    print("#" * 80, file=out)
    print("RESCALING WITH NEW TARGET CELL", file=out)
    print("  average cell: %g %g %g %g %g %g" % \
      work_params.target_unit_cell.parameters(), file=out)
    print("", file=out)
    scaler.reset()
    scaler.scale_all(frame_files)
    scaler.show_unit_cell_histograms()
  if False : #(work_params.output.show_plots) :
    try :
      plot_overall_completeness(completeness)
    except Exception as e :
      print("ERROR: can't show plots")
      print("  %s" % str(e))
  print("\n", file=out)

  # Sum the observations of I and I/sig(I) for each reflection.
  sum_I, sum_I_SIGI = scaler.sum_intensities()

  miller_set_avg = miller_set.customized_copy(
    unit_cell=work_params.target_unit_cell)
  table1 = show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.completeness,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for all reflections",
    out=out,
    work_params=work_params)
  print("", file=out)
  if work_params.model is not None:
    n_refl, corr = scaler.get_overall_correlation(sum_I)
  else:
    n_refl, corr = ((scaler.completeness > 0).count(True), 0)
  print("\n", file=out)
  table2 = show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.summed_N,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for reflections where I > 0",
    out=out,
    work_params=work_params)
  #from libtbx import easy_pickle
  #easy_pickle.dump(file_name="stats.pickle", obj=stats)
  #stats.report(plot=work_params.plot)
  #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
  #  stats.counts != 0)
  #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
  #  file_name="nobs.mtz")
  if work_params.data_subsubsets.subsubset is not None and work_params.data_subsubsets.subsubset_total is not None:
    easy_pickle.dump("scaler_%d.pickle"%work_params.data_subsubsets.subsubset, scaler)
  explanation = """
Explanation:
Completeness       = # unique Miller indices present in data / # Miller indices theoretical in asymmetric unit
Asu. Multiplicity  = # measurements / # Miller indices theoretical in asymmetric unit
Obs. Multiplicity  = # measurements / # unique Miller indices present in data
Pred. Multiplicity = # predictions on all accepted images / # Miller indices theoretical in asymmetric unit"""
  print(explanation, file=out)
  mtz_file, miller_array = scaler.finalize_and_save_data()
  #table_pickle_file = "%s_graphs.pkl" % work_params.output.prefix
  #easy_pickle.dump(table_pickle_file, [table1, table2])
  loggraph_file = os.path.abspath("%s_graphs.log" % work_params.output.prefix)
  f = open(loggraph_file, "w")
  f.write(table1.format_loggraph())
  f.write("\n")
  f.write(table2.format_loggraph())
  f.close()
  result = scaling_result(
    miller_array=miller_array,
    plots=scaler.get_plot_statistics(),
    mtz_file=mtz_file,
    loggraph_file=loggraph_file,
    obs_table=table1,
    all_obs_table=table2,
    n_reflections=n_refl,
    overall_correlation=corr)
  easy_pickle.dump("%s.pkl" % work_params.output.prefix, result)
  return result

def show_overall_observations(
  obs, redundancy, summed_wt_I, summed_weight, ISIGI,
  n_bins=15, out=None, title=None, work_params=None, redundancy_to_edge=None):
  if out is None:
    out = sys.stdout
  obs.setup_binner(d_max=100000, d_min=work_params.d_min, n_bins=n_bins)
  result = []

  cumulative_unique = 0
  cumulative_meas   = 0
  cumulative_n_pred = 0
  cumulative_pred   = 0
  cumulative_theor  = 0
  cumulative_In     = 0
  cumulative_I      = 0.0
  cumulative_Isigma = 0.0

  for i_bin in obs.binner().range_used():
    sel_w = obs.binner().selection(i_bin)
    sel_fo_all = obs.select(sel_w)
    obs.binner()._have_format_strings=True
    obs.binner().fmt_bin_range_used ="%6.3f -%6.3f"
    d_range = obs.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)
    sel_redundancy = redundancy.select(sel_w)
    if redundancy_to_edge is not None:
      sel_redundancy_pred = redundancy_to_edge.select(sel_w)
    else:
      sel_redundancy_pred = None
    sel_absent = sel_redundancy.count(0)
    n_present = sel_redundancy.size() - sel_absent
    sel_complete_tag = "[%d/%d]" % (n_present, sel_redundancy.size())
    sel_measurements = flex.sum(sel_redundancy)

    # Alternatively, redundancy (or multiplicity) is calculated as the
    # average number of observations for the observed
    # reflections--missing reflections do not affect the redundancy
    # adversely, and the reported value becomes
    # completeness-independent.
    val_redundancy_obs = 0
    if n_present > 0:
      val_redundancy_obs = flex.sum(sel_redundancy) / n_present
    # Repeat on the full set of predictions, calculated w.r.t. the asymmetric unit
    val_redundancy_pred = 0
    if redundancy_to_edge is not None and n_present > 0:
      val_redundancy_pred = flex.sum(sel_redundancy_pred) / sel_redundancy.size()

    # Per-bin sum of I and I/sig(I).  For any reflection, the weight
    # of the merged intensity must be positive for this to make sense.
    sel_o = (sel_w & (summed_weight > 0))
    intensity = summed_wt_I.select(sel_o) / summed_weight.select(sel_o)
    I_sum = flex.sum(intensity)
    I_sigI_sum = flex.sum(intensity * flex.sqrt(summed_weight.select(sel_o)))
    I_n = sel_o.count(True)

    if work_params is not None and \
     (work_params.plot_single_index_histograms and \
      N >= 30 and \
      work_params.data_subset == 0):
      scaling_manager.single_reflection_histograms(obs, ISIGI)

    if sel_measurements > 0:
      mean_I = mean_I_sigI = 0
      if I_n > 0:
        mean_I = I_sum / I_n
        mean_I_sigI = I_sigI_sum / I_n
      bin = resolution_bin(
        i_bin=i_bin,
        d_range=d_range,
        d_min=obs.binner().bin_d_min(i_bin),
        redundancy_asu=flex.mean(sel_redundancy.as_double()),
        redundancy_obs=val_redundancy_obs,
        redundancy_to_edge=val_redundancy_pred,
        complete_tag=sel_complete_tag,
        completeness=n_present / sel_redundancy.size(),
        measurements=sel_measurements,
        predictions=sel_redundancy_pred,
        mean_I=mean_I,
        mean_I_sigI=mean_I_sigI)
      result.append(bin)
    cumulative_unique += n_present
    cumulative_meas   += sel_measurements
    if redundancy_to_edge is not None:
      cumulative_n_pred += flex.sum(sel_redundancy_pred)
      cumulative_pred   += redundancy_to_edge
    cumulative_theor  += sel_redundancy.size()
    cumulative_In     += I_n
    cumulative_I      += I_sum
    cumulative_Isigma += I_sigI_sum

  if (title is not None) :
    print(title, file=out)
  from libtbx import table_utils
  table_header = ["","","","","<asu","<obs","<pred","","","",""]
  table_header2 = ["Bin","Resolution Range","Completeness","%","multi>","multi>","multi>",
    "n_meas", "n_pred","<I>","<I/sig(I)>"]
  use_preds = (redundancy_to_edge is not None and flex.sum(redundancy_to_edge) > 0)
  include_columns = [True, True, True, True, True, True, use_preds, True, use_preds, True, True]
  table_data = []
  table_data.append(table_header)
  table_data.append(table_header2)
  for bin in result:
    table_row = []
    table_row.append("%3d" % bin.i_bin)
    table_row.append("%-13s" % bin.d_range)
    table_row.append("%13s" % bin.complete_tag)
    table_row.append("%5.2f" % (100*bin.completeness))
    table_row.append("%6.2f" % bin.redundancy_asu)
    table_row.append("%6.2f" % bin.redundancy_obs)
    table_row.append("%6.2f" % (0 if redundancy_to_edge is None else bin.redundancy_to_edge))
    table_row.append("%6d" % bin.measurements)
    table_row.append("%6d" % (0 if redundancy_to_edge is None else flex.sum(bin.predictions)))
    table_row.append("%8.0f" % bin.mean_I)
    table_row.append("%8.3f" % bin.mean_I_sigI)
    table_data.append(table_row)
  if len(table_data) <= 2:
    print("Table could not be constructed -- no bins accepted.", file=out)
    return
  table_data.append([""]*len(table_header))
  table_data.append(  [
      format_value("%3s",   "All"),
      format_value("%-13s", "                 "),
      format_value("%13s",  "[%d/%d]"%(cumulative_unique,cumulative_theor)),
      format_value("%5.2f", 100*(cumulative_unique/cumulative_theor)),
      format_value("%6.2f", cumulative_meas/cumulative_theor),
      format_value("%6.2f", cumulative_meas/cumulative_unique),
      format_value("%6.2f", (0 if redundancy_to_edge is None else cumulative_n_pred/cumulative_theor)),
      format_value("%6d",   cumulative_meas),
      format_value("%6d",   (0 if redundancy_to_edge is None else flex.sum(redundancy_to_edge))),
      format_value("%8.0f", cumulative_I/cumulative_In),
      format_value("%8.3f", cumulative_Isigma/cumulative_In),
  ])
  table_data = table_utils.manage_columns(table_data, include_columns)

  print()
  print(table_utils.format(table_data,has_header=2,justify='center',delim=" "), file=out)

  # XXX generate table object for displaying plots
  if (title is None) :
    title = "Data statistics by resolution"
  table = data_plots.table_data(
    title=title,
    x_is_inverse_d_min=True,
    force_exact_x_labels=True)
  table.add_column(
    column=[1 / bin.d_min**2 for bin in result],
    column_name="d_min",
    column_label="Resolution")
  table.add_column(
    column=[bin.redundancy_asu for bin in result],
    column_name="redundancy",
    column_label="Redundancy")
  table.add_column(
    column=[bin.completeness for bin in result],
    column_name="completeness",
    column_label="Completeness")
  table.add_column(
    column=[bin.mean_I_sigI for bin in result],
    column_name="mean_i_over_sigI",
    column_label="<I/sig(I)>")
  table.add_graph(
    name="Redundancy vs. resolution",
    type="GRAPH",
    columns=[0,1])
  table.add_graph(
    name="Completeness vs. resolution",
    type="GRAPH",
    columns=[0,2])
  table.add_graph(
    name="<I/sig(I)> vs. resolution",
    type="GRAPH",
    columns=[0,3])
  return table

class resolution_bin(object):
  def __init__(self,
               i_bin=None,
               d_range=None,
               d_min=None,
               redundancy_asu=None,
               redundancy_obs=None,
               redundancy_to_edge=None,
               absent=None,
               complete_tag=None,
               completeness=None,
               measurements=None,
               predictions=None,
               mean_I=None,
               mean_I_sigI=None,
               sigmaa=None):
    adopt_init_args(self, locals())

class scaling_result (group_args) :
  """
  Container for any objects that might need to be saved for future use (e.g.
  in a GUI).  Must be pickle-able!
  """
  pass

#-----------------------------------------------------------------------
# graphical goodies
def plot_overall_completeness(completeness):
  completeness_range = range(-1,flex.max(completeness)+1)
  completeness_counts = [completeness.count(n) for n in completeness_range]
  from matplotlib import pyplot as plt
  plt.plot(completeness_range,completeness_counts,"r+")
  plt.show()

class plot_statistics (object) :
  """
  Container for assorted histograms of frame statistics.  The resolution bin
  plots are stored separately, since they can be displayed using the loggraph
  viewer.
  """
  def __init__ (self,
                prefix,
                unit_cell_statistics,
                reference_cell,
                correlations,
                min_corr,
                rejected_fractions,
                frame_d_min) :
    adopt_init_args(self, locals())

  def show_all_pyplot (self, n_slots=20) :
    """
    Display histograms using pyplot.  For use in a wxPython GUI the figure
    should be created separately in a wx.Frame.
    """
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(9,12))
    self.plot_unit_cell_histograms(
      figure=fig,
      a_values=self.unit_cell_statistics.all_uc_a_values,
      b_values=self.unit_cell_statistics.all_uc_b_values,
      c_values=self.unit_cell_statistics.all_uc_c_values,
      n_slots=n_slots,
      title=\
        "Unit cell length distribution (all frames with compatible indexing): %s" % self.prefix)
    plt.show()
    fig = plt.figure(figsize=(9,12))
    self.plot_unit_cell_histograms(
      figure=fig,
      a_values=self.unit_cell_statistics.uc_a_values,
      b_values=self.unit_cell_statistics.uc_b_values,
      c_values=self.unit_cell_statistics.uc_c_values,
      n_slots=n_slots,
      title=\
        "Unit cell length distribution (frames with acceptable correlation): %s" % self.prefix)
    plt.show()
    fig = plt.figure(figsize=(9,12))
    self.plot_statistics_histograms(
      figure=fig,
      n_slots=n_slots)
    plt.show()

  def plot_unit_cell_histograms (self,
      figure,
      a_values,
      b_values,
      c_values,
      n_slots=20,
      title="Distribution of unit cell edge lengths") :
    [a0,b0,c0,alpha0,beta0,gamma0] = self.reference_cell.parameters()
    ax1 = figure.add_axes([0.1, 0.1, 0.8, 0.25])
    ax2 = figure.add_axes([0.1, 0.4, 0.8, 0.25])
    ax3 = figure.add_axes([0.1, 0.7, 0.8, 0.25])
    ax1.hist(c_values, n_slots, color=[1.0,0.0,0.0])
    ax2.hist(b_values, n_slots, color=[0.0,1.0,0.0])
    ax3.hist(a_values, n_slots, color=[0.0,0.5,1.0])
    ax1.axvline(c0, linestyle='.', linewidth=2, color='k')
    ax2.axvline(b0, linestyle='.', linewidth=2, color='k')
    ax3.axvline(a0, linestyle='.', linewidth=2, color='k')
    ax1.set_xlabel("c edge")
    ax2.set_xlabel("b edge")
    ax3.set_xlabel("a edge")
    ax3.set_title("%s: %s" % (title, self.prefix))

  def plot_statistics_histograms (self,
      figure,
      n_slots=20) :
    ax1 = figure.add_axes([0.1, 0.1, 0.8, 0.25])
    ax2 = figure.add_axes([0.1, 0.4, 0.8, 0.25])
    ax3 = figure.add_axes([0.1, 0.7, 0.8, 0.25])
    ax1.hist(self.correlations, n_slots, color=[1.0,0.0,0.0])
    ax2.hist(self.rejected_fractions, n_slots, color=[0.0,1.0,0.0])
    ax3.hist(self.frame_d_min, n_slots, color=[0.0,0.5,1.0])
    ax1.axvline(self.min_corr, linestyle='.', linewidth=2, color='k')
    ax1.set_xlabel("Correlation to reference dataset")
    ax2.set_xlabel("Fraction of rejected zero or negative intensities")
    ax3.set_xlabel("Integrated resolution limit")
    ax1.set_title("Correlation by frame (%s)" % self.prefix)
    ax2.set_title("Rejected reflections by frame (%s)" % self.prefix)
    ax3.set_title("Resolution by frame (%s)" % self.prefix)

if (__name__ == "__main__"):
  show_plots = False
  if ("--plots" in sys.argv) :
    sys.argv.remove("--plots")
    show_plots = True
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)
  if (show_plots) :
    try :
      result.plots.show_all_pyplot()
      from wxtbx.command_line import loggraph
      loggraph.run([result.loggraph_file])
    except Exception as e :
      print("Can't display plots")
      print("You should be able to view them by running this command:")
      print("  wxtbx.loggraph %s" % result.loggraph_file)
      raise e
