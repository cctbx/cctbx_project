from __future__ import division, print_function

import iotbx.phil

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


if (__name__ == "__main__"):
  pass
