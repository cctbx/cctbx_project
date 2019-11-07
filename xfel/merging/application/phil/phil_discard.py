from __future__ import absolute_import, division, print_function

from xfel.merging.database.merging_database import mysql_master_phil
master_phil="""
short_circuit = False
  .type = bool
  .help = Assess per-image resolution limits and exit early.

data_subset = 0
  .type = int
  .help = 0: use all data / 1: use odd-numbered frames / 2: use even-numbered frames

wavelength = None
  .type = float
pixel_size = None
  .type = float
  .help = Detector-specific parameter for pixel size in mm
include_bulk_solvent = True
  .type = bool
elements = None
  .type = str
  .multiple = True


include_negatives = False
  .type = bool
  .help = Whether to include negative intensities during scaling and merging
include_negatives_fix_27May2018 = True
  .type = bool
  .help = Bugfix for include negatives. Affects cxi.xmerge.

rescale_with_average_cell = False
  .type = bool
  NOTE:  deprecating this option in the new program.  If the unit cell is not known
         to precision, run the cluster algorithm first (without scaling and merging)
         then run the program again with the cluster-averaged unit cell.

merging {
  refine_G_Imodel = False
    .type = bool
    .help = "Refine per-frame scaling factors and model intensities
             after initial merging"
}
scaling {
  algorithm = levmar
    .type = choice
    .help = "Sparse-matrix Levenberg-Marquardt scaling and merging.  Under development, not available in cxi.merge"
  simulation = None
    .type = str
    .help = To test scaling, use simulated data from the model instead of the actual observations
    .help = String value governs how the sim data are calculated
  simulation_data = None
    .type = floats
    .help = Extra parameters for the simulation, exact meaning depends on calculation method
}

error {
  propagate_errors = False
    .type = bool
    .help = Propagate errors from estimated parameters
}

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

data_subsubsets {
  subsubset = None
    .type = int
  subsubset_total = None
    .type = int
}



""" + mysql_master_phil


if (__name__ == "__main__"):
  pass
