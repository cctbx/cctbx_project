# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.merge
#
# $Id$

from rstbx.apps.stills.simple_integration import show_observations
import iotbx.phil
from iotbx import data_plots
from cctbx.array_family import flex
from cctbx import miller
from cctbx.crystal import symmetry
from cctbx.sgtbx.bravais_types import bravais_lattice
from cctbx import uctbx
from libtbx.str_utils import format_value
from libtbx.utils import Usage, multi_out
from libtbx import easy_pickle
from libtbx import adopt_init_args, group_args, Auto
from cStringIO import StringIO
import os
import math
import time
import sys
op = os.path

master_phil="""
data = None
  .type = path
  .multiple = True
  .help = Directory containing integrated data in pickle format.  Repeat to \
    specify additional directories.
model = None
  .type = str
  .help = PDB filename containing atomic coordinates & isomorphous cryst1 record
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
k_sol = 0.35
  .type = float
  .help = bulk solvent scale factor - approximate mean value in PDB \
    (according to Pavel)
b_sol = 46.00
  .type = float
  .help = bulk solvent B-factor - approximate mean value in PDB \
    (according to Pavel)
include_bulk_solvent = True
  .type = bool
wavelength = None
  .type = float
min_corr = 0.1
  .type = float
  .help = Correlation cutoff for rejecting individual frames
unit_cell_length_tolerance = 0.1
  .type = float
  .help = Fractional change in unit cell dimensions allowed (versus target \
    cell).
unit_cell_angle_tolerance = 2.
  .type = float
nproc = None
  .type = int
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
"""

def get_observations (data_dirs):
  file_names = []
  for dir_name in data_dirs :
    for file_name in os.listdir(dir_name):
      if (file_name.endswith(".pickle")):
        file_names.append(os.path.join(dir_name, file_name))
  print "Number of pickle files found:", len(file_names)
  print
  return file_names

class WrongBravaisError (Exception) :
  pass

class OutlierCellError (Exception) :
  pass

def load_result (file_name,
                 ref_bravais_type,
                 reference_cell,
                 params,
                 out) :
  """
  Take a pickle file, confirm that it contains the appropriate data, and
  check the lattice type and unit cell against the reference settings - if
  rejected, raises an exception (for tracking statistics).
  """
  if (file_name.endswith("stats.pickle")) :
    return None
  obj = easy_pickle.load(file_name=file_name)
  if (not obj.has_key("observations")) :
    return None
  result_array = obj["observations"][0]
  unit_cell = result_array.unit_cell()
  sg_info = result_array.space_group_info()
  print >> out, ""
  print >> out, "-" * 80
  print >> out, file_name
  print >> out, sg_info
  print >> out, unit_cell

  HARD_CODED_PIXEL_SZ_MM = 0.11
  assert obj['mapped_predictions'][0].size() == obj["observations"][0].size()
  mm_predictions = HARD_CODED_PIXEL_SZ_MM*(obj['mapped_predictions'][0])
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
  print >> out, "Integrated data:"
  result_array.show_summary(f=out, prefix="  ")
  # XXX don't force reference setting here, it will be done later, after the
  # original unit cell is recorded
  return obj

class intensity_data (object) :
  """
  Container for scaled intensity data.
  """
  def __init__ (self, n_refl) :
    self.n_refl = n_refl
    self.initialize()

  def initialize (self) :
    self.ISIGI        = {}
    self.completeness = flex.int(self.n_refl, 0)
    self.summed_N     = flex.int(self.n_refl, 0)
    self.summed_weight= flex.double(self.n_refl, 0.)
    self.summed_wt_I  = flex.double(self.n_refl, 0.)

class frame_data (intensity_data) :
  """
  Intensity data for a single frame.
  """
  def __init__ (self, n_refl, file_name) :
    intensity_data.__init__(self, n_refl)
    self.file_name = file_name
    self.n_obs = 0
    self.n_rejected = 0
    self.corr = 0
    self.d_min = -1
    self.accept = False
    self.indexed_cell = None
    self.log_out = file_name
    self.wavelength = None

  def set_indexed_cell (self, unit_cell) :
    self.indexed_cell = unit_cell

  def set_log_out (self, out_str) :
    self.log_out = out_str

  def show_log_out (self, out) :
    print >> out, self.log_out

class null_data (object) :
  """
  Stand-in for a frame rejected due to conflicting symmetry.  (No flex arrays
  included, to save pickling time during multiprocessing.)
  """
  def __init__ (self, file_name, log_out, wrong_cell, wrong_bravais) :
    adopt_init_args(self, locals())

  def show_log_out (self, out) :
    print >> out, self.log_out

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
    print >> out, ""
    labels = ["a","b","c"]
    ref_edges = [a0,b0,c0]
    def _show_each (edges) :
      for edge, ref_edge, label in zip(edges, ref_edges, labels) :
        h = flex.histogram(edge, n_slots=n_slots)
        smin, smax = flex.min(edge), flex.max(edge)
        stats = flex.mean_and_variance(edge)
        print >> out, "  %s edge" % label
        print >> out, "     range:     %6.2f - %.2f" % (smin, smax)
        print >> out, "     mean:      %6.2f +/- %6.2 on N = %d" % (
          flex.mean(edge), stats.unweighted_sample_standard_deviation(), edge.size())
        print >> out, "     reference: %6.2f" % ref_edge
        h.show(f=out, prefix="    ", format_cutoffs="%6.2f")
        print >> out, ""
    edges = [self.all_uc_a_values, self.all_uc_b_values, self.all_uc_c_values]
    print >> out, \
      "Unit cell length distribution (all frames with compatible indexing):"
    _show_each(edges)
    edges = [self.uc_a_values, self.uc_b_values, self.uc_c_values]
    print >> out, \
      "Unit cell length distribution (frames with acceptable correlation):"
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
    intensity_data.__init__(self, i_model.size())
    self.reset()

  def reset (self) :
    self.n_processed = 0
    self.n_accepted = 0
    self.n_wrong_bravais = 0
    self.n_wrong_cell = 0
    self.n_low_corr = 0
    self.observations = flex.int()
    self.corr_values = flex.double()
    self.rejected_fractions = flex.double()
    self.uc_values = unit_cell_distribution()
    self.d_min_values = flex.double()
    self.wavelength = flex.double()
    self.initialize()

  def scale_all (self, file_names) :
    t1 = time.time()
    try :
      import multiprocessing
    except ImportError, e :
      print >> self.log, \
        "multiprocessing module not available (requires Python >= 2.6)"
      print >> self.log, "will scale frames serially"
      self._scale_all_serial(file_names)
    else :
      if (self.params.nproc == 1) :
        self._scale_all_serial(file_names)
      else :
        self._scale_all_parallel(file_names)
    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(file_names))
    print >> self.log, "  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais
    print >> self.log, "  %d rejected for unit cell outliers" % \
      self.n_wrong_cell
    print >> self.log, "  %d rejected due to poor correlation" % \
      self.n_low_corr

  def _scale_all_parallel (self, file_names) :
    import multiprocessing
    import libtbx.introspection
    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto) :
      nproc = libtbx.introspection.number_of_processors()
    pool = multiprocessing.Pool(processes=nproc)
    # Round-robin the frames through the process pool.  Each process
    # accumulates its own statistics in serial, and the grand total is
    # eventually collected by the main process' _add_all_frames()
    # function.
    for i in xrange(nproc) :
      sm = scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[[file_names[j] for j in xrange(i, len(file_names), nproc)]],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

  def _scale_all_serial (self, file_names) :
    """
    Scale frames sequentially (single-process)
    """
    for file_name in file_names :
      scaled = self.scale_frame(file_name)
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
      return
    #data.show_log_out(self.log)
    #self.log.flush()
    if (isinstance(data, null_data)) :
      if (data.wrong_bravais) :
        self.n_wrong_bravais += 1
      elif (data.wrong_cell) :
        self.n_wrong_cell += 1
      return
    if (data.accept) :
      self.n_accepted    += 1
      self.completeness  += data.completeness
      self.summed_N      += data.summed_N
      self.summed_weight += data.summed_weight
      self.summed_wt_I   += data.summed_wt_I
      for index, isigi in data.ISIGI.iteritems() :
        if (index in self.ISIGI):
          self.ISIGI[index] += isigi
        else:
          self.ISIGI[index] = isigi
    else :
      self.n_low_corr += 1
    self.uc_values.add_cell(data.indexed_cell,
      rejected=(not data.accept))
    self.observations.append(data.n_obs)
    if (data.n_obs > 0) :
      frac_rejected = data.n_rejected / data.n_obs
      self.rejected_fractions.append(frac_rejected)
      self.d_min_values.append(data.d_min)
    self.corr_values.append(data.corr)
    self.wavelength.append(data.wavelength)
    del data

  def _add_all_frames (self, data) :
    """The _add_all_frames() function collects the statistics
    accumulated in @p data by the individual scaling processes in
    process pool.  XXX Sure this does not need a lock?
    """
    self.n_accepted += data.n_accepted
    self.n_low_corr += data.n_low_corr
    self.n_processed += data.n_processed
    self.n_wrong_bravais += data.n_wrong_bravais
    self.n_wrong_cell += data.n_wrong_cell

    for index, isigi in data.ISIGI.iteritems() :
      if (index in self.ISIGI):
        self.ISIGI[index] += isigi
      else:
        self.ISIGI[index] = isigi

    self.completeness += data.completeness
    self.summed_N += data.summed_N
    self.summed_weight += data.summed_weight
    self.summed_wt_I += data.summed_wt_I

    self.corr_values.extend(data.corr_values)
    self.d_min_values.extend(data.d_min_values)
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
    for i in xrange(len(self.summed_N)):
      if (self.summed_N[i] <= 0):
        continue
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
    print >> self.log, \
      "SUMMARY: For %d reflections, got slope %f, correlation %f" \
        % (N, slope, corr)
    return N, corr

  def finalize_and_save_data (self) :
    """
    Assemble a Miller array with the summed data, setting the unit cell to
    the consensus average if desired, and write to an MTZ file (including
    merged/non-anomalous data too).
    """
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "OUTPUT FILES"
    Iobs_all = flex.double(self.i_model.size())
    SigI_all = flex.double(self.i_model.size())
    for i in xrange(len(Iobs_all)):
      if (self.summed_weight[i] > 0.):
        Iobs_all[i] = self.summed_wt_I[i] / self.summed_weight[i]
        SigI_all[i] = math.sqrt(1. / self.summed_weight[i])
    if (self.params.set_average_unit_cell) :
      # XXX since XFEL crystallography runs at room temperature, it may not
      # be appropriate to use the cell dimensions from a cryo structure.
      # also, some runs seem to have huge variance in the indexed cell
      # dimensions, so downstream programs (MR, refinement) may run better
      # with the cell set to the mean edge lengths.
      abc = self.uc_values.get_average_cell_dimensions()
      print >> self.log, "  (will use final unit cell edges %g %g %g)" % abc
      angles = self.miller_set.unit_cell().parameters()[3:]
      unit_cell = uctbx.unit_cell(list(abc) + list(angles))
      final_symm = symmetry(
        unit_cell=unit_cell,
        space_group_info=self.miller_set.space_group_info())
    else :
      final_symm = self.miller_set
    all_obs = self.i_model.customized_copy(data=Iobs_all,
      sigmas=SigI_all,
      crystal_symmetry=final_symm).set_observation_type_xray_intensity()
    mtz_file = "%s.mtz" % self.params.output.prefix
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
    print >> self.log, "  Anomalous and mean data:\n    %s" % \
      os.path.abspath(mtz_file)
    print >> self.log, ""
    print >> self.log, "Final data:"
    all_obs.show_summary(self.log, prefix="  ")
    return mtz_file, all_obs

  def __call__ (self, file_names) :
    try :
      return self._scale_all_serial(file_names)
    except Exception, e :
      print >> self.log, str(e)
      return None

  def scale_frame (self, file_name) :
    """
    Scales the data from a single frame against the reference dataset, and
    returns an intensity_data object.  Can be called either serially or
    via a multiprocessing map() function.
    """
    # XXX VERY IMPORTANT: this method must not modify any internal data or
    # the parallelization will not yield usable results!
    out = StringIO()
    wrong_cell = wrong_bravais = False
    try :
      result = load_result(
        file_name=file_name,
        reference_cell=self.params.target_unit_cell,
        ref_bravais_type=self.ref_bravais_type,
        params=self.params,
        out=out)
    except OutlierCellError, e :
      print >> out, str(e)
      result = None
      wrong_cell = True
    except WrongBravaisError, e :
      print >> out, str(e)
      result = None
      wrong_bravais = True
    if (result is None) :
      null = null_data(
        file_name=file_name,
        log_out=out.getvalue(),
        wrong_bravais=wrong_bravais,
        wrong_cell=wrong_cell)
      return null

    # The wavelength from the command line overrides the wavelength in
    # the pickled integration file.  XXX The wavelength parameter
    # should probably be removed from master_phil once all integration
    # pickle files contain it.
    if (self.params.wavelength is not None):
      wavelength = self.params.wavelength
    elif (result.has_key("wavelength")):
      wavelength = result["wavelength"]
    else:
      # XXX Give error, or raise exception?
      return None
    assert (wavelength > 0)

    observations = result["observations"][0]
    cos_two_polar_angle = result["cos_two_polar_angle"]

    assert observations.size() == cos_two_polar_angle.size()
    tt_vec = observations.two_theta(wavelength)
    cos_tt_vec = flex.cos( tt_vec.data() )
    sin_tt_vec = flex.sin( tt_vec.data() )
    cos_sq_tt_vec = cos_tt_vec * cos_tt_vec
    sin_sq_tt_vec = sin_tt_vec * sin_tt_vec
    P_nought_vec = 0.5 * (1. + cos_sq_tt_vec)

    F_prime = -1.0 # Hard-coded value defines the incident polarization axis
    P_prime = 0.5 * F_prime * cos_two_polar_angle * sin_sq_tt_vec
    observations = observations / ( P_nought_vec - P_prime )
    # This corrects observations for polarization assuming 100% polarization on
    # one axis (thus the F_prime = -1.0 rather than the perpendicular axis, 1.0)
    # Polarization model as described by Kahn, Fourme, Gadet, Janin, Dumas & Andre
    # (1982) J. Appl. Cryst. 15, 330-337, equations 13 - 15.

    indexed_cell = observations.unit_cell()
    # Now do manipulate the data to conform to unit cell, asu, and space group
    # of reference
    # Only works if there is NOT an indexing ambiguity!
    observations = observations.customized_copy(
      crystal_symmetry=self.miller_set.crystal_symmetry()
      ).resolution_filter(d_min=self.params.d_min).map_to_asu()
    print >> out, "Data in reference setting:"
    #observations.show_summary(f=out, prefix="  ")
    show_observations(observations, out=out)
    # Match up the observed intensities against the reference data
    # set, i_model, instead of the pre-generated miller set,
    # miller_set.
    matches = miller.match_indices(
      self.i_model.indices(),
      observations.indices())
    data = frame_data(self.n_refl, file_name)
    data.set_indexed_cell(indexed_cell)
    data.d_min = observations.d_min()
    # Update the count for each matched reflection.  This counts
    # reflections with negative intensities, too.
    data.completeness +=  (~matches.single_selection(0)).as_int()

    data.wavelength = wavelength

    # Initialise first- and second-order statistics.
    N = 0
    sum_xx = 0
    sum_xy = 0
    sum_yy = 0
    sum_x = 0
    sum_y = 0
    for pair in matches.pairs():
      data.n_obs += 1
      if (observations.data()[pair[1]] <= 0):
        data.n_rejected += 1
        continue
      # Update statistics using reference intensities (I_r), and
      # observed intensities (I_o).
      I_r = self.i_model.data()[pair[0]]
      I_o = observations.data()[pair[1]]
      N      += 1
      sum_xx += I_r**2
      sum_yy += I_o**2
      sum_xy += I_r * I_o
      sum_x  += I_r
      sum_y  += I_o
    # Linearly fit I_r to I_o, i.e. find slope and offset such that
    # I_o = slope * I_r + offset, optimal in a least-squares sense.
    # XXX This is backwards, really.
    slope = (N * sum_xy - sum_x * sum_y) / (N * sum_xx - sum_x**2)
    offset = (sum_xx * sum_y - sum_x * sum_xy) / (N * sum_xx - sum_x**2)
    corr  = (N * sum_xy - sum_x * sum_y) / (math.sqrt(N * sum_xx - sum_x**2) *
             math.sqrt(N * sum_yy - sum_y**2))
    data.corr = corr
    print >> out, "For %d reflections, got slope %f, correlation %f" % \
        (N, slope, corr)
    print >> out, "average obs",sum_y/N, "average calc",sum_x/N
    print >> out, "Rejected %d reflections with negative intensities" % \
        (len(matches.pairs()) - N)
    if (corr > self.params.min_corr) :
      data.accept = True
      for pair in matches.pairs():
        if (observations.data()[pair[1]] <= 0) :
          continue
        Intensity = observations.data()[pair[1]] /  slope

        # Add the reflection as a two-tuple of intensity and I/sig(I)
        # to the dictionary of observations.
        index = self.i_model.indices()[pair[0]]
        isigi = (Intensity,
                 observations.data()[pair[1]] / observations.sigmas()[pair[1]])
        if (index in data.ISIGI):
          data.ISIGI[index].append(isigi)
        else:
          data.ISIGI[index] = [isigi]

        sigma = observations.sigmas()[pair[1]] /  slope
        variance = sigma * sigma
        data.summed_N[pair[0]] += 1
        data.summed_wt_I[pair[0]] += Intensity / variance
        data.summed_weight[pair[0]] += 1. / variance
    else :
      print >> out, "Skipping these data - correlation too low."
    data.set_log_out(out.getvalue())
    if corr > 0.5:
      print "Selected file %s"%file_name.replace("integration","out").replace("int","idx")
      print "Selected distance %6.2f mm"%float(result["distance"])
      data.show_log_out(sys.stdout)
    return data

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
    return

  if ((work_params.d_min is None) or
      (work_params.data is None) or
      (work_params.model is None)) :
    raise Usage("cxi.merge "
                "d_min=4.0 "
                "data=~/scratch/r0220/006/strong/ "
                "model=3bz1_3bz2_core.pdb")
  if ((work_params.rescale_with_average_cell) and
      (not work_params.set_average_unit_cell)) :
    raise Usage("If rescale_with_average_cell=True, you must also specify "+
      "set_average_unit_cell=True.")
  # Read Nat's reference model from an MTZ file.  XXX The observation
  # type is given as F, not I--should they be squared?  Check with Nat!
  log = open("%s.log" % work_params.output.prefix, "w")
  out = multi_out()
  out.register("log", log, atexit_send_to=None)
  out.register("stdout", sys.stdout)
  print >> out, "I model"
  from xfel.cxi.merging.general_fcalc import run
  i_model = run(work_params)
  i_model.show_summary()
  if (work_params.target_unit_cell is None) :
    work_params.target_unit_cell = i_model.unit_cell()
  if (work_params.target_space_group is None) :
    work_params.target_space_group = i_model.space_group_info()

  print >> out, "Target unit cell and space group:"
  print >> out, "  ", work_params.target_unit_cell
  print >> out, "  ", work_params.target_space_group

  miller_set = symmetry(
      unit_cell=work_params.target_unit_cell,
      space_group_info=work_params.target_space_group
    ).build_miller_set(
      anomalous_flag=True,
      d_min=work_params.d_min)

  #miller_set.show_summary()
  #from iotbx import mtz
  #i_model = mtz.object(file_name = "data.mtz").as_miller_arrays()[0]
  #i_model.show_summary()
  #exit()
  # reality check
  #recip_cell_volume = work_params.target_unit_cell.reciprocal().volume()
  #recip_sphere_volume = (4/3)*math.pi*math.pow(1./work_params.d_min,3)
  #resolution_cells = recip_sphere_volume/recip_cell_volume
  #print "Number of asu's in sphere=",resolution_cells/miller_set.size()

  frame_files = get_observations(work_params.data)
  scaler = scaling_manager(
    miller_set=miller_set,
    i_model=i_model,
    params=work_params,
    log=out)
  scaler.scale_all(frame_files)
  scaler.show_unit_cell_histograms()
  if (work_params.rescale_with_average_cell) :
    average_cell_abc = scaler.uc_values.get_average_cell_dimensions()
    average_cell = uctbx.unit_cell(list(average_cell_abc) +
      list(work_params.target_unit_cell.parameters()[3:]))
    work_params.target_unit_cell = average_cell
    print >> out, ""
    print >> out, "#" * 80
    print >> out, "RESCALING WITH NEW TARGET CELL"
    print >> out, "  average cell: %g %g %g %g %g %g" % \
      work_params.target_unit_cell.parameters()
    print >> out, ""
    scaler.reset()
    scaler.scale_all(frame_files)
    scaler.show_unit_cell_histograms()
  if False : #(work_params.output.show_plots) :
    try :
      plot_overall_completeness(completeness)
    except Exception, e :
      print "ERROR: can't show plots"
      print "  %s" % str(e)
  print >> out, "\n"

  # Sum the observations of I and I/sig(I) for each reflection.
  sum_I = flex.double(i_model.size(), 0.)
  sum_I_SIGI = flex.double(i_model.size(), 0.)
  for i in xrange(i_model.size()) :
    index = i_model.indices()[i]
    if index in scaler.ISIGI :
      for t in scaler.ISIGI[index]:
        sum_I[i] += t[0]
        sum_I_SIGI[i] += t[1]

  table1 = show_overall_observations(
    obs=i_model,
    redundancy=scaler.completeness,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for all reflections",
    out=out)
  print >> out, ""
  n_refl, corr = scaler.get_overall_correlation(sum_I)
  print >> out, "\n"
  table2 = show_overall_observations(
    obs=i_model,
    redundancy=scaler.summed_N,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for reflections where I > 0",
    out=out)
  #from libtbx import easy_pickle
  #easy_pickle.dump(file_name="stats.pickle", obj=stats)
  #stats.report(plot=work_params.plot)
  #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
  #  stats.counts != 0)
  #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
  #  file_name="nobs.mtz")
  print >> out, ""
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
  obs, redundancy, ISIGI, n_bins=15, out=None, title=None):
  if out==None:
    out = sys.stdout
  obs.setup_binner(n_bins=n_bins)
  result = []

  # R_iso_tot and R_merge_tot are two-membered lists, holding the
  # numerator and the denominator.
  R_iso_tot = [0, 0]
  R_merge_tot = [0, 0]

  cumulative_unique = 0
  cumulative_meas   = 0
  cumulative_theor  = 0
  cumulative_Isigma = 0.0

  for i_bin in obs.binner().range_used():
    sel_w = obs.binner().selection(i_bin)
    sel_fo_all = obs.select(sel_w)
    d_max_,d_min_ = sel_fo_all.d_max_min()
    d_range = obs.binner().bin_legend(
      i_bin=i_bin, show_bin_number=False, show_counts=False)
    sel_redundancy = redundancy.select(sel_w)
    sel_absent = sel_redundancy.count(0)
    n_present = sel_redundancy.size() - sel_absent
    sel_complete_tag = "[%d/%d]" % (n_present, sel_redundancy.size())
    sel_measurements = flex.sum(sel_redundancy)

    # Per-bin sum of I and I/sig(I) for each observation.  Accumulate
    # numerators for R_merge (Stout & Jensen, 1968) and R_iso (Chapman
    # et al., 2011).
    I_sum = 0
    I_sigI_sum = 0
    R_iso = [0, 0]
    R_merge = [0, 0]
    for i in obs.binner().array_indices(i_bin) :
      index = obs.indices()[i]
      if (index in ISIGI) :
        # Compute m, the "merged" intensity, as the average intensity
        # of all observations of the reflection with the given index.
        N = 0
        m = 0
        for t in ISIGI[index] :
          I_sigI_sum += t[1]
          N += 1
          m += t[0]
        I_sum += m
        if (N > 0):
          m /= N
          R_iso[0] += abs(m - obs.data()[i])
          R_iso[1] += m
          for t in ISIGI[index] :
            R_merge[0] += abs(t[0] - m)
            R_merge[1] += t[0]

    # Keep track of total sums for global statistics.
    for i in xrange(2) :
      R_iso_tot[i] += R_iso[i]
      R_merge_tot[i] += R_merge[i]

    if (sel_measurements > 0 and R_iso[1] > 0 and R_merge[1] > 0):
      bin = resolution_bin(
        i_bin        = i_bin,
        d_range      = d_range,
        d_min        = obs.binner().bin_d_min(i_bin),
        redundancy   = flex.mean(sel_redundancy.as_double()),
        complete_tag = sel_complete_tag,
        completeness = n_present / sel_redundancy.size(),
        measurements = sel_measurements,
        mean_I       = I_sum / sel_measurements,
        mean_I_sigI  = I_sigI_sum / sel_measurements,
        R_iso        = R_iso[0] / R_iso[1],
        R_merge      = R_merge[0] / R_merge[1],
        )
      result.append(bin)
    cumulative_unique += n_present
    cumulative_meas   += sel_measurements
    cumulative_theor  += sel_redundancy.size()
    cumulative_Isigma += I_sigI_sum

  if (title is not None) :
    print >> out, title
  print >>out, "\n Bin  Resolution Range  Completeness <Redundancy>  n_meas      <I> <I/sig(I)>    R_iso  R_merge"
  for bin in result:
    fmt = " %s %s %s       %s  %s %s   %s %s %s"
    print >>out,fmt%(
      format_value("%3d",   bin.i_bin),
      format_value("%-13s", bin.d_range), # XXX This doesn't always work, compare thermolysin runs with PSII with spots out to 181 A.
      format_value("%13s",  bin.complete_tag),
      format_value("%6.2f", bin.redundancy),
      format_value("%6d",   bin.measurements),
      format_value("%8.0f", bin.mean_I),
      format_value("%8.3f", bin.mean_I_sigI),
      format_value("%8.3f", bin.R_iso),
      format_value("%8.3f", bin.R_merge)
    )
  print >>out,fmt%(
      format_value("%3s",   "All"),
      format_value("%-13s", "                 "),
      format_value("%13s",  "[%d/%d]"%(cumulative_unique,cumulative_theor)),
      format_value("%6.2f", cumulative_meas/cumulative_theor),
      format_value("%6d",   cumulative_meas),
      format_value("%8.0f", 0.),
      format_value("%8.3f", cumulative_Isigma/cumulative_meas),
      format_value("%8.3f", 0.),
      format_value("%8.3f", 0.)
    )
  print "Global R_iso   ", R_iso_tot[0] / R_iso_tot[1]
  print "Global R_merge ", R_merge_tot[0] / R_merge_tot[1]
  # XXX generate table object for displaying plots
  if (title is None) :
    title = "Data statistics by resolution"
  table = data_plots.table_data(
    title=title,
    x_is_inverse_d_min=True,
    force_exact_x_labels=True)
  table.add_column(
    column=[ 1/bin.d_min**2 for bin in result ],
    column_name="d_min",
    column_label="Resolution")
  table.add_column(
    column=[ bin.redundancy for bin in result ],
    column_name="redundancy",
    column_label="Redundancy")
  table.add_column(
    column=[ bin.completeness for bin in result ],
    column_name="completeness",
    column_label="Completeness")
  table.add_column(
    column=[ bin.mean_I_sigI for bin in result ],
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
               i_bin         = None,
               d_range       = None,
               d_min         = None,
               redundancy    = None,
               absent        = None,
               complete_tag  = None,
               completeness  = None,
               measurements  = None,
               mean_I        = None,
               mean_I_sigI   = None,
               R_iso         = None,
               R_merge       = None,
               sigmaa        = None):
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
  completeness_range = xrange(-1,flex.max(completeness)+1)
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
    ax3.hist(self.d_min_values, n_slots, color=[0.0,0.5,1.0])
    ax1.axvline(self.min_corr, linestyle='.', linewidth=2, color='k')
    ax1.set_xlabel("Correlation to reference dataset")
    ax2.set_xlabel("Fraction of rejected zero or negative intensities")
    ax3.set_xlabel("Integrated resolution limit")
    ax1.set_title("Correlation by frame (%s)" % self.prefix)
    ax2.set_title("Rejected reflections by frame (%s)" % self.prefix)
    ax3.set_title("Resolution by frame (%s)" % self.prefix)
    plt.show()

if (__name__ == "__main__"):
  show_plots = False
  if ("--plots" in sys.argv) :
    sys.argv.remove("--plots")
    show_plots = True
  result = run(args=sys.argv[1:])
  if (show_plots) :
    try :
      result.plots.show_all_pyplot()
      from wxtbx.command_line import loggraph
      loggraph.run([result.loggraph_file])
    except Exception, e :
      print "Can't display plots"
      print "You should be able to view them by running this command:"
      print "  wxtbx.loggraph %s" % result.loggraph_file
      raise e
