from __future__ import absolute_import, division, print_function
from six.moves import range

from rstbx.dials_core.integration_core import show_observations
import iotbx.phil
from cctbx.array_family import flex
from cctbx import miller
from cctbx import uctbx
from libtbx.utils import Usage, multi_out
from libtbx import easy_pickle
from libtbx import adopt_init_args, group_args, Auto
import os
import math
import time
import sys
from scitbx import matrix
import six
from six.moves import zip
op = os.path

from xfel.command_line.cxi_merge import master_phil
from xfel.command_line.cxi_merge import get_observations
from xfel.command_line.cxi_merge import frame_data, null_data
from xfel.command_line.cxi_merge import consistent_set_and_model

class unit_cell_distribution(object):
  """
  Container for collecting unit cell edge length statistics - for frames
  included in the final dataset.
  (Frames with incompatible indexing solutions will not be included.)
  """
  # TODO make this more general - currently assumes that angles are fixed,
  # which is true for the systems studied so far
  def __init__(self):
    self.all_uc_a_values = flex.double()
    self.all_uc_b_values = flex.double()
    self.all_uc_c_values = flex.double()

  def add_cell(self, unit_cell, rejected=False):
    if (unit_cell is None):
      return
    (a,b,c,alpha,beta,gamma) = unit_cell.parameters()
    self.all_uc_a_values.append(a)
    self.all_uc_b_values.append(b)
    self.all_uc_c_values.append(c)

  def add_cells(self, uc):
    """Addition operation for unit cell statistics."""
    self.all_uc_a_values.extend(uc.all_uc_a_values)
    self.all_uc_b_values.extend(uc.all_uc_b_values)
    self.all_uc_c_values.extend(uc.all_uc_c_values)

  def show_histograms(self, reference, out, n_slots=20):
    [a0,b0,c0,alpha0,beta0,gamma0] = reference.parameters()
    print("", file=out)
    labels = ["a","b","c"]
    ref_edges = [a0,b0,c0]
    def _show_each(edges):
      for edge, ref_edge, label in zip(edges, ref_edges, labels):
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

  def get_average_cell_dimensions(self):
    a = flex.mean(self.all_uc_a_values)
    b = flex.mean(self.all_uc_b_values)
    c = flex.mean(self.all_uc_c_values)
    return a,b,c

#-----------------------------------------------------------------------
from xfel.command_line.cxi_merge import scaling_manager as scaling_manager_base
class scaling_manager(scaling_manager_base):
  def __init__(self, miller_set, i_model, params, log=None):
    scaling_manager_base.__init__(self, miller_set, i_model, params, log)

  def reset(self):
    self.n_processed = 0
    self.n_accepted = 0
    self.n_file_error = 0
    self.n_low_signal = 0
    self.n_wrong_bravais = 0
    self.n_wrong_cell = 0
    self.observations = flex.int()
    self.corr_values = flex.double()
    self.rejected_fractions = flex.double()
    self.uc_values = unit_cell_distribution()
    self.d_min_values = flex.double()
    self.wavelength = flex.double()
    self.initialize()

  def scale_all(self, file_names):
    t1 = time.time()
    if self.params.backend == 'MySQL':
      from xfel.merging.database.merging_database import manager
    elif self.params.backend == 'SQLite':
      from xfel.merging.database.merging_database_sqlite3 import manager
    elif self.params.backend == 'FS':
      from xfel.merging.database.merging_database_fs import manager
    elif self.params.backend == 'Flex':
      from xfel.merging.database.merging_database_flex import manager

    import multiprocessing
    print("Allocating intensities")
    intensity_proxy = multiprocessing.Array('d',self.params.memory.shared_array_allocation,lock=True)
    print("Allocating sigmas")
    sigma_proxy = multiprocessing.Array('d',self.params.memory.shared_array_allocation,lock=True)
    print("Allocating frame_id")
    frame_proxy = multiprocessing.Array('l',self.params.memory.shared_array_allocation,lock=True)
    print("Allocating miller_id")
    miller_proxy = multiprocessing.Array('l',self.params.memory.shared_array_allocation,lock=True)
    H_proxy = multiprocessing.Array('i',self.params.memory.shared_array_allocation,lock=True)
    K_proxy = multiprocessing.Array('i',self.params.memory.shared_array_allocation,lock=True)
    L_proxy = multiprocessing.Array('i',self.params.memory.shared_array_allocation,lock=True)
    xtal_proxy = multiprocessing.Array('c',self.params.memory.shared_array_allocation,lock=True)
    print("Finished allocating")
    rows_observation = multiprocessing.Array('l',[0],lock=True)

    data_dict = dict(intensity_proxy=intensity_proxy,
                     sigma_proxy=sigma_proxy,
                     frame_proxy=frame_proxy,
                     miller_proxy=miller_proxy,
                     H_proxy=H_proxy,
                     K_proxy=K_proxy,
                     L_proxy=L_proxy,
                     rows=rows_observation,
                     xtal_proxy=xtal_proxy
                     )
    db_mgr = manager(self.params,data_dict)
    db_mgr.initialize_db(self.miller_set)
    # Unless the number of requested processes is greater than one,
    # try parallel multiprocessing on a parallel host.  Block until
    # all database commands have been processed.
    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto):
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
    pickled_data = db_mgr.join(data_dict)

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
    print("  %d rejected for low signal" % \
      self.n_low_signal, file=self.log)
    print("  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error, file=self.log)
    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell
    assert checksum == len(file_names)

    self.join_obs = pickled_data

  def _scale_all_parallel(self, file_names, db_mgr):
    import multiprocessing
    import libtbx.introspection

    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto):
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
    for i in range(nproc):
      sm = scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[input_queue, db_mgr],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

    # Block until the input queue has been emptied.
    input_queue.join()


  def _scale_all_serial(self, file_names, db_mgr):
    """
    Scale frames sequentially (single-process).  The return value is
    picked up by the callback.
    """
    for file_name in file_names :
      scaled = self.scale_frame(file_name, db_mgr)
      if (scaled is not None):
        self.add_frame(scaled)
    return (self)

  def add_frame(self, data):
    """
    Combine the scaled data from a frame with the current overall dataset.
    Also accepts None or null_data objects, when data are unusable but we
    want to record the file as processed.
    """
    self.n_processed += 1
    if (data is None):
      return None
    #data.show_log_out(self.log)
    #self.log.flush()
    if (isinstance(data, null_data)):
      if (data.file_error):
        self.n_file_error += 1
      elif (data.low_signal):
        self.n_low_signal += 1
      elif (data.wrong_bravais):
        self.n_wrong_bravais += 1
      elif (data.wrong_cell):
        self.n_wrong_cell += 1
      return
    if (data.accept):
      self.n_accepted    += 1
      self.completeness  += data.completeness
      self.summed_N      += data.summed_N
      self.summed_weight += data.summed_weight
      self.summed_wt_I   += data.summed_wt_I
      for index, isigi in six.iteritems(data.ISIGI):
        if (index in self.ISIGI):
          self.ISIGI[index] += isigi
        else:
          self.ISIGI[index] = isigi

    self.uc_values.add_cell(data.indexed_cell,
      rejected=(not data.accept))
    self.observations.append(data.n_obs)
    if (data.n_obs > 0):
      frac_rejected = data.n_rejected / data.n_obs
      self.rejected_fractions.append(frac_rejected)
      self.d_min_values.append(data.d_min)
    self.corr_values.append(data.corr)
    self.wavelength.append(data.wavelength)

  def _add_all_frames(self, data):
    """The _add_all_frames() function collects the statistics accumulated
    in @p data by the individual scaling processes in the process
    pool.  This callback function is run in serial, so it does not
    need a lock.
    """
    self.n_accepted += data.n_accepted
    self.n_file_error += data.n_file_error
    self.n_low_signal += data.n_low_signal
    self.n_processed += data.n_processed
    self.n_wrong_bravais += data.n_wrong_bravais
    self.n_wrong_cell += data.n_wrong_cell

    for index, isigi in six.iteritems(data.ISIGI):
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

  def get_plot_statistics(self):
    return plot_statistics(
      prefix=self.params.output.prefix,
      unit_cell_statistics=self.uc_values,
      reference_cell=self.miller_set.unit_cell(),
      correlations=self.corr_values,
      rejected_fractions=self.rejected_fractions,
      frame_d_min=self.d_min_values)

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
    cos_two_polar_angle = result["cos_two_polar_angle"]

    assert observations.size() == cos_two_polar_angle.size()
    tt_vec = observations.two_theta(wavelength)
    #print "mean tt degrees",180.*flex.mean(tt_vec.data())/math.pi
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
    print("mean diff is",flex.mean(diff2), "range",flex.min(diff2), flex.max(diff2))
    # XXX done
    observations = observations / ( P_nought_vec - P_prime )
    # This corrects observations for polarization assuming 100% polarization on
    # one axis (thus the F_prime = -1.0 rather than the perpendicular axis, 1.0)
    # Polarization model as described by Kahn, Fourme, Gadet, Janin, Dumas & Andre
    # (1982) J. Appl. Cryst. 15, 330-337, equations 13 - 15.

    print("Step 3. Correct for polarization.")
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
    observations = observations.customized_copy(
      anomalous_flag=not self.params.merge_anomalous,
      crystal_symmetry=self.miller_set.crystal_symmetry()
      ).map_to_asu()

    observations_original_index = observations_original_index.customized_copy(
      anomalous_flag=not self.params.merge_anomalous,
      crystal_symmetry=self.miller_set.crystal_symmetry()
      )
    print("Step 4. Filter on global resolution and map to asu")
    print("Data in reference setting:", file=out)
    #observations.show_summary(f=out, prefix="  ")
    show_observations(observations, out=out)

    #if self.params.significance_filter.apply is True:
    #  raise Exception("significance filter not implemented in samosa")
    if self.params.significance_filter.apply is True: #------------------------------------
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
      print("Total obs %d Choose n bins = %d"%(N_obs_pre_filter,N_bins))
      bin_results = show_observations(observations, out=out, n_bins=N_bins)
      #show_observations(observations, out=sys.stdout, n_bins=N_bins)
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
        print("New resolution filter at %7.2f"%imposed_res_filter,file_name)
      print("N acceptable bins",N_acceptable_bins)
      print("Old n_obs: %d, new n_obs: %d"%(N_obs_pre_filter,observations.size()))
      print("Step 5. Frame by frame resolution filter")
      # Finished applying the binwise I/sigma filter---------------------------------------

    if self.params.raw_data.sdfac_auto is True:
      raise Exception("sdfac auto not implemented in samosa.")

    print("Step 6.  Match to reference intensities, filter by correlation, filter out negative intensities.")
    assert len(observations_original_index.indices()) \
      ==   len(observations.indices())

    data = frame_data(self.n_refl, file_name)
    data.set_indexed_cell(indexed_cell)
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

    # Early return if there are no positive reflections on the frame.
    if data.n_obs <= data.n_rejected:
      return null_data(
        file_name=file_name, log_out=out.getvalue(), low_signal=True)

    # Update the count for each matched reflection.  This counts
    # reflections with non-positive intensities, too.
    data.completeness += matches.number_of_matches(0).as_int()
    data.wavelength = wavelength

    if not self.params.scaling.enable: # Do not scale anything
      print("Scale factor to an isomorphous reference PDB will NOT be applied.")
      slope = 1.0
      offset = 0.0

    observations_original_index_indices = observations_original_index.indices()
    if db_mgr is None: return unpack(MINI.x) # special exit for two-color indexing

    kwargs = {'wavelength': wavelength,
              'beam_x': result['xbeam'],
              'beam_y': result['ybeam'],
              'distance': result['distance'],
              'unique_file_name': data.file_name}

    ORI = result["current_orientation"][0]
    Astar = matrix.sqr(ORI.reciprocal_matrix())

    kwargs['res_ori_1'] = Astar[0]
    kwargs['res_ori_2'] = Astar[1]
    kwargs['res_ori_3'] = Astar[2]
    kwargs['res_ori_4'] = Astar[3]
    kwargs['res_ori_5'] = Astar[4]
    kwargs['res_ori_6'] = Astar[5]
    kwargs['res_ori_7'] = Astar[6]
    kwargs['res_ori_8'] = Astar[7]
    kwargs['res_ori_9'] = Astar[8]
    assert self.params.scaling.report_ML is True
    kwargs['half_mosaicity_deg'] = result["ML_half_mosaicity_deg"][0]
    kwargs['domain_size_ang'] = result["ML_domain_size_ang"][0]

    frame_id_0_base = db_mgr.insert_frame(**kwargs)

    xypred = result["mapped_predictions"][0]
    indices = flex.size_t([pair[1] for pair in matches.pairs()])

    sel_observations = flex.intersection(
      size=observations.data().size(),
      iselections=[indices])
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

    print("Lattice: %d reflections" % (data.n_obs - data.n_rejected), file=out)
    print("average obs", sum_y / (data.n_obs - data.n_rejected), \
      "average calc", sum_x / (data.n_obs - data.n_rejected), file=out)
    print("Rejected %d reflections with negative intensities" % \
        data.n_rejected, file=out)

    data.accept = True
    for pair in matches.pairs():
      if not self.params.include_negatives and (observations.data()[pair[1]] <= 0):
        continue
      Intensity = observations.data()[pair[1]]
      # Super-rare exception. If saved sigmas instead of I/sigmas in the ISIGI dict, this wouldn't be needed.
      if Intensity == 0:
        continue

      # Add the reflection as a two-tuple of intensity and I/sig(I)
      # to the dictionary of observations.
      index = self.miller_set.indices()[pair[0]]
      isigi = (Intensity,
               observations.data()[pair[1]] / observations.sigmas()[pair[1]],
               1.0)
      if index in data.ISIGI:
        data.ISIGI[index].append(isigi)
      else:
        data.ISIGI[index] = [isigi]

      sigma = observations.sigmas()[pair[1]]
      variance = sigma * sigma
      data.summed_N[pair[0]] += 1
      data.summed_wt_I[pair[0]] += Intensity / variance
      data.summed_weight[pair[0]] += 1 / variance


    data.set_log_out(out.getvalue())
    return data

#-----------------------------------------------------------------------
def run(args):
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application,samosa
  application(work_params)
  samosa(work_params)
  if ("--help" in args):
    libtbx.phil.parse(master_phil.show())
    return

  if ((work_params.d_min is None) or
      (work_params.data is None) ):
    command_name = os.environ["LIBTBX_DISPATCHER_NAME"]
    raise Usage(command_name + " "
                "d_min=4.0 "
                "data=~/scratch/r0220/006/strong/ "
                "model=3bz1_3bz2_core.pdb")
  if ((work_params.rescale_with_average_cell) and
      (not work_params.set_average_unit_cell)):
    raise Usage("If rescale_with_average_cell=True, you must also specify "+
      "set_average_unit_cell=True.")
  if work_params.raw_data.sdfac_auto and work_params.raw_data.sdfac_refine:
    raise Usage("Cannot specify both sdfac_auto and sdfac_refine")

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
  if (work_params.rescale_with_average_cell):
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
  if False : #(work_params.output.show_plots):
    try :
      plot_overall_completeness(completeness)
    except Exception as e :
      print("ERROR: can't show plots")
      print("  %s" % str(e))
  print("\n", file=out)

  # Sum the observations of I and I/sig(I) for each reflection.
  sum_I = flex.double(miller_set.size(), 0.)
  sum_I_SIGI = flex.double(miller_set.size(), 0.)
  for i in range(miller_set.size()):
    index = miller_set.indices()[i]
    if index in scaler.ISIGI :
      for t in scaler.ISIGI[index]:
        sum_I[i] += t[0]
        sum_I_SIGI[i] += t[1]

  miller_set_avg = miller_set.customized_copy(
    unit_cell=work_params.target_unit_cell)
  table1 = show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.completeness,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for all reflections",
    out=out,
    work_params=work_params)
  print("", file=out)
  n_refl, corr = ((scaler.completeness > 0).count(True), 0)
  print("\n", file=out)
  table2 = show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.summed_N,
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
  print("", file=out)
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

from xfel.command_line.cxi_merge import show_overall_observations

class scaling_result(group_args):
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

from xfel.command_line.cxi_merge import plot_statistics as plot_base
class plot_statistics(plot_base):
  """
  Container for assorted histograms of frame statistics.  The resolution bin
  plots are stored separately, since they can be displayed using the loggraph
  viewer.
  """
  def __init__(self,
                prefix,
                unit_cell_statistics,
                reference_cell,
                correlations,
                rejected_fractions,
                frame_d_min):
    adopt_init_args(self, locals())

  def show_all_pyplot(self, n_slots=20):
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
    self.plot_statistics_histograms(
      figure=fig,
      n_slots=n_slots)
    plt.show()

  def plot_statistics_histograms(self,
      figure,
      n_slots=20):
    ax3 = figure.add_axes([0.1, 0.7, 0.8, 0.25])
    ax3.hist(self.frame_d_min, n_slots, color=[0.0,0.5,1.0])
    ax3.set_xlabel("Integrated resolution limit")
    ax3.set_title("Resolution by frame (%s)" % self.prefix)

if (__name__ == "__main__"):
  show_plots = False
  if ("--plots" in sys.argv):
    sys.argv.remove("--plots")
    show_plots = True
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)
  if (show_plots):
    result.plots.show_all_pyplot()
    from wxtbx.command_line import loggraph
    loggraph.run([result.loggraph_file])
