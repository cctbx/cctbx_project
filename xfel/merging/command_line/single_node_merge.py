# LIBTBX_SET_DISPATCHER_NAME singlenode.merge
from __future__ import division
from xfel.command_line.cxi_merge import run
from xfel.command_line.cxi_merge import OutlierCellError, WrongBravaisError
import sys,glob,os,math,time
from xfel.command_line import cxi_merge
from cctbx.array_family import flex
from cctbx.sgtbx.bravais_types import bravais_lattice
from libtbx import Auto
from cStringIO import StringIO
from libtbx.utils import Sorry
from scitbx import matrix
from xfel.cxi.merging_utils import intensity_data, frame_data, null_data

def get_observations (work_params):
  try:
    data_dirs = work_params.data
    data_subset = work_params.data_subset
    subsubset = work_params.data_subsubsets.subsubset
    subsubset_total = work_params.data_subsubsets.subsubset_total
    extension = work_params.filename_extension
  except Exception,e:
    exit("Changed the interface for get_observations, please contact authors "+str(e))

  print "Step 1.  Get a list of all files"
  assert work_params.a_list is None
  integration_pickle_names = []
  assert work_params.targlob is not None
  tar_list = glob.glob(work_params.targlob)
  for tarname in tar_list:
      import tarfile
      T = tarfile.open(name=tarname, mode='r')
      K = T.getmembers()
      NT = len(K)
      for nt in xrange(NT):
        k = os.path.basename(K[nt].path)
        integration_pickle_names.append("%s;member%05d;timestamp%s"%(tarname,nt,k))
      print tarname,NT
  print "Number of pickle files found:", len(integration_pickle_names)
  print
  return tar_list,integration_pickle_names
cxi_merge.get_observations = get_observations

def load_result (file_hash,
                 frame_obj,
                 ref_bravais_type,
                 reference_cell,
                 params,
                 reindex_op,
                 out,
                 get_predictions_to_edge=False,
                 image_info=None,
                 exclude_CSPAD_sensor=None) :

  print >> out, "-" * 80
  print >> out, "Step 2.  Load pickle file into dictionary obj and filter on lattice & cell with",reindex_op
  print >> out, file_hash
  """
  Take a pickle file, confirm that it contains the appropriate data, and
  check the lattice type and unit cell against the reference settings - if
  rejected, raises an exception (for tracking statistics).
  """

  obj = frame_obj
  if (not obj.has_key("observations")) :
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
    for idx in xrange(len(obj["observations"])):
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
  print >> out, ""

  print >> out, sg_info
  print >> out, unit_cell

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
  print >> out, "Integrated data:"
  result_array.show_summary(f=out, prefix="  ")
  # XXX don't force reference setting here, it will be done later, after the
  # original unit cell is recorded

  #Remove observations on a selected sensor if requested
  if exclude_CSPAD_sensor is not None:
    print >> out, "excluding CSPAD sensor %d" % exclude_CSPAD_sensor
    fast_min1, slow_min1, fast_max1, slow_max1, \
    fast_min2, slow_min2, fast_max2, slow_max2 = \
      get_boundaries_from_sensor_ID(exclude_CSPAD_sensor)
    fast_min = min(fast_min1, fast_min2)
    slow_min = min(slow_min1, slow_min2)
    fast_max = max(fast_max1, fast_max2)
    slow_max = max(slow_max1, slow_max2)
    accepted = flex.bool()
    px_preds = obj['mapped_predictions'][0]
    for idx in xrange(px_preds.size()):
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

from xfel.command_line.cxi_merge import scaling_manager as scaling_manager_base
class scaling_manager(scaling_manager_base):

  def tar_to_scale_frame_adapter(self, tar_list, db_mgr):
    for tarname in tar_list:
      import tarfile
      T = tarfile.open(name=tarname, mode='r')
      K = T.getmembers()
      NT = len(K)
      for nt in xrange(NT):
        k = os.path.basename(K[nt].path)
        file_hash = "%s;member%05d;timestamp%s"%(tarname,nt,k)
        this_member = K[nt]
        fileIO = T.extractfile(member=this_member)
        import cPickle as pickle
        frame_obj = pickle.load(fileIO)
        scaled = self.scale_frame(file_hash, frame_obj, db_mgr)
        if (scaled is not None) :
          self.add_frame(scaled)

  def scale_all (self, file_names) :
    tar_list,integration_pickle_names = file_names
    t1 = time.time()
    if self.params.backend == 'MySQL':
      from xfel.cxi.merging_database import manager
    elif self.params.backend == 'SQLite':
      from xfel.cxi.merging_database_sqlite3 import manager
    else:
      from xfel.cxi.merging_database_fs import manager

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
        self._scale_all_parallel(tar_list, db_mgr)
      except ImportError, e :
        print >> self.log, \
          "multiprocessing module not available (requires Python >= 2.6)\n" \
          "will scale frames serially"
        self._scale_all_serial(tar_list, db_mgr)
    else:
      self._scale_all_serial(tar_list, db_mgr)
    db_mgr.join()

    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(integration_pickle_names))
    print >> self.log, "  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais
    print >> self.log, "  %d rejected for unit cell outliers" % \
      self.n_wrong_cell
    print >> self.log, "  %d rejected for low signal" % \
      self.n_low_signal
    print >> self.log, "  %d rejected due to up-front poor correlation under min_corr parameter" % \
      self.n_low_corr
    print >> self.log, "  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error
    for key in self.failure_modes.keys():
      print >>self.log, "  %d rejected due to %s"%(self.failure_modes[key], key)

    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_corr + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell \
               + sum([val for val in self.failure_modes.itervalues()])
    assert checksum == len(integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print >> self.log, "Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min)

    if self.params.raw_data.sdfac_refine:
      self.scale_errors()

    if self.params.raw_data.errors_from_sample_residuals:
      self.errors_from_residuals()

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
      print file_name
      input_queue.put(file_name)
    pool = multiprocessing.Pool(processes=nproc)
    # Each process accumulates its own statistics in serial, and the
    # grand total is eventually collected by the main process'
    # _add_all_frames() function.
    for i in xrange(nproc) :
      sm = scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[input_queue, db_mgr],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

    # Block until the input queue has been emptied.
    input_queue.join()

  def _scale_all_serial (self, tar_list, db_mgr) :
    """
    Scale frames sequentially (single-process).  The return value is
    picked up by the callback.
    """
    self.tar_to_scale_frame_adapter(tar_list, db_mgr)
    return (self)

  def __call__ (self, input_queue, db_mgr) :
    """ Scale frames sequentially within the current process.  The
    return value is picked up by the callback.  See also
    self.scale_all_serial()"""
    from Queue import Empty

    try :
      while True:
        try:
          file_name = input_queue.get_nowait()
        except Empty:
          return self
        self.tar_to_scale_frame_adapter(tar_list=[file_name,], db_mgr=db_mgr)
        input_queue.task_done()

    except Exception, e :
      print >> self.log, str(e)
      return None

  def scale_frame (self, file_hash, frame_obj, db_mgr) :
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
      reindex_op = self.reverse_lookup.get(file_hash, None)
      if reindex_op is None:
        return null_data(file_name=file_hash, log_out=out.getvalue(), file_error=True)
    if (self.params.predictions_to_edge.apply) and (self.params.predictions_to_edge.image is not None):
      from xfel.merging.predictions_to_edges import ImageInfo
      image_info = ImageInfo(self.params.predictions_to_edge.image, detector_phil=self.params.predictions_to_edge.detector_phil)
    else:
      image_info = None
    try :
      result = load_result(
        file_hash,
        frame_obj,
        reference_cell=self.params.target_unit_cell,
        ref_bravais_type=self.ref_bravais_type,
        params=self.params,
        reindex_op = reindex_op,
        out=out,
        get_predictions_to_edge=self.params.predictions_to_edge.apply,
        image_info=image_info,
        exclude_CSPAD_sensor=self.params.validation.exclude_CSPAD_sensor)
      if result is None:
        return null_data(
          file_name=file_hash, log_out=out.getvalue(), file_error=True)
    except OutlierCellError, e :
      print >> out, str(e)
      return null_data(
        file_name=file_hash, log_out=out.getvalue(), wrong_cell=True)
    except WrongBravaisError, e :
      print >> out, str(e)
      return null_data(
        file_name=file_hash, log_out=out.getvalue(), wrong_bravais=True)
    return self.scale_frame_detail(result, file_hash, db_mgr, out)

cxi_merge.scaling_manager = scaling_manager

if (__name__ == "__main__"):
  result = run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

