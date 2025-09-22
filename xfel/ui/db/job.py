from __future__ import absolute_import, division, print_function
from xfel.ui import settings_dir
from xfel.ui.db import db_proxy, get_run_path, write_xtc_locator, get_image_mode
import os, shutil, copy

known_job_statuses = ["DONE", "ERR", "PEND", "RUN", "SUSP", "PSUSP", "SSUSP", "UNKWN", "EXIT", "DONE", "ZOMBI", "DELETED", "SUBMIT_FAIL", "SUBMITTED", "HOLD", "TIMEOUT"]
finished_job_statuses = ["DONE", "EXIT", "DELETED", "UNKWN", "ERR", "SUBMIT_FAIL", "TIMEOUT"]

class JobFactory(object):
  @staticmethod
  def from_job(job, task_type = None):
    if job.task_id is None:
      return IndexingJob(job.app, job.id, **job._db_dict)

    if task_type is None:
      task_type = job.app.get_task(job.task_id).type
    if task_type == "indexing":
      return IndexingJob(job.app, job.id, **job._db_dict)
    if task_type == "ensemble_refinement":
      return EnsembleRefinementJob(job.app, job.id, **job._db_dict)
    if task_type == "scaling":
      return ScalingJob(job.app, job.id, **job._db_dict)
    if task_type == "merging":
      return MergingJob(job.app, job.id, **job._db_dict)
    if task_type == "phenix":
      return PhenixJob(job.app, job.id, **job._db_dict)

  @staticmethod
  def from_args(app, job_id = None, **kwargs):
    if 'task_type' in kwargs:
      task_type = kwargs.pop('task_type')
    else:
      task_type = None
    return JobFactory.from_job(Job(app, job_id, **kwargs), task_type=task_type)

class Job(db_proxy):
  def __init__(self, app, job_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_job" % app.params.experiment_tag, id = job_id, **kwargs)
    self.job_id = self.id
    self._run = None
    self._rungroup = None
    self._trial = None
    self._task = None
    self._dataset = None
    self._dataset_version = None

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name in ["run", "rungroup", "trial", "task", "dataset", "dataset_version"]:
      _name = "_" + name
      name_id = name + "_id"
      if getattr(self, _name) is None:
        if name == "dataset_version":
          if self.dataset_id is not None:
            self._dataset_version = self.app.get_job_dataset_version(self.id)
        elif getattr(self, name_id) is not None:
          setattr(self, _name, getattr(self.app, "get_" + name)(**{name_id:getattr(self, name_id)}))
      return getattr(self, _name)
    elif name == "scope":
      return task_scope[task_types.index(self.type)]
    else:
      return super(Job, self).__getattr__(name)

  def __setattr__(self, name, value):
    if name in ["run", "rungroup", "trial", "task", "dataset", "dataset_version"]:
      setattr(self, "_"+name, value)
    else:
      super(Job, self).__setattr__(name, value)

  def get_log_path(self):
    run_path = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run)
    return os.path.join(run_path, "stdout", "out.log")

  def submit(self, previous_job = None):
    raise NotImplementedError("Override me!")

  def delete(self, output_only=False):
    raise NotImplementedError("Override me!")

  def get_output_files(self):
    # Retrun folder and experiment and reflection table suffixes
    raise NotImplementedError("Override me!")

  def remove_from_db(self):
    assert self.status == "DELETED"

    print("Removing job %d from the db"%self.id, end=' ')
    tag = self.app.params.experiment_tag
    query = """DELETE job FROM `%s_job` job
               WHERE job.id = %d""" % (
               tag, self.id)
    cursor = self.app.execute_query(query, commit=True)
    print("(%d)"%cursor.rowcount)

  def get_identifier_string(self):
    if self.app.params.facility.name == 'lcls':
      s =  "%s_%s_r%04d_t%03d_rg%03d"% \
        (self.app.params.facility.lcls.experiment, self.app.params.experiment_tag, int(self.run.run), self.trial.trial, self.rungroup.id)
    else:
      s =  "%s_%s_t%03d_rg%03d"% \
        (self.app.params.experiment_tag, self.run.run, self.trial.trial, self.rungroup.id)
    if self.task is not None:
      s += "_task%03d"%self.task.id
    return s

class AveragingJob(Job):
  def get_identifier_string(self):
    # Override this function because rungroups are not used for averaging
    if self.app.params.facility.name == 'lcls':
      s = "%s_%s_r%04d"% \
        (self.app.params.facility.lcls.experiment, self.app.params.experiment_tag, int(self.run.run))
    else:
      s = "%s_%s"% \
        (self.app.params.experiment_tag, self.run.run)
    return s

  def submit(self, previous_job=None):
    from xfel.command_line.cxi_mpi_submit import Script as submit_script
    params = copy.deepcopy(self.app.params)
    params.dispatcher = 'dxtbx.image_average'
    configs_dir = os.path.join(settings_dir, "cfgs")
    if not os.path.exists(configs_dir):
      os.makedirs(configs_dir)
    identifier_string = self.get_identifier_string()

    # Make an argument list that can be submitted to cxi_mpi_submit.
    # dxtbx.image_average does not use phil files.
    extra_args = "-a <output_dir>/avg.cbf -m <output_dir>/max.cbf -s <output_dir>/std.cbf"
    if self.skip_images > 0:
      extra_args += f' --skip-images={self.skip_images}'
    if self.num_images > 0:
      extra_args += f' --num-images={self.num_images}'
    self.args = [
      f'input.run_num = {self.run.run}',
      'input.dispatcher = dxtbx.image_average',
      'output.output_dir = {0}'.format(os.path.join(params.output_folder, 'averages')),
      'output.split_logs = False',
      'output.add_output_dir_option = False',
      f'mp.extra_args = {extra_args}',
      f'mp.method = {params.mp.method}',
    ]
    for opt in params.mp.extra_options:
      self.args.append(f'mp.extra_options = {opt}')

    if params.mp.method != 'local' or (params.mp.method == 'local' and params.facility.name == 'lcls'):
      mp_args = [
        f'mp.use_mpi = {params.mp.use_mpi}',
        f'mp.mpi_command = {params.mp.mpi_command}',
        f'mp.mpi_option = "--mpi=True"',
        f'mp.nnodes = {params.mp.nnodes}',
        f'mp.nproc = {params.mp.nproc}',
        f'mp.nproc_per_node = {params.mp.nproc_per_node}',
        f'mp.queue = {params.mp.queue}',
        f'mp.env_script = {params.mp.env_script[0]}',
        f'mp.wall_time = {params.mp.wall_time}',
        f'mp.htcondor.executable_path = {params.mp.htcondor.executable_path}',
      ]
      for arg in mp_args:
        self.args.append(arg)
    if params.mp.shifter.shifter_image is not None:
      shifter_args = [
        f'mp.shifter.nersc_shifter_image = {params.mp.shifter.shifter_image}',
        f'mp.shifter.sbatch_script_template = {params.mp.shifter.sbatch_script_template}',
        f'mp.shifter.srun_script_template = {params.mp.shifter.srun_script_template}',
        f'mp.shifter.nersc_partition = {params.mp.shifter.partition}',
        f'mp.shifter.nersc_jobname = {params.mp.shifter.jobname}',
        f'mp.shifter.nersc_project = {params.mp.shifter.project}',
        f'mp.shifter.nersc_constraint = {params.mp.shifter.constraint}',
        f'mp.shifter.nersc_reservation = {params.mp.shifter.reservation}',
        f'mp.shifter.staging = {params.mp.shifter.staging}',
      ]
      for arg in shifter_args:
        self.args.append(arg)

    if params.facility.name == 'lcls':
      locator_path = os.path.join(configs_dir, identifier_string + ".loc")
      self.args.append(f'input.locator = {locator_path}')
      write_xtc_locator(locator_path, params, self.run, self.rungroup)
    else:
      self.args.append(self.run.path)
    result = submit_script().run(self.args)
    return result

class IndexingJob(Job):
  def get_output_files(self):
    run_path = str(get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run))
    return os.path.join(run_path, 'out'), '_integrated.expt', '_integrated.refl', None, None

  def submit(self, previous_job = None):
    import libtbx.load_env
    configs_dir = os.path.join(settings_dir, "cfgs")
    if not os.path.exists(configs_dir):
      os.makedirs(configs_dir)
    identifier_string = self.get_identifier_string()

    target_phil_path = os.path.join(configs_dir, identifier_string + "_params.phil")
    dispatcher = self.app.params.dispatcher
    phil_str = self.trial.target_phil_str
    if phil_str is None: phil_str = ""
    if self.rungroup.extra_phil_str is not None:
      phil_str += "\n" + self.rungroup.extra_phil_str

    from xfel.ui import load_phil_scope_from_dispatcher
    if dispatcher == "cxi.xtc_process":
      image_format = 'pickle'
    else:
      orig_phil_scope = load_phil_scope_from_dispatcher(dispatcher)
      if os.path.isfile(dispatcher):
        dispatcher = 'libtbx.python ' + dispatcher
      from iotbx.phil import parse
      if self.rungroup.two_theta_low is not None or self.rungroup.two_theta_high is not None:
        override_str = """
        radial_average {
          enable = True
          show_plots = False
          verbose = False
          output_bins = False
          mask = %s
        }
        """%(self.rungroup.untrusted_pixel_mask_path)
        phil_scope = orig_phil_scope.fetch(parse(override_str))
      else:
        phil_scope = orig_phil_scope

      trial_params = phil_scope.fetch(parse(phil_str)).extract()

      image_format = self.rungroup.format
      mode = get_image_mode(self.rungroup)

      if hasattr(trial_params, 'format'):
        trial_params.format.file_format = image_format
        trial_params.format.cbf.mode = mode

      if hasattr(trial_params.indexing.stills, 'known_orientations') and \
         len(trial_params.indexing.stills.known_orientations) == 1:
        try:
          ko_trial, ko_rungroup = trial_params.indexing.stills.known_orientations[0].split('_')
          ko_trial = self.app.get_trial(trial_number=int(ko_trial))
          ko_rungroup = self.app.get_rungroup(int(ko_rungroup.lstrip('rg')))
        except (IndexError, ValueError):
          pass
        else:
          ko_run_path = get_run_path(self.app.params.output_folder, ko_trial, ko_rungroup, self.run)
          ko_wildcard = trial_params.output.refined_experiments_filename.replace('%s', '*')
          trial_params.indexing.stills.known_orientations[0] = os.path.join(ko_run_path, 'out', ko_wildcard)

    if self.rungroup.calib_dir is not None or self.rungroup.config_str is not None or dispatcher == 'cxi.xtc_process' or image_format == 'pickle':
      config_path = os.path.join(configs_dir, identifier_string + ".cfg")
    else:
      config_path = None

    if hasattr(trial_params.dispatch, 'process_percent'):
      trial_params.dispatch.process_percent = self.trial.process_percent

    # Dictionary for formating the submit phil and, if used, the labelit cfg file
    d = dict(
      # Generally for the LABELIT backend or image pickles
      address                   = self.rungroup.detector_address,
      default_calib_dir         = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0"),
      dark_avg_path             = self.rungroup.dark_avg_path,
      dark_stddev_path          = self.rungroup.dark_stddev_path,
      untrusted_pixel_mask_path = self.rungroup.untrusted_pixel_mask_path,
      detz_parameter            = self.rungroup.detz_parameter,
      gain_map_path             = self.rungroup.gain_map_path,
      gain_mask_level           = self.rungroup.gain_mask_level,
      beamx                     = self.rungroup.beamx,
      beamy                     = self.rungroup.beamy,
      energy                    = self.rungroup.energy,
      binning                   = self.rungroup.binning,
      two_theta_low             = self.rungroup.two_theta_low,
      two_theta_high            = self.rungroup.two_theta_high,
      # Generally for job submission
      dry_run                   = self.app.params.dry_run,
      dispatcher                = dispatcher,
      cfg                       = config_path,
      experiment                = self.app.params.facility.lcls.experiment, # LCLS specific parameter
      run_num                   = self.run.run,
      output_dir                = self.app.params.output_folder,
      use_ffb                   = self.app.params.facility.lcls.use_ffb, # LCLS specific parameter
      # Generally for both
      trial                     = self.trial.trial,
      rungroup                  = self.rungroup.rungroup_id,
      experiment_tag            = self.app.params.experiment_tag,
      calib_dir                 = self.rungroup.calib_dir,
      nproc                     = self.app.params.mp.nproc,
      nnodes                    = self.app.params.mp.nnodes_index or self.app.params.mp.nnodes,
      nproc_per_node            = self.app.params.mp.nproc_per_node,
      queue                     = self.app.params.mp.queue or None,
      env_script                = self.app.params.mp.env_script[0] if self.app.params.mp.env_script is not None and len(self.app.params.mp.env_script) > 0 and len(self.app.params.mp.env_script[0]) > 0 else None,
      phenix_script             = self.app.params.mp.phenix_script[0] if self.app.params.mp.phenix_script is not None and len(self.app.params.mp.phenix_script) > 0 and len(self.app.params.mp.phenix_script[0]) > 0 else None,
      method                    = self.app.params.mp.method,
      wall_time                 = self.app.params.mp.wall_time,
      htcondor_executable_path  = self.app.params.mp.htcondor.executable_path,
      nersc_shifter_image       = self.app.params.mp.shifter.shifter_image,
      sbatch_script_template    = self.app.params.mp.shifter.sbatch_script_template,
      srun_script_template      = self.app.params.mp.shifter.srun_script_template,
      nersc_partition           = self.app.params.mp.shifter.partition,
      nersc_jobname             = self.app.params.mp.shifter.jobname,
      nersc_project             = self.app.params.mp.shifter.project,
      nersc_constraint          = self.app.params.mp.shifter.constraint,
      nersc_reservation         = self.app.params.mp.shifter.reservation,
      nersc_staging             = self.app.params.mp.shifter.staging,
      target                    = target_phil_path,
      host                      = self.app.params.db.host,
      dbname                    = self.app.params.db.name,
      user                      = self.app.params.db.user,
      port                      = self.app.params.db.port,
      # always use mpi for 'lcls'
      use_mpi                   = self.app.params.mp.method != 'local' or (self.app.params.mp.method == 'local' and self.app.params.facility.name == 'lcls'),
      mpi_command               = self.app.params.mp.mpi_command,
      extra_options             = "\n".join(["extra_options = %s"%opt for opt in self.app.params.mp.extra_options]),
    )
    if self.app.params.mp.method == 'sge':
      d['use_mpi'] = False
    if self.app.params.db.password is not None and len(self.app.params.db.password) == 0:
      d['password'] = None
    else:
      d['password'] = self.app.params.db.password

    phil = open(target_phil_path, "w")

    if dispatcher == 'cxi.xtc_process':
      phil.write(phil_str)
    else:
      extra_scope = None
      if hasattr(trial_params, 'format'):
        if image_format == "cbf":
          trial_params.input.address = self.rungroup.detector_address
          trial_params.format.cbf.detz_offset = self.rungroup.detz_parameter
          trial_params.format.cbf.override_energy = self.rungroup.energy
          trial_params.format.cbf.invalid_pixel_mask = self.rungroup.untrusted_pixel_mask_path
          if mode == 'cspad':
            trial_params.format.cbf.cspad.gain_mask_value = self.rungroup.gain_mask_level
          elif mode == 'rayonix':
            trial_params.format.cbf.rayonix.bin_size = self.rungroup.binning
            trial_params.format.cbf.rayonix.override_beam_x = self.rungroup.beamx
            trial_params.format.cbf.rayonix.override_beam_y = self.rungroup.beamy

        if trial_params.input.known_orientations_folder is not None:
          trial_params.input.known_orientations_folder = trial_params.input.known_orientations_folder.format(run=self.run.run)
      else:
        if trial_params.spotfinder.lookup.mask is None:
          trial_params.spotfinder.lookup.mask = self.rungroup.untrusted_pixel_mask_path
        if trial_params.integration.lookup.mask is None:
          trial_params.integration.lookup.mask = self.rungroup.untrusted_pixel_mask_path

        if self.app.params.facility.name == 'lcls':
          locator_path = os.path.join(configs_dir, identifier_string + ".loc")
          write_xtc_locator(locator_path, self.app.params, self.run, self.rungroup)
          if mode == 'rayonix':
            from xfel.cxi.cspad_ana import rayonix_tbx
            pixel_size = rayonix_tbx.get_rayonix_pixel_size(self.rungroup.binning)
            extra_scope = parse("geometry { detector { panel { origin = (%f, %f, %f) } } }"%(-self.rungroup.beamx * pixel_size,
                                                                                              self.rungroup.beamy * pixel_size,
                                                                                             -self.rungroup.detz_parameter))
          d['locator'] = locator_path
        else:
          d['locator'] = None

      if self.rungroup.two_theta_low is not None or self.rungroup.two_theta_high is not None:
        try:
          trial_params.radial_average.two_theta_low = self.rungroup.two_theta_low
          trial_params.radial_average.two_theta_high = self.rungroup.two_theta_high
        except AttributeError:
          pass # not all dispatchers support radial averaging

      working_phil = phil_scope.format(python_object=trial_params)
      if extra_scope:
        working_phil = working_phil.fetch(extra_scope)
      diff_phil = orig_phil_scope.fetch_diff(source=working_phil)

      phil.write(diff_phil.as_str())
    phil.close()

    if config_path is not None:
      if dispatcher != 'cxi.xtc_process':
        d['untrusted_pixel_mask_path'] = None # Don't pass a pixel mask to mod_image_dict as it will
                                              # will be used during dials processing directly

      config_str = "[psana]\n"
      if self.rungroup.calib_dir is not None:
        config_str += "calib-dir=%s\n"%self.rungroup.calib_dir
      modules = []
      if self.rungroup.config_str is not None:
        for line in self.rungroup.config_str.split("\n"):
          if line.startswith('['):
            modules.append(line.lstrip('[').rstrip(']'))
      if dispatcher == 'cxi.xtc_process':
        modules.insert(0, 'my_ana_pkg.mod_radial_average')
        modules.extend(['my_ana_pkg.mod_hitfind:index','my_ana_pkg.mod_dump:index'])
      elif image_format == 'pickle':
        modules.insert(0, 'my_ana_pkg.mod_radial_average')
        modules.extend(['my_ana_pkg.mod_image_dict'])
      if self.app.params.facility.lcls.dump_shots:
        modules.insert(0, 'my_ana_pkg.mod_dump:shot')

      if len(modules) > 0:
        config_str += "modules = %s\n"%(" ".join(modules))

      if self.rungroup.config_str is not None:
        config_str += self.rungroup.config_str + "\n"

      if dispatcher == 'cxi.xtc_process' or image_format == 'pickle':
        d['address'] = d['address'].replace('.','-').replace(':','|') # old style address
        if dispatcher == 'cxi.xtc_process':
          template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "index_all.cfg"))
        elif image_format == 'pickle':
          template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "image_dict.cfg"))
        for line in template.readlines():
          config_str += line.format(**d)
        template.close()
        d['address'] = self.rungroup.detector_address

      cfg = open(config_path, 'w')
      cfg.write(config_str)
      cfg.close()

      if dispatcher != 'cxi.xtc_process':
        d['untrusted_pixel_mask_path'] = self.rungroup.untrusted_pixel_mask_path

    submit_phil_path = os.path.join(configs_dir, identifier_string + "_submit.phil")
    submit_root = libtbx.env.find_in_repositories("xfel/ui/db/cfgs")
    if dispatcher in ['cxi.xtc_process', 'cctbx.xfel.xtc_process']:
      template = open(os.path.join(submit_root, "submit_xtc_process.phil"))
    else:
      test_root = os.path.join(submit_root, "submit_" + dispatcher + ".phil")
      if os.path.exists(test_root):
        template = open(test_root)
      else:
        if hasattr(trial_params, 'format'):
          template = open(os.path.join(submit_root, "submit_xtc_process.phil"))
        else:
          template = open(os.path.join(submit_root, "submit_xfel_process.phil"))
    phil = open(submit_phil_path, "w")

    if dispatcher == 'cxi.xtc_process':
      d['target'] = None # any target phil will be in mod_hitfind

    for line in template.readlines():
      phil.write(line.format(**d))

    d['target'] = target_phil_path

    template.close()
    phil.close()

    from xfel.command_line.cxi_mpi_submit import Script as submit_script
    args = [submit_phil_path]
    if self.app.params.facility.name not in ['lcls']:
      args.append(self.run.path)
    return submit_script().run(args)

  def delete(self, output_only=False):
    if self.status not in finished_job_statuses:
      print("Job is not finished (status = %s)"%self.status)
      return

    if self.status == "DELETED":
      return

    job_folder = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run)
    if os.path.exists(job_folder):
      print("Deleting job folder for job", self.id)
      shutil.rmtree(job_folder)
    else:
      print("Cannot find job folder (%s)"%job_folder)

    # Have to be careful to delete from the tables in the right order
    tag = self.app.params.experiment_tag

    def delete_and_commit(query):
      cursor = self.app.execute_query(query, commit=True)
      print("(%d)"%cursor.rowcount)

    print("Deleting cell_bin entries", end=' ')
    query = """DELETE cell_bin FROM `%s_cell_bin` cell_bin
               JOIN `%s_crystal` crystal ON crystal.id = cell_bin.crystal_id
               JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    ids = {}
    for item in "crystal", "beam", "detector":
      print("Listing %s ids"%item, end=' ')
      query = """SELECT %s.id FROM `%s_%s` %s
                 JOIN `%s_experiment` expr ON expr.%s_id = %s.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
                 item, tag, item, item, tag, item, item, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print("(%d)"%len(item_ids))
      ids[item] = ",".join(item_ids)

    if len(self.trial.isoforms) == 0:
      print("Listing bin entries", end=' ')
      query = """SELECT bin.id FROM `%s_bin` bin
                 JOIN `%s_cell` cell ON bin.cell_id = cell.id
                 JOIN `%s_crystal` crystal ON crystal.cell_id = cell.id
                 JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d
                 AND cell.trial_id is NULL""" % (
                 tag, tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print("(%d)"%len(item_ids))
      bin_ids = ",".join(item_ids)

      print("Listing cell entries", end=' ')
      query = """SELECT cell.id FROM `%s_cell` cell
                 JOIN `%s_crystal` crystal ON crystal.cell_id = cell.id
                 JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d
                 AND cell.trial_id IS NULL""" % (
                 tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print("(%d)"%len(item_ids))
      cell_ids = ",".join(item_ids)

    print("Deleting experiment entries", end=' ')
    query = """DELETE expr FROM `%s_experiment` expr
               JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    for item in "crystal", "beam", "detector":
      if len(ids[item]) > 0:
        print("Deleting %s entries"%item, end=' ')
        query = """DELETE %s FROM `%s_%s` %s
                   WHERE %s.id IN (%s)""" % (
                   item, tag, item, item, item, ids[item])
        delete_and_commit(query)

    if len(self.trial.isoforms) == 0 and len(bin_ids) > 0:
      print("Deleting bin entries", end=' ')
      query = """DELETE bin FROM `%s_bin` bin
                 WHERE bin.id IN (%s)""" % (
                 tag, bin_ids)
      delete_and_commit(query)

    if len(self.trial.isoforms) == 0 and len(cell_ids) > 0:
      print("Deleting cell entries", end=' ')
      query = """DELETE cell FROM `%s_cell` cell
                 WHERE cell.id IN (%s)""" % (
                 tag, cell_ids)
      delete_and_commit(query)

    print("Listing imageset entries", end=' ')
    query = """SELECT imgset.id FROM `%s_imageset` imgset
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    cursor = self.app.execute_query(query)
    item_ids = ["%d"%i[0] for i in cursor.fetchall()]
    print("(%d)"%len(item_ids))
    imageset_ids = ",".join(item_ids)

    print("Deleting imageset_event entries", end=' ')
    query = """DELETE is_e FROM `%s_imageset_event` is_e
               JOIN `%s_event` evt ON evt.id = is_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    if len(imageset_ids) > 0:
      print("Deleting imageset entries", end=' ')
      query = """DELETE imgset FROM `%s_imageset` imgset
                 WHERE imgset.id IN (%s)""" % (
                 tag, imageset_ids)
      delete_and_commit(query)

    print("Deleting event entries", end=' ')
    query = """DELETE evt FROM `%s_event` evt
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    self.status = "DELETED"

class EnsembleRefinementJob(Job):
  def delete(self, output_only=False):
    job_folder = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    if os.path.exists(job_folder):
      print("Deleting job folder for job", self.id)
      shutil.rmtree(job_folder)
    else:
      print("Cannot find job folder (%s)"%job_folder)
    self.status = "DELETED"

  def get_output_files(self):
    run_path = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    return os.path.join(run_path, 'combine_experiments_t%03d'%self.trial.trial, 'intermediates', "*reintegrated*"), '.expt', '.refl', None, None

  def get_log_path(self):
    run_path = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    return os.path.join(run_path, 'combine_experiments_t%03d'%self.trial.trial, 'intermediates',
      "combine_t%03d_rg%03d_chunk000.out"%(self.trial.trial, self.rungroup.id)) # XXX there can be multiple chunks or multiple clusters

  def submit(self, previous_job = None):
    from xfel.command_line.striping import Script
    from xfel.command_line.cxi_mpi_submit import get_submission_id
    from libtbx import easy_run
    configs_dir = os.path.join(settings_dir, "cfgs")
    identifier_string = self.get_identifier_string()
    target_phil_path = os.path.join(configs_dir, identifier_string + "_params.phil")
    with open(target_phil_path, 'w') as f:
      if self.task.parameters:
        f.write(self.task.parameters)

    path = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    os.mkdir(path)

    arguments = """
    mp.queue={}
    mp.nnodes={}
    mp.nproc_per_node={}
    mp.method={}
    {}
    {}
    mp.wall_time={}
    mp.use_mpi=False
    mp.mpi_command={}
    {}
    mp.shifter.submit_command={}
    mp.shifter.shifter_image={}
    mp.shifter.sbatch_script_template={}
    mp.shifter.srun_script_template={}
    mp.shifter.partition={}
    mp.shifter.jobname={}
    mp.shifter.project={}
    mp.shifter.reservation={}
    mp.shifter.constraint={}
    mp.shifter.staging={}
    striping.results_dir={}
    striping.trial={}
    striping.rungroup={}
    striping.run={}
    {}
    striping.chunk_size=64
    striping.stripe=False
    striping.dry_run=True
    striping.output_folder={}
    reintegration.integration.lookup.mask={}
    mp.local.include_mp_in_command=False
    """.format(self.app.params.mp.queue if len(self.app.params.mp.queue) > 0 else None,
               self.app.params.mp.nnodes_tder or self.app.params.mp.nnodes,
               self.app.params.mp.nproc_per_node,
               self.app.params.mp.method,
               '\n'.join(['mp.env_script={}'.format(p) for p in self.app.params.mp.env_script if p]),
               '\n'.join(['mp.phenix_script={}'.format(p) for p in self.app.params.mp.phenix_script if p]),
               self.app.params.mp.wall_time,
               self.app.params.mp.mpi_command,
               "\n".join(["mp.extra_options={}".format(opt) for opt in self.app.params.mp.extra_options]),
               self.app.params.mp.shifter.submit_command,
               self.app.params.mp.shifter.shifter_image,
               self.app.params.mp.shifter.sbatch_script_template,
               self.app.params.mp.shifter.srun_script_template,
               self.app.params.mp.shifter.partition,
               self.app.params.mp.shifter.jobname,
               self.app.params.mp.shifter.project,
               self.app.params.mp.shifter.reservation,
               self.app.params.mp.shifter.constraint,
               self.app.params.mp.shifter.staging,
               self.app.params.output_folder,
               self.trial.trial,
               self.rungroup.id,
               self.run.run,
               target_phil_path,
               path,
               self.rungroup.untrusted_pixel_mask_path,
               ).split('\n')
    arguments = [arg.strip() for arg in arguments]

    try:
      commands = Script(arguments).run()
    except Exception as e:
      if 'no DIALS integration results found' in str(e):
        print("No DIALS integration results found")
        self.status = "EXIT"
        return
      else: raise
    submission_ids = []
    if self.app.params.mp.method == 'local':
      self.status = "RUNNING"
    for command in commands:
      try:
        result = easy_run.fully_buffered(command=command)
        result.raise_if_errors()
      except Exception as e:
        if not "Warning: job being submitted without an AFS token." in str(e):
          raise e
      submission_ids.append(get_submission_id(result, self.app.params.mp.method))
    if self.app.params.mp.method == 'local':
      self.status = "DONE"
    else:
      return ",".join(submission_ids)

class ScalingJob(Job):
  def delete(self, output_only=False):
    job_folder = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    if os.path.exists(job_folder):
      print("Deleting job folder for job", self.id)
      shutil.rmtree(job_folder)
    else:
      print("Cannot find job folder (%s)"%job_folder)
    self.status = "DELETED"

  def get_output_files(self):
    run_path = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task)
    return os.path.join(run_path, 'out'), ".expt", ".refl", None, None

  def write_submit_phil(self, submit_phil_path, target_phil_path):
    import libtbx.load_env
    from xfel.ui.db.task import task_types, task_dispatchers

    submit_root = libtbx.env.find_in_repositories("xfel/ui/db/cfgs")
    d =  dict(
      dry_run                   = self.app.params.dry_run,
      dispatcher                = task_dispatchers[task_types.index(self.task.type)],
      run_num                   = self.run.run,
      output_dir                = self.app.params.output_folder,
      trial                     = self.trial.trial,
      rungroup                  = self.rungroup.rungroup_id,
      task                      = self.task.id,
      nproc                     = self.app.params.mp.nproc,
      nproc_per_node            = self.app.params.mp.nproc_per_node,
      queue                     = self.app.params.mp.queue or None,
      env_script                = self.app.params.mp.env_script[0] if len(self.app.params.mp.env_script) > 0 and len(self.app.params.mp.env_script[0]) > 0 else None,
      phenix_script                = self.app.params.mp.phenix_script[0] if len(self.app.params.mp.phenix_script) > 0 and len(self.app.params.mp.phenix_script[0]) > 0 else None,
      method                    = self.app.params.mp.method,
      htcondor_executable_path  = self.app.params.mp.htcondor.executable_path,
      nersc_shifter_image       = self.app.params.mp.shifter.shifter_image,
      sbatch_script_template    = self.app.params.mp.shifter.sbatch_script_template,
      srun_script_template      = self.app.params.mp.shifter.srun_script_template,
      nersc_partition           = self.app.params.mp.shifter.partition,
      nersc_jobname             = self.app.params.mp.shifter.jobname,
      nersc_project             = self.app.params.mp.shifter.project,
      nersc_constraint          = self.app.params.mp.shifter.constraint,
      nersc_reservation         = self.app.params.mp.shifter.reservation,
      nersc_staging             = self.app.params.mp.shifter.staging,
      target                    = target_phil_path,
      # always use mpi for 'lcls'
      use_mpi                   = self.app.params.mp.method != 'local' or (self.app.params.mp.method == 'local' and self.app.params.facility.name == 'lcls'),
      mpi_command               = self.app.params.mp.mpi_command,
      nnodes                    = self.app.params.mp.nnodes_scale or self.app.params.mp.nnodes,
      wall_time                 = self.app.params.mp.wall_time,
      extra_options             = "\n".join(["extra_options = %s"%opt for opt in self.app.params.mp.extra_options]),
    )

    with open(submit_phil_path, "w") as phil:
      for line in open(os.path.join(submit_root, "submit_xfel_merge.phil")).readlines():
        phil.write(line.format(**d))

  def submit(self, previous_job = None):
    from xfel.command_line.cxi_mpi_submit import Script as submit_script

    output_path = os.path.join(get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run, self.task), 'out')

    configs_dir = os.path.join(settings_dir, "cfgs")
    if not os.path.exists(configs_dir):
      os.makedirs(configs_dir)
    identifier_string = self.get_identifier_string()
    submit_phil_path = os.path.join(configs_dir, identifier_string + "_submit.phil")

    target_phil_path = os.path.join(configs_dir, identifier_string + "_params.phil")
    input_folder, expt_suffix, refl_suffix, _, _ = previous_job.get_output_files()

    with open(target_phil_path, 'w') as f:
      f.write("input.path=%s\n"%input_folder)
      f.write("input.experiments_suffix=%s\n"%expt_suffix)
      f.write("input.reflections_suffix=%s\n"%refl_suffix)
      f.write("output.output_dir=%s\n"%output_path)
      f.write("output.prefix=%s_%d\n"%(self.task.type, self.task.id))
      f.write(self.task.parameters)

    self.write_submit_phil(submit_phil_path, target_phil_path)

    args = [submit_phil_path]
    return submit_script().run(args)

class MergingJob(Job):
  def get_global_path(self):
    if self.dataset_version is None:
      return None
    return self.dataset_version.output_path()

  def get_log_path(self):
    return self.get_global_path()

  def get_identifier_string(self):
    return "%s_%s%03d_v%03d"%(self.dataset.name, self.task.type, self.task.id, self.dataset_version.version)

  def delete(self, output_only=False):
    job_folder = self.get_global_path()
    if job_folder and os.path.exists(job_folder):
      print("Deleting job folder for job", self.id)
      shutil.rmtree(job_folder)
    else:
      print("Cannot find job folder (%s)"%job_folder)
    self.status = "DELETED"

  def get_output_files(self):
    path = self.get_global_path()
    return path, None, None, "%s_v%03d_all.mtz"%(self.dataset.name, self.dataset_version.version), None

  def submit(self, previous_job = None):
    from xfel.command_line.cxi_mpi_submit import do_submit

    output_path = self.get_global_path()
    if not os.path.exists(output_path):
      os.makedirs(output_path)
    identifier_string = self.get_identifier_string()
    target_phil_path = os.path.join(output_path, identifier_string + "_params.phil")

    with open(target_phil_path, 'w') as f:
      expt_suffix = refl_suffix = None
      for job in self.dataset_version.jobs:
        input_folder, _expt_suffix, _refl_suffix, _, _ = job.get_output_files()
        if expt_suffix is None: expt_suffix = _expt_suffix
        else: assert expt_suffix == _expt_suffix
        if refl_suffix is None: refl_suffix = _refl_suffix
        else: assert refl_suffix == _refl_suffix
        f.write("input.path=%s\n"%input_folder)

      f.write("input.experiments_suffix=%s\n"%expt_suffix)
      f.write("input.reflections_suffix=%s\n"%refl_suffix)
      f.write("output.output_dir=%s\n"%output_path)
      f.write("output.prefix=%s_v%03d\n"%(self.dataset.name, self.dataset_version.version))
      f.write(self.task.parameters)

    command = "cctbx.xfel.merge %s"%target_phil_path
    submit_path = os.path.join(output_path, identifier_string + "_submit.sh")

    params = self.app.params.mp
    if params.nnodes_merge:
      params = copy.deepcopy(params)
      params.nnodes = params.nnodes_merge

    return do_submit(command, submit_path, output_path, params, log_name="out.log", err_name="err.log", job_name=identifier_string)

class PhenixJob(Job):
  def get_global_path(self):
    return os.path.join(self.dataset_version.output_path(), self.get_identifier_string())

  def get_log_path(self):
    return self.get_global_path()

  def get_identifier_string(self):
    return "%s_%s%03d_v%03d"%(self.dataset.name, self.task.type, self.task.id, self.dataset_version.version)

  def delete(self, output_only=False):
    try:
      job_folder = self.get_global_path()
    except AttributeError: # If job fails to submit for example
      job_folder = None
    if job_folder and os.path.exists(job_folder):
      print("Deleting job folder for job", self.id)
      shutil.rmtree(job_folder)
    else:
      print("Cannot find job folder (%s)"%job_folder)
    self.status = "DELETED"

  def get_output_files(self):
    path = self.get_global_path()
    return path, None, None, ".mtz", ".pdb"

  def submit(self, previous_job = None):
    from xfel.command_line.cxi_mpi_submit import do_submit

    output_path = self.get_global_path()
    if not os.path.exists(output_path):
      os.makedirs(output_path)
    identifier_string = self.get_identifier_string()
    target_phil_path = os.path.join(output_path, identifier_string + "_params.phil")
    input_folder, _, _, input_mtz, _ = previous_job.get_output_files()

    def replace_keywords(input_str):
      input_str = input_str.replace('<PREVIOUS_TASK_MTZ>', os.path.join(input_folder, input_mtz))
      input_str = input_str.replace('<PREVIOUS_TASK_FOLDER>', input_folder)
      input_str = input_str.replace('<DATASET_NAME>', self.dataset.name)
      input_str = input_str.replace('<DATASET_VERSION>', str(self.dataset_version.version))
      return input_str

    command = replace_keywords(self.task.parameters.split('\n')[0])
    phil_params = replace_keywords('\n'.join(self.task.parameters.split('\n')[1:]))

    with open(target_phil_path, 'w') as f:
      f.write(phil_params)

    command = "%s %s"%(command, target_phil_path)
    submit_path = os.path.join(output_path, identifier_string + "_submit.sh")

    params = copy.deepcopy(self.app.params.mp)
    if params.nnodes_merge:
      params.nnodes = params.nnodes_merge
    params.use_mpi = False
    params.shifter.staging = None
    if not('upload' in command or 'evaluate_anom' in command):
      params.env_script = params.phenix_script
      if not params.extra_options:
        params.extra_options = ""
      params.extra_options.append("--export=None")
    params.nnodes = 1
    params.nproc_per_node = 1

    if params.method == 'shifter' and 'upload' not in command:
       import libtbx.load_env
       params.shifter.sbatch_script_template = os.path.join( \
         libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "phenix_sbatch.sh")
       params.shifter.srun_script_template = os.path.join( \
         libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "phenix_srun.sh")

    return do_submit(command, submit_path, output_path, params, log_name="out.log", err_name="err.log", job_name=identifier_string)

# Support classes and functions for job submission

class _job(object):
  """Used to represent a job that may not have been submitted into the cluster or database yet"""
  def __init__(self, trial, rungroup, run, task=None, dataset=None, dataset_version=None):
    self.trial = trial
    self.rungroup = rungroup
    self.run = run
    self.task = task
    self.dataset = dataset
    self.dataset_version = dataset_version

  def __str__(self):
    s = "Job: Trial %s, rg %s, run %s"%(str(self.trial.trial) if self.trial else 'None',
                                        str(self.rungroup.id) if self.rungroup else 'None',
                                        self.run.run if self.run else 'None')
    if self.task:
      s += ", task %d %s"%(self.task.id, self.task.type)
    if self.dataset:
      s += ", dataset %d %s"%(self.dataset.id, self.dataset.name)
    if self.dataset_version:
      s += ", dataset_version %d %d"%(self.dataset_version.id, self.dataset_version.version)
    return s

  @staticmethod
  def job_hash(job):
    ret = []
    check = ['trial', 'rungroup', 'run', 'task', 'dataset', 'dataset_version']
    for subitem_name in check:
      subitem = getattr(job, subitem_name)
      if subitem is None:
        ret.append(None)
      else:
        ret.append(subitem.id)
    return tuple(ret)

  def __eq__(self, other):
    return _job.job_hash(self) == _job.job_hash(other)

def submit_all_jobs(app):
  submitted_jobs = {_job.job_hash(j):j for j in app.get_all_jobs()}
  if app.params.mp.method == 'local': # only run one job at a time
    for job in submitted_jobs.values():
      if job.status in ['RUN', 'UNKWN', 'SUBMITTED']: return

  if app.params.mp.max_queued is not None:
    running_jobs = sum([1 for job in submitted_jobs.values() if job.status in ['RUN', 'SUBMITTED', 'PEND']])
    if running_jobs >= app.params.mp.max_queued:
      print("Waiting for space in the queue to submit next job")
      return

  runs = app.get_all_runs()
  trials = app.get_all_trials(only_active = True)

  needed_jobs = []
  for trial in trials:
    for rungroup in trial.rungroups:
      assert rungroup.active
      for run in rungroup.runs:
        needed_jobs.append(_job(trial, rungroup, run))

  for job in needed_jobs:
    if _job.job_hash(job) in submitted_jobs:
      continue

    print("Submitting job: trial %d, rungroup %d, run %s"%(job.trial.trial, job.rungroup.id, job.run.run))

    j = JobFactory.from_args(app,
                             trial_id = job.trial.id,
                             rungroup_id = job.rungroup.id,
                             run_id = job.run.id,
                             status = "SUBMITTED")
    j.trial = job.trial; j.rungroup = job.rungroup; j.run = job.run
    try:
      j.submission_id = j.submit()
    except Exception as e:
      print("Couldn't submit job:", str(e))
      j.status = "SUBMIT_FAIL"
      raise

    if app.params.mp.method == 'local': # only run one job at a time
      return

    if app.params.mp.max_queued is not None:
      running_jobs += 1
      if running_jobs >= app.params.mp.max_queued:
        print("Waiting for space in the queue to submit next job")
        return

  datasets = app.get_all_datasets()
  for dataset_idx, dataset in enumerate(datasets):
    if not dataset.active: continue

    # one of the tasks will have a trial, otherwise we don't know where to save the data
    trial = None
    tasks = dataset.tasks
    for task in tasks:
      if task.trial is not None:
        if trial is None:
          trial = task.trial
        else:
          assert trial.id == task.trial.id, "Found multiple trials, don't know where to save the results"
    assert trial, "No trial found in task list, don't know where to save the results"
    trial_tags_ids = [t.id for t in trial.tags]
    dataset_tags = [t for t in dataset.tags if t.id in trial_tags_ids]
    if not dataset_tags or len(dataset_tags) < len(dataset.tags): continue
    runs_rungroups = []
    for rungroup in trial.rungroups:
      runs_rungroups.extend([(run, rungroup) for run in app.get_rungroup_runs_by_tags(rungroup, dataset_tags, dataset.tag_operator)])

    # Datasets always start with indexing
    global_tasks = {}
    for run, rungroup in runs_rungroups:
      submit_next_task = False
      last_task_status = ""
      previous_job = None
      for task_idx, task in enumerate(tasks):
        if task.scope == 'global':
          if previous_job.status in ["DONE", "EXIT"]:
            key = (dataset_idx, task_idx)
            if key not in global_tasks:
              global_tasks[key] = []
            global_tasks[key].append(previous_job)
          continue
        if task.type == 'indexing':
          job = _job(trial, rungroup, run)
        else:
          job = _job(trial, rungroup, run, task)
        try:
          submitted_job = submitted_jobs[_job.job_hash(job)]
        except KeyError:
          if not submit_next_task:
            print("Warning, expected to find submitted %s job: trial %d, rungroup %d, run %s, task %d"% \
              (task.type, trial.trial, rungroup.id, run.run, task.id))
            break
        else:
          if not task_idx+1 < len(tasks): break # no more tasks to do after this one
          next_task = tasks[task_idx+1]
          if submitted_job.status not in finished_job_statuses or submitted_job.status == "UNKWN":
            print ("Task %s waiting on job %d (%s) for trial %d, rungroup %d, run %s, task %d" % \
              (next_task.type, submitted_job.id, submitted_job.status, trial.trial, rungroup.id, run.run, next_task.id))
            break
          if submitted_job.status not in ["DONE"]:
            if submitted_job.status != "EXIT":
              print ("Task %s cannot start due to unexpected status for job %d (%s) for trial %d, rungroup %d, run %s, task %d" % \
                (next_task.type, submitted_job.id, submitted_job.status, trial.trial, rungroup.id, run.run, next_task.id))
            break
          if submitted_job.status in ("SUBMIT_FAIL", "DELETED", "UNKWN") and job.task and job.task.type == "ensemble_refinement":
            break # XXX need a better way to indicate that a job has failed and shouldn't go through the pipeline due to no data
          submit_next_task = True
          previous_job = submitted_job
          continue

        print("Submitting %s job: trial %d, rungroup %d, run %s, task %d"% \
          (task.type, trial.trial, rungroup.id, run.run, task.id))

        j = JobFactory.from_args(app,
                                 trial_id = trial.id,
                                 rungroup_id = rungroup.id,
                                 run_id = run.id,
                                 task_id = task.id,
                                 status = "SUBMITTED")
        j.trial = job.trial; j.rungroup = job.rungroup; j.run = job.run; j.task = job.task
        try:
          j.submission_id = j.submit(previous_job)
        except Exception as e:
          print("Couldn't submit job:", str(e))
          j.status = "SUBMIT_FAIL"
          raise

        previous_job = j

        if app.params.mp.method == 'local': # only run one job at a time
          return
        if app.params.mp.max_queued is not None:
          running_jobs += 1
          if running_jobs >= app.params.mp.max_queued:
            print("Waiting for space in the queue to submit next job")
            return
        break # job submitted so don't look for more in this run for this dataset

    versions = dataset.versions
    for task_idx, task in enumerate(tasks):
      if task.scope == 'local':
        # only global tasks follow global tasks
        if task_idx: assert tasks[task_idx-1].scope != 'global'
        continue
      assert task.scope == 'global' # only two task scopes
      assert task_idx # first task cannot be global
      prev_task = tasks[task_idx-1]
      if prev_task.scope == 'global':
        # Submit a job for this task for any versions where it has not been
        for version in versions:
          prev_j = _job(None, None, None, prev_task, dataset, version)
          test_j = _job(None, None, None, task, dataset, version)
          prev_job = this_job = None
          for j in version.jobs:
            if prev_j == j:
              prev_job = j
              continue
            elif test_j == j:
              this_job = j
              continue
            if prev_job and this_job: break

          if not this_job and prev_job and prev_job.status == 'DONE':
            j = JobFactory.from_args(app,
                                     task_id = task.id,
                                     dataset_id = dataset.id,
                                     status = "SUBMITTED")
            j.task = task; j.dataset = dataset; j.dataset_version = version

            try:
              j.submission_id = j.submit(prev_job)
            except Exception as e:
              print("Couldn't submit job:", str(e))
              j.status = "SUBMIT_FAIL"
              raise
            version.add_job(j)

            if app.params.mp.method == 'local': # only run one job at a time
              return
            if app.params.mp.max_queued is not None:
              running_jobs += 1
              if running_jobs >= app.params.mp.max_queued:
                print("Waiting for space in the queue to submit next job")
                return

      key = dataset_idx, task_idx
      if key not in global_tasks: continue # no jobs ready yet
      latest_version = dataset.latest_version
      next_version = None
      if latest_version is None:
        next_version = 0
      else:
        latest_version_local_jobs = [j.id for j in latest_version.jobs if j.task and j.task.scope == 'local']
        new_jobs = [j for j in global_tasks[key] if j.id not in latest_version_local_jobs]
        if new_jobs:
          next_version = latest_version.version + 1

      if next_version is not None:
        latest_version = app.create_dataset_version(dataset_id = dataset.id, version=next_version)
        print("Created a new dataset version", next_version, "with jobs", [j.id for j in global_tasks[key]])
        for job in global_tasks[key]:
          latest_version.add_job(job)

        j = JobFactory.from_args(app,
                                 task_id = task.id,
                                 dataset_id = dataset.id,
                                 status = "SUBMITTED")
        j.task = task; j.dataset = dataset; j.dataset_version = latest_version

        try:
          j.submission_id = j.submit()
        except Exception as e:
          print("Couldn't submit job:", str(e))
          j.status = "SUBMIT_FAIL"
          raise
        latest_version.add_job(j)

        if app.params.mp.method == 'local': # only run one job at a time
          return
        if app.params.mp.max_queued is not None:
          running_jobs += 1
          if running_jobs >= app.params.mp.max_queued:
            print("Waiting for space in the queue to submit next job")
            return
