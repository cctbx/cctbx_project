from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 09/18/2019
Description : IOTA GUI Threads and PostEvents
'''

import os
import wx
import shutil
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep
from libtbx import easy_run

from dxtbx.model.experiment_list import ExperimentListFactory
import multiprocessing

import subprocess
import sys
import signal

from xfel.clustering.cluster import Cluster
from cctbx import crystal
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import lattice_symmetry

from iota.components.iota_utils import InputFinder, Capturing, IOTATermination
from iota.components.iota_analysis import Analyzer

ginp = InputFinder()

# for testing
import time
assert time

# -------------------------------- Processing -------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

tp_EVT_IMGDONE = wx.NewEventType()
EVT_IMGDONE = wx.PyEventBinder(tp_EVT_IMGDONE, 1)

tp_EVT_OBJDONE = wx.NewEventType()
EVT_OBJDONE = wx.PyEventBinder(tp_EVT_OBJDONE, 1)

class ImageFinderAllDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, input_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.input_list = input_list
  def GetValue(self):
    return self.input_list

class ObjectFinderAllDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, obj_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.obj_list = obj_list
  def GetValue(self):
    return self.obj_list

class AllDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, img_objects=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_objects = img_objects
  def GetValue(self):
    return self.image_objects

class JobSubmitThread(Thread):
  """ Thread for easy_run submissions so that they don't block GUI """
  def __init__(self, parent, params, out_type='gui_silent'):
    Thread.__init__(self)
    self.parent = parent
    self.params = params
    self.job_id = None

    # Enforce silent run if on mp queue (unless need debug statements)
    if self.params.mp.method != 'multiprocessing' and out_type != 'gui_debug':
      self.out_type = 'gui_silent'
    else:
      self.out_type = out_type

  def submit(self):
    run_path = self.parent.info.int_base
    iota_cmd = 'iota.run --run_path {} -o {}'.format(run_path, self.out_type)
    if self.params.mp.submit_command:
      command = self.params.mp.submit_command.replace('<iota_command>', iota_cmd)
    else:
      command = iota_cmd

    if command is not None:
      print(command)
      self.job = CustomRun(command=str(command), join_stdout_stderr=True)
      self.job.run()
      if self.job_id is not None:
        print('JOB NAME = ', self.job_id)
      return
    else:
      print('IOTA ERROR: COMMAND NOT ISSUED!')
      return

  def abort(self):
    if self.params.mp.kill_command:
      CustomRun(command=self.params.mp.kill_command)
    else:
      try:
        self.job.kill_thread()
      except Exception as e:
        print ('IOTA JOB ERROR: Cannot kill job thread! {}'.format(e))

  def run(self):
    return self.submit()

class ObjectReader():
  def __init__(self):
    pass

  def update_info(self, info):
    # # Determine chunk of image objects to check for existing
    # if self.n_proc:
    #   chunk = self.n_proc
    # else:
    #   chunk = 25
    # if len(self.info.unread_files) < chunk:
    #   chunk = len(self.info.unread_files)
    #
    # # Create filelist of written image objects
    # filelist = []
    # for fp in self.info.unread_files[:chunk]:
    #   if os.path.isfile(fp):
    #     filelist.append(fp)
    #     self.info.unread_files.pop(fp)
    #
    # # Perform stat extraction
    # if filelist:
    #   from iota.components.iota_analysis import Analyzer
    #   analyzer = Analyzer(info=self.info)
    #   stats_OK = analyzer.get_results(filelist=filelist)
    #   if stats_OK:
    #     self.info = analyzer.info
    #     self.obs = analyzer.obs

    finished_objects = info.get_finished_objects_from_file()
    if finished_objects:
      analyzer = Analyzer(info=info, gui_mode=True)
      stats_OK = analyzer.run_get_results(finished_objects=finished_objects)
      if stats_OK:
        return analyzer.info
    return None

  def run(self, info):
    return self.update_info(info)

class ObjectReaderThread(Thread):
  """ Thread for reading processed objects and making all calculations for
      plotting in main GUI thread """
  def __init__(self, parent, info=None):
    Thread.__init__(self, name='object_reader')
    self.parent = parent
    self.info = info
    self.obj_worker = ObjectReader()

  def run(self):
    info = self.obj_worker.run(self.info)
    if info:
      info.export_json()
      evt = ObjectFinderAllDone(tp_EVT_OBJDONE, -1, obj_list=info)
      wx.PostEvent(self.parent, evt)


class ImageFinderThread(Thread):
  """ Worker thread generated to poll filesystem on timer. Will check to see
  if any new images have been found. Put on a thread to run in background """
  def __init__(self,
               parent,
               input,
               input_list,
               min_back=None,
               last_file=None,
               back_to_thread=False):
    Thread.__init__(self)
    self.parent = parent
    self.input = input
    self.input_list = input_list
    self.min_back = min_back
    self.last_file = last_file
    self.back_to_thread = back_to_thread

  def run(self):
    # Poll filesystem and determine which files are new (if any)

    ext_file_list = ginp.make_input_list(self.input,
                                         filter_results=True,
                                         filter_type='image',
                                         min_back=self.min_back,
                                         last=self.last_file)
    new_input_list = list(set(ext_file_list) - set(self.input_list))

    if self.back_to_thread:
      wx.CallAfter(self.parent.onImageFinderDone, new_input_list)
    else:
      evt = ImageFinderAllDone(tp_EVT_IMGDONE, -1, input_list=new_input_list)
      wx.PostEvent(self.parent, evt)

class ObjectFinderThread(Thread):
  """ Worker thread that polls filesystem on timer for image objects. Will
  collect and extract info on images processed so far"""
  def __init__(self,
               parent,
               object_folder,
               last_object=None,
               new_fin_base = None):
    Thread.__init__(self)
    self.parent = parent
    self.object_folder = object_folder
    self.new_fin_base = new_fin_base
    self.last_object = last_object

  def run(self):
    if self.last_object is not None:
      last = self.last_object.obj_file
    else:
      last = None
    object_files = ginp.get_file_list(self.object_folder,
                                      ext_only='int',
                                      last=last)
    new_objects = [self.read_object_file(i) for i in object_files]
    new_finished_objects = [i for i in new_objects if i is not None]

    evt = ObjectFinderAllDone(tp_EVT_OBJDONE, -1, obj_list=new_finished_objects)
    wx.PostEvent(self.parent, evt)

  def read_object_file(self, filepath):
    try:
      object = ep.load(filepath)
      return object
    except EOFError as e:
      print ('OBJECT_IMPORT_ERROR: ', e)
      return None

class ImageViewerThread(Thread):
  """ Worker thread that will move the image viewer launch away from the GUI
  and hopefully will prevent the image selection dialog freezing on MacOS"""
  def __init__(self,
               parent,
               file_string=None,
               viewer='dials.image_viewer',
               img_type=None):
    Thread.__init__(self)
    self.parent = parent
    self.file_string = file_string
    self.viewer = viewer
    self.img_type = img_type

  def run(self):
    command = '{} {}'.format(self.viewer, self.file_string)
    easy_run.fully_buffered(command)

# ---------------------------------- PRIME ----------------------------------- #

tp_EVT_PRIMEDONE = wx.NewEventType()
EVT_PRIMEDONE = wx.PyEventBinder(tp_EVT_PRIMEDONE, 1)

class PRIMEAllDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class PRIMEThread(Thread):
  """ Thread for running all PRIME calculations; will prepare info for plotting
      in main GUI thread """
  def __init__(self, parent, info, params, best_pg=None, best_uc=None):
    Thread.__init__(self, name='live_prime')
    self.parent = parent
    self.params = params
    self.info = info
    self.best_pg = best_pg
    self.best_uc = best_uc

  def prepare_PRIME_input(self):
    """ Prepare the list of integrated pickles as well as pertinent
    parameters; create a PRIME input file """

    # Check if any pg/uc information is available; pg/uc info from UI
    # supercedes pg/uc from Cluster

    if self.best_pg is None:
      self.best_pg = self.info.best_pg

    if self.best_uc is None:
      self.best_uc = self.info.best_uc
    else:
      uc_params = str(self.best_uc).rsplit()
      self.best_uc = [float(i) for i in uc_params]

    # Only run PRIME if both pg and uc are provided
    if self.best_pg and self.best_uc:
      from iota.components.iota_analysis import Analyzer

      # Create a file with list of integrated pickles (overwrite old file)
      int_pickles_file = os.path.join(self.info.int_base, 'int_pickles.lst')
      with open(int_pickles_file, 'w') as ff:
        ff.write('\n'.join(self.info.categories['integrated'][0]))

      # Create PRIME input file
      analyzer = Analyzer(info=self.info, params=self.params)
      analyzer.prime_data_path = int_pickles_file
      analyzer.best_pg = self.best_pg
      analyzer.best_uc = self.best_uc

      prime_phil = analyzer.make_prime_input(filename='live_prime.phil',
                                             run_zero=True)
      self.pparams = prime_phil.extract()

      # Modify specific options based in IOTA settings
      # Queue options
      if (
              self.params.mp.method == 'lsf' and
              self.params.mp.queue is not None
      ):
        self.pparams.queue.mode = 'bsub'
        self.pparams.queue.qname = self.params.mp.queue

      # Number of processors (automatically, 1/2 of IOTA procs)
      self.pparams.n_processors = int(self.params.mp.n_processors / 2)

      # Generate command args
      cmd_args_list = ['n_postref_cycle=0',
                       'queue.mode={}'.format(self.pparams.queue.mode),
                       'queue.qname={}'.format(self.pparams.queue.qname),
                       'n_processors={}'.format(self.pparams.n_processors)]
      if self.pparams.queue.mode == 'bsub':
        cmd_args_list.append('timeout_seconds=120')
      cmd_args = ' '.join(cmd_args_list)

      return cmd_args
    else:
      return None

  def get_prime_stats(self):
    stats_folder = os.path.join(self.pparams.run_no, 'stats')
    prime_info = None
    if os.path.isdir(stats_folder):
      stat_files = [os.path.join(stats_folder, i) for i in
                    os.listdir(stats_folder) if i.endswith('stat')]
      if stat_files:
        assert len(stat_files) == 1
        stat_file = stat_files[0]
        if os.path.isfile(stat_file):
          prime_info = ep.load(stat_file)
          live_prime_info_file = os.path.join(self.info.int_base,
                                              'life_prime_info.pickle')
          shutil.copyfile(stat_file, live_prime_info_file)

    # Convert space_group_info object to str, to make compatible with JSON
    if prime_info and 'space_group_info' in prime_info:
      for i in prime_info['space_group_info']:
        sg = str(i).split('(')[0].rstrip(' ')
        idx = prime_info['space_group_info'].index(i)
        prime_info['space_group_info'][idx] = sg

    return prime_info

  def abort(self):
    # TODO: put in an LSF kill command
    if hasattr(self, 'job'):
      try:
        self.job.kill_thread()
      except Exception as e:
        print ('PRIME THREAD ERROR: Cannot terminate thread! {}'.format(e))

  def run(self):
    # Generate PRIME input
    cmd_args = self.prepare_PRIME_input()

    # Run PRIME
    if cmd_args:
      # remove previous run to avoid conflict
      prime_dir = os.path.join(self.info.int_base, 'prime/000')
      if os.path.isdir(prime_dir):
        shutil.rmtree(prime_dir)

      # Launch PRIME
      prime_file = os.path.join(self.info.int_base, 'live_prime.phil')
      cmd = 'prime.run {} {}'.format(prime_file, cmd_args)

      try:
        # easy_run.fully_buffered(cmd, join_stdout_stderr=True).show_stdout()
        self.job = CustomRun(command=cmd,
                             join_stdout_stderr=True)
        self.job.run()
        prime_info = self.get_prime_stats()
      except Exception as e:
        print ("LIVE PRIME ERROR: ", e)
        prime_info = None

    else:
      prime_info = None

    # Send signal to UI
    evt = PRIMEAllDone(tp_EVT_PRIMEDONE, -1, info=prime_info)
    wx.PostEvent(self.parent, evt)

#------------------------------ IMAGE TRACKING ------------------------------ #

tp_EVT_SPFDONE = wx.NewEventType()
EVT_SPFDONE = wx.PyEventBinder(tp_EVT_SPFDONE, 1)

tp_EVT_SPFALLDONE = wx.NewEventType()
EVT_SPFALLDONE = wx.PyEventBinder(tp_EVT_SPFALLDONE, 1)

tp_EVT_SLICEDONE = wx.NewEventType()
EVT_SLICEDONE = wx.PyEventBinder(tp_EVT_SLICEDONE)

tp_EVT_SPFTERM = wx.NewEventType()
EVT_SPFTERM = wx.PyEventBinder(tp_EVT_SPFTERM)

class SpotFinderAllDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderOneDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderTerminated(wx.PyCommandEvent):
  """ Send event when spotfinder terminated """
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)
  def GetValue(self):
    return None

class SpotFinderDIALSThread():
  def __init__(self, parent, processor, term_file,
               run_indexing=False, run_integration=False):
    self.meta_parent = parent.parent
    self.processor = processor
    self.term_file = term_file
    self.run_indexing = run_indexing
    self.run_integration = run_integration

  def run(self, idx, img):
    if self.meta_parent.terminated:
      raise IOTATermination('IOTA_TRACKER: SPF Termination signal received!')
    else:
      with Capturing() as junk_output:
        fail = False
        sg = None
        uc = None
        try:
          experiments = ExperimentListFactory.from_filenames([img])[0]
          observed = self.processor.find_spots(experiments=experiments)
        except Exception:
          fail = True
          observed = []
          pass

        # TODO: Indexing / lattice determination very slow (how to speed up?)
        if self.run_indexing:
          if not fail:
            try:
              experiments, indexed = self.processor.index(
                experiments=experiments, reflections=observed)
            except Exception:
              fail = True
              pass

          if not fail:
            try:
              solution = self.processor.refine_bravais_settings(
                reflections=indexed, experiments=experiments)

              # Only reindex if higher-symmetry solution found
              if solution is not None:
                experiments, indexed = self.processor.reindex(
                  reflections=indexed,
                  experiments=experiments,
                  solution=solution)
              lat = experiments[0].crystal.get_space_group().info()
              sg = str(lat).replace(' ', '')
            except Exception:
              fail = True
              pass

          if not fail:
            unit_cell = experiments[0].crystal.get_unit_cell().parameters()
            uc = ' '.join(['{:.4f}'.format(i) for i in unit_cell])

          if self.run_integration:
            if not fail:
              try:
                # Run refinement
                experiments, indexed = self.processor.refine(
                  experiments=experiments, centroids=indexed)
              except Exception:
                fail = True
                pass

            if not fail:
              try:
                print (experiments)
                print (indexed)
                integrated = self.processor.integrate(experiments=experiments,
                                                      indexed=indexed)
              except Exception:
                pass

      return [idx, int(len(observed)), img, sg, uc]

class SpotFinderMosflmThread():
  def __init__(self, parent, term_file):
    self.meta_parent = parent.parent
    self.term_file = term_file

  def run(self, idx, img):
    if os.path.isfile(self.term_file):
      raise IOTATermination('IOTA_TRACKER: Termination signal received!')
    else:
      # First, parse filepath to create Mosflm template
      directory = os.path.dirname(img)
      filepath = os.path.basename(img).split('.')
      fname = filepath[0]
      extension = filepath[1]
      if '_' in fname:
        suffix = fname.split('_')[-1]
      elif '-' in fname:
        suffix = fname.split('-')[-1]
      elif '.' in fname:
        suffix = fname.split('.')[-1]
      else:
        suffix = fname
      img_number = int(''.join(n if n.isdigit() else '' for n in suffix))
      prefix = fname.replace(suffix, '')
      n_suffix = ''.join("#" if c.isdigit() else c for c in suffix)
      template = '{}{}.{}'.format(prefix, n_suffix, extension)

      # Create autoindex.com w/ Mosflm script
      # Write to temporary file and change permissions to run
      autoindex = ['#! /bin/tcsh -fe',
                   'ipmosflm << eof-ipmosflm'.format(fname),
                   'NEWMATRIX {0}.mat'.format(fname),
                   'DIRECTORY {}'.format(directory),
                   'TEMPLATE {}'.format(template),
                   'AUTOINDEX DPS THRESH 0.1 IMAGE {} PHI 0 0.01'.format(
                     img_number),
                   'GO',
                   'eof-ipmosflm'
                   ]
      autoindex_string = '\n'.join(autoindex)
      autoindex_filename = 'autoindex_{}.com'.format(idx)

      with open(autoindex_filename, 'w') as af:
        af.write(autoindex_string)
      os.chmod(autoindex_filename, 0o755)

      # Run Mosflm autoindexing
      command = './{}'.format(autoindex_filename)
      out = easy_run.fully_buffered(command, join_stdout_stderr=True)

      # Scrub text output
      final_spots = [l for l in out.stdout_lines if
                     'spots written for image' in l]
      final_cell_line = [l for l in out.stdout_lines if 'Final cell' in l]
      final_sg_line = [l for l in out.stdout_lines if 'space group' in l]

      if final_spots:
        spots = final_spots[0].rsplit()[0]
      else:
        spots = 0
      if final_cell_line:
        cell = final_cell_line[0].replace('Final cell (after refinement) is',
                                          '').rsplit()
      else:
        cell = None
      if final_sg_line:
        sg = final_sg_line[0].rsplit()[6]
      else:
        sg = None

      # Temp file cleanup
      try:
        os.remove('{}.mat'.format(fname))
      except Exception:
        pass
      try:
        os.remove('{}.spt'.format(prefix[:-1]))
      except Exception:
        pass
      try:
        os.remove('SUMMARY')
      except Exception:
        pass
      try:
        os.remove(autoindex_filename)
      except Exception:
        pass

      return [idx, spots, img, sg, cell]

class SpotFinderThread(Thread):
  """ Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected """
  def __init__(self,
               parent,
               data_list=None,
               term_file=None,
               proc_params=None,
               backend='dials',
               n_proc=0,
               run_indexing=False,
               run_integration=False):
    Thread.__init__(self)
    self.parent = parent
    self.data_list = data_list
    self.term_file = term_file
    self.terminated = False
    self.backend = backend
    self.run_indexing = run_indexing
    self.run_integration = run_integration
    if n_proc > 0:
      self.n_proc = n_proc
    else:
      self.n_proc = multiprocessing.cpu_count() - 2

    if self.backend == 'dials':
      # Modify default DIALS parameters
      # These parameters will be set no matter what
      proc_params.output.experiments_filename = None
      proc_params.output.indexed_filename = None
      proc_params.output.strong_filename = None
      proc_params.output.refined_experiments_filename = None
      proc_params.output.integrated_filename = None
      proc_params.output.integrated_experiments_filename = None
      proc_params.output.profile_filename = None
      proc_params.output.integration_pickle = None

      from iota.components.iota_processing import IOTAImageProcessor
      self.processor = IOTAImageProcessor(phil=proc_params)

  def run(self):
    try:
      parallel_map(iterable=self.data_list,
                   func=self.spf_wrapper,
                   callback=self.callback,
                   processes=self.n_proc)
    except IOTATermination as e:
      self.terminated = True
      print (e)

    # Signal that this batch is finished
    try:
      if self.terminated:
        print ('RUN TERMINATED!')
        evt = SpotFinderTerminated(tp_EVT_SPFTERM, -1)
        wx.PostEvent(self.parent, evt)

      wx.CallAfter(self.parent.onSpfAllDone, self.data_list)

      # info = self.data_list
      # evt = SpotFinderAllDone(tp_EVT_SPFALLDONE, -1, info=info)
      # wx.PostEvent(self.parent, evt)
      return
    except TypeError as e:
      print (e)
      return

  def spf_wrapper(self, img):
    try:
      if os.path.isfile(img):
        if self.backend == 'dials':
          spf_worker = SpotFinderDIALSThread(self,
                                             processor=self.processor,
                                             term_file=self.term_file,
                                             run_indexing=self.run_indexing,
                                             run_integration=self.run_integration
                                             )
          result = spf_worker.run(idx=int(self.data_list.index(img)), img=img)
        elif self.backend == 'mosflm':
          spf_worker = SpotFinderMosflmThread(self, self.term_file)
          result = spf_worker.run(idx=int(self.data_list.index(img)), img=img)
        else:
          result = [int(self.data_list.index(img)), 0, img, None, None]
        return result
      else:
        return [int(self.data_list.index(img)), 0, img, None, None]
    except IOTATermination as e:
      raise e

  def callback(self, info):
    try:
      wx.CallAfter(self.parent.onSpfOneDone, info)
      # evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
      # wx.PostEvent(self.parent.parent, evt)
    except TypeError:
      pass

  # def terminate_thread(self):
  #   raise IOTATermination('IOTA_TRACKER: SPF THREAD Terminated!')

class InterceptorFileThread(Thread):
  def __init__(self,
               parent,
               results_file,
               reorder=False):
    Thread.__init__(self)
    self.parent = parent
    self.results_file = results_file
    self.reorder = reorder

    self.bookmark = 0
    self.msg = ''
    self.spotfinding_info = []
    self.cluster_info = None

    self.prc_timer = wx.Timer()
    self.cls_timer = wx.Timer()

    #Bindings
    self.prc_timer.Bind(wx.EVT_TIMER, self.onProcTimer)
    self.cls_timer.Bind(wx.EVT_TIMER, self.onClusterTimer)


  def run(self):
    # self.timer.Start(1000)
    pass

  def onProcTimer(self, e):
    if os.path.isfile(self.results_file):
      with open(self.results_file, 'r') as rf:
        rf.seek(self.bookmark)
        split_info = [i.replace('\n', '').split(' ') for i in rf.readlines()]
        self.bookmark = rf.tell()

      if self.reorder:
        idx_offset = len(self.spotfinding_info)
        new_info = [
          [split_info.index(i) + idx_offset,
           int(i[1]), i[2], i[3], tuple(i[4:10])] if len(i) > 5 else
          [split_info.index(i) + idx_offset,
           int(i[1]), i[2], util.makenone(i[3]), util.makenone(i[4])]
          for i in split_info]
      else:
        new_info = [
          [int(i[0]), int(i[1]), i[2], i[3], tuple(i[4:10])] if len(i) > 5 else
          [int(i[0]), int(i[1]), i[2], util.makenone(i[3]), util.makenone(i[4])]
          for i in split_info]

      if len(new_info) > 0:
        self.spotfinding_info.extend(new_info)

        if len(self.spotfinding_info) > 0:
          self.msg = 'Tracking new images in {} ...'.format(self.results_file)
      else:
        self.msg = 'Waiting for new images in {} ...'.format(self.results_file)

    else:
      self.msg = 'Waiting for new run to initiate...'

    info = [self.msg, self.spotfinding_info, self.cluster_info]
    evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)

  def onClusterTimer(self, e):
    input = []
    for item in self.spotfinding_info:
      if item[4] is not None:
        try:
          if type(item[4]) in (tuple, list):
            uc = item[4]
          else:
            uc = item[4].rsplit()
          info_line = [float(i) for i in uc]
          info_line.append(item[3])
          input.append(info_line)
        except ValueError as e:
          print ('CLUSTER ERROR: ', e)
          pass

    if len(input) > 0:
      self.running_clustering = True
      cluster_thread = ClusterWorkThread(self)
      self.cluster_info = cluster_thread.run(iterable=input)

  def terminate_thread(self):
    raise IOTATermination('IOTA_TRACKER: Termination signal received!')

class InterceptorThread(Thread):
  """ Thread for the full Interceptor image processing process; will also
   house the processing timer, which will update the UI front end and initiate
   plotting """
  def __init__(self,
               parent,
               data_folder=None,
               term_file=None,
               proc_params=None,
               backend='dials',
               n_proc=0,
               min_back=None,
               run_indexing=False,
               run_integration=False):
    Thread.__init__(self)
    self.parent = parent
    self.data_folder = data_folder
    self.term_file = term_file
    self.terminated = False
    self.backend = backend
    self.run_indexing = run_indexing
    self.run_integration = run_integration
    self.min_back = min_back
    self.submit_new_images = True

    self.spotfinding_info = []
    self.cluster_info = None
    self.msg = None
    self.done_list = []
    self.data_list = []
    self.new_data = []

    self.spf_thread = None

    if n_proc > 0:
      self.n_proc = n_proc
    else:
      self.n_proc = multiprocessing.cpu_count() - 2

    if self.backend == 'dials':
      # Modify default DIALS parameters
      # These parameters will be set no matter what
      proc_params.output.experiments_filename = None
      proc_params.output.indexed_filename = None
      proc_params.output.strong_filename = None
      proc_params.output.refined_experiments_filename = None
      proc_params.output.integrated_filename = None
      proc_params.output.integrated_experiments_filename = None
      proc_params.output.profile_filename = None
      proc_params.output.integration_pickle = None

      self.proc_params = proc_params

      # from iota.components.iota_dials import IOTADialsProcessor
      # self.processor = IOTADialsProcessor(params=proc_params)

    self.prc_timer = wx.Timer()
    self.cls_timer = wx.Timer()

    #Bindings
    self.prc_timer.Bind(wx.EVT_TIMER, self.onProcTimer)
    self.cls_timer.Bind(wx.EVT_TIMER, self.onClusterTimer)


  def run(self):
    pass

  def onProcTimer(self, e):
    """ Main timer (1 sec) will send data to UI, find new images, and submit
    new processing run """

    # Send current data to UI
    info = [self.msg, self.spotfinding_info, self.cluster_info]
    evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)

    # Find new images
    if self.data_list:
      last_file = self.data_list[-1]
    else:
      last_file = None
    self.find_new_images(last_file=last_file, min_back=self.min_back)

    if self.spf_thread is not None:
      if not self.spf_thread.isAlive():
        self.submit_new_images = True

    # Submit new images (if found)
    if self.submit_new_images and len(self.new_data) > 0:
      self.submit_new_images = False
      self.run_processing()

  def onClusterTimer(self, e):
    input = []
    if len(self.spotfinding_info) > 0:
      for item in self.spotfinding_info:
        if item[4] is not None:
          try:
            if type(item[4]) in (tuple, list):
              uc = item[4]
            else:
              uc = item[4].rsplit()
            info_line = [float(i) for i in uc]
            info_line.append(item[3])
            input.append(info_line)
          except ValueError as e:
            print ('CLUSTER ERROR: ', e)
            pass

      if len(input) > 0:
        self.running_clustering = True
        cluster_thread = ClusterWorkThread(self)
        self.cluster_info = cluster_thread.run(iterable=input)

  def find_new_images(self, min_back=None, last_file=None):
    found_files = ginp.make_input_list([self.data_folder],
                                       filter_results=True,
                                       filter_type='image',
                                       last=last_file,
                                       min_back=min_back)

    # Sometimes duplicate files are found anyway; clean that up
    found_files = list(set(found_files) - set(self.data_list))

    # Add new files to the data list & clean up
    self.new_data.extend(found_files)
    self.new_data = sorted(self.new_data, key=lambda i:i)
    self.data_list.extend(self.new_data)

  def run_processing(self):
    self.spf_thread = SpotFinderThread(self,
                                       data_list=self.new_data,
                                       term_file=self.term_file,
                                       proc_params=self.proc_params,
                                       backend=self.backend,
                                       n_proc=self.n_proc,
                                       run_indexing=self.run_indexing,
                                       run_integration=self.run_integration)
    self.new_data = []
    self.spf_thread.start()


  def onSpfOneDone(self, info):
    info[0] = int(info[0]) + len(self.done_list)
    self.spotfinding_info.append(info)

  def onSpfAllDone(self, done_list):
    self.done_list.extend(done_list)

  def terminate_thread(self):
    self.terminated = True

# ------------------------------ UC CLUSTERING ------------------------------- #

tp_EVT_CLUSTERDONE = wx.NewEventType()
EVT_CLUSTERDONE = wx.PyEventBinder(tp_EVT_CLUSTERDONE, 1)

class ClusteringDone(wx.PyCommandEvent):
  """ Send event when finished all cycles  """
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class ClusterWorkThread():
  def __init__(self, parent):
    self.parent = parent
    self.abort = False

  def run(self, iterable):

    # with Capturing() as junk_output:
    errors = []
    try:
      ucs = Cluster.from_iterable(iterable=iterable)
      clusters, _ = ucs.ab_cluster(5000, log=False, write_file_lists=False,
                                   schnell=True, doplot=False)
    except Exception as e:
      print ('IOTA ERROR (CLUSTERING): ', e)
      clusters = []
      errors.append(e)

    info = []
    if clusters:
      for cluster in clusters:
        uc_init = unit_cell(cluster.medians)
        symmetry = crystal.symmetry(unit_cell=uc_init, space_group_symbol='P1')
        groups = lattice_symmetry.metric_subgroups(input_symmetry=symmetry,
                                                   max_delta=3)
        top_group = groups.result_groups[0]
        best_sg = str(groups.lattice_group_info()).split('(')[0]
        best_uc = top_group['best_subsym'].unit_cell().parameters()
        uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "".format(best_uc[0], best_uc[1], best_uc[2],
                                best_uc[3], best_uc[4], best_uc[5])
        cluster_info = {'number': len(cluster.members),
                        'pg': str(best_sg),
                        'uc': uc_no_stdev}
        info.append(cluster_info)

    return info, errors


class ClusterThread(Thread):
  """ Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected """
  def __init__(self,
               parent,
               iterable):
    Thread.__init__(self, name='live_cluster')
    self.parent = parent
    self.iterable = iterable
    self.clustering = ClusterWorkThread(self)

  def abort(self):
    self.clustering.abort = True

  def run(self):
    clusters, errors = self.clustering.run(iterable=self.iterable)

    if clusters:
      clusters = sorted(clusters, key=lambda i: i['number'], reverse=True)
    evt = SpotFinderOneDone(tp_EVT_CLUSTERDONE, -1, info=[clusters, errors])
    wx.PostEvent(self.parent, evt)


# ------------------------------- SUBPROCESS --------------------------------- #

class CustomRun(easy_run.fully_buffered_base):
  ''' A subclass from easy_run with a "kill switch" for easy process
      termination from UI; took out timeout, since that won't be used,
      and doesn't work on all systems, anyway

      Tested on Mac OS X 10.13.6 and CentOS 6 so far.
  '''
  def __init__(self,
        command,
        timeout=None,
        stdin_lines=None,
        join_stdout_stderr=False,
        stdout_splitlines=True,
        bufsize=-1):
    self.command = command
    self.join_stdout_stderr = join_stdout_stderr
    self.timeout = timeout
    self.stdin_lines = stdin_lines
    self.stdout_splitlines=stdout_splitlines
    self.bufsize = bufsize
    self.thread = None

  def target(self, process, lines, result):
    o, e = process.communicate(input=lines)
    result[0] = o
    result[1] = e

  def kill_thread(self):
    os.killpg(os.getpgid(self.p.pid), signal.SIGTERM)

  def run(self):
    if (not isinstance(self.command, str)):
      self.command = subprocess.list2cmdline(self.command)
    if (sys.platform == 'darwin'):   # bypass SIP on OS X 10.11
      self.command = ("DYLD_LIBRARY_PATH=%s exec "%
                      os.environ.get("DYLD_LIBRARY_PATH","")) + self.command
    if (self.stdin_lines is not None):
      if (not isinstance(self.stdin_lines, str)):
        self.stdin_lines = os.linesep.join(self.stdin_lines)
        if (len(self.stdin_lines) != 0):
          self.stdin_lines += os.linesep
    if (self.join_stdout_stderr):
      stderr = subprocess.STDOUT
    else:
      stderr = subprocess.PIPE

    self.p = subprocess.Popen(
      args=self.command,
      shell=True,
      bufsize=self.bufsize,
      stdin=subprocess.PIPE,
      stdout=subprocess.PIPE,
      stderr=stderr,
      universal_newlines=True,
      close_fds=(sys.platform != 'win32'),
      preexec_fn=os.setsid if sys.platform != 'win32' else None)
    o, e = self.p.communicate(input=self.stdin_lines)
    if (self.stdout_splitlines):
      self.stdout_buffer = None
      self.stdout_lines = o.splitlines()
    else:
      self.stdout_buffer = o
      self.stdout_lines = None
    if (self.join_stdout_stderr):
      self.stderr_lines = []
    else:
      self.stderr_lines = e.splitlines()
    self.return_code = self.p.returncode
