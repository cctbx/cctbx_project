from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 08/31/2018
Description : IOTA GUI Threads and PostEvents
'''

import os
import wx
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep
from libtbx import easy_run

from dxtbx.datablock import DataBlockFactory
import multiprocessing

from xfel.clustering.cluster import Cluster
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import lattice_symmetry
from cctbx import crystal

from iota.components.iota_utils import InputFinder
from iota.components.iota_misc import Capturing
import iota.components.iota_image as img
import iota.components.iota_misc as misc

ginp = InputFinder()

# -------------------------------- Threading --------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

tp_EVT_IMGDONE = wx.NewEventType()
EVT_IMGDONE = wx.PyEventBinder(tp_EVT_IMGDONE, 1)

tp_EVT_OBJDONE = wx.NewEventType()
EVT_OBJDONE = wx.PyEventBinder(tp_EVT_OBJDONE, 1)

class ImageFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, image_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_list = image_list
  def GetValue(self):
    return self.image_list

class ObjectFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, obj_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.obj_list = obj_list
  def GetValue(self):
    return self.obj_list

class AllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, img_objects=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_objects = img_objects
  def GetValue(self):
    return self.image_objects

class ProcessImage():
  ''' Wrapper class to do full processing of an image '''
  def __init__(self, init, input_entry, input_type = 'image', abort=False):
    self.init = init
    self.input_entry = input_entry
    self.input_type = input_type
    self.abort = abort

  def run(self):
    if self.abort:
      raise IOTATermination('IOTA: Run aborted by user')
    else:
      if self.input_type == 'image':
        img_object = img.SingleImage(self.input_entry, self.init)
        img_object.import_image()
      elif self.input_type == 'object':
        img_object = self.input_entry[2]
        img_object.import_int_file(self.init)
      else:
        img_object = None

      if self.init.params.image_conversion.convert_only:
        return img_object
      else:
        img_object.process()
        return img_object

class ProcThread(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''
  def __init__(self,
               parent,
               init,
               iterable,
               term_file,
               input_type='image'):
    Thread.__init__(self)
    self.parent = parent
    self.init = init
    self.iterable = iterable
    self.type = input_type
    self.term_file = term_file
    self.aborted = False

  def run(self):
    try:
      img_objects = parallel_map(iterable=self.iterable,
                                 func = self.full_proc_wrapper,
                                 processes=self.init.params.n_processors)
    except IOTATermination, e:
      self.aborted = True
      print e
      return

    # Send "all done" event to GUI
    try:
      evt = AllDone(tp_EVT_ALLDONE, -1, img_objects=img_objects)
      wx.PostEvent(self.parent, evt)
    except Exception, e:
      pass

  def full_proc_wrapper(self, input_entry):
    abort = os.path.isfile(self.term_file)
    if abort:
      os.remove(self.term_file)
    try:
      proc_image_instance = ProcessImage(init=self.init,
                                         input_entry=input_entry,
                                         input_type=self.type,
                                         abort=abort)
      proc_image = proc_image_instance.run()
      return proc_image
    except IOTATermination, e:
      raise e
    except Exception, e:
      pass

class ImageFinderThread(Thread):
  ''' Worker thread generated to poll filesystem on timer. Will check to see
  if any new images have been found. Put on a thread to run in background '''
  def __init__(self,
               parent,
               image_paths,
               image_list):
    Thread.__init__(self)
    self.parent = parent
    self.image_paths = image_paths
    self.image_list = image_list

  def run(self):
    # Poll filesystem and determine which files are new (if any)

    ext_file_list = ginp.make_input_list(self.image_paths,
                                         filter=True,
                                         filter_type='image')
    old_file_list = [i[2] for i in self.image_list]
    new_file_list = [i for i in ext_file_list if i not in old_file_list]

    # Generate list of new images
    new_img = [[i, len(ext_file_list) + 1, j] for i, j in enumerate(
      new_file_list, len(old_file_list) + 1)]

    evt = ImageFinderAllDone(tp_EVT_IMGDONE, -1, image_list=new_img)
    wx.PostEvent(self.parent, evt)

class ObjectFinderThread(Thread):
  ''' Worker thread that polls filesystem on timer for image objects. Will
  collect and extract info on images processed so far'''
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
    except EOFError, e:
      print 'OBJECT_IMPORT_ERROR: ', e
      return None

class ImageViewerThread(Thread):
  ''' Worker thread that will move the image viewer launch away from the GUI
  and hopefully will prevent the image selection dialog freezing on MacOS'''
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

# -------------------------------- UI Thread --------------------------------- #

# Set up events for finishing one timer cycle

tp_EVT_PROCTIMER = wx.NewEventType()
EVT_PROCTIMER = wx.PyEventBinder(tp_EVT_PROCTIMER, 1)

class ProcTimerDone(wx.PyCommandEvent):
  ''' Send event at every ProcTimer ping  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class ProcessingInfo(object):

  def __init__(self):
    ''' constructor '''

    self.image_objects = None
    self.nref_list = None
    self.res_list = None
    self.img_list = None
    self.indices = None
    self.b_factors = None
    self.cluster_info = None
    self.prime_info = None

class IOTAUIThread(Thread):
  ''' Main thread for IOTA UI; will contain all times and call all the other
  threads - processing, object finding, etc. - separately; will use
  PostEvents to send data to the main UI thread, which will plot only. The
  idea is to prevent UI blocking as much as possible '''

  def __init__(self,
               parent,
               gparams,
               target_phil,
               tmp_aborted_file=None,
               proc_info_file=None,
               recover=False):
    Thread.__init__(self)
    self.parent = parent
    self.gparams = gparams
    self.target_phil = target_phil
    self.tmp_aborted_file = tmp_aborted_file
    self.recover = recover

    # Instantiate info object
    if proc_info_file is None:
      self.info = ProcessingInfo()
    else:
      self.info = ep.load(proc_info_file)

    # Timers
    self.plt_timer = wx.Timer()
    self.anl_timer = wx.Timer()

    # Bindings
    self.plt_timer.Bind(wx.EVT_TIMER, self.onPlotTimer)
    self.anl_timer.Bind(wx.EVT_TIMER, self.onAnalysisTimer)

  def run(self):
    if self.recover:
      self.recover_info()
    else:
      pass

  def onPlotTimer(self, e):
    ''' One second timer for status check and plotting '''

    # Send info to UI
    info = []
    evt = SpotFinderOneDone(tp_EVT_PROCTIMER, -1, info=info)
    wx.PostEvent(self.parent, evt)

  def onAnalysisTimer(self, e):
    pass

  def recover_info(self):
    ''' Recover information from previous run (is here only for
    backwards compatibility reasons; normally all info will come from file) '''
    pass

  def run_clustering_thread(self):
    # Run clustering
    pass

  def onFinishedCluster(self, e):
    pass

  def run_prime_thread(self):
    # Run PRIME (basic merge only)
    pass

  def onFinishedPRIME(self, e):
    pass

  def get_prime_stats(self):
    pass

  def process_images(self):
    ''' One-fell-swoop importing / triaging / integration of images '''
    pass

  def analyze_results(self, analysis=None):
    pass


  def find_objects(self, find_old=False):
    pass

  def read_object_file(self, filepath):
    pass

  def populate_data_points(self, objects=None):
    pass


#------------------------------ IMAGE TRACKING ------------------------------ #

tp_EVT_SPFDONE = wx.NewEventType()
EVT_SPFDONE = wx.PyEventBinder(tp_EVT_SPFDONE, 1)

tp_EVT_SPFALLDONE = wx.NewEventType()
EVT_SPFALLDONE = wx.PyEventBinder(tp_EVT_SPFALLDONE, 1)

tp_EVT_SLICEDONE = wx.NewEventType()
EVT_SLICEDONE = wx.PyEventBinder(tp_EVT_SLICEDONE)

tp_EVT_SPFTERM = wx.NewEventType()
EVT_SPFTERM = wx.PyEventBinder(tp_EVT_SPFTERM)

class IOTATermination(Exception):
  def __init__(self, termination):
    Exception.__init__(self, termination)

class SpotFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderOneDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderTerminated(wx.PyCommandEvent):
  ''' Send event when spotfinder terminated '''
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
          datablock = DataBlockFactory.from_filenames([img])[0]
          observed = self.processor.find_spots(datablock=datablock)
        except Exception, e:
          fail = True
          observed = []
          pass

        # TODO: Indexing / lattice determination very slow (how to speed up?)
        if self.run_indexing:
          if not fail:
            try:
              experiments, indexed = self.processor.index(
                datablock=datablock, reflections=observed)
            except Exception, e:
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
                  experiments=experiments,
                  centroids=indexed)
              except Exception, e:
                fail = True
                pass

            if not fail:
              try:
                print experiments
                print indexed
                integrated = self.processor.integrate(experiments=experiments,
                                                      indexed=indexed)
              except Exception, e:
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
      os.chmod(autoindex_filename, 0755)

      # Run Mosflm autoindexing
      command = './{}'.format(autoindex_filename)
      out = easy_run.fully_buffered(command, join_stdout_stderr=True)

      # Scrub text output
      final_spots = [l for l in out.stdout_lines if
                     'spots written for image' in l]
      final_cell_line = [l for l in out.stdout_lines if 'Final cell' in l]
      final_sg_line = [l for l in out.stdout_lines if 'space group' in l]

      if final_spots != []:
        spots = final_spots[0].rsplit()[0]
      else:
        spots = 0
      if final_cell_line != []:
        cell = final_cell_line[0].replace('Final cell (after refinement) is',
                                          '').rsplit()
      else:
        cell = None
      if final_sg_line != []:
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
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
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
      proc_params.output.datablock_filename = None
      proc_params.output.indexed_filename = None
      proc_params.output.strong_filename = None
      proc_params.output.refined_experiments_filename = None
      proc_params.output.integrated_filename = None
      proc_params.output.integrated_experiments_filename = None
      proc_params.output.profile_filename = None
      proc_params.output.integration_pickle = None

      from iota.components.iota_dials import IOTADialsProcessor
      self.processor = IOTADialsProcessor(params=proc_params)

  def run(self):
    try:
      parallel_map(iterable=self.data_list,
                   func=self.spf_wrapper,
                   callback=self.callback,
                   processes=self.n_proc)
    except IOTATermination, e:
      self.terminated = True
      print e

    # Signal that this batch is finished
    try:
      if self.terminated:
        print 'RUN TERMINATED!'
        evt = SpotFinderTerminated(tp_EVT_SPFTERM, -1)
        wx.PostEvent(self.parent, evt)

      wx.CallAfter(self.parent.onSpfAllDone, self.data_list)

      # info = self.data_list
      # evt = SpotFinderAllDone(tp_EVT_SPFALLDONE, -1, info=info)
      # wx.PostEvent(self.parent, evt)
      return
    except TypeError, e:
      print e
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
    except IOTATermination, e:
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
           int(i[1]), i[2], misc.makenone(i[3]), misc.makenone(i[4])]
          for i in split_info]
      else:
        new_info = [
          [int(i[0]), int(i[1]), i[2], i[3], tuple(i[4:10])] if len(i) > 5 else
          [int(i[0]), int(i[1]), i[2], misc.makenone(i[3]), misc.makenone(i[4])]
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
        except ValueError, e:
          print 'CLUSTER ERROR: ', e
          pass

    if len(input) > 0:
      self.running_clustering = True
      cluster_thread = ClusterWorkThread(self)
      self.cluster_info = cluster_thread.run(iterable=input)

  def terminate_thread(self):
    raise IOTATermination('IOTA_TRACKER: Termination signal received!')

class InterceptorThread(Thread):
  ''' Thread for the full Interceptor image processing process; will also
   house the processing timer, which will update the UI front end and initiate
   plotting '''
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
      proc_params.output.datablock_filename = None
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
    ''' Main timer (1 sec) will send data to UI, find new images, and submit
    new processing run '''

    # Send current data to UI
    info = [self.msg, self.spotfinding_info, self.cluster_info]
    evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)

    # Find new images
    if self.data_list != []:
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
          except ValueError, e:
            print 'CLUSTER ERROR: ', e
            pass

      if len(input) > 0:
        self.running_clustering = True
        cluster_thread = ClusterWorkThread(self)
        self.cluster_info = cluster_thread.run(iterable=input)

  def find_new_images(self, min_back=None, last_file=None):
    found_files = ginp.make_input_list([self.data_folder],
                                       filter=True,
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
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class ClusterWorkThread():
  def __init__(self, parent):
    self.parent = parent

  def run(self, iterable):

    with Capturing() as junk_output:
      try:
        ucs = Cluster.from_iterable(iterable=iterable)
        clusters, _ = ucs.ab_cluster(5000,
                                     log=False, write_file_lists=False,
                                     schnell=True, doplot=False)
      except Exception:
        clusters = []

    if len(clusters) > 0:
      info = []
      for cluster in clusters:
        uc_init = unit_cell(cluster.medians)
        symmetry = crystal.symmetry(unit_cell=uc_init, space_group_symbol='P1')
        groups = lattice_symmetry.metric_subgroups(input_symmetry=symmetry,
                                                   max_delta=3)
        top_group = groups.result_groups[0]
        best_uc = top_group['best_subsym'].unit_cell().parameters()
        best_sg = top_group['best_subsym'].space_group_info()

        uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "".format(best_uc[0], best_uc[1], best_uc[2],
                                best_uc[3], best_uc[4], best_uc[5])
        cluster_info = {'number': len(cluster.members),
                        'pg': str(best_sg),
                        'uc': uc_no_stdev}

        info.append(cluster_info)

    else:
      info = None

    return info

class ClusterThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
  def __init__(self,
               parent,
               iterable):
    Thread.__init__(self)
    self.parent = parent
    self.iterable = iterable
    self.clustering = ClusterWorkThread(self)

  def run(self):
    info = self.clustering.run(iterable=self.iterable)
    evt = SpotFinderOneDone(tp_EVT_CLUSTERDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)
