from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 02/12/2018
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

from iota.components.iota_utils import InputFinder
import iota.components.iota_image as img

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
    # if self.init.params.mp_method == 'multiprocessing':
    try:
      img_objects = parallel_map(iterable=self.iterable,
                                 func = self.full_proc_wrapper,
                                 processes=self.init.params.n_processors)
    except IOTATermination, e:
      self.aborted = True
      print e
      return
    # else:
    #   # write iterable
    #   img_objects = None
    #   queue = self.init.params.mp_queue
    #   iter_path = os.path.join(self.init.int_base, 'iter.cfg')
    #   init_path = os.path.join(self.init.int_base, 'init.cfg')
    #   nproc = self.init.params.n_processors
    #   ep.dump(iter_path, self.iterable)
    #   ep.dump(init_path, self.init)
    #
    #   if self.init.params.mp_method == 'lsf':
    #     logfile = os.path.join(self.init.int_base, 'bsub.log')
    #     command = 'bsub -o {} -q {} -n {} ' \
    #               'iota.process {} --files {} --type {} --stopfile {}' \
    #               ''.format(logfile, queue, nproc,
    #                         init_path, iter_path, self.type, self.term_file)
    #   elif self.init.params.mp_method == 'torq':
    #     params = '{} --files {} --type {} --stopfile {}' \
    #              ''.format(init_path, iter_path, self.type, self.term_file)
    #     command = 'qsub -e /dev/null -o /dev/null -d {} iota.process -F "{}"' \
    #               ''.format(self.init.params.output, params)
    #   else:
    #     command = None
    #   if command is not None:
    #     try:
    #       print command
    #       easy_run.fully_buffered(command, join_stdout_stderr=True).show_stdout()
    #     except IOTATermination, e:
    #       print e
    #   else:
    #     print 'IOTA ERROR: COMMAND NOT ISSUED!'
    #     return

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

    ext_file_list = ginp.make_input_list(self.image_paths)
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
               fix_paths = False,
               new_fin_base = None):
    Thread.__init__(self)
    self.parent = parent
    self.object_folder = object_folder
    self.fix_paths = fix_paths
    self.new_fin_base = new_fin_base

  def run(self):
    object_files = ginp.get_file_list(self.object_folder, ext_only='int')
    new_objects = [self.read_object_file(i) for i in object_files]
    new_finished_objects = [i for i in new_objects if
                            (i is not None and i.status == 'final')]

    # # If recovering and need to fix paths of images (final only):
    # if self.fix_paths:
    #   for obj in new_finished_objects:
    #     print 'Changing ', obj.final['final']
    #     if obj.final['final'] is not None:
    #       filename = os.path.basename(obj.final['final'])
    #       obj.final['final'] = os.path.join(self.new_fin_base, filename)
    #       print 'Generating filename ', obj.final['final']

    evt = ObjectFinderAllDone(tp_EVT_OBJDONE, -1, obj_list=new_finished_objects)
    wx.PostEvent(self.parent, evt)

  def read_object_file(self, filepath):
    try:
      object = ep.load(filepath)
      return object
    except EOFError:
      pass

class ImageViewerThread(Thread):
  ''' Worker thread that will move the image viewer launch away from the GUI
  and hopefully will prevent the image selection dialog freezing on MacOS'''
  def __init__(self,
               parent,
               file_string,
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


# ------------------------------ IMAGE TRACKING ------------------------------ #

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

class SpotFinderOneThread():
  def __init__(self, parent, processor, term_file):
    self.meta_parent = parent.parent
    self.processor = processor
    self.term_file = term_file

  def run(self, idx, img):
    if os.path.isfile(self.term_file):
      raise IOTATermination('IOTA_TRACKER: Termination signal received!')
    else:
      datablock = DataBlockFactory.from_filenames([img])[0]
      observed = self.processor.find_spots(datablock=datablock)
      return [idx, int(len(observed)), img]

class SpotFinderThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
  def __init__(self,
               parent,
               data_list=None,
               term_file=None,
               processor=None):
    Thread.__init__(self)
    self.parent = parent
    self.data_list = data_list
    self.term_file = term_file
    self.processor = processor
    self.terminated = False

  def run(self):
    total_procs = multiprocessing.cpu_count()
    try:
      parallel_map(iterable=self.data_list,
                   func=self.spf_wrapper,
                   callback=self.callback,
                   processes=total_procs-5)
    except IOTATermination, e:
      self.terminated = True
      print e

    # Signal that this batch is finished
    try:
      if self.terminated:
        print 'RUN TERMINATED!'
        evt = SpotFinderTerminated(tp_EVT_SPFTERM, -1)
        wx.PostEvent(self.parent, evt)
      info = self.data_list
      evt = SpotFinderOneDone(tp_EVT_SPFALLDONE, -1, info=info)
      wx.PostEvent(self.parent, evt)
      return
    except TypeError:
      return

  def spf_wrapper(self, img):
    # It appears that having a separate thread for each process prevents
    # blocking of the GUI and makes navigation faster
    try:
      if os.path.isfile(img):
        spf_worker = SpotFinderOneThread(self, self.processor, self.term_file)
        result = spf_worker.run(idx=int(self.data_list.index(img)), img=img)
        return result
      else:
        return [int(self.data_list.index(img)), 0, img]
    except IOTATermination, e:
      raise e


  # def spf_wrapper(self, img):
  #   if os.path.isfile(self.term_file):
  #     os.remove(self.term_file)
  #     raise IOTATermination('IOTA_TRACKER: Termination signal received!')
  #   else:
  #     if os.path.isfile(img):
  #       datablock = DataBlockFactory.from_filenames([img])[0]
  #       observed = self.processor.find_spots(datablock=datablock)
  #       return [int(self.data_list.index(img)), int(len(observed)), img]
  #     else:
  #       return [int(self.data_list.index(img)), 0, img]

  def callback(self, info):
    try:
      evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
      wx.PostEvent(self.parent, evt)
    except TypeError:
      pass
