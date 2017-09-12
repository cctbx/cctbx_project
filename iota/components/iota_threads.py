from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 09/12/2017
Description : IOTA GUI Threads and PostEvents
'''

import os
import wx
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep
from libtbx import easy_run

from dxtbx.datablock import DataBlockFactory

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
      raise Exception('IOTA: Run aborted by user')
    else:
      if self.input_type == 'image':
        img_object = img.SingleImage(self.input_entry, self.init)
        img_object.import_image()
      elif self.input_type == 'object':
        img_object = self.input_entry[2]
        img_object.import_int_file(self.init)

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

  def run(self):
    if self.init.params.mp_method == 'multiprocessing':
      try:
        img_objects = parallel_map(iterable=self.iterable,
                                   func = self.full_proc_wrapper,
                                   processes=self.init.params.n_processors)
      except Exception, e:
        print e
        return
    else:
      # write iterable
      img_objects = None
      queue = self.init.params.mp_queue
      iter_path = os.path.join(self.init.int_base, 'iter.cfg')
      init_path = os.path.join(self.init.int_base, 'init.cfg')
      nproc = self.init.params.n_processors
      ep.dump(iter_path, self.iterable)
      ep.dump(init_path, self.init)
      try:
        if self.init.params.mp_method == 'lsf':
          command = 'bsub -q {} -n {} iota.process {} --files {} --type {} ' \
                    '--stopfile {}'.format(queue, nproc, init_path, iter_path,
                                             self.type, self.term_file)
        if self.init.params.mp_method == 'torq':
          params = '{} --files {} --type {} --stopfile {}' \
                   ''.format(init_path, iter_path, self.type, self.term_file)
          command = 'qsub -e /dev/null -o /dev/null -d {} iota.process -F "{}"' \
                    ''.format(self.init.params.output, params)
        print command
        easy_run.fully_buffered(command, join_stdout_stderr=True)
      except Exception, e:
        print e

    # Send "all done" event to GUI
    try:
      evt = AllDone(tp_EVT_ALLDONE, -1, img_objects)
      wx.PostEvent(self.parent, evt)
    except TypeError, e:
      pass

  def full_proc_wrapper(self, input_entry):
    abort = os.path.isfile(self.term_file)
    try:
      proc_image_instance = ProcessImage(self.init, input_entry, self.type, abort)
      proc_image = proc_image_instance.run()
      return proc_image
    except Exception, e:
      raise e

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

class SpotFinderOneThread():
  def __init__(self, parent, processor, term_file):
    self.meta_parent = parent.parent
    self.processor = processor
    self.term_file = term_file

  def run(self, idx, datablock, img):
    if os.path.isfile(self.term_file):
      raise Exception('IOTA_TRACKER: Termination signal received!')
    else:
      observed = self.processor.find_spots(datablock=datablock)
      return [idx, len(observed), img]


class SpotFinderThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
  def __init__(self,
               parent,
               data_list,
               term_file,
               processor):
    Thread.__init__(self)
    self.parent = parent
    self.data_list = data_list
    self.term_file = term_file
    self.processor = processor
    #self.spotfinder = SpotFinderOneThread(self, processor, term_file)

  def run(self):
    try:
      parallel_map(iterable=self.data_list,
                   func=self.spf_wrapper,
                   callback=self.callback,
                   preserve_exception_message=True,
                   processes=None)
    except Exception, e:
      print 'SPOTFINDING THREAD:', e

    # Signal that this batch is finished
    try:
      if os.path.isfile(self.term_file):
        info = []
      else:
        info = self.data_list
      evt = SpotFinderOneDone(tp_EVT_SPFALLDONE, -1, info=info)
      wx.PostEvent(self.parent, evt)
    except TypeError:
      pass

  def spf_wrapper(self, img):
    if os.path.isfile(self.term_file):
      print 'TERMINATING {}, image {} of {}' \
            ''.format(img, self.data_list.index(img), len(self.data_list))
      raise Exception('Termination signal received!')
    else:
      if os.path.isfile(img):
        datablock = DataBlockFactory.from_filenames([img])[0]
        #info = self.spotfinder.run(self.data_list.index(img), datablock, img)
        #return info
        observed = self.processor.find_spots(datablock=datablock)
        return [int(self.data_list.index(img)), int(len(observed)), img]
      else:
        print 'DEBUG: FILE {} DOES NOT EXIST'.format(img)
        return [int(self.data_list.index(img)), 0, img]

  def callback(self, info):
    try:
      evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
      wx.PostEvent(self.parent, evt)
    except TypeError:
      pass
