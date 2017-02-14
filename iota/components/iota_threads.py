from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 01/06/2017
Description : IOTA GUI Threads and PostEvents
'''

import os
import wx
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep
from libtbx import easy_run

from iota.components.iota_utils import GenerateInput
import iota.components.iota_image as img

# -------------------------------- Threading --------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

tp_EVT_IMGDONE = wx.NewEventType()
EVT_IMGDONE = wx.PyEventBinder(tp_EVT_IMGDONE, 1)

class ImageFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, image_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_list = image_list
  def GetValue(self):
    return self.image_list

class AllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, img_objects=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_objects = img_objects
  def GetValue(self):
    return self.image_objects

class ProcessImage():
  ''' Wrapper class to do full processing of an image '''
  def __init__(self, init, input_entry, input_type = 'image'):
    self.init = init
    self.input_entry = input_entry
    self.input_type = input_type
  def run(self):
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
      #result_file = os.path.splitext(img_object.obj_file)[0] + '.fin'
      #ep.dump(result_file, img_object)
      return img_object

class ProcThread(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''
  def __init__(self,
               parent,
               init,
               iterable,
               input_type='image'):
    Thread.__init__(self)
    self.parent = parent
    self.init = init
    self.iterable = iterable
    self.type = input_type

  def run(self):
    if self.init.params.mp_method == 'multiprocessing':
      img_objects = parallel_map(iterable=self.iterable,
                                 func = self.full_proc_wrapper,
                                 processes=self.init.params.n_processors)
    else:
      # write iterable
      img_objects = None
      queue = self.init.params.mp_queue
      iter_path = os.path.join(self.init.params.output, 'iter.cfg')
      init_path = os.path.join(self.init.params.output, 'init.cfg')
      nproc = self.init.params.n_processors
      ep.dump(iter_path, self.iterable)
      ep.dump(init_path, self.init)
      if self.init.params.mp_method == 'lsf':
        try:
          command = 'bsub -q {} -n {} iota.process {} --files {} --type {}' \
                    ''.format(queue, nproc, init_path, iter_path, self.type)
        except Exception, e:
          print e
      if self.init.params.mp_method == 'torq':



        params = '{} --files {} --type {}'.format(init_path, iter_path, self.type)
        try:
          command = 'qsub -e /dev/null -o /dev/null -d {} iota.process -F "{}"' \
                    ''.format(self.init.params.output, params)
        except Exception, e:
          print e

      print command
      easy_run.fully_buffered(command, join_stdout_stderr=True)

    # Send "all done" event to GUI
    evt = AllDone(tp_EVT_ALLDONE, -1, img_objects)
    wx.PostEvent(self.parent, evt)

  def full_proc_wrapper(self, input_entry):
    proc_image_instance = ProcessImage(self.init, input_entry, self.type)
    proc_image = proc_image_instance.run()
    return proc_image

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
    ginp = GenerateInput()

    input_entries = []
    for path in self.image_paths:
      if os.path.isdir(path):
        input_entries.extend([i for i in path if i != None])
    ext_file_list = ginp.make_input_list(input_entries)
    old_file_list = [i[2] for i in self.image_list]
    new_file_list = [i for i in ext_file_list if i not in old_file_list]

    # Test if file is an image pickle or raw image (may be slow!)
    tested_file_list = [i for i in new_file_list if ginp.get_file_type(i) !=
                        'not image']

    # Generate list of new images
    new_img = [[i, len(ext_file_list) + 1, j] for i, j in enumerate(
      tested_file_list, len(old_file_list) + 1)]

    evt = ImageFinderAllDone(tp_EVT_IMGDONE, -1, image_list=new_img)
    wx.PostEvent(self.parent, evt)