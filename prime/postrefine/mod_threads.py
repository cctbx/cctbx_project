from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 06/20/2019
Description : PRIME GUI Threading module
'''

import os
import wx
from threading import Thread

from iota.threads.iota_threads import CustomRun

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = 'python'
elif wx.Platform == '__WXMAC__':
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif (wx.Platform == '__WXMSW__'):
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9

def str_split(string, delimiters=(' ', ','), maxsplit=0):
  import re
  rexp = '|'.join(map(re.escape, delimiters))
  return re.split(rexp, string, maxsplit)


# -------------------------------- Threading --------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

class AllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class PRIMEThread(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               prime_file,
               out_file=None,
               command=None,
               cmd_args=None,
               signal_finished=False,
               debug=False,
               verbose=False):
    Thread.__init__(self)
    self.parent = parent
    self.prime_file = prime_file
    self.out_file = out_file
    self.command = command
    self.cmd_args = cmd_args
    self.signal_finished = signal_finished
    self.verbose = verbose
    self.debug = debug
    self.job = None

  def run(self):
    if os.path.isfile(self.out_file):
      os.remove(self.out_file)
    if self.command is None:
      if self.cmd_args is None:
        args = ''
      else:
        args = self.cmd_args

      cmd = 'prime.run {} {}'.format(self.prime_file, args)
    else:
      cmd = self.command

    if self.debug:
      from libtbx import easy_run
      if self.verbose:
        print (cmd)
        easy_run.fully_buffered(cmd, join_stdout_stderr=True).show_stdout()
      else:
        easy_run.fully_buffered(cmd, join_stdout_stderr=True)
    else:
      try:
        self.job = CustomRun(command=cmd, join_stdout_stderr=True)
        self.job.run()
      except Exception as e:
        print ("PRIME ERROR: ", e)

    if self.signal_finished:
      evt = AllDone(tp_EVT_ALLDONE, -1)
      wx.PostEvent(self.parent, evt)

  def abort(self):
    # TODO: put in an LSF kill command
    if self.job:
      try:
        self.job.kill_thread()
      except Exception as e:
        print ('PRIME THREAD ERROR: Cannot terminate thread! {}'.format(e))
