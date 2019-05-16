'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 08/31/2018
Description : PRIME GUI Threading module
'''
from __future__ import absolute_import, division, print_function

import os
import wx
from threading import Thread

from libtbx import easy_run
from six.moves import map

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

# user = os.getlogin()
# icons = os.path.join(os.path.dirname(os.path.abspath(ct.__file__)), 'icons/')

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
               verbose=False):
    Thread.__init__(self)
    self.parent = parent
    self.prime_file = prime_file
    self.out_file = out_file
    self.command = command
    self.cmd_args = cmd_args
    self.signal_finished = signal_finished
    self.verbose = verbose

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

    if self.verbose:
      print(cmd)
      easy_run.fully_buffered(cmd, join_stdout_stderr=True).show_stdout()
    else:
      easy_run.fully_buffered(cmd, join_stdout_stderr=True)

    if self.signal_finished:
      evt = AllDone(tp_EVT_ALLDONE, -1)
      wx.PostEvent(self.parent, evt)
