from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 11/03/2017
Description : PRIME GUI dialogs module
'''

import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wxtbx import bitmaps

import iota.components.iota_controls as ct
from iota.components.iota_misc import WxFlags

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

f = WxFlags()


class BaseDialog(wx.Dialog):
  def __init__(self, parent, style=wx.DEFAULT_DIALOG_STYLE,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, style=style, *args, **kwargs)

    self.envelope = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.envelope.Add(self.main_sizer, 1, flag=wx.EXPAND | wx.ALL, border=5)
    self.SetSizer(self.envelope)

    if label_style == 'normal':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    elif label_style == 'bold':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    elif label_style == 'italic':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.NORMAL)
    elif label_style == 'italic_bold':
      self.font = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.BOLD)

    if content_style == 'normal':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
    elif content_style == 'bold':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.NORMAL, wx.BOLD)
    elif content_style == 'italic':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.NORMAL)
    elif content_style == 'italic_bold':
      self.cfont = wx.Font(norm_font_size, wx.DEFAULT, wx.ITALIC, wx.BOLD)

user = os.getlogin()
icons = os.path.join(os.path.dirname(os.path.abspath(ct.__file__)), 'icons/')

def str_split(string, delimiters=(' ', ','), maxsplit=0):
  import re
  rexp = '|'.join(map(re.escape, delimiters))
  return re.split(rexp, string, maxsplit)


class PRIMEAdvancedOptions(BaseDialog):
  ''' Advanced Options Dialog'''

  def __init__(self, parent, phil=None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 300),
                        *args, **kwargs)

    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self, size=(-1, 300))
    options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(options_sizer)

    phil_box = wx.StaticBox(self, label='PRIME Input')
    phil_box_sizer = wx.StaticBoxSizer(phil_box, wx.VERTICAL)

    # Target file input
    self.phil = ct.PHILBox(self,
                           btn_import=False,
                           btn_export=False,
                           btn_default=False,
                           ctr_size=(500, 200),
                           ctr_value='')
    phil_box_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # PRIME Options (there is some redundancy with the PHIL textbox)
    opt_box = wx.StaticBox(self.options, label='Advanced Options')
    opt_box_sizer = wx.StaticBoxSizer(opt_box, wx.VERTICAL)

    # Resolution
    self.res_override = wx.CheckBox(self.options, label='Resolution override')
    self.res = ct.OptionCtrl(self.options,
                             label='Resolution: ',
                             label_size=(120, -1),
                             label_style='normal',
                             ctrl_size=wx.DefaultSize,
                             items=[('high', 50),
                                    ('low', 1.5)])
    opt_box_sizer.Add(self.res_override, flag=f.stack, border=10)
    opt_box_sizer.Add(self.res, flag=f.stack, border=10)

    # Target space group
    self.sg = ct.OptionCtrl(self.options,
                            label='Space Group: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('spacegroup','P212121')])
    opt_box_sizer.Add(self.sg, flag=f.stack, border=10)

    # Target unit cell
    self.uc = ct.OptionCtrl(self.options,
                            label='Unit Cell: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(300, -1),
                            items=[('unit_cell', '72 120 134 90 90 90')])
    self.anom = wx.CheckBox(self.options, label='Anomalous')
    self.anom.SetValue(False)
    opt_box_sizer.Add(self.uc, flag=f.stack, border=10)
    opt_box_sizer.Add(self.anom, flag=f.stack, border=10)

    # CC cutoff
    self.cc = ct.OptionCtrl(self.options,
                            label='CC cutoff: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('cc_cutoff', 0.25)])
    opt_box_sizer.Add(self.cc, flag=f.stack, border=10)

    # Pixel size
    self.pix = ct.OptionCtrl(self.options,
                             label='Pixel size: ',
                             label_size=(120, -1),
                             label_style='normal',
                             ctrl_size=(100, -1),
                             items=[('pixel_size', 0.172)])
    opt_box_sizer.Add(self.pix, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    self.cycles = ct.SpinCtrl(self.options,
                              label='No. of Cycles:',
                              label_size=(120, -1),
                              label_style='normal',
                              ctrl_size=(60, -1))
    opt_box_sizer.Add(self.cycles, flag=wx.ALL, border=10)

    # Dialog controls
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    options_sizer.Add(opt_box_sizer, 1, flag=f.expand, border=10)
    self.main_sizer.Add(phil_box_sizer, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.options, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.options.SetupScrolling()


class PRIMEPreferences(wx.Dialog):
  def __init__(self, *args, **kwargs):
    super(PRIMEPreferences, self).__init__(*args, **kwargs)

    self.method = None
    self.queue = None
    main_sizer = wx.BoxSizer(wx.VERTICAL)

    main_box = wx.StaticBox(self, label='PRIME Preferences')
    vbox = wx.StaticBoxSizer(main_box, wx.VERTICAL)

    self.SetSizer(main_sizer)

    q_choices = ['psanaq', 'psnehq', 'psfehq'] + ['custom']
    self.queues = ct.ChoiceCtrl(self,
                                label='Queue:',
                                label_size=(120, -1),
                                label_style='bold',
                                ctrl_size=wx.DefaultSize,
                                choices=q_choices)
    vbox.Add(self.queues, flag=wx.ALL, border=10)

    self.custom_queue = ct.OptionCtrl(self,
                                      items=[('cqueue', '')],
                                      label='Custom Queue:',
                                      label_size=(120, -1),
                                      label_style='normal',
                                      ctrl_size=(150, -1))
    self.custom_queue.Disable()
    vbox.Add(self.custom_queue, flag=wx.ALL, border=10)

    mp_choices = ['multiprocessing', 'bsub']
    self.mp_methods = ct.ChoiceCtrl(self,
                                    label='Method:',
                                    label_size=(120, -1),
                                    label_style='bold',
                                    ctrl_size=wx.DefaultSize,
                                    choices=mp_choices)
    vbox.Add(self.mp_methods, flag=wx.ALL, border=10)

    main_sizer.Add(vbox, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.Bind(wx.EVT_CHOICE, self.onQueue, self.queues.ctr)
    self.Bind(wx.EVT_CHOICE, self.onMethod, self.mp_methods.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def set_choices(self, method, queue):
    # Set queue to default value
    if queue is None:
      queue = 'None'
    inp_queue = self.queues.ctr.FindString(queue)
    if inp_queue != wx.NOT_FOUND:
      self.queues.ctr.SetSelection(inp_queue)
    else:
      self.custom_queue.Enable()
      self.custom_queue.cqueue.SetValue(queue)
      custom_queue = self.queues.ctr.FindString('custom')
      self.queues.ctr.SetSelection(custom_queue)

    # Set method to default value
    inp_method = self.mp_methods.ctr.FindString(str(method))
    if inp_method != wx.NOT_FOUND:
      self.mp_methods.ctr.SetSelection(inp_method)

    self.check_method()

  def onMethod(self, e):
    self.check_method()

  def check_method(self):
    choice = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    if choice == 'multiprocessing':
      self.queues.Disable()
      self.custom_queue.Disable()
    else:
      self.queues.Enable()
      queue = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
      if queue == 'custom':
        self.custom_queue.Enable()

  def onQueue(self, e):
    choice = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
    if choice == 'custom':
      self.custom_queue.Enable()
    else:
      self.custom_queue.Disable()

  def onOK(self, e):
    self.method = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    queue_selection = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
    if queue_selection == 'custom':
      if self.custom_queue.cqueue.GetValue() == '':
        wx.MessageBox('Please choose or enter a queue', wx.OK)
      else:
        self.queue = self.custom_queue.cqueue.GetValue()
        e.Skip()
    else:
      self.queue = queue_selection
      e.Skip()


class TextFileView(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               contents=None,
               *args, **kwargs):

    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP

    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)

    self.txt_panel = ScrolledPanel(self)
    self.txt_sizer = wx.BoxSizer(wx.VERTICAL)
    self.txt_panel.SetSizer(self.txt_sizer)

    self.txt = wx.StaticText(self.txt_panel, label=contents)
    self.txt_sizer.Add(self.txt)

    self.txt_panel.SetupScrolling()
    self.main_sizer.Add(self.txt_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)

class RecoveryDialog(BaseDialog):

  def __init__(self,
               parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP
    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(400, 400),
                        title='Recover Completed Run',
                        *args, **kwargs)


    self.pathlist = wx.ListCtrl(self, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    self.selected = None

    self.pathlist.InsertColumn(0, "")
    self.pathlist.InsertColumn(1, "Status")
    self.pathlist.InsertColumn(2, "Path")

    bmps = wx.ImageList(width=16, height=16)
    finished_bmp = bitmaps.fetch_icon_bitmap('actions', 'button_ok', size=16)
    aborted_bmp = bitmaps.fetch_icon_bitmap('actions', 'cancel', size=16)
    unknown_bmp = bitmaps.fetch_icon_bitmap('actions', 'status_unknown',
                                            size=16)
    bmps.Add(finished_bmp)
    bmps.Add(aborted_bmp)
    bmps.Add(unknown_bmp)
    self.pathlist.AssignImageList(bmps, which=1)

    self.main_sizer.Add(self.pathlist, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)


  def insert_paths(self, pathlist):
    for i in range(len(pathlist)):
      logfile = os.path.join(pathlist[i], 'log.txt')
      if os.path.isfile(logfile):
        with open(logfile, 'r') as lf:
          log_contents = lf.readlines()
      else:
        log_contents = []
        img_id = 2
        status = "Unknown"

      stats_folder = os.path.join(pathlist[i], 'stats')
      stat_files = [os.path.join(stats_folder, f) for f in
                    os.listdir(stats_folder) if f.endswith('stat')]
      if stat_files != []:
        assert len(stat_files) == 1
        stat_file = stat_files[0]
        if log_contents != [] and 'Finished' in log_contents[-1]:
          img_id = 0
          status = 'Finished'
        elif log_contents != [] and os.path.isfile(stat_file):
          img_id = 1
          status = 'Aborted'
      else:
        img_id = 2
        status = "Unknown"

      idx = self.pathlist.InsertImageItem(i, img_id)
      self.pathlist.SetStringItem(idx, 1, status)
      self.pathlist.SetStringItem(idx, 2, pathlist[i])
      self.pathlist.SetColumnWidth(0, width=-1)
      self.pathlist.SetColumnWidth(1, width=-1)
      self.pathlist.SetColumnWidth(2, width=-1)


  def onOK(self, e):
    for i in range(self.pathlist.GetItemCount()):
      if self.pathlist.IsSelected(i):

        self.selected = [self.pathlist.GetItemText(i, col=1),
                         self.pathlist.GetItemText(i, col=2)]
    e.Skip()
