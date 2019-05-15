from __future__ import absolute_import, division, print_function
from six.moves import range

'''
Author      : Lyubimov, A.Y.
Created     : 05/01/2016
Last Changed: 10/21/2018
Description : PRIME GUI dialogs module
'''

import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wxtbx import bitmaps

from iotbx import phil as ip

from iota.components.iota_utils import WxFlags, Capturing
from iota.components.iota_ui_base import BaseDialog, BaseBackendDialog
import iota.components.iota_ui_controls as ct

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
user = os.getlogin()


class PRIMEBaseBackendDialog(BaseBackendDialog):
  def __init__(self, parent,
               content_style='normal',
               label_style='bold',
               *args, **kwargs):
    BaseBackendDialog.__init__(self, parent,
                               content_style=content_style,
                               label_style=label_style,
                               *args, **kwargs)

  def get_target_file(self):
    dlg = wx.FileDialog(
      self, message="Select PRIME settings file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.FD_OPEN | wx.FD_CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]

      with open(filepath, 'r') as phil_file:
        phil_content = phil_file.read()
      return phil_content
    else:
      return None

def str_split(string, delimiters=(' ', ','), maxsplit=0):
  import re
  rexp = '|'.join(map(re.escape, delimiters))
  return re.split(rexp, string, maxsplit)


class PRIMEAdvancedOptions(PRIMEBaseBackendDialog):
  ''' Advanced Options Dialog'''

  def __init__(self, parent, phil=None, *args, **kwargs):
    PRIMEBaseBackendDialog.__init__(self, parent,
                                    backend_name='PRIME',
                                    phil=phil,
                                    phil_size=(500, 500),
                                    opt_size=(500, 500),
                                    *args, **kwargs)

    self.prime_phil = phil
    self.new_prime_phil = None
    self.pparams = self.prime_phil.extract()

    self.splitter.SplitVertically(self.options, self.phil_panel)

    # Target file input
    self.phil = ct.PHILBox(self.phil_panel,
                           btn_import=True,
                           btn_import_label='Import PHIL',
                           btn_export=False,
                           btn_default=False,
                           btn_pos='bottom',
                           ctr_size=(-1, 300),
                           ctr_value='')
    self.phil_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=5)

    # PRIME Options (there is some redundancy with the PHIL textbox)
    self.prm_options = wx.Panel(self.options)
    opt_box = wx.StaticBox(self.prm_options, label='Advanced Options')
    opt_box_sizer = wx.StaticBoxSizer(opt_box, wx.VERTICAL)
    self.prm_options.SetSizer(opt_box_sizer)

    # Resolution
    self.scale_res = ct.OptionCtrl(self.prm_options,
                                   label='Scale Resolution: ',
                                   label_size=(120, -1),
                                   label_style='normal',
                                   ctrl_size=wx.DefaultSize,
                                   items=[('low', 50),
                                          ('high', 1.5)])
    opt_box_sizer.Add(self.scale_res, flag=f.stack, border=10)

    self.merge_res = ct.OptionCtrl(self.prm_options,
                                   label='Merge Resolution: ',
                                   label_size=(120, -1),
                                   label_style='normal',
                                   ctrl_size=wx.DefaultSize,
                                   items=[('low', 50),
                                          ('high', 1.5)])
    opt_box_sizer.Add(self.merge_res, flag=f.stack, border=10)

    self.p_scale_res = ct.OptionCtrl(self.prm_options,
                                     label='Postref Scale Res.: ',
                                     label_size=(120, -1),
                                     label_style='normal',
                                     ctrl_size=wx.DefaultSize,
                                     items=[('low', 50),
                                            ('high', 1.5)])
    opt_box_sizer.Add(self.p_scale_res, flag=f.stack, border=10)

    self.p_cryst_res = ct.OptionCtrl(self.prm_options,
                                     checkbox=True,
                                     checkbox_label='Crystal Orientation: ',
                                     label_size=(120, -1),
                                     label_style='normal',
                                     ctrl_size=wx.DefaultSize,
                                     items=[('low', 50),
                                            ('high', 1.5)])
    opt_box_sizer.Add(self.p_cryst_res, flag=f.stack, border=10)

    self.p_rrange_res = ct.OptionCtrl(self.prm_options,
                                      checkbox=True,
                                      checkbox_label='Reflecting Range: ',
                                      label_size=(120, -1),
                                      label_style='normal',
                                      ctrl_size=wx.DefaultSize,
                                      items=[('low', 50),
                                             ('high', 1.5)])
    opt_box_sizer.Add(self.p_rrange_res, flag=f.stack, border=10)

    self.p_uc_res = ct.OptionCtrl(self.prm_options,
                                  checkbox=True,
                                  checkbox_label='Unit Cell: ',
                                  label_size=(120, -1),
                                  label_style='normal',
                                  ctrl_size=wx.DefaultSize,
                                  items=[('low', 50),
                                         ('high', 1.5)])
    opt_box_sizer.Add(self.p_uc_res, flag=f.stack, border=10)

    self.p_all_res = ct.OptionCtrl(self.prm_options,
                                   checkbox=True,
                                   checkbox_label='All Parameters: ',
                                   label_size=(120, -1),
                                   label_style='normal',
                                   ctrl_size=wx.DefaultSize,
                                   items=[('low', 50),
                                          ('high', 1.5)])
    opt_box_sizer.Add(self.p_all_res, flag=f.stack, border=10)

    self.btn_synch_res = wx.Button(self.prm_options,
                                   label="Set Res. Limits to Postref. Scale")
    opt_box_sizer.Add(self.btn_synch_res, flag=f.stack, border=10)


    # Target space group
    self.sg = ct.OptionCtrl(self.prm_options,
                            label='Space Group: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('spacegroup','P212121')])
    opt_box_sizer.Add(self.sg, flag=f.stack, border=10)

    # Target unit cell
    self.uc = ct.OptionCtrl(self.prm_options,
                            label='Unit Cell: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(300, -1),
                            items=[('unit_cell', '72 120 134 90 90 90')])
    self.anom = wx.CheckBox(self.prm_options, label='Anomalous')
    self.anom.SetValue(False)
    opt_box_sizer.Add(self.uc, flag=f.stack, border=10)
    opt_box_sizer.Add(self.anom, flag=f.stack, border=10)

    # CC cutoff
    self.cc = ct.OptionCtrl(self.prm_options,
                            label='CC cutoff: ',
                            label_size=(120, -1),
                            label_style='normal',
                            ctrl_size=(100, -1),
                            items=[('cc_cutoff', 0.25)])
    opt_box_sizer.Add(self.cc, flag=f.stack, border=10)

    # Pixel size
    self.pix = ct.OptionCtrl(self.prm_options,
                             label='Pixel size: ',
                             label_size=(120, -1),
                             label_style='normal',
                             ctrl_size=(100, -1),
                             items=[('pixel_size', 0.172)])
    opt_box_sizer.Add(self.pix, flag=wx.RIGHT | wx.LEFT | wx.TOP, border=10)

    self.cycles = ct.SpinCtrl(self.prm_options,
                              label='No. of Cycles:',
                              label_size=(120, -1),
                              label_style='normal',
                              ctrl_size=(60, -1))
    opt_box_sizer.Add(self.cycles, flag=wx.ALL, border=10)

    # self.options_sizer.Add(self.phil_sizer, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.options_sizer.Add(self.prm_options, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onImportPHIL, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onHideScript, self.btn_hide_script)
    self.Bind(wx.EVT_BUTTON, self.onResSynchronize, self.btn_synch_res)
    self.Bind(wx.EVT_CHOICE, self.onAdvanced, self.dlg_ctr.choice)

    self.show_hide_script()
    self.show_hide_advanced(show=False)
    self.Layout()
    self.options.SetupScrolling()
    self.read_param_phil()

  def onAdvanced(self, e):
    mode = self.dlg_ctr.choice.GetSelection()
    if mode == 0:
      self.show_hide_advanced(show=False)
    else:
      self.show_hide_advanced(show=True)

  def show_hide_advanced(self, show=False):
    if show:
      self.p_all_res.Show()
      self.p_cryst_res.Show()
      self.p_rrange_res.Show()
      self.p_scale_res.Show()
      self.p_uc_res.Show()
      self.btn_synch_res.Show()
    else:
      self.p_all_res.Hide()
      self.p_cryst_res.Hide()
      self.p_rrange_res.Hide()
      self.p_scale_res.Hide()
      self.p_uc_res.Hide()
      self.btn_synch_res.Hide()

    self.options.Layout()
    self.options.SetupScrolling()

  def read_param_phil(self):
    # Put in PHIL script
    self.generate_phil_string()
    self.phil.ctr.SetValue(self.phil_string)

    # Set values to default parameters
    self.scale_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.scale.d_max))
    self.scale_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.scale.d_min))
    self.merge_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.merge.d_max))
    self.merge_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.merge.d_min))
    self.p_scale_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.postref.scale.d_max))
    self.p_scale_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.postref.scale.d_min))
    self.p_cryst_res.toggle_boxes(flag_on=self.pparams.postref.crystal_orientation.flag_on)
    self.p_cryst_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.postref.crystal_orientation.d_max))
    self.p_cryst_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.postref.crystal_orientation.d_min))
    self.p_rrange_res.toggle_boxes(flag_on=self.pparams.postref.reflecting_range.flag_on)
    self.p_rrange_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.postref.reflecting_range.d_max))
    self.p_rrange_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.postref.reflecting_range.d_min))
    self.p_uc_res.toggle_boxes(flag_on=self.pparams.postref.unit_cell.flag_on)
    self.p_uc_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.postref.unit_cell.d_max))
    self.p_uc_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.postref.unit_cell.d_min))
    self.p_all_res.toggle_boxes(flag_on=self.pparams.postref.allparams.flag_on)
    self.p_all_res.low.SetValue(
      '{:4.2f}'.format(self.pparams.postref.allparams.d_max))
    self.p_all_res.high.SetValue(
      '{:4.2f}'.format(self.pparams.postref.allparams.d_min))

    self.sg.spacegroup.SetValue(str(self.pparams.target_space_group))
    if str(self.pparams.target_unit_cell).lower() != 'none':
      uc = ' '.join(list(map(str, self.pparams.target_unit_cell.parameters())))
      self.uc.unit_cell.SetValue(uc)
    else:
      self.uc.unit_cell.SetValue(str(self.pparams.target_unit_cell))
    self.anom.SetValue(self.pparams.target_anomalous_flag)
    self.cc.cc_cutoff.SetValue(str(self.pparams.frame_accept_min_cc))
    self.pix.pixel_size.SetValue(str(self.pparams.pixel_size_mm))
    self.cycles.ctr.SetValue(int(self.pparams.n_postref_cycle))

  def generate_phil_string(self):
    with Capturing() as txt_output:
      self.prime_phil.show()
    self.phil_string = ''
    for one_output in txt_output:
      self.phil_string += one_output + '\n'

  def onHideScript(self, e):
    self.opt_size = self.options.GetSize()
    self.phil_size = self.phil_panel.GetSize()
    self.sash_position = self.opt_size[0]
    self.show_hide_script(initialized=True)

  def onImportPHIL(self, e):
    phil_content = self.get_target_file()
    if phil_content is not None:
      self.phil.ctr.SetValue(phil_content)

  def onResSynchronize(self, e):
    self.synchronize_resolution()

  def synchronize_resolution(self):
    # Set all activated resolution ranges to Merge Resolution
    scale_low = self.merge_res.low.GetValue()
    scale_high = self.merge_res.high.GetValue()
    self.p_scale_res.low.SetValue(scale_low)
    self.p_scale_res.high.SetValue(scale_high)
    if self.p_cryst_res.toggle.GetValue():
      self.p_cryst_res.low.SetValue(scale_low)
      self.p_cryst_res.high.SetValue(scale_high)
    if self.p_rrange_res.toggle.GetValue():
      self.p_rrange_res.low.SetValue(scale_low)
      self.p_rrange_res.high.SetValue(scale_high)
    if self.p_uc_res.toggle.GetValue():
      self.p_uc_res.low.SetValue(scale_low)
      self.p_uc_res.high.SetValue(scale_high)
    if self.p_all_res.toggle.GetValue():
      self.p_all_res.low.SetValue(scale_low)
      self.p_all_res.high.SetValue(scale_high)

  def onOK(self, e):
    self.phil_script = ip.parse(self.phil.ctr.GetValue())

    prime_phil_text = '\n'.join([
      'scale',
      '{',
      '  d_min = {}'.format(float(self.scale_res.high.GetValue())),
      '  d_max = {}'.format(float(self.scale_res.low.GetValue())),
      '}',
      'postref',
      '{',
      '  scale',
      '  {',
      '    d_min = {}'.format(float(self.p_scale_res.high.GetValue())),
      '    d_max = {}'.format(float(self.p_scale_res.low.GetValue())),
      '  }',
      '  crystal_orientation',
      '  {',
      '    flag_on = {}'.format(self.p_cryst_res.toggle.GetValue()),
      '    d_min = {}'.format(float(self.p_cryst_res.high.GetValue())),
      '    d_max = {}'.format(float(self.p_cryst_res.low.GetValue())),
      '  }',
      '  reflecting_range',
      '  {',
      '    flag_on = {}'.format(self.p_rrange_res.toggle.GetValue()),
      '    d_min = {}'.format(float(self.p_rrange_res.high.GetValue())),
      '    d_max = {}'.format(float(self.p_rrange_res.low.GetValue())),
      '  }',
      '  unit_cell',
      '  {',
      '    flag_on = {}'.format(self.p_uc_res.toggle.GetValue()),
      '    d_min = {}'.format(float(self.p_uc_res.high.GetValue())),
      '    d_max = {}'.format(float(self.p_uc_res.low.GetValue())),
      '  }',
      '  allparams',
      '  {',
      '    flag_on = {}'.format(self.p_all_res.toggle.GetValue()),
      '    d_min = {}'.format(float(self.p_all_res.high.GetValue())),
      '    d_max = {}'.format(float(self.p_all_res.low.GetValue())),
      '  }',
      '}',
      'merge',
      '{',
      '  d_min = {}'.format(float(self.merge_res.high.GetValue())),
      '  d_max = {}'.format(float(self.merge_res.low.GetValue())),
      '}',
      'target_unit_cell = {}'.format(self.uc.unit_cell.GetValue()),
      'target_space_group = {}'.format(self.sg.spacegroup.GetValue()),
      'target_anomalous_flag = {}'.format(self.anom.GetValue()),
      'pixel_size_mm = {}'.format(self.pix.pixel_size.GetValue()),
      'frame_accept_min_cc = {}'.format(self.cc.cc_cutoff.GetValue()),
      'n_postref_cycle = {}'.format(self.cycles.ctr.GetValue())
    ])

    self.new_prime_phil = ip.parse(prime_phil_text)
    e.Skip()

class PRIMEPreferences(BaseDialog):
  def __init__(self, parent, phil=None, *args, **kwargs):
    BaseDialog.__init__(self, parent, *args, **kwargs)

    self.pparams = phil.extract()
    self.pref_phil = None

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

    self.nodes = ct.SpinCtrl(self,
                             label='No. of nodes',
                             label_size=(120, -1),
                             ctrl_value=12,
                             ctrl_size=(80, -1),
                             ctrl_min=1,
                             ctrl_max=5000)
    self.nodes.Disable()
    vbox.Add(self.nodes, flag=wx.ALL, border=10)

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

    self.set_choices()

  def set_choices(self):
    # Set queue to default value
    queue = self.pparams.queue.qname
    method = self.pparams.queue.mode
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

    self.nodes.ctr.SetValue(self.pparams.queue.n_nodes)
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
      self.nodes.Enable()
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

    if self.method != 'multiprocessing':
      n_nodes = self.nodes.ctr.GetValue()
    else:
      n_nodes = 12

    prime_prefs_string = '\n'.join([
      'queue {',
      '  mode = {}'.format(self.method),
      '  qname = {}'.format(self.queue),
      '  n_nodes = {}'.format(n_nodes),
      '}'
    ])

    self.pref_phil = ip.parse(prime_prefs_string)
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
    self.recovery_mode = 0

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
    self.dlg_ctr = ct.DialogButtonsCtrl(self, preset='OK_CANCEL',
                                        choices=['everything', 'settings only'],
                                        choice_label='Recover: ')
    self.main_sizer.Add(self.dlg_ctr, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

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
        self.recovery_mode = self.dlg_ctr.choice.GetSelection()
    e.Skip()
