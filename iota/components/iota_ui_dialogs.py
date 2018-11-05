from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 01/17/2017
Last Changed: 11/05/2018
Description : IOTA GUI Dialogs
'''

import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.buttons import GenToggleButton
from wxtbx import bitmaps

from iotbx import phil as ip

import iota.components.iota_ui_controls as ct
from iota.components.iota_input import master_phil
from iota.components.iota_utils import UnicodeCharacters, WxFlags, noneset, \
  norm_font_size

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif wx.Platform == '__WXMAC__':
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
elif (wx.Platform == '__WXMSW__'):
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9

# Initialize unicode font and wx flags
u = UnicodeCharacters()
f = WxFlags()


# ------------------------------ Base Classes -------------------------------- #

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


class BaseBackendDialog(BaseDialog):
  def __init__(self, parent, phil,
               backend_name = 'BACKEND',
               target=None,
               content_style='normal',
               label_style='bold',
               opt_size=(500, 500),
               phil_size=(500, 500),
               *args, **kwargs):
    BaseDialog.__init__(self, parent,
                        content_style=content_style,
                        label_style=label_style,
                        *args, **kwargs)

    self.parent = parent
    self.target_phil = target
    self.backend = backend_name
    self.params = phil.extract()
    self.opt_size = opt_size
    self.phil_size = phil_size
    self.sash_position = None

    self.splitter = wx.SplitterWindow(self, style=wx.SP_LIVE_UPDATE |
                                                  wx.SP_3DSASH |
                                                  wx.SP_NOBORDER)

    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self.splitter, size=self.opt_size)
    self.options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(self.options_sizer)

    # Create PHIL panel
    phil_label = "{} Target Settings".format(backend_name)
    self.phil_panel = wx.Panel(self.splitter, size=self.opt_size)
    phil_box = wx.StaticBox(self.phil_panel, label=phil_label)
    self.phil_sizer = wx.StaticBoxSizer(phil_box, wx.VERTICAL)
    self.phil_panel.SetSizer(self.phil_sizer)

    # Dialog control
    self.dlg_ctr = ct.DialogButtonsCtrl(self, preset='PROC_DIALOG')

    # Splitter button
    self.btn_hide_script = GenToggleButton(self, label='Show Script >>>')
    self.show_hide_script()
    self.btn_hide_script.SetValue(False)

    self.main_sizer.Add(self.btn_hide_script, flag=wx.ALIGN_RIGHT | wx.ALL,
                        border=10)
    self.main_sizer.Add(self.splitter, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.dlg_ctr,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.RIGHT,
                        border=10)

  def show_hide_script(self, initialized=False):
    if self.btn_hide_script.GetValue():
      if initialized:
        h = self.GetSize()[1]
        w = self.GetSize()[0] + self.phil_size[0]
        self.SetSize((w, h))
      self.splitter.SplitVertically(self.options, self.phil_panel)
      self.splitter.SetSashPosition(self.sash_position)
      self.phil_panel.SetSize(self.phil_size)
      self.options.SetSize(self.opt_size)
      self.btn_hide_script.SetLabel('<<< Hide Script')
    else:
      h = self.GetSize()[1]
      w = self.GetSize()[0] - self.phil_size[0]
      self.SetSize((w, h))
      self.splitter.Unsplit()
      self.phil_panel.SetSize(self.phil_size)
      self.options.SetSize(self.opt_size)
      self.btn_hide_script.SetLabel('Show Script >>>')
    self.splitter.SizeWindows()

  def get_target_file(self):
    dlg = wx.FileDialog(
      self, message="Select CCTBX.XFEL target file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.OPEN | wx.CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]

      with open(filepath, 'r') as phil_file:
        phil_content = phil_file.read()
      return phil_content
    else:
      return None

  def write_default_phil(self):
    if str.lower(self.backend) in ('cctbx.xfel', 'dials'):
      method = 'cctbx.xfel'
    elif str.lower(self.backend) in ('cctbx', 'ha14', 'labelit'):
      method = 'ha14'
    else:
      method = 'current'
    from iota.components.iota_input import write_defaults
    default_phil, _ = write_defaults(method=method, write_target_file=False,
                                     write_param_file=False)
    self.target_phil = '\n'.join(default_phil)


class BasePanel(wx.Panel):
  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(800, 500))

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

# ---------------------------------------------------------------------------- #

class IOTAPreferences(BaseDialog):
  """ Class for dialog that houses IOTA interface preferences, e.g.:
        - multiprocessing / queue settings
        - monitor mode settings
        - miscellaneous interface-only settings """

  def __init__(self, parent, phil=None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP

    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)

    # Import current PHIL and set params to current values
    if phil is None:
      self.params = master_phil.extract()
    else:
      self.params = phil.extract()
    self.method = self.params.mp.method
    self.queue = self.params.mp.queue
    self.monitor_mode = self.params.gui.monitor_mode
    self.mm_timeout = self.params.gui.monitor_mode_timeout
    self.mm_timeout_len = self.params.gui.monitor_mode_timeout_length
    self.random_subset = self.params.advanced.random_sample.flag_on
    self.random_subset_number = self.params.advanced.random_sample.number
    self.image_range = self.params.advanced.image_range.flag_on
    self.image_range_string = self.params.advanced.image_range.range

    # Queue Preferences
    queue_box = wx.StaticBox(self, label='Multiprocessing Preferences')
    queue_sizer = wx.StaticBoxSizer(queue_box, wx.VERTICAL)

    mp_choices = ['multiprocessing', 'lsf', 'torq']
    self.mp_methods = ct.ChoiceCtrl(self,
                                    label='Method:',
                                    label_size=(120, -1),
                                    label_style='bold',
                                    ctrl_size=wx.DefaultSize,
                                    choices=mp_choices)
    queue_sizer.Add(self.mp_methods, flag=wx.ALL,border=10)

    q_choices = ['psanaq', 'psnehq', 'psfehq'] + ['custom']
    self.queues = ct.ChoiceCtrl(self,
                                label='Queue:',
                                label_size=(120, -1),
                                label_style='bold',
                                ctrl_size=wx.DefaultSize,
                                choices=q_choices)
    queue_sizer.Add(self.queues, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM,
                    border=10)

    self.custom_queue = ct.OptionCtrl(self, items=[('cqueue', '')],
                                      label='Custom Queue:',
                                      label_size=(120, -1), ctrl_size=(150, -1))
    queue_sizer.Add(self.custom_queue, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM,
                    border=10)

    # Advanced Preferences
    adv_box = wx.StaticBox(self, label='Advanced Preferences')
    adv_sizer = wx.StaticBoxSizer(adv_box, wx.VERTICAL)

    # Backend choice (NOT RECOMMENDED!)
    self.proc_backend = ct.ChoiceCtrl(self,
                                 label='Processing Backend:',
                                 label_size=(150, -1),
                                 choices=['cctbx.xfel', 'cctbx.xfel HA14'],
                                 ctrl_size=(200, -1))
    self.proc_backend.ctr.SetSelection(0)
    adv_sizer.Add(self.proc_backend, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM,
                  border=10)

    # Viewer preferences
    v_choices = ['dials.image_viewer',
                 'cctbx.image_viewer',
                 'distl.image_viewer',
                 'cxi.view']
    self.viewers = ct.ChoiceCtrl(self,
                                 label='Image Viewer:',
                                 label_size=(120, -1),
                                 label_style='bold',
                                 ctrl_size=wx.DefaultSize,
                                 choices=v_choices)
    adv_sizer.Add(self.viewers, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM, border=10)

    # Monitor Mode preferences
    self.chk_cont_mode = wx.CheckBox(self, label='Process in Monitor Mode')
    self.chk_cont_mode.SetValue(False)
    self.chk_mm_timeout = wx.CheckBox(self, label='Monitor Mode Timeout')
    self.chk_mm_timeout.SetValue(False)
    self.opt_timeout = ct.OptionCtrl(self, items=[('timeout', '')],
                                     label='Timeout (sec):',
                                     label_size=(120, -1), ctrl_size=(150, -1))
    self.chk_mm_timeout.Disable()
    self.opt_timeout.Disable()
    adv_sizer.Add(self.chk_cont_mode, flag=f.stack, border=10)
    adv_sizer.Add(self.chk_mm_timeout, flag=f.stack,  border=10)
    adv_sizer.Add(self.opt_timeout, flag=f.stack, border=10)

    # Sub-sample preferences
    self.chk_image_range = ct.OptionCtrl(self, items=[('range', '')],
                                         label_size=wx.DefaultSize,
                                         checkbox=True,
                                         checkbox_label='Image Range: ',
                                         ctrl_size=(200, -1))
    adv_sizer.Add(self.chk_image_range, flag=f.stack, border=10)

    self.chk_random_sample = ct.SpinCtrl(self, checkbox=True,
                                         checkbox_label='Random subset: ',
                                         ctrl_value=100, ctrl_size=(80, -1),
                                         ctrl_max=5000)
    adv_sizer.Add(self.chk_random_sample, flag=f.stack, border=10)

    # PRIME prefix
    self.prime_prefix = ct.OptionCtrl(self,
                                      label='PRIME prefix',
                                      items=[('prefix', 'prime')],
                                      label_size=(160, -1),
                                      ctrl_size=(200, -1))
    adv_sizer.Add(self.prime_prefix, flag=wx.ALL | wx.EXPAND, border=10)

    # Temp folder
    self.temp_folder = ct.InputCtrl(self,
                                    label='Temp folder:',
                                    buttons=True)
    adv_sizer.Add(self.temp_folder, flag=wx.ALL | wx.EXPAND, border=10)

    self.main_sizer.Add(queue_sizer, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(adv_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_CHOICE, self.onQueue, self.queues.ctr)
    self.Bind(wx.EVT_CHOICE, self.onMethod, self.mp_methods.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_CHECKBOX, self.onMonitor, self.chk_cont_mode)
    self.Bind(wx.EVT_CHECKBOX, self.onTimeout, self.chk_mm_timeout)
    self.Bind(wx.EVT_BUTTON, self.onTempBrowse, self.temp_folder.btn_browse)

    self.Fit()

  def onTempBrowse(self, e):
    """ On clicking the Browse button: show the DirDialog and populate 'Output'
        box w/ selection """
    dlg = wx.DirDialog(self, "Choose the temporary output directory:",
                       style=wx.DD_DEFAULT_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
      self.temp_folder.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

  def set_choices(self):
    # Set queue to default value
    if self.queue is None:
      self.queue = 'None'
    inp_queue = self.queues.ctr.FindString(self.queue)
    if inp_queue != wx.NOT_FOUND:
      self.queues.ctr.SetSelection(inp_queue)
    else:
      self.custom_queue.Enable()
      self.custom_queue.cqueue.SetValue(self.queue)
      custom_sel = self.queues.ctr.FindString('custom')
      self.queues.ctr.SetSelection(custom_sel)

    # Set method to default value
    inp_method = self.mp_methods.ctr.FindString(noneset(self.method))
    if inp_method != wx.NOT_FOUND:
      self.mp_methods.ctr.SetSelection(inp_method)
    self.check_method()

    # Set backend values

    if self.params.advanced.processing_backend == 'cctbx.xfel':
      pos = 0
    elif self.params.advanced.processing_backend == 'ha14':
      pos = 1
    else:
      pos = 0
    self.proc_backend.ctr.SetSelection(pos)

    # Set viewer values
    try:   # Need this for backwards compatibility with old scripts
      viewer = self.params.gui.image_viewer
    except Exception:
      viewer = 'dials.image_viewer'
    v_selection = self.viewers.ctr.FindString(viewer)
    if v_selection != wx.NOT_FOUND:
      self.viewers.ctr.SetSelection(v_selection)

    # Set Monitor Mode values
    if self.monitor_mode:
      self.chk_cont_mode.SetValue(True)
      self.chk_mm_timeout.Enable()
      if self.mm_timeout:
        self.chk_mm_timeout.SetValue(True)
        self.opt_timeout.Enable()
        self.opt_timeout.timeout.SetValue(str(self.mm_timeout_len))

    # Set subset values
    if self.random_subset:
      self.chk_random_sample.toggle_boxes()
      self.chk_random_sample.ctr.SetValue(self.random_subset_number)

    if self.image_range:
      self.chk_image_range.toggle_boxes(flag_on=self.image_range)
      self.chk_image_range.range.SetValue(noneset(self.image_range_string))

    # Set PRIME prefix
    self.prime_prefix.prefix.SetValue(self.params.advanced.prime_prefix)

    # Set temp folder
    if self.params.advanced.temporary_output_folder is not None:
      self.temp_folder.ctr.SetValue(self.params.advanced.temporary_output_folder)


  def onMethod(self, e):
    self.check_method()

  def check_method(self):
    choice = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    if choice == 'lsf':
      self.queues.Enable()
      queue = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
      if queue == 'custom':
        self.custom_queue.Enable()
    else:
      self.queues.Disable()
      self.custom_queue.Disable()

  def onQueue(self, e):
    choice = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
    if choice == 'custom':
      self.custom_queue.Enable()
    else:
      self.custom_queue.Disable()

  def onMonitor(self, e):
    if self.chk_cont_mode.GetValue():
      self.chk_mm_timeout.Enable()
    else:
      self.chk_mm_timeout.Disable()
      self.chk_mm_timeout.SetValue(False)
      self.opt_timeout.Disable()
      self.opt_timeout.timeout.SetValue('')

  def onTimeout(self, e):
    if self.chk_mm_timeout.GetValue():
      self.opt_timeout.Enable()
      self.opt_timeout.timeout.SetValue('30')
    else:
      self.opt_timeout.Disable()
      self.opt_timeout.timeout.SetValue('')

  def onOK(self, e):
    # Get continous mode settings
    self.monitor_mode = self.chk_cont_mode.GetValue()
    self.mm_timeout = self.chk_mm_timeout.GetValue()

    if self.opt_timeout.timeout.GetValue() == '':
      self.mm_timeout_len = 0
    else:
      self.mm_timeout_len = int(self.opt_timeout.timeout.GetValue())

    self.method = self.mp_methods.ctr.GetString(self.mp_methods.ctr.GetSelection())
    if self.method == 'lsf':
      queue_selection = self.queues.ctr.GetString(self.queues.ctr.GetSelection())
      if queue_selection == 'custom':
        if self.custom_queue.cqueue.GetValue() == '':
          wx.MessageBox('Please choose or enter a queue', wx.OK)
        else:
          self.queue = self.custom_queue.cqueue.GetValue()
      else:
        self.queue = queue_selection
    else:
      self.queue = None

    n_bck = self.proc_backend.ctr.GetSelection()
    if n_bck == 0:
      backend = 'cctbx.xfel'
    elif n_bck == 1:
      backend = 'ha14'
    else:
      backend = 'cctbx.xfel'

    viewer = self.viewers.ctr.GetString(self.viewers.ctr.GetSelection())
    temp_folder = noneset(self.temp_folder.ctr.GetValue())

    if self.chk_random_sample.toggle.IsChecked():
      random_number = self.chk_random_sample.ctr.GetValue()
    else:
      random_number = 0

    if self.chk_image_range.toggle.IsChecked():
      img_range = self.chk_image_range.range.GetValue()
    else:
      img_range = None


    # test generation of PHIL settings
    prefs_text = '\n'.join([
      'mp {',
      'method = {}'.format(str(self.method)),
      'queue = {}'.format(str(self.queue)),
      '}',
      'gui {',
      '  image_viewer = {}'.format(viewer),
      '  monitor_mode = {}'.format(self.monitor_mode),
      '  monitor_mode_timeout = {}'.format(self.mm_timeout),
      '  monitor_mode_timeout_length = {}'.format(int(self.mm_timeout_len)),
      '}',
      'advanced {',
      '  processing_backend = {}'.format(backend),
      '  prime_prefix = {}'.format(self.prime_prefix.prefix.GetValue()),
      '  temporary_output_folder = {}'.format(temp_folder),
      '  image_range',
      '  {',
      '    flag_on = {}'.format(self.chk_image_range.toggle.GetValue()),
      '    range = {}'.format(img_range),
      '  }',
      '  random_sample',
      '  {',
      '    flag_on = {}'.format(self.chk_random_sample.toggle.GetValue()),
      '    number = {}'.format(random_number),
      '  }',
      '}'
    ])
    self.prefs_phil = ip.parse(prefs_text)
    e.Skip()


class ImportWindow(BaseDialog):
  # Import window - image import, modification and triage

  def __init__(self, parent, phil,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)

    self.params = phil.extract()
    self.import_phil = None

    conv_box = wx.StaticBox(self, label='Image Conversion Options')
    conv_box_sizer = wx.StaticBoxSizer(conv_box, wx.VERTICAL)

    # Image conversion options
    self.conv_rename = ct.ChoiceCtrl(self,
                                     label='Rename pickle:',
                                     choices=['keep_file_structure',
                                              'auto_filename',
                                              'custom_filename'],
                                     custom_choices=['custom_filename'],
                                     ctrl_size=(150, -1))
    conv_box_sizer.Add(self.conv_rename, flag=wx.ALL | wx.EXPAND, border=10)

    self.mod_square = ct.ChoiceCtrl(self,
                                    label='Image Modification:',
                                    choices=['no_modification', 'crop', 'pad'],
                                    ctrl_size=(150, -1))
    conv_box_sizer.Add(self.mod_square, 1, flag=wx.ALL | wx.EXPAND, border=10)

    self.flip_beamxy = wx.CheckBox(self, label='Flip BeamXY for crop or pad')
    conv_box_sizer.Add(self.flip_beamxy, flag=wx.ALL | wx.EXPAND, border=10)
    self.flip_beamxy.SetValue(False)

    self.mod_beamstop = ct.OptionCtrl(self,
                                      items=[('threshold', 0.0)],
                                      label='Beamstop shadow threshold',
                                      label_size=(200, -1),
                                      ctrl_size=(100, -1))
    conv_box_sizer.Add(self.mod_beamstop, flag=wx.ALL | wx.EXPAND, border=10)

    self.mod_detZ = ct.OptionCtrl(self,
                                  items=[('detZ', 0.0)],
                                  label='Detector distance (mm):',
                                  label_size=(200, -1),
                                  ctrl_size=(100, -1))
    conv_box_sizer.Add(self.mod_detZ, flag=wx.ALL | wx.EXPAND, border=10)

    self.mod_beamXY = ct.OptionCtrl(self,
                                    items=[('X', 0), ('Y', 0)],
                                    label='Direct beam XY (pixels)',
                                    label_size=(200, -1),
                                    ctrl_size=(100, -1))
    conv_box_sizer.Add(self.mod_beamXY, flag=wx.ALL | wx.EXPAND, border=10)

    self.mod_mask = ct.InputCtrl(self,
                                 label='Mask',
                                 label_size=wx.DefaultSize,
                                 buttons=True)
    conv_box_sizer.Add(self.mod_mask, 1, flag=wx.ALL | wx.EXPAND, border=10)

    # Image triage options
    trg_box = wx.StaticBox(self, label='Diffraction Triage Options')
    trg_box_sizer = wx.StaticBoxSizer(trg_box, wx.VERTICAL)

    self.img_triage = ct.ChoiceCtrl(self,
                                    label='Diffraction triage:',
                                    choices=['no_triage', 'simple',
                                             'grid_search'],
                                    ctrl_size=(150, -1))
    trg_box_sizer.Add(self.img_triage, flag=wx.ALL | wx.EXPAND, border=10)

    self.min_bragg_peaks = ct.OptionCtrl(self,
                                         items=[('n_bragg', 10)],
                                         label='Minimum Bragg peaks:',
                                         label_size=(200, -1),
                                         ctrl_size=(100, -1))
    trg_box_sizer.Add(self.min_bragg_peaks, flag=wx.ALL | wx.EXPAND, border=10)

    self.triage_spot_height = ct.OptionCtrl(self,
                                            items=[('min', ''), ('max', '')],
                                            label='Triage spot height:',
                                            label_size=(200, -1),
                                            ctrl_size=(50, -1))
    trg_box_sizer.Add(self.triage_spot_height, flag=wx.ALL | wx.EXPAND,
                      border=10)

    self.triage_spot_area = ct.OptionCtrl(self,
                                            items=[('min', ''), ('max', '')],
                                            label='Triage spot area:',
                                            label_size=(200, -1),
                                            ctrl_size=(50, -1))
    trg_box_sizer.Add(self.triage_spot_area, flag=wx.ALL | wx.EXPAND, border=10)

    self.triage_step_size = ct.OptionCtrl(self,
                                          items=[('step', '')],
                                          label='Grid search step:',
                                          label_size=(200, -1),
                                          ctrl_size=(50, -1))
    trg_box_sizer.Add(self.triage_step_size, flag=wx.ALL | wx.EXPAND, border=10)


    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    # Add all to main sizer
    self.main_sizer.Add(conv_box_sizer, flag=wx.ALL | wx.EXPAND, border=15)
    self.main_sizer.Add(trg_box_sizer, flag=wx.ALL | wx.EXPAND, border=15)
    self.main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.Bind(wx.EVT_CHOICE, self.onImageModChoice, self.mod_square.ctr)
    self.Bind(wx.EVT_CHOICE, self.onTriageChoice, self.img_triage.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onMaskBrowse, self.mod_mask.btn_browse)
    self.Bind(wx.EVT_BUTTON, self.onViewMask, self.mod_mask.btn_mag)

    self.read_phil()

  def onImageModChoice(self, e):
    selection = self.mod_square.ctr.GetSelection()
    if selection > 0:
      self.flip_beamxy.Enable()
    else:
      self.flip_beamxy.Disable()

  def onTriageChoice(self, e):
    selection = self.img_triage.ctr.GetString(self.img_triage.ctr.GetSelection())
    self.triage_choice(selection=selection)

  def onMaskBrowse(self, e):
    dlg = wx.FileDialog(
      self, message="Select mask file",
      defaultDir=os.curdir,
      defaultFile="*.pickle",
      wildcard="*.pickle",
      style=wx.OPEN | wx.CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]
      self.mod_mask.ctr.SetValue(filepath)

  def onViewMask(self, e):
    import iota.components.iota_threads as thr
    filepath = self.mod_mask.ctr.GetValue()
    if os.path.isfile(filepath):
      viewer = thr.ImageViewerThread(self,
                                     viewer=self.params.gui.image_viewer,
                                     file_string=filepath)
      viewer.start()

  def triage_choice(self, selection):

    if str(self.params.cctbx_ha14.image_triage.type).lower() == 'none':
      selection = 'no_triage'
    self.img_triage.ctr.SetSelection(self.img_triage.ctr.FindString(selection))
    if selection.lower() == 'no_triage':
      self.min_bragg_peaks.Disable()
      self.triage_spot_area.Disable()
      self.triage_spot_height.Disable()
      self.triage_step_size.Disable()
    elif selection.lower() == 'simple':
      self.min_bragg_peaks.Enable()
      self.min_bragg_peaks.n_bragg.SetValue(
        str(self.params.cctbx_ha14.image_triage.min_Bragg_peaks))
      self.triage_spot_area.Disable()
      self.triage_spot_height.Disable()
      self.triage_step_size.Disable()
    elif selection.lower() == 'grid_search':
      self.min_bragg_peaks.Enable()
      self.min_bragg_peaks.n_bragg.SetValue(
        str(self.params.cctbx_ha14.image_triage.min_Bragg_peaks))
      self.triage_spot_area.Enable()
      self.triage_spot_area.min.SetValue(
        str(self.params.cctbx_ha14.image_triage.grid_search.area_min))
      self.triage_spot_area.max.SetValue(
        str(self.params.cctbx_ha14.image_triage.grid_search.area_max))
      self.triage_spot_height.Enable()
      self.triage_spot_height.min.SetValue(
        str(self.params.cctbx_ha14.image_triage.grid_search.height_min))
      self.triage_spot_height.max.SetValue(
        str(self.params.cctbx_ha14.image_triage.grid_search.height_max))
      self.triage_step_size.Enable()
      self.triage_step_size.step.SetValue(
        str(self.params.cctbx_ha14.image_triage.grid_search.step_size))

  def read_phil(self):
    """ TODO: make PHIL reading more automated! """

    # Rename pickle prefix
    conv_prefix = str(self.params.cctbx_ha14.image_conversion.rename_pickle).lower()
    custom_filename = str(self.params.cctbx_ha14.image_conversion.rename_pickle_prefix).lower()

    if conv_prefix in self.conv_rename.choices:
      idx = self.conv_rename.ctr.FindString(conv_prefix)
      self.conv_rename.ctr.SetSelection(idx)
    else:
      self.conv_rename.ctr.SetSelection(0)
    self.conv_rename.set_choice(custom=custom_filename)

    # Image modification
    square_mode = self.params.cctbx_ha14.image_conversion.square_mode
    if square_mode in self.mod_square.choices:
      idx = self.mod_square.ctr.FindString(square_mode)
    else:
      idx = 1
    self.mod_square.ctr.SetSelection(idx)
    if idx > 0:
      self.flip_beamxy.Enable()
    else:
      self.flip_beamxy.Disable()

    self.flip_beamxy.SetValue(self.params.advanced.flip_beamXY)

    if str(self.params.image_import.mask).lower() == 'none':
      self.mod_mask.ctr.SetValue('')
    else:
      self.mod_mask.ctr.SetValue(str(self.params.image_import.mask))
    self.mod_beamstop.threshold.SetValue(
      str(self.params.image_import.beamstop))
    self.mod_detZ.detZ.SetValue(str(self.params.image_import.distance))
    self.mod_beamXY.X.SetValue(str(self.params.image_import.beam_center.x))
    self.mod_beamXY.Y.SetValue(str(self.params.image_import.beam_center.y))

    # Image modification
    self.triage_choice(self.params.cctbx_ha14.image_triage.type)

  def onOK(self, e):
    """ Accept changes and populate the PHIL scope """

    if self.mod_mask.ctr.GetValue() == '':
      maskpath = None
    else:
      maskpath = self.mod_mask.ctr.GetValue()
    if self.triage_spot_area.min.GetValue() == '':
      triage_spot_area_min = None
    else:
      triage_spot_area_min = self.triage_spot_area.min.GetValue()
    if self.triage_spot_area.max.GetValue() == '':
      triage_spot_area_max = None
    else:
      triage_spot_area_max = self.triage_spot_area.max.GetValue()
    if self.triage_spot_height.min.GetValue() == '':
      triage_spot_height_min = None
    else:
      triage_spot_height_min = self.triage_spot_height.min.GetValue()
    if self.triage_spot_height.max.GetValue() == '':
      triage_spot_height_max = None
    else:
      triage_spot_height_max = self.triage_spot_height.max.GetValue()
    if self.triage_step_size.step.GetValue() == '':
      triage_step_size = None
    else:
      triage_step_size = self.triage_step_size.step.GetValue()

    if self.conv_rename.custom.GetValue() == '':
      conv_pickle_prefix = None
    else:
      conv_pickle_prefix = self.conv_rename.custom.GetValue()

    if self.flip_beamxy.Enabled:
      if self.flip_beamxy.GetValue():
        flip_beamXY = True
      else:
        flip_beamXY = False
    else:
      flip_beamXY = False

    self.phil_text = '\n'.join([
    'image_conversion',
    '{',
    '  rename_pickle = {}'.format(self.conv_rename.ctr.GetString(
      self.conv_rename.ctr.GetSelection())),
    '  rename_pickle_prefix = {}'.format(conv_pickle_prefix),
    '  square_mode = {}'.format(self.mod_square.ctr.GetString(
      self.mod_square.ctr.GetSelection())),
    '  mask = {}'.format(maskpath),
    '  beamstop = {}'.format(self.mod_beamstop.threshold.GetValue()),
    '  distance = {}'.format(self.mod_detZ.detZ.GetValue()),
    '  beam_center',
    '  {',
    '    x = {}'.format(self.mod_beamXY.X.GetValue()),
    '    y = {}'.format(self.mod_beamXY.Y.GetValue()),
    '  }',
    '}',
    'image_triage',
    '{',
    '  type = {}'.format(self.img_triage.ctr.GetString(
      self.img_triage.ctr.GetSelection())),
    '    .type = choice',
    '  min_Bragg_peaks = {}'.format(self.min_bragg_peaks.n_bragg.GetValue()),
    '  grid_search',
    '  {',
    '    area_min = {}'.format(triage_spot_area_min),
    '    area_max = {}'.format(triage_spot_area_max),
    '    height_min = {}'.format(triage_spot_height_min),
    '    height_max = {}'.format(triage_spot_height_max),
    '    step_size = {}'.format(triage_step_size),
    '  }',
    '}',
    'cctbx_xfel {',
    '  min_Bragg_peaks = {}'.format(self.min_bragg_peaks.n_bragg.GetValue()),
    '}',
    'advanced',
    '{',
    '  flip_beamXY = {}'.format(flip_beamXY),
    '}'])
    self.import_phil = ip.parse(self.phil_text)
    e.Skip()


class OldBackendOptions(BaseBackendDialog):
  # CCTBX.XFEL options

  def __init__(self, parent,
               phil, target,
               *args, **kwargs):

    BaseBackendDialog.__init__(self, parent,
                               backend_name='HA14',
                               target=target,
                               phil=phil,
                               opt_size=(500, 500),
                               phil_size=(500, 500),
                               *args, **kwargs)

    self.params = phil.extract()
    self.proc_phil = None

    self.splitter.SplitVertically(self.options, self.phil_panel)

    # Target file input
    self.phil = ct.PHILBox(self.phil_panel, btn_pos='bottom',
                           ctr_size=(-1, 300))
    self.phil_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=5)

    # Grid search options
    # Type selection
    self.xtal_options = wx.Panel(self.options)
    xtal_box = wx.StaticBox(self.xtal_options, label='Processing Options')
    xtal_box_sizer = wx.StaticBoxSizer(xtal_box, wx.VERTICAL)
    self.xtal_options.SetSizer(xtal_box_sizer)

    self.res_limits = ct.OptionCtrl(self.xtal_options,
                                    label='Resolution Limits: ',
                                    label_size=(120, -1),
                                    items=[('lowres', '50.0'),
                                           ('hires', '1.5')],
                                    sub_labels=['low', 'high'],
                                    sub_label_justify=wx.ALIGN_RIGHT,
                                    ctrl_size=(100, -1))
    xtal_box_sizer.Add(self.res_limits, flag=wx.ALL, border=10)

    self.target_uc = ct.OptionCtrl(self.xtal_options,
                                   label='Target Unit Cell: ',
                                   label_size=(120, -1),
                                   items = [('cell', '')])
    xtal_box_sizer.Add(self.target_uc, flag=wx.ALL, border=10)

    lattice_types = ["None", "triclinic", "monoclinic", "orthorhombic",
                     "tetragonal", "hexagonal", "rhombohedral",  "cubic"]
    self.target_lattice = ct.ChoiceCtrl(self.xtal_options,
                                        label='Target Lattice Type:',
                                        label_size=(120, -1),
                                        choices=lattice_types,
                                        ctrl_size=(150, -1))
    xtal_box_sizer.Add(self.target_lattice, flag=wx.ALL, border=10)

    centering_types = ["None", "P - Primitive", "C - Base centered",
                       "I - Body centered", "R - Rhombohedral",
                       "F- Face centered"]
    self.target_centering = ct.ChoiceCtrl(self.xtal_options,
                                          label='Target Centering Type:',
                                          label_size=(120, -1),
                                          choices=centering_types,
                                          ctrl_size=(150, -1))
    xtal_box_sizer.Add(self.target_centering, flag=wx.ALL, border=10)

    self.gs_options = wx.Panel(self.options)
    gs_box = wx.StaticBox(self.gs_options, label='Grid Search Options')
    gs_box_sizer = wx.StaticBoxSizer(gs_box, wx.VERTICAL)

    self.gs_type = ct.ChoiceCtrl(self.gs_options,
                                 label='Grid search method:',
                                 choices=['no_grid_search',
                                          'brute_force',
                                          'smart'],
                                 ctrl_size=(150, -1),
                                 label_size=(200, -1))
    gs_box_sizer.Add(self.gs_type, flag=f.stack, border=10)

    self.signal_search = wx.CheckBox(self.gs_options, label='Signal height search')
    gs_box_sizer.Add(self.signal_search, flag=f.stack, border=10)

    self.spot_height = ct.OptionCtrl(self.gs_options,
                                     items=[('median', 0), ('range', 0)],
                                     label='Minimum spot height:',
                                     sub_labels=['', '+/-'],
                                     label_size=(200, -1),
                                     ctrl_size=(50, -1))
    gs_box_sizer.Add(self.spot_height, flag=f.stack, border=10)

    self.spot_area = ct.OptionCtrl(self.gs_options,
                                   items=[('median', 0), ('range', 0)],
                                   label='Minimum spot area:',
                                   sub_labels=['', '+/-'],
                                   label_size=(200, -1),
                                   ctrl_size=(50, -1))
    gs_box_sizer.Add(self.spot_area, flag=wx.ALL, border=10)

    self.gs_options.SetSizer(gs_box_sizer)


    # Selection options
    self.sel_options = wx.Panel(self.options)
    sel_box = wx.StaticBox(self.sel_options,
                           label='Grid Search Selection Options')
    sel_box_sizer = wx.StaticBoxSizer(sel_box, wx.VERTICAL)

    # self.select_only = wx.CheckBox(self.sel_options, label="Select only")
    # sel_box_sizer.Add(self.select_only, flag=f.stack, border=10)
    #
    # self.img_objects_path = ct.InputCtrl(self.sel_options,
    #                                      label='Image objects:',
    #                                      label_size=(120, -1),
    #                                      buttons=True)
    # sel_box_sizer.Add(self.img_objects_path, 1, flag=f.expand, border=10)

    self.select_by = ct.ChoiceCtrl(self.sel_options,
                                   label='Select by:',
                                   choices=['epv', 'mosaicity'],
                                   label_size=(120, -1),
                                   ctrl_size=(150, -1))
    sel_box_sizer.Add(self.select_by, flag=f.stack, border=10)

    self.min_sigma = ct.OptionCtrl(self.sel_options,
                                   label='Minimum I / {}(I)'.format(u.sigma),
                                   items=[('sigma', '5.0')],
                                   label_size=(120, -1),
                                   ctrl_size=(80, -1))
    sel_box_sizer.Add(self.min_sigma, flag=wx.ALL, border=10)

    self.sel_options.SetSizer(sel_box_sizer)

    # Filters
    self.filter_options = wx.Panel(self.options)
    filter_box = wx.StaticBox(self.filter_options, label='Selection Filters')
    filter_box_sizer = wx.StaticBoxSizer(filter_box, wx.VERTICAL)

    self.filt_lattice = ct.OptionCtrl(self.filter_options,
                                      items=[('lattice', 'P4')],
                                      checkbox=True,
                                      checkbox_label='Bravais Lattice:',
                                      label_size=(160, -1),
                                      ctrl_size=(150, -1))
    filter_box_sizer.Add(self.filt_lattice, flag=f.stack, border=10)

    self.filt_uc = ct.OptionCtrl(self.filter_options,
                                 items=[('a', '79.4'), ('b', '79.4'),
                                        ('c', '38.1'), ('alpha', '90'),
                                        ('beta', '90'), ('gamma', '90'),
                                        ('tolerance', '0.05')],
                                 sub_labels=['a =', 'b =', 'c =',
                                             '{} ='.format(u.alpha),
                                             '{} ='.format(u.beta),
                                             '{} ='.format(u.gamma),
                                             '{} ='.format(u.sigma)],
                                 sub_label_justify=wx.ALIGN_RIGHT,
                                 label_size=(160, -1),
                                 grid=(3, 6),
                                 checkbox=True,
                                 checkbox_label='Unit Cell',
                                 ctrl_size=(50, -1))
    filter_box_sizer.Add(self.filt_uc, flag=f.stack, border=10)

    self.filt_res = ct.OptionCtrl(self.filter_options,
                                  items=[('res', '2.5')],
                                  checkbox=True,
                                  checkbox_label='Resolution:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_res, flag=f.stack, border=10)

    self.filt_ref = ct.OptionCtrl(self.filter_options,
                                  items=[('ref', '100')],
                                  checkbox=True,
                                  checkbox_label='Num. of reflections:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_ref, flag=wx.ALL, border=10)

    self.f_spacer = filter_box_sizer.AddSpacer((-1, 10))

    self.filter_options.SetSizer(filter_box_sizer)

    # Add everything to sizer
    self.options_sizer.Add(self.xtal_options, flag=wx.BOTTOM | wx.EXPAND,
                           border=10)
    self.options_sizer.Add(self.gs_options, flag=wx.BOTTOM | wx.EXPAND,
                           border=10)
    self.options_sizer.Add(self.sel_options, flag=wx.BOTTOM | wx.EXPAND,
                           border=10)
    self.options_sizer.Add(self.filter_options, flag=wx.BOTTOM | wx.EXPAND,
                           border=10)

    self.show_hide_script()
    self.read_param_phil()
    self.show_hide_advanced()

    self.Layout()
    self.options.SetupScrolling()

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onImportPHIL, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDefaultPHIL, self.phil.btn_default)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_CHOICE, self.onGSChoice, self.gs_type.ctr)
    self.Bind(wx.EVT_CHOICE, self.onLatChoice, self.target_lattice.ctr)
    self.Bind(wx.EVT_CHOICE, self.onCentChoice, self.target_centering.ctr)
    # self.Bind(wx.EVT_CHECKBOX, self.onSelCheck, self.select_only)
    self.Bind(wx.EVT_BUTTON, self.onHideScript, self.btn_hide_script)
    self.Bind(wx.EVT_CHOICE, self.onAdvanced, self.dlg_ctr.choice)

  def onAdvanced(self, e):
    mode = self.dlg_ctr.choice.GetSelection()
    if mode == 0:
      self.show_hide_advanced()
    else:
      self.show_hide_advanced(show=True)

  def show_hide_advanced(self, show=False):
    if show:
      self.sel_options.Show()
      self.signal_search.Show()
      self.filt_ref.Show()
      self.filt_res.Show()
      self.f_spacer.Show(False)
    else:
      self.sel_options.Hide()
      self.signal_search.Hide()
      self.filt_ref.Hide()
      self.filt_res.Hide()
      self.f_spacer.Show(True)

    self.options.Layout()
    self.options.SetupScrolling()

  def onHideScript(self, e):
    self.opt_size = self.options.GetSize()
    self.phil_size = self.phil_panel.GetSize()
    self.sash_position = self.opt_size[0]
    self.show_hide_script(initialized=True)

  def onImportPHIL(self, e):
    phil_content = self.get_target_file()
    if phil_content is not None:
      self.phil.ctr.SetValue(phil_content)

  def onDefaultPHIL(self, e):
    self.write_default_phil()
    self.phil.ctr.SetValue(self.target_phil)

  # def onSelCheck(self, e):
  #   self.img_objects_path.Enable(self.select_only.GetValue())

  def onGSChoice(self, e):
    self.set_grid_search(self.gs_type.ctr.GetSelection())

  def onLatChoice(self, e):
    lat_choice = self.target_lattice.ctr.GetSelection()
    if lat_choice == 0:
      self.target_centering.ctr.SetSelection(0)
    else:
      if self.target_centering.ctr.GetSelection() == 0:
        self.target_centering.ctr.SetSelection(1)

  def onCentChoice(self, e):
    cent_choice = self.target_centering.ctr.GetSelection()
    if cent_choice == 0:
      self.target_lattice.ctr.SetSelection(0)
    else:
      if self.target_lattice.ctr.GetSelection() == 0:
        self.target_lattice.ctr.SetSelection(1)

  def set_grid_search(self, idx=1):
    self.gs_type.ctr.SetSelection(idx)

    if idx == 0:
      self.signal_search.Disable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Disable()
      self.spot_height.median.SetValue(str(self.params.cctbx_ha14.grid_search.height_median))
      self.spot_height.range.SetValue('0')
      self.spot_area.range.Disable()
      self.spot_area.median.SetValue(str(self.params.cctbx_ha14.grid_search.area_median))
      self.spot_area.range.SetValue('0')
    elif idx == 1:
      self.signal_search.Enable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Enable()
      self.spot_height.median.SetValue(
        str(self.params.cctbx_ha14.grid_search.height_median))
      self.spot_height.range.SetValue(
        str(self.params.cctbx_ha14.grid_search.height_range))
      self.spot_area.range.Enable()
      self.spot_area.median.SetValue(
        str(self.params.cctbx_ha14.grid_search.area_median))
      self.spot_area.range.SetValue(
        str(self.params.cctbx_ha14.grid_search.area_range))
    elif idx == 2:
      self.signal_search.Enable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Disable()
      self.spot_height.median.SetValue(
        str(self.params.cctbx_ha14.grid_search.height_median))
      self.spot_height.range.SetValue('1')
      self.spot_area.range.Disable()
      self.spot_area.median.SetValue(
        str(self.params.cctbx_ha14.grid_search.area_median))
      self.spot_area.range.SetValue('1')

  def read_param_phil(self):
    """ Reads parameters in the IOTA param PHIL and populates option controls """

    # LABELIT target file settings
    if self.target_phil is None:
      self.write_default_phil()
    self.phil.ctr.SetValue(self.target_phil)

    # Resolution limits
    # "Try/except" for backwards compatibility
    try:
      lowres = self.params.cctbx_ha14.resolution_limits.low
      hires = self.params.cctbx_ha14.resolution_limits.high
      self.res_limits.lowres.SetValue(str(lowres))
      self.res_limits.hires.SetValue(str(hires))
    except AttributeError:
      pass

    # Target options
    # "Try/except" for backwards compatibility
    try:
      t_uc = self.params.cctbx_ha14.target_unit_cell
      t_lat = self.params.cctbx_ha14.target_lattice_type
      l_idx = self.target_lattice.ctr.FindString(str(t_lat))
      t_ctype = self.params.cctbx_ha14.target_centering_type
      if t_ctype == 'P':
        c_idx = 1
      elif t_ctype == 'C':
        c_idx = 2
      elif t_ctype == 'I':
        c_idx = 3
      elif t_ctype == 'R':
        c_idx = 4
      elif t_ctype == 'F':
        c_idx = 5
      else:
        c_idx = 0
      if t_uc is not None:
        uc_str = [str(i) for i in t_uc.parameters()]
        self.target_uc.cell.SetValue(' '.join(uc_str))
      self.target_lattice.ctr.SetSelection(l_idx)
      self.target_centering.ctr.SetSelection(c_idx)
    except AttributeError:
      pass

    # Grid search options
    idx = self.gs_type.ctr.FindString(self.params.cctbx_ha14.grid_search.type)
    self.set_grid_search(idx=idx)
    self.signal_search.SetValue(self.params.cctbx_ha14.grid_search.sig_height_search)

    # # Selection options
    # self.select_only.SetValue(self.params.cctbx_ha14.selection.select_only.flag_on)
    # self.img_objects_path.Enable(self.select_only.GetValue())

    idx = self.select_by.ctr.FindString(self.params.cctbx_ha14.selection.select_by)
    self.select_by.ctr.SetSelection(idx)

    self.min_sigma.sigma.SetValue(str(self.params.cctbx_ha14.selection.min_sigma))

    # Selection filters
    if self.params.cctbx_ha14.selection.prefilter.flag_on:
      pg = self.params.cctbx_ha14.selection.prefilter.target_pointgroup
      ut = self.params.cctbx_ha14.selection.prefilter.target_uc_tolerance
      rs = self.params.cctbx_ha14.selection.prefilter.min_resolution
      rf = self.params.cctbx_ha14.selection.prefilter.min_reflections
      if self.params.cctbx_ha14.selection.prefilter.target_unit_cell is not None:
        try:
          uc = self.params.cctbx_ha14.selection.prefilter.target_unit_cell.parameters()
        except AttributeError:
          uc = None
      else:
        uc = None

      if str(pg).lower() != 'none':
        self.filt_lattice.toggle_boxes()
        self.filt_lattice.lattice.SetValue(str(pg))
      if str(uc).lower() != 'none':
        self.filt_uc.toggle_boxes()
        self.filt_uc.a.SetValue(str(uc[0]))
        self.filt_uc.b.SetValue(str(uc[1]))
        self.filt_uc.c.SetValue(str(uc[2]))
        self.filt_uc.alpha.SetValue(str(uc[3]))
        self.filt_uc.beta.SetValue(str(uc[4]))
        self.filt_uc.gamma.SetValue(str(uc[5]))
        self.filt_uc.tolerance.SetValue(str(ut))
      if str(rs).lower() != 'none':
        self.filt_res.toggle_boxes()
        self.filt_res.res.SetValue(str(rs))
      if str(rf).lower() != 'none':
        self.filt_ref.toggle_boxes()
        self.filt_ref.ref.SetValue(str(rf))

  def onOK(self, e):
    """ Output PHIL settings & save target file """

    # Read cctbx.xfel PHIL string
    if self.phil.ctr.GetValue() == '':
      trg_warning = wx.MessageDialog(None,
                                     message='No target parameters specified! Generate defaults?',
                                     caption='No Target Parameters',
                                     style=wx.YES_NO | wx.ICON_EXCLAMATION)
      if trg_warning.ShowModal() == wx.ID_YES:
        self.write_default_phil()
      trg_warning.Destroy()
      return
    else:
      self.target_phil = self.phil.ctr.GetValue()

    # Resolution limits
    # Resolution limits
    lowres = noneset(self.res_limits.lowres.GetValue())
    hires = noneset(self.res_limits.hires.GetValue())

    # Target crystal parameters
    t_uc = self.target_uc.cell.GetValue()
    t_lat = self.target_lattice.ctr.GetString(
      self.target_lattice.ctr.GetSelection())
    t_ctype = self.target_centering.ctr.GetSelection()
    if noneset(t_uc) != "None":
      target_uc = str(t_uc)
    else:
      target_uc = None

    if noneset(t_lat) == "None":
      target_lattice = None
    else:
      target_lattice = str(t_lat)

    if t_ctype == 1:
      target_centering = 'P'
    elif t_ctype == 2:
      target_centering = 'C'
    elif t_ctype == 3:
      target_centering = 'I'
    elif t_ctype == 4:
      target_centering = 'R'
    elif t_ctype == 5:
      target_centering = 'F'
    else:
      target_centering = None


    # # Grid search path (for select-only option)
    # grid_search_path = noneset(self.img_objects_path.ctr.GetValue())

    # Filter params
    filter_on = bool(self.filt_lattice.toggle.GetValue() +
                     self.filt_uc.toggle.GetValue() +
                     self.filt_ref.toggle.GetValue() +
                     self.filt_res.toggle.GetValue()
                     )
    if self.filt_lattice.toggle.GetValue():
      lattice = self.filt_lattice.lattice.GetValue()
    else:
      lattice = None

    if self.filt_uc.toggle.GetValue():
      uc = ', '.join([
        self.filt_uc.a.GetValue(),
        self.filt_uc.b.GetValue(),
        self.filt_uc.c.GetValue(),
        self.filt_uc.alpha.GetValue(),
        self.filt_uc.beta.GetValue(),
        (self.filt_uc.gamma.GetValue())
      ])
      tolerance = self.filt_uc.tolerance.GetValue()
    else:
      uc = None
      tolerance = None

    if self.filt_ref.toggle.GetValue():
      ref = noneset(self.filt_ref.ref.GetValue())
    else:
      ref = None

    if self.filt_res.toggle.GetValue():
      res = noneset(self.filt_res.res.GetValue())
    else:
      res = None

    # Populate IOTA settings
    proc_phil_text = '\n'.join([
      'cctbx_ha14',
      '{',
      '  resolution_limits',
      '  {',
      '    low = {}'.format(lowres),
      '    high = {}'.format(hires),
      '  }',
      '  target_unit_cell = {}'.format(target_uc),
      '  target_lattice_type = {}'.format(target_lattice),
      '  target_centering_type = {}'.format(target_centering),
      '  grid_search',
      '  {',
      '    type = {}'.format(self.gs_type.ctr.GetString(
        self.gs_type.ctr.GetSelection())),
      '    area_median = {}'.format(self.spot_area.median.GetValue()),
      '    area_range = {}'.format(self.spot_area.range.GetValue()),
      '    height_median = {}'.format(self.spot_height.median.GetValue()),
      '    height_range = {}'.format(self.spot_height.range.GetValue()),
      '    sig_height_search = {}'.format(self.signal_search.GetValue()),
      '    }',
      '  selection',
      '  {',
      # '    select_only',
      # '    {',
      # '      flag_on = {}'.format(self.select_only.GetValue()),
      # '      grid_search_path = {}'.format(grid_search_path),
      # '    }',
      '    min_sigma = {}'.format(self.min_sigma.sigma.GetValue()),
      '    select_by = {}'.format(self.select_by.ctr.GetString(
        self.select_by.ctr.GetSelection())),
      '    prefilter',
      '    {',
      '      flag_on = {}'.format(filter_on),
      '      target_pointgroup = {}'.format(lattice),
      '      target_unit_cell = {}'.format(uc),
      '      target_uc_tolerance = {}'.format(tolerance),
      '      min_reflections = {}'.format(ref),
      '      min_resolution = {}'.format(res),
      '    }',
      '  }',
      '}'
    ])

    self.proc_phil = ip.parse(proc_phil_text)
    e.Skip()


class BackendOptions(BaseBackendDialog):
  # DIALS options

  def __init__(self, parent,
               phil, target,
               *args, **kwargs):

    BaseBackendDialog.__init__(self, parent,
                               backend_name='CCTBX.XFEL',
                               phil=phil,
                               target=target,
                               phil_size=(500, 500),
                               opt_size=(500, 500),
                               *args, **kwargs)

    self.params = phil.extract()
    self.proc_phil = None

    self.splitter.SplitVertically(self.options, self.phil_panel)

    # Target file input
    self.phil = ct.PHILBox(self.phil_panel, btn_pos='bottom',
                           ctr_size=(-1, 300))
    self.phil_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=5)

    # Target parameters
    self.trg_options = wx.Panel(self.options)
    target_box = wx.StaticBox(self.trg_options, label='Target Parameters')
    target_box_sizer = wx.StaticBoxSizer(target_box, wx.VERTICAL)
    self.trg_options.SetSizer(target_box_sizer)

    self.target_sg = ct.OptionCtrl(self.trg_options,
                                   label='Target Space Group: ',
                                   label_size=(150, -1),
                                   ctrl_size=(150, -1),
                                   items = [('sg', '')])
    target_box_sizer.Add(self.target_sg, flag=f.stack, border=10)

    self.target_uc = ct.OptionCtrl(self.trg_options,
                                   label='Target Unit Cell: ',
                                   label_size=(150, -1),
                                   items = [('cell', '')])
    target_box_sizer.Add(self.target_uc, flag=f.stack, border=10)

    self.use_fft3d = wx.CheckBox(self.trg_options,
                                 label='Use FFT3D for indexing')
    self.use_fft3d.SetValue(False)
    target_box_sizer.Add(self.use_fft3d, flag=f.stack, border=10)

    self.sig_filter = ct.OptionCtrl(self.trg_options,
                                    items=[('sigma', '1.0')],
                                    checkbox=True,
                                    checkbox_state=True,
                                    checkbox_label='Significance Filter:',
                                    label_size=(150, -1),
                                    ctrl_size=(50, -1))
    self.sig_filter.toggle_boxes(False)
    target_box_sizer.Add(self.sig_filter, flag=wx.ALL, border=10)

    self.t_spacer = target_box_sizer.AddSpacer((-1, 10))

    # Optimization options
    self.opt_options = wx.Panel(self.options)
    optz_box = wx.StaticBox(self.opt_options, label='Optimization Options')
    optz_box_sizer = wx.StaticBoxSizer(optz_box, wx.VERTICAL)
    self.opt_options.SetSizer(optz_box_sizer)

    self.reindex = wx.CheckBox(self.opt_options,
                               label='Determine space group and reindex')
    self.reindex.SetValue(True)
    optz_box_sizer.Add(self.reindex, flag=wx.ALL, border=10)

    self.estimate_gain = wx.CheckBox(self.opt_options,
                                     label='Estimate gain for each image')
    self.estimate_gain.SetValue(True)
    optz_box_sizer.Add(self.estimate_gain, flag=wx.ALL, border=10)

    self.auto_threshold = wx.CheckBox(self.opt_options,
                                      label='Estimate threshold for each image')
    self.auto_threshold.SetValue(True)
    optz_box_sizer.Add(self.auto_threshold, flag=wx.ALL, border=10)

    # Filters
    self.filt_options = wx.Panel(self.options)
    filter_box = wx.StaticBox(self.filt_options, label='Filters')
    filter_box_sizer = wx.StaticBoxSizer(filter_box, wx.VERTICAL)
    self.filt_options.SetSizer(filter_box_sizer)

    self.filt_lattice = ct.OptionCtrl(self.filt_options,
                                      items=[('lattice', 'P4')],
                                      checkbox=True,
                                      checkbox_label='Bravais Lattice:',
                                      label_size=(160, -1),
                                      ctrl_size=(150, -1))
    filter_box_sizer.Add(self.filt_lattice, flag=f.stack, border=10)

    self.filt_uc = ct.OptionCtrl(self.filt_options,
                                 items=[('a', '79.4'), ('b', '79.4'),
                                        ('c', '38.1'), ('alpha', '90'),
                                        ('beta', '90'), ('gamma', '90'),
                                        ('tolerance', '0.05')],
                                 sub_labels=['a =', 'b =', 'c =',
                                             '{} ='.format(u.alpha),
                                             '{} ='.format(u.beta),
                                             '{} ='.format(u.gamma),
                                             '{} ='.format(u.sigma)],
                                 sub_label_justify=wx.ALIGN_RIGHT,
                                 label_size=(160, -1),
                                 grid=(3, 6),
                                 checkbox=True,
                                 checkbox_label='Unit Cell',
                                 ctrl_size=(50, -1))
    filter_box_sizer.Add(self.filt_uc, flag=f.stack, border=10)

    self.filt_res = ct.OptionCtrl(self.filt_options,
                                  items=[('res', '2.5')],
                                  checkbox=True,
                                  checkbox_label='Resolution:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_res, flag=f.stack, border=10)

    self.filt_ref = ct.OptionCtrl(self.filt_options,
                                  items=[('ref', '100')],
                                  checkbox=True,
                                  checkbox_label='Num. of reflections:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_ref, flag=wx.ALL, border=10)

    self.f_spacer = filter_box_sizer.AddSpacer((-1, 10))

    # Add all to sizers
    self.options_sizer.Add(self.trg_options, flag=wx.ALL | wx.EXPAND, border=10)
    self.options_sizer.Add(self.opt_options, flag=wx.ALL | wx.EXPAND, border=10)
    self.options_sizer.Add(self.filt_options, flag=wx.ALL | wx.EXPAND, border=10)

    self.show_hide_script()
    self.show_hide_advanced()
    self.Layout()
    self.options.SetupScrolling()
    self.read_param_phil()

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onImportPHIL, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDefaultPHIL, self.phil.btn_default)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onHideScript, self.btn_hide_script)
    self.Bind(wx.EVT_CHOICE, self.onAdvanced, self.dlg_ctr.choice)

  def onAdvanced(self, e):
    mode = self.dlg_ctr.choice.GetSelection()
    if mode == 0:
      self.show_hide_advanced()
    else:
      self.show_hide_advanced(show=True)

  def show_hide_advanced(self, show=False):
    if show:
      self.use_fft3d.Show()
      self.estimate_gain.Show()
      self.auto_threshold.Show()
      self.filt_ref.Show()
      self.filt_res.Show()
      self.t_spacer.Show(False)
      self.f_spacer.Show(False)
    else:
      self.use_fft3d.Hide()
      self.estimate_gain.Hide()
      self.auto_threshold.Hide()
      self.filt_ref.Hide()
      self.filt_res.Hide()
      self.t_spacer.Show(True)
      self.f_spacer.Show(True)

    self.options.Layout()
    self.options.SetupScrolling()


  def onHideScript(self, e):
    self.opt_size = self.options.GetSize()
    self.phil_size = self.phil_panel.GetSize()
    self.sash_position = self.opt_size[0]
    self.show_hide_script(initialized=True)

  def onImportPHIL(self, e):
    phil_content = self.get_target_file()
    if phil_content is not None:
      self.phil.ctr.SetValue(phil_content)

  def onDefaultPHIL(self, e):
    self.write_default_phil()
    self.phil.ctr.SetValue(self.target_phil)

  def read_param_phil(self):
    # DIALS target file settings
    if self.target_phil is None:
      self.write_default_phil()
    self.phil.ctr.SetValue(self.target_phil)

    # Target options
    # "Try/except" for backwards compatibility
    try:
      if self.params.cctbx_xfel.target_unit_cell is not None:
        uc_str = [str(i) for i in self.params.cctbx_xfel.target_unit_cell.parameters()]
        self.target_uc.cell.SetValue(' '.join(uc_str))
      if self.params.cctbx_xfel.target_space_group is not None:
        sg_info = str(self.params.cctbx_xfel.target_space_group).replace(' ', '')
        self.target_sg.sg.SetValue(sg_info)
      self.use_fft3d.SetValue(self.params.cctbx_xfel.use_fft3d)
      if self.params.cctbx_xfel.significance_filter.flag_on:
        sigma = self.params.cctbx_xfel.significance_filter.sigma
        self.sig_filter.toggle_boxes()
        self.sig_filter.sigma.SetValue(str(sigma))
    except AttributeError:
      pass

    # Optimization options
    self.reindex.SetValue(self.params.cctbx_xfel.determine_sg_and_reindex)
    self.estimate_gain.SetValue(self.params.advanced.estimate_gain)
    self.auto_threshold.SetValue(self.params.cctbx_xfel.auto_threshold)

   # Selection filters
    try:
      if self.params.cctbx_xfel.filter.flag_on:
        pg = self.params.cctbx_xfel.filter.target_pointgroup
        ut = self.params.cctbx_xfel.filter.target_uc_tolerance
        rs = self.params.cctbx_xfel.filter.min_resolution
        rf = self.params.cctbx_xfel.filter.min_reflections
        if self.params.cctbx_xfel.filter.target_unit_cell is not None:
          try:
            uc = self.params.cctbx_xfel.filter.target_unit_cell.parameters()
          except AttributeError:
            uc = None
        else:
          uc = None

        if str(pg).lower() != 'none':
          self.filt_lattice.toggle_boxes()
          self.filt_lattice.lattice.SetValue(str(pg))
        if str(uc).lower() != 'none':
          self.filt_uc.toggle_boxes()
          self.filt_uc.a.SetValue(str(uc[0]))
          self.filt_uc.b.SetValue(str(uc[1]))
          self.filt_uc.c.SetValue(str(uc[2]))
          self.filt_uc.alpha.SetValue(str(uc[3]))
          self.filt_uc.beta.SetValue(str(uc[4]))
          self.filt_uc.gamma.SetValue(str(uc[5]))
          self.filt_uc.tolerance.SetValue(str(ut))
        if str(rs).lower() != 'none':
          self.filt_res.toggle_boxes()
          self.filt_res.res.SetValue(str(rs))
        if str(rf).lower() != 'none':
          self.filt_ref.toggle_boxes()
          self.filt_ref.ref.SetValue(str(rf))
    except AttributeError:
      pass


  def onOK(self, e):
    """ Populate param PHIL file for DIALS options """

    if self.phil.ctr.GetValue() == '':
      trg_warning = wx.MessageDialog(None,
                                     message='No target parameters specified! Generate defaults?',
                                     caption='No Target Parameters',
                                     style=wx.YES_NO | wx.ICON_EXCLAMATION)
      if trg_warning.ShowModal() == wx.ID_YES:
        self.write_default_phil()
      trg_warning.Destroy()
      return
    else:
      self.target_phil = self.phil.ctr.GetValue()

    # Target params
    if self.target_sg.sg.GetValue() not in (None, ''):
      t_sg = self.target_sg.sg.GetValue()
    else:
      t_sg = None

    if self.target_uc.cell.GetValue() not in (None, ''):
      t_uc = self.target_uc.cell.GetValue()
    else:
      t_uc = None

    # Filter params
    filter_on = bool(self.filt_lattice.toggle.GetValue() +
                     self.filt_uc.toggle.GetValue() +
                     self.filt_ref.toggle.GetValue() +
                     self.filt_res.toggle.GetValue()
                     )
    if self.filt_lattice.toggle.GetValue():
      lattice = self.filt_lattice.lattice.GetValue()
    else:
      lattice = None

    if self.filt_uc.toggle.GetValue():
      uc = ', '.join([
                      self.filt_uc.a.GetValue(),
                      self.filt_uc.b.GetValue(),
                      self.filt_uc.c.GetValue(),
                      self.filt_uc.alpha.GetValue(),
                      self.filt_uc.beta.GetValue(),
                      (self.filt_uc.gamma.GetValue())
                      ])
      tolerance = self.filt_uc.tolerance.GetValue()
    else:
      uc = None
      tolerance = None

    if self.filt_ref.toggle.GetValue():
      ref = noneset(self.filt_ref.ref.GetValue())
    else:
      ref = None

    if self.filt_res.toggle.GetValue():
      res = noneset(self.filt_res.res.GetValue())
    else:
      res = None

    dials_phil_text = '\n'.join([
      'cctbx_xfel',
      '{',
      '  target_space_group = {}'.format(t_sg),
      '  target_unit_cell = {}'.format(t_uc),
      '  use_fft3d = {}'.format(str(self.use_fft3d.GetValue())),
      '  significance_filter',
      '  {',
      '    flag_on = {}'.format(self.sig_filter.toggle.GetValue()),
      '    sigma = {}'.format(noneset(self.sig_filter.sigma.GetValue())),
      '  }',
      '  determine_sg_and_reindex = {}'.format(self.reindex.GetValue()),
      '  auto_threshold = {}'.format(self.auto_threshold.GetValue()),
      '  filter',
      '    {',
      '      flag_on = {}'.format(filter_on),
      '      target_pointgroup = {}'.format(lattice),
      '      target_unit_cell = {}'.format(uc),
      '      target_uc_tolerance = {}'.format(tolerance),
      '      min_reflections = {}'.format(ref),
      '      min_resolution = {}'.format(res),
      '    }',
      '}',
      'advanced',
      '{',
      '  estimate_gain = {}'.format(self.estimate_gain.GetValue()),
      '}'
    ])

    self.proc_phil = ip.parse(dials_phil_text)
    e.Skip()


class AnalysisWindow(BaseDialog):
  # Import window - image import, modification and triage

  def __init__(self, parent, phil,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        *args, **kwargs)

    self.params = phil.extract()
    self.viz_phil = None

    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self, size=(-1, 250))
    options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(options_sizer)

    viz_box = wx.StaticBox(self.options, label='Analysis Options')
    viz_box_sizer = wx.StaticBoxSizer(viz_box, wx.VERTICAL)

    # Unit cell clustering options
    self.clustering = ct.OptionCtrl(self.options,
                                    items=[('threshold', '5000'),
                                           ('limit', '5'),
                                           ('n_images', '0')],
                                    sub_labels=['Threshold', 'Cluster limit',
                                                'No. images'],
                                    checkbox=True,
                                    checkbox_label='Unit Cell Clustering',
                                    #sub_label_vertical=wx.ALIGN_TOP,
                                    grid=(4, 2),
                                    label_size=(160, -1),
                                    ctrl_size=(100, -1))
    viz_box_sizer.Add(self.clustering, flag=f.stack, border=10)

    # Visualization options
    self.visualization = ct.ChoiceCtrl(self.options,
                                       choices=['no_visualization',
                                                'integration',
                                                'cv_vectors'],
                                       label='Visualization:',
                                       label_size=(160, -1),
                                       ctrl_size=(160, -1))
    viz_box_sizer.Add(self.visualization, flag=f.stack, border=10)

    self.proc_charts = wx.CheckBox(self.options,
                                   label='Output processing charts')
    viz_box_sizer.Add(self.proc_charts, flag=f.stack, border=10)

    self.summary_graphs = wx.CheckBox(self.options,
                                      label='Output run summary graphs')
    viz_box_sizer.Add(self.summary_graphs, flag=wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    # Add to sizer
    options_sizer.Add(viz_box_sizer, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.options, 1, flag=wx.EXPAND| wx.ALL, border=10)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Layout()
    self.options.SetupScrolling()

    self.read_param_phil()

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def read_param_phil(self):
    """ Read in parameters from IOTA parameter PHIL"""

    if self.params.analysis.run_clustering:
      self.clustering.toggle_boxes()
      self.clustering.threshold.SetValue(
        str(self.params.analysis.cluster_threshold))
      self.clustering.limit.SetValue(
        str(self.params.analysis.cluster_limit))
      self.clustering.n_images.SetValue(str(
        self.params.analysis.cluster_n_images))

    viz_idx = self.visualization.ctr.FindString(str(self.params.analysis.viz))
    if str(self.params.analysis.viz).lower() == 'none':
      viz_idx = 0
    self.visualization.ctr.SetSelection(viz_idx)

    self.proc_charts.SetValue(self.params.analysis.charts)
    self.summary_graphs.SetValue(self.params.analysis.summary_graphs)

  def onOK(self, e):
    """ Populate param PHIL and pass on to main """

    viz = self.visualization.ctr.GetString(self.visualization.ctr.GetSelection())

    analysis_phil_txt = '\n'.join([
      'analysis',
      '{',
      ' run_clustering = {}'.format(self.clustering.toggle.GetValue()),
      ' cluster_threshold = {}'.format(noneset(
          self.clustering.threshold.GetValue())),
      ' cluster_limit = {}'.format(noneset(
          self.clustering.limit.GetValue())),
      ' cluster_n_images = {}'.format(noneset(
          self.clustering.n_images.GetValue())),
      ' viz = {}'.format(viz),
      ' charts = {}'.format(self.proc_charts.GetValue()),
      ' summary_graphs = {}'.format(self.summary_graphs.GetValue()),
      '}'
    ])

    self.viz_phil = ip.parse(analysis_phil_txt)
    e.Skip()


class WatchModeTimeOut(wx.Dialog):
  def __init__(self, *args, **kwargs):
    super(WatchModeTimeOut, self).__init__(*args, **kwargs)

    self.timeout_length = None
    main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(main_sizer)

    option_box = wx.StaticBox(self, label='Timeout Options')
    option_sizer = wx.StaticBoxSizer(option_box, wx.VERTICAL)

    rb_sizer = wx.FlexGridSizer(4, 3, 10, 5)
    self.rb_30sec = wx.RadioButton(self, label='30 seconds', style=wx.RB_GROUP)
    self.rb_60sec = wx.RadioButton(self, label='60 seconds')
    self.rb_90sec = wx.RadioButton(self, label='90 seconds')
    self.rb_custom = wx.RadioButton(self, label='Custom: ')
    self.opt_custom = wx.TextCtrl(self, size=(100, -1))
    self.opt_custom.Disable()
    self.txt_custom = wx.StaticText(self, label='seconds')
    self.txt_custom.Disable()

    rb_sizer.AddMany([self.rb_30sec, (0, 0), (0, 0),
                     self.rb_60sec, (0, 0), (0, 0),
                     self.rb_90sec, (0, 0), (0, 0),
                     self.rb_custom, self.opt_custom, self.txt_custom])
    option_sizer.Add(rb_sizer, flag=wx.ALL | wx.EXPAND, border=10)
    main_sizer.Add(option_sizer, flag=wx.ALL | wx.EXPAND)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.rb_30sec.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb_60sec.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb_90sec.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb_custom.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onCustom(self, e):
    if self.rb_custom.GetValue():
      self.txt_custom.Enable()
      self.opt_custom.Enable()
      self.opt_custom.SetValue('120')
    else:
      self.txt_custom.Disable()
      self.opt_custom.Disable()
      self.opt_custom.SetValue('')

  def onOK(self, e):
    if self.rb_30sec.GetValue():
      self.timeout_length = 30
    elif self.rb_60sec.GetValue():
      self.timeout_length = 60
    elif self.rb_90sec.GetValue():
      self.timeout_length = 90
    elif self.rb_custom.GetValue():
      self.timeout_length = int(self.opt_custom.GetValue())
    e.Skip()


class DirView(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):


    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP

    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)

    self.dir_ctrl = wx.GenericDirCtrl(self, id=wx.ID_ANY,
                                      dir=os.path.abspath(os.curdir))
    self.main_sizer.Add(self.dir_ctrl, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)


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


    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)


class ViewerWarning(BaseDialog):
  def __init__(self, parent,
               img_list_length = None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    dlg_style = wx.CAPTION | wx.CLOSE_BOX | wx.RESIZE_BORDER | wx.STAY_ON_TOP
    BaseDialog.__init__(self, parent, style=dlg_style,
                        label_style=label_style,
                        content_style=content_style,
                        size=(400, 400),
                        *args, **kwargs)

    self.img_list_length = img_list_length
    self.no_images = 0

    self.opt_sizer = wx.FlexGridSizer(6, 3, 10, 10)

    self.rb1_img_view = wx.RadioButton(self, label='First 1 image',
                                       style=wx.RB_GROUP)
    self.rb2_img_view = wx.RadioButton(self, label='First 10 images')
    self.rb3_img_view = wx.RadioButton(self, label='First 50 images')
    self.rb4_img_view = wx.RadioButton(self, label='First 100 images')
    self.rb5_img_view = wx.RadioButton(self, label='All {} images'
                                      ''.format(self.img_list_length))
    self.rb_custom = wx.RadioButton(self, label='Other: ')
    self.opt_custom = wx.TextCtrl(self, size=(100, -1))
    self.opt_custom.Disable()
    self.txt_custom = wx.StaticText(self, label='images')
    self.txt_custom.Disable()

    self.opt_sizer.AddMany([self.rb1_img_view, (0, 0), (0, 0),
                            self.rb2_img_view, (0, 0), (0, 0),
                            self.rb3_img_view, (0, 0), (0, 0),
                            self.rb4_img_view, (0, 0), (0, 0),
                            self.rb5_img_view, (0, 0), (0, 0),
                            self.rb_custom, self.opt_custom, self.txt_custom])

    # Grey out irrelevant radio buttons
    if self.img_list_length < 100:
      self.rb4_img_view.Disable()
    if self.img_list_length < 50:
      self.rb3_img_view.Disable()
    if self.img_list_length < 10:
      self.rb3_img_view.Disable()

    self.main_sizer.Add(self.opt_sizer, flag=wx.ALL, border=10)

    # Dialog control
    self.main_sizer.Add(self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL),
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL, border=10)

    self.rb1_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb2_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb3_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb4_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb5_img_view.Bind(wx.EVT_RADIOBUTTON, self.onCustom)
    self.rb_custom.Bind(wx.EVT_RADIOBUTTON, self.onCustom)

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onCustom(self, e):
    if self.rb_custom.GetValue():
      self.txt_custom.Enable()
      self.opt_custom.Enable()
      if self.img_list_length < 25:
        value = str(self.img_list_length)
      else:
        value = '25'
      self.opt_custom.SetValue(value)
    else:
      self.txt_custom.Disable()
      self.opt_custom.Disable()
      self.opt_custom.SetValue('')

  def onOK(self, e):
    if self.rb1_img_view.GetValue():
      self.no_images = '1'
    elif self.rb2_img_view.GetValue():
      self.no_images = '1-10'
    elif self.rb3_img_view.GetValue():
      self.no_images = '1-50'
    elif self.rb4_img_view.GetValue():
      self.no_images = '1-100'
    elif self.rb5_img_view.GetValue():
      self.no_images = '1-{}'.format(self.img_list_length)
    elif self.rb_custom.GetValue():
      self.no_images = self.opt_custom.GetValue()
    self.EndModal(wx.ID_OK)


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
                        title='Recover Previous Run',
                        *args, **kwargs)


    self.pathlist = wx.ListCtrl(self, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    self.selected = None
    self.recovery_mode = 0

    self.pathlist.InsertColumn(0, "")
    self.pathlist.InsertColumn(1, "#")
    self.pathlist.InsertColumn(2, "Status")
    self.pathlist.InsertColumn(3, "Integration Path")

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
    pathlist = sorted(pathlist, key=lambda i:i.lower())
    for i in range(len(pathlist)):
      final_file = os.path.join(pathlist[i], 'integrated.lst')
      run_no = int(os.path.basename(pathlist[i]))
      if os.path.isfile(final_file):
        img_id = 0       # "finished" icon
        status = 'Finished'
      elif os.path.isfile(os.path.join(pathlist[i], '.abort.tmp')):
        img_id = 1       # "aborted" icon
        status = 'Aborted'
      else:
        img_id = 2       # "unknown" icon
        status = 'Unknown'
      idx = self.pathlist.InsertImageItem(i, img_id)
      self.pathlist.SetStringItem(idx, 1, str(run_no))
      self.pathlist.SetStringItem(idx, 2, status)
      self.pathlist.SetStringItem(idx, 3, pathlist[i])
      self.pathlist.SetColumnWidth(0, width=-1)
      self.pathlist.SetColumnWidth(1, width=-1)
      self.pathlist.SetColumnWidth(2, width=-1)
      self.pathlist.SetColumnWidth(3, width=-1)

      self.Fit()

  def onOK(self, e):
    for i in range(self.pathlist.GetItemCount()):
      if self.pathlist.IsSelected(idx=i):

        self.selected = [self.pathlist.GetItemText(i, col=2),
                         self.pathlist.GetItemText(i, col=3)]
        self.recovery_mode = self.dlg_ctr.choice.GetSelection()
    e.Skip()


class DIALSSpfDialog(BaseDialog):
  def __init__(self, parent,
               phil,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        *args, **kwargs)

    self.params = phil.extract()
    self.spf_phil = None

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    common_size = (200, -1)

    # Spotfinding settings
    self.options = wx.Panel(self, size=(500, -1))
    self.options_box = wx.StaticBox(self.options, label='Spotfinding Options')
    self.options_sizer = wx.StaticBoxSizer(self.options_box, wx.VERTICAL)
    self.options.SetSizer(self.options_sizer)
    self.sigma_background = ct.OptionCtrl(self.options,
                                          items=[('s_bkg', 6)],
                                          label='Sigma background',
                                          label_size=common_size,
                                          ctrl_size=(100, -1))
    self.options_sizer.Add(self.sigma_background, flag=wx.ALL, border=5)

    self.sigma_strong = ct.OptionCtrl(self.options,
                                      items=[('s_strong', 3)],
                                      label='Sigma strong',
                                      label_size=common_size,
                                      ctrl_size=(100, -1))
    self.options_sizer.Add(self.sigma_strong, flag=wx.ALL, border=5)

    self.global_threshold = ct.OptionCtrl(self.options,
                                          items=[('threshold', 0)],
                                          label='Global threshold',
                                          label_size=common_size,
                                          ctrl_size=(100, -1))
    self.options_sizer.Add(self.global_threshold, flag=wx.ALL, border=5)

    self.min_local = ct.OptionCtrl(self.options,
                                   items=[('min_local', 2)],
                                   label='Min. local',
                                   label_size=common_size,
                                   ctrl_size=(100, -1))
    self.options_sizer.Add(self.min_local, flag=wx.ALL, border=5)

    self.gain = ct.OptionCtrl(self.options,
                              items=[('gain', 1.0)],
                              label='Detector gain',
                              label_size=common_size,
                              ctrl_size=(100, -1))
    self.options_sizer.Add(self.gain, flag=wx.ALL, border=5)

    self.kernel_size = ct.OptionCtrl(self.options,
                                     items=[('kernel', '3 3')],
                                     label='Kernel size',
                                     label_size=common_size,
                                     ctrl_size=(100, -1))
    self.options_sizer.Add(self.kernel_size, flag=wx.ALL, border=5)

    self.mod_mask = ct.InputCtrl(self.options,
                                 label='Mask',
                                 label_size=wx.DefaultSize,
                                 buttons=True)
    self.options_sizer.Add(self.mod_mask, flag=wx.EXPAND | wx.ALL, border=5)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    self.main_sizer.Add(self.options, 1, wx.EXPAND)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    # Bindings:
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onMaskBrowse,
              self.mod_mask.btn_browse)
    self.Bind(wx.EVT_BUTTON, self.onViewMask,
              self.mod_mask.btn_mag)

    self.Fit()
    self.read_param_phil()

  def onMaskBrowse(self, e):
    dlg = wx.FileDialog(self,
                        message="Select mask file",
                        defaultDir=os.curdir,
                        defaultFile="*.pickle",
                        wildcard="*.pickle",
                        style=wx.OPEN | wx.CHANGE_DIR)
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]
      self.mod_mask.ctr.SetValue(filepath)

  def onViewMask(self, e):
    import iota.components.iota_threads as thr
    filepath = self.mod_mask.ctr.GetValue()
    if os.path.isfile(filepath):
      viewer = thr.ImageViewerThread(self,
                                     file_string=filepath,
                                     viewer='dials')
      viewer.start()

  def read_param_phil(self):

    # DIALS Spotfinder options
    self.sigma_background.s_bkg.SetValue(
      str(self.params.spotfinder.threshold.dispersion.sigma_background))
    self.sigma_strong.s_strong.SetValue(
      str(self.params.spotfinder.threshold.dispersion.sigma_strong))
    self.global_threshold.threshold.SetValue(
      str(self.params.spotfinder.threshold.dispersion.global_threshold))
    self.min_local.min_local.SetValue(
      str(self.params.spotfinder.threshold.dispersion.min_local))
    self.gain.gain.SetValue(str(self.params.spotfinder.threshold.dispersion.gain))
    self.kernel_size.kernel.SetValue('{} {}'.format(
      self.params.spotfinder.threshold.dispersion.kernel_size[0],
      self.params.spotfinder.threshold.dispersion.kernel_size[1]))
    if str(self.params.spotfinder.lookup.mask).lower() != 'none':
      self.mod_mask.ctr.SetValue(str(self.params.spotfinder.lookup.mask))
    else:
      self.mod_mask.ctr.SetValue('')

  def onOK(self, e):
    """ Populate param PHIL file for DIALS options """

    s_bkg = self.sigma_background.s_bkg.GetValue()
    s_str = self.sigma_strong.s_strong.GetValue()
    thresh = self.global_threshold.threshold.GetValue()
    min_local = self.min_local.min_local.GetValue()
    gain = self.gain.gain.GetValue()
    kernel = self.kernel_size.kernel.GetValue()
    mask = self.mod_mask.ctr.GetValue()
    if mask == '':
      mask = None

    phil_string = '\n'.join(['spotfinder {',
                             '  lookup.mask = {}'.format(mask),
                             '  threshold {',
                             '    dispersion {',
                             '      gain = {}'.format(gain),
                             '      kernel_size = {}'.format(kernel),
                             '      sigma_background = {}'.format(s_bkg),
                             '      sigma_strong = {}'.format(s_str),
                             '      min_local = {}'.format(min_local),
                             '      global_threshold = {}'.format(thresh),
                             '    }',
                             '  }',
                             '}'
                             ])
    self.spf_phil = ip.parse(phil_string)
    e.Skip()


class ClusterDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

    # Clustering parameters
    self.cluster_options = wx.Panel(self)
    cluster_box = wx.StaticBox(self.cluster_options, label='Cluster Parameters')
    cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_options.SetSizer(cluster_box_sizer)

    self.write_files = wx.CheckBox(self.cluster_options,
                                   label='Write Cluster Files')
    cluster_box_sizer.Add(self.write_files, flag=f.stack, border=10)

    self.cluster_n_images = ct.SpinCtrl(self.cluster_options,
                                        checkbox_label='No. images',
                                        checkbox=True, ctrl_size=(100, -1),
                                        ctrl_value=1000)
    cluster_box_sizer.Add(self.cluster_n_images, flag=f.expand, border=10)

    self.cluster_threshold = ct.SpinCtrl(self.cluster_options,
                                         label='Threshold: ',
                                         ctrl_size=(100, -1), ctrl_value=5000)
    cluster_box_sizer.Add(self.cluster_threshold, flag=f.expand, border=10)

    self.cluster_limit = ct.SpinCtrl(self.cluster_options,
                                     label='Minimum cluster size: ',
                                     ctrl_size=(100, -1), ctrl_value=10)
    cluster_box_sizer.Add(self.cluster_limit, flag=wx.EXPAND | wx.ALL,
                          border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    self.main_sizer.Add(self.cluster_options, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    # Bindings:
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

    self.Fit()

  def onOK(self, e):
    e.Skip()
