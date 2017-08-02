from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 01/17/2017
Last Changed: 08/01/2017
Description : IOTA GUI Dialogs
'''

import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wxtbx import bitmaps

from iotbx import phil as ip

import iota.components.iota_controls as ct
from iota.components.iota_input import master_phil
from iota.components.iota_misc import UnicodeCharacters, WxFlags, noneset


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


# ---------------------------------------------------------------------------- #

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

class IOTAPreferences(BaseDialog):
  ''' Class for dialog that houses IOTA interface preferences, e.g.:
        - multiprocessing / queue settings
        - monitor mode settings
        - miscellaneous interface-only settings '''

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
    self.method = self.params.mp_method
    self.queue = self.params.mp_queue
    self.monitor_mode = self.params.advanced.monitor_mode
    self.mm_timeout = self.params.advanced.monitor_mode_timeout
    self.mm_timeout_len = self.params.advanced.monitor_mode_timeout_length
    self.random_subset = self.params.advanced.random_sample.flag_on
    self.random_subset_number = self.params.advanced.random_sample.number

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

    self.custom_queue = ct.OptionCtrl(self,
                                      items=[('cqueue', '')],
                                      label='Custom Queue:',
                                      label_size=(120, -1),
                                      label_style='normal',
                                      ctrl_size=(150, -1))
    # self.custom_queue.Disable()
    queue_sizer.Add(self.custom_queue, flag=wx.LEFT | wx.RIGHT | wx.BOTTOM,
                    border=10)

    # Advanced Preferences
    adv_box = wx.StaticBox(self, label='Advanced Preferences')
    adv_sizer = wx.StaticBoxSizer(adv_box, wx.VERTICAL)

    # Monitor Mode preferences
    self.chk_cont_mode = wx.CheckBox(self, label='Process in Monitor Mode')
    self.chk_cont_mode.SetValue(False)
    self.chk_mm_timeout = wx.CheckBox(self, label='Monitor Mode Timeout')
    self.chk_mm_timeout.SetValue(False)
    self.opt_timeout = ct.OptionCtrl(self,
                                     items=[('timeout', '')],
                                     label='Timeout (sec):',
                                     label_size=(120, -1),
                                     label_style='normal',
                                     ctrl_size=(150, -1))
    self.chk_mm_timeout.Disable()
    self.opt_timeout.Disable()
    adv_sizer.Add(self.chk_cont_mode, flag=f.stack, border=10)
    adv_sizer.Add(self.chk_mm_timeout, flag=f.stack,  border=10)
    adv_sizer.Add(self.opt_timeout, flag=f.stack, border=10)

    # Random sample preferences
    self.chk_random_sample = wx.CheckBox(self, label='Process a random subset')
    self.random_number = ct.SpinCtrl(self,
                                     label='No. images in subset:',
                                     label_size=(160, -1),
                                     label_style='normal',
                                     ctrl_size=(80, -1),
                                     ctrl_min=0,
                                     ctrl_max=5000)
    self.random_number.Disable()
    adv_sizer.Add(self.chk_random_sample, flag=f.stack, border=10)
    adv_sizer.Add(self.random_number, flag=f.stack, border=10)

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
    self.Bind(wx.EVT_CHECKBOX, self.onRandom, self.chk_random_sample)
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

    # Set method to default value
    inp_method = self.mp_methods.ctr.FindString(noneset(self.method))
    if inp_method != wx.NOT_FOUND:
      self.mp_methods.ctr.SetSelection(inp_method)
    self.check_method()

    # Set Monitor Mode values
    if self.monitor_mode:
      self.chk_cont_mode.SetValue(True)
      self.chk_mm_timeout.Enable()
      if self.mm_timeout:
        self.chk_mm_timeout.SetValue(True)
        self.opt_timeout.Enable()
        self.opt_timeout.timeout.SetValue(str(self.mm_timeout_len))

    # Set random subset values
    if self.random_subset:
      self.chk_random_sample.SetValue(True)
      self.random_number.Enable()
    else:
      self.chk_random_sample.SetValue(False)
      self.random_number.Disable()
    self.random_number.ctr.SetValue(self.random_subset_number)

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

  def onRandom(self, e):
    if self.chk_random_sample.GetValue():
      self.random_number.Enable()
    else:
      self.random_number.Disable()

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

    temp_folder = noneset(self.temp_folder.ctr.GetValue())

    # test generation of PHIL settings
    prefs_text = '\n'.join([
      'mp_method = {}'.format(str(self.method)),
      'mp_queue = {}'.format(str(self.queue)),
      'advanced {',
      '  monitor_mode = {}'.format(self.monitor_mode),
      '  monitor_mode_timeout = {}'.format(self.mm_timeout),
      '  monitor_mode_timeout_length = {}'.format(int(self.mm_timeout_len)),
      '  prime_prefix = {}'.format(self.prime_prefix.prefix.GetValue()),
      '  temporary_output_folder = {}'.format(temp_folder),
      '  random_sample',
      '  {',
      '    flag_on = {}'.format(self.chk_random_sample.GetValue()),
      '    number = {}'.format(self.random_number.ctr.GetValue()),
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

    self.mask_invert = wx.CheckBox(self,
                                   label='Invert boolean mask')
    self.mask_invert.SetValue(False)
    self.mask_invert.Disable()
    conv_box_sizer.Add(self.mask_invert, flag=wx.ALL, border=10)


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

    self.Bind(wx.EVT_CHOICE, self.onTriageChoice, self.img_triage.ctr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onMaskBrowse, self.mod_mask.btn_browse)
    self.Bind(wx.EVT_BUTTON, self.onViewMask, self.mod_mask.btn_mag)

    self.read_phil()

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
      self.mask_invert.Enable()

  def onViewMask(self, e):
    import iota.components.iota_threads as thr
    filepath = self.mod_mask.ctr.GetValue()
    if os.path.isfile(filepath):
      viewer = thr.ImageViewerThread(self,
                                     backend=self.params.advanced.integrate_with,
                                     file_string=filepath)
      viewer.start()

  def triage_choice(self, selection):

    if str(self.params.image_triage.type).lower() == 'none':
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
        str(self.params.image_triage.min_Bragg_peaks))
      self.triage_spot_area.Disable()
      self.triage_spot_height.Disable()
      self.triage_step_size.Disable()
    elif selection.lower() == 'grid_search':
      self.min_bragg_peaks.Enable()
      self.min_bragg_peaks.n_bragg.SetValue(
        str(self.params.image_triage.min_Bragg_peaks))
      self.triage_spot_area.Enable()
      self.triage_spot_area.min.SetValue(
        str(self.params.image_triage.grid_search.area_min))
      self.triage_spot_area.max.SetValue(
        str(self.params.image_triage.grid_search.area_max))
      self.triage_spot_height.Enable()
      self.triage_spot_height.min.SetValue(
        str(self.params.image_triage.grid_search.height_min))
      self.triage_spot_height.max.SetValue(
        str(self.params.image_triage.grid_search.height_max))
      self.triage_step_size.Enable()
      self.triage_step_size.step.SetValue(
        str(self.params.image_triage.grid_search.step_size))

  def read_phil(self):
    ''' TODO: make PHIL reading more automated! '''

    # Rename pickle prefix
    conv_prefix = str(self.params.image_conversion.rename_pickle).lower()
    custom_filename = str(self.params.image_conversion.rename_pickle_prefix).lower()

    if conv_prefix in self.conv_rename.choices:
      idx = self.conv_rename.ctr.FindString(conv_prefix)
      self.conv_rename.ctr.SetSelection(idx)
    else:
      self.conv_rename.ctr.SetSelection(0)
    self.conv_rename.set_choice(custom=custom_filename)

    # Image modification
    square_mode = self.params.image_conversion.square_mode
    if square_mode in self.mod_square.choices:
      idx = self.mod_square.ctr.FindString(square_mode)
    else:
      idx = 1
    self.mod_square.ctr.SetSelection(idx)
    if str(self.params.image_conversion.mask).lower() == 'none':
      self.mod_mask.ctr.SetValue('')
      self.mask_invert.Disable()
    else:
      self.mod_mask.ctr.SetValue(str(self.params.image_conversion.mask))
      self.mask_invert.Enable()
      self.mask_invert.SetValue(self.params.image_conversion.invert_boolean_mask)
    self.mod_beamstop.threshold.SetValue(
      str(self.params.image_conversion.beamstop))
    self.mod_detZ.detZ.SetValue(str(self.params.image_conversion.distance))
    self.mod_beamXY.X.SetValue(str(self.params.image_conversion.beam_center.x))
    self.mod_beamXY.Y.SetValue(str(self.params.image_conversion.beam_center.y))

    # Image modification
    self.triage_choice(self.params.image_triage.type)

  def onOK(self, e):
    ''' Accept changes and populate the PHIL scope '''

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

    self.phil_text = '\n'.join([
    'image_conversion',
    '{',
    '  rename_pickle = {}'.format(self.conv_rename.ctr.GetString(
      self.conv_rename.ctr.GetSelection())),
    '  rename_pickle_prefix = {}'.format(conv_pickle_prefix),
    '  square_mode = {}'.format(self.mod_square.ctr.GetString(
      self.mod_square.ctr.GetSelection())),
    '  mask = {}'.format(maskpath),
    '  invert_boolean_mask = {}'.format(self.mask_invert.GetValue()),
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
    '}'])
    self.import_phil = ip.parse(self.phil_text)
    e.Skip()


class CCTBXOptions(BaseDialog):
  # CCTBX.XFEL options

  def __init__(self, parent,
               phil, target,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 500),
                        *args, **kwargs)

    self.target_phil = target
    self.params = phil.extract()
    self.proc_phil = None

    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self, size=(-1, 300))
    options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(options_sizer)

    phil_box = wx.StaticBox(self, label='LABELIT Target Settings')
    phil_box_sizer = wx.StaticBoxSizer(phil_box, wx.VERTICAL)

    # Target file input
    self.phil = ct.PHILBox(self,
                             btn_import=True,
                             btn_import_label='Import PHIL',
                             btn_export=False,
                             btn_default=True,
                             btn_default_label='Default PHIL',
                             ctr_size=(-1, 300),
                             ctr_value='')
    phil_box_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Grid search options
    # Type selection
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

    self.select_only = wx.CheckBox(self.sel_options, label="Select only")
    sel_box_sizer.Add(self.select_only, flag=f.stack, border=10)

    self.img_objects_path = ct.InputCtrl(self.sel_options,
                                         label='Image objects:',
                                         label_size=(120, -1),
                                         buttons=True)
    sel_box_sizer.Add(self.img_objects_path, 1, flag=f.expand, border=10)

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

    self.filter_options.SetSizer(filter_box_sizer)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    # Add everything to sizer
    self.main_sizer.Add(phil_box_sizer, 1, flag=wx.ALL | wx.EXPAND, border=10)
    options_sizer.Add(self.gs_options, flag=wx.ALL | wx.EXPAND, border=10)
    options_sizer.Add(self.sel_options, flag=wx.ALL | wx.EXPAND, border=10)
    options_sizer.Add(self.filter_options, flag=wx.ALL | wx.EXPAND,
                      border=10)

    self.main_sizer.Add(self.options, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(dialog_box, flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onImportPHIL, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDefaultPHIL, self.phil.btn_default)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_CHOICE, self.onGSChoice, self.gs_type.ctr)
    self.Bind(wx.EVT_CHECKBOX, self.onSelCheck, self.select_only)

    self.Layout()
    self.options.SetupScrolling()

    self.read_param_phil()

  def onImportPHIL(self, e):
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
      self.phil.ctr.SetValue(phil_content)

  def onDefaultPHIL(self, e):
    self.write_default_phil()

  def write_default_phil(self):
    from iota.components.iota_input import write_defaults
    default_phil, _ = write_defaults(current_path=None,
                                     txt_out=None,
                                     method='cctbx',
                                     write_target_file=False,
                                     write_param_file=False)
    self.phil.ctr.SetValue('\n'.join(default_phil))

  def onSelCheck(self, e):
    self.img_objects_path.Enable(self.select_only.GetValue())

  def onGSChoice(self, e):
    self.set_grid_search(self.gs_type.ctr.GetSelection())

  def set_grid_search(self, idx=1):

    self.gs_type.ctr.SetSelection(idx)

    if idx == 0:
      self.signal_search.Disable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Disable()
      self.spot_height.median.SetValue(str(self.params.cctbx.grid_search.height_median))
      self.spot_height.range.SetValue('0')
      self.spot_area.range.Disable()
      self.spot_area.median.SetValue(str(self.params.cctbx.grid_search.area_median))
      self.spot_area.range.SetValue('0')
    elif idx == 1:
      self.signal_search.Enable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Enable()
      self.spot_height.median.SetValue(
        str(self.params.cctbx.grid_search.height_median))
      self.spot_height.range.SetValue(
        str(self.params.cctbx.grid_search.height_range))
      self.spot_area.range.Enable()
      self.spot_area.median.SetValue(
        str(self.params.cctbx.grid_search.area_median))
      self.spot_area.range.SetValue(
        str(self.params.cctbx.grid_search.area_range))
    elif idx == 2:
      self.signal_search.Enable()
      self.signal_search.SetValue(False)
      self.spot_height.range.Disable()
      self.spot_height.median.SetValue(
        str(self.params.cctbx.grid_search.height_median))
      self.spot_height.range.SetValue('1')
      self.spot_area.range.Disable()
      self.spot_area.median.SetValue(
        str(self.params.cctbx.grid_search.area_median))
      self.spot_area.range.SetValue('1')

  def read_param_phil(self):
    ''' Reads parameters in the IOTA param PHIL and populates option controls '''

    # LABELIT target file settings
    if self.target_phil is not None:
      self.phil.ctr.SetValue(self.target_phil)
    else:
      self.phil.ctr.SetValue('')

    # Grid search options
    idx = self.gs_type.ctr.FindString(self.params.cctbx.grid_search.type)
    self.set_grid_search(idx=idx)
    self.signal_search.SetValue(self.params.cctbx.grid_search.sig_height_search)

    # Selection options
    self.select_only.SetValue(self.params.cctbx.selection.select_only.flag_on)
    self.img_objects_path.Enable(self.select_only.GetValue())

    idx = self.select_by.ctr.FindString(self.params.cctbx.selection.select_by)
    self.select_by.ctr.SetSelection(idx)

    self.min_sigma.sigma.SetValue(str(self.params.cctbx.selection.min_sigma))

    # Selection filters
    if self.params.cctbx.selection.prefilter.flag_on:
      pg = self.params.cctbx.selection.prefilter.target_pointgroup
      ut = self.params.cctbx.selection.prefilter.target_uc_tolerance
      rs = self.params.cctbx.selection.prefilter.min_resolution
      rf = self.params.cctbx.selection.prefilter.min_reflections
      if self.params.cctbx.selection.prefilter.target_unit_cell is not None:
        try:
          uc = self.params.cctbx.selection.prefilter.target_unit_cell.parameters()
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
    ''' Output PHIL settings & save target file '''

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

    grid_search_path = noneset(self.img_objects_path.ctr.GetValue())

    filter_on = bool(self.filt_lattice.toggle.GetValue() +
                     self.filt_uc.toggle.GetValue() +
                     self.filt_ref.toggle.GetValue() +
                     self.filt_res.toggle.GetValue()
                     )
    lattice = noneset(self.filt_lattice.lattice.GetValue())
    uc = noneset(', '.join([noneset(self.filt_uc.a.GetValue()),
                            noneset(self.filt_uc.b.GetValue()),
                            noneset(self.filt_uc.c.GetValue()),
                            noneset(self.filt_uc.alpha.GetValue()),
                            noneset(self.filt_uc.beta.GetValue()),
                            noneset((self.filt_uc.gamma.GetValue()))
                            ]
                           )
                 )
    tolerance = noneset(self.filt_uc.tolerance.GetValue())
    ref = noneset(self.filt_ref.ref.GetValue())
    res = noneset(self.filt_res.res.GetValue())

    proc_phil_text = '\n'.join([
      'cctbx',
      '{',
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
      '    select_only',
      '    {',
      '      flag_on = {}'.format(self.select_only.GetValue()),
      '      grid_search_path = {}'.format(grid_search_path),
      '    }',
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


class DIALSOptions(BaseDialog):
  # DIALS options

  def __init__(self, parent,
               phil, target,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        *args, **kwargs)

    self.target_phil = target
    self.params = phil.extract()
    self.proc_phil = None


    # Create options panel (all objects should be called as self.options.object)
    self.options = ScrolledPanel(self, size=(200, 200))
    options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(options_sizer)

    phil_box = wx.StaticBox(self, label='DIALS Target Settings')
    phil_box_sizer = wx.StaticBoxSizer(phil_box, wx.VERTICAL)
    self.phil = ct.PHILBox(self,
                           btn_import=True,
                           btn_import_label='Import PHIL',
                           btn_export=False,
                           btn_default=True,
                           btn_default_label='Default PHIL',
                           ctr_size=(500, 300),
                           ctr_value='')
    phil_box_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=10)

    dials_box = wx.StaticBox(self.options, label='DIALS Options')
    dials_box_sizer = wx.StaticBoxSizer(dials_box, wx.VERTICAL)

    # DIALS options
    self.reindex = wx.CheckBox(self.options,
                               label='Determine space group and reindex')
    self.reindex.SetValue(True)
    dials_box_sizer.Add(self.reindex, flag=wx.ALL, border=10)

    self.estimate_gain = wx.CheckBox(self.options,
                                     label='Estimate gain for each image')
    self.estimate_gain.SetValue(True)
    dials_box_sizer.Add(self.estimate_gain, flag=wx.ALL, border=10)

    self.auto_threshold = wx.CheckBox(self.options,
                                      label='Estimate threshold for each image')
    self.auto_threshold.SetValue(True)
    dials_box_sizer.Add(self.auto_threshold, flag=wx.ALL, border=10)

    # Filters
    filter_box = wx.StaticBox(self.options, label='Filters')
    filter_box_sizer = wx.StaticBoxSizer(filter_box, wx.VERTICAL)

    self.filt_lattice = ct.OptionCtrl(self.options,
                                      items=[('lattice', 'P4')],
                                      checkbox=True,
                                      checkbox_label='Bravais Lattice:',
                                      label_size=(160, -1),
                                      ctrl_size=(150, -1))
    filter_box_sizer.Add(self.filt_lattice, flag=f.stack, border=10)

    self.filt_uc = ct.OptionCtrl(self.options,
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

    self.filt_res = ct.OptionCtrl(self.options,
                                  items=[('res', '2.5')],
                                  checkbox=True,
                                  checkbox_label='Resolution:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_res, flag=f.stack, border=10)

    self.filt_ref = ct.OptionCtrl(self.options,
                                  items=[('ref', '100')],
                                  checkbox=True,
                                  checkbox_label='Num. of reflections:',
                                  label_size=(160, -1),
                                  ctrl_size=(100, -1))
    filter_box_sizer.Add(self.filt_ref, flag=wx.ALL, border=10)


    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)

    # Add all to sizers
    options_sizer.Add(dials_box_sizer, flag=wx.ALL | wx.EXPAND, border=10)
    options_sizer.Add(filter_box_sizer, flag=wx.ALL | wx.EXPAND, border=10)
    self.main_sizer.Add(phil_box_sizer, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.options, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onImportPHIL, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDefaultPHIL, self.phil.btn_default)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

    self.Layout()
    self.options.SetupScrolling()

    self.read_param_phil()

  def onImportPHIL(self, e):
    dlg = wx.FileDialog(
      self, message="Select DIALS target file",
      defaultDir=os.curdir,
      defaultFile="*.phil",
      wildcard="*",
      style=wx.OPEN | wx.CHANGE_DIR
    )
    if dlg.ShowModal() == wx.ID_OK:
      filepath = dlg.GetPaths()[0]

      with open(filepath, 'r') as phil_file:
        phil_content = phil_file.read()
      self.phil.ctr.SetValue(phil_content)


  def onDefaultPHIL(self, e):
    self.write_default_phil()

  def write_default_phil(self):
    from iota.components.iota_input import write_defaults
    default_phil, _ = write_defaults(current_path=None,
                                     txt_out=None,
                                     method='dials',
                                     write_target_file=False,
                                     write_param_file=False)
    self.phil.ctr.SetValue('\n'.join(default_phil))

  def read_param_phil(self):

    # DIALS target file settings
    if self.target_phil is not None:
      self.phil.ctr.SetValue(self.target_phil)
    else:
      self.phil.ctr.SetValue('')

    # DIALS options
    self.reindex.SetValue(self.params.dials.determine_sg_and_reindex)
    self.estimate_gain.SetValue(self.params.advanced.estimate_gain)
    self.auto_threshold.SetValue(self.params.dials.auto_threshold)

   # Selection filters
    try:
      if self.params.dials.filter.flag_on:
        pg = self.params.dials.filter.target_pointgroup
        ut = self.params.dials.filter.target_uc_tolerance
        rs = self.params.dials.filter.min_resolution
        rf = self.params.dials.filter.min_reflections
        if self.params.dials.filter.target_unit_cell is not None:
          try:
            uc = self.params.dials.filter.target_unit_cell.parameters()
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
    ''' Populate param PHIL file for DIALS options '''

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

    filter_on = bool(self.filt_lattice.toggle.GetValue() +
                     self.filt_uc.toggle.GetValue() +
                     self.filt_ref.toggle.GetValue() +
                     self.filt_res.toggle.GetValue()
                     )
    lattice = noneset(self.filt_lattice.lattice.GetValue())
    uc = noneset(', '.join([noneset(self.filt_uc.a.GetValue()),
                            noneset(self.filt_uc.b.GetValue()),
                            noneset(self.filt_uc.c.GetValue()),
                            noneset(self.filt_uc.alpha.GetValue()),
                            noneset(self.filt_uc.beta.GetValue()),
                            noneset((self.filt_uc.gamma.GetValue()))
                            ]
                           )
                 )
    tolerance = noneset(self.filt_uc.tolerance.GetValue())
    ref = noneset(self.filt_ref.ref.GetValue())
    res = noneset(self.filt_res.res.GetValue())

    dials_phil_text = '\n'.join([
      'dials',
      '{',
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
    self.options = ScrolledPanel(self, size=(-1, 200))
    options_sizer = wx.BoxSizer(wx.VERTICAL)
    self.options.SetSizer(options_sizer)

    viz_box = wx.StaticBox(self.options, label='Analysis Options')
    viz_box_sizer = wx.StaticBoxSizer(viz_box, wx.VERTICAL)

    # Unit cell clustering options
    self.clustering = ct.OptionCtrl(self.options,
                                    items=[('threshold', '5000')],
                                    sub_labels=['Threshold'],
                                    checkbox=True,
                                    checkbox_label='Unit Cell Clustering',
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
    ''' Read in parameters from IOTA parameter PHIL'''

    if self.params.analysis.run_clustering:
      self.clustering.toggle_boxes(flag_on=True)
      self.clustering.threshold.SetValue(
        str(self.params.analysis.cluster_threshold))

    viz_idx = self.visualization.ctr.FindString(str(self.params.analysis.viz))
    if str(self.params.analysis.viz).lower() == 'none':
      viz_idx = 0
    self.visualization.ctr.SetSelection(viz_idx)

    self.proc_charts.SetValue(self.params.analysis.charts)
    self.summary_graphs.SetValue(self.params.analysis.summary_graphs)

  def onOK(self, e):
    ''' Populate param PHIL and pass on to main '''

    viz = self.visualization.ctr.GetString(self.visualization.ctr.GetSelection())

    analysis_phil_txt = '\n'.join([
      'analysis',
      '{',
      ' run_clustering = {}'.format(self.clustering.toggle.GetValue()),
      ' cluster_threshold = {}'.format(noneset(
          self.clustering.threshold.GetValue())),
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
    self.rb5_img_view = wx.RadioButton(self, label='All images')
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
                        size=(400, 400),
                        title='Recover Previous Run',
                        *args, **kwargs)


    self.pathlist = wx.ListCtrl(self, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    self.selected = None

    self.pathlist.InsertColumn(0, "")
    self.pathlist.InsertColumn(1, "Status")
    self.pathlist.InsertColumn(2, "Integration Path")

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
      final_file = os.path.join(pathlist[i], 'integrated.lst')
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
                                     backend='dials',
                                     file_string=filepath)
      viewer.start()

  def read_param_phil(self):

    # DIALS Spotfinder options
    self.sigma_background.s_bkg.SetValue(
      str(self.params.spotfinder.threshold.xds.sigma_background))
    self.sigma_strong.s_strong.SetValue(
      str(self.params.spotfinder.threshold.xds.sigma_strong))
    self.global_threshold.threshold.SetValue(
      str(self.params.spotfinder.threshold.xds.global_threshold))
    self.min_local.min_local.SetValue(
      str(self.params.spotfinder.threshold.xds.min_local))
    self.gain.gain.SetValue(str(self.params.spotfinder.threshold.xds.gain))
    self.kernel_size.kernel.SetValue('{} {}'.format(
      self.params.spotfinder.threshold.xds.kernel_size[0],
      self.params.spotfinder.threshold.xds.kernel_size[1]))
    if str(self.params.spotfinder.lookup.mask).lower() != 'none':
      self.mod_mask.ctr.SetValue(str(self.params.spotfinder.lookup.mask))
    else:
      self.mod_mask.ctr.SetValue('')

  def onOK(self, e):
    ''' Populate param PHIL file for DIALS options '''

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
                             '    xds {',
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
