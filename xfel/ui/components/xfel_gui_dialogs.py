from __future__ import absolute_import, division, print_function
from six.moves import range
import six

'''
Author      : Lyubimov, A.Y.
Created     : 06/04/2016
Last Changed: 06/04/2016
Description : XFEL UI Custom Dialogs
'''

import os
import wx
import wx.richtext as rt
from wx.lib.mixins.listctrl import TextEditMixin, getListCtrlSelection
from wx.lib.scrolledpanel import ScrolledPanel

import xfel.ui.components.xfel_gui_controls as gctr

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')

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


class EdListCtrl(wx.ListCtrl, TextEditMixin):
  ''' TextEditMixin allows any column to be edited. '''

  # ----------------------------------------------------------------------
  def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
               size=wx.DefaultSize, style=0):
    """Constructor"""
    wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
    TextEditMixin.__init__(self)
    self.curRow = -1


class BaseDialog(wx.Dialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, *args, **kwargs)

    self.envelope = wx.BoxSizer(wx.VERTICAL)
    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.envelope.Add(self.main_sizer, flag=wx.EXPAND | wx.ALL, border=5)
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


# --------------------------------- Dialogs ---------------------------------- #

class SettingsDialog(BaseDialog):
  ''' Initial settings for cctbx.xfel; accepts DB credentials (may be
  separate dialog), populates experiment name / tag, separate dialog for
  multiprocessing and other settings; starts sentinels on OK '''

  def __init__(self, parent, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, -1),
                        *args, **kwargs)

    self.params = params
    self.drop_tables = False

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Experiment tag and DB Credentials button
    self.db_cred = gctr.TextButtonCtrl(self,
                                       label='Experiment Tag',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button=True,
                                       big_button_label='DB Credentials...',
                                       big_button_size=(130, -1),
                                       value=self.params.experiment_tag if self.params.experiment_tag is not None else "")
    self.main_sizer.Add(self.db_cred,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Facility control
    self.facility_sizer = wx.BoxSizer(wx.HORIZONTAL)

    choices = ['LCLS', 'Standalone']
    lower_choices = [f.lower() for f in choices]
    self.facility = gctr.ChoiceCtrl(self,
                                    label='Facility',
                                    label_size=(150, -1),
                                    ctrl_size=(180, -1),
                                    label_style='bold',
                                    choices=choices)
    self.facility_sizer.Add(self.facility, flag=wx.EXPAND | wx.ALL, border=10)
    try:
      self.facility.ctr.SetSelection(lower_choices.index(params.facility.name))
    except ValueError:
      pass

    self.btn_facility_options = wx.Button(self, label='Options...')
    self.facility_sizer.Add(self.btn_facility_options, flag=wx.EXPAND | wx.ALL, border=10)

    self.main_sizer.Add(self.facility_sizer, flag=wx.EXPAND | wx.ALL)

    # Experiment name control
    self.experiment = gctr.TextButtonCtrl(self,
                                          label='Experiment',
                                          label_style='bold',
                                          label_size=(150, -1),
                                          big_button_size=(130, -1),
                                          value=self.params.facility.lcls.experiment)
    self.main_sizer.Add(self.experiment,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Output folder text control w/ Browse / magnifying glass button
    if self.params.output_folder == '':
      current_folder = os.path.abspath(os.curdir)
    else:
      current_folder = self.params.output_folder
    self.output = gctr.TextButtonCtrl(self,
                                      label='Output',
                                      label_style='bold',
                                      label_size=(150, -1),
                                      big_button=True,
                                      big_button_label='Browse...',
                                      big_button_size=(120, -1),
                                      value=current_folder)
    self.main_sizer.Add(self.output,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    #self.btn_mp = wx.Button(self, label='Multiprocessing...')
    self.btn_op = wx.Button(self, label='Advanced Settings...')
    self.btn_OK = wx.Button(self, label="OK", id=wx.ID_OK)
    self.btn_cancel = wx.Button(self, label="Cancel", id=wx.ID_CANCEL)

    button_sizer = wx.FlexGridSizer(1, 4, 0, 10)
    button_sizer.AddMany([#(self.btn_mp),
                          (self.btn_op),
                          (0,0),
                          (self.btn_OK),
                          (self.btn_cancel)])

    button_sizer.AddGrowableCol(1)
    self.main_sizer.Add(button_sizer,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)
    self.SetSizer(self.main_sizer)

    self.Bind(wx.EVT_BUTTON, self.onDBCredentialsButton, id=self.db_cred.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onAdvanced, id=self.btn_op.GetId())
    self.Bind(wx.EVT_BUTTON, self.onOK, id=self.btn_OK.GetId())
    self.Bind(wx.EVT_BUTTON, self.onBrowse, id=self.output.btn_big.GetId())
    self.Bind(wx.EVT_CHOICE, self.onFacilityChoice)
    self.Bind(wx.EVT_BUTTON, self.onFacilityOptions, id=self.btn_facility_options.GetId())

    self.setup_facility_options()

  def onBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose the input directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.output.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

  def onFacilityChoice(self, e):
    self.params.facility.name = self.facility.ctr.GetStringSelection().lower()
    self.setup_facility_options()

  def setup_facility_options(self):
    if self.params.facility.name == 'lcls':
      self.experiment.Enable()
    else:
      self.experiment.Disable()

  def onFacilityOptions(self, e):
    if self.params.facility.name == 'lcls':
      opts = LCLSFacilityOptions(self, self.params)
      opts.Fit()
      opts.Center()
      opts.ShowModal()
    elif self.params.facility.name == 'standalone':
      opts = StandaloneOptions(self, self.params)
      opts.Fit()
      opts.Center()
      opts.ShowModal()

  def onAdvanced(self, e):
    adv = AdvancedSettingsDialog(self, self.params)
    adv.Fit()
    adv.Center()

    adv.ShowModal()

  def onDBCredentialsButton(self, e):
    creds = DBCredentialsDialog(self, self.params)
    creds.Center()
    if (creds.ShowModal() == wx.ID_OK):
      self.params.db.host     = creds.db_host.ctr.GetValue()
      self.params.db.port     = int(creds.db_port.ctr.GetValue())
      self.params.db.name     = creds.db_name.ctr.GetValue()
      self.params.db.user     = creds.db_user.ctr.GetValue()
      self.params.db.password = creds.db_password.ctr.GetValue()
      if self.params.facility.name == 'lcls':
        self.params.facility.lcls.web.user     = creds.web_user.ctr.GetValue()
        self.params.facility.lcls.web.password = creds.web_password.ctr.GetValue()

      self.drop_tables = creds.chk_drop_tables.GetValue()


  def onOK(self, e):
    self.params.facility.name = self.facility.ctr.GetStringSelection().lower()
    self.params.experiment_tag = self.db_cred.ctr.GetValue()
    self.params.facility.lcls.experiment = self.experiment.ctr.GetValue()
    self.params.output_folder = self.output.ctr.GetValue()
    e.Skip()


class DBCredentialsDialog(BaseDialog):
  ''' DB credentials entry '''

  def __init__(self, parent, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.params = params
    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        style=wx.NO_BORDER,
                        *args, **kwargs)

    # Host name
    self.db_host = gctr.TextButtonCtrl(self,
                                       label='DB Host name',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button_size=(130, -1),
                                       value=params.db.host)
    self.main_sizer.Add(self.db_host, flag=wx.EXPAND | wx.ALL, border=10)

    # Host name
    self.db_port = gctr.TextButtonCtrl(self,
                                       label='DB Port number',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button_size=(130, -1),
                                       value=str(params.db.port))
    self.main_sizer.Add(self.db_port, flag=wx.EXPAND | wx.ALL, border=10)

    # Database name
    self.db_name = gctr.TextButtonCtrl(self,
                                       label='DB name',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button_size=(130, -1),
                                       value=params.db.name)
    self.main_sizer.Add(self.db_name, flag=wx.EXPAND | wx.ALL, border=10)

    # User name
    self.db_user = gctr.TextButtonCtrl(self,
                                       label='DB user name',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button_size=(130, -1),
                                       value=params.db.user)
    self.main_sizer.Add(self.db_user, flag=wx.EXPAND | wx.ALL, border=10)

    # Password
    self.db_password = gctr.TextButtonCtrl(self,
                                           label='DB Password',
                                           label_style='bold',
                                           label_size=(150, -1),
                                           text_style=wx.TE_PASSWORD,
                                           big_button_size=(130, -1),
                                           value=params.db.password)
    self.main_sizer.Add(self.db_password, flag=wx.EXPAND | wx.ALL, border=10)

    # Drop tables button
    self.chk_drop_tables = wx.CheckBox(self,
                                       label='Delete and regenerate all tables')
    self.main_sizer.Add(self.chk_drop_tables, flag=wx.ALL, border=10)

    if params.facility.name == 'lcls':
      self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
      # LCLS user name
      self.web_user = gctr.TextButtonCtrl(self,
                                         label='LCLS user name',
                                         label_style='bold',
                                         label_size=(150, -1),
                                         big_button_size=(130, -1),
                                         value=params.facility.lcls.web.user)
      self.main_sizer.Add(self.web_user, flag=wx.EXPAND | wx.ALL, border=10)

      # LCLS password
      self.web_password = gctr.TextButtonCtrl(self,
                                             label='LCLS Password',
                                             label_style='bold',
                                             label_size=(150, -1),
                                             text_style=wx.TE_PASSWORD,
                                             big_button_size=(130, -1),
                                             value=params.facility.lcls.web.password)
      self.main_sizer.Add(self.web_password, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_CHECKBOX, self.onDropTables, self.chk_drop_tables)

    self.Fit()
    self.SetTitle('Database Credentials')

  def onDropTables(self, e):
    if self.chk_drop_tables.GetValue():
      msg = wx.MessageDialog(self,
                             message='Are you sure?',
                             caption='Warning',
                             style=wx.YES_NO |  wx.ICON_EXCLAMATION)

      if (msg.ShowModal() == wx.ID_NO):
        self.chk_drop_tables.SetValue(False)
    e.Skip()

class LCLSFacilityOptions(BaseDialog):
  ''' Options settings specific to LCLS'''
  def __init__(self, parent, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.params = params
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    self.chk_use_ffb = wx.CheckBox(self,
                                   label='Use ffb (fast feedback) file system. Active experiment only, on hiprio or prio queues')
    self.chk_use_ffb.SetValue(params.facility.lcls.use_ffb)

    self.chk_dump_shots = wx.CheckBox(self,
                                      label='Dump all images to disk. Useful for tuning spotfinding and indexing parameters')
    self.chk_dump_shots.SetValue(params.facility.lcls.dump_shots)

    self.chk_enforce80 = wx.CheckBox(self,
                                     label='Require stream 80 (FEE spectrometer) before processing')
    self.chk_enforce80.SetValue(params.facility.lcls.web.enforce80)


    self.chk_enforce81 = wx.CheckBox(self,
                                     label='Require stream 81 (FEE spectrometer) before processing')
    self.chk_enforce81.SetValue(params.facility.lcls.web.enforce81)

    self.main_sizer.Add(self.chk_use_ffb, flag=wx.ALL, border=10)
    self.main_sizer.Add(self.chk_dump_shots, flag=wx.ALL, border=10)
    self.main_sizer.Add(self.chk_enforce80, flag=wx.ALL, border=10)
    self.main_sizer.Add(self.chk_enforce81, flag=wx.ALL, border=10)
    self.SetSizer(self.main_sizer)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.SetTitle('LCLS Settings')

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onOK(self, e):
    self.params.facility.lcls.use_ffb = bool(self.chk_use_ffb.GetValue())
    self.params.facility.lcls.dump_shots = bool(self.chk_dump_shots.GetValue())
    self.params.facility.lcls.web.enforce80 = bool(self.chk_enforce80.GetValue())
    self.params.facility.lcls.web.enforce81 = bool(self.chk_enforce81.GetValue())
    e.Skip()

class StandaloneOptions(BaseDialog):
  ''' Options settings specific to standalone GUI '''
  def __init__(self, parent, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.params = params
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    # Output folder text control w/ Browse / magnifying glass button
    if self.params.facility.standalone.data_dir is None:
      current_folder = os.path.abspath(os.curdir)
    else:
      current_folder = self.params.facility.standalone.data_dir
    self.data_dir = gctr.TextButtonCtrl(self,
                                        label='Folder to monitor',
                                        label_style='bold',
                                        label_size=(300, -1),
                                        big_button=True,
                                        big_button_label='Browse...',
                                        big_button_size=(120, -1),
                                        value=current_folder)
    self.main_sizer.Add(self.data_dir,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    self.SetSizer(self.main_sizer)

    # Raw image option
    self.monitor_for = gctr.RadioCtrl(self,
                                      label='Monitor for',
                                      label_style='bold',
                                      label_size=(-1, -1),
                                      direction='horizontal',
                                      items={'files':'files',
                                            'folders':'folders'})
    getattr(self.monitor_for, self.params.facility.standalone.monitor_for).SetValue(1)

    self.main_sizer.Add(self.monitor_for, flag=wx.EXPAND | wx.ALL, border=10)

    # File matching template control
    if self.params.facility.standalone.template is None:
      self.params.facility.standalone.template = ''
    self.template = gctr.TextButtonCtrl(self,
                                          label='File matching template (example *.h5)',
                                          label_style='bold',
                                          label_size=(300, -1),
                                          value=self.params.facility.standalone.template)
    self.main_sizer.Add(self.template,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Composite check
    self.chk_composite = wx.CheckBox(self,
                                   label='Files are composite (like HDF5, files are submitted as individual runs.\nOtherwise, groups of files are submitted as single runs)')
    self.chk_composite.SetValue(params.facility.standalone.composite_files)
    self.main_sizer.Add(self.chk_composite,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    self.SetTitle('Standalone settings')

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onBrowse, id=self.data_dir.btn_big.GetId())

  def onOK(self, e):
    self.params.facility.standalone.data_dir = self.data_dir.ctr.GetValue()
    if self.monitor_for.files.GetValue():
      self.params.facility.standalone.monitor_for = 'files'
    else:
      self.params.facility.standalone.monitor_for = 'folders'
    self.params.facility.standalone.template = self.template.ctr.GetValue()
    self.params.facility.standalone.composite_files = self.chk_composite.GetValue()
    e.Skip()

  def onBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose the input directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.data_dir.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

class AdvancedSettingsDialog(BaseDialog):
  ''' Advanced settings for the cctbx.xfel front end '''
  def __init__(self, parent, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.params = params
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    # Multiprocessing settings
    mp_box = wx.StaticBox(self, label='Multiprocessing Options')
    self.mp_sizer = wx.StaticBoxSizer(mp_box, wx.VERTICAL)

    choices = ['python', 'lsf', 'mpi', 'sge', 'pbs', 'custom']
    self.mp_option = gctr.ChoiceCtrl(self,
                                     label='Multiprocessing:',
                                     label_size=(200, -1),
                                     label_style='bold',
                                     choices=choices)
    self.mp_sizer.Add(self.mp_option, flag=wx.EXPAND | wx.ALL, border=10)
    try:
      self.mp_option.ctr.SetSelection(choices.index(params.mp.method))
    except ValueError:
      pass

    if params.facility.name == 'lcls':
      # Queue
      queues = ['psanaq', 'psanaq', 'psdebugq','psanaidleq', 'psnehhiprioq',
                'psnehprioq', 'psnehq', 'psfehhiprioq', 'psfehprioq', 'psfehq']
      self.queue = gctr.ChoiceCtrl(self,
                                   label='Queue:',
                                   label_size=(200, -1),
                                   label_style='bold',
                                   choices=queues)
      self.Bind(wx.EVT_CHOICE, self.onQueueChoice, self.queue.ctr)
      self.mp_sizer.Add(self.queue, flag=wx.EXPAND | wx.ALL, border=10)
      try:
        self.queue.ctr.SetSelection(queues.index(params.mp.queue))
      except ValueError:
        pass

      self.nproc = gctr.SpinCtrl(self,
                                 label='Number of processors:',
                                 label_size=(200, -1),
                                 label_style='normal',
                                 ctrl_size=(100, -1),
                                 ctrl_value='%d'%params.mp.nproc,
                                 ctrl_min=1,
                                 ctrl_max=1000)
      self.mp_sizer.Add(self.nproc, flag=wx.EXPAND | wx.ALL, border=10)
    else:
      # Queue
      self.queue = gctr.TextButtonCtrl(self,
                                       label='Queue:',
                                       label_style='bold',
                                       label_size=(200, -1),
                                       value=self.params.mp.queue \
                                             if params.mp.queue is not None else '')
      self.mp_sizer.Add(self.queue, flag=wx.EXPAND | wx.ALL, border=10)

      self.nproc = gctr.SpinCtrl(self,
                                 label='Total number of processors:',
                                 label_size=(240, -1),
                                 label_style='normal',
                                 ctrl_size=(100, -1),
                                 ctrl_value='%d'%params.mp.nproc,
                                 ctrl_min=1,
                                 ctrl_max=1000)
      self.mp_sizer.Add(self.nproc, flag=wx.EXPAND | wx.ALL, border=10)
      self.nproc_per_node = gctr.SpinCtrl(self,
                                          label='Number of processors per node:',
                                          label_size=(240, -1),
                                          label_style='normal',
                                          ctrl_size=(100, -1),
                                          ctrl_value='%d'%params.mp.nproc_per_node,
                                          ctrl_min=1,
                                          ctrl_max=1000)
      self.mp_sizer.Add(self.nproc_per_node, flag=wx.EXPAND | wx.ALL, border=10)

      self.env_script = gctr.TextButtonCtrl(self,
                                       label='Environment setup script:',
                                       label_style='bold',
                                       label_size=(200, -1),
                                       value=self.params.mp.env_script[0] \
                                             if len(params.mp.env_script) > 0 else '')
      self.mp_sizer.Add(self.env_script, flag=wx.EXPAND | wx.ALL, border=10)

    self.main_sizer.Add(self.mp_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Data analysis settings
    analysis_box = wx.StaticBox(self, label='Data Analysis Options')
    self.analysis_sizer = wx.StaticBoxSizer(analysis_box, wx.VERTICAL)

    # Processing back-ends
    self.dispatchers_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.back_ends = ['cctbx.xfel (XTC+CBF mode)', 'cctbx.xfel (XTC mode)', 'Ha14', 'Small cell', 'custom']
    self.dispatchers = ['cctbx.xfel.xtc_process', 'cctbx.xfel.process', 'cxi.xtc_process', 'cctbx.xfel.small_cell_process', 'custom']
    self.dispatcher_descriptions = [
      'Process the data according to Brewster 2018, using DIALS for indexing, refinement and integration, with stills-specific defaults. Converts XTC into CBF in memory and optionally provides dumping of CBFs',
      'Process the data according to Brewster 2018, using DIALS for indexing, refinement and integration, with stills-specific defaults. Reads XTC directly.',
      'Process the data according to Hattne 2014, using LABELIT for initial indexing and stills-specific refinement and integration code implemented in the package cctbx.rstbx.',
      'Process the data according to Brewster 2015, using small cell for initial indexing and using DIALS for refinement and integration, with stills-specific defaults.',
      'Provide a custom program. See authors for details.']

    self.back_end = gctr.ChoiceCtrl(self,
                                    label='Processing back end:',
                                    label_size=(180, -1),
                                    label_style='bold',
                                    ctrl_size=(220, -1),
                                    choices=self.back_ends)
    self.Bind(wx.EVT_CHOICE, self.onBackendChoice)
    self.dispatchers_sizer.Add(self.back_end, flag=wx.ALIGN_LEFT)

    self.custom_dispatcher = gctr.TextCtrl(self,
                                           ctrl_size=(300, -1),
                                           value="")
    self.dispatchers_sizer.Add(self.custom_dispatcher, flag=wx.EXPAND | wx.ALL)

    try:
      self.back_end.ctr.SetSelection(self.dispatchers.index(params.dispatcher))
      self.custom_dispatcher.Hide()
    except ValueError:
      self.back_end.ctr.SetSelection(len(self.dispatchers)-1)
      self.custom_dispatcher.ctr.SetValue(params.dispatcher)

    self.analysis_sizer.Add(self.dispatchers_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    self.dispatcher_help = wx.StaticText(self, label=self.dispatcher_descriptions[self.back_end.ctr.GetSelection()], size=(600,80))
    self.dispatcher_help.Wrap(600)
    self.analysis_sizer.Add(self.dispatcher_help, flag=wx.EXPAND | wx.ALL, border=10)

    img_types = ['corrected', 'raw']
    self.avg_img_type = gctr.ChoiceCtrl(self,
                                        label='Avg. Image Type:',
                                        label_size=(200, -1),
                                        label_style='bold',
                                        ctrl_size=(200, -1),
                                        choices=img_types)
    if params.average_raw_data:
      i = img_types.index('raw')
    else:
      i = img_types.index('corrected')
    self.avg_img_type.ctr.SetSelection(i)

    self.analysis_sizer.Add(self.avg_img_type, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.analysis_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.SetTitle('Advanced Settings')

    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onQueueChoice(self, e):
    queue = self.queue.ctr.GetString(self.queue.ctr.GetSelection())
    if 'neh' in queue or 'feh' in queue:
      self.nproc.ctr.SetValue(16)
      self.nproc.ctr.SetIncrement(16)
    elif 'psana' in queue or 'debug' in queue:
      self.nproc.ctr.SetValue(12)
      self.nproc.ctr.SetIncrement(12)
    else:
      self.nproc.ctr.SetValue(1)
      self.nproc.ctr.SetIncrement(1)

  def onBackendChoice(self, e):
    self.params.dispatcher = self.dispatchers[self.back_end.ctr.GetSelection()]
    self.dispatcher_help.SetLabel(self.dispatcher_descriptions[self.back_end.ctr.GetSelection()])
    self.dispatcher_help.Wrap(600)
    if self.params.dispatcher == 'custom':
      self.custom_dispatcher.Show()
      self.Layout()
    else:
      self.custom_dispatcher.Hide()
      self.Layout()

  def onOK(self, e):
    self.params.dispatcher = self.dispatchers[self.back_end.ctr.GetSelection()]
    if self.params.dispatcher == 'custom':
      self.params.dispatcher = self.custom_dispatcher.ctr.GetValue()
    self.params.mp.method = self.mp_option.ctr.GetStringSelection()
    if self.params.facility.name == 'lcls':
      self.params.mp.queue = self.queue.ctr.GetStringSelection()
    else:
      self.params.mp.queue = self.queue.ctr.GetValue()
      self.params.mp.nproc_per_node = int(self.nproc_per_node.ctr.GetValue())
      self.params.mp.env_script = [self.env_script.ctr.GetValue()]
    self.params.mp.nproc = int(self.nproc.ctr.GetValue())
    self.params.average_raw_data = self.avg_img_type.ctr.GetStringSelection() == 'raw'
    e.Skip()

class CalibrationDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               db=None,
               *args, **kwargs):

    self.parent = parent
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.db = db

    # Metrology version name
    self.version_name = gctr.TextButtonCtrl(self,
                                            label='Metrology Version Name:',
                                            label_size=(180, -1),
                                            label_style='bold',
                                            value='cspad')
    self.main_sizer.Add(self.version_name, flag=wx.EXPAND | wx.ALL, border=10)

    # Reflection type and number of images in subset
    self.top_sizer = wx.FlexGridSizer(1, 2, 0, 20)
    self.top_sizer.AddGrowableCol(1, 1)
    choices = ['indexed', 'reindexed strong', 'integrated']
    self.reflections = gctr.ChoiceCtrl(self,
                                       label='Reflections:',
                                       label_size=(100, -1),
                                       label_style='bold',
                                       choices=choices)
    self.n_subset = gctr.SpinCtrl(self,
                                  label='Images in subset:',
                                  label_size=(120, -1),
                                  label_style='normal',
                                  ctrl_value='1000',
                                  ctrl_max=10000,
                                  ctrl_min=10)
    self.top_sizer.Add(self.reflections)
    self.top_sizer.Add(self.n_subset, wx.ALIGN_RIGHT)

    #Trial & runs
    self.trial_sizer = wx.FlexGridSizer(1, 2, 0, 20)
    self.trial_sizer.AddGrowableCol(1, 1)
    trials = [str(i.trial) for i in self.db.get_all_trials()]
    self.trial_number = gctr.ChoiceCtrl(self,
                                        label='Trial:',
                                        label_size=(40, -1),
                                        label_style='normal',
                                        ctrl_size=(200, -1),
                                        choices=trials)
    self.trial_runs = gctr.CheckListCtrl(self,
                                         label='Runs:',
                                         label_size=(40, -1),
                                         label_style='normal',
                                         ctrl_size=(200, -1),
                                         choices=[])
    self.trial_sizer.Add(self.trial_number)
    self.trial_sizer.Add(self.trial_runs, wx.ALIGN_RIGHT)

    #Phil blob
    self.phil_text = rt.RichTextCtrl(self, size=(550, 300), style=wx.VSCROLL)
    self.phil_path = gctr.TwoButtonCtrl(self,
                                        label='PHIL path:',
                                        label_size=(80, -1),
                                        label_style='normal',
                                        button1=True,
                                        button1_label='Browse...',
                                        button2=True,
                                        button2_label='Default PHIL')

    self.main_sizer.Add(self.top_sizer, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.trial_sizer, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.phil_text, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.phil_path, flag=wx.EXPAND | wx.ALL, border=10)

    self.chk_split_dataset = wx.CheckBox(self,
                                         label='Split dataset into 2 halves (outputs statistics, double runtime, uses 2x number of images')
    self.chk_split_dataset.SetValue(True)
    self.main_sizer.Add(self.chk_split_dataset, flag=wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)
    # Bindings
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice,
              id=self.trial_number.ctr.GetId())
    self.Bind(wx.EVT_BUTTON, self.onBrowse, self.phil_path.button1)
    self.Bind(wx.EVT_BUTTON, self.onDefault, self.phil_path.button2)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    self.Bind(wx.EVT_BUTTON, self.onCancel, id=wx.ID_CANCEL)

    self.SetTitle('Calibration Settings')
    self.find_runs()

    self.frames_sentinel = None
    self.start_frames_sentinel()

  def start_frames_sentinel(self):
    if self.frames_sentinel is None:
      from xfel.ui.components.xfel_gui_init import FramesSentinel
      self.frames_sentinel = FramesSentinel(self)
    self.frames_sentinel.active = True
    self.frames_sentinel.start()

  def stop_frames_sentinel(self):
    if self.frames_sentinel is not None:
      self.frames_sentinel.active = False
      self.frames_sentinel.join()

  def onOK(self, e):
    from xfel.ui.db import get_run_path
    from xfel.util.mp import get_lsf_submit_command
    from xfel.ui import settings_dir
    from libtbx import easy_run
    import copy
    params = copy.deepcopy(self.parent.params)
    params.mp.nproc = 1

    version_str = self.version_name.ctr.GetValue()
    working_dir = os.path.join(params.output_folder, "metrology", version_str)
    if os.path.exists(working_dir):
      wx.MessageBox('Version name %s already used'%version_str, 'Warning!',
                    wx.ICON_EXCLAMATION)
      return

    self.stop_frames_sentinel()

    runs = [self.db.get_run(run_number=int(r)) for r in self.trial_runs.ctr.GetCheckedStrings()]
    run_ids = [r.id for r in runs]
    command = "cspad.cbf_metrology reflections=%s tag=%s split_dataset=%s n_subset=%d "% (
      self.reflections.ctr.GetStringSelection(), self.version_name.ctr.GetValue(),
      self.chk_split_dataset.GetValue(), self.n_subset.ctr.GetValue())

    phil_file = os.path.join(settings_dir, "cfgs", "%s.phil"%version_str)
    f = open(phil_file, 'w')
    f.write(self.phil_text.GetValue())
    f.close()
    command += phil_file + " "

    trial = self.db.get_trial(trial_number=int(self.trial_number.ctr.GetStringSelection()))
    runs_found = []
    run_paths = []
    for rungroup in trial.rungroups:
      for run in rungroup.runs:
        if run.id in runs_found:
          continue
        if run.id in run_ids:
          run_ids.pop(run_ids.index(run.id))
          runs_found.append(run.id)
          run_paths.append(os.path.join(get_run_path(params.output_folder, trial, rungroup, run), 'out'))
    assert len(run_ids) == 0

    command += " ".join(run_paths)

    submit_path = os.path.join(settings_dir, "%s.sh"%version_str)
    command = str(get_lsf_submit_command(command, submit_path, working_dir, params.mp)()) # comes back as unicode which throws off easy_run
    cwd = os.getcwd()
    os.makedirs(working_dir)
    os.chdir(working_dir)

    print("Submitting metrology refinement. Command:")
    print(command)
    try:
      results = easy_run.fully_buffered(command=command)
      results.show_stdout()
      results.show_stderr()
      results.raise_if_errors()
    except Exception as exc:
      if not "Warning: job being submitted without an AFS token." in str(exc):
        raise exc
    os.chdir(cwd)
    print("Output will be in", working_dir)

    e.Skip()

  def onTrialChoice(self, e):
    self.find_runs()

  def find_runs(self):
    self.trial_runs.ctr.Clear()
    trial = self.db.get_all_trials()[self.trial_number.ctr.GetSelection()]
    runs = [str(i.run) for i in trial.runs]
    self.trial_runs.ctr.InsertItems(items=runs, pos=0)

  def onBrowse(self, e):
    ''' Open dialog for selecting PHIL file '''
    load_dlg = wx.FileDialog(self,
                             message="Load PHIL file",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      target_file = load_dlg.GetPaths()[0]
      with open(target_file, 'r') as phil_file:
        phil_file_contents = phil_file.read()
      self.phil_text.SetValue(phil_file_contents)
      self.phil_path.ctr.SetValue(target_file)
    load_dlg.Destroy()

  def onDefault(self, e):
    # TODO: Generate default PHIL parameters
    pass

  def onCancel(self, e):
    self.stop_frames_sentinel()
    e.Skip()

class AveragingDialog(BaseDialog):
  def __init__(self, parent, run, params,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.run = run
    self.params = params

    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    # Raw image option
    self.raw_toggle = gctr.RadioCtrl(self,
                                     label='',
                                     label_style='normal',
                                     label_size=(-1, -1),
                                     direction='horizontal',
                                     items={'corrected':'corrected',
                                            'raw':'raw'})
    self.raw_toggle.corrected.SetValue(1)
    self.main_sizer.Add(self.raw_toggle, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onOK(self, e):
    from xfel.ui.components.averaging import AveragingCommand
    from libtbx import easy_run
    raw = self.raw_toggle.raw.GetValue() == 1
    average_command = AveragingCommand(self.run, self.params, raw)()
    print("executing", average_command)
    result = easy_run.fully_buffered(command=average_command)
    result.show_stdout()
    e.Skip()

class TrialTagSelectionDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               db=None,
               *args, **kwargs):
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)
    self.db = db
    self.parent = parent
    self.run_paths = []
    self.tags = []
    self.tag_names = []

    self.trials = db.get_all_trials()
    trial_numbers = ["%d"%trial.trial for trial in self.trials]
    self.trials_ctrl = gctr.ChoiceCtrl(self,
                                       label='Trial number:',
                                       label_size=(120, -1),
                                       label_style='bold',
                                       choices=trial_numbers)
    self.main_sizer.Add(self.trials_ctrl, flag=wx.EXPAND | wx.ALL, border=10)


    self.trial_tags = gctr.CheckListCtrl(self,
                                         label='Tags (optional):',
                                         label_size=(40, -1),
                                         label_style='normal',
                                         ctrl_size=(200, -1),
                                         choices=[])
    self.main_sizer.Add(self.trial_tags, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.SetTitle('Pick a trial and optionally some tags from that trial')

    self.Bind(wx.EVT_CHOICE, self.onTrialSelect, id=self.trials_ctrl.ctr.GetId())
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)
    wx.CallAfter(self.refresh_tags)

  def onTrialSelect(self, e):
    self.refresh_tags()

  def refresh_tags(self):
    self.trial_tags.ctr.Clear()
    trial = self.trials[self.trials_ctrl.ctr.GetSelection()]
    self.tags = []
    self.tag_names = []
    tag_ids = []
    for run in trial.runs:
      for tag in run.tags:
        if tag.id not in tag_ids:
          self.tags.append(tag)
          tag_ids.append(tag.id)
          self.tag_names.append(tag.name)

    self.trial_tags.ctr.InsertItems(items=self.tag_names, pos=0)

  def onOK(self, e):
    from xfel.ui.db import get_run_path
    self.run_paths = []

    trial = self.trials[self.trials_ctrl.ctr.GetSelection()]
    tags = [self.tags[self.tag_names.index(t)] for t in self.trial_tags.ctr.GetCheckedStrings()]
    tag_ids = [t.id for t in tags]

    run_ids = []
    for rungroup in trial.rungroups:
      for run in rungroup.runs:
        if run.id not in run_ids:
          if len(tags) == 0:
            self.run_paths.append(os.path.join(
              get_run_path(self.parent.main.params.output_folder, trial, rungroup, run),
              'out'))
          else:
            run_tag_ids = [t.id for t in run.tags]
            for tag_id in tag_ids:
              if tag_id in run_tag_ids:
                run_ids.append(run.id)
                self.run_paths.append(os.path.join(
                  get_run_path(self.parent.main.params.output_folder, trial, rungroup, run),
                  'out'))
                break
    e.Skip()

class MultiRunTagDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               db=None,
               *args, **kwargs):
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)
    self.parent = parent
    self.db = db
    self.db_tags = self.db.get_all_tags()
    self.db_tag_names = [t.name for t in self.db_tags]
    self.db_tags_selected_sort = self.db_tag_names
    self.db_runs = self.db.get_all_runs()
    self.db_run_numbers = [str(r.run) for r in self.db_runs]
    self.runs_selected = []
    self.buttons_for_runs_selected = []

    self.multiruntag_sizer = wx.BoxSizer(wx.VERTICAL)

    self.select_runs_panel = wx.Panel(self)
    self.select_runs_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.select_runs =  gctr.CheckListCtrl(self,
                                           label='Selected runs:',
                                           label_size=(200, -1),
                                           label_style='normal',
                                           ctrl_size=(240, 300),
                                           direction='vertical',
                                           choices=self.db_run_numbers)
    self.select_runs_sizer.Add(self.select_runs,
                               flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                               border=10)
    self.select_runs_panel.SetSizer(self.select_runs_sizer)

    self.sort_type_panel = wx.Panel(self)
    self.sort_type_sizer = wx.BoxSizer(wx.VERTICAL)
    self.radio_sort = gctr.RadioCtrl(self.sort_type_panel,
                                     label='Sort tags by',
                                     label_style='normal',
                                     label_size=(80, -1),
                                     direction='horizontal',
                                     items={'id':'DB ID',
                                            'alpha':'name'})
    self.sort_type_sizer.Add(self.radio_sort,
                             flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                             border=5)
    self.sort_type_panel.SetSizer(self.sort_type_sizer)

    self.button_panel = wx.Panel(self)
    self.button_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.btn_add_tags    = wx.Button(self.button_panel, size=(115, -1),
                                     label='Add tags')
    self.btn_remove_tags = wx.Button(self.button_panel, size=(115, -1),
                                     label='Remove tags')
    self.button_sizer.Add(self.btn_add_tags,
                          flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                          border=5)
    self.button_sizer.Add(self.btn_remove_tags,
                          flag=wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP,
                          border=5)
    self.button_panel.SetSizer(self.button_sizer)

    # Add panels to main sizer
    self.multiruntag_sizer.Add(self.select_runs_panel,
                               flag=wx.EXPAND | wx.LEFT | wx.RIGHT)
    self.multiruntag_sizer.Add(self.sort_type_panel,
                               flag=wx.EXPAND | wx.BOTTOM | wx.TOP)
    self.multiruntag_sizer.Add(self.button_panel,
                               flag=wx.EXPAND | wx.BOTTOM | wx.TOP)
    self.main_sizer.Add(self.multiruntag_sizer,
                        flag=wx.EXPAND | wx.RIGHT | wx.LEFT)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Layout()
    self.SetTitle('Manage multiple runs')

    self.Bind(wx.EVT_CHECKLISTBOX, self.onRunChoice, self.select_runs.ctr)
    self.Bind(wx.EVT_RADIOBUTTON, self.onSortDefault, self.radio_sort.id)
    self.Bind(wx.EVT_RADIOBUTTON, self.onSortAlphanum, self.radio_sort.alpha)
    self.Bind(wx.EVT_BUTTON, self.onAddTags, self.btn_add_tags)
    self.Bind(wx.EVT_BUTTON, self.onRemoveTags, self.btn_remove_tags)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onRunChoice(self, e):
    run_numbers_selected = self.select_runs.ctr.GetCheckedStrings()
    self.selected = {}
    for r in self.db_runs:
      if r.run in run_numbers_selected:
        self.selected[r.run] = [r]
    for b in self.parent.all_tag_buttons:
      if b.run.run in self.selected:
        self.selected[b.run.run].append(b)

  def onSortDefault(self, e):
    self.db_tags_selected_sort = self.db_tag_names

  def onSortAlphanum(self, e):
    self.db_tags_selected_sort = sorted(self.db_tag_names)

  def onAddTags(self, e):
    tag_dlg = wx.MultiChoiceDialog(self,
                                   message='Add these tags to selected runs',
                                   caption='Add tags to multiple runs',
                                   choices=self.db_tags_selected_sort)
    tag_dlg.Fit()

    if (tag_dlg.ShowModal() == wx.ID_OK):
      tag_indices = tag_dlg.GetSelections()
      add_tags = [t for t in self.db_tags if self.db_tag_names.index(t.name) in tag_indices]

      runs_with_tags_buttons = self.parent.all_tag_buttons

      for tag in add_tags:
        for run_number in self.selected.keys():
          run, button = self.selected[run_number]
          run_tags = run.app.get_run_tags(run.run_id)
          run_tag_ids = [t.tag_id for t in run_tags]
          if tag.tag_id not in run_tag_ids:
            run.add_tag(tag)
            button.tags.append(tag)
            button.update_label()


  def onRemoveTags(self, e):
    tag_dlg = wx.MultiChoiceDialog(self,
                                   message='Remove these tags from selected runs',
                                   caption='Remove tags from multiple runs',
                                   choices=self.db_tags_selected_sort)
    tag_dlg.Fit()

    if (tag_dlg.ShowModal() == wx.ID_OK):
      tag_indices = tag_dlg.GetSelections()
      remove_tags = [t for t in self.db_tags if self.db_tag_names.index(t.name) in tag_indices]

      for tag in remove_tags:
        for run_number in self.selected.keys():
          run, button = self.selected[run_number]
          run_tags = run.app.get_run_tags(run.run_id)
          run_tag_ids = [t.tag_id for t in run_tags]
          if tag.tag_id in run_tag_ids:
            run.remove_tag(tag)
            tag_on_button = [t for t in button.tags if t.tag_id == tag.tag_id][0] # proxy tag =/= other proxy tag
            button.tags.remove(tag_on_button)
            button.update_label()

  def onOK(self, e):
    e.Skip()

class TagDialog(BaseDialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               db=None,
               *args, **kwargs):
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.db = db
    self.db_tags = self.db.get_all_tags()
    self.deleted_tags = []
    self.new_tags = []
    self.edited_tags =[]
    self.index = 0

    self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.button_panel = wx.Panel(self)
    self.button_sizer = wx.BoxSizer(wx.VERTICAL)
    self.btn_add = wx.Button(self.button_panel, size=(120, -1),
                             label='Add Tag')
    self.btn_rmv = wx.Button(self.button_panel, size=(120, -1),
                             label='Remove Tags')
    self.btn_clr = wx.Button(self.button_panel, size=(120, -1),
                             label='Clear All')
    self.button_sizer.Add(self.btn_add)
    self.button_sizer.Add(self.btn_rmv)
    self.button_sizer.Add(self.btn_clr)
    self.button_panel.SetSizer(self.button_sizer)

    self.tag_panel = ScrolledPanel(self, size=(500, 400))
    self.tag_list = EdListCtrl(self.tag_panel,
                               style=wx.LC_REPORT | wx.SUNKEN_BORDER)
    self.tag_sizer = wx.BoxSizer(wx.VERTICAL)
    self.tag_panel.SetSizer(self.tag_sizer)

    self.tag_list.InsertColumn(0, 'Sample Tag', width=200)
    self.tag_list.InsertColumn(1, 'Comments', width=300)

    # Populate tags with current values from db
    if len(self.db_tags) > 0:
      for tag in self.db_tags:
        self.tag_list.InsertStringItem(self.index, str(tag.name))
        self.tag_list.SetStringItem(self.index, 1, str(tag.comment))
        self.tag_list.SetItemData(self.index, tag.tag_id)
        self.index += 1

    self.tag_sizer.Add(self.tag_list, 1, flag=wx.EXPAND)

    # Add panels to main sizer
    self.top_sizer.Add(self.button_panel,
                       flag=wx.LEFT, border=10)
    self.top_sizer.Add(self.tag_panel,
                       flag=wx.EXPAND | wx.RIGHT | wx.LEFT, border=10)
    self.main_sizer.Add(self.top_sizer,
                        flag=wx.EXPAND| wx.TOP | wx.BOTTOM, border=10)
    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.Layout()
    self.SetTitle('Manage Tags')

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onRemove, self.btn_rmv)
    self.Bind(wx.EVT_BUTTON, self.onClear, self.btn_clr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onAdd(self, e):
    ''' Add a string item to list; focus on item & provide default tag name'''
    new_tag = ('default tag {}'.format(self.index), '', self.index)
    self.new_tags.append(new_tag)
    self.tag_list.InsertStringItem(self.index, new_tag[0])
    self.tag_list.SetStringItem(self.index, 1, new_tag[1])
    self.tag_list.SetItemData(self.index, -1)
    #self.tag_list.Select(self.index)
    #self.tag_list.Focus(self.index)
    self.index += 1

  def onRemove(self, e):
    selected_indices = getListCtrlSelection(self.tag_list)
    tag_ids = [self.tag_list.GetItemData(i) for i in selected_indices]
    self.deleted_tags = [i for i in self.db_tags if i.tag_id in tag_ids]

    for i in range(len(selected_indices)):
      to_delete = getListCtrlSelection(self.tag_list)
      self.tag_list.DeleteItem(to_delete[0])

  def onClear(self, e):
    warning = wx.MessageDialog(self,
                               message='Are you sure you want to delete all '
                                       'tags? \n This cannot be reversed!',
                               caption='Warning',
                               style=wx.YES_NO | wx.ICON_EXCLAMATION)
    if (warning.ShowModal() == wx.ID_YES):
      self.deleted_tags = self.db_tags
      self.tag_list.DeleteAllItems()
      self.index = 0
      self.db_tags = []

  def onOK(self, e):
    try:
      # Delete tags from DB
      for tag in self.deleted_tags:
        self.db.delete_tag(tag=tag)

      # Update names for edited tags
      all_items = [(self.tag_list.GetItemData(i),
                    self.tag_list.GetItem(itemId=i, col=0),
                    self.tag_list.GetItem(itemId=i, col=1))
                    for i in range(self.tag_list.GetItemCount())]

      self.db_tags = self.db.get_all_tags()
      tag_ids = [i.tag_id for i in self.db_tags]
      new_tag_names = [i[1].m_text for i in all_items]

      if len([i for i in new_tag_names if new_tag_names.count(i) > 1]) != 0:
        wx.MessageBox('Need a unique tag name!', 'Warning',
                       wx.ICON_EXCLAMATION)
      else:
        for item in all_items:
          if item[0] in tag_ids:
            tag = self.db.get_tag(tag_id=item[0])
            tag.name = item[1].m_text
            tag.comment = item[2].m_text
          elif item[0] == -1:
            self.db.create_tag(name=item[1].m_text, comment=item[2].m_text)

    except Exception as exception:
      print(str(exception))

    e.Skip()


class RunBlockDialog(BaseDialog):
  ''' Comes up when individual run block button is clicked; allows for run
  block settings to be manipulated by user '''

  def __init__(self, parent, db,
               block=None, trial=None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.parent = parent
    self.block = block
    self.all_blocks = []
    self.db = db
    self.use_ids = db.params.facility.name != 'lcls'

    all_runs = db.get_all_runs()
    if self.use_ids:
      runs_available = sorted([i.id for i in all_runs])
    else:
      runs_available = sorted([int(i.run) for i in all_runs])
    self.first_avail = min(runs_available)
    self.last_avail = max(runs_available)

    if block is None:
      runs = self.db.get_all_runs()
      if self.use_ids:
        run_numbers = [r.id for r in runs]
      else:
        run_numbers = [int(r.run) for r in runs]
      assert len(set(run_numbers)) == len(run_numbers)

      if trial is not None:
        self.all_blocks = trial.rungroups

      if len(runs) == 0:
        wx.MessageBox("No runs found", "Error", wx.OK | wx.ICON_EXCLAMATION)
        assert False # Close and destroy dialog properly here

      self.first_run = min(run_numbers)
      self.last_run = None

      class defaults(object):
        def __getattr__(self, item):
          if item == "detector_address":
            return 'CxiDs2.0:Cspad.0'
          elif item == "detz_parameter":
            return 580
          elif item == "format":
            return "cbf"
          elif item == "two_theta_low":
            return 12.5 # Defaults are from kapton tape experiments (this is kapton ring)
          elif item == "two_theta_high":
            return 22.8 # Defaults are from kapton tape experiments (this is water ring)
          elif item in ["extra_phil_str", "calib_dir", "dark_avg_path", "dark_stddev_path",
            "gain_map_path", "beamx", "beamy", "gain_mask_level", "untrusted_pixel_mask_path",
            "binning", "energy", "comment", "config_str"]:
            return None
          else:
            raise AttributeError(item)
      block = defaults()

    else:
      db = block.app
      self.first_run, self.last_run = block.get_first_and_last_runs()
      if self.first_run is None:
        self.first_run = self.first_avail
      else:
        if self.use_ids:
          if self.first_run is not None: self.first_run = self.first_run.id
          if self.last_run is not None: self.last_run = self.last_run.id
        else:
          if self.first_run is not None: self.first_run = int(self.first_run.run)
          if self.last_run is not None: self.last_run = int(self.last_run.run)

    self.orig_first_run = self.first_run
    self.orig_last_run = self.last_run

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 600),
                        *args, **kwargs)

    # Run block start / end points (choice widgets)

    self.config_panel = wx.Panel(self)
    config_box = wx.StaticBox(self.config_panel, label='Configuration')
    self.config_sizer = wx.StaticBoxSizer(config_box)
    self.config_panel.SetSizer(self.config_sizer)

    self.phil_panel = wx.Panel(self)
    phil_box = wx.StaticBox(self.phil_panel, label='Extra phil parameters')
    self.phil_sizer = wx.StaticBoxSizer(phil_box)
    self.phil_panel.SetSizer(self.phil_sizer)

    runblock_box = wx.StaticBox(self, label='Options')
    self.runblock_box_sizer = wx.StaticBoxSizer(runblock_box, wx.VERTICAL)
    self.runblock_panel = ScrolledPanel(self, size=(550, 225))
    self.runblock_sizer = wx.BoxSizer(wx.VERTICAL)
    self.runblock_panel.SetSizer(self.runblock_sizer)

    # Configuration text ctrl (user can put in anything they want)
    self.config = gctr.PHILBox(self.config_panel,
                               btn_import=True,
                               btn_import_label='Import Config',
                               btn_export=False,
                               btn_default=True,
                               btn_default_label='Default Config',
                               ctr_size=(-1, 100),
                               ctr_value=str(block.config_str))
    self.config_sizer.Add(self.config, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Extra phil
    self.phil = gctr.PHILBox(self.phil_panel,
                             btn_import=True,
                             btn_import_label='Import PHIL',
                             btn_export=False,
                             btn_default=True,
                             btn_default_label='Default PHIL',
                             ctr_size=(-1, 100),
                             ctr_value=str(block.extra_phil_str))
    self.phil_sizer.Add(self.phil, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Image format choice
    if self.parent.trial.app.params.dispatcher == "cxi.xtc_process":
      image_choices = ['pickle']
    else:
      image_choices = ['cbf','pickle']
    self.img_format = gctr.ChoiceCtrl(self.runblock_panel,
                                      label='Image Format:',
                                      label_size=(100, -1),
                                      ctrl_size=(150, -1),
                                      choices=image_choices)
    try:
      self.img_format.ctr.SetSelection(image_choices.index(block.format))
    except Exception:
      pass #in case of selecting an unavailable default
    self.runblock_sizer.Add(self.img_format, flag=wx.TOP | wx.LEFT, border=10)

    self.start_stop_sizer = wx.FlexGridSizer(1, 3, 60, 20)

    self.runblocks_start = gctr.SpinCtrl(self.runblock_panel,
                                   label='Start run:',
                                   label_style='bold',
                                   label_size=(100, -1),
                                   ctrl_value=(self.first_run or self.first_avail),
                                   ctrl_min=self.first_avail,
                                   ctrl_max=self.last_avail)
    self.runblocks_end = gctr.SpinCtrl(self.runblock_panel,
                                   label='End run:',
                                   label_style='bold',
                                   label_size=(100, -1),
                                   ctrl_value=(self.last_run or self.last_avail),
                                   ctrl_min=self.first_avail,
                                   ctrl_max=self.last_avail)
    self.end_type = gctr.RadioCtrl(self.runblock_panel,
                                   label='',
                                   label_style='normal',
                                   label_size=(100, -1),
                                   direction='vertical',
                                   items={'auto':'Auto add runs',
                                          'specify':'Specify end run'})
    self.start_stop_sizer.AddMany([(self.runblocks_start),
                                   (self.runblocks_end),
                                   (self.end_type)])
    self.runblock_sizer.Add(self.start_stop_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Detector address
    self.address = gctr.TextButtonCtrl(self.runblock_panel,
                                       label='Detector Address:',
                                       label_style='bold',
                                       label_size=(100, -1),
                                       value=block.detector_address)
    self.runblock_sizer.Add(self.address, flag=wx.EXPAND | wx.ALL, border=10)


    # Beam XYZ (X, Y - pickle only)
    self.beam_xyz = gctr.OptionCtrl(self.runblock_panel,
                                    label='Beam:',
                                    label_style='bold',
                                    label_size=(100, -1),
                                    ctrl_size=(60, -1),
                                    items=[('X', block.beamx),
                                           ('Y', block.beamy),
                                           ('DetZ', block.detz_parameter)])
    self.runblock_sizer.Add(self.beam_xyz, flag=wx.EXPAND | wx.ALL, border=10)

    # Binning, energy, gain mask level
    self.bin_nrg_gain = gctr.OptionCtrl(self.runblock_panel,
                                        ctrl_size=(80, -1),
                                        items=[('binning', block.binning),
                                               ('energy', block.energy),
                                               ('gain_mask_level', block.gain_mask_level)])
    self.runblock_sizer.Add(self.bin_nrg_gain, flag=wx.EXPAND | wx.ALL, border=10)

    # Two theta values for droplet hit finding
    self.two_thetas = gctr.OptionCtrl(self.runblock_panel,
                                      ctrl_size=(80, -1),
                                      items=[('two_theta_low', block.two_theta_low),
                                             ('two_theta_high', block.two_theta_high)])
    self.runblock_sizer.Add(self.two_thetas, flag=wx.EXPAND | wx.ALL, border=10)


    # Untrusted pixel mask path
    self.untrusted_path = gctr.TextButtonCtrl(self.runblock_panel,
                                              label='Untrusted Pixel Mask:',
                                              label_style='normal',
                                              label_size=(100, -1),
                                              big_button=True,
                                              value=str(block.untrusted_pixel_mask_path))
    self.runblock_sizer.Add(self.untrusted_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Calibration folder
    self.calib_dir = gctr.TextButtonCtrl(self.runblock_panel,
                                         label='Calibration:',
                                         label_style='normal',
                                         label_size=(100, -1),
                                         big_button=True,
                                         value=str(block.calib_dir))
    self.runblock_sizer.Add(self.calib_dir, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark average path (pickle only)
    self.dark_avg_path = gctr.TextButtonCtrl(self.runblock_panel,
                                             label='Dark Average:',
                                             label_style='normal',
                                             label_size=(100, -1),
                                             big_button=True,
                                             value=str(block.dark_avg_path))
    self.runblock_sizer.Add(self.dark_avg_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark stddev path (pickle only)
    self.dark_stddev_path = gctr.TextButtonCtrl(self.runblock_panel,
                                                label='Dark StdDev:',
                                                label_style='normal',
                                                label_size=(100, -1),
                                                big_button=True,
                                                value=str(block.dark_stddev_path))
    self.runblock_sizer.Add(self.dark_stddev_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark map path (pickle only)
    self.gain_map_path = gctr.TextButtonCtrl(self.runblock_panel,
                                             label='Gain Map:',
                                             label_style='normal',
                                             label_size=(100, -1),
                                             big_button=True,
                                             value=str(block.gain_map_path))
    self.runblock_sizer.Add(self.gain_map_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Comment
    self.comment = gctr.TextButtonCtrl(self.runblock_panel,
                                       label='Comment:',
                                       label_style='normal',
                                       label_size=(100, -1))
    self.comment.ctr.SetValue(str(block.comment))
    self.runblock_sizer.Add(self.comment, flag=wx.EXPAND | wx.ALL,
                            border=10)

    self.main_sizer.Add(self.phil_panel, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.config_panel, flag=wx.EXPAND | wx.ALL, border=10)
    self.runblock_box_sizer.Add(self.runblock_panel)
    self.main_sizer.Add(self.runblock_box_sizer, flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_RADIOBUTTON, self.onAutoEnd, self.end_type.auto)
    self.Bind(wx.EVT_RADIOBUTTON, self.onSpecifyEnd, self.end_type.specify)
    self.Bind(wx.EVT_BUTTON, self.onDarkAvgBrowse,
              id=self.dark_avg_path.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onImportConfig, self.config.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDefaultConfig, self.config.btn_default)
    self.Bind(wx.EVT_BUTTON, self.onImportPhil, self.phil.btn_import)
    self.Bind(wx.EVT_BUTTON, self.onDarkMapBrowse,
              id=self.gain_map_path.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onDarkStdBrowse,
              id=self.dark_stddev_path.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onCalibDirBrowse,
              id=self.calib_dir.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onUntrustedBrowse,
              id=self.untrusted_path.btn_big.GetId())
    self.Bind(wx.EVT_CHOICE, self.onImageFormat, id=self.img_format.ctr.GetId())
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)


    self.fill_in_fields()
    self.configure_controls()
    self.Layout()
    self.SetTitle('Run Block Settings')

  def onAutoEnd(self, e):
    self.last_run = None
    self.runblocks_end.Disable()

  def onSpecifyEnd(self, e):
    self.runblocks_end.Enable()

  def onImportConfig(self, e):
    cfg_dlg = wx.FileDialog(self,
                            message="Load configuration file",
                            defaultDir=os.curdir,
                            defaultFile="*",
                            wildcard="*",
                            style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                            )
    if cfg_dlg.ShowModal() == wx.ID_OK:
      config_file = cfg_dlg.GetPaths()[0]
      with open(config_file, 'r') as cfg:
        cfg_contents = cfg.read()
      self.config.ctr.SetValue(cfg_contents)
    cfg_dlg.Destroy()

  def onImportPhil(self, e):
    phil_dlg = wx.FileDialog(self,
                             message="Load phil file",
                             defaultDir=os.curdir,
                             defaultFile="*",
                             wildcard="*.phil",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if phil_dlg.ShowModal() == wx.ID_OK:
      phil_file = phil_dlg.GetPaths()[0]
      with open(phil_file, 'r') as phil:
        phil_contents = phil.read()
      self.phil.ctr.SetValue(phil_contents)
    phil_dlg.Destroy()

  def onDefaultConfig(self, e):
    # TODO: Generate default config parameters (re-do based on pickle / CBF)
    pass

  def onOK(self, e):
    try:
      first = int(self.runblocks_start.ctr.GetValue())
      assert first > 0 and first >= self.first_avail
      self.first_run = first
    except (ValueError, AssertionError) as e:
      print("Please select a run between %d and %d." % (self.first_avail, self.last_avail))
      raise e
    if self.end_type.specify.GetValue() == 1:
      try:
        last = int(self.runblocks_end.ctr.GetValue())
        assert last > 0 and last <= self.last_avail and last >= first
        self.last_run = last
      except (ValueError, AssertionError) as e:
        print("Please select a run between %d and %d." % (self.first_run, self.last_avail))
        raise e
    elif self.end_type.specify.GetValue() == 0:
      self.last_run = None
    else:
      assert False
    rg_open = self.last_run is None

    rg_dict = dict(active=True,
                   open=rg_open,
                   format=self.img_format.ctr.GetStringSelection(),
                   config_str=self.config.ctr.GetValue(),
                   extra_phil_str=self.phil.ctr.GetValue(),
                   detector_address=self.address.ctr.GetValue(),
                   detz_parameter=self.beam_xyz.DetZ.GetValue(),
                   beamx=self.beam_xyz.X.GetValue(),
                   beamy=self.beam_xyz.Y.GetValue(),
                   binning=self.bin_nrg_gain.binning.GetValue(),
                   energy=self.bin_nrg_gain.energy.GetValue(),
                   untrusted_pixel_mask_path=self.untrusted_path.ctr.GetValue(),
                   dark_avg_path=self.dark_avg_path.ctr.GetValue(),
                   dark_stddev_path=self.dark_stddev_path.ctr.GetValue(),
                   gain_map_path=self.gain_map_path.ctr.GetValue(),
                   gain_mask_level=self.bin_nrg_gain.gain_mask_level.GetValue(),
                   calib_dir=self.calib_dir.ctr.GetValue(),
                   two_theta_low=self.two_thetas.two_theta_low.GetValue(),
                   two_theta_high=self.two_thetas.two_theta_high.GetValue(),
                   comment=self.comment.ctr.GetValue())
    for key, value in six.iteritems(rg_dict):
      if str(value) == 'None' or str(value) == '':
        rg_dict[key] = None
      elif type(value) == bool:
        rg_dict[key] = int(value)

    if self.block is None:
      self.block = self.db.create_rungroup(**rg_dict)
      self.block.sync_runs(self.first_run, self.last_run, self.use_ids)
      self.parent.trial.add_rungroup(self.block)
    else:
      # if all the parameters are unchanged, do nothing
      all_the_same = [str(rg_dict[key]) == str(getattr(self.block, key)) for key in rg_dict].count(False) == 0
      all_the_same &= self.first_run == self.orig_first_run and self.last_run == self.orig_last_run
      if not all_the_same:
        # if all the parameters except open and comment are the same,
        # only update those fields
        keep_old_run_group = [str(rg_dict[key]) == str(getattr(self.block, key)) for key in rg_dict \
                              if key not in ['open', 'comment']].count(False) == 0
        if keep_old_run_group:
          main = self.parent.parent.GetParent().main
          running = main.job_sentinel is not None and main.job_sentinel.active
          if running:
            main.stop_job_sentinel()

          self.block.open = rg_open
          self.block.comment = rg_dict['comment']
          self.block.sync_runs(self.first_run, self.last_run, self.use_ids)

          if running:
            main.start_job_sentinel()

        else:
          # enough parameters have changed to warrant creating a new run group
          self.block.active = False
          self.block = self.db.create_rungroup(**rg_dict)
          self.block.sync_runs(self.first_run, self.last_run, self.use_ids)
          self.parent.trial.add_rungroup(self.block)

    e.Skip()

  def fill_in_fields(self):
    ''' If previous rungroups exist in trial, fill in fields in nascent block '''
    if len(self.all_blocks) > 0:
      last = self.all_blocks[-1]
      self.config.ctr.SetValue(str(last.config_str))
      self.phil.ctr.SetValue(str(last.extra_phil_str))
      self.address.ctr.SetValue(str(last.detector_address))
      self.beam_xyz.DetZ.SetValue(str(last.detz_parameter))
      self.beam_xyz.X.SetValue(str(last.beamx))
      self.beam_xyz.Y.SetValue(str(last.beamy))
      self.bin_nrg_gain.binning.SetValue(str(last.binning))
      self.bin_nrg_gain.energy.SetValue(str(last.energy))
      self.bin_nrg_gain.gain_mask_level.SetValue(str(last.gain_mask_level))
      self.two_thetas.two_theta_low.SetValue(str(last.two_theta_low))
      self.two_thetas.two_theta_high.SetValue(str(last.two_theta_high))
      self.untrusted_path.ctr.SetValue(str(last.untrusted_pixel_mask_path))
      self.dark_avg_path.ctr.SetValue(str(last.dark_avg_path))
      self.dark_stddev_path.ctr.SetValue(str(last.dark_stddev_path))
      self.gain_map_path.ctr.SetValue(str(last.gain_map_path))
      self.calib_dir.ctr.SetValue(str(last.calib_dir))
      self.comment.ctr.SetValue(str(last.comment))

  def onImageFormat(self, e):
    sel= self.img_format.ctr.GetString(self.img_format.ctr.GetSelection())
    if 'cbf' in sel:
      #self.beam_xyz.X.Disable()
      #self.beam_xyz.Y.Disable()
      #self.bin_nrg_gain.binning.Disable()
      #self.two_thetas.two_theta_low.Disable()
      #self.two_thetas.two_theta_high.Disable()
      #self.dark_avg_path.Hide()
      #self.dark_stddev_path.Hide()
      #self.gain_map_path.Hide()
      self.runblock_panel.SetupScrolling()
    elif 'pickle' in sel:
      #self.beam_xyz.X.Enable()
      #self.beam_xyz.Y.Enable()
      #self.bin_nrg_gain.binning.Enable()
      #self.two_thetas.two_theta_low.Enable()
      #self.two_thetas.two_theta_high.Enable()
      #self.dark_avg_path.Show()
      #self.dark_stddev_path.Show()
      #self.gain_map_path.Show()
      self.runblock_panel.Layout()
      self.runblock_panel.SetupScrolling()

  def configure_controls(self):
    self.onImageFormat(None)
    if self.last_run is None:
      self.runblocks_end.Disable()
      self.end_type.auto.SetValue(1)
    else:
      self.runblocks_end.Enable()
      self.end_type.specify.SetValue(1)

  def onDarkAvgBrowse(self, e):
    dark_dlg = wx.FileDialog(self,
                             message="Load dark average file",
                             defaultDir=os.curdir,
                             defaultFile="*.cbf",
                             wildcard="*.cbf",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )

    if dark_dlg.ShowModal() == wx.ID_OK:
      self.dark_avg_path.ctr.SetValue(dark_dlg.GetPaths()[0])
    dark_dlg.Destroy()

  def onDarkMapBrowse(self, e):
    dark_dlg = wx.FileDialog(self,
                             message="Load dark map file",
                             defaultDir=os.curdir,
                             defaultFile="*.cbf",
                             wildcard="*.cbf",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )

    if dark_dlg.ShowModal() == wx.ID_OK:
      self.gain_map_path.ctr.SetValue(dark_dlg.GetPaths()[0])
    dark_dlg.Destroy()

  def onDarkStdBrowse(self, e):
    dark_dlg = wx.FileDialog(self,
                             message="Load dark stddev file",
                             defaultDir=os.curdir,
                             defaultFile="*.cbf",
                             wildcard="*.cbf",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )

    if dark_dlg.ShowModal() == wx.ID_OK:
      self.dark_stddev_path.ctr.SetValue(dark_dlg.GetPaths()[0])
    dark_dlg.Destroy()

  def onUntrustedBrowse(self, e):
    dlg = wx.FileDialog(self,
                             message="Load untrusted pixel mask",
                             defaultDir=os.curdir,
                             defaultFile="*.pickle",
                             wildcard="*.pickle",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )

    if dlg.ShowModal() == wx.ID_OK:
      self.untrusted_path.ctr.SetValue(dlg.GetPaths()[0])
    dlg.Destroy()

  def onCalibDirBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose calibration directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.calib_dir.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

class SelectRunBlocksDialog(BaseDialog):
  def __init__(self, parent, trial,
               label_style='bold',
               content_style='normal',
               db=None,
               *args, **kwargs):
    BaseDialog.__init__(self, parent, label_style=label_style,
                        content_style=content_style, *args, **kwargs)

    self.db = db
    self.trial = trial

    self.top_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.button_panel = wx.Panel(self)
    self.button_sizer = wx.BoxSizer(wx.VERTICAL)
    self.button_panel.SetSizer(self.button_sizer)

    self.runblocks_panel = ScrolledPanel(self, size=(500, 400))

    # Populate rungroups with current values from db
    self.trial_rungroups = [t.id for t in trial.rungroups]
    self.all_rungroups = self.db.get_all_rungroups()
    choices = []
    selected = []
    for rungroup in self.all_rungroups:
      selected.append(rungroup.id in self.trial_rungroups)
      first_run, last_run = rungroup.get_first_and_last_runs()
      if last_run is None:
        if first_run is None:
          desc = "[%d]"%(rungroup.id)
        else:
          desc = "[%d] %d+"%(rungroup.id, int(first_run.id))
      else:
        desc = "[%d] %d-%d"%(rungroup.id, int(first_run.id), int(last_run.id))
      if rungroup.comment is not None:
        desc += " " + rungroup.comment

      choices.append(desc)

    self.runblocks_list = gctr.CheckListCtrl(self.runblocks_panel,
                                             label='Select runblocks',
                                             label_size=(40, -1),
                                             label_style='normal',
                                             ctrl_size=(450, 350),
                                             direction='vertical',
                                             choices=choices)
    for i in range(len(selected)):
      self.runblocks_list.ctr.Check(i, selected[i])

    self.runblocks_sizer = wx.BoxSizer(wx.VERTICAL)
    self.runblocks_panel.SetSizer(self.runblocks_sizer)

    self.runblocks_sizer.Add(self.runblocks_list, 1, flag=wx.EXPAND)

    # Add panels to main sizer
    self.top_sizer.Add(self.button_panel,
                       flag=wx.LEFT, border=10)
    self.top_sizer.Add(self.runblocks_panel,
                       flag=wx.EXPAND | wx.RIGHT | wx.LEFT, border=10)
    self.main_sizer.Add(self.top_sizer,
                        flag=wx.EXPAND| wx.TOP | wx.BOTTOM, border=10)
    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                   flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                   border=10)

    self.Layout()
    self.SetTitle('Select run blocks')

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onOK(self, e):
    trial_rungroups = [t.id for t in self.db.get_trial_rungroups(self.trial.id, only_active=False)]
    trials = self.db.get_all_trials()

    for i, rungroup in enumerate(self.all_rungroups):
      if self.runblocks_list.ctr.IsChecked(i):
        if not rungroup.id in trial_rungroups:
          self.trial.add_rungroup(rungroup)
        rungroup.active = True
      else:
        if rungroup.id in trial_rungroups:
          self.trial.remove_rungroup(rungroup)
        if not any([rungroup.id in [rg.id for rg in t.rungroups] for t in trials]):
          rungroup.active = False

    e.Skip()

class TrialDialog(BaseDialog):
  def __init__(self, parent, db,
               new=True,
               trial=None,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.db = db
    self.new = new
    self.trial = trial
    self.all_trials = db.get_all_trials()

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 600),
                        *args, **kwargs)

    if trial is None:
      trials = [t.trial for t in self.all_trials]
      if len(trials) == 0:
        trial_number = 0
      else:
        trial_number = max(trials) + 1
      d_min = 1.5
    else:
      trial_number = trial.trial
      d_min = trial.d_min if trial.d_min is not None else 1.5

    self.trial_info = gctr.TwoButtonCtrl(self,
                                         label='Trial number:',
                                         label_size=(100, -1),
                                         label_style='bold',
                                         button1=True,
                                         button1_label='Import PHIL',
                                         button1_size=(120, -1),
                                         button2=True,
                                         button2_label='Default PHIL',
                                         button2_size=(120, -1),
                                         value="{}".format(trial_number))
    self.trial_comment = gctr.TextButtonCtrl(self,
                                             label='Comment:',
                                             label_size=(100, -1),
                                             label_style='bold',
                                             ghost_button=False)

    self.phil_box = rt.RichTextCtrl(self, style=wx.VSCROLL, size=(-1, 400))

    choices = [('None', None)] + \
              [('Trial {}'.format(t.trial), t.trial) for t in self.all_trials]
    self.copy_runblocks = gctr.ChoiceCtrl(self,
                                          label='Copy runblocks from',
                                          label_style='normal',
                                          label_size=(180, -1),
                                          ctrl_size=(80, -1),
                                          choices=choices)
    self.copy_runblocks.ctr.SetSelection(0)
    self.throttle = gctr.SpinCtrl(self,
                                  label='Percent events processed:',
                                  label_size=(180, -1),
                                  label_style='bold',
                                  ctrl_size=(80, -1),
                                  ctrl_value='100',
                                  ctrl_min=1,
                                  ctrl_max=100)
    self.num_bins = gctr.SpinCtrl(self,
                                  label='Number of bins:',
                                  label_size=(180, -1),
                                  label_style='bold',
                                  ctrl_size=(80, -1),
                                  ctrl_value='20',
                                  ctrl_min=1,
                                  ctrl_max=100,
                                  ctrl_step=1)
    self.d_min = gctr.SpinCtrl(self,
                               label='High res. limit ({}):'
                               ''.format(u'\N{ANGSTROM SIGN}'.encode('utf-8')),
                               label_size=(180, -1),
                               label_style='bold',
                               ctrl_size=(80, -1),
                               ctrl_value=str(d_min),
                               ctrl_min=0.1,
                               ctrl_max=100,
                               ctrl_step=0.1,
                               ctrl_digits=1)

    self.option_sizer = wx.FlexGridSizer(3, 2, 10, 20)
    self.option_sizer.AddMany([(self.copy_runblocks),
                               (0, 0),
                               (self.throttle),
                               (0, 0),
                               (self.num_bins),
                               (self.d_min)])

    self.main_sizer.Add(self.trial_info,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.trial_comment,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.phil_box, 1,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.option_sizer, flag=wx.EXPAND | wx.ALL, border=10)


    # Dialog control
    if self.new:
      dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    else:
      dialog_box = self.CreateSeparatedButtonSizer(wx.OK)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Layout()

    if self.new:
      self.SetTitle('New Trial Settings')

      # If previous trials exist, propagate previous settings from them
      if len(self.all_trials) > 0:
        self.last_trial = self.all_trials[-1]
        target_phil_str = self.last_trial.target_phil_str
        self.trial_comment.ctr.SetValue(str(self.last_trial.comment))
        process_percent = self.last_trial.process_percent

      else:
        target_phil_str = ''
        process_percent = 100

    else:
      self.SetTitle('Trial {}'.format(self.trial.trial))
      target_phil_str = self.trial.target_phil_str
      self.trial_comment.ctr.SetValue(str(self.trial.comment))
      process_percent = self.trial.process_percent

      # Disable controls for viewing
      self.trial_info.button1.Disable()
      self.trial_info.button2.Disable()
      self.trial_info.ctr.SetEditable(False)
      self.copy_runblocks.Hide()
      self.phil_box.SetEditable(False)
      self.throttle.ctr.Disable()
      self.num_bins.ctr.Disable()
      self.d_min.ctr.Disable()

    if target_phil_str is None:
      target_phil_str = ""
    if process_percent is None:
      process_percent = 100
    self.phil_box.SetValue(target_phil_str)
    self.throttle.ctr.SetValue(process_percent)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onBrowse, self.trial_info.button1)
    self.Bind(wx.EVT_BUTTON, self.onDefault, self.trial_info.button2)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onBrowse(self, e):
    ''' Open dialog for selecting PHIL file '''

    load_dlg = wx.FileDialog(self,
                             message="Load PHIL file",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      target_file = load_dlg.GetPaths()[0]
      with open(target_file, 'r') as phil_file:
        phil_file_contents = phil_file.read()
      self.phil_box.SetValue(phil_file_contents)
    load_dlg.Destroy()

  def onDefault(self, e):
    # TODO: Generate default PHIL parameters
    pass

  def onOK(self, e):
    if self.new:
      target_phil_str = self.phil_box.GetValue()

      # Parameter validation
      dispatcher = self.db.params.dispatcher
      if dispatcher == 'cxi.xtc_process': #LABELIT
        from spotfinder.applications.xfel import cxi_phil
        phil_scope = cxi_phil.cxi_versioned_extract().persist.phil_scope
      else:
        from xfel.ui import load_phil_scope_from_dispatcher
        phil_scope = load_phil_scope_from_dispatcher(dispatcher)

      from iotbx.phil import parse
      msg = None
      try:
        trial_params, unused = phil_scope.fetch(parse(target_phil_str), track_unused_definitions = True)
      except Exception as e:
        msg = '\nParameters incompatible with %s dispatcher:\n%s\n' % (dispatcher, str(e))
      else:
        if len(unused) > 0:
          msg = [str(item) for item in unused]
          msg = '\n'.join(['  %s' % line for line in msg])
          msg = 'The following definitions were not recognized:\n%s\n' % msg

        try:
          params = trial_params.extract()
        except Exception as e:
          if msg is None: msg = ""
          msg += '\nOne or more values could not be parsed:\n%s\n' % str(e)

      if msg is not None:
        msg += '\nFix the parameters and press OK again'
        msgdlg = wx.MessageDialog(self,
                                  message=msg,
                                  caption='Warning',
                                  style=wx.OK |  wx.ICON_EXCLAMATION)
        msgdlg.ShowModal()
        return

      comment = self.trial_comment.ctr.GetValue()
      process_percent = int(self.throttle.ctr.GetValue())
      d_min = float(self.d_min.ctr.GetValue())
      n_bins = int(self.num_bins.ctr.GetValue())
      if process_percent == 100:
        process_percent = None

      if self.trial is None:
        self.db.create_trial(
          trial = int(self.trial_info.ctr.GetValue()),
          active = False,
          target_phil_str = target_phil_str,
          comment = comment,
          process_percent = process_percent,
          d_min = d_min,
          n_bins = n_bins)
        self.trial = self.db.get_trial(trial_number=int(self.trial_info.ctr.GetValue()))

        template_trial_number = self.copy_runblocks.ctr.GetClientData(
          self.copy_runblocks.ctr.GetSelection())

        if template_trial_number is not None:
          template_trial = self.db.get_trial(trial_number=template_trial_number)
          for block in template_trial.rungroups:
            self.trial.add_rungroup(block)

      else:
        self.trial.target_phil_str = target_phil_str
        self.trial.comment = comment
    else:
      if self.trial.comment != self.trial_comment.ctr.GetValue():
        self.trial.comment = self.trial_comment.ctr.GetValue()
    e.Skip()
