from __future__ import division

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
                        size=(600, 200),
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
                                       value=self.params.experiment_tag)
    self.main_sizer.Add(self.db_cred,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)

    # Experiment name control
    self.experiment = gctr.TextButtonCtrl(self,
                                          label='Experiment',
                                          label_style='bold',
                                          label_size=(150, -1),
                                          big_button_size=(130, -1),
                                          value=self.params.experiment)
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

  def onBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose the input directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.output.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

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
      self.params.db.name     = creds.db_name.ctr.GetValue()
      self.params.db.user     = creds.db_user.ctr.GetValue()
      self.params.db.password = creds.db_password.ctr.GetValue()
      self.params.web.user     = creds.web_user.ctr.GetValue()
      self.params.web.password = creds.web_password.ctr.GetValue()

      self.drop_tables = creds.chk_drop_tables.GetValue()


  def onOK(self, e):
    self.params.experiment_tag = self.db_cred.ctr.GetValue()
    self.params.experiment = self.experiment.ctr.GetValue()
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
                                           label='Password',
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
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)

    # LCLS user name
    self.web_user = gctr.TextButtonCtrl(self,
                                       label='LCLS user name',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       big_button_size=(130, -1),
                                       value=params.web.user)
    self.main_sizer.Add(self.web_user, flag=wx.EXPAND | wx.ALL, border=10)

    # LCLS password
    self.web_password = gctr.TextButtonCtrl(self,
                                           label='LCLS Password',
                                           label_style='bold',
                                           label_size=(150, -1),
                                           text_style=wx.TE_PASSWORD,
                                           big_button_size=(130, -1),
                                           value=params.web.password)
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

    choices = ['python', 'lsf', 'mpi', 'sge', 'pbi', 'custom']
    self.mp_option = gctr.ChoiceCtrl(self,
                                     label='Multiprocessing:',
                                     label_size=(120, -1),
                                     label_style='bold',
                                     choices=choices)
    self.mp_sizer.Add(self.mp_option, flag=wx.EXPAND | wx.ALL, border=10)
    try:
      self.mp_option.ctr.SetSelection(choices.index(params.mp.method))
    except ValueError:
      pass

    # Queue
    queues = ['psanaq', 'psanaq', 'psdebugq','psanaidleq', 'psnehhiprioq',
              'psnehprioq', 'psnehq', 'psfehhiprioq', 'psfehprioq', 'psfehq']
    self.queue = gctr.ChoiceCtrl(self,
                                 label='Queue:',
                                 label_size=(120, -1),
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
                               label_size=(120, -1),
                               label_style='normal',
                               ctrl_size=(100, -1),
                               ctrl_value='%d'%params.mp.nproc,
                               ctrl_min=1,
                               ctrl_max=1000)
    self.mp_sizer.Add(self.nproc, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.mp_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Data analysis settings
    analysis_box = wx.StaticBox(self, label='Data Analysis Options')
    self.analysis_sizer = wx.StaticBoxSizer(analysis_box, wx.VERTICAL)

    # Processing back-ends
    self.back_ends = ['DIALS', 'LABELIT']
    self.dispatchers = ['cctbx.xfel.xtc_process', 'cxi.xtc_process']
    self.back_end = gctr.ChoiceCtrl(self,
                                    label='Processing back end:',
                                    label_size=(120, -1),
                                    label_style='bold',
                                    choices=self.back_ends)
    self.analysis_sizer.Add(self.back_end, flag=wx.EXPAND | wx.ALL, border=10)
    try:
      self.back_end.ctr.SetSelection(self.dispatchers.index(params.dispatcher))
    except ValueError:
      pass

    img_types = ['corrected', 'raw']
    self.avg_img_type = gctr.ChoiceCtrl(self,
                                        label='Avg. Image Type:',
                                        label_size=(120, -1),
                                        label_style='bold',
                                        choices=img_types)
    if params.average_raw_data:
      i = img_types.index('raw')
    else:
      i = img_types.index('corrected')
    self.avg_img_type.ctr.SetSelection(i)

    self.analysis_sizer.Add(self.avg_img_type, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.analysis_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    self.chk_enforce80 = wx.CheckBox(self,
                                     label='Require stream 80 (FEE spectrometer) before processing')
    self.chk_enforce80.SetValue(params.web.enforce80)

    self.main_sizer.Add(self.chk_enforce80, flag=wx.ALL, border=10)

    self.chk_enforce81 = wx.CheckBox(self,
                                     label='Require stream 81 (FEE spectrometer) before processing')
    self.chk_enforce81.SetValue(params.web.enforce81)

    self.main_sizer.Add(self.chk_enforce81, flag=wx.ALL, border=10)

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

  def onOK(self, e):
    self.params.dispatcher = self.dispatchers[self.back_end.ctr.GetSelection()]
    self.params.mp.method = self.mp_option.ctr.GetStringSelection()
    self.params.mp.queue = self.queue.ctr.GetStringSelection()
    self.params.mp.nproc = int(self.nproc.ctr.GetValue())
    self.params.average_raw_data = self.avg_img_type.ctr.GetStringSelection() == 'raw'
    self.params.web.enforce80 = bool(self.chk_enforce80.GetValue())
    self.params.web.enforce81 = bool(self.chk_enforce81.GetValue())
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
    from xfel.command_line.cxi_mpi_submit import get_submit_command
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
    command = str(get_submit_command(command, submit_path, working_dir, params.mp)) # comes back as unicode which throws off easy_run
    cwd = os.getcwd()
    os.makedirs(working_dir)
    os.chdir(working_dir)

    print "Submitting metrology refinement. Command:"
    print command
    try:
      results = easy_run.fully_buffered(command=command)
      results.show_stdout()
      results.show_stderr()
      results.raise_if_errors()
    except Exception, exc:
      if not "Warning: job being submitted without an AFS token." in str(exc):
        raise exc
    os.chdir(cwd)
    print "Output will be in", working_dir

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

    # Calib folder
    self.calib_dir = gctr.TextButtonCtrl(self,
                                         label='Calibration:',
                                         label_style='normal',
                                         label_size=(150, -1),
                                         ctrl_size=(350, -1),
                                         big_button=True)
    self.main_sizer.Add(self.calib_dir, flag=wx.EXPAND | wx.ALL, border=10)
    self.Bind(wx.EVT_BUTTON, self.onCalibDirBrowse, self.calib_dir.btn_big)

    # Detector Address
    self.address = gctr.TextButtonCtrl(self,
                                       label='Detector Address:',
                                       label_style='bold',
                                       label_size=(150, -1),
                                       value='')
    self.main_sizer.Add(self.address, flag=wx.EXPAND | wx.ALL, border=10)

    # Raw vs. corrected
    img_types = ['raw', 'corrected']
    self.avg_img_type = gctr.ChoiceCtrl(self,
                                        label='Avg. Image Type:',
                                        label_size=(150, -1),
                                        label_style='bold',
                                        choices=img_types)
    self.main_sizer.Add(self.avg_img_type, flag=wx.EXPAND | wx.ALL, border=10)

    # DetZ
    self.detz = gctr.OptionCtrl(self,
                                label='DetZ:',
                                label_style='bold',
                                label_size=(150, -1),
                                ctrl_size=(100, -1),
                                items=[('DetZ', 580)])
    self.main_sizer.Add(self.detz, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

  def onCalibDirBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose calibration directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.calib_dir.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()
    e.Skip()

    self.SetTitle('Averaging Settings')

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

    except Exception, exception:
      print str(exception)

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

    if block is None:
      runs = self.db.get_all_runs()
      run_numbers = [r.run for r in runs]
      assert len(set(run_numbers)) == len(run_numbers)

      if trial is not None:
        self.all_blocks = trial.rungroups

      if len(runs) == 0:
        wx.MessageBox("No runs found", "Error", wx.OK | wx.ICON_EXCLAMATION)
        assert False # Close and destroy dialog properly here

      min_run = runs[run_numbers.index(min(run_numbers))]
      self.first_run = min_run.run
      self.last_run = None

      class defaults(object):
        def __getattr__(self, item):
          if item == "detector_address":
            return 'CxiDs2.0:Cspad.0'
          elif item == "detz_parameter":
            return 580
          else:
            return None
      block = defaults()

    else:
      db = block.app
      self.first_run = db.get_run(run_id=block.startrun).run
      if block.endrun is None:
        self.last_run = None
      else:
        self.last_run = db.get_run(run_id=block.endrun).run

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 600),
                        *args, **kwargs)

    # Run block start / end points (choice widgets)

    start = [str(i.run) for i in db.get_all_runs()]
    stop = start + ['None']
    firstidx = start.index(str(self.first_run))
    lastidx = stop.index(str(self.last_run))

    self.option_sizer = wx.FlexGridSizer(1, 3, 0, 10)
    self.option_sizer.AddGrowableCol(2)

    self.config_panel = wx.Panel(self)
    config_box = wx.StaticBox(self.config_panel, label='Configuration')
    self.config_sizer = wx.StaticBoxSizer(config_box)
    self.config_panel.SetSizer(self.config_sizer)

    runblock_box = wx.StaticBox(self, label='Options')
    self.runblock_box_sizer = wx.StaticBoxSizer(runblock_box, wx.VERTICAL)
    self.runblock_panel = ScrolledPanel(self, size=(550, 350))
    self.runlbock_sizer = wx.BoxSizer(wx.VERTICAL)
    self.runblock_panel.SetSizer(self.runlbock_sizer)

    # Image format choice
    self.img_format = gctr.ChoiceCtrl(self,
                                      label='Image Format:',
                                      label_size=(100, -1),
                                      ctrl_size=(150, -1),
                                      choices=['CBF', 'pickle'])
    self.btn_import_cfg = wx.Button(self, label='Import Config', size=(120, -1))
    self.btn_default_cfg = wx.Button(self, label='Default Config', size=(120, -1))
    self.option_sizer.Add(self.btn_import_cfg, flag=wx.ALIGN_RIGHT)
    self.option_sizer.Add(self.btn_default_cfg, flag=wx.ALIGN_RIGHT)
    self.option_sizer.Add(self.img_format, flag=wx.ALIGN_RIGHT)

    # Configuration text ctrl (user can put in anything they want)
    self.config = rt.RichTextCtrl(self.config_panel,
                                  size=(-1, 250),
                                  style=wx.VSCROLL,
                                  value=str(block.config_str))
    self.config_sizer.Add(self.config, 1, flag=wx.EXPAND | wx.ALL, border=10)

    # Runblock runs choice
    self.runblocks = gctr.MultiChoiceCtrl(self.runblock_panel,
                                          label='Run block:',
                                          label_style='bold',
                                          label_size=(100, -1),
                                          ctrl_size=(150, -1),
                                          items={'start':start, 'end':stop})
    self.runblocks.start.SetSelection(firstidx)
    self.runblocks.end.SetSelection(lastidx)
    self.runlbock_sizer.Add(self.runblocks, flag=wx.EXPAND | wx.ALL, border=10)

    # Detector address
    self.address = gctr.TextButtonCtrl(self.runblock_panel,
                                       label='Detector Address:',
                                       label_style='bold',
                                       label_size=(100, -1),
                                       value=block.detector_address)
    self.runlbock_sizer.Add(self.address, flag=wx.EXPAND | wx.ALL, border=10)


    # Beam XYZ (X, Y - pickle only)
    self.beam_xyz = gctr.OptionCtrl(self.runblock_panel,
                                    label='Beam:',
                                    label_style='bold',
                                    label_size=(100, -1),
                                    ctrl_size=(60, -1),
                                    items=[('X', block.beamx),
                                           ('Y', block.beamy),
                                           ('DetZ', block.detz_parameter)])
    self.runlbock_sizer.Add(self.beam_xyz, flag=wx.EXPAND | wx.ALL, border=10)

    # Binning, energy, gain mask level
    self.bin_nrg_gain = gctr.OptionCtrl(self.runblock_panel,
                                        ctrl_size=(80, -1),
                                        items=[('binning', block.binning),
                                               ('energy', block.energy),
                                               ('gain_mask_level', block.gain_mask_level)])
    self.runlbock_sizer.Add(self.bin_nrg_gain, flag=wx.EXPAND | wx.ALL, border=10)

    # Untrusted pixel mask path
    self.untrusted_path = gctr.TextButtonCtrl(self.runblock_panel,
                                              label='Untrusted Pixel Mask:',
                                              label_style='normal',
                                              label_size=(100, -1),
                                              big_button=True,
                                              value=str(block.untrusted_pixel_mask_path))
    self.runlbock_sizer.Add(self.untrusted_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Calibration folder
    self.calib_dir = gctr.TextButtonCtrl(self.runblock_panel,
                                         label='Calibration:',
                                         label_style='normal',
                                         label_size=(100, -1),
                                         big_button=True,
                                         value=str(block.calib_dir))
    self.runlbock_sizer.Add(self.calib_dir, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark average path (pickle only)
    self.dark_avg_path = gctr.TextButtonCtrl(self.runblock_panel,
                                             label='Dark Average:',
                                             label_style='normal',
                                             label_size=(100, -1),
                                             big_button=True,
                                             value=str(block.dark_avg_path))
    self.runlbock_sizer.Add(self.dark_avg_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark stddev path (pickle only)
    self.dark_stddev_path = gctr.TextButtonCtrl(self.runblock_panel,
                                                label='Dark StdDev:',
                                                label_style='normal',
                                                label_size=(100, -1),
                                                big_button=True,
                                                value=str(block.dark_stddev_path))
    self.runlbock_sizer.Add(self.dark_stddev_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Dark map path (pickle only)
    self.gain_map_path = gctr.TextButtonCtrl(self.runblock_panel,
                                             label='Gain Map:',
                                             label_style='normal',
                                             label_size=(100, -1),
                                             big_button=True,
                                             value=str(block.gain_map_path))
    self.runlbock_sizer.Add(self.gain_map_path, flag=wx.EXPAND | wx.ALL,
                            border=10)

    # Comment
    self.comment = gctr.TextButtonCtrl(self.runblock_panel,
                                       label='Comment:',
                                       label_style='normal',
                                       label_size=(100, -1))
    self.runlbock_sizer.Add(self.comment, flag=wx.EXPAND | wx.ALL,
                            border=10)

    self.main_sizer.Add(self.option_sizer, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.config_panel, flag=wx.EXPAND | wx.ALL, border=10)

    self.runblock_box_sizer.Add(self.runblock_panel)
    self.main_sizer.Add(self.runblock_box_sizer, flag=wx.EXPAND | wx.ALL,
                        border=10)



    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_BUTTON, self.onDarkAvgBrowse,
              id=self.dark_avg_path.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onImportConfig, self.btn_import_cfg)
    self.Bind(wx.EVT_BUTTON, self.onDefaultConfig, self.btn_default_cfg)
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
      self.config.SetValue(cfg_contents)
    cfg_dlg.Destroy()

  def onDefaultConfig(self, e):
    # TODO: Generate default config parameters (re-do based on pickle / CBF)
    pass

  def onOK(self, e):

    startrun_number = int(self.runblocks.start.GetString(self.runblocks.start.GetSelection()))
    startrun = self.db.get_run(run_number=startrun_number).id

    endrun_number = self.runblocks.end.GetString(self.runblocks.end.GetSelection())
    if endrun_number == "None":
      endrun = None
    else:
      endrun = self.db.get_run(run_number=int(endrun_number)).id

    rg_dict = dict(startrun=startrun,
                   endrun=endrun,
                   active=True,
                   config_str=self.config.GetValue(),
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
                   comment=self.comment.ctr.GetValue())
    for key, value in rg_dict.iteritems():
      if str(value) == 'None' or str(value) == '':
        rg_dict[key] = None
      elif type(value) == bool:
        rg_dict[key] = int(value)

    if self.block is None:
      self.block = self.db.create_rungroup(**rg_dict)
      self.parent.trial.add_rungroup(self.block)
    else:
      # if all the parameters are unchanged, do nothing
      all_the_same = [str(rg_dict[key]) == str(getattr(self.block, key)) for key in rg_dict].count(False) == 0
      if not all_the_same:
        # if all the parameters except startrun, endrun and comment are the same,
        # only update those fields
        keep_old_run_group = [str(rg_dict[key]) == str(getattr(self.block, key)) for key in rg_dict \
                              if key not in ['startrun', 'endrun', 'comment']].count(False) == 0
        if keep_old_run_group:
          main = self.parent.parent.GetParent().main
          running = main.job_sentinel is not None and main.job_sentinel.active
          if running:
            main.stop_job_sentinel()

          self.block.startrun = startrun
          self.block.endrun = endrun
          self.block.comment = rg_dict['comment']

          if running:
            main.start_job_sentinel()

        else:
          # enough parameters have changed to warrant creating a new run group
          self.block.active = False
          self.block = self.db.create_rungroup(**rg_dict)
          self.parent.trial.add_rungroup(self.block)

    e.Skip()

  def fill_in_fields(self):
    ''' If previous rungroups exist in trial, fill in fields in nascent block '''
    if len(self.all_blocks) > 0:
      last = self.all_blocks[-1]
      self.config.SetValue(str(last.config_str))
      self.address.ctr.SetValue(str(last.detector_address))
      self.beam_xyz.DetZ.SetValue(str(last.detz_parameter))
      self.beam_xyz.X.SetValue(str(last.beamx))
      self.beam_xyz.Y.SetValue(str(last.beamy))
      self.bin_nrg_gain.binning.SetValue(str(last.binning))
      self.bin_nrg_gain.energy.SetValue(str(last.energy))
      self.bin_nrg_gain.gain_mask_level.SetValue(str(last.gain_mask_level))
      self.untrusted_path.ctr.SetValue(str(last.untrusted_pixel_mask_path))
      self.dark_avg_path.ctr.SetValue(str(last.dark_avg_path))
      self.dark_stddev_path.ctr.SetValue(str(last.dark_stddev_path))
      self.gain_map_path.ctr.SetValue(str(last.gain_map_path))
      self.calib_dir.ctr.SetValue(str(last.calib_dir))
      self.comment.ctr.SetValue(str(last.comment))

  def onImageFormat(self, e):
    self.configure_controls()

  def configure_controls(self):
    sel= self.img_format.ctr.GetString(self.img_format.ctr.GetSelection())
    if 'CBF' in sel:
      self.beam_xyz.X.Disable()
      self.beam_xyz.Y.Disable()
      self.bin_nrg_gain.binning.Disable()
      self.dark_avg_path.Hide()
      self.dark_stddev_path.Hide()
      self.gain_map_path.Hide()
      self.runblock_panel.SetupScrolling()
    elif 'pickle' in sel:
      self.beam_xyz.X.Enable()
      self.beam_xyz.Y.Enable()
      self.bin_nrg_gain.binning.Enable()
      self.dark_avg_path.Show()
      self.dark_stddev_path.Show()
      self.gain_map_path.Show()
      self.runblock_panel.Layout()
      self.runblock_panel.SetupScrolling()

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
    else:
      trial_number = trial.trial

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
                                          ctrl_size=(100, -1),
                                          choices=choices)
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
                               label_size=(120, -1),
                               label_style='bold',
                               ctrl_size=(80, -1),
                               ctrl_value='1.5',
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
