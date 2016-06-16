from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/04/2016
Last Changed: 06/04/2016
Description : XFEL UI Custom Dialogs
'''

import wx
import os
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

class BaseDialog(wx.Dialog):
  def __init__(self, parent,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):
    wx.Dialog.__init__(self, parent, *args, **kwargs)

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)

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

    self.btn_mp = wx.Button(self, label='Multiprocessing...')
    self.btn_op = wx.Button(self, label='Advanced Settings...')
    self.btn_OK = wx.Button(self, label="OK", id=wx.ID_OK)
    self.btn_cancel = wx.Button(self, label="Cancel", id=wx.ID_CANCEL)

    button_sizer = wx.FlexGridSizer(1, 5, 0, 10)
    button_sizer.AddMany([(self.btn_mp),
                          (self.btn_op),
                          (0,0),
                          (self.btn_OK),
                          (self.btn_cancel)])

    button_sizer.AddGrowableCol(2)
    self.main_sizer.Add(button_sizer,
                        flag=wx.EXPAND | wx.ALL,
                        border=10)
    self.SetSizer(self.main_sizer)

    self.Bind(wx.EVT_BUTTON, self.onDBCredentialsButton, id=self.db_cred.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onOK, id=self.btn_OK.GetId())
    self.Bind(wx.EVT_BUTTON, self.onBrowse, id=self.output.btn_big.GetId())

  def onBrowse(self, e):
    dlg = wx.DirDialog(self, "Choose the input directory:",
                       style=wx.DD_DEFAULT_STYLE)

    if dlg.ShowModal() == wx.ID_OK:
      self.output.ctr.SetValue(dlg.GetPath())
    dlg.Destroy()

  def onDBCredentialsButton(self, e):
    creds = DBCredentialsDialog(self, self.params)
    creds.Center()
    if (creds.ShowModal() == wx.ID_OK):
      self.params.db.host     = creds.db_host.ctr.GetValue()
      self.params.db.user     = creds.db_user.ctr.GetValue()
      self.params.db.password = creds.db_password.ctr.GetValue()
      self.params.web.user     = creds.web_user.ctr.GetValue()
      self.params.web.password = creds.web_password.ctr.GetValue()

  def onOK(self, e):
    self.params.experiment_tag = self.db_cred.ctr.GetValue()
    self.params.experiment = self.params.db.name = self.experiment.ctr.GetValue()
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
                        size=(600, 350),
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
    self.main_sizer.Add(self.db_password, flag=wx.EXPAND | wx.ALL,border=10)
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

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
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
        self.tag_list.InsertStringItem(self.index, tag.name)
        self.tag_list.SetStringItem(self.index, 1, tag.comment)
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

    self.SetSizer(self.main_sizer)
    self.Layout()

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onAdd, self.btn_add)
    self.Bind(wx.EVT_BUTTON, self.onRemove, self.btn_rmv)
    self.Bind(wx.EVT_BUTTON, self.onClear, self.btn_clr)
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

  def onAdd(self, e):
    ''' Add a string item to list; focus on item & provide default tag name'''
    new_tag = ('default tag {}'.format(self.index), '')
    self.new_tags.append(new_tag)
    self.tag_list.InsertStringItem(self.index, new_tag[0])
    self.tag_list.SetStringItem(self.index, 1, new_tag[1])
    self.tag_list.SetItemData(self.index, -1)
    self.tag_list.Select(self.index)
    self.tag_list.Focus(self.index)
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
      self.tag_list.DeleteAllItems()
      self.index = 0
      self.db_tags = []

  def onOK(self, e):
    count = self.tag_list.GetItemCount()
    if count == 0:
      warning = wx.MessageDialog(self,
                                 message='Must have at least one tag!',
                                 caption='Warning',
                                 style=wx.OK | wx.ICON_EXCLAMATION)
      warning.ShowModal()
    else:
      try:
        # Update names for edited tags
        all_items = [(self.tag_list.GetItemData(i),
                    self.tag_list.GetItem(itemId=i, col=0),
                    self.tag_list.GetItem(itemId=i, col=1))
                   for i in range(count)]
        edited_items = [i for i in all_items if i[0] != -1]
        for tag in self.db_tags:
          for item in edited_items:
            if tag.tag_id == item[0]:
              if tag.name != item[1].m_text:
                tag.name = item[1].m_text
              if tag.comment != item[2].m_text:
                tag.comment = item[2].m_text

        # Delete tags from DB
        for tag in self.deleted_tags:
          self.db.delete_tag(tag)

        # Add new tags to DB
        for tag in self.new_tags:
          self.db.create_tag(name=tag[0], comment=tag[1])
      except Exception:
        pass

      e.Skip()

class RunBlockDialog(BaseDialog):
  ''' Comes up when individual run block button is clicked; allows for run
  block settings to be manipulated by user '''

  def __init__(self, parent, block,
               label_style='bold',
               content_style='normal',
               *args, **kwargs):

    self.block = block
    db = block.app
    self.first_run = db.get_run(run_id=block.startrun).run
    if block.endrun is None:
      self.last_run = None
    else:
      self.last_run = db.get_run(run_id=block.endrun).run

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 200),
                        *args, **kwargs)

    # Run block start / end points (choice widgets)

    start = [str(i.run) for i in db.get_all_runs()]
    stop = start + ['None']
    firstidx = start.index(str(self.first_run))
    lastidx = stop.index(str(self.last_run))

    self.runblocks = gctr.MultiChoiceCtrl(self,
                                          label='Run block:',
                                          label_style='bold',
                                          label_size=(100, -1),
                                          ctrl_size=(150, -1),
                                          items={'start':start, 'end':stop})
    self.runblocks.start.SetSelection(firstidx)
    self.runblocks.end.SetSelection(lastidx)
    self.main_sizer.Add(self.runblocks, flag=wx.EXPAND | wx.ALL, border=10)

    # Beam XYZ
    self.beam_xyz = gctr.OptionCtrl(self,
                                    label='Beam:',
                                    label_style='bold',
                                    label_size=(100, -1),
                                    ctrl_size=(100, -1),
                                    items={'X':self.block.beamx,
                                           'Y':self.block.beamy,
                                           'DetZ':self.block.detz_parameter})
    self.main_sizer.Add(self.beam_xyz, flag=wx.EXPAND | wx.ALL, border=10)

    # Dark path
    self.dark_path = gctr.TextButtonCtrl(self,
                                         label='Dark CBF Path:',
                                         label_style='bold',
                                         label_size=(100, -1),
                                         big_button=True,
                          value=str(self.block.dark_avg_path))
    self.main_sizer.Add(self.dark_path, flag=wx.EXPAND | wx.ALL, border=10)

    # Dialog control
    dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Bind(wx.EVT_BUTTON, self.onDarkBrowse,
              id=self.dark_path.btn_big.GetId())
    self.Bind(wx.EVT_BUTTON, self.onOK, id=wx.ID_OK)

    self.Layout()

  def onOK(self, e):
    startrun_number = int(self.runblocks.start.GetString(self.runblocks.start.GetSelection()))
    startrun = self.block.app.get_run(run_number=startrun_number).id

    endrun_number = self.runblocks.end.GetString(self.runblocks.end.GetSelection())
    if endrun_number == "None":
      endrun = None
    else:
      endrun = self.block.app.get_run(run_number=int(endrun_number)).id

    self.block.startrun = startrun
    self.block.endrun = endrun
    e.Skip()

  def onDarkBrowse(self, e):
    dark_dlg = wx.FileDialog(self,
                             message="Load dark file",
                             defaultDir=os.curdir,
                             defaultFile="*.cbf",
                             wildcard="*.cbf",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )

    if dark_dlg.ShowModal() == wx.ID_OK:
      self.dark_file = dark_dlg.GetPaths()[0]
    dark_dlg.Destroy()

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

    BaseDialog.__init__(self, parent,
                        label_style=label_style,
                        content_style=content_style,
                        size=(600, 600),
                        style=wx.RESIZE_BORDER,
                        *args, **kwargs)

    if trial is None:
      trials = [t.trial for t in db.get_all_trials()]
      if len(trials) == 0:
        trial_number = 0
      else:
        trial_number = max(trials) + 1
    else:
      trial_number = trial.trial

    self.trial_number = gctr.TextButtonCtrl(self,
                                            label='Trial number:',
                                            label_size=(150, -1),
                                            label_style='bold',
                                            big_button=True,
                                            big_button_label='Default PHIL',
                                            value="{}".format(trial_number))
    self.trial_comment = gctr.TextButtonCtrl(self,
                                             label='Comment:',
                                             label_size=(150, -1),
                                             label_style='bold')
    self.trial_phil = gctr.TextButtonCtrl(self,
                                          label='PHIL Path:',
                                          label_size=(150, -1),
                                          label_style='bold',
                                          big_button=True)
    self.phil_box = wx.TextCtrl(self, style=wx.TE_MULTILINE)


    # TODO: show trial's PHIL blob & inactivate everything
    # if not self.new:
    #   self.trial_number.ctr.SetStyle(wx.TE_READONLY)
    #   self.trial_number.ctr.Disable()
    #   self.trial_comment.ctr.SetStyle(wx.TE_READONLY)
    #   self.trial_phil.ctr.SetStyle(wx.TE_READONLY)
    #   self.phil_box.SetStyle(wx.TE_READONLY)

    self.main_sizer.Add(self.trial_number,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.phil_box, 1,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.trial_phil,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)
    self.main_sizer.Add(self.trial_comment,
                        flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT,
                        border=10)

    # Dialog control
    if self.new:
      dialog_box = self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL)
    else:
      dialog_box = self.CreateSeparatedButtonSizer(wx.OK)
    self.main_sizer.Add(dialog_box,
                        flag=wx.EXPAND | wx.ALIGN_RIGHT | wx.ALL,
                        border=10)

    self.Layout()

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onBrowse, self.trial_phil.btn_big)
    self.Bind(wx.EVT_BUTTON, self.onDefault, self.trial_number.btn_big)
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

  def onReadIn(self, e):
    ''' When PHIL file is selected, check for
          a) existence,
          b) readability,
        and populate the text box with the PHIL parameters '''

    pass

  def onDefault(self, e):
    # TODO: Generate default PHIL parameters
    pass

  def onOK(self, e):
    target_phil_str = self.phil_box.GetValue()
    comment = self.trial_comment.ctr.GetValue()

    if self.trial is None:
      self.db.create_trial(
        trial = int(self.trial_number.ctr.GetValue()),
        active = True,
        target_phil_str = target_phil_str,
        comment = comment)
    else:
      self.trial.target_phil_str = target_phil_str
      self.trial.comment = comment
    e.Skip()
