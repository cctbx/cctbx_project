from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/02/2016
Last Changed: 06/02/2016
Description : XFEL UI Initialization module
'''

import os
import wx
import time
from threading import Thread
from wx.lib.scrolledpanel import ScrolledPanel

import xfel.ui.components.xfel_gui_controls as gctr
import xfel.ui.components.xfel_gui_dialogs as dlg
from xfel.ui import load_cached_settings, save_cached_settings

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')

license = 'cctbx.xfel and cctbx.xfel UI are developed under the open source ' \
          'license'

description = 'The cctbx.xfel UI is developed for use during data collection ' \
              'and initial processing at LCSL XFEL beamlines.'

# ------------------------------- Threading ---------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_REFRESH = wx.NewEventType()
EVT_REFRESH = wx.PyEventBinder(tp_EVT_REFRESH, 1)

class RefreshRuns(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class RunSentinel(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def run(self):
    while self.active:
      self.parent.run_window.runs_tab.runs = self.parent.db.get_all_runs()
      self.parent.run_window.runs_tab.all_tags = self.parent.db.get_all_tags()
      evt = RefreshRuns(tp_EVT_REFRESH, -1)
      wx.PostEvent(self.parent.run_window.runs_tab, evt)
      wx.PostEvent(self.parent.run_window.trials_tab, evt)
      time.sleep(1)


# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(800, 500))

    self.params = load_cached_settings()
    self.db = None

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT,
                       label='Quit',
                       bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                       shortHelp='Quit',
                       longHelp='Exit iXFEL')
    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY,
                      label='Run Jobs',
                      bitmap=wx.Bitmap('{}/32x32/play.png'.format(icons)),
                      shortHelp='Run All Jobs',
                      longHelp='Activate all pending jobs')
    self.tb_btn_pause = self.toolbar.AddLabelTool(wx.ID_ANY,
                        label='Pause Jobs',
                        bitmap=wx.Bitmap('{}/32x32/pause.png'.format(icons)),
                        shortHelp='Pause All Jobs',
                        longHelp='Pause all pending jobs')
    self.toolbar.AddSeparator()
    self.tb_btn_options = self.toolbar.AddLabelTool(wx.ID_ANY,
                        label='Settings',
                        bitmap=wx.Bitmap('{}/32x32/settings.png'.format(icons)),
                        shortHelp='Settings',
                        longHelp='Database, user and experiment settings')

    self.toolbar.Realize()

    # Status bar
    self.sb = self.CreateStatusBar()
    self.sb.SetFieldsCount(3)
    self.sb.SetStatusWidths([320, 200, -2])

    # Menu bar
    menubar = wx.MenuBar()
    m_help = wx.Menu()
    self.mb_about = m_help.Append(wx.ID_ANY, '&About')
    menubar.Append(m_help, '&Help')
    self.SetMenuBar(menubar)

    # Place elements in main PRIME window
    main_box = wx.BoxSizer(wx.VERTICAL)

    # Instantiate windows
    self.run_window = RunWindow(self)

    # Single input window
    main_box.Add(self.run_window, 1, flag=wx.ALL | wx.EXPAND, border=10)
    main_box.Add((-1, 20))

    # Menubar button bindings
    self.Bind(wx.EVT_MENU, self.OnAboutBox, self.mb_about)

    # Bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)

    # Draw the main window sizer
    self.SetSizer(main_box)

  def connect_to_db(self):
    from xfel.ui.db import get_db_connection
    from xfel.ui.db.xfel_db import xfel_db_application
    try:
      conn = get_db_connection(self.params)
      print 'Connecting to database...'
    except Exception, e:
      print "Couldn't connect to database"
      return False
    self.db = xfel_db_application(conn, self.params)

    if not self.db.verify_tables():
      self.db.create_tables()
      print 'Creating experiment tables...'
      if not self.db.verify_tables():
        from libtbx.utils import Sorry
        raise Sorry("Couldn't create experiment tables")
        return False

    return True

  def start_sentinels(self):
    self.start_run_sentinel()
    self.start_job_sentinel()

  def stop_sentinels(self):
    self.stop_run_sentinel()
    self.stop_job_sentinel()

  def start_run_sentinel(self):
    self.run_sentinel = RunSentinel(self, active=True)
    self.run_sentinel.start()

  def stop_run_sentinel(self):
    self.run_sentinel.active = False
    self.run_sentinel.join()

  def start_job_sentinel(self):
    pass

  def stop_job_sentinel(self):
    pass

  def OnAboutBox(self, e):
    ''' About dialog '''
    info = wx.AboutDialogInfo()
    info.SetName('iXFEL')
    info.SetLicense(license)
    info.SetDescription(description)
    info.AddDeveloper('Artem Lyubimov')
    info.AddDeveloper('Aaron Brewster')
    info.AddDeveloper('Axel Brunger')
    info.AddDeveloper('Nicholas Sauter')
    wx.AboutBox(info)

  def onRun(self, e):
    ''' All the sentinels will be activated here '''
    pass

  def onQuit(self, e):
    self.stop_sentinels()
    save_cached_settings(self.params)
    self.Destroy()


class RunWindow(wx.Panel):
  ''' Window panel that will house all the run tabs '''

  def __init__(self, parent):
    self.parent = parent
    super(RunWindow, self).__init__(self.parent)

    self.main_panel = wx.Panel(self)
    self.main_nbook = wx.Notebook(self.main_panel, style=0)
    self.runs_tab = RunTab(self.main_nbook, main=self.parent)
    self.trials_tab = TrialsTab(self.main_nbook, main=self.parent)
    self.jobs_tab = JobsTab(self.main_nbook)
    self.status_tab = StatusTab(self.main_nbook)
    self.merge_tab = MergeTab(self.main_nbook)
    self.main_nbook.AddPage(self.runs_tab, 'Runs')
    self.main_nbook.AddPage(self.trials_tab, 'Trials')
    self.main_nbook.AddPage(self.jobs_tab, 'Jobs')
    self.main_nbook.AddPage(self.status_tab, 'Status')
    self.main_nbook.AddPage(self.merge_tab, 'Merge')

    nb_sizer = wx.BoxSizer(wx.VERTICAL)
    nb_sizer.Add(self.main_nbook, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.main_panel.SetSizer(nb_sizer)

    main_sizer = wx.BoxSizer(wx.VERTICAL)
    main_sizer.Add(self.main_panel, 1, flag=wx.EXPAND | wx.ALL, border=3)
    self.SetSizer(main_sizer)


# --------------------------------- UI Tabs ---------------------------------- #

class BaseTab(wx.Panel):
  ''' Base class for runtime tab '''

  def __init__(self, parent):
    wx.Panel.__init__(self, parent=parent, id=wx.ID_ANY, size=(300, 400))

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.main_sizer)


class RunTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)
    self.main = main
    self.last_run = 0
    self.runs = []
    self.all_tags = []
    self.all_tag_buttons = []

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.run_panel = ScrolledPanel(self)
    self.run_sizer = wx.BoxSizer(wx.VERTICAL)
    self.run_panel.SetSizer(self.run_sizer)

    self.colname_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    run_label = wx.StaticText(self, label='Run', size=(60, -1))
    tag_label = wx.StaticText(self, label='Sample Tags', size=(620, -1))
    self.colname_sizer.Add(run_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(tag_label, flag=wx.ALIGN_RIGHT | wx.EXPAND)
    self.colname_sizer.AddGrowableCol(1, 1)
    self.main_sizer.Add(self.colname_sizer, flag=wx.ALL | wx.EXPAND, border=10)

    self.btn_manage_tags = wx.Button(self, label='Manage Tags', size=(120, -1))
    self.main_sizer.Add(self.run_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_manage_tags,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM | wx.ALIGN_RIGHT,
                        border=10)
    self.SetSizer(self.main_sizer)

    # Bindings
    self.Bind(EVT_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_BUTTON, self.onManageTags, self.btn_manage_tags)

    self.refresh_rows()

  def onRefresh(self, e):
    self.refresh_rows()

  def onManageTags(self, e):
    ''' User can add / remove / edit sample tags '''
    mtag_dlg = dlg.TagDialog(self, tags=self.all_tags, db=self.main.db)
    mtag_dlg.Fit()
    mtag_dlg.ShowModal()
    mtag_dlg.Destroy()

    # Update tags on all tag buttons
    for btn in self.all_tag_buttons:
      btn.update_label()

  def refresh_rows(self):

    # Check for new runs and create if necessary
    tag_buttons = [i.run.run_id for i in self.all_tag_buttons]
    for run in self.runs:
      if run.run_id not in tag_buttons:
        self.add_row(run)

    # Update labels on all new tag buttons
    for button in self.all_tag_buttons:
      button.all_tags = self.all_tags
      button.update_label()

    self.run_panel.SetupScrolling()
    self.run_panel.Refresh()

  def add_row(self, run):
    ''' Adds run row to table, matching colname_sizer '''
    row_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    run_no = wx.StaticText(self.run_panel, label=str(run.run),
                           size=(60, -1))
    tag_button = gctr.TagButton(self.run_panel, run=run,
                                all_tags = self.all_tags)
    self.Bind(wx.EVT_BUTTON, self.onTrialButton, id=tag_button.GetId())
    self.all_tag_buttons.append(tag_button)
    row_sizer.Add(run_no, flag=wx.EXPAND | wx.ALIGN_CENTRE)
    row_sizer.Add(tag_button, flag=wx.EXPAND)
    row_sizer.AddGrowableCol(1)
    self.run_sizer.Add(row_sizer, flag=wx.ALL | wx.EXPAND, border=0)

  def onTrialButton(self, e):
    e.GetEventObject().change_tags()


class TrialsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)

    self.trial_number = 1
    self.main = main

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    self.trial_panel = ScrolledPanel(self, size=(300, 350))
    self.trial_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.btn_add_trial = wx.Button(self, label='New Trial',
                                   size=(120, -1))

    self.main_sizer.Add(self.trial_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_add_trial,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM | wx.ALIGN_RIGHT,
                        border=10)

    self.trial_panel.SetSizer(self.trial_sizer)
    self.SetSizer(self.main_sizer)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddTrial, self.btn_add_trial)
    self.Bind(EVT_REFRESH, self.onRefresh)

  def onRefresh(self, e):
    self.db = self.main.db
    self.all_trials = self.db.get_all_trials()
    self.all_runs = self.db.get_all_runs()
    self.first_run = self.all_runs[0]
    self.last_run = self.all_runs[-1]

  def onAddTrial(self, e):
    self.add_new_trial()
    self.trial_panel.Layout()
    self.trial_panel.SetupScrolling()

  def add_new_trial(self):
    self.trial_number = len(self.all_trials) + 1
    self.db.create_trial(trial_id=self.trial_number)
    new_trial = TrialPanel(self.trial_panel,
                           db = self.db,
                           trial=self.db.get_trial(trial_id=self.trial_number),
                           box_label='Trial {}'.format(self.trial_number))
    new_trial.tgl_active.SetValue(True)
    self.trial_sizer.Add(new_trial, flag=wx.EXPAND | wx.ALL, border=10)


class JobsTab(BaseTab):
  def __init__(self, parent):
    BaseTab.__init__(self, parent=parent)


class StatusTab(BaseTab):
  def __init__(self, parent):
    BaseTab.__init__(self, parent=parent)


class MergeTab(BaseTab):
  def __init__(self, parent):
    BaseTab.__init__(self, parent=parent)

# ------------------------------- UI Elements -------------------------------- #

class TrialPanel(wx.Panel):
  def __init__(self, parent, db, trial, box_label=None):
    wx.Panel.__init__(self, parent=parent, size=(200, 300))

    self.db = db
    self.trial = trial
    self.first_run = self.db.get_all_runs()[0].run_id

    trial_box = wx.StaticBox(self, label=box_label)
    self.main_sizer = wx.StaticBoxSizer(trial_box, wx.VERTICAL)

    self.block_panel = ScrolledPanel(self, size=(150, 200))
    self.block_sizer = wx.BoxSizer(wx.VERTICAL)
    self.block_panel.SetSizer(self.block_sizer)
    self.add_panel = wx.Panel(self)
    self.add_sizer = wx.BoxSizer(wx.VERTICAL)
    self.add_panel.SetSizer(self.add_sizer)

    self.add_new_block()

    # Add "New Block" button to a separate sizer (so it is always on bottom)
    self.btn_add_block = wx.Button(self.add_panel, label='New Block',
                                   size=(120, -1))
    self.tgl_active = wx.ToggleButton(self.add_panel, label='De-Activate',
                                      size=(120, -1))
    self.add_sizer.Add(self.tgl_active,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.btn_add_block,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border = 10)

    self.main_sizer.Add(self.block_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.add_panel, flag=wx.ALL | wx.ALIGN_BOTTOM, border=5)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddBlock, self.btn_add_block)
    self.tgl_active.Bind(wx.EVT_TOGGLEBUTTON, self.onToggleActivity)

    self.SetSizer(self.main_sizer)

  def onToggleActivity(self, e):
    if self.tgl_active.GetValue():
      self.trial.active = True
      self.tgl_active.SetLabel('De-Activate')
    else:
      self.trial.active = False
      self.tgl_active.SetLabel('Activate')

  def onAddBlock(self, e):
    ''' Add new block button '''

    # Get list of all buttons in block sizer
    all_children = self.block_sizer.GetChildren()
    run_buttons = [i for i in all_children if isinstance(i.GetWindow(),
                                                         gctr.RunBlockButton)]

    # Get last button in block sizer, update label with latest run
    self.last_run = self.db.get_all_runs()[-1]
    last_block = run_buttons[-1].GetWindow()
    last_block.last_run = self.last_run.run_id
    last_block.update_label()
    self.first_run = self.last_run.run_id + 1

    self.add_new_block()
    self.block_panel.Layout()
    self.block_panel.SetupScrolling()

  def add_new_block(self):
    ''' Add new run block button '''
    new_block = gctr.RunBlockButton(self.block_panel, size=(120, -1))
    new_block.first_run = self.first_run
    new_block.update_label()
    self.Bind(wx.EVT_BUTTON, self.onRunBlockOptions, id=new_block.GetId())
    self.block_sizer.Add(new_block,
                         flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                         border=5)

  def onRunBlockOptions(self, e):
    ''' Open dialog and change run_block options '''
    run_block = e.GetEventObject()
    label = run_block.block_label
    first_run = run_block.first_run
    last_run = run_block.last_run

    rblock_dlg = dlg.RunBlockDialog(self, self.db)
    rblock_dlg.Fit()



    rblock_dlg.ShowModal()
