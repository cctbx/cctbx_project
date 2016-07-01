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
import numpy as np
from threading import Thread
from wx.lib.scrolledpanel import ScrolledPanel
from libtbx import easy_run

try:
  from MySQLdb import OperationalError
except ImportError:
  from libtbx.utils import Sorry
  raise Sorry('Mysql not available')

from xfel.clustering.cluster import Cluster
import xfel.ui.components.xfel_gui_controls as gctr
import xfel.ui.components.xfel_gui_dialogs as dlg
import xfel.ui.components.xfel_gui_plotter as pltr

from xfel.ui import load_cached_settings, save_cached_settings
from xfel.ui.db import get_run_path
from xfel.ui.db.xfel_db import xfel_db_application

from prime.postrefine.mod_gui_init import PRIMEInputWindow, PRIMERunWindow
from prime.postrefine.mod_input import master_phil
from iota.components import iota_misc as misc

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')
user = os.getlogin()

license = 'cctbx.xfel and cctbx.xfel UI are developed under the open source ' \
          'license'

description = 'The cctbx.xfel UI is developed for use during data collection ' \
              'and initial processing at LCSL XFEL beamlines.'


# ------------------------------- Run Sentinel ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_RUN_REFRESH = wx.NewEventType()
EVT_RUN_REFRESH = wx.PyEventBinder(tp_EVT_RUN_REFRESH, 1)

class RefreshRuns(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class RunSentinel(Thread):
  ''' Worker thread for runs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    evt = RefreshRuns(tp_EVT_RUN_REFRESH, -1)
    wx.PostEvent(self.parent.run_window.runs_tab, evt)
    wx.PostEvent(self.parent.run_window.trials_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()

    while self.active:
      try:
        db = xfel_db_application(self.parent.params)
        self.parent.run_window.run_light.change_status('on')
      except OperationalError:
        self.parent.run_window.run_light.change_status('alert')

      # Find the delta
      known_runs = [r.run for r in db.get_all_runs()]
      unknown_runs = [run['run'] for run in db.list_lcls_runs() if
                      run['run'] not in known_runs]

      if len(unknown_runs) > 0:
        for run_number in unknown_runs:
          db.create_run(run = run_number)
        print "%d new runs" % len(unknown_runs)
        self.post_refresh()
      else:
        pass #print "No new data..."
      time.sleep(1)

# ------------------------------- Job Monitor ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_JOB_MONITOR = wx.NewEventType()
EVT_JOB_MONITOR = wx.PyEventBinder(tp_EVT_JOB_MONITOR, 1)

class MonitorJobs(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class JobMonitor(Thread):
  ''' Monitor thread for jobs; generated so that the GUI does not lock up when
      monitoring is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    evt = MonitorJobs(tp_EVT_JOB_MONITOR, -1)
    wx.PostEvent(self.parent.run_window.jobs_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()

    #from xfel.ui.db.job import submit_all_jobs
    #try:
    #  db = xfel_db_application(self.parent.params)
    #  self.parent.run_window.job_light.change_status('on')
    #except OperationalError:
    #  self.parent.run_window.job_light.change_status('alert')

    while self.active:
      self.post_refresh()
      time.sleep(1)


# ------------------------------- Job Sentinel ------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_JOB_REFRESH = wx.NewEventType()
EVT_JOB_REFRESH = wx.PyEventBinder(tp_EVT_JOB_REFRESH, 1)

class RefreshJobs(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)

class JobSentinel(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active

  def post_refresh(self):
    pass

  def run(self):
    # one time post for an initial update
    self.post_refresh()

    from xfel.ui.db.job import submit_all_jobs
    try:
      db = xfel_db_application(self.parent.params)
      self.parent.run_window.job_light.change_status('on')
    except OperationalError:
      self.parent.run_window.job_light.change_status('alert')

    while self.active:
      submit_all_jobs(db)
      self.post_refresh()
      time.sleep(1)

# ----------------------------- Progress Sentinel ---------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_PRG_REFRESH = wx.NewEventType()
EVT_PRG_REFRESH = wx.PyEventBinder(tp_EVT_PRG_REFRESH, 1)

class RefreshStats(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class ProgressSentinel(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               active=True):
    Thread.__init__(self)
    self.parent = parent
    self.active = active
    self.output = self.parent.params.output_folder
    self.number_of_pickles = 0
    self.info = {}

  def post_refresh(self):
    evt = RefreshStats(tp_EVT_PRG_REFRESH, -1, self.info)
    wx.PostEvent(self.parent.run_window.status_tab, evt)

  def run(self):
    # one time post for an initial update
    self.post_refresh()

    while self.active:
      try:
        db = xfel_db_application(self.parent.params)
        self.parent.run_window.prg_light.change_status('idle')
      except OperationalError:
        self.parent.run_window.prg_light.change_status('alert')
        time.sleep(5)
        continue

      if len(db.get_all_trials()) > 0:
        trial = db.get_trial(
          trial_number=self.parent.run_window.status_tab.trial_no)

        tags = self.parent.run_window.status_tab.selected_tags
        cells = db.get_stats(trial=trial, tags=tags)()

        if self.parent.run_window.status_tab.tag_trial_changed:
          self.parent.run_window.status_tab.redraw_windows = True
          self.parent.run_window.status_tab.tag_trial_changed = False

        for cell in cells:

          # Check to see if different high-res bin has been selected
          current_rows = self.parent.run_window.status_tab.rows
          if current_rows != {}:
            bins = cell.bins[:int(current_rows[cell.isoform.name]['high_bin'])]
          else:
            bins = cell.bins

          # Check for cell isoform
          if cell.isoform is None:
            self.info[cells.index(cell)] = {'a':cell.cell_a,
                                            'b':cell.cell_b,
                                            'c':cell.cell_c,
                                            'alpha':cell.cell_alpha,
                                            'beta':cell.cell_beta,
                                            'gamma':cell.cell_gamma}
          else:
            counts = [int(i.count) for i in bins]
            totals = [int(i.total_hkl) for i in bins]

            # Apply throttle to multiplicity calculation
            if trial.process_percent is None:
              process_percent = 100
            else:
              process_percent = trial.process_percent

            # Generate multiplicity graph for isoforms
            mult = sum(counts) / sum(totals) * (process_percent / 100)
            self.info[cell.isoform.name] = {'multiplicity':mult,
                                            'bins':bins,
                                            'isoform':cell.isoform.name,
                                            'a':cell.cell_a,
                                            'b':cell.cell_b,
                                            'c':cell.cell_c,
                                            'alpha':cell.cell_alpha,
                                            'beta':cell.cell_beta,
                                            'gamma':cell.cell_gamma}
      self.post_refresh()
      self.info = {}
      self.parent.run_window.prg_light.change_status('on')
      time.sleep(5)

# ------------------------------- Clustering --------------------------------- #

# Set up events for and for finishing all cycles
tp_EVT_CLUSTERING = wx.NewEventType()
EVT_CLUSTERING = wx.PyEventBinder(tp_EVT_CLUSTERING, 1)

class ClusteringResult(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid, result=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.result = result
  def GetValue(self):
    return self.result

class Clusterer():
  def __init__(self, trial, runblocks, output, sample_size, threshold,
               histogram=False):
    self.trial = trial
    self.runblocks = runblocks
    self.output = output
    self.sample_size = sample_size
    self.threshold = threshold
    self.histogram = histogram

  def unit_cell_clustering(self):

    # 1. Get all pickle files, check if new ones arrived
    run_numbers = []
    rb_paths = []
    for rb in self.runblocks:
      for run in rb.runs:
        if run.run not in run_numbers:
          run_numbers.append(run.run)
          rb_paths.append(os.path.join(get_run_path(self.output, self.trial,
                                                    rb, run), "out"))
    all_pickles = []
    for path in rb_paths:
      try:
        pickles = [os.path.join(path, i) for i in os.listdir(path) if
                   i.endswith('pickle') and 'int-' in i]
        all_pickles = all_pickles + pickles
      except OSError, error:
        print 'Folder not found!'
        print error

    if len(all_pickles) == 0:
      print 'No images integrated (yet)'
      return

    if self.histogram:
      # If histogram button was pressed, return set of unit cells
      # 2. Load pickle files, extract and return unit cell tuples
      from libtbx import easy_pickle as ep
      all_cells = []
      for pickle in all_pickles:
        info = ep.load(pickle)
        cell = info['observations'][0].unit_cell().parameters()
        all_cells.append(cell)

      return all_cells

    else:
      # If clustering button was pressed, do clustering
      # 2. Pick subset
      print len(all_pickles), self.sample_size
      subset = list(np.random.choice(all_pickles, size=int(self.sample_size)))

      # 3. Do clustering
      ucs = Cluster.from_files(subset, use_b=True)
      clusters, _ = ucs.ab_cluster(self.threshold,
                                   log=False, write_file_lists=False,
                                   schnell=False, doplot=False)

      return clusters

class ClusteringWorker(Thread):
  ''' Worker thread for jobs; generated so that the GUI does not lock up when
      processing is running '''

  def __init__(self,
               parent,
               trial,
               runblocks,
               output,
               sample_size=1000,
               threshold=250,):
    Thread.__init__(self)
    self.parent = parent
    self.trial = trial
    self.runblocks = runblocks
    self.output = output
    self.sample_size = sample_size
    self.threshold = threshold

  def run(self):
    clusterer = Clusterer(self.trial, self.runblocks, self.output,
                          self.sample_size, self.threshold)
    self.clusters = clusterer.unit_cell_clustering()

    evt = ClusteringResult(tp_EVT_CLUSTERING, -1, self.clusters)
    wx.PostEvent(self.parent, evt)

  def histogram(self):
    clusterer = Clusterer(self.trial, self.runblocks, self.output,
                          self.sample_size, self.threshold, histogram=True)

    self.all_cells = clusterer.unit_cell_clustering()

    evt = ClusteringResult(tp_EVT_CLUSTERING, -1, self.all_cells)
    wx.PostEvent(self.parent, evt)

# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(800, 500))

    self.run_sentinel = None
    self.job_sentinel = None
    self.prg_sentinel = None

    self.params = load_cached_settings()
    self.db = None

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT,
                       label='Quit',
                       bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                       shortHelp='Quit',
                       longHelp='Exit CCTBX.XFEL')
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
    self.tb_btn_calibrate = self.toolbar.AddLabelTool(wx.ID_ANY,
                        label='Calibration',
                        bitmap=wx.Bitmap('{}/32x32/calib.png'.format(icons)),
                        shortHelp='Calibration',
                        longHelp='Detector geometry calibration')
    self.toolbar.AddSeparator()
    self.tb_btn_settings = self.toolbar.AddLabelTool(wx.ID_ANY,
                        label='Settings',
                        bitmap=wx.Bitmap('{}/32x32/settings.png'.format(icons)),
                        shortHelp='Settings',
                        longHelp='Database, user and experiment settings')

    self.toolbar.Realize()
    self.toolbar.EnableTool(self.tb_btn_pause.GetId(), False)

    # Status bar
    self.sb = self.CreateStatusBar()

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
    self.Bind(wx.EVT_TOOL, self.onPause, self.tb_btn_pause)
    self.Bind(wx.EVT_TOOL, self.onCalibration, self.tb_btn_calibrate)
    self.Bind(wx.EVT_TOOL, self.onSettings, self.tb_btn_settings)

    # Draw the main window sizer
    self.SetSizer(main_box)

  def connect_to_db(self, drop_tables = False):
    self.db = xfel_db_application(self.params)

    if drop_tables:
      self.db.drop_tables()

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
    self.start_job_monitor()
    self.start_prg_sentinel()

  def stop_sentinels(self):
    self.stop_run_sentinel()
    self.stop_job_sentinel()
    self.stop_job_monitor()
    self.stop_prg_sentinel()

  def start_run_sentinel(self):
    self.run_sentinel = RunSentinel(self, active=True)
    self.run_sentinel.start()
    self.run_window.run_light.change_status('on')

  def stop_run_sentinel(self):
    self.run_window.run_light.change_status('off')
    self.run_sentinel.active = False
    self.run_sentinel.join()

  def start_job_monitor(self):
    self.job_monitor = JobMonitor(self, active=True)
    self.job_monitor.start()

  def stop_job_monitor(self):
    self.job_monitor.active = False
    self.job_monitor.join()

  def start_job_sentinel(self):
    self.job_sentinel = JobSentinel(self, active=True)
    self.job_sentinel.start()
    self.run_window.job_light.change_status('on')

    self.toolbar.EnableTool(self.tb_btn_run.GetId(), False)
    self.toolbar.EnableTool(self.tb_btn_pause.GetId(), True)

  def stop_job_sentinel(self):
    if self.job_sentinel is not None:
      self.run_window.job_light.change_status('off')
      self.job_sentinel.active = False
      self.job_sentinel.join()

    self.toolbar.EnableTool(self.tb_btn_run.GetId(), True)
    self.toolbar.EnableTool(self.tb_btn_pause.GetId(), False)

  def start_prg_sentinel(self):
    self.prg_sentinel = ProgressSentinel(self, active=True)
    self.prg_sentinel.start()
    self.run_window.prg_light.change_status('on')

  def stop_prg_sentinel(self):
    self.run_window.prg_light.change_status('off')
    self.prg_sentinel.active = False
    self.prg_sentinel.join()

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

  def onSettings(self, e):
    settings_dlg = dlg.SettingsDialog(self,
                                      params=self.params)
    settings_dlg.db_cred.btn_big.Disable()
    settings_dlg.SetTitle('Settings')

    if (settings_dlg.ShowModal() == wx.ID_OK):
      self.params = settings_dlg.params
      self.title = 'CCTBX.XFEL | {} | {}'.format(self.params.experiment_tag,
                                                 self.params.experiment)

  def onCalibration(self, e):
    calib_dlg = dlg.CalibrationDialog(self, db=self.db)
    calib_dlg.Fit()

    calib_dlg.ShowModal()

  def onRun(self, e):
    ''' All the jobs will be activated here '''
    self.start_job_sentinel()

  def onPause(self, e):
    ''' Pause all jobs '''
    self.stop_job_sentinel()

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
    self.jobs_tab = JobsTab(self.main_nbook, main=self.parent)
    self.status_tab = StatusTab(self.main_nbook, main=self.parent)
    self.merge_tab = MergeTab(self.main_nbook)
    self.main_nbook.AddPage(self.runs_tab, 'Runs')
    self.main_nbook.AddPage(self.trials_tab, 'Trials')
    self.main_nbook.AddPage(self.jobs_tab, 'Jobs')
    self.main_nbook.AddPage(self.status_tab, 'Status')
    self.main_nbook.AddPage(self.merge_tab, 'Merge')

    self.sentinel_box = wx.BoxSizer(wx.HORIZONTAL)
    self.run_light = gctr.SentinelStatus(self.main_panel, label='Run Sentinel')
    self.job_light = gctr.SentinelStatus(self.main_panel, label='Job Sentinel')
    self.prg_light = gctr.SentinelStatus(self.main_panel, label='Progress Sentinel')
    self.sentinel_box.Add(self.run_light, flag=wx.LEFT, border=10)
    self.sentinel_box.Add(self.job_light, flag=wx.LEFT, border=20)
    self.sentinel_box.Add(self.prg_light, flag=wx.LEFT, border=20)

    nb_sizer = wx.BoxSizer(wx.VERTICAL)
    nb_sizer.Add(self.main_nbook, 1, flag=wx.EXPAND | wx.ALL, border=3)
    nb_sizer.Add(self.sentinel_box, flag=wx.EXPAND | wx.ALL, border=3)
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
    self.all_runs = []
    self.all_tags = []
    self.all_tag_buttons = []

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

    # Bindings
    self.Bind(EVT_RUN_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_BUTTON, self.onManageTags, self.btn_manage_tags)

  def onRefresh(self, e):
      self.refresh_rows()

  def onManageTags(self, e):
    ''' User can add / remove / edit sample tags '''
    mtag_dlg = dlg.TagDialog(self, db=self.main.db)
    mtag_dlg.Fit()
    mtag_dlg.ShowModal()
    mtag_dlg.Destroy()

    # Update tags on all tag buttons
    for btn in self.all_tag_buttons:
      btn.tags = btn.run.tags
      btn.update_label()

  def refresh_rows(self, all=False):

    # Get new runs
    old_run_numbers = [run.run for run in self.all_runs]
    all_runs = self.main.db.get_all_runs()
    new_runs = [run for run in all_runs if run.run not in old_run_numbers]

    # Update either all or only new runs
    if all:
      runs = self.all_runs
    else:
      runs = new_runs
    for run in runs:
      self.add_row(run)

    self.all_runs = all_runs

    # Update labels on all new tag buttons
    self.all_tags = self.main.db.get_all_tags()
    for button in self.all_tag_buttons:
      button.all_tags = self.all_tags
      button.update_label()

    self.run_panel.SetupScrolling()
    self.run_panel.Refresh()

  def add_row(self, run):
    ''' Adds run row to table, matching colname_sizer '''
    run_row = RunEntry(self.run_panel, run=run, params=self.main.params)
    self.all_tag_buttons.append(run_row.tag_button)
    self.run_sizer.Add(run_row, flag=wx.ALL | wx.EXPAND, border=0)


class TrialsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)

    self.main = main
    self.show_active_only = False

    self.trial_panel = ScrolledPanel(self, size=(300, 350))
    self.trial_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.trial_panel.SetSizer(self.trial_sizer)

    self.btn_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.btn_sizer.AddGrowableCol(0)
    self.btn_add_trial = wx.Button(self, label='New Trial', size=(120, -1))
    self.btn_active_only = wx.ToggleButton(self,
                                           label='Show Only Active Trials',
                                    size=(180, self.btn_add_trial.GetSize()[1]))
    self.btn_sizer.Add(self.btn_active_only, flag=wx.ALIGN_RIGHT)
    self.btn_sizer.Add(self.btn_add_trial)

    self.main_sizer.Add(self.trial_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_sizer, flag=wx.EXPAND | wx.ALL, border=10)

    # Bindings
    self.Bind(EVT_RUN_REFRESH, self.onRefresh)
    self.Bind(wx.EVT_BUTTON, self.onAddTrial, self.btn_add_trial)
    self.Bind(wx.EVT_TOGGLEBUTTON, self.onActiveOnly, self.btn_active_only)

  def onRefresh(self, e):
    self.refresh_trials()

  def refresh_trials(self):
    self.trial_sizer.Clear(deleteWindows=True)
    self.all_trials = self.main.db.get_all_trials()
    for trial in self.all_trials:
      if self.show_active_only:
        if trial.active:
          self.add_trial(trial=trial)
      else:
        self.add_trial(trial=trial)

    self.trial_panel.SetSizer(self.trial_sizer)
    self.trial_panel.Layout()
    self.trial_panel.SetupScrolling()

  def add_trial(self, trial):
    new_trial = TrialPanel(self.trial_panel,
                           db=self.main.db,
                           trial=trial,
                           box_label='Trial {}'.format(trial.trial))
    new_trial.chk_active.SetValue(trial.active)
    new_trial.refresh_trial()
    self.trial_sizer.Add(new_trial, flag=wx.EXPAND | wx.ALL, border=10)

  def onAddTrial(self, e):
    new_trial_dlg = dlg.TrialDialog(self, db=self.main.db)
    new_trial_dlg.Fit()

    if new_trial_dlg.ShowModal() == wx.ID_OK:
      self.refresh_trials()

  def onActiveOnly(self, e):
    self.show_active_only = self.btn_active_only.GetValue()
    self.refresh_trials()


class JobsTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)

    self.main = main
    self.db = None
    self.all_trials = []
    self.filter = 'All jobs'

    self.job_panel = ScrolledPanel(self)
    self.job_sizer = wx.BoxSizer(wx.VERTICAL)
    self.job_panel.SetSizer(self.job_sizer)

    self.colname_sizer = wx.FlexGridSizer(1, 3, 0, 10)
    trial_label = wx.StaticText(self, label='Trial', size=(60, -1))
    run_label = wx.StaticText(self, label='Run', size=(60, -1))
    tag_label = wx.StaticText(self, label='Status', size=(560, -1))
    self.colname_sizer.Add(trial_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(run_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(tag_label, flag=wx.ALIGN_RIGHT | wx.EXPAND)
    self.colname_sizer.AddGrowableCol(2, 1)
    self.main_sizer.Add(self.colname_sizer, flag=wx.ALL | wx.EXPAND, border=10)

    self.trial_choice = gctr.ChoiceCtrl(self,
                                        label='Filter by:',
                                        label_size=(80, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])
    self.btn_kill_all = wx.Button(self, label='Kill All', size=(120, -1))
    self.option_sizer = wx.FlexGridSizer(1, 2, 0, 20)
    self.option_sizer.AddMany([(self.trial_choice), (self.btn_kill_all)])

    self.main_sizer.Add(self.job_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.option_sizer, flag=wx.EXPAND | wx.ALL, border=10)


    self.Bind(wx.EVT_BUTTON, self.onKillAll, self.btn_kill_all)
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_choice.ctr)
    self.Bind(EVT_JOB_MONITOR, self.onMonitorJobs)

  def onTrialChoice(self, e):
    self.filter = self.trial_choice.ctr.GetString(
                  self.trial_choice.ctr.GetSelection())

  def onKillAll(self, e):
    # TODO: make method to kill a job, apply to all jobs here
    pass

  def onMonitorJobs(self, e):
    self.db = xfel_db_application(self.main.params)

    # Find new trials
    if self.db is not None:
      all_db_trials = [str(i.trial) for i in self.db.get_all_trials()]
      new_trials = [i for i in all_db_trials if i not in self.all_trials]
      if len(new_trials) > 0:
        self.find_trials()
        self.all_trials = all_db_trials

    self.monitor_jobs()
    self.job_panel.SetupScrolling()
    self.job_panel.Refresh()

  def find_trials(self):
    if self.db is not None:
      choices = ['All jobs'] + \
                ['trial {}'.format(i.trial) for i in self.db.get_all_trials()]
      self.trial_choice.ctr.Clear()
      for choice in choices:
        self.trial_choice.ctr.Append(choice)

      if self.filter == 'All jobs':
        self.trial_choice.ctr.SetSelection(0)
      else:
        self.trial_choice.ctr.SetSelection(int(self.filter[-1]))

  def monitor_jobs(self):
    if self.db is not None:
      jobs = self.db.get_all_jobs()
      if str(self.filter).lower() != 'all jobs':
        jobs = [i for i in jobs if i.trial_id == int(self.filter[-1])]

      self.job_sizer.DeleteWindows()
      for job in jobs:
        self.add_job_row(job)

  def add_job_row(self, job):
    job_row_sizer = wx.FlexGridSizer(1, 3, 0, 10)
    job_row_sizer.Add(wx.StaticText(self.job_panel, label=str(job.trial.trial),
                                    size=(60, -1)))
    job_row_sizer.Add(wx.StaticText(self.job_panel, label=str(job.run.run),
                                    size=(60, -1)))
    job_row_sizer.Add(wx.StaticText(self.job_panel, label=str(job.status),
                                    size=(560, -1)))
    job_row_sizer.AddGrowableCol(2, 1)
    self.job_sizer.Add(job_row_sizer, flag=wx.EXPAND)


class StatusTab(BaseTab):
  def __init__(self, parent, main):
    BaseTab.__init__(self, parent=parent)

    self.main = main
    self.trial_no = 0
    self.tags = None
    self.all_trials = []
    self.all_tags = []
    self.selected_tags = []
    self.rows = {}
    self.tag_trial_changed = False
    self.redraw_windows = True
    self.multiplicity_goal = 10

    self.status_panel = ScrolledPanel(self, size=(900, 120))
    self.status_box = wx.StaticBox(self.status_panel, label='Data Statistics')
    self.status_sizer = wx.StaticBoxSizer(self.status_box, wx.VERTICAL)
    self.status_panel.SetSizer(self.status_sizer)

    self.iso_panel = ScrolledPanel(self, size=(-1, 100))
    self.iso_box = wx.StaticBox(self.iso_panel, label='Isoforms')
    self.iso_box_sizer = wx.StaticBoxSizer(self.iso_box, wx.VERTICAL)
    self.iso_panel.SetSizer(self.iso_box_sizer)

    self.trial_number = gctr.ChoiceCtrl(self,
                                        label='Trial:',
                                        label_size=(120, -1),
                                        label_style='normal',
                                        ctrl_size=(100, -1),
                                        choices=[])
    self.tag_list = gctr.CheckListCtrl(self,
                                       ctrl_size=(200, 100),
                                       choices=[],
                                       direction='vertical')

    self.opt_multi = gctr.OptionCtrl(self,
                                     label='Goal multiplicity:',
                                     label_size=(120, -1),
                                     ctrl_size=(100, -1),
                                     items=[('goal', 10)])

    self.bottom_sizer = wx.FlexGridSizer(1, 2, 0, 10)

    multi_box = wx.StaticBox(self, label='Statistics Options')
    self.multi_box_sizer = wx.StaticBoxSizer(multi_box, wx.VERTICAL)
    self.multi_opt_sizer = wx.GridBagSizer(2, 2)

    self.multi_opt_sizer.Add(self.opt_multi, pos=(0, 0),
                             flag=wx.LEFT | wx.TOP | wx.RIGHT, border=10)
    self.multi_opt_sizer.Add(self.trial_number, pos=(1, 0),
                             flag=wx.ALL, border=10)
    self.multi_opt_sizer.Add(self.tag_list, pos=(0, 1), span=(2, 1),
                             flag=wx.BOTTOM | wx.TOP | wx.RIGHT, border=10)
    self.multi_box_sizer.Add(self.multi_opt_sizer, flag=wx.EXPAND)
    self.bottom_sizer.Add(self.multi_box_sizer)

    self.opt_cluster = gctr.VerticalOptionCtrl(self,
                                               ctrl_size=(100, -1),
                                               sub_labels=['No. of images',
                                                           'Threshold'],
                                               items=[('num_images', 1000),
                                                      ('threshold', 250)])

    self.cluster_btn_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.cluster_btn_sizer.AddGrowableCol(0)
    self.cluster_btn_sizer.AddGrowableCol(1)
    self.btn_cluster = wx.Button(self, label='Clustering')
    self.btn_histogram = wx.Button(self, label='Histogram')
    self.cluster_btn_sizer.Add(self.btn_cluster, flag=wx.EXPAND)
    self.cluster_btn_sizer.Add(self.btn_histogram, flag=wx.EXPAND)

    cluster_box = wx.StaticBox(self, label='Unit Cell Clustering')
    self.cluster_box_sizer = wx.StaticBoxSizer(cluster_box, wx.VERTICAL)
    self.cluster_box_sizer.Add(self.opt_cluster, flag=wx.ALL, border=10)
    self.cluster_box_sizer.Add(self.cluster_btn_sizer,
                               flag=wx.ALL | wx.EXPAND, border=10)
    self.bottom_sizer.Add(self.cluster_box_sizer)

    self.main_sizer.Add(self.status_panel, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.iso_panel, 1,
                        flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.bottom_sizer,
                        flag=wx.EXPAND | wx.ALL, border=10)


    # Bindings
    self.Bind(wx.EVT_CHOICE, self.onTrialChoice, self.trial_number.ctr)
    self.Bind(wx.EVT_CHECKLISTBOX, self.onTagCheck, self.tag_list.ctr)
    self.Bind(wx.EVT_TEXT_ENTER, self.onMultiplicityGoal, self.opt_multi.goal)
    self.Bind(wx.EVT_BUTTON, self.onClustering, self.btn_cluster)
    self.Bind(wx.EVT_BUTTON, self.onHistogram, self.btn_histogram)
    self.Bind(EVT_PRG_REFRESH, self.onRefresh)
    self.Bind(EVT_CLUSTERING, self.onClusteringResult)


  def onHistogram(self, e):
    trial = self.main.db.get_trial(trial_number=self.trial_no)
    runblocks = trial.rungroups

    clustering = ClusteringWorker(self, trial=trial, runblocks=runblocks,
                                  output=self.main.params.output_folder,
                                  threshold=self.opt_cluster.threshold.GetValue(),
                                  sample_size=self.opt_cluster.num_images.GetValue())
    clustering.histogram()

  def onClustering(self, e):
    trial = self.main.db.get_trial(trial_number=self.trial_no)
    runblocks = trial.rungroups

    clustering = ClusteringWorker(self, trial=trial, runblocks=runblocks,
                                  output=self.main.params.output_folder,
                                  threshold=self.opt_cluster.threshold.GetValue(),
                                  sample_size=self.opt_cluster.num_images.GetValue())
    clustering.run()

  def onClusteringResult(self, e):
    from string import ascii_uppercase
    self.iso_box_sizer.DeleteWindows()

    if e.GetValue() is None:                       # nothing returned
      print 'Nothing to cluster!'
    elif type(e.GetValue()) in (list, tuple):      # uc param list returned
      plotter = pltr.PopUpCharts()
      plotter.plot_uc_histogram(info=e.GetValue())
    else:                                          # uc clusters returned
      counter = 0
      clusters = sorted(e.GetValue(), key=lambda x: x.members, reverse=True)
      for cluster in clusters:
        sorted_pg_comp = sorted(cluster.pg_composition.items(),
                                key=lambda x: -1 * x[1])
        pg_nums = [pg[1] for pg in sorted_pg_comp]
        cons_pg = sorted_pg_comp[np.argmax(pg_nums)]

        if len(cluster.members) > 1:

          # format and record output
          uc_line = "{:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}" \
                    "".format(cluster.medians[0], cluster.medians[1],
                              cluster.medians[2], cluster.medians[3],
                              cluster.medians[4], cluster.medians[5])

          iso = gctr.IsoformInfoCtrl(self.iso_panel)
          iso.ctr_iso.SetValue(ascii_uppercase[counter])
          iso.ctr_pg.SetValue(cons_pg[0])
          iso.ctr_num.SetValue(str(len(cluster.members)))
          iso.ctr_uc.SetValue(uc_line)
          iso.uc_values = [i.uc for i in cluster.members]
          self.iso_box_sizer.Add(iso,
                                 flag=wx.EXPAND| wx.TOP | wx.LEFT | wx.RIGHT,
                                 border=10)

          counter += 1

    self.iso_panel.SetSizer(self.iso_box_sizer)
    self.iso_panel.Layout()
    self.iso_panel.SetupScrolling()

  def onTagCheck(self, e):
    checked_items = self.tag_list.ctr.GetCheckedStrings()
    self.selected_tags = [i for i in self.main.db.get_all_tags() if i.name
                          in checked_items]
    self.main.run_window.prg_light.change_status('idle')
    self.tag_trial_changed = True
    self.rows = {}

  def onTrialChoice(self, e):
    self.trial_no = self.trial_number.ctr.GetSelection()
    self.main.run_window.prg_light.change_status('idle')
    self.tag_trial_changed = True
    self.rows = {}

  def find_tags(self):
    self.tag_list.ctr.Clear()
    self.all_tags = [str(i.name) for i in self.main.db.get_all_tags()]
    self.tag_list.ctr.InsertItems(items=self.all_tags, pos=0)

  def find_trials(self):
    self.all_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    self.trial_number.ctr.Clear()
    for trial in self.all_trials:
      self.trial_number.ctr.Append(trial)
    self.trial_number.ctr.SetSelection(self.trial_no)

  def onMultiplicityGoal(self, e):
    goal = self.opt_multi.goal.GetValue()
    if goal.isdigit() and goal != '':
      self.multiplicity_goal = int(goal)

  def onRefresh(self, e):
    # Find new tags
    all_db_tags = [str(i.name) for i in self.main.db.get_all_tags()]
    new_tags = [i for i in all_db_tags if i not in self.all_tags]
    if len(new_tags) > 0:
      self.find_tags()

    # Find new trials
    all_db_trials = [str(i.trial) for i in self.main.db.get_all_trials()]
    new_trials = [i for i in all_db_trials if i not in self.all_trials]
    if len(new_trials) > 0:
      self.find_trials()

    # Show info
    self.refresh_rows(info=e.GetValue())
    self.status_panel.SetSizer(self.status_sizer)
    #self.status_panel.Layout()
    self.status_sizer.Layout()
    self.status_panel.SetupScrolling()
    self.status_box.SetLabel('Data Statistics - Trial {}'.format(self.trial_no))

  def refresh_rows(self, info={}):
    ''' Refresh status data '''

    if info != {}:
      max_multiplicity = max([info[i]['multiplicity'] for i in info])
      if max_multiplicity > self.multiplicity_goal:
        xmax = max_multiplicity
      else:
        xmax = self.multiplicity_goal

      if self.redraw_windows:
        self.status_sizer.Clear(deleteWindows=True)
        self.rows = {}
        row = None

      for iso, values in info.iteritems():
        if not self.redraw_windows:
          row = self.rows[values['isoform']]['row']

        self.update_row(row=row,
                        name=values['isoform'],
                        value=values['multiplicity'],
                        bins=values['bins'],
                        xmax=xmax)
      self.redraw_windows = False

    else:
      self.status_sizer.Clear(deleteWindows=True)

  def update_row(self, row, name, value, bins, xmax):
    ''' Add new row, or update existing '''

    bin_choices = [("Bin {}:  {:3.2f} - {:3.2f}" \
                    "".format(b.number, float(b.d_max), float(b.d_min)),
                    bins.index(b)) for b in bins]
    if row is None:
      self.rows[name] = {}
      row = pltr.SingleBarPlot(self.status_panel,
                               label='Isoform {}'.format(name),
                               gauge_size=(250, 15),
                               choice_label='High res. limit:',
                               choice_size=(160, -1),
                               choices = bin_choices)
      self.status_sizer.Add(row, flag=wx.EXPAND | wx.ALL, border=10)
      row.bins.ctr.SetSelection(row.bins.ctr.GetCount() - 1)

    row.redraw_axes(value=value, goal=self.multiplicity_goal, xmax=xmax)
    self.rows[name]['row'] = row
    self.rows[name]['high_bin'] = row.bins.ctr.GetClientData(
      row.bins.ctr.GetSelection())

class MergeTab(BaseTab):
  def __init__(self, parent, prefix='prime'):
    BaseTab.__init__(self, parent=parent)

    #TODO: This seems like a placeholder. Need to have a proper window that
    #TODO: accepts multiple input entries, whether they are files or folders

    #TODO: alternatively, concatenate any combo of tags into an input file

    self.prefix = prefix
    self.prime_filename = '{}.phil'.format(self.prefix)

    self.prime_panel = PRIMEInputWindow(self)
    self.btn_get_tags = wx.Button(self, label='Select Tags...', size=(120, -1))
    self.btn_load_phil = wx.Button(self, label='Load PHIL', size=(120, -1))
    self.btn_save_phil = wx.Button(self, label='Save PHIL', size=(120, -1))
    self.btn_run_prime = wx.Button(self, label='Run PRIME', size=(120, -1))
    self.btn_run_prime.Disable()
    self.btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.btn_sizer.Add(self.btn_get_tags)
    self.btn_sizer.Add(self.btn_load_phil, flag=wx.LEFT, border=5)
    self.btn_sizer.Add(self.btn_save_phil, flag=wx.LEFT, border=5)
    self.btn_sizer.Add(self.btn_run_prime, flag=wx.LEFT, border=5)

    self.main_sizer.Add(self.prime_panel, flag=wx.ALL | wx.EXPAND, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_sizer,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM,
                        border=10)

    self.Bind(wx.EVT_BUTTON, self.onInput, self.prime_panel.inp_box.btn_browse)
    self.Bind(wx.EVT_TEXT, self.onInput, self.prime_panel.inp_box.ctr)
    self.Bind(wx.EVT_BUTTON, self.onIsoRef, self.prime_panel.ref_box.btn_browse)
    self.Bind(wx.EVT_TEXT, self.onIsoRef, self.prime_panel.ref_box.ctr)
    self.Bind(wx.EVT_BUTTON, self.onRun, self.btn_run_prime)
    self.Bind(wx.EVT_BUTTON, self.onLoad, self.btn_load_phil)
    self.Bind(wx.EVT_BUTTON, self.onSave, self.btn_save_phil)

  def onInput(self, e):
    if self.prime_panel.inp_box.ctr.GetValue() != '':
      self.btn_run_prime.Enable()
    else:
      self.btn_run_prime.Disable()

  def onLoad(self, e):
    # Extract params from file
    load_dlg = wx.FileDialog(self,
                             message="Load script file",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.OPEN | wx.FD_FILE_MUST_EXIST,
                             )
    if load_dlg.ShowModal() == wx.ID_OK:
      script = load_dlg.GetPaths()[0]
      out_dir = os.path.dirname(script)
      self.prime_filename = os.path.basename(script)
      self.load_script(out_dir=out_dir)
    load_dlg.Destroy()

  def load_script(self, out_dir):
    ''' Loads PRIME script '''
    import iotbx.phil as ip

    script = os.path.join(out_dir, self.prime_filename)
    user_phil = ip.parse(open(script).read())
    self.pparams = master_phil.fetch(sources=[user_phil]).extract()
    self.prime_panel.pparams = self.pparams

    self.prime_panel.inp_box.ctr.SetValue(str(self.pparams.data[0]))
    current_dir = os.path.dirname(self.pparams.run_no)
    self.prime_panel.out_box.ctr.SetValue(str(current_dir))
    if str(self.prime_panel.out_box.ctr.GetValue).lower() == '':
      self.prime_panel.out_box.ctr.SetValue(self.out_dir)
    if str(self.pparams.title).lower() != 'none':
      self.prime_panel.title_box.ctr.SetValue(str(self.pparams.title))
    if str(self.pparams.hklisoin).lower() != 'none':
      self.prime_panel.ref_box.ctr.SetValue(str(self.pparams.hklisoin))
    elif str(self.pparams.hklrefin).lower() != 'none':
      self.prime_panel.ref_box.ctr.SetValue(str(self.pparams.hklrefin))
      self.prime_panel.opt_chk_useref.SetValue(True)
    if str(self.pparams.n_residues).lower() == 'none':
      self.prime_panel.opt_spc_nres.SetValue(500)
    else:
      self.prime_panel.opt_spc_nres.SetValue(int(self.pparams.n_residues))
    self.prime_panel.opt_spc_nproc.SetValue(int(self.pparams.n_processors))

  def onSave(self, e):
    self.init_settings()

    # Generate text of params
    final_phil = master_phil.format(python_object=self.pparams)
    with misc.Capturing() as txt_output:
      final_phil.show()
    txt_out = ''
    for one_output in txt_output:
      txt_out += one_output + '\n'

    # Save param file
    save_dlg = wx.FileDialog(self,
                             message="Save PRIME Script",
                             defaultDir=os.curdir,
                             defaultFile="*.phil",
                             wildcard="*.phil",
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT
                             )
    if save_dlg.ShowModal() == wx.ID_OK:
      with open(save_dlg.GetPath(), 'w') as savefile:
        savefile.write(txt_out)

  def onIsoRef(self, e):
    if self.prime_panel.ref_box.ctr.GetValue() != '':
      self.prime_panel.opt_chk_useref.Enable()
    else:
      self.prime_panel.opt_chk_useref.Disable()

  def init_settings(self):
    self.pparams = self.prime_panel.pparams
    self.pparams.data = [self.prime_panel.inp_box.ctr.GetValue()]
    self.out_dir = self.prime_panel.out_box.ctr.GetValue()
    self.pparams.run_no = misc.set_base_dir(out_dir=self.out_dir)
    self.pparams.title = self.prime_panel.title_box.ctr.GetValue()
    if str(self.prime_panel.ref_box.ctr.GetValue()).lower() != '':
      self.pparams.hklisoin = self.prime_panel.ref_box.ctr.GetValue()
      if self.prime_panel.opt_chk_useref.GetValue():
        self.pparams.hklrefin = self.prime_panel.ref_box.ctr.GetValue()
    self.pparams.n_residues = self.prime_panel.opt_spc_nres.GetValue()
    self.pparams.n_processors = self.prime_panel.opt_spc_nproc.GetValue()

  def onRun(self, e):
    # Run full processing

    self.init_settings()
    prime_phil = master_phil.format(python_object=self.pparams)

    with misc.Capturing() as output:
      prime_phil.show()

    txt_out = ''
    for one_output in output:
      txt_out += one_output + '\n'

    prime_file = os.path.join(self.out_dir, self.prime_filename)
    out_file = os.path.join(self.out_dir, 'stdout.log')
    with open(prime_file, 'w') as pf:
      pf.write(txt_out)

    self.prime_run_window = PRIMERunWindow(self, -1,
                                           title='PRIME Output',
                                           params=self.pparams,
                                           prime_file=prime_file,
                                           out_file=out_file)
    self.prime_run_window.prev_pids = easy_run.fully_buffered('pgrep -u {} {}'
                                          ''.format(user, 'python')).stdout_lines
    self.prime_run_window.Show(True)


# ------------------------------- UI Elements -------------------------------- #

class TrialPanel(wx.Panel):
  ''' A scrolled panel that contains run blocks and trial controls '''

  def __init__(self, parent, db, trial, box_label=None):
    wx.Panel.__init__(self, parent=parent, size=(250, 300))

    self.db = db
    self.trial = trial
    self.parent = parent

    trial_box = wx.StaticBox(self, label=box_label)
    self.main_sizer = wx.StaticBoxSizer(trial_box, wx.VERTICAL)

    self.block_panel = ScrolledPanel(self, size=(150, 180))
    self.block_sizer = wx.BoxSizer(wx.VERTICAL)
    self.block_panel.SetSizer(self.block_sizer)
    self.one_block_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.add_panel = wx.Panel(self)
    self.add_sizer = wx.BoxSizer(wx.VERTICAL)
    self.add_panel.SetSizer(self.add_sizer)

    # Add "New Block" button to a separate sizer (so it is always on bottom)
    self.btn_add_block = wx.Button(self.add_panel, label='New Block',
                                   size=(200, -1))
    self.btn_view_phil = wx.BitmapButton(self.add_panel,
                        bitmap=wx.Bitmap('{}/16x16/viewmag.png'.format(icons)))
    self.chk_active = wx.CheckBox(self.add_panel, label='Active Trial')
    self.view_sizer = wx.FlexGridSizer(1, 2, 0, 10)
    self.view_sizer.Add(self.btn_view_phil)
    self.view_sizer.Add(self.chk_active, flag=wx.EXPAND)

    self.add_sizer.Add(self.btn_add_block,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                       border=10)
    self.add_sizer.Add(self.view_sizer,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT,
                       border=10)

    self.main_sizer.Add(self.block_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.add_panel, flag=wx.ALL | wx.ALIGN_BOTTOM, border=5)

    # Bindings
    self.Bind(wx.EVT_BUTTON, self.onAddBlock, self.btn_add_block)
    self.Bind(wx.EVT_BUTTON, self.onViewPHIL, self.btn_view_phil)
    self.chk_active.Bind(wx.EVT_CHECKBOX, self.onToggleActivity)

    self.SetSizer(self.main_sizer)

  def onViewPHIL(self, e):
    view_dlg = dlg.TrialDialog(self, db=self.db, trial=self.trial, new=False)
    view_dlg.ShowModal()
    view_dlg.Destroy()

  def onToggleActivity(self, e):
    if self.chk_active.GetValue():
      self.trial.active = True
    else:
      self.trial.active = False

  def onAddBlock(self, e):
    rblock_dlg = dlg.RunBlockDialog(self, trial=self.trial,
                                    db=self.db)
    rblock_dlg.Fit()

    if (rblock_dlg.ShowModal() == wx.ID_OK):
      self.refresh_trial()
    rblock_dlg.Destroy()

  def refresh_trial(self):
    self.block_sizer.DeleteWindows()
    self.active_blocks = self.trial.rungroups
    for block in self.active_blocks:
      self.draw_block_button(block)
    self.block_panel.Layout()
    self.block_panel.SetupScrolling()

  def draw_block_button(self, block):
    ''' Add new run block button '''
    new_block = gctr.RunBlock(self.block_panel, block=block)
    self.Bind(wx.EVT_BUTTON, self.onRunBlockOptions, new_block.new_runblock)
    self.block_sizer.Add(new_block,
                         flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                         border=5)

  def onRunBlockOptions(self, e):
    ''' Open dialog and change run_block options '''
    run_block = e.GetEventObject().block
    rblock_dlg = dlg.RunBlockDialog(self, block=run_block,
                                    db=self.db)
    rblock_dlg.Fit()

    if (rblock_dlg.ShowModal() == wx.ID_OK):
      wx.CallAfter(self.refresh_trial)
    rblock_dlg.Destroy()

class RunEntry(wx.Panel):
  ''' Adds run row to table, with average and view buttons'''
  def __init__(self, parent, run, params):
    self.run = run
    self.params = params

    wx.Panel.__init__(self, parent=parent)

    self.sizer = wx.FlexGridSizer(1, 4, 0, 10)
    run_no = wx.StaticText(self, label=str(run.run),
                           size=(60, -1))
    self.tag_button = gctr.TagButton(self, run=run)
    self.avg_button = wx.Button(self, label='Average')
    self.view_button = wx.Button(self, label='View')
    self.view_button.Hide()

    self.sizer.Add(run_no, flag=wx.EXPAND | wx.ALIGN_CENTRE)
    self.sizer.Add(self.tag_button, flag=wx.EXPAND)
    self.sizer.AddGrowableCol(1)
    self.sizer.Add(self.avg_button)
    self.sizer.Add(self.view_button, flag=wx.RESERVE_SPACE_EVEN_IF_HIDDEN)

    # Button Bindings
    self.Bind(wx.EVT_BUTTON, self.onTagButton, self.tag_button)
    self.Bind(wx.EVT_BUTTON, self.onAvgButton, self.avg_button)
    self.Bind(wx.EVT_BUTTON, self.onViewButton, self.view_button)

    self.SetSizer(self.sizer)

  def onTagButton(self, e):
    self.tag_button.change_tags()

  def onAvgButton(self, e):
    avg = dlg.AveragingDialog(self, self.run, self.params)
    avg.Fit()
    avg.Center()

    if (avg.ShowModal() == wx.ID_OK):
      e.GetEventObject().SetLabel('Running')
      e.GetEventObject().Disable()
      self.view_button.Show()
      # TODO: hook up the calibration app

  def onViewButton(self, e):
    # TODO: hook up view function
    pass
