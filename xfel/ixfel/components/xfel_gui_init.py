from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/02/2016
Last Changed: 06/02/2016
Description : XFEL UI Initialization module
'''

import os
import wx
import random
from wx.lib.scrolledpanel import ScrolledPanel

import xfel.ui.components.xfel_gui_controls as gctr
import xfel.ui.components.xfel_gui_dialogs as dlg

icons = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'icons/')

license = 'cctbx.xfel and cctbx.xfel UI are developed under the open source ' \
          'license'

description = 'The cctbx.xfel UI is developed for use during data collection ' \
              'and initial processing at LCSL XFEL beamlines.'

# ------------------------------- Main Window -------------------------------- #

class MainWindow(wx.Frame):

  def __init__(self, parent, id, title):
    wx.Frame.__init__(self, parent, id, title, size=(800, 500))

    # Toolbar
    self.toolbar = self.CreateToolBar(wx.TB_TEXT)
    self.tb_btn_quit = self.toolbar.AddLabelTool(wx.ID_EXIT, label='Quit',
                                                 bitmap=wx.Bitmap('{}/32x32/exit.png'.format(icons)),
                                                 shortHelp='Quit',
                                                 longHelp='Exit iXFEL')

    self.toolbar.AddSeparator()
    self.tb_btn_run = self.toolbar.AddLabelTool(wx.ID_ANY, label='Run',
                                                bitmap=wx.Bitmap('{}/32x32/run.png'.format(icons)),
                                                shortHelp='Run',
                                                longHelp='Activate all sentinels')
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

    # Toolbar button bindings
    self.Bind(wx.EVT_TOOL, self.onQuit, self.tb_btn_quit)
    self.Bind(wx.EVT_TOOL, self.onRun, self.tb_btn_run)

    # Draw the main window sizer
    self.SetSizer(main_box)

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
    self.Destroy()


class RunWindow(wx.Panel):
  ''' Window panel that will house all the run tabs '''

  def __init__(self, parent):
    self.parent = parent
    super(RunWindow, self).__init__(self.parent)

    self.main_panel = wx.Panel(self)
    self.main_nbook = wx.Notebook(self.main_panel, style=0)
    self.runs_tab = RunTab(self.main_nbook)
    self.trials_tab = TrialsTab(self.main_nbook)
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
  def __init__(self, parent):
    BaseTab.__init__(self, parent=parent)
    self.last_run = 0
    self.all_tag_buttons = []

    # TODO: Need to figure out how to set up initial tags
    self.all_tags = [gctr.Sample(os.getlogin(), 'Login name'),
                     gctr.Sample('test', 'Tags for testing'),
                     gctr.Sample('initial sample', 'WT luciferase'),
                     gctr.Sample('first exposure', '40 ns at 400 nm')]

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.run_panel = ScrolledPanel(self)
    self.run_sizer = wx.BoxSizer(wx.VERTICAL)
    self.run_panel.SetSizer(self.run_sizer)

    self.colname_sizer = wx.FlexGridSizer(1, 4, 0, 10)
    run_label = wx.StaticText(self, label='Run', size=(60, -1))
    tag_label = wx.StaticText(self, label='Sample Tags', size=(300, -1))
    img_label = wx.StaticText(self, label='# of Images', size=(160, -1))
    dur_label = wx.StaticText(self, label='Duration', size=(160, -1))
    self.colname_sizer.Add(run_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(tag_label, flag=wx.ALIGN_RIGHT | wx.EXPAND)
    self.colname_sizer.Add(img_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.Add(dur_label, flag=wx.ALIGN_RIGHT)
    self.colname_sizer.AddGrowableCol(1, 1)
    self.main_sizer.Add(self.colname_sizer,
                        flag=wx.ALL | wx.EXPAND, border=10)

    self.btn_manage_tags = wx.Button(self, label='Manage Tags', size=(120, -1))
    self.main_sizer.Add(self.run_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_manage_tags,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM | wx.ALIGN_RIGHT,
                        border=10)

    # Timer to simulate run acquisition
    # TODO: replace timer with actual run updates
    self.last_run += 1
    self.add_row()
    self.run_panel.SetupScrolling()
    self.run_panel.Refresh()
    self.timer = wx.Timer(self)
    self.timer.Start(25000)

    self.SetSizer(self.main_sizer)

    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())
    self.Bind(wx.EVT_BUTTON, self.onManageTags, self.btn_manage_tags)

  def onTimer(self, e):
    self.last_run += 1
    self.add_row()
    self.run_panel.SetupScrolling()
    self.run_panel.Refresh()

  def onManageTags(self, e):
    ''' User can add / remove / edit sample tags '''
    mtag_dlg = dlg.TagDialog(self,
                             tags=self.all_tags)
    mtag_dlg.Fit()

    # Get new tag objects
    if (mtag_dlg.ShowModal() == wx.ID_OK):
      self.all_tags = mtag_dlg.tags
    mtag_dlg.Destroy()

    # Update tags on all tag buttons
    for btn in self.all_tag_buttons:
      tags_in_button = [i.tag for i in btn.tags]
      btn.tags = [j for j in self.all_tags if j.tag in tags_in_button]
      btn.update_label()

  def add_row(self):
    ''' Adds run row to table, matching colname_sizer '''
    row_sizer = wx.FlexGridSizer(1, 4, 0, 10)
    run_no = wx.StaticText(self.run_panel, label=str(self.last_run),
                           size=(60, -1))
    btn_tags = gctr.TagButton(self.run_panel)
    btn_tags.tags = [self.all_tags[0]]
    btn_tags.update_label()
    self.Bind(wx.EVT_BUTTON, self.onTrialButton, id=btn_tags.GetId())
    self.all_tag_buttons.append(btn_tags)
    images = wx.StaticText(self.run_panel,
                           label=str(random.randrange(100, 20000)),
                           size=(160, -1),
                           style=wx.ALIGN_CENTER)
    duration = wx.StaticText(self.run_panel,
                             label='{:02d} min, {:02d} sec'
                                   ''.format(random.randrange(0, 20),
                                             random.randrange(0, 59)),
                             size=(160, -1),
                             style=wx.ALIGN_CENTER)
    row_sizer.Add(run_no, flag=wx.EXPAND | wx.ALIGN_CENTRE)
    row_sizer.Add(btn_tags, flag=wx.EXPAND)
    row_sizer.Add(images, flag=wx.EXPAND)
    row_sizer.Add(duration, flag=wx.EXPAND)
    row_sizer.AddGrowableCol(1)
    self.run_sizer.Add(row_sizer, flag=wx.ALL | wx.EXPAND, border=0)

  def onTrialButton(self, e):
    e.GetEventObject().change_tags(self.all_tags)


class TrialsTab(BaseTab):
  def __init__(self, parent):
    BaseTab.__init__(self, parent=parent)

    self.trial_number = 1

    self.main_sizer = wx.BoxSizer(wx.VERTICAL)

    self.trial_panel = ScrolledPanel(self, size=(300, 350))
    self.trial_sizer = wx.BoxSizer(wx.HORIZONTAL)

    self.add_new_trial()

    self.btn_add_trial = wx.Button(self, label='New Trial',
                                   size=(120, -1))

    self.main_sizer.Add(self.trial_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(wx.StaticLine(self), flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.btn_add_trial,
                        flag=wx.RIGHT | wx.LEFT | wx.BOTTOM | wx.ALIGN_RIGHT,
                        border=10)

    self.trial_panel.SetSizer(self.trial_sizer)
    #self.new_trial_panel.SetSizer(self.new_trial_sizer)
    self.SetSizer(self.main_sizer)

    self.Bind(wx.EVT_BUTTON, self.onAddTrial, self.btn_add_trial)

  def onAddTrial(self, e):
    all_children = self.trial_sizer.GetChildren()
    trials = [i for i in all_children if isinstance(i.GetWindow(),
                                                    TrialPanel)]
    self.trial_number = len(trials) + 1
    self.add_new_trial()
    self.trial_panel.Layout()
    self.trial_panel.SetupScrolling()

  def add_new_trial(self):
    new_trial = TrialPanel(self.trial_panel,
                           box_label='Trial {}'.format(self.trial_number))
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
  def __init__(self, parent, box_label=None):
    wx.Panel.__init__(self, parent=parent, size=(200, 300))

    # Run variables
    self.first_run = 1
    self.last_run = 0

    trial_box = wx.StaticBox(self, label=box_label)
    self.main_sizer = wx.StaticBoxSizer(trial_box, wx.VERTICAL)

    self.block_panel = ScrolledPanel(self, size=(150, 200))
    self.block_sizer = wx.BoxSizer(wx.VERTICAL)
    self.block_panel.SetSizer(self.block_sizer)
    self.add_panel = wx.Panel(self)
    self.add_sizer = wx.BoxSizer(wx.VERTICAL)
    self.add_panel.SetSizer(self.add_sizer)

    self.chk_active = wx.CheckBox(self, label='Active')
    self.main_sizer.Add(self.chk_active, flag=wx.ALIGN_CENTER)

    self.add_new_block()

    # Add "New Block" button to a separate sizer (so it is always on bottom)
    self.btn_add_block = wx.Button(self.add_panel, label='New Block',
                                   size=(120, -1))
    self.add_sizer.Add(self.btn_add_block,
                       flag=wx.TOP | wx.LEFT | wx.RIGHT | wx.ALIGN_CENTER,
                    border = 10)

    self.main_sizer.Add(self.block_panel, 1, flag=wx.EXPAND | wx.ALL, border=10)
    self.main_sizer.Add(self.add_panel, flag=wx.ALL | wx.ALIGN_BOTTOM, border=5)

    # Create and start timer
    # TODO: Replace timer with actual run updates
    self.timer = wx.Timer(self)
    self.timer.Start(5000)

    # Button bindings
    self.Bind(wx.EVT_BUTTON, self.onAddBlock, self.btn_add_block)
    self.Bind(wx.EVT_TIMER, self.onTimer, id=self.timer.GetId())

    self.SetSizer(self.main_sizer)

  def onAddBlock(self, e):
    ''' Add new block button '''

    # Get list of all buttons in block sizer
    all_children = self.block_sizer.GetChildren()
    run_buttons = [i for i in all_children if isinstance(i.GetWindow(),
                                                         gctr.RunBlockButton)]

    # Get last button in block sizer, update label with latest run
    last_block = run_buttons[-1].GetWindow()
    last_block.last_run = self.last_run
    last_block.update_label()
    self.first_run = self.last_run + 1

    self.add_new_block()
    self.block_panel.Layout()
    self.block_panel.SetupScrolling()

  def add_new_block(self):
    ''' Add new run block button '''
    new_block = gctr.RunBlockButton(self.block_panel, size=(120, -1))
    new_block.first_run = self.first_run
    new_block.update_label()
    self.block_sizer.Add(new_block, flag=wx.TOP | wx.LEFT | wx.RIGHT |
                                         wx.ALIGN_CENTER,
                         border=5)

  def onTimer(self, e):
    self.last_run += 1
