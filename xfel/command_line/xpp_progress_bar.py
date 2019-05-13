# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME xpp.progress_bar
#
'''
This is a wxpython GUI application for reading a database and
displaying a plot showing multiplicity by trial.
Currently it features:
* contains the matplotlib navigation toolbar
* allows user input of desired experiment trial number
* finds the multiplicity of data currently in the trial
* Displays this information graphically
* automatically updates the plot on a timer
* The plot can be saved to a file from the file menu

'''
from __future__ import division
from __future__ import print_function
from six.moves import range
import os
import wx
import numpy as np
import sys
import time
import threading
# The recommended way to use wx with mpl is with the WXAgg
# backend.
#
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from xfel.command_line.xpp_progress_detail import master_phil
import iotbx.phil

myEVT_DB_STATS = wx.NewEventType()
EVT_DB_STATS = wx.PyEventBinder(myEVT_DB_STATS, 1)
class StatsEvent(wx.PyCommandEvent):
  """Event to signal that a count value is ready"""
  def __init__(self, etype, eid, value=None):
      """Creates the event object"""
      wx.PyCommandEvent.__init__(self, etype, eid)
      self._value = value

  def GetValue(self):
      """Returns the value from the event.
      @return: the value of this event
      """
      return self._value

class DB_thread(threading.Thread):
    def __init__(self, parent):
       """
         @param parent: The gui object that should recieve the value
         @param value: value to 'calculate' to

       """
       threading.Thread.__init__(self)
       self.parent = parent

    def run(self):
       while True:
         from xfel.xpp.progress_utils import application
         stats = application(self.parent.params, loop = False)
         evt = StatsEvent(myEVT_DB_STATS, -1, stats)
         wx.PostEvent(self.parent, evt)
         import time
         time.sleep(5)

class BarsFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Multiplicity Progress'

    def __init__(self, params):
        wx.Frame.__init__(self, None, -1, self.title)
        self.params = params
        self.stats = None

        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()


        #set up the thread
        self.Bind(EVT_DB_STATS, self.on_stats_update)
        self.worker = DB_thread(self)
        self.worker.start()

    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

        menu_help = wx.Menu()
        m_help = menu_help.Append(-1, "&About\tF1", "About the data")
        self.Bind(wx.EVT_MENU, self.on_about, m_help)

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.
        # 10.5x5.5 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((10.5, 5.5), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)

        # set up trial text box
        self.textbox = wx.TextCtrl(self.panel, -1, '')
        if self.params.trial is not None:
          self.textbox.SetValue(str(self.params.trial))
        self.fnameLabel = wx.StaticText(self.panel, wx.ID_ANY, 'Trial Number ')
        self.fnameLabel.SetFont(wx.Font (12, wx.SWISS, wx.NORMAL, wx.BOLD))

        # set up resolution text box
        self.restextbox = wx.TextCtrl(self.panel, -1, '')
        if self.params.resolution is not None:
          self.restextbox.SetValue(str(self.params.resolution))
        self.resLabel = wx.StaticText(self.panel, wx.ID_ANY, 'Resolution Shell')
        self.resLabel.SetFont(wx.Font (12, wx.SWISS, wx.NORMAL, wx.BOLD))

        # set up tags text box
        # this determines which crystal isoforms in which illumination stats are plotted
        self.tagstextbox = wx.TextCtrl(self.panel, -1, '')
        if self.params.run_tags is not None:
          self.tagstextbox.SetValue(self.params.run_tags)
        self.tagsLabel = wx.StaticText(self.panel, wx.ID_ANY, 'Tags')
        self.tagsLabel.SetFont(wx.Font (12, wx.SWISS, wx.NORMAL, wx.BOLD))

        # set up the submit button
        self.SubmitButton = wx.Button(self.panel, wx.ID_ANY, 'Submit')
        self.Bind(wx.EVT_BUTTON, self.on_submit, self.SubmitButton)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        #

        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)

        #
        # Layout with box sizers
        #
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW )
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.DataSizer = wx.BoxSizer(wx.HORIZONTAL)

        # create a box for entering the trial number
        self.DataSizer.Add(self.fnameLabel, 0 , wx.ALL | wx.ALIGN_CENTER, 5)
        self.DataSizer.Add(self.textbox, 1, wx.ALL | wx.ALIGN_CENTER , 5)

        # create a box for entering the resolution shell
        self.DataSizer.Add(self.resLabel, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        self.DataSizer.Add(self.restextbox, 1, wx.ALIGN_CENTER| wx.ALL, 5)

        self.DataSizer.Add(self.tagsLabel, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        self.DataSizer.Add(self.tagstextbox, 1, wx.ALIGN_CENTER| wx.ALL, 5)

        self.DataSizer.Add(self.SubmitButton, 0, wx.ALL | wx.EXPAND , 5)

        self.vbox.Add(self.DataSizer, 0 , wx.ALL |wx.CENTER, 5)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_figure(self):
        """ Redraws the figure
        """
        stats = self.stats
        if stats is None: return
        if len(stats) == 0:
          return

        print(stats)
        res = self.restextbox.GetValue()
        trial = self.textbox.GetValue()
        self.mult_highest = [stats[key]['multiplicity_highest'] for key in stats.keys()]
        self.mult = [stats[key]['multiplicity'] for key in stats.keys()]
        plot_max = max(max(self.mult), max(self.mult_highest)) + 1
        pos = np.arange(len(self.mult))+0.5 # the bar centers on the y-axis
        labels = [stats.keys()[i] for i in range(len(stats.keys()))]
        n = len(labels)
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()
        self.axes.barh(
            pos,
            self.mult,
            align='center',
            color='lightskyblue')
        self.axes.barh(
            pos,
            self.mult_highest,
            align='center',
            color='blue')
        title = 'Trial: %(tr)s Overall multplicity (light blue) and multiplicity at %(ang)s angstroms (blue).\n' \
          % {"tr": trial, "ang": res} \
          + 'Updated %(HH)02d:%(MM)02d:%(SS)02d' \
          % {"HH":time.localtime().tm_hour, "MM":time.localtime().tm_min, "SS":time.localtime().tm_sec}
        self.axes.set_title(title)
        self.axes.set_yticks(pos)
        self.axes.set_yticklabels(labels)
        self.axes.set_xlabel('Multiplicity')
        self.axes.set_xlim(0.0,plot_max)
        self.canvas.draw()

    def on_submit(self, event):
        self.params.trial = int(self.textbox.GetValue())
        self.params.resolution = float(self.restextbox.GetValue())
        self.params.run_tags = self.tagstextbox.GetValue()

    def on_stats_update(self,event):
        self.stats = event.GetValue()
        self.draw_figure()

    def on_save_plot(self, event):
        file_choices = "PNG (*.png)|*.png"

        dlg = wx.FileDialog(
            self,
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="plot.png",
            wildcard=file_choices,
            style=wx.SAVE)

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path, dpi=self.dpi)
            self.flash_status_message("Saved to %s" % path)

    def on_exit(self, event):
        self.redraw_timer.Stop()
        self.Destroy()

    def on_about(self, event):
        msg = """ A progress plot using wxPython with matplotlib:

         * Use the text box to enter experiment trial number
         * Save the plot to a file using the File menu

        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()

    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER,
            self.on_flash_status_off,
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')

if __name__ == '__main__':
  args = sys.argv[1:]
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.xpp.progress_utils import phil_validation
  phil_validation(work_params)
  if ("--help" in args) :
    libtbx.phil.parse(master_phil.show())
  else:
    app = wx.PySimpleApp()
    app.frame = BarsFrame(work_params)
    app.frame.Show()
    app.MainLoop()
