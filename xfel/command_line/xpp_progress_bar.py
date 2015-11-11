# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME xpp.progress_bar
#
'''
This is a wxpython GUI application for reading a database and
displaying a plot showing completeness and multiplicity by trial.
Currently it features:
* contains the matplotlib navigation toolbar
* allows user input of desired experiment trial number
* finds the multiplicity and completeness of data currently in the trial
* Displays this information graphically
* automatically updates the plot on a timer
* The plot can be saved to a file from the file menu

'''
from __future__ import division
import os
import random
import wx
import numpy as np
# The recommended way to use wx with mpl is with the WXAgg
# backend.
#
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

class BarsFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Completeness and Multiplicity at 2.5 Angstrom Resolution Progress'

    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)

        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()

        #set up the timer for testing
        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)
        self.redraw_timer.Start(1000)
        self.trial = None

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
        self.fnameLabel = wx.StaticText(self.panel, wx.ID_ANY, 'Trial Number ')
        self.fnameLabel.SetFont(wx.Font (12, wx.SWISS, wx.NORMAL, wx.BOLD))
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
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        # create a box for entering the trial number
        self.vbox.Add(self.fnameLabel, 0 , wx.ALL, 5)
        self.vbox.Add(self.textbox, 0, wx.ALL , 5)
        self.vbox.Add(self.SubmitButton, 0, wx.ALL , 5)


        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_figure(self):
        """ Redraws the figure
        """
        # multiplicity at 2.5 Angstroms
        self.mult = [random.uniform(0.0,4.0) for i in xrange(0,4)]
        self.plot_max = max(self.mult) + 1
        # completeness at 2.5 Angstroms
        self.data = [random.uniform(0.0,1.0)*self.plot_max for i in xrange(0,4)]
        pos = np.arange(len(self.data))+0.5 # the bar centers on the y-axis
        labels = ['Dark Isoform A', 'Dark Isoform B', '2F Isoform A', '2F Isoform B']
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()
        self.axes.barh(
            pos,
            self.data,
            align='center')
        self.axes.set_title('Trial: %s Completeness and Multiplicity at 2.5 Angstroms' % self.textbox.GetValue())
        self.axes.set_yticks(pos)
        self.axes.set_yticklabels(labels)
        self.axes.set_xlabel('Multiplicity')
        self.axes.set_xlim(0.0,self.plot_max)
        self.axes.axvline(x=self.mult[0], ymin=0.0, ymax=0.25, color='k', linestyle='dashed', linewidth=4)
        self.axes.axvline(x=self.mult[1], ymin=0.25, ymax=0.50, color='k', linestyle='dashed', linewidth=4)
        self.axes.axvline(x=self.mult[2], ymin=0.5, ymax=0.75, color='k', linestyle='dashed', linewidth=4)
        self.axes.axvline(x=self.mult[3], ymin=0.75, ymax=1.0, color='k', linestyle='dashed', linewidth=4)
        self.canvas.draw()

    def on_submit(self, event):
        print 'on_submit handler'
        self.trial = self.textbox.GetValue()
        print 'Textbox value is: %s'% self.trial
        self.on_redraw_timer()

    def on_redraw_timer(self,event=None):
        if self.trial is not None:
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
    app = wx.PySimpleApp()
    app.frame = BarsFrame()
    app.frame.Show()
    app.MainLoop()
