
########################################################################
# (excerpted and paraphrased) email from JH 11-23-11:
#
# fields:
# - DAQ det-z in mm (floating point in range [-400, 10])
# - mm of Si-foil in the beam (integer, in range [0, 20460], typically
#   less than 1000)
# - Wavelength, should be tightly distributed (floating point)
#
# The det-z from the data acquisition system (field 5) is related to the
# sample--detector distance.  Because we may want to change the exact
# relationship during the experiment, plotting the det-z value from the DAQ
# should give more consistent results (and I'd think we'd only like to know
# whether it changes or not).  If it changes, it could (semi-)continuously
# vary over the entire permitted range.
#
# The thickness of the attenuation Si-foil is integer (actually we can only
# ever see certain values, ask again if you'd like to know how this works).
# These aren't continuous changes, and they'll typically jump back and forth
# as the sheets of foil are inserted and removed.
#
# I expect the wavelength to have a sigma of around 0.001 A (even though we
# could see deviations as large as 0.1 A).  Basically, whatever wavelength we
# see for the first shout shouldn't change much.
########################################################################

import wxtbx.plots
import wx
import time

def format_time(x, pos=None):
  lt = time.localtime(x)
  return time.strftime("%H:%M:%S", lt)

def draw_plot (figure, t, det_z, laser01, laser04, laser04_power, si_foil,
               wavelength, bragg=None, run_id=None) :
  assert (len(t) == len(si_foil))
  import matplotlib.ticker as ticker
  left, width = 0.1, 0.8
  rect1 = [left, 0.65, width, 0.15] # detector Z
  rect2 = [left, 0.50, width, 0.15] # mm of Si-foil
  rect3 = [left, 0.35, width, 0.15] # wavelength
  #rect4 = [left, 0.05, width, 0.4] # bragg spots
  rect5 = [left, 0.20, width, 0.15] # Laser #1 & #4
  rect6 = [left, 0.05, width, 0.15] # Laser #4 power
  ax1 = figure.add_axes(rect1)
  ax2 = figure.add_axes(rect2, sharex=ax1)
  ax3 = figure.add_axes(rect3, sharex=ax1)
  #ax4 = figure.add_axes(rect4, sharex=ax1)
  ax51 = figure.add_axes(rect5, sharex=ax1)
  ax52 = ax51.twinx()
  ax6 = figure.add_axes(rect6, sharex=ax1)
  ax1.grid(True, color="0.75")
  ax2.grid(True, color="0.75")
  ax3.grid(True, color="0.75")
  #ax4.grid(True, color="0.75")
  ax51.grid(True, color="0.75")
  ax52.grid(True, color="0.75")
  ax6.grid(True, color="0.75")
  try :
    ax1.plot(t, det_z, linewidth=1, color=[1.0, 0.0, 0.0])
  except ValueError, e :
    print "error plotting det_z:", e
    print len(t), len(det_z)
  try :
    ax2.plot(t, si_foil, linewidth=1, color=[0.2, 0.4, 0.6])
  except ValueError, e :
    print "error plotting si_foil:", e
    print len(t), len(si_foil)
  ax3.plot(t, wavelength, linewidth=1, color=[0.0,1.0,0.5])
  #ax4.bar(t, bragg, facecolor=[0.0,0.5,1.0], edgecolor=[0.0,0.5,1.0])
  #ax4.plot(t, bragg, linewidth=1, color=[0.0,0.5,1.0])
  ax51.plot(t, laser01, linewidth=1, color=[1.0,0.0,0.5])
  ax52.plot(t, laser04, linewidth=1, color=[0.0,1.0,0.5])
  ax6.plot(t, laser04_power, linewidth=1, color=[0.0,1.0,0.5])
  ax1.set_ylabel("Detector Z (mm)")
  ax2.set_ylabel("Si foil (mm)")
  ax3.set_ylabel("Wavelength (A)")
  #ax4.set_ylabel("Bragg spots")
  ax51.set_ylabel("Laser #1", color=[1.0,0.0,0.5])
  ax52.set_ylabel("Laser #4", color=[0.0,1.0,0.5])
  ax6.set_ylabel("Laser #4 power")
  ax1.set_xlabel("Time")
  for ax in ax1, ax2, ax3, ax51, ax52, ax6: #, ax4:
    if (ax is not ax6) :
      for label in ax.get_xticklabels():
        label.set_visible(False)
      ax.get_yticklabels()[0].set_visible(False)
  ax6.xaxis.set_major_formatter(ticker.FuncFormatter(format_time))
  if (run_id is not None) :
    ax1.set_title("CXI experiment status - run %d" % run_id)
  else :
    ax1.set_title("CXI experiment status - RUN NUMBER NOT SET")
  figure.autofmt_xdate()

UPDATE_ID = wx.NewId()
RUN_NUMBER_ID = wx.NewId()
SAVE_IMG_ID = wx.NewId()

class UpdateEvent (wx.PyEvent) :
  """
  Contains new data for plotting (one x series and three y series).
  """
  def __init__ (self, t, det_z, laser01, laser04, laser04_power, si_foil,
                wavelength) :
    self.data = (t, det_z, laser01, laser04, laser04_power, si_foil, wavelength)
    wx.PyEvent.__init__(self)
    self.SetEventType(UPDATE_ID)

class RunNumberEvent (wx.PyEvent) :
  """
  Sets the current run ID, which will appear in the frame title bar and the
  plot.
  """
  def __init__ (self, run_id) :
    self.run_id = run_id
    wx.PyEvent.__init__(self)
    self.SetEventType(RUN_NUMBER_ID)

class SaveImageEvent (wx.PyEvent) :
  """
  Signals that an image of the current plot should be saved to disk (the file
  name will be something like rID_status.png, where ID is the run number.
  """
  def __init__ (self) :
    wx.PyEvent.__init__(self)
    self.SetEventType(SAVE_IMG_ID)

class StatusFrame (wxtbx.plots.plot_frame) :
  show_controls_default = False
  def __init__ (self, *args, **kwds) :
    self.run_id = None
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)
    self.Connect(-1, -1, UPDATE_ID, self.OnUpdate)
    self.Connect(-1, -1, RUN_NUMBER_ID, self.OnSetRunNumber)
    self.Connect(-1, -1, SAVE_IMG_ID, self.OnSaveImage)

  def create_plot_panel (self) :
    return StatusPlot(
      parent=self,
      figure_size=(16,10))

  def OnUpdate (self, event) :
    t, det_z, laser01, laser04, laser04_power, si_foil, wavelength = event.data
    self.plot_panel.show_plot(t, det_z, laser01, laser04, laser04_power,
                              si_foil, wavelength, run_id=self.run_id)

  def OnSetRunNumber (self, event) :
    self.run_id = event.run_id
    self.SetTitle("CXI experiment status - run %d" % self.run_id)

  def OnSaveImage (self, event) :
    print "StatusFrame.OnSaveImage(event)"
    assert (self.run_id is not None)
    file_name = "r%d_status.png" % self.run_id
    self.plot_panel.figure.savefig(file_name, format="png")
    print "SAVED PLOT: %s" % file_name

class StatusPlot (wxtbx.plots.plot_container) :
  def show_plot (self, *args, **kwds) :
    self.figure.clear()
    draw_plot(self.figure, *args, **kwds)
    self.canvas.draw()
    self.parent.Refresh()

def draw_plot_to_file (file_name, *args, **kwds) :
  """
  For offline plotting, which turns out to be necessary when running the
  status plot through pyana.
  """
  import matplotlib
  import matplotlib.figure
  from matplotlib.backends.backend_agg import FigureCanvasAgg
  figure = matplotlib.figure.Figure((16,10))
  canvas = FigureCanvasAgg(figure)
  draw_plot(figure, *args, **kwds)
  canvas.draw()
  figure.savefig(file_name, format="png")
