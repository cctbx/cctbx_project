from __future__ import division

#-----------------------------------------------------------------------
# more-or-less real-time plotting of Bragg peak count and XES detector
# skewness.
#-----------------------------------------------------------------------

from xfel.cxi.gfx import status_plot
import wxtbx.plots
from scitbx.array_family import flex
import libtbx.phil
from libtbx.utils import Usage, Sorry
import wx
import matplotlib.ticker as ticker
import time
import os

master_phil = libtbx.phil.parse("""
  trial_id = None
    .type = int
  t_wait = 8000
    .type = int
  hit_cutoff = 12
    .type = int
  average_window = 1000
    .type = int
""")

class TrialsPlotFrame (wxtbx.plots.plot_frame) :
  show_controls_default = False
  def __init__ (self, *args, **kwds) :
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)

  def create_plot_panel (self) :
    return TrialsPlot(
      parent=self,
      figure_size=(16,10))

  def show_plot(self):
    self.plot_panel.show_plot()

  def set_trial_id(self, trial_id):
    self.trial_id = trial_id
    self.plot_panel.trial_id = trial_id

class TrialsPlot (wxtbx.plots.plot_container) :
  def show_plot (self):
    assert (self.trial_id is not None)
    self.figure.clear()
    xmin = xmax = None

    import cxi_xdr_xes.cftbx.cspad_ana.db as cxidb
    db=cxidb.dbconnect()
    assert(db is not None and db.open)

    cursor = db.cursor()
    cursor.execute("SELECT eventstamp, data FROM cxi_braggs_front WHERE trial = %s"%self.trial_id)

    t1 = flex.double()
    bragg = flex.int()

    for eventstamp, data in cursor.fetchall():
      t1.append(float(eventstamp))
      bragg.append(int(data))

  
    xmin, xmax = min(t1), max(t1) #when add more axis, use the logic below
#    if (len(t1) > 0) :
#      xmin, xmax = min(t1), max(t1)
#    if (len(t2) > 0) :
#      if (xmin is not None) :
#        xmin, xmax = min(min(t2), xmin), max(max(t2), xmax)
#      else :
#        xmin, xmax = min(t2), max(t2)
#    perm = flex.sort_permutation(t3)
#    t3 = t3.select(perm)
#    hit_rate = hit_rate.select(perm)
    ax1 = self.figure.add_axes([0.1, 0.05, 0.8, 0.4])
#    ax2 = self.figure.add_axes([0.1, 0.45, 0.8, 0.15], sharex=ax1)
#    ax3 = self.figure.add_axes([0.1, 0.6, 0.8, 0.25], sharex=ax1)
    ax1.grid(True, color="0.75")
#    ax2.grid(True, color="0.75")
#    ax3.grid(True, color="0.75")
    ax1.plot(t1, bragg, 'd', color=[0.0,0.5,1.0])
#    ax2.plot(t3, hit_rate, 'o-', color=[0.0,1.0,0.0])
#    ax3.plot(t4, photon_counts, '^', color=[0.8,0.0,0.2])
    ax1.set_ylabel("# of Bragg spots")
#    ax2.set_ylabel("Hit rate (%)")
#    ax3.set_ylabel("XES photon count")
#    if (len(photon_counts) > 0) :
#      ax3.set_ylim(-1, max(photon_counts))
    ax1.set_xlim(xmin, xmax)
    ax1.set_xlabel("Time")
#    for ax in ax1, ax2, ax3:
#      if (ax is not ax1) :
#        for label in ax.get_xticklabels():
#          label.set_visible(False)
#      ax.get_yticklabels()[0].set_visible(False)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(status_plot.format_time))
    ax1.set_title("Detector analysis for trial %d" % self.trial_id) #was axis 3
    self.figure.autofmt_xdate()
    self.canvas.draw()
    self.parent.Refresh()

def run (args) :
  user_phil = []
  # TODO: replace this stuff with iotbx.phil.process_command_line_with_files
  # as soon as I can safely modify it
  for arg in args :
    #if (os.path.isdir(arg)) :
      #user_phil.append(libtbx.phil.parse("""status_dir=\"%s\"""" % arg))
    #elif (not "=" in arg) :
    if (not "=" in arg) :
      try :
        user_phil.append(libtbx.phil.parse("""trial_id=%d""" % int(arg)))
      except ValueError, e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = master_phil.fetch(sources=user_phil).extract()
  if (params.trial_id is None) :
    master_phil.show()
    raise Usage("trial_id must be defined (either trial_id=XXX, or the integer "+
      "ID alone).")
  #if (params.status_dir is None) :
    #master_phil.show()
    #raise Usage("status_dir must be defined!")
  #elif (not os.path.isdir(params.status_dir)) :
    #raise Sorry("%s does not exist or is not a directory!" % params.status_dir)
  assert (params.t_wait is not None) and (params.t_wait > 0)
  assert (params.hit_cutoff is not None) and (params.hit_cutoff > 0)
  assert (params.average_window is not None) and (params.average_window > 0)
  app = wx.App(0)
  frame = TrialsPlotFrame(None, -1, "Detector status for trial %d" %
    params.trial_id)
  frame.set_trial_id(params.trial_id)
  frame.show_plot()
  frame.Show()
  app.MainLoop()
