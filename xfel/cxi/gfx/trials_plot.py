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
  n_points = 1000
    .type = int
""")

class Run (object):
  def __init__(self, runId):
    self.runId = runId
    self.bragg_times = flex.double()
    self.braggs = flex.double()
    self.culled_braggs = flex.double()
    self.culled_bragg_times = flex.double()
    self.hit_rates_times = flex.double()
    self.hit_rates = flex.double()

  def width(self):
    return max(self.bragg_times)-min(self.bragg_times)

  def min(self):
    return min(self.bragg_times)

  def max(self):
    return max(self.bragg_times)

  def cull_braggs(self, count):
    if count <= 0:
      self.culled_braggs = self.braggs
      self.culled_bragg_times = self.bragg_times
      return

    self.culled_braggs = flex.double()
    self.culled_bragg_times = flex.double()

    window = int(len(self.bragg_times)/count)
    for i in range(count):
      value = self.braggs[i*window]
      time = self.bragg_times[i*window]
      for j in range(window):
        idx = (i*window)+j
        if self.braggs[idx] > value:
          value = self.braggs[idx]
          time = self.bragg_times[idx]
      self.culled_braggs.append(value)
      self.culled_bragg_times.append(time)

class TrialsPlotFrame (wxtbx.plots.plot_frame) :
  show_controls_default = False
  def __init__ (self, *args, **kwds) :
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)

    label = wx.StaticText(self.toolbar, -1, " Zoom: ", style=wx.ALIGN_CENTER)
    self.toolbar.AddControl(label)

    self.zoomSlider = wx.Slider(self.toolbar, size= (250, -1), minValue=0, maxValue=99)
    self.toolbar.AddControl(self.zoomSlider)
    self.Bind(wx.EVT_SCROLL, self.OnZoom, self.zoomSlider)

    label = wx.StaticText(self.toolbar, -1, " Pan: ", style=wx.ALIGN_CENTER)
    self.toolbar.AddControl(label)

    self.panSlider = wx.Slider(self.toolbar, size= (250, -1), minValue=1, maxValue=100)
    self.panSlider.SetValue(50)
    self.panSlider.Disable()
    self.toolbar.AddControl(self.panSlider)
    self.Bind(wx.EVT_SCROLL, self.OnPan, self.panSlider)

    self.timerCheck = wx.CheckBox(self.toolbar, -1, "Auto load new data")
    self.timerCheck.SetValue(True)
    self.toolbar.AddControl(self.timerCheck)

    self.zoom = 100
    self.pan = 50

    self.full_data_load = True

  def create_plot_panel (self) :
    return TrialsPlot(
      parent=self,
      figure_size=(16,10))

  def show_plot(self):
    self.plot_panel.show_plot()

  def set_params(self, params):
    self.params = params
    self.trial_id = params.trial_id

    self._timer = wx.Timer(owner=self)
    self.Bind(wx.EVT_TIMER, self.OnTimer)
    self._timer.Start(params.t_wait)

  def OnTimer(self, event):
    if(self.timerCheck.GetValue()):
      #t1 = time.time()
      if self.load_data(): pass
      self.show_plot()
      #t2 = time.time()
      #print "Updated in %.2fs" % (t2 - t1)

  def OnZoom(self, event):
    self.zoom = 100-event.GetPosition()
    if self.zoom < 100:
      self.panSlider.Enable()
    else:
      self.panSlider.Disable()

    self.cull_braggs()

    self.plot_panel.show_plot()

  def OnPan(self, event):
    self.pan = event.GetPosition()
    self.plot_panel.show_plot()

  # Returns true if new data was loaded, otherwise false
  def load_data (self):
    assert (self.trial_id is not None)

    import cxi_xdr_xes.cftbx.cspad_ana.db as cxidb
    db=cxidb.dbconnect()
    assert(db is not None and db.open)

    # retrieve the run IDs in this trial
    cursor = db.cursor()
    cursor.execute("SELECT DISTINCT(run) FROM cxi_braggs_front WHERE trial = %s ORDER BY run"%self.trial_id)
    if(self.full_data_load):
      self.runs = []

    new_data = False

    for runId in cursor.fetchall():
      if self.full_data_load:
        run = Run(int(runId[0]))
        self.runs.append(run)
      else:
        for runtest in self.runs:
          foundit = False
          if runtest.runId == int(runId[0]):
            foundit = True
            run = runtest
            break
        if not foundit:
          run = Run(int(runId[0]))

      if self.full_data_load or not hasattr(run, "latest_entry_id"):
        cursor.execute("SELECT id, eventstamp, hitcount FROM cxi_braggs_front WHERE trial = %s AND run = %s ORDER BY eventstamp"%(self.trial_id,run.runId))
      else:
        cursor.execute("SELECT id, eventstamp, hitcount FROM cxi_braggs_front WHERE trial = %s AND run = %s AND id > %s ORDER BY eventstamp"%(self.trial_id,run.runId,run.latest_entry_id ))

      hit_counter = self.params.average_window

      ids = flex.int()

      for id, eventstamp, data in cursor.fetchall():
        run.bragg_times.append(float(eventstamp))
        run.braggs.append(int(data))
        ids.append(id)

        if (len(run.bragg_times) >= self.params.average_window) and hit_counter >= self.params.average_window:
          start = len(run.bragg_times) - self.params.average_window
          window = run.braggs[start:]
          isel = (window > self.params.hit_cutoff).iselection()
          ratio = float(len(isel)) / float(self.params.average_window)
          run.hit_rates_times.append(run.bragg_times[-1])
          run.hit_rates.append(ratio*100)
          hit_counter = 0
        else :
          hit_counter += 1

      if len(ids) > 0:
        run.latest_entry_id = max(ids)
        new_data = True


#    self.xmin, self.xmax = min(t1), max(t1) #when add more axis, use the logic below
#    if (len(t1) > 0) :
#      xmin, xmax = min(t1), max(t1)
#    if (len(t2) > 0) :
#      if (xmin is not None) :
#        xmin, xmax = min(min(t2), xmin), max(max(t2), xmax)
#      else :
#        xmin, xmax = min(t2), max(t2)

    self.total_width = 0
    for run in self.runs:
      perm = flex.sort_permutation(run.hit_rates_times)
      run.hit_rates_times = run.hit_rates_times.select(perm)
      run.hit_rates = run.hit_rates.select(perm)

      self.total_width += run.width()

    self.cull_braggs()

    self.full_data_load = False
    return new_data

  def cull_braggs(self):
    for run in self.runs:
      run.cull_braggs(int(self.params.n_points*run.width()/self.total_width))


class TrialsPlot (wxtbx.plots.plot_container) :
  def show_plot(self):
    """
    We assume that the runs are sorted in ascending order and that they do not overlap in time

    The final graph will have a number of time units displayed equal to a percentage
    of the total number of time units in all the runs, not including gaps.  Each run
    may or may not be displayed depending on the zoom and pan.

    total_width: how many time units total in these runs, not including gaps
    xmin: smallest time value in the earliest run
    xmax: largest time value in the lastest run
    newwidth: number of time units to display on xaxis total over all runs
    newmid: a point on the total time units scale within the zoom and pan limits. Usually
     the center, but since it could cause x to fall out of bounds it is adjusted
    newxmin: how many time units on the left of the graph are not displayed, OR, the first
     time unit to be displayed
    newxmax: the first time unit not to be displayed
    """
    total_width = self.GetParent().total_width

    newwidth = total_width * (self.GetParent().zoom / 100)
    newmid = total_width * (self.GetParent().pan/100)
    newxmin = newmid - (newwidth/2)
    newxmax = newxmin + newwidth

    if newxmin < 0:
      newxmin = 0
      newxmax = newwidth
    elif newxmax > total_width:
      newxmax = total_width
      newxmin = newxmax - newwidth

    assert newxmin >= 0 and newxmin <= total_width

    #print "**** Zoom: %s, pan: %s, total_width: %s, newwidth: %s, newmid: %s, newxmin: %s, newxmax: %s" \
    #  %(self.GetParent().zoom,self.GetParent().pan,total_width,newwidth,newmid,newxmin,newxmax)

    left = 0
    width_so_far = 0
    self.figure.clear()
    for run in self.GetParent().runs:
      right = left + run.width()

      if right < newxmin or left > newxmax:
        left += run.width()
        #print "Not showing run %s"%run.runId
        continue

      if left < newxmin:
        xmin = run.min() + (newxmin - left)
      else:
        xmin = run.min()

      if right > newxmax:
        xmax = run.min() + (newxmax - left)
      else:
        xmax = run.max()

      #print "Run: %s, run.width(): %s, left: %s, right: %s, run.min(): %s, run.max(): %s, xmin: %s, xmax: %s, width_so_far: %s, xmax-xmin: %s" \
        #%(run.runId,run.width(),left,right,run.min(),run.max(),xmin,xmax,width_so_far,xmax-xmin)

      ax1 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.05, (xmax-xmin)/newwidth, 0.4])
      ax2 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.45, (xmax-xmin)/newwidth, 0.15], sharex=ax1)
      left += run.width()
      width_so_far += (xmax-xmin)

#    ax3 = self.figure.add_axes([0.1, 0.6, 0.8, 0.25], sharex=ax1)
      ax1.grid(True, color="0.75")
      ax2.grid(True, color="0.75")
#    ax3.grid(True, color="0.75")
      ax1.plot(run.culled_bragg_times, run.culled_braggs, 'd', color=[0.0,0.5,1.0])
      ax2.plot(run.hit_rates_times, run.hit_rates, 'o-', color=[0.0,1.0,0.0])
#    ax3.plot(t4, photon_counts, '^', color=[0.8,0.0,0.2])
      ax1.set_ylabel("# of Bragg spots")
      ax2.set_ylabel("Hit rate (%)")
#    ax3.set_ylabel("XES photon count")
#    if (len(photon_counts) > 0) :
#      ax3.set_ylim(-1, max(photon_counts))
      ax1.set_xlim(xmin, xmax)
      ax1.set_xlabel("Time")
#    for ax in ax1, ax2, ax3:
      for ax in ax1, ax2:
        if (ax is not ax1) :
          for label in ax.get_xticklabels():
            label.set_visible(False)
        ax.get_yticklabels()[0].set_visible(False)

      ax2.xaxis.set_major_formatter(ticker.FuncFormatter(status_plot.format_time))
      ax2.set_title("Run number: %d" % run.runId) #was axis 3

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
  assert (params.t_wait is not None) and (params.t_wait > 0)
  assert (params.hit_cutoff is not None) and (params.hit_cutoff > 0)
  assert (params.average_window is not None) and (params.average_window > 0)
  assert (params.n_points is not None) # zero or less means display all the points
  app = wx.App(0)
  frame = TrialsPlotFrame(None, -1, "Detector status for trial %d" %
      params.trial_id)
  frame.set_params(params)
  frame.load_data()
  frame.show_plot()
  frame.Show()
  app.MainLoop()
