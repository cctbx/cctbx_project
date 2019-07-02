from __future__ import absolute_import, division, print_function

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
import operator
import math
from six.moves import range

master_phil = libtbx.phil.parse("""
  trial_id = None
    .type = int
  t_wait = 8000
    .type = int
  hit_cutoff = 16
    .type = int
  average_window = 1000
    .type = int
  n_points = 1000
    .type = int
  display_time = 1800
    .type = int
  run_num = None
    .type = int
  run_min = None
    .type = int
  run_max = None
    .type = int
""")

class Run (object):
  def __init__(self, runId):
    self.runId = runId
    self.hits_count = 0
    self.bragg_times = flex.double()
    self.braggs = flex.double()
    self.culled_braggs = flex.double()
    self.culled_bragg_times = flex.double()
    self.hit_rates_times = flex.double()
    self.hit_rates = flex.double()

    self.distances = flex.double()
    self.culled_distances = flex.double()

    self.sifoils = flex.double()
    self.culled_sifoils = flex.double()

    self.wavelengths = flex.double()
    self.culled_wavelengths = flex.double()

    self.indexed = flex.bool()
    self.culled_indexed = flex.bool()

  def width(self):
    return max(self.bragg_times)-min(self.bragg_times)

  def min(self):
    return min(self.bragg_times)

  def max(self):
    return max(self.bragg_times)

  def cull_braggs(self, count):
    if count == 0:
      window = 0
    else:
      window = len(self.bragg_times)/count

    if count <= 0 or window < 1:
      self.culled_braggs = self.braggs
      self.culled_bragg_times = self.bragg_times
      self.culled_distances = self.distances
      self.culled_sifoils = self.sifoils
      self.culled_wavelengths = self.wavelengths
      self.culled_indexed = self.indexed
      return

    self.culled_braggs = flex.double()
    self.culled_bragg_times = flex.double()
    self.culled_distances = flex.double()
    self.culled_sifoils = flex.double()
    self.culled_wavelengths = flex.double()
    self.culled_indexed = flex.bool()

    for i in range(count):
      braggs =  self.braggs[i*int(window):(i+1)*int(window)]
      idx = int(i*window) + flex.max_index(braggs)
      self.culled_braggs     .append(self.braggs[idx])
      self.culled_bragg_times.append(self.bragg_times[idx])
      self.culled_distances  .append(self.distances[idx])
      self.culled_sifoils    .append(self.sifoils[idx])
      self.culled_wavelengths.append(self.wavelengths[idx])
      self.culled_indexed    .append(self.indexed[idx])

  def recalc_hits(self, windowLen, hit_cutoff):

    self.hit_rates_times = flex.double()
    self.hit_rates = flex.double()
    self.hits_count = 0
    if len(self.braggs) <= 0 or windowLen <= 0: return

    for i in range(int(math.floor(len(self.braggs)/windowLen))):
      window = self.braggs[i*windowLen:(i+1)*windowLen]
      isel = (window >= hit_cutoff).iselection()
      ratio = float(len(isel)) / float(windowLen)
      self.hit_rates.append(ratio*100)
      self.hit_rates_times.append(self.bragg_times[(i*windowLen)+int(math.floor((windowLen/2)))])
      self.hits_count += len(isel)

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

    self.timelockCheck = wx.CheckBox(self.toolbar, -1, "Display last %s minutes: ")
    self.timelockCheck.SetValue(False)
    self.timelockCheck.newly_checked = False
    self.toolbar.AddControl(self.timelockCheck)
    self.timelockCheck.Bind(wx.EVT_CHECKBOX, self.OnTimeLockCheck, self.timelockCheck)

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

    self.timelockCheck.SetLabel("Display last %s minutes: "%(self.params.display_time/60))

  def OnTimeLockCheck(self, event):
    self.timelockCheck.newly_checked = True

  def OnTimer(self, event):
    if(self.timerCheck.GetValue()):
      if self.load_data() or self.timelockCheck.newly_checked:
        if self.timelockCheck.GetValue():
          self.panSlider.Disable()
          self.zoomSlider.Disable()
          if self.total_width > 0:
            self.zoom = 100 * self.params.display_time / self.total_width
            self.pan = 100

            if self.zoom > 100: self.zoom = 100

        else:
          self.zoomSlider.Enable()

        self.timelockCheck.newly_checked = False
        self.show_plot()
      else:
        print("No new data")

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
    ttop = time.time()
    print("Loading data...")
    assert (self.trial_id is not None)

    import cxi_xdr_xes.cftbx.cspad_ana.db as cxidb
    db=cxidb.dbconnect()
    assert(db is not None and db.open)

    # retrieve the run IDs in this trial
    #t1 = time.time()
    cursor = db.cursor()
    #cursor.execute("SELECT DISTINCT(run) FROM %s WHERE trial = %s"%(cxidb.table_name,self.trial_id))
    cmd = "SELECT DISTINCT(run) FROM %s WHERE trial = %s"
    if self.params.run_num is not None:
      extra = " AND run = %s"%self.params.run_num
    elif self.params.run_min is not None and self.params.run_max is not None:
      extra = " AND run >= %s AND run <= %s"%(self.params.run_min, self.params.run_max)
    else:
      extra = " ORDER BY run DESC LIMIT 5"
    cursor.execute(cmd%(cxidb.table_name,self.trial_id) + extra)
    #t2 = time.time()
    #print "Runs queried in %.2fs" % (t2 - t1)

    if(self.full_data_load):
      self.runs = []
    if(len(self.runs) > 5):
      self.runs = self.runs[-5:]

    new_data = False

    for runId in cursor.fetchall():
      if self.full_data_load:
        run = Run(int(runId[0]))
        self.runs.append(run)
      else:
        foundit=False
        for runtest in self.runs:
          if runtest.runId == int(runId[0]):
            foundit = True
            run = runtest
            break
        if not foundit:
          print("New run: %s"%runId)
          run = Run(int(runId[0]))
          self.runs.append(run)

      #t1 = time.time()
      #print "Loading data from run %s" % (run.runId)
      if self.full_data_load or not hasattr(run, "latest_entry_id"):
        print("Full load")
        cursor.execute("SELECT id, eventstamp, hitcount, distance, sifoil, wavelength, indexed FROM %s \
          WHERE trial = %s AND run = %s ORDER BY eventstamp"%(cxidb.table_name,self.trial_id,run.runId))
      else:
        print("Partial load")
        cursor.execute("SELECT id, eventstamp, hitcount, distance, sifoil, wavelength, indexed FROM %s \
          WHERE trial = %s AND run = %s AND id > %s ORDER BY eventstamp"%(cxidb.table_name,self.trial_id,run.runId,run.latest_entry_id ))

      #t2 = time.time()
      #print "Query ran in %.2fs" % (t2 - t1)

      ids = flex.int()

      for id, eventstamp, hitcount, distance, sifoil, wavelength, indexed in cursor.fetchall():
        run.bragg_times.append(float(eventstamp))
        run.braggs.append(int(hitcount))
        ids.append(id)

        run.distances.append(float(distance))
        run.sifoils.append(float(sifoil))
        run.wavelengths.append(float(wavelength))
        run.indexed.append(bool(indexed))

      if len(ids) > 0:
        run.latest_entry_id = max(ids)
        new_data = True
        run.recalc_hits(self.params.average_window, self.params.hit_cutoff)


    self.total_width = 0
    for run in self.runs:
      perm = flex.sort_permutation(run.hit_rates_times)
      run.hit_rates_times = run.hit_rates_times.select(perm)
      run.hit_rates = run.hit_rates.select(perm)

      self.total_width += run.width()

    self.cull_braggs()

    #self.full_data_load = False #always do a full load
    self.runs.sort(key=operator.attrgetter('runId'))
    tbot = time.time()
    print("Data loaded in %.2fs" % (tbot - ttop))
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
    runs = self.GetParent().runs
    if len(runs) <= 0: return

    t1 = time.time()
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
    braggsmax = max(flex.max(r.culled_braggs) for r in runs)
    braggsmin = min(flex.min(r.culled_braggs) for r in runs)
    distsmax = max(flex.max(r.culled_distances) for r in runs)
    distsmin = min(flex.min(r.culled_distances) for r in runs)
    sifomax = max(flex.max(r.culled_sifoils) for r in runs)
    sifomin = min(flex.min(r.culled_sifoils) for r in runs)
    wavemax = max(flex.max(r.culled_wavelengths) for r in runs)
    wavemin = min(flex.min(r.culled_wavelengths) for r in runs)

    #above tricks don't work for hit rates as they can be empty if the run is new
    goodruns = []
    for run in runs:
      if len(run.hit_rates) > 0: goodruns.append(run)
    if len(goodruns) > 0:
      hitsmax = max(flex.max(r.hit_rates) for r in goodruns)
      hitsmin = min(flex.min(r.hit_rates) for r in goodruns)
    else:
      hitsmax = hitsmin = 0

    first_run = True
    for run in runs:
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

      ax1 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.05, 0.9*(xmax-xmin)/newwidth, 0.4])
      ax2 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.45, 0.9*(xmax-xmin)/newwidth, 0.2], sharex=ax1)
      ax3 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.65, 0.9*(xmax-xmin)/newwidth, 0.1], sharex=ax1)
      ax4 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.75, 0.9*(xmax-xmin)/newwidth, 0.1], sharex=ax1)
      ax5 = self.figure.add_axes([0.05+(0.9*width_so_far/newwidth), 0.85, 0.9*(xmax-xmin)/newwidth, 0.1], sharex=ax1)
      left += run.width()
      width_so_far += (xmax-xmin)

      ax1.grid(True, color="0.75")
      ax2.grid(True, color="0.75")
      ax3.grid(True, color="0.75")
      ax4.grid(True, color="0.75")
      ax5.grid(True, color="0.75")
      ax1.plot(run.culled_bragg_times.select(run.culled_indexed),
               run.culled_braggs.select(run.culled_indexed), 'd', color=[0.0,1.0,0.0])
      ax1.plot(run.culled_bragg_times.select(~run.culled_indexed),
               run.culled_braggs.select(~run.culled_indexed), 'd', color=[0.0,0.5,1.0])
      ax2.plot(run.hit_rates_times, run.hit_rates, 'o-', color=[0.0,1.0,0.0])
      ax3.plot(run.culled_bragg_times, run.culled_wavelengths, '^', color=[0.8,0.0,0.2])
      ax4.plot(run.culled_bragg_times, run.culled_sifoils, '<', color=[0.8,0.0,0.2])
      ax5.plot(run.culled_bragg_times, run.culled_distances, '>', color=[0.8,0.0,0.2])
      ax1.set_ylabel("# of Bragg spots")
      ax2.set_ylabel("Hit rate (%)")
      ax3.set_ylabel("WaveL")
      ax4.set_ylabel("SiFoils(mm)")
      ax5.set_ylabel("Dist (mm)")
      ax1.set_xlim(xmin, xmax)
      ax1.set_ylim(braggsmin, braggsmax)
      ax2.set_ylim(hitsmin, hitsmax)
      ax3.set_ylim(wavemin, wavemax)
      ax4.set_ylim(sifomin-10, sifomax+10)
      ax5.set_ylim(distsmin-3, distsmax+3)
      ax1.set_xlabel("Time")
      for ax in ax1, ax2, ax3, ax4, ax5:
        if (ax is not ax1) :
          for label in ax.get_xticklabels():
            label.set_visible(False)
        ax.get_yticklabels()[0].set_visible(False)
        if not first_run:
          ax.get_yaxis().set_visible(False)

      ax1.xaxis.set_major_formatter(ticker.FuncFormatter(status_plot.format_time))
      ax3.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))
      ax5.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.0f"))
      ax5.set_title("%d:%d/%d:%.1f%% I:%d"%(run.runId, run.hits_count, len(run.braggs), 100*run.hits_count/len(run.braggs),run.indexed.count(True)))

      labels = ax1.get_xticklabels()
      for label in labels:
        label.set_rotation(30)

      first_run = False

    self.figure.autofmt_xdate()
    self.canvas.draw()
    self.parent.Refresh()

    t2 = time.time()
    print("Plotted in %.2fs" % (t2 - t1))

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
      except ValueError as e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
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
  assert (params.display_time is not None) and (params.display_time > 0)

  app = wx.App(0)
  frame = TrialsPlotFrame(None, -1, "Detector status for trial %d" %
      params.trial_id)
  frame.set_params(params)
  frame.load_data()
  frame.show_plot()
  frame.Show()
  app.MainLoop()
