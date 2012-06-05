
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
  status_dir = None
    .type = path
  run_id = None
    .type = int
  t_wait = 8000
    .type = int
  hit_cutoff = 12
    .type = int
  average_window = 1000
    .type = int
""")

class DetectorPlotFrame (wxtbx.plots.plot_frame) :
  show_controls_default = False
  def __init__ (self, *args, **kwds) :
    self.params = None
    self._watch_dir = None
    self._watch_files = []
    self._line_offsets = {}
    self._t1 = flex.double()
    self._t2 = flex.double()
    self._t3 = flex.double()
    self._t4 = flex.double()
    self._hit_sample_last = 0
    self._bragg = flex.int()
    self._skewness = flex.double()
    self._photon_counts = flex.double()
    self._hit_ratio = flex.double()
    wxtbx.plots.plot_frame.__init__(self, *args, **kwds)

  def create_plot_panel (self) :
    return DetectorPlot(
      parent=self,
      figure_size=(16,10))

  def set_run_id (self, run_id) :
    self.plot_panel.set_run_id(run_id)

  def set_watch_dir (self, dir_name, params) :
    assert os.path.isdir(dir_name)
    self.params = params
    self._watch_dir = dir_name
    self.update_from_logs()
    self._timer = wx.Timer(owner=self)
    self.Bind(wx.EVT_TIMER, self.OnTimer)
    self._timer.Start(self.params.t_wait)

  def find_logs (self) :
    assert (self._watch_dir is not None)
    current_dir = os.path.join(self._watch_dir, "stdout")
    print "Current directory: %s" % current_dir
    for file_name in os.listdir(current_dir) :
      if (file_name.endswith(".out")) :
        full_path = os.path.join(current_dir, file_name)
        print "Adding %s to list of files to monitor" % full_path
        f = open(full_path)
        self._watch_files.append(f)

  def update_from_logs (self, force_update_hits=False) :
    if (len(self._watch_files) == 0) :
      self.find_logs()
    if (len(self._watch_files) > 0) :
      #print "Collecting new data @ %s" % time.strftime("%H:%M:%S",
      #  time.localtime())
      for fh in self._watch_files :
        for line in fh.readlines() :
          if ("BRAGG" in line) :
            fields1 = line.split(":")
            fields2 = fields1[-1].strip().split()
            self._t1.append(float(fields2[1]))
            self._bragg.append(int(fields2[2]))
            hit_point_min = self._hit_sample_last+self.params.average_window
            if (not force_update_hits) and (len(self._t1) > hit_point_min) :
              self.update_hit_rate()
              self._hit_sample_last = len(self._t1)
          elif ("SKEW" in line) :
            fields1 = line.split(":")
            fields2 = fields1[-1].strip().split()
            self._t2.append(float(fields2[1]))
            self._skewness.append(float(fields2[2]))
          elif ("N_PHOTONS" in line) :
            fields1 = line.split(":")
            fields2 = fields1[-1].strip().split()
            self._t4.append(float(fields2[1]))
            self._photon_counts.append(float(fields2[2]))
      if (force_update_hits) :
        self.update_hit_rate()
      self.plot_panel.show_plot(
        t1=self._t1,
        bragg=self._bragg,
        t2=self._t2,
        skewness=self._skewness,
        t3=self._t3,
        hit_rate=self._hit_ratio,
        t4=self._t4,
        photon_counts=self._photon_counts,
      )

  def update_hit_rate (self) :
    if (len(self._t1) >= self.params.average_window) :
      start = len(self._t1) - self.params.average_window
      window = self._bragg[start:]
      isel = (window > self.params.hit_cutoff).iselection()
      ratio = float(len(isel)) / float(self.params.average_window)
      self._t3.append(self._t1[-1])
      self._hit_ratio.append(ratio*100)

  def OnTimer (self, event) :
    t1 = time.time()
    self.update_from_logs(True)
    t2 = time.time()
    #print "Updated in %.2fs" % (t2 - t1)

  def OnSave (self, event) :
    self.plot_panel.save_png()

class DetectorPlot (wxtbx.plots.plot_container) :
  def set_run_id (self, run_id) :
    self.run_id = run_id

  def show_plot (self, t1, bragg, t2, skewness, t3, hit_rate,
                 t4, photon_counts) :
    assert (self.run_id is not None)
    self.figure.clear()
    xmin = xmax = None
    if (len(t1) > 0) :
      xmin, xmax = min(t1), max(t1)
    if (len(t2) > 0) :
      if (xmin is not None) :
        xmin, xmax = min(min(t2), xmin), max(max(t2), xmax)
      else :
        xmin, xmax = min(t2), max(t2)
    perm = flex.sort_permutation(t3)
    t3 = t3.select(perm)
    hit_rate = hit_rate.select(perm)
    ax1 = self.figure.add_axes([0.1, 0.05, 0.8, 0.4])
    ax2 = self.figure.add_axes([0.1, 0.45, 0.8, 0.15], sharex=ax1)
    ax3 = self.figure.add_axes([0.1, 0.6, 0.8, 0.25], sharex=ax1)
    ax1.grid(True, color="0.75")
    ax2.grid(True, color="0.75")
    ax3.grid(True, color="0.75")
    ax1.plot(t1, bragg, 'd', color=[0.0,0.5,1.0])
    ax2.plot(t3, hit_rate, 'o-', color=[0.0,1.0,0.0])
    ax3.plot(t4, photon_counts, '^', color=[0.8,0.0,0.2])
    ax1.set_ylabel("# of Bragg spots")
    ax2.set_ylabel("Hit rate (%)")
    ax3.set_ylabel("XES photon count")
    if (len(photon_counts) > 0) :
      ax3.set_ylim(-1, max(photon_counts))
    ax1.set_xlim(xmin, xmax)
    ax1.set_xlabel("Time")
    for ax in ax1, ax2, ax3:
      if (ax is not ax1) :
        for label in ax.get_xticklabels():
          label.set_visible(False)
      ax.get_yticklabels()[0].set_visible(False)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(status_plot.format_time))
    ax3.set_title("Detector analysis for run %d" % self.run_id)
    self.figure.autofmt_xdate()
    self.canvas.draw()
    self.parent.Refresh()

  def save_png (self) :
    if (getattr(self, "run_id", None) is not None) :
      file_name = "run%d_detector_status.png" % self.run_id
      self.figure.savefig(file_name, format="png")
      print "Saved image to %s" % os.path.abspath(file_name)
    else :
      print "Can't save an image until run ID is set"

def run (args) :
  user_phil = []
  # TODO: replace this stuff with iotbx.phil.process_command_line_with_files
  # as soon as I can safely modify it
  for arg in args :
    if (os.path.isdir(arg)) :
      user_phil.append(libtbx.phil.parse("""status_dir=\"%s\"""" % arg))
    elif (not "=" in arg) :
      try :
        user_phil.append(libtbx.phil.parse("""run_id=%d""" % int(arg)))
      except ValueError, e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = master_phil.fetch(sources=user_phil).extract()
  if (params.run_id is None) :
    master_phil.show()
    raise Usage("run_id must be defined (either run_id=XXX, or the integer "+
      "ID alone).")
  if (params.status_dir is None) :
    master_phil.show()
    raise Usage("status_dir must be defined!")
  elif (not os.path.isdir(params.status_dir)) :
    raise Sorry("%s does not exist or is not a directory!" % params.status_dir)
  assert (params.t_wait is not None) and (params.t_wait > 0)
  assert (params.hit_cutoff is not None) and (params.hit_cutoff > 0)
  assert (params.average_window is not None) and (params.average_window > 0)
  app = wx.App(0)
  frame = DetectorPlotFrame(None, -1, "Detector status for run %d" %
    params.run_id)
  frame.set_run_id(params.run_id)
  frame.set_watch_dir(params.status_dir, params)
  frame.Show()
  app.MainLoop()
