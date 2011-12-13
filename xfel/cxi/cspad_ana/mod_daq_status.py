import logging
import threading
import wx

from xfel.cxi.gfx import status_plot
from xfel.cxi.cspad_ana import cspad_tbx

class StatusFrame_thread(threading.Thread):
  """The XrayFrame_thread class allows Run MainLoop() to be run as a
  thread, which is necessary because all calls to wxPython must be
  made from the same thread that originally imported wxPython.

  This is all based on "Running MainLoop in a separate thread",
  http://wiki.wxpython.org/MainLoopAsThread.
  """
  def __init__(self):
    threading.Thread.__init__(self) # XXX super()?
    self.setDaemon(1)
    self.start_orig = self.start
    self.start      = self.start_local
    self.frame      = None
    self.lock       = threading.Lock()
    self.lock.acquire()
    self.start()

  def run(self):
    import wx
    app   = wx.App(0)
    frame = status_plot.StatusFrame(None, -1, "CXI experiment status")
    frame.Show()
    self.frame = frame
    self.lock.release()
    app.MainLoop()

  def start_local(self):
    """The start_local() function calls the run() function through
    self.start_orig, and exists only after self.lock has been
    released.  This eliminates a race condition which could cause
    updates to be sent to non-existent frame."""
    self.start_orig()
    self.lock.acquire()

class mod_daq_status (object) :
  def __init__ (self, address) :
    self.address = address
    self.initialize()
    self.logger = logging.getLogger(__name__)
    self.logger.setLevel(logging.INFO)
    self.display_thread = StatusFrame_thread()
    self.window = self.display_thread.frame
    self.run_id = None

  def initialize (self) :
    self.nfail = 0
    self.nshots = 0
    self._t = []
    self._wavelength = []
    self._det_z = []
    self._si_foil = []

  def beginjob(self, evt, env):
    env.update(evt)
    self.initialize()
    self.run_id = evt.run()
    event = status_plot.RunNumberEvent(self.run_id)
    wx.PostEvent(self.window, event)

  def event (self, evt, env) :
    if (evt.get("skip_event")) :
      return
    det_z = cspad_tbx.env_detz(env)
    si_foil = cspad_tbx.env_sifoil(env)
    wavelength = cspad_tbx.evt_wavelength(evt)
    t = evt.getTime()
    s = None
    self.nshots += 1
    if (t is not None):
      s = t.seconds() + (t.nanoseconds() / 1000000000)
    else :
      self.nfail += 1
      self.logger.warn("event(): no timestamp, shot skipped")
      evt.put(True, "skip_event")
      return
    if (det_z is None):
      self.nfail += 1
      self.logger.warn("event(): no distance, shot skipped")
      evt.put(True, "skip_event")
      return
    if (si_foil is None):
      self.nfail += 1
      self.logger.warn("event(): no Si-foil thickness, shot skipped")
      evt.put(True, "skip_event")
      return
    if (wavelength is None):
      self.nfail += 1
      self.logger.warn("event(): no wavelength, shot skipped")
      evt.put(True, "skip_event")
      return
    if (not isinstance(s, float)) :
      raise RuntimeError("Wrong type for 's': %s" % type(s).__name__)
    if (not (isinstance(si_foil, float) or isinstance(si_foil, int))) :
      raise RuntimeError("Wrong type for 'si_foil': %s"% type(si_foil).__name__)
    self._t.append(s)
    self._si_foil.append(si_foil)
    self._wavelength.append(wavelength)
    self._det_z.append(det_z)
    if (self.nshots % 120 == 0) :
      self.update_plot()

  def update_plot (self) :
    """
    Post an update event with current plot values to redraw the window.
    """
    event = status_plot.UpdateEvent(self._t, self._det_z, self._si_foil,
      self._wavelength)
    wx.PostEvent(self.window, event)

  def endjob (self, env) :
    print "END OF RUN"
    wx.PostEvent(self.window, status_plot.SaveImageEvent())

    # Uncomment to close the frame immediately.  Otherwise, wouldn't
    # it be nice if the window's title bar indicated that the run has
    # ended, so that one wouldn't have to watch the controlling
    # terminal?  XXX It may be safer to post an event than to call
    # Close() directly.
#    self.window.Close()

    self.display_thread.join()
