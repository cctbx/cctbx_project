# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$
"""Interactive image viewer for CSPad images

XXX Known issues, wishlist:

  * Is it slow?  Yes!

  * Radial distribution plot, requested by Jan F. Kern

  * Coordinate and resolution tool-tips in sub-pixel zoom.  Can we
    choose the colour automagically?
"""
from __future__ import absolute_import, division, print_function

__version__ = "$Revision$"

import multiprocessing
import thread
import threading
import time

from dxtbx.format.FormatPYunspecified import FormatPYunspecified
from iotbx.detectors.npy import NpyImage
from rstbx.viewer.frame import XrayFrame
from xfel.cxi.cspad_ana import common_mode, cspad_tbx
from xfel.cxi.cspad_ana import skip_event_flag

class _ImgDict(NpyImage):
  """Minimal iotbx detector class for in-memory dictionary
  representation of NpyImage images.  This class must be defined at
  the top level of the module so that it can be pickled.
  """

  def __init__(self, data, parameters):
    # self.vendortype is required to guess the beam centre convention.

    super(_ImgDict, self).__init__('')
    self.parameters = parameters
    self.vendortype = 'npy_raw'
    self.bin_safe_set_data(data)

  def readHeader(self):
    pass

  def read(self):
    pass


class _Format(FormatPYunspecified):
  """Minimal Format class for in-memory dxtbx representation of
  FormatPYunspecified images.  This class must be defined at the top
  level of the module so that it can be pickled.
  """

  def __init__(self, **kwargs):
    """The _Format constructor builds a dxtbx Format instance from the
    supplied keyworded arguments.  It should be equivalent to the
    FormatPYunspecified constructor.
    """

    from copy import deepcopy

    from dxtbx.model.beam import Beam, BeamFactory
    from dxtbx.model.detector import Detector, DetectorFactory
    from dxtbx.model.goniometer import Goniometer, GoniometerFactory
    from dxtbx.model import Scan, ScanFactory

    from spotfinder.applications.xfel import cxi_phil
    from iotbx.detectors.cspad_detector_formats import detector_format_version

    # From Format.__init__().
    self._goniometer_factory = GoniometerFactory
    self._detector_factory = DetectorFactory
    self._beam_factory = BeamFactory
    self._scan_factory = ScanFactory

    # From FormatPYunspecified._start().  Split the keyworded
    # arguments into a parameter dictionary suitable for
    # iotbx.detectors.npy.NpyImage and a separate image data object.
    parameters = dict(
      BEAM_CENTER_X=kwargs['PIXEL_SIZE'] * kwargs['BEAM_CENTER'][0],
      BEAM_CENTER_Y=kwargs['PIXEL_SIZE'] * kwargs['BEAM_CENTER'][1],
      CCD_IMAGE_SATURATION=kwargs['SATURATED_VALUE'],
      DISTANCE=kwargs['DISTANCE'],
      OSC_RANGE=0,
      OSC_START=0,
      PIXEL_SIZE=kwargs['PIXEL_SIZE'],
      SATURATED_VALUE=kwargs['SATURATED_VALUE'],
      SIZE1=kwargs['DATA'].focus()[0],
      SIZE2=kwargs['DATA'].focus()[1],
      WAVELENGTH=kwargs['WAVELENGTH'])
    self.detectorbase = _ImgDict(kwargs['DATA'], parameters)

    # Attempt to apply tile translations only for detectors that
    # support it.
    version_lookup = detector_format_version(
      kwargs['DETECTOR_ADDRESS'],
      kwargs['TIME_TUPLE'][0])
    if version_lookup is not None:
      params = cxi_phil.cxi_versioned_extract(
        "distl.detector_format_version=" + version_lookup)

      # Necessary to keep the phil parameters for subsequent calls to
      # get_tile_manager().
      horizons_phil = params.persist.commands

      self.detectorbase.translate_tiles(horizons_phil)
      self.detectorbase.horizons_phil_cache = deepcopy(horizons_phil)

    # From Format.setup().
    goniometer_instance = self._goniometer()
    assert(isinstance(goniometer_instance, Goniometer))
    self._goniometer_instance = goniometer_instance

    detector_instance = self._detector()
    assert(isinstance(detector_instance, Detector))
    self._detector_instance = detector_instance

    beam_instance = self._beam()
    assert(isinstance(beam_instance, Beam))
    self._beam_instance = beam_instance

    scan_instance = self._scan()
    assert(isinstance(scan_instance, Scan))
    self._scan_instance = scan_instance


class _XrayFrameThread(threading.Thread):
  """The _XrayFrameThread class allows MainLoop() to be run as a
  thread, which is necessary because all calls to wxPython must be
  made from the same thread that originally imported wxPython.

  This is all based on "Running MainLoop in a separate thread",
  http://wiki.wxpython.org/MainLoopAsThread.
  """

  def __init__(self):
    """The thread is started automatically on initialisation.
    self.run() will initialise self.frame and release self._init_lock.
    """

    super(_XrayFrameThread, self).__init__()
    self.setDaemon(1)
    self._init_lock = threading.Lock()
    self._next_semaphore = threading.Semaphore()
    self._start_orig = self.start
    self._frame = None
    self.start = self._start_local

    self._init_lock.acquire()
    self.start()


  def _start_local(self):
    """The _start_local() function calls the run() function through
    self._start_orig, and exists only after self._init_lock has been
    released.  This eliminates a race condition which could cause
    updates to be sent to a non-existent frame.
    """

    self._start_orig()
    self._init_lock.acquire()


  def run(self):
    """The run() function defines the frame and starts the main loop.
    self._init_lock is released only when all initialisation is done.

    Whatever thread is the current one when wxWindows is initialised
    is what it will consider the "main thread." For wxPython 2.4 that
    happens when wxPython.wx is imported the first time.  For 2.5 it
    will be when the wx.App object is created.
    """

    import wx
    from wxtbx import bitmaps
    app = wx.App(0)
    self._bitmap_pause = bitmaps.fetch_icon_bitmap('actions', 'stop')
    self._bitmap_run = bitmaps.fetch_icon_bitmap('actions', 'runit')
    self._frame = XrayFrame(None, -1, "X-ray image display", size=(800, 720))

    self._frame.Bind(wx.EVT_IDLE, self.OnIdle)

    self.setup_toolbar(self._frame.toolbar)
    self._frame.Show()

    self._init_lock.release()
    app.MainLoop()

    # Avoid deadlock where the send_data() function is waiting for the
    # semaphore after the frame has closed.
    self._next_semaphore.release()


  def send_data(self, img, title):
    """The send_data() function updates the wxPython application with
    @p img and @p title by sending it an ExternalUpdateEvent().  The
    function blocks until the event is processed."""

    from rstbx.viewer.frame import ExternalUpdateEvent

    event = ExternalUpdateEvent()
    event.img = img
    event.title = title
    if self.isAlive():
      try:
        # Saturating the event queue makes the whole caboodle
        # uselessly unresponsive.  Therefore, block until idle events
        # are processed.
        while self.isAlive() and not self._run_pause.IsToggled():
          pass
        self._frame.AddPendingEvent(event)
        self._is_idle = False
        while self.isAlive() and not self._is_idle:
          pass
      except Exception:
        pass


  def setup_toolbar(self, toolbar):
    import wx
    from wxtbx import icons

    toolbar.ClearTools()

    btn = toolbar.AddLabelTool(
      id=wx.ID_ANY,
      label="Settings",
      bitmap=icons.advancedsettings.GetBitmap(),
      shortHelp="Settings",
      kind=wx.ITEM_NORMAL)
    self._frame.Bind(wx.EVT_MENU, self._frame.OnShowSettings, btn)

    btn = toolbar.AddLabelTool(
      id=wx.ID_ANY,
      label="Zoom",
      bitmap=icons.search.GetBitmap(),
      shortHelp="Zoom",
      kind=wx.ITEM_NORMAL)
    self._frame.Bind(wx.EVT_MENU, self._frame.OnZoom, btn)

    # Reset the normal bitmap after the tool has been created, so that
    # it will update on the next event.  See also OnPauseRun()
    self._run_pause = toolbar.AddCheckLabelTool(
      id=wx.ID_ANY,
      label="Run/Pause",
      bitmap=self._bitmap_run,
      shortHelp="Run/Pause")
    self._run_pause.SetNormalBitmap(self._bitmap_pause)
    self._frame.Bind(wx.EVT_MENU, self.OnPauseRun, self._run_pause)


  def OnIdle(self, event):
    self._is_idle = True
    event.RequestMore()


  def OnPauseRun(self, event):
    if self._run_pause.IsToggled():
      self._run_pause.SetNormalBitmap(self._bitmap_run)
    else:
      self._run_pause.SetNormalBitmap(self._bitmap_pause)


  def stop(self):
    from wx import CloseEvent
    self._frame.AddPendingEvent(CloseEvent())


def _xray_frame_process(queue, linger=True, wait=None):
  """The _xray_frame_process() function starts the viewer in a
  separate thread.  It then continuously reads data from @p queue and
  dispatches update events to the viewer.  The function returns when
  it reads a @c None object from @p queue or when the viewer thread
  has exited.
  """

  from Queue import Empty
  import rstbx.viewer

  # Start the viewer's main loop in its own thread, and get the
  # interface for sending updates to the frame.
  thread = _XrayFrameThread()
  send_data = thread.send_data

  while True:
    try:
      payload = queue.get(timeout=1)

      if payload is None:
        if linger:
          thread.join()
        else:
          thread.stop()
        return

      if not thread.isAlive():
        thread.join()
        return

      if wait is not None:
        time.sleep(wait)

      # All kinds of exceptions--not just PyDeadObjectError--may occur
      # if the viewer process exits during this call.  XXX This may be
      # dangerous!
      try:
        send_data(rstbx.viewer.image(payload[0]), payload[1])
      except Exception:
        pass
    except Empty:
      pass


class mod_view(common_mode.common_mode_correction):
  """XXX
  """

  def __init__(self,
               address,
               n_collate   = None,
               n_update    = 120,
               common_mode_correction = "none",
               wait=None,
               photon_counting=False,
               sigma_scaling=False,
               **kwds):
    """The mod_view class constructor XXX.

    @param address         Full data source address of the DAQ device
    @param calib_dir       Directory with calibration information
    @param common_mode_correction The type of common mode correction to apply
    @param dark_path       Path to input average dark image
    @param dark_stddev     Path to input standard deviation dark
                           image, required if @p dark_path is given
    @param wait            Minimum time (in seconds) to wait on the current
                           image before moving on to the next
    @param n_collate       Number of shots to average, or <= 0 to
                           average all shots
    @param n_update        Number of shots between updates
    """

    super(mod_view, self).__init__(
      address=address,
      common_mode_correction=common_mode_correction,
      **kwds)

    self.detector = cspad_tbx.address_split(address)[0]
    self.nvalid   = 0
    self.ncollate = cspad_tbx.getOptInteger(n_collate)
    self.nupdate  = cspad_tbx.getOptInteger(n_update)
    self.photon_counting = cspad_tbx.getOptBool(photon_counting)
    self.sigma_scaling = cspad_tbx.getOptBool(sigma_scaling)
    if (self.ncollate is None):
      self.ncollate = self.nupdate
    if (self.ncollate > self.nupdate):
      self.ncollate = self.nupdate
      self.logger.warning("n_collate capped to %d" % self.nupdate)

    linger = True # XXX Make configurable
    wait = cspad_tbx.getOptFloat(wait)

    # Create a managed FIFO queue shared between the viewer and the
    # current process.  The current process will produce images, while
    # the viewer process will consume them.
    manager = multiprocessing.Manager()
    self._queue = manager.Queue()
    self._proc = multiprocessing.Process(
      target=_xray_frame_process, args=(self._queue, linger, wait))
    self._proc.start()

    self.n_shots = 0


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    XXX Since the viewer is now running in a parallel process, the
    averaging here is now the bottleneck.

    @param evt Event data object, a configure object
    @param env Environment object
    """
    from pyana.event import Event

    self.n_shots += 1

    super(mod_view, self).event(evt, env)
    if evt.status() != Event.Normal or evt.get('skip_event'): # XXX transition
      return

    # Get the distance for the detectors that should have it, and set
    # it to NaN for those that should not.
    if self.detector == 'CxiDs1' or \
       self.detector == 'CxiDsd' or \
       self.detector == 'XppGon':
      distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
      if distance is None:
        self.nfail += 1
        self.logger.warning("event(): no distance, shot skipped")
        evt.put(skip_event_flag(), "skip_event")
        return
    else:
      distance = float('nan')

    if not self._proc.is_alive():
      evt.setStatus(Event.Stop)

    # Early return if the next update to the viewer is more than
    # self.ncollate shots away.  XXX Since the common_mode.event()
    # function does quite a bit of processing, the savings are
    # probably not so big.
    next_update = (self.nupdate - 1) - (self.nshots - 1) % self.nupdate
    if (self.ncollate > 0 and next_update >= self.ncollate):
      return

    if self.sigma_scaling:
      self.do_sigma_scaling()

    if self.photon_counting:
      self.do_photon_counting()

    # Trim the disabled section from the Sc1 detector image.  XXX This
    # is a bit of a kludge, really.
#    if (self.address == "CxiSc1-0|Cspad2x2-0"):
#      self.cspad_img = self.cspad_img[185:2 * 185, :]

    # Update the sum of the valid images, starting a new collation if
    # appropriate.  This guarantees self.nvalid > 0.
    if (self.nvalid == 0 or self.ncollate > 0 and self.nvalid >= self.ncollate):
      self.img_sum = self.cspad_img
      self.nvalid  = 1
    else:
      self.img_sum += self.cspad_img
      self.nvalid  += 1

    # Update the viewer to display the current average image, and
    # start a new collation, if appropriate.
    if (next_update == 0):
      from time import localtime, strftime

      time_str = strftime("%H:%M:%S", localtime(evt.getTime().seconds()))
      title = "r%04d@%s: average of %d last images on %s" \
          % (evt.run(), time_str, self.nvalid, self.address)

      # See also mod_average.py.
      device = cspad_tbx.address_split(self.address)[2]
      if device == 'Cspad':
        beam_center = self.beam_center
        pixel_size = cspad_tbx.pixel_size
        saturated_value = cspad_tbx.cspad_saturated_value

      elif device == 'marccd':
        beam_center = tuple(t // 2 for t in self.img_sum.focus())
        pixel_size = 0.079346
        saturated_value = 2**16 - 1

      # Wait for the viewer process to empty the queue before feeding
      # it a new image, and ensure not to hang if the viewer process
      # exits.  Because of multithreading/multiprocessing semantics,
      # self._queue.empty() is unreliable.
      fmt = _Format(BEAM_CENTER=beam_center,
                      DATA=self.img_sum / self.nvalid,
                      DETECTOR_ADDRESS=self.address,
                      DISTANCE=distance,
                      PIXEL_SIZE=pixel_size,
                      SATURATED_VALUE=saturated_value,
                      TIME_TUPLE=cspad_tbx.evt_time(evt),
                      WAVELENGTH=self.wavelength)

      while not self._queue.empty():
        if not self._proc.is_alive():
          evt.setStatus(Event.Stop)
          return
      while True:
        try:
          self._queue.put((fmt, title), timeout=1)
          break
        except Exception:
          pass

      if (self.ncollate > 0):
        self.nvalid = 0

  #signature for pyana:
  #def endjob(self, env):

  #signature for psana:
  #def endjob(self, evt, env):

  def endjob(self, obj1, obj2=None):
    """The endjob() function terminates the viewer process by sending
    it a @c None object, and waiting for it to finish.

    @param evt Event object (psana only)
    @param env Environment object
    """

    if obj2 is None:
      env = obj1
    else:
      evt = obj1
      env = obj2

    super(mod_view, self).endjob(env)
    try:
      self._queue.put(None)
    except Exception:
      pass
    self.logger.info("endjob(): end of stream")
    self._proc.join()
