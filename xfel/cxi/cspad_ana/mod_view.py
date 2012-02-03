# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$
"""Interactive image viewer for CSPad images

XXX Known issues, wishlist:

  * Is it slow?  Even changing between two cached images takes a lot
    of time.

  * Radial distribution plot, requested by Jan F. Kern

  * Sub-pixel zoom like in xdisp, maybe with coordinate tool-tips on
    mouse-over.  Can we choose the colour automagically?  Zoom window
    not updated with hold!

XXX
"""

__version__ = "$Revision$"

import multiprocessing
import thread
import threading
import sys
import time

from iotbx.detectors.detectorbase import DetectorImageBase
from rstbx.viewer.frame import XrayFrame
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx


class _CxiDict(DetectorImageBase):
  """Minimal detector class for in-memory dictionary representation of
  images.
  """

  def __init__(self, data, parameters):
    """self.vendortype is required to guess the beam centre
    convention.
    """

    super(_CxiDict, self).__init__("")
    self.parameters = parameters
    self.vendortype = "npy_raw"
    self.bin_safe_set_data(data)


  def readHeader(self):
    pass


  def read(self):
    pass


def _image_factory(img_dict):
  """The ImageFactory() function splits @p img_dict into a parameter
  dictionary suitable for iotbx.detectors.DetectorImageBase, and a
  separate image data object.  The function relies on constants from
  cspad_tbx.  See also iotbx.detectors.NpyImage.
  """

  data       = img_dict["DATA"]
  parameters = dict(
    BEAM_CENTER_X        = cspad_tbx.pixel_size * img_dict["BEAM_CENTER"][0],
    BEAM_CENTER_Y        = cspad_tbx.pixel_size * img_dict["BEAM_CENTER"][1],
    CCD_IMAGE_SATURATION = cspad_tbx.dynamic_range,
    DISTANCE             = img_dict["DISTANCE"],
    OSC_RANGE            = 0,
    OSC_START            = 0,
    PIXEL_SIZE           = cspad_tbx.pixel_size,
    SATURATED_VALUE      = cspad_tbx.dynamic_range,
    SIZE1                = data.focus()[0],
    SIZE2                = data.focus()[1],
    WAVELENGTH           = img_dict["WAVELENGTH"])
  return (_CxiDict(data, parameters))

from iotbx import detectors
detectors.ImageFactory = _image_factory


class _XrayFrameThread(threading.Thread):
  """The _XrayFrameThread class allows MainLoop() to be run as a
  thread, which is necessary because all calls to wxPython must be
  made from the same thread that originally imported wxPython.

  This is all based on "Running MainLoop in a separate thread",
  http://wiki.wxpython.org/MainLoopAsThread.
  """

  def __init__(self, hold=False):
    """The thread is started automatically on initialisation.
    self.run() will initialise self.frame and release self._init_lock.
    """

    super(_XrayFrameThread, self).__init__()
    self.setDaemon(1)
    self._hold = hold
    self._init_lock = threading.Lock()
    self._next_lock = threading.Semaphore()
    self._start_orig = self.start
    self._frame = None
    self.start = self._start_local

    self._init_lock.acquire()
    self.start()


  def _forward (self, event):
    """The _forward() function implements an alternate event handler
    for the viewer's "Next"-button.  It increments the semaphore by
    one.
    """

    self._next_lock.release()


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
    app = wx.App(0)
    self._frame = XrayFrame(None, -1, "X-ray image display", size=(800, 720))
    if (self._hold):
      self._frame.Bind(wx.EVT_MENU, self._forward, id=wx.ID_FORWARD)
    self._frame.Show()

    self._init_lock.release()
    app.MainLoop()

    # Avoid deadlock where the send_data() function is waiting for the
    # semaphore after the frame has closed.
    self._next_lock.release()


  def send_data(self, img, title):
    """The send_data() function updates the wxPython application with
    @p img and @p title by sending it an ExternalUpdateEvent()."""

    from rstbx.viewer.frame import ExternalUpdateEvent
    from wx import PyDeadObjectError

    event = ExternalUpdateEvent()
    event.img = img
    event.title = title
    if (self._hold):
      # Decrement the counter by one and return immediately if the
      # counter is larger than zero.
      self._next_lock.acquire()
    if (self.isAlive()):
      try:
        self._frame.AddPendingEvent(event)
      except PyDeadObjectError:
        pass


def _xray_frame_process(pipe, hold=False, wait=None):
  """The _xray_frame_process() function starts the viewer in a
  separate thread.  It then continuously reads data from @p pipe and
  dispatches update events to the viewer.  The function returns when
  it reads a @c None object from @p pipe or when the viewer thread has
  exited.
  """

  import rstbx.viewer

  # Start the viewer's main loop in its own thread, and get the
  # interface for sending updates to the frame.
  thread = _XrayFrameThread(hold)
  send_data = thread.send_data

  while (True):
    payload = pipe.recv()
    if (payload is None or not thread.isAlive()):
      break
    if wait is not None:
      time.sleep(wait)
    send_data(rstbx.viewer.image(payload[0]), payload[1])


class mod_view(common_mode.common_mode_correction):
  """XXX
  """

  def __init__(self,
               address,
               n_collate   = None,
               n_update    = 120,
               common_mode_correction = "none",
               hold=False,
               wait=None,
               photon_counting=False,
               sigma_scaling=False,
               **kwds):
    """The mod_view class constructor XXX.

    @param address         Address string XXX Que?!
    @param calib_dir       Directory with calibration information
    @param common_mode_correction The type of common mode correction to apply
    @param dark_path       Path to input average dark image
    @param dark_stddev     Path to input standard deviation dark
                           image, required if @p dark_path is given
    @param hold            Whether to wait for user input after each
                           displayed image or not
    @param wait            Length of time (in seconds) to wait on the current
                           image before moving on to the next
    @param n_collate       Number of shots to average, or <= 0 to
                           average all shots
    @param n_update        Number of shots between updates
    """

    super(mod_view, self).__init__(
      address=address,
      common_mode_correction=common_mode_correction,
      **kwds)

    self.nvalid   = 0
    self.ncollate = cspad_tbx.getOptInteger(n_collate)
    self.nupdate  = cspad_tbx.getOptInteger(n_update)
    self.photon_counting = cspad_tbx.getOptBool(photon_counting)
    self.sigma_scaling = cspad_tbx.getOptBool(sigma_scaling)
    if (self.ncollate is None):
      self.ncollate = self.nupdate
    if (self.ncollate > self.nupdate):
      self.ncollate = self.nupdate
      self.logger.warn("n_collate capped to %d" % self.nupdate)

    wait = cspad_tbx.getOptFloat(wait)
    hold = cspad_tbx.getOptBool(hold)
    # Create a unidirectional pipe and hand its read end to the viewer
    # process.  The write end is kept for sending updates.
    pipe_recv, self._pipe = multiprocessing.Pipe(False)
    self._proc = multiprocessing.Process(
      target=_xray_frame_process, args=(pipe_recv, hold, wait))
    self._proc.start()


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    XXX Since the viewer is now running in a parallel process, the
    averaging here is now the bottleneck.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_view, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    if (not self._proc.is_alive()):
      sys.exit(self._proc.exitcode)

    # Early exit if the next update to the viewer is more than
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
      title = "r%04d@%s: average of %d last images" \
          % (evt.run(), time_str, self.nvalid)

      self._pipe.send((dict(
          BEAM_CENTER = self.beam_center,
          DATA = self.img_sum / self.nvalid,
          DISTANCE = self.distance,
          WAVELENGTH = self.wavelength), title))

      if (self.ncollate > 0):
        self.nvalid = 0


  def endjob(self, env):
    """The endjob() terminates the viewer process by sending it a @c
    None object, and waiting for it to finish.

    @param env Environment object
    """

    super(mod_view, self).endjob(env)
    self._pipe.send(None)
    self._proc.join()
