# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id$
"""Interactive image viewer for CSPad images

XXX Known issues, wishlist:

  * Is it slow?  Even changing between two cached images takes a lot
    of time.

  * Radial distribution plot, requested by Jan F. Kern

  * Sub-pixel zoom like in xdisp, maybe with coordinate tool-tips on
    mouse-over

  * Should use Python's multiprocessing module instead of thread and
    threading.

XXX
"""

__version__ = "$Revision$"

import thread
import threading

from iotbx.detectors.detectorbase import DetectorImageBase
from rstbx.viewer.frame import XrayFrame, ExternalUpdateEvent
from xfel.cxi.cspad_ana import common_mode
from xfel.cxi.cspad_ana import cspad_tbx


class cxi_dict(DetectorImageBase):
  """Minimal detector class for in-memory dictionary representation of
  images.
  """

  def __init__(self, data, parameters):
    """self.vendortype is required to guess the beam centre
    convention.
    """

    super(cxi_dict, self).__init__("")
    self.parameters = parameters
    self.vendortype = "npy_raw"
    self.bin_safe_set_data(data)


  def readHeader(self):
    pass


  def read(self):
    pass


def ImageFactory(img_dict):
  """The ImageFactory() function splits @p img_dict into a parameter
  dictionary suitable for iotbx.detectors.DetectorImageBase, and a
  separate image data object.  The function relies on constants from
  cspad_tbx.  See also iotbx.detectors.NpyImage"""
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
  return (cxi_dict(data, parameters))

from iotbx import detectors
detectors.ImageFactory = ImageFactory


class DataDispatch(object):
  """Interface for sending images to the wxPython application."""

  def __init__(self, parent):
    self.parent = parent


  def send_data(self, img, title):
    """The send_data() function creates and sends the event.
    """
    event       = ExternalUpdateEvent()
    event.img   = img
    event.title = title
    self.parent.AddPendingEvent(event)


class XrayFrame_thread(threading.Thread):
  """The XrayFrame_thread class allows Run MainLoop() to be run as a
  thread, which is necessary because all calls to wxPython must be
  made from the same thread that originally imported wxPython.

  This is all based on "Running MainLoop in a separate thread",
  http://wiki.wxpython.org/MainLoopAsThread.
  """

  def __init__(self, n_frames):
    """The thread is started automatically on initialisation.
    self.run() will populate self.frames and release self.lock.
    """

    super(XrayFrame_thread, self).__init__()
    self.setDaemon(1)
    self.start_orig = self.start
    self.start      = self.start_local
    self.frames     = []
    self.n_frames   = n_frames
    self.lock       = threading.Lock()
    self.lock.acquire()
    self.start()


  def run(self):
    """The run() function defines the frame and starts the main loop.
    self.lock is released only when all initialisation is done.

    Whatever thread is the current one when wxWindows is initialised
    is what it will consider the "main thread." For wxPython 2.4 that
    happens when wxPython.wx is imported the first time.  For 2.5 it
    will be when the wx.App object is created.
    """

    import wx
    app = wx.App(0)
    for i in xrange(self.n_frames):
      self.frames.append(
        XrayFrame(None,
                  -1,
                  "X-ray image display [%d/%d]" % (i + 1, self.n_frames),
                  size = (800, 720)))
      self.frames[i].Show()

    self.lock.release()
    app.MainLoop()


  def start_local(self):
    """The start_local() function calls the run() function through
    self.start_orig, and exists only after self.lock has been
    released.  This eliminates a race condition which could cause
    updates to be sent to non-existent frame."""
    self.start_orig()
    self.lock.acquire()


class mod_view(common_mode.common_mode_correction):
  """XXX
  """

  def __init__(self,
               address,
               n_collate   = None,
               n_update    = 120,
               n_view      = "1/1",
               common_mode_correction = "none",
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
    @param n_collate       Number of shots to average, or <= 0 to
                           average all shots
    @param n_update        Number of shots between updates
    @param n_view          Index of current viewer and number of
                           simultaneous viewers in configuration file,
                           on the form "n/m"
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

    # Start the viewer's main loop in its own thread, and get the
    # interface for sending updates to the frame.  If there already is
    # a running thread, use the appropriate frame from that.
    n_view         = cspad_tbx.getOptString(n_view).split("/")
    self.send_data = None
    for t in threading.enumerate():
      if (isinstance(t, XrayFrame_thread)):
        self.send_data = DataDispatch(t.frames[int(n_view[0]) - 1]).send_data
        break
    if (self.send_data is None):
      self.send_data = DataDispatch(
        XrayFrame_thread(int(n_view[1])).frames[0]).send_data


  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    super(mod_view, self).event(evt, env)
    if (evt.get("skip_event")):
      return

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
      import rstbx.viewer

      img = rstbx.viewer.image(dict(
        BEAM_CENTER = self.beam_center,
        DATA        = self.img_sum / self.nvalid,
        DISTANCE    = self.distance,
        WAVELENGTH  = self.wavelength))
      time_str  = strftime("%H:%M:%S", localtime(evt.getTime().seconds()))
      title = "r%04d@%s: average of %d last images" \
          % (evt.run(), time_str, self.nvalid)

      self.send_data(img, title)
      if (self.ncollate > 0):
        self.nvalid = 0
