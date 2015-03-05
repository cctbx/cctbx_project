from __future__ import division
import wx
import wx.lib.plot as plot
import multiprocessing
import thread
import threading
import time
from pyana.event import Event
EVT_EXTERNAL_UPDATE = wx.PyEventBinder(wx.NewEventType(), 0)
import time


class ExternalUpdateEvent (wx.PyCommandEvent) :
  """XXX This class, along with the EVT_EXTERNAL_UPDATE instance
  should perhaps move into its own file?
  """

  def __init__ (self, eventType=EVT_EXTERNAL_UPDATE.evtType[0], id=0) :
    wx.PyCommandEvent.__init__(self, eventType, id)
    self.data = None
    self.number = None
    self.calib= None

class ShotFrame (wx.Frame) :

  def __init__ (self, *args, **kwds) :

    super(ShotFrame, self).__init__(*args, **kwds)
    self.panel = wx.Panel(self)
#    self.panel.SetBackgroundColour("yellow")
    self.mb = wx.MenuBar()
    self.Bind(EVT_EXTERNAL_UPDATE, self.OnExternalUpdate)
    menubar = wx.MenuBar()
    file = wx.Menu()
    file.Append(22, '&Quit', 'Exit Calculator')
    menubar.Append(file, '&File')
    self.SetMenuBar(menubar)
    wx.EVT_MENU(self, 20, self.OnClose)
    self.counter=0
    self.plotter = plot.PlotCanvas(self.panel)
    self.plotter.SetInitialSize(size=(500, 500))
    self.plotter.SetPosition((150,150))

  def OnExternalUpdate (self, event) :
    """The OnExternalUpdate() function updates the image and the title
    from @p event.
    """
    self.draw(event.data, event.number, event.calib)
    return
  def draw(self,data,number, calib):
        max=0
        min=10000000
        b=[]
        if calib is not None:
          x=calib[0]
          y=calib[1]
          z=calib[2]
        else:
          x=0;y=1;z=0
        for i in range(len(data)):
         E=i*i*x+i*y+z
         a=(E,data[i])
         b.append(a)
         if E>max: max=E
         if E<min: min=E

        # draw points as a line
        line = plot.PolyLine(b, colour='blue', width=1)
        # set up text, axis and draw
        gc = plot.PlotGraphics([line], "Event Number "+str(number), 'Energy', 'Intensity')
        self.plotter.Draw(gc,xAxis=(min,max))

  def OnClose(self, event):
    self.Close()

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
    check=False
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
    app = wx.App(0)
    self._frame = ShotFrame( None, -1, "Shot Spectrum", size=(800, 720))

    self._frame.Bind(wx.EVT_IDLE, self.OnIdle)

    self._frame.Show()

    self._init_lock.release()
    app.MainLoop()

    # Avoid deadlock where the send_data() function is waiting for the
    # semaphore after the frame has closed.
    self._next_semaphore.release()


  def send_data(self, data):
    """The send_data() function updates the wxPython application with
    @p img and @p title by sending it an ExternalUpdateEvent().  The
    function blocks until the event is processed."""
    event = ExternalUpdateEvent()
    event.data = data[0]
    event.number = data[1]
    event.calib= data[2]
    if self.isAlive():
      try:
        # Saturating the event queue makes the whole caboodle
        # uselessly unresponsive.  Therefore, block until idle events
        # are processed.
#         while self.isAlive(): #and not self._run_pause.IsToggled():
#          pass
#        print "event sent"
          self._frame.AddPendingEvent(event)
          self._is_idle = False
          while self.isAlive() and not self._is_idle:
           pass
      except Exception:
        pass

  def OnIdle(self, event):
    self._is_idle = True
    event.RequestMore()

def _xray_frame_process(queue):
  """The _xray_frame_process() function starts the viewer in a
  separate thread.  It then continuously reads data from @p queue and
  dispatches update events to the viewer.  The function returns when
  it reads a @c None object from @p queue or when the viewer thread
  has exited.
  """

  from Queue import Empty

  # Start the viewer's main loop in its own thread, and get the
  # interface for sending updates to the frame.
  thread = _XrayFrameThread()
  send_data = thread.send_data

  while True:
    try:
      data = queue.get(timeout=1)

      if data is None:
        if True:
          thread.join()
        else:
          thread.stop()
        return

      if not thread.isAlive():
        thread.join()
        return

#      if wait is not None:
#        time.sleep(wait)

      # All kinds of exceptions--not just PyDeadObjectError--may occur
      # if the viewer process exits during this call.  XXX This may be
      # dangerous!
      try:
        send_data(data)
      except Exception:
        pass
    except Empty:
      pass



class mod_spectra_view:

  def __init__(self, calib=None,
               **kwds):

    """

    @param address         Address string XXX Que?!
    @param dirname         Directory portion of output pickle file
    @param basename        Filename prefix of output pickle file
    @param dark_path       Path to input dark image
    @param dark_stddev     Path to input dark standard deviation
    """

#    linger = True # XXX Make configurable

    # Create a managed FIFO queue shared between the viewer and the
    # current process.  The current process will produce images, while
    # the viewer process will consume them.

    manager = multiprocessing.Manager()
    self._queue = manager.Queue()
    self._proc = multiprocessing.Process(target=_xray_frame_process,args=(self._queue,))
    self._proc.start()
    self.nv = 0
    if calib is not None:
      self.calib=calib.split(",")
      self.calib=map(float,self.calib)

  def beginjob(self, evt, env):
    print "Start"

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.
    Once self.nshots shots are accumulated, this function turns into
    a nop.
    @param evt Event data object, a configure object
    @param env Environment object
    """
    if (evt.get("skip_event")):
      return

    self.nv=self.nv+1
    sum1=evt.get("cctbx_spectra")
    while not self._queue.empty():
        if not self._proc.is_alive():
          evt.setStatus(Event.Stop)
          return
    while True:
        try:
          self._queue.put((sum1, self.nv, self.calib), timeout=1)
          break
        except Exception:
          pass

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

    try:
      self._queue.put(None)
    except Exception:
      pass
    print "endjob(): end of stream"
    self._proc.join()
