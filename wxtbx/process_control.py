
from __future__ import absolute_import
from wx.lib.agw import pyprogress
import wx
from libtbx import thread_utils
from libtbx.utils import Sorry
import threading

JOB_START_ID = wx.NewId()
LOG_UPDATE_ID = wx.NewId()
CALLBACK_ID = wx.NewId()
JOB_EXCEPTION_ID = wx.NewId()
JOB_KILLED_ID = wx.NewId()
JOB_COMPLETE_ID = wx.NewId()
DOWNLOAD_COMPLETE_ID = wx.NewId()

class SubprocessEvent (wx.PyEvent) :
  event_id = None

  def __init__ (self, data, **kwds) :
    self.data = data
    self.__dict__.update(kwds)
    wx.PyEvent.__init__(self)
    self.SetEventType(self.event_id)

class JobStartEvent (SubprocessEvent) :
  event_id = JOB_START_ID

class LogEvent (SubprocessEvent) :
  event_id = LOG_UPDATE_ID

class JobExceptionEvent (SubprocessEvent) :
  event_id = JOB_EXCEPTION_ID

class JobKilledEvent (SubprocessEvent) :
  event_id = JOB_KILLED_ID

class JobCompleteEvent (SubprocessEvent) :
  event_id = JOB_COMPLETE_ID

class CallbackEvent (SubprocessEvent) :
  event_id = CALLBACK_ID

class DownloadCompleteEvent (SubprocessEvent) :
  event_id = DOWNLOAD_COMPLETE_ID

def setup_stdout_logging_event (window, OnPrint) :
  window.Connect(-1, -1, LOG_UPDATE_ID, OnPrint)

def setup_process_gui_events (
    window,
    OnStart=None,
    OnPrint=None,
    OnUpdate=None,
    OnExcept=None,
    OnAbort=None,
    OnComplete=None) :
  if OnStart is not None :
    assert hasattr(OnStart, "__call__")
    window.Connect(-1, -1, JOB_START_ID, OnStart)
  if OnPrint is not None :
    assert hasattr(OnPrint, "__call__")
    window.Connect(-1, -1, LOG_UPDATE_ID, OnPrint)
  if OnUpdate is not None :
    assert hasattr(OnUpdate, "__call__")
    window.Connect(-1, -1, CALLBACK_ID, OnUpdate)
  if OnExcept is not None :
    assert hasattr(OnExcept, "__call__")
    window.Connect(-1, -1, JOB_EXCEPTION_ID, OnExcept)
  if OnAbort is not None :
    assert hasattr(OnAbort, "__call__")
    window.Connect(-1, -1, JOB_KILLED_ID, OnAbort)
  if OnComplete is not None :
    assert hasattr(OnComplete, "__call__")
    window.Connect(-1, -1, JOB_COMPLETE_ID, OnComplete)

class event_agent (object) :
  def __init__ (self, window, **kwds) :
    self.window = window
    self.__dict__.update(kwds)

  def get_kwds (self) :
    return {}

  def callback_start (self, data) :
    kwds = self.get_kwds()
    event = JobStartEvent(data, **kwds)
    wx.PostEvent(self.window, event)

  def callback_stdout (self, data) :
    kwds = self.get_kwds()
    event = LogEvent(data, **kwds)
    wx.PostEvent(self.window, event)

  def callback_error (self, error, traceback_info) :
    kwds = self.get_kwds()
    event = JobExceptionEvent((error, traceback_info), **kwds)
    wx.PostEvent(self.window, event)

  def callback_abort (self) :
    kwds = self.get_kwds()
    event = JobKilledEvent(None, **kwds)
    wx.PostEvent(self.window, event)

  def callback_final (self, result) :
    kwds = self.get_kwds()
    event = JobCompleteEvent(result, **kwds)
    wx.PostEvent(self.window, event)

  def callback_other (self, data) :
    kwds = self.get_kwds()
    event = CallbackEvent(data, **kwds)
    wx.PostEvent(self.window, event)

# this just adds event posting callbacks to the original class
class process_with_gui_callbacks (thread_utils.process_with_callbacks) :
  def __init__ (self, proxy, target, args=(), kwargs={}, buffer_stdout=True) :
    thread_utils.process_with_callbacks.__init__(self,
      target = target,
      args=args,
      kwargs=kwargs,
      callback_stdout = proxy.callback_stdout,
      callback_final  = proxy.callback_final,
      callback_err    = proxy.callback_error,
      callback_abort  = proxy.callback_abort,
      callback_other  = proxy.callback_other,
      buffer_stdout   = buffer_stdout)

  def set_job (self, job) :
    pass

class simple_gui_process (process_with_gui_callbacks) :
  def __init__ (self, window, target, args=(), kwargs={}) :
    proxy = event_agent(window, None, None)
    process_with_gui_callbacks.__init__(self,
      proxy=proxy,
      target=target,
      args=args,
      kwargs=kwargs,
      buffer_stdout=True)

class ThreadProgressDialog (pyprogress.PyProgress) :
  def __init__ (self, parent, title, message) :
    pyprogress.PyProgress.__init__(self, parent, -1, title, message,
      agwStyle=wx.PD_ELAPSED_TIME|wx.PD_APP_MODAL)
    self.SetGaugeProportion(0.15)
    self.SetGaugeSteps(50)
    self.SetGaugeBackground(wx.Colour(235, 235, 235))
    self.SetFirstGradientColour(wx.Colour(235,235,235))
    self.SetSecondGradientColour(wx.Colour(120, 200, 255))

class download_file_basic (object) :
  def __init__ (self, window, dl_func, args) :
    assert isinstance(window, wx.EvtHandler)
    assert hasattr(dl_func, "__call__")
    assert (isinstance(args, list) or isinstance(args, tuple))
    self.window = window
    window.Connect(-1, -1, DOWNLOAD_COMPLETE_ID, self.OnComplete)
    self.dl_func = dl_func
    self.args = args
    self.t = threading.Thread(target=self.run)
    self.t.start()

  def run (self) :
    try :
      result = self.dl_func(self.args)
    except Exception, e :
      result = (None, str(e))
    finally :
      wx.PostEvent(self.window, DownloadCompleteEvent(result))
    return result

  def OnComplete (self, event) :
    if isinstance(event.data, str) :
      wx.MessageBox(message="File downloaded to %s" % event.data)
    else :
      wx.MessageBox(message="Error downloading file: %s" % event.data[1],
        caption="Download error", style=wx.ICON_ERROR)
    self.t.join()

def run_function_as_thread_in_dialog (parent, thread_function, title, message) :
  dlg = ThreadProgressDialog(None, title, message)
  t = thread_utils.simple_task_thread(thread_function, dlg)
  t.start()
  while True :
    if t.is_complete() or t.exception_raised() :
      #dlg.Destroy()
      dlg.Hide()
      break
    else :
      dlg.UpdatePulse()
    wx.MilliSleep(30)
  dlg.Destroy()
  wx.SafeYield()
  if t.exception_raised() :
    raise RuntimeError("An exception occurred while running this process: %s" %
      t.get_error())
  return t.return_value

# TODO
class ProcessDialog (wx.Dialog) :
  def __init__ (self, parent, message, caption, callback=None) :
    wx.Dialog.__init__(self,
      parent=parent,
      title=caption,
      style=wx.RAISED_BORDER|wx.CAPTION)
    self.callback = callback
    self.process = None
    self._error = None
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL, 5)
    msg_txt = wx.StaticText(self, -1, message)
    msg_txt.Wrap(400)
    szr2.Add(msg_txt, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)
    self.gauge = wx.Gauge(parent=self)
    self.gauge.SetRange(100)
    szr2.Add(self.gauge, 1, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5)
    abort_btn = wx.Button(parent=self,
      label="Abort")
    self.Bind(wx.EVT_BUTTON, self.OnAbort, abort_btn)
    szr2.Add(abort_btn, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5)
    szr.Fit(self)
    self.Centre(wx.BOTH)

  def run (self, process) :
    self.process = process
    self.process.start()
    self.gauge.Pulse()
    return self.ShowModal()

  def OnAbort (self, event) :
    self.process.abort()
    self.EndModal(wx.ID_CANCEL)

  def OnError (self, event) :
    self._error = event.data
    self.EndModal(wx.ID_CANCEL)

  def exception_raised (self) :
    return (self._error is not None)

  def handle_error (self) :
    if isinstance(self._error, Exception) :
      raise event.data
    elif isinstance(self._error, tuple) :
      exception, traceback = self._error
      if (isinstance(exception, Sorry)) :
        raise Sorry(str(exception))
      raise RuntimeError("""\
Error in subprocess!
 Original error: %s
 Original traceback:
%s""" % (str(exception), traceback))
    else :
      raise Sorry("error in child process: %s" % str(self._error))
   # finally :
   #   self.EndModal(wx.ID_CANCEL)

  def OnComplete (self, event) :
    try :
      if (self.callback is not None) :
        self.callback(event.data)
    finally :
      self._result = event.data
      self.EndModal(wx.ID_OK)

  def get_result (self) :
    return getattr(self, "_result", None)

def run_function_as_process_in_dialog (
    parent,
    thread_function,
    title,
    message,
    callback=None) :
  dlg = ProcessDialog(
    parent=parent,
    message=message,
    caption=title,
    callback=callback)
  setup_process_gui_events(
    window=dlg,
    OnExcept=dlg.OnError,
    OnComplete=dlg.OnComplete)
  cb = event_agent(dlg)
  p = thread_utils.process_with_callbacks(
    target=thread_function,
    callback_final=cb.callback_final,
    callback_err=cb.callback_error,
    buffer_stdout=True,
    sleep_after_start=1)
  result = None
  if (dlg.run(p) == wx.ID_OK) :
    result = dlg.get_result()
  elif dlg.exception_raised() :
    dlg.handle_error()
  wx.CallAfter(dlg.Destroy)
  return result

if (__name__ == "__main__") :
  from libtbx.test_utils import approx_equal, Exception_expected
  import math
  import sys
  def test_function_1 (*args, **kwds) :
    n = 0
    for i in range(25000) :
      x = math.sqrt(i)
      print x
      n += x
    return n
  def test_function_2 (*args, **kwds) :
    n = 0
    for i in range(100000) :
      x = math.sqrt(i)
      n += x
    return n
  def test_function_3 (*args, **kwds) :
    raise RuntimeError("This is a test!")
  def excepthook (*args, **kwds) :
    pass
  sys._excepthook = excepthook
  app = wx.App(0)
  result = run_function_as_process_in_dialog(
    parent=None,
    thread_function=test_function_1,
    title="Test subprocess",
    message="Running test function as separate process...",
    callback=None)
  if (result is not None) :
    assert approx_equal(result, 2635152.11891, eps=0.0001)
  #result2 = run_function_as_thread_in_dialog(
  #  parent=None,
  #  thread_function=test_function_2,
  #  title="Test subprocess",
  #  message="Running test function in Python thread...")
  #assert approx_equal(result2, 21081692.7462, eps=0.0001)
  try :
    result = run_function_as_process_in_dialog(
      parent=None,
      thread_function=test_function_3,
      title="Test subprocess",
      message="Running test function as separate process...",
      callback=None)
  except RuntimeError :
    pass
  else :
    raise Exception_expected
  wx.Yield()
  print "OK"
