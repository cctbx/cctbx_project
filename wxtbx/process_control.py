
# TODO more comprehensive tests

from __future__ import absolute_import, division, print_function
from wx.lib.agw import pyprogress
import wx
from libtbx import thread_utils
from libtbx import runtime_utils
from libtbx import easy_pickle
from libtbx import easy_run
from libtbx.utils import Sorry, Abort, download_progress
import threading
import random
import locale
import math
import os
from six.moves import range

JOB_START_ID = wx.NewId()
LOG_UPDATE_ID = wx.NewId()
CALLBACK_ID = wx.NewId()
JOB_EXCEPTION_ID = wx.NewId()
JOB_KILLED_ID = wx.NewId()
JOB_COMPLETE_ID = wx.NewId()
JOB_PAUSE_ID = wx.NewId()
JOB_RESUME_ID = wx.NewId()
DOWNLOAD_COMPLETE_ID = wx.NewId()
DOWNLOAD_INCREMENT_ID = wx.NewId()

class SubprocessEvent(wx.PyEvent):
  event_id = None

  def __init__(self, data, **kwds):
    super(SubprocessEvent, self).__init__()
    self.data = data
    if hasattr(self, '_getAttrDict'):
      self._getAttrDict().update(kwds)
    else:
      self.__dict__.update(kwds)
    self.SetEventType(self.event_id)

class JobStartEvent(SubprocessEvent):
  event_id = JOB_START_ID

  def __init__(self, data, **kwds):
    super(JobStartEvent, self).__init__(data, **kwds)

class LogEvent(SubprocessEvent):
  event_id = LOG_UPDATE_ID

  def __init__(self, data, **kwds):
    super(LogEvent, self).__init__(data, **kwds)


class JobExceptionEvent(SubprocessEvent):
  event_id = JOB_EXCEPTION_ID

  def __init__(self, data, **kwds):
    super(JobExceptionEvent, self).__init__(data, **kwds)

class JobKilledEvent(SubprocessEvent):
  event_id = JOB_KILLED_ID

  def __init__(self, data, **kwds):
    super(JobKilledEvent, self).__init__(data, **kwds)

class JobCompleteEvent(SubprocessEvent):
  event_id = JOB_COMPLETE_ID

  def __init__(self, data, **kwds):
    super(JobCompleteEvent, self).__init__(data, **kwds)

class CallbackEvent(SubprocessEvent):
  event_id = CALLBACK_ID

  def __init__(self, data, **kwds):
    super(CallbackEvent, self).__init__(data, **kwds)

class JobPauseEvent(SubprocessEvent):
  event_id = JOB_PAUSE_ID

  def __init__(self, data, **kwds):
    super(JobPauseEvent, self).__init__(data, **kwds)

class JobResumeEvent(SubprocessEvent):
  event_id = JOB_RESUME_ID

  def __init__(self, data, **kwds):
    super(JobResumeEvent, self).__init__(data, **kwds)

class DownloadCompleteEvent(SubprocessEvent):
  event_id = DOWNLOAD_COMPLETE_ID

  def __init__(self, data, **kwds):
    super(DownloadCompleteEvent, self).__init__(data, **kwds)

class DownloadIncrementEvent(SubprocessEvent):
  event_id = DOWNLOAD_INCREMENT_ID

  def __init__(self, data, **kwds):
    super(DownloadIncrementEvent, self).__init__(data, **kwds)

def setup_stdout_logging_event(window, OnPrint):
  window.Connect(-1, -1, LOG_UPDATE_ID, OnPrint)

def setup_process_gui_events(
    window,
    OnStart=None,
    OnPrint=None,
    OnUpdate=None,
    OnExcept=None,
    OnAbort=None,
    OnComplete=None,
    OnPause=None,
    OnResume=None):
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
  if OnPause is not None :
    assert hasattr(OnPause, "__call__")
    window.Connect(-1, -1, JOB_PAUSE_ID, OnPause)
  if OnResume is not None :
    assert hasattr(OnResume, "__call__")
    window.Connect(-1, -1, JOB_RESUME_ID, OnResume)

class event_agent(object):
  def __init__(self, window, **kwds):
    self.window = window
    self._kwds = dict(kwds)
    self.__dict__.update(kwds)

  def get_kwds(self):
    return self._kwds

  def callback_start(self, data):
    kwds = self.get_kwds()
    event = JobStartEvent(data, **kwds)
    wx.PostEvent(self.window, event)

  def callback_stdout(self, data):
    kwds = self.get_kwds()
    event = LogEvent(data, **kwds)
    wx.PostEvent(self.window, event)

  def callback_error(self, error, traceback_info):
    kwds = self.get_kwds()
    event = JobExceptionEvent((error, traceback_info), **kwds)
    wx.PostEvent(self.window, event)

  def callback_abort(self):
    kwds = self.get_kwds()
    event = JobKilledEvent(None, **kwds)
    wx.PostEvent(self.window, event)

  def callback_final(self, result):
    kwds = self.get_kwds()
    event = JobCompleteEvent(result, **kwds)
    wx.PostEvent(self.window, event)

  def callback_other(self, data):
    kwds = self.get_kwds()
    event = CallbackEvent(data, **kwds)
    wx.PostEvent(self.window, event)

  def callback_pause(self):
    kwds = self.get_kwds()
    event = JobPauseEvent(None, **kwds)
    wx.PostEvent(self.window, event)

  def callback_resume(self):
    kwds = self.get_kwds()
    event = JobResumeEvent(None, **kwds)
    wx.PostEvent(self.window, event)

# simplified for when the window is really the app object
class background_event_agent(event_agent):
  def callback_stdout(self, data):
    pass

  def callback_other(self, data):
    pass

class detached_process(runtime_utils.detached_process_client):
  def __init__(self, params, proxy):
    runtime_utils.detached_process_client.__init__(self, params)
    self.proxy = proxy

  def callback_start(self, data):
    self.proxy.callback_start(data)

  def callback_stdout(self, data):
    self.proxy.callback_stdout(data)

  def callback_other(self, data):
    self.proxy.callback_other(data)

  def callback_abort(self):
    self.proxy.callback_abort()

  def callback_final(self, result):
    self.proxy.callback_final(result)

  def callback_error(self, error, traceback_info):
    self.proxy.callback_error(error, traceback_info)

  def callback_pause(self):
    self.proxy.callback_pause()

  def callback_resume(self):
    self.proxy.callback_resume()

  def start(self):
    pass

# this just adds event posting callbacks to the original class
class process_with_gui_callbacks(thread_utils.process_with_callbacks):
  def __init__(self, proxy, target, args=(), kwargs={}, buffer_stdout=True):
    thread_utils.process_with_callbacks.__init__(self,
      target = target,
      args=args,
      kwargs=kwargs,
      callback_stdout = proxy.callback_stdout,
      callback_final  = proxy.callback_final,
      callback_err    = proxy.callback_error,
      callback_abort  = proxy.callback_abort,
      callback_other  = proxy.callback_other,
      callback_pause  = proxy.callback_pause,
      callback_resume = proxy.callback_resume,
      buffer_stdout   = buffer_stdout)

  def set_job(self, job):
    pass

  def purge_files(self):
    pass

class simple_gui_process(process_with_gui_callbacks):
  def __init__(self, window, target, args=(), kwargs={}):
    # XXX fix for phenix gui - is this necessary?
    proxy = event_agent(window, project_id=None, job_id=None)
    process_with_gui_callbacks.__init__(self,
      proxy=proxy,
      target=target,
      args=args,
      kwargs=kwargs,
      buffer_stdout=True)

class ThreadProgressDialog(pyprogress.PyProgress):
  def __init__(self, parent, title, message):
    pyprogress.PyProgress.__init__(self, parent, -1, title, message,
      agwStyle=wx.PD_ELAPSED_TIME|wx.PD_APP_MODAL)
    self.SetGaugeProportion(0.15)
    self.SetGaugeSteps(50)
    self.SetGaugeBackground(wx.Colour(235, 235, 235))
    self.SetFirstGradientColour(wx.Colour(235,235,235))
    self.SetSecondGradientColour(wx.Colour(120, 200, 255))

class download_file_basic(object):
  def __init__(self, window, dl_func, args):
    assert isinstance(window, wx.EvtHandler)
    assert hasattr(dl_func, "__call__")
    assert (isinstance(args, list) or isinstance(args, tuple))
    self.window = window
    window.Connect(-1, -1, DOWNLOAD_COMPLETE_ID, self.OnComplete)
    self.dl_func = dl_func
    self.args = args
    self.t = threading.Thread(target=self.run)
    self.t.start()

  def run(self):
    try :
      result = self.dl_func(self.args)
    except Exception as e :
      result = (None, str(e))
    finally :
      wx.PostEvent(self.window, DownloadCompleteEvent(result))
    return result

  def OnComplete(self, event):
    from six import string_types
    if isinstance(event.data, string_types):
      wx.MessageBox(message="File downloaded to %s" % event.data)
    else :
      wx.MessageBox(message="Error downloading file: %s" % event.data[1],
        caption="Download error", style=wx.ICON_ERROR)
    self.t.join()

class DownloadProgressDialog(wx.ProgressDialog, download_progress):
  """
  Dialog for displaying download progress.  The actual download (not
  implemented here) should be run in a separate thread, with a reasonable
  chunk size, and call download_progress.increment() as each new chunk is
  downloaded.
  """
  def __init__(self, parent, title, message):
    download_progress.__init__(self)
    wx.ProgressDialog.__init__(self, parent=parent,
      title=title,
      message=message,
      style=wx.PD_ELAPSED_TIME|wx.PD_CAN_ABORT|wx.PD_AUTO_HIDE,
      maximum=100)
    self.Connect(-1, -1, DOWNLOAD_INCREMENT_ID, self.OnIncrement)
    self.Connect(-1, -1, DOWNLOAD_COMPLETE_ID, self.OnComplete)
    self._continue = True

  def show_progress(self):
    if (not self._continue):
      return False
    locale.setlocale(locale.LC_ALL, 'en_US')
    pct = self.percent_finished()
    msg = "%s/%s KB downloaded" % (
      locale.format("%d", self.n_kb_elapsed, grouping=True),
      locale.format("%d", self.n_kb_total, grouping=True))
    evt = DownloadIncrementEvent(data=(pct, msg))
    wx.PostEvent(self, evt)
    return self._continue

  def OnIncrement(self, event):
    (cont, skip) = self.Update(value=event.data[0], newmsg=event.data[1])
    self._continue = cont

  def OnComplete(self, event):
    self.Hide()
    self.Close()
    # FIXME destroying the dialog crashes wxPython 2.9.5/osx-coocoa

  def complete(self):
    evt = DownloadCompleteEvent(data=None)
    wx.PostEvent(self, evt)

class BackgroundDownloadDialog(pyprogress.PyProgress, download_progress):
  """
  Placeholder for downloads which block the child thread; will pulse
  continuously but not show changing status.
  """
  def __init__(self, parent, title, message):
    download_progress.__init__(self)
    pyprogress.PyProgress.__init__(self, parent, -1, title, message,
      agwStyle=wx.PD_ELAPSED_TIME|wx.PD_CAN_ABORT|wx.PD_AUTO_HIDE)
    self.SetGaugeProportion(0.15)
    self.SetGaugeSteps(100)
    self.SetGaugeBackground(wx.Colour(235, 235, 235))
    self.SetFirstGradientColour(wx.Colour(235,235,235))
    self.SetSecondGradientColour(wx.Colour(120, 200, 255))
    self.Connect(-1, -1, DOWNLOAD_COMPLETE_ID, self.OnComplete)
    self._continue = True

  def show_progress(self):
    if (not self._continue):
      return False
    return self._continue

  def OnComplete(self, event):
    self.Hide()
    self.Close()

  def complete(self):
    evt = DownloadCompleteEvent(data=None)
    wx.PostEvent(self, evt)

def run_function_as_thread_in_dialog(parent, thread_function, title, message):
  dlg = ThreadProgressDialog(None, title, message)
  t = thread_utils.simple_task_thread(thread_function, dlg)
  t.start()
  while True :
    if t.is_complete() or t.exception_raised():
      #dlg.Destroy()
      dlg.Hide()
      break
    else :
      dlg.UpdatePulse()
    wx.MilliSleep(30)
  dlg.Destroy()
  wx.SafeYield()
  if t.exception_raised():
    raise RuntimeError("An exception occurred while running this process: %s" %
      t.get_error())
  return t.return_value

# TODO
class ProcessDialog(wx.Dialog):
  def __init__(self, parent, message, caption, callback=None):
    wx.Dialog.__init__(self,
      parent=parent,
      title=caption,
      style=wx.RAISED_BORDER|wx.CAPTION)
    self.callback = callback
    self.process = None
    self._error = None
    self._aborted = False
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL, 5)
    msg_txt = wx.StaticText(self, -1, message)
    msg_txt.Wrap(400)
    szr2.Add(msg_txt, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALL, 5)
    self.gauge = wx.Gauge(parent=self, size=(300,-1))
    self.gauge.SetRange(100)

    # XXX  TEMPORARY PYTHON 3 FIX TT
    #szr2.Add(self.gauge, 1, wx.ALL|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL, 5)
    szr2.Add(self.gauge, 1, wx.ALL, 5)
    abort_btn = wx.Button(parent=self,
      label="Abort")
    self.Bind(wx.EVT_BUTTON, self.OnAbort, abort_btn)
    # XXX  TEMPORARY PYTHON 3 FIX TT
    #szr2.Add(abort_btn, 0, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5)
    szr2.Add(abort_btn, 0, wx.ALL, 5)
    self.SetMinSize((300,100))
    szr.Fit(self)
    self.Centre(wx.BOTH)

  def run(self, process):
    self.process = process
    self._timer = wx.Timer(owner=self)
    self.Bind(wx.EVT_TIMER, self.OnTimer)
    self._timer.Start(100)
    self.process.start()
    self.gauge.Pulse()
    return self.ShowModal()

  def OnTimer(self, event):
    if hasattr(self.process,'update'):
      self.process.update()
    self.gauge.Pulse()

  def OnAbort(self, event):
    self.process.abort()
    self._aborted = True
    self._timer.Stop()
    try:
      self.EndModal(wx.ID_CANCEL)
    except Exception as e:
      pass # C++ was deleted

  def OnError(self, event):
    self._error = event.data
    self._timer.Stop()
    try:
      self.EndModal(wx.ID_CANCEL)
    except Exception as e:
      pass # C++ was deleted

  def exception_raised(self):
    return (self._error is not None)

  def was_aborted(self):
    return (self._aborted)

  def handle_error(self):
    if isinstance(self._error, Exception):
      raise event.data
    elif isinstance(self._error, tuple):
      exception, traceback = self._error
      if (isinstance(exception, Sorry)):
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

  def OnComplete(self, event):
    try :
      if (self.callback is not None):
        self.callback(event.data)
    finally :
      self._result = event.data
      self._timer.Stop()
      self.EndModal(wx.ID_OK)

  def get_result(self):
    return getattr(self, "_result", None)

def run_function_as_process_in_dialog(
    parent,
    thread_function,
    title,
    message,
    callback=None,
    project_id=None,
    job_id=None):
  dlg = ProcessDialog(
    parent=parent,
    message=message,
    caption=title,
    callback=callback)
  setup_process_gui_events(
    window=dlg,
    OnExcept=dlg.OnError,
    OnComplete=dlg.OnComplete)
  cb = event_agent(dlg, project_id=project_id, job_id=job_id)
  p = thread_utils.process_with_callbacks(
    target=thread_function,
    callback_final=cb.callback_final,
    callback_err=cb.callback_error,
    buffer_stdout=True,
    sleep_after_start=1)
  result = None
  abort = False
  if (dlg.run(p) == wx.ID_OK):
    result = dlg.get_result()
  elif dlg.exception_raised():
    dlg.handle_error()
  elif (dlg.was_aborted()):
    abort = True
  wx.CallAfter(dlg.Destroy)
  if (abort):
    raise Abort()
  return result

# TODO this is awful, needs to be re-thought
def run_function_as_detached_process_in_dialog(
    parent,
    thread_function,
    title,
    message,
    tmp_dir,
    callback=None,
    project_id=None,
    job_id=None):
  if (tmp_dir is None):
    tmp_dir = os.getcwd()
  params = runtime_utils.process_master_phil.extract()
  params.tmp_dir = tmp_dir
  if (job_id is None):
    job_id = str(os.getpid()) + "_" + str(int(random.random() * 1000))
  params.prefix = str(job_id)
  target = runtime_utils.detached_process_driver(target=thread_function)
  run_file = os.path.join(tmp_dir, "libtbx_run_%s.pkl" % job_id)
  easy_pickle.dump(run_file, target)
  params.run_file = run_file
  eff_file = os.path.join(tmp_dir, "libtbx_run_%s.eff" % job_id)
  runtime_utils.write_params(params, eff_file)
  dlg = ProcessDialog(
    parent=parent,
    message=message,
    caption=title,
    callback=callback)
  setup_process_gui_events(
    window=dlg,
    OnExcept=dlg.OnError,
    OnAbort=dlg.OnAbort,
    OnComplete=dlg.OnComplete)
  agent = event_agent(
    window=dlg,
    project_id=project_id,
    job_id=job_id)
  process = detached_process(params, proxy=agent)
  cb = event_agent(dlg, project_id=project_id, job_id=job_id)
  easy_run.call("libtbx.start_process \"%s\" &" % eff_file)
  result = None
  abort = False
  if (dlg.run(process) == wx.ID_OK):
    result = dlg.get_result()
  elif dlg.exception_raised():
    dlg.handle_error()
  elif (dlg.was_aborted()):
    abort = True
  wx.CallAfter(dlg.Destroy)
  if (abort):
    raise Abort()
  return result

########################################################################
# XXX regression testing utilities
def test_function_1(*args, **kwds):
  n = 0
  for i in range(25000):
    x = math.sqrt(i)
    print(x)
    n += x
  return n
def test_function_2(*args, **kwds):
  n = 0
  for i in range(100000):
    x = math.sqrt(i)
    n += x
  return n
def test_function_3(*args, **kwds):
  raise RuntimeError("This is a test!")
