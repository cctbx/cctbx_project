
from wx.lib.agw import pyprogress
import wx
import threading
import sys

class task_thread (threading.Thread) :
  def __init__ (self, thread_function, parent_window=None) :
    threading.Thread.__init__(self)
    self.f = thread_function
    self.return_value = None
    self._parent_window = parent_window
    self.is_complete = False
    self.exception_raised = None

  def Launch (self) :
    self.start()

  def run (self) :
    try :
      self.return_value = self.f()
    except Exception, e :
      sys.stderr.write(str(e))
      self.exception_raised = str(e)
    self.is_complete = True

class ThreadProgressDialog (pyprogress.PyProgress) :
  def __init__ (self, parent, title, message) :
    pyprogress.PyProgress.__init__(self, parent, -1, title, message,
      style=wx.PD_ELAPSED_TIME|wx.PD_APP_MODAL)
    self.SetGaugeProportion(0.15)
    self.SetGaugeSteps(50)
    self.SetGaugeBackground(wx.Colour(235, 235, 235))
    self.SetFirstGradientColour(wx.Colour(235,235,235))
    self.SetSecondGradientColour(wx.Colour(120, 200, 255))

def run_function_as_thread_in_dialog (thread_function, title, message) :
  dlg = ThreadProgressDialog(None, title, message)
  t = task_thread(thread_function, dlg)
  t.Launch()
  while True :
    if t.is_complete or t.exception_raised :
      #dlg.Destroy()
      dlg.Hide()
      break
    else :
      dlg.UpdatePulse()
    wx.MilliSleep(30)
  dlg.Destroy()
  wx.SafeYield()
  if t.exception_raised :
    raise RuntimeError("An exception occurred while running this process: %s")
  return t.return_value
