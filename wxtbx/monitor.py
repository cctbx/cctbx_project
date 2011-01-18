
import wxtbx.bitmaps
from wxtbx import metallicbutton
from wx.lib.agw import flatnotebook as fnb
import wx

class StatusPanel (wx.Panel) :
  def __init__ (self, *args, **kwds) :
    wx.Panel.__init__(self, *args, **kwds)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.content_sizer = wx.BoxSizer(wx.VERTICAL)
    self._errors = []
    self._btns = {}
    self._current_status = None
    self._status_warn = False
    self._std_font = wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL,wx.NORMAL)
    self._warn_font = wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.BOLD)
    self._btn_font = wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.NORMAL,wx.NORMAL)
    self.btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.sizer.Add(self.btn_sizer)
    s2 = wx.BoxSizer(wx.HORIZONTAL)
    self.sizer.Add(s2)
    status_label = wx.StaticText(self, -1, "Status:")
    status_label.SetFont(self._warn_font)
    s2.Add(status_label, 0, wx.TOP|wx.LEFT|wx.RIGHT|wx.ALIGN_LEFT, 5)
    self.status = wx.StaticText(self, -1, "N/A")
    self.status.SetFont(self._std_font)
    s2.Add(self.status, 0, wx.TOP|wx.LEFT|wx.RIGHT|wx.ALIGN_LEFT, 5)
    self.sizer.Add(self.content_sizer)

  def AddActionButton (self, label, icon=None) :
    bmp = None
    if (icon is not None) :
      bmp = wxtbx.bitmaps.fetch_icon_bitmap("actions", icon, 16)
    btn = metallicbutton.MetallicButton(
      parent=self,
      label=label,
      bmp=bmp,
      highlight_color=(200,220,240))
    btn.SetFont(self._std_font)
    if (len(self._btns) == 0) :
      self.btn_sizer.Add((5,1))
    self.btn_sizer.Add(btn, 0, wx.TOP|wx.RIGHT|wx.ALIGN_LEFT, 5)
    self._btns[label] = btn
    return btn

  def GetButton (self, label) :
    return self._btns.get(label, None)

  def CreateErrorBox (self) :
    error_txt = wx.StaticText(self, -1, "Errors:")
    error_txt.SetFont(wx.Font(9, wx.FONTFAMILY_DEFAULT, wx.NORMAL, wx.BOLD))
    self.sizer.Add(error_txt, 0, wx.TOP|wx.LEFT|wx.RIGHT|wx.ALIGN_LEFT, 5)
    self.errors = wx.TextCtrl(self, -1, size=(400,120),
      style=wx.TE_MULTILINE|wx.TE_READONLY|wx.NO_BORDER)
    self.errors.SetFont(wx.Font(9, wx.MODERN, wx.NORMAL, wx.NORMAL))
    self.sizer.Add(self.errors, 1, wx.ALL|wx.EXPAND|wx.ALIGN_LEFT, 5)

  def UpdateErrors (self, errors) :
    if (len(self._errors) != len(errors)) :
      self.errors.SetValue("\n".join(errors))

  def UpdateStatus (self, status, warn=False) :
    if (self._current_status != status) :
      self.status.SetLabel(status)
      if (warn != self._status_warn) :
        self._status_warn = warn
        if warn :
          self.status.SetFont(self._warn_font)
          self.status.SetForegroundColour((200,0,0))
        else :
          self.status.SetFont(self._std_font)
          self.status.SetForegroundColour((0,0,0))
      self.sizer.Layout()

class MonitorWindow (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    wx.MiniFrame.__init__(self, *args, **kwds)
    #self.panel = wx.Panel(self)
    self.nb = fnb.FlatNotebook(
      parent=self,
      agwStyle=fnb.FNB_FF2|fnb.FNB_NO_NAV_BUTTONS|fnb.FNB_NO_X_BUTTON|fnb.FNB_NODRAG)
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.sizer.Add(self.nb, 1)
    self.SetSizer(self.sizer)
    #self.panel_sizer = wx.BoxSizer(wx.VERTICAL)
    #self.panel.SetSizer(self.panel_sizer)
    self._panels = {}
    self._errors = {}
    self.timer = None
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def AddPanel (self, title, panel_name) :
    p = StatusPanel(self.nb, -1, style=wx.NO_BORDER)
    p.CreateErrorBox()
    self._panels[panel_name] = p
    self.nb.AddPage(p, title)

  def GetPanel (self, panel_name) :
    return self._panels[panel_name]

  def UpdateStatus (self, panel_name, status, errors, warn=False) :
    p = self._panels[panel_name]
    p.UpdateErrors(errors)
    p.UpdateStatus(status, warn=warn)
    self.Refresh()

  def OnClose (self, event) :
    if (self.timer is not None) :
      self.timer.Stop()
    self.Destroy()

  def OnDestroy (self, event) : # XXX override this in subclasses
    pass

  def OnTimer (self, event) :
    pass

  def StartTimer (self, wait=30) :
    self.timer = wx.Timer(owner=self)
    self.Bind(wx.EVT_TIMER, self.OnTimer)
    self.timer.Start(wait * 1000)

if (__name__ == "__main__") :
  app = wx.App(0)
  f = MonitorWindow(None, -1, "Test monitor window",
    style=wx.CAPTION|wx.CLOSE_BOX)
  f.AddPanel("Coot", "coot")
  f.GetPanel("coot").AddActionButton("Quit", "stop")
  #f.GetPanel("coot").
  f.AddPanel("PyMOL", "pymol")
  f.Fit()
  f.Show()
  f.UpdateStatus("coot", status="7 socket timeouts",
    errors=["12:14:01 :: connection refused",
            "12:14:04 :: broken pipe",
            "12:14:05 :: broken pipe",],
    warn=True)
  f.UpdateStatus("pymol", status="running",
    errors=[],
    warn=False)
  app.MainLoop()
