
import wxtbx.bitmaps
_import_error = None
try :
  from wx import webkit
except ImportError, e :
  webkit = None
  _import_error = str(e)
import wx

class webkit_frame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.SetupToolbar()
    self._was_shown = False
    if (webkit is None) :
      raise ImportError(("Sorry, this module is only supported on Macintosh "+
        "at present.  (Original error: %s)") % _import_error)
    self.viewer = webkit.WebKitCtrl(self, -1)
    szr.Add(self.viewer, 1, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.SetInitialSize((800,640))
    #self.Bind(wx.EVT_WINDOW_CREATE, self.OnShow)
    self.Bind(webkit.EVT_WEBKIT_STATE_CHANGED, self.OnChangeState, self.viewer)
    self.Bind(webkit.EVT_WEBKIT_BEFORE_LOAD, self.OnBeforeLoad, self.viewer)

  def SetHomepage (self, url) :
    self.home_url = url

  def SetupToolbar (self) :
    if (wxtbx.bitmaps.icon_lib is None) :
      return
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    commands = [
      ("filesystems", "folder_home", "OnHome", "Home"),
      ("actions", "back", "OnBack", "Back"),
      ("actions", "forward", "OnForward", "Forward"),
      ("actions", "reload", "OnReload", "Reload"),
      ("actions", "stop", "OnStop", "Stop"),
    ]
    for (icon_class, icon_name, fname, label) in commands :
      bmp = wxtbx.bitmaps.fetch_icon_bitmap(icon_class, icon_name, 32)
      tool_button = self.toolbar.AddLabelTool(-1, label, bmp,
        shortHelp=label, kind=wx.ITEM_NORMAL)
      self.Bind(wx.EVT_MENU, getattr(self, fname), tool_button)
    self.toolbar.AddSeparator()
    phenix_bmp = wxtbx.bitmaps.fetch_custom_icon_bitmap("phenix.refine")
    phenix_btn = self.toolbar.AddLabelTool(-1, "PHENIX homepage", phenix_bmp,
      shortHelp="PHENIX homepage", kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnPhenixWeb, phenix_btn)
    self.toolbar.Realize()

  def OnShow (self, event) :
    if (not self._was_shown) :
      self.viewer.LoadURL(self.home_url)
      self._was_shown = True

  def OnHome (self, event) :
    self.viewer.LoadURL(self.home_url)

  def OnBack (self, event) :
    if self.viewer.CanGoBack() :
      self.viewer.GoBack()

  def OnForward (self, event) :
    if self.viewer.CanGoForward() :
      self.viewer.GoForward()

  def OnReload (self, event) :
    self.viewer.Reload()

  def OnStop (self, event) :
    self.viewer.Stop()

  def Open (self, url) :
    self.viewer.LoadURL(url)

  def OnPhenixWeb (self, event) :
    self.viewer.LoadURL("http://www.phenix-online.org")

  def OnChangeState (self, event) :
    state = event.GetState()
    url = event.GetURL()
    if (state == wx.webkit.WEBKIT_STATE_START) :
      self.statusbar.SetStatusText("Opening %s" % url)
    elif (state == wx.webkit.WEBKIT_STATE_TRANSFERRING) :
      self.statusbar.SetStatusText("Loading %s" % url)
    elif (state == wx.webkit.WEBKIT_STATE_STOP) :
      self.statusbar.SetStatusText("Viewing %s" % url)
    elif (state == wx.webkit.WEBKIT_STATE_FAILED) :
      self.statusbar.SetStatusText("Failed loading %s" % url)
    else :
      self.statusbar.SetStatusText("")

  def OnBeforeLoad (self, event) :
    pass
    #print event.GetNavigationType()

if __name__ == "__main__" :
  app = wx.App(0)
  frame = webkit_frame(None, -1, "wxtbx.browser example")
  #  size=(800,600))
  frame.SetHomepage("http://cci.lbl.gov")
  frame.Show()
  app.MainLoop()
