
import wxtbx.xray_viewer
import wx
import os

class XrayView (wx.ScrolledWindow) :
  def __init__ (self, *args, **kwds) :
    self._img = None
    super(XrayView, self).__init__(*args, **kwds)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.Bind(wx.EVT_SIZE, self.OnSize)
    self.Bind(wx.EVT_MOTION, self.OnMove)

  def DoGetBestSize (self) :
    if (self._img is None) :
      return (640, 640)
    else :
      x, y = self._img.get_size()
      return (x+20, y+20)

  def OnPaint (self, event) :
    dc = wx.AutoBufferedPaintDCFactory(self)
    x_offset, y_offset = self.GetViewStart()
    bitmap = self._img.get_bitmap()
    dc.DrawBitmap(bitmap, 10 - x_offset, 10 - y_offset)
    center_x, center_y = self._img.get_beam_center()
    dc.SetPen(wx.Pen('red'))
    x0 = center_x - x_offset + 10
    y0 = center_y - y_offset + 10
    dc.DrawLine(x0 - 10, y0, x0 + 10, y0)
    dc.DrawLine(x0, y0 - 10, x0, y0 + 10)

  def OnSize (self, event) :
    w, h = self.DoGetBestSize()
    self.SetVirtualSize((w,h))
    self.SetScrollbars(1, 1, w, h)

  def OnMove (self, event) :
    x, y = event.GetPositionTuple()
    x_offset, y_offset = self.GetViewStart()
    img_x = x + x_offset - 10
    img_y = y + y_offset - 10
    img_w, img_h = self._img.get_size()
    if (img_x < 0) or (img_x > img_w) or (img_y < 0) or (img_y > img_h) :
      self.GetParent().display_resolution(None)
    else :
      center_x, center_y = self._img.get_beam_center()
      d_min = self._img.get_d_min_at_point(img_x, img_y)
      self.GetParent().display_resolution(d_min)

  def load_image (self, file_name) :
    wx.GetApp().Yield(True)
    self._img = wxtbx.xray_viewer.image(file_name)
    self.OnSize(None)
    self.Refresh()

class XrayFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    super(XrayFrame, self).__init__(*args, **kwds)
    self._panel = XrayView(self, -1, size=(640,640))
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(self.sizer)
    self.sizer.Add(self._panel, 1, wx.EXPAND)
    self.statusbar = self.CreateStatusBar()
    self.Fit()

  def load_image (self, file_name) :
    file_name = os.path.abspath(file_name)
    self._panel.load_image(file_name)
    self.SetTitle(file_name)
    self.display_resolution(None)
    self.Layout()

  def display_resolution (self, d_min) :
    if (d_min is None) :
      self.statusbar.SetStatusText("Move mouse over image to show resolution")
    else :
      self.statusbar.SetStatusText("d_min = %.2f A" % d_min)
