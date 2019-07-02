
from __future__ import absolute_import, division, print_function
from gltbx import wx_viewer
from gltbx.gl import *
from gltbx.glu import *
import gltbx.util
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
from math import log
import wx
from six.moves import range

def glbox(xyz1, xyz2):
  x1 = min(xyz1[0], xyz2[0])
  x2 = max(xyz1[0], xyz2[0])
  y1 = min(xyz1[1], xyz2[1])
  y2 = max(xyz1[1], xyz2[1])
  z1 = min(xyz1[2], xyz2[2])
  z2 = max(xyz1[2], xyz2[2])
  # top face (red)
  glBegin(GL_QUADS)
  glNormal3f(0, 0, 1)
  glVertex3f(x1, y1, z2)
  glVertex3f(x2, y1, z2)
  glVertex3f(x2, y2, z2)
  glVertex3f(x1, y2, z2)
  glEnd()
  # bottom face (blue)
  glBegin(GL_QUADS)
  glNormal3f(0, 0, -1)
  glVertex3f(x1, y1, z1)
  glVertex3f(x1, y2, z1)
  glVertex3f(x2, y2, z1)
  glVertex3f(x2, y1, z1)
  glEnd()
  # top side (y = y)
  glBegin(GL_QUADS)
  glNormal3f(0, -1, 0)
  glVertex3f(x1, y1, z1)
  glVertex3f(x2, y1, z1)
  glVertex3f(x2, y1, z2)
  glVertex3f(x1, y1, z2)
  glEnd()
  # bottom side (y = y + 1)
  glBegin(GL_QUADS)
  glNormal3f(0, 1, 0)
  glVertex3f(x1, y2, z2)
  glVertex3f(x2, y2, z2)
  glVertex3f(x2, y2, z1)
  glVertex3f(x1, y2, z1)
  glEnd()
  # left side (x' = x)
  glBegin(GL_QUADS)
  glNormal3f(-1, 0, 0)
  glVertex3f(x1, y2, z2)
  glVertex3f(x1, y2, z1)
  glVertex3f(x1, y1, z1)
  glVertex3f(x1, y1, z2)
  glEnd()
  # right side (x' = x+1)
  glBegin(GL_QUADS)
  glNormal3f(1, 0, 0)
  glVertex3f(x2, y1, z1)
  glVertex3f(x2, y2, z1)
  glVertex3f(x2, y2, z2)
  glVertex3f(x2, y1, z2)
  glEnd()

class pixel_viewer_3d(wx_viewer.wxGLWindow):
  def __init__(self, *args, **kwds):
    wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    self.SetSize((400,400))
    self.SetMinSize((400,400))
    self.flag_draw_boxes = False
    self.flag_log_scale = False
    self._img = None
    self.x_center = None
    self.y_center = None
    self.buffer_factor = 2.0
    self.scale_factor = 100
    self.rotation_center = (20,-20,0)
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([[0,0,0],[40,-40,40],[40,0,0],[0,-40,40]]))
    self.zoom_level = 16

  def InitGL(self):
    gltbx.util.handle_error()
    b = self.background_rgb
    glClearColor(b[0], b[1], b[2], 0.0)
    glDisable(GL_LIGHT0)
    glDisable(GL_LIGHTING)
    glDisable(GL_BLEND)
    glEnable(GL_LINE_SMOOTH)
    glEnable(GL_CULL_FACE)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_POLYGON_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
    self.initialize_modelview()
    gltbx.util.handle_error()

  def process_pick_points(self):
    pass

  def set_image(self, img):
    self._img = img

  def recenter(self, x, y):
    self.x_center = x
    self.y_center = y
    self.Refresh()

  def DrawGL(self):
    if (None in [self._img, self.x_center, self.y_center]) : return
    if self.flag_draw_boxes :
      r, g, b = self._img.get_opengl_background()
      glClearColor(r, g, b, 0.)
      self._draw_boxes()
    else :
      glClearColor(0, 0, 0, 0.)
      self._draw_lines()

  def _draw_lines(self):
    values = self._img.get_intensities_in_box(
      x=self.x_center,
      y=self.y_center,
      boxsize=40,
      mag=1)#self.zoom_level)
    scale_factor = self.scale_factor
    glColor3f(1, 1, 1)
    w_range = [0.0,0.0]
    glGetFloatv(GL_LINE_WIDTH_RANGE, w_range)
    line_width = 0.1
    if (w_range[0] > 0.1):
      line_width = w_range[0]
    glLineWidth(line_width)
    for i, row in enumerate(values):
      y = -i
      glBegin(GL_LINES)
      for j, value in enumerate(row):
        if (j > 0):
          x = j
          val1 = row[j-1] / scale_factor
          val2 = value / scale_factor
          if (self.flag_log_scale):
            val1 = log(val1)
            val2 = log(val2)
          glVertex3f(x-1, y, val1)
          glVertex3f(x, y, val2)
      glEnd()
    for i in range(40):
      x = i
      glBegin(GL_LINES)
      for j, row in enumerate(values):
        if (j > 0):
          y = -j
          val1 = values[j-1][i] / scale_factor
          val2 = row[i] / scale_factor
          if (self.flag_log_scale):
            val1 = log(val1)
            val2 = log(val2)
          glVertex3f(x, y+1, val1)
          glVertex3f(x, y, val2)
      glEnd()

  def _draw_boxes(self):
    values = self._img.get_intensities_in_box(
      x=self.x_center,
      y=self.y_center,
      boxsize=40,
      mag=1)#self.zoom_level)
    wx_image = self._img.get_zoomed_bitmap(
      x=self.x_center,
      y=self.y_center,
      boxsize=40,
      mag=1)
    if (wx_image is None) : return
    glPolygonMode(GL_FRONT, GL_FILL)
    scale_factor = self.scale_factor
    for y_, row in enumerate(values):
      y = -y_
      for x, value in enumerate(row):
        z = value / scale_factor
        if (self.flag_log_scale):
          z = log(z)
        R = wx_image.GetRed(x,y_)
        G = wx_image.GetGreen(x,y_)
        B = wx_image.GetBlue(x,y_)
        glColor3f(R/255., G/255., B/255.)
        glbox((x,y,0),(x+1,y+1,z))

class pixel_viewer_3d_frame(wx.MiniFrame):
  def __init__(self, *args, **kwds):
    wx.MiniFrame.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    ctrls = wx.Panel(self, style=wx.RAISED_BORDER)
    szr.Add(ctrls, 0, wx.EXPAND)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    ctrls.SetSizer(szr2)
    box1 = wx.CheckBox(ctrls, label="Wireframe rendering")
    box1.SetValue(True)
    self.Bind(wx.EVT_CHECKBOX, self.OnChangeRenderStyle, box1)
    szr2.Add(box1, 0, wx.ALL, 5)
    box2 = wx.CheckBox(ctrls, label="Log scaling")
    self.Bind(wx.EVT_CHECKBOX, self.OnChangeScale, box2)
    szr2.Add(box2, 0, wx.ALL, 5)
    szr2.Add(wx.StaticText(ctrls, label="Scale:"), 0, wx.ALL, 5)
    slider = wx.Slider(ctrls, size=(120,-1), style=wx.SL_AUTOTICKS)
    slider.SetMin(10)
    slider.SetMax(100)
    slider.SetValue(100)
    slider.SetTickFreq(10,1)
    self.Bind(wx.EVT_SLIDER, self.OnSetScale, slider)
    szr2.Add(slider, 0, wx.ALL, 5)
    szr2.Fit(ctrls)
    self._viewer = pixel_viewer_3d(parent=self)
    szr.Add(self._viewer, 1, wx.EXPAND)
    szr.Fit(self._viewer)
    self.Fit()

  def __getattr__(self, name):
    return getattr(self._viewer, name)

  def OnChangeRenderStyle(self, event):
    wires = event.GetEventObject().GetValue()
    self._viewer.flag_draw_boxes = not wires
    self._viewer.Refresh()

  def OnChangeScale(self, event):
    scale = event.GetEventObject().GetValue()
    self._viewer.flag_log_scale = scale
    self._viewer.Refresh()

  def OnSetScale(self, event):
    scale = event.GetEventObject().GetValue()
    self._viewer.scale_factor = scale
    self._viewer.Refresh()
