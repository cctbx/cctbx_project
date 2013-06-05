
from __future__ import division
from gltbx import wx_viewer
from gltbx.gl import *
from gltbx.glu import *
import gltbx.util
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
import wx

class pixel_viewer_3d (wx_viewer.wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    self.SetSize((400,400))
    self.SetMinSize((400,400))
    self._img = None
    self.x_center = None
    self.y_center = None
    self.buffer_factor = 2.0
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([[0,0,0],[40,40,40],[40,0,0],[0,40,40]]))
    self.zoom_level = 16

  def InitGL(self):
    gltbx.util.handle_error()
    b = self.background_rgb
    glClearColor(b[0], b[1], b[2], 0.0)
    glDisable(GL_LIGHT0)
    glDisable(GL_LIGHTING)
    glDisable(GL_BLEND)
    glEnable(GL_LINE_SMOOTH)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
    self.initialize_modelview()
    gltbx.util.handle_error()

  def process_pick_point (self) :
    pass

  def set_image (self, img) :
    self._img = img

  def recenter (self, x, y) :
    self.x_center = x
    self.y_center = y
    self.Refresh()

  def DrawGL (self) :
    if (None in [self._img, self.x_center, self.y_center]) : return
    scale_factor = 1000
    values = self._img.get_intensities_in_box(
      x=self.x_center,
      y=self.y_center,
      boxsize=40,
      mag=1)#self.zoom_level)
    glColor3f(1, 1, 1)
    w_range = [0.0,0.0]
    glGetFloatv(GL_LINE_WIDTH_RANGE, w_range)
    line_width = 0.1
    if (w_range[0] > 0.1) :
      line_width = w_range[0]
    glLineWidth(line_width)
    for i, row in enumerate(values) :
      glBegin(GL_LINES)
      for j, value in enumerate(row) :
        if (j > 0) :
          glVertex3f(i, j-1, row[j-1] / scale_factor)
          glVertex3f(i, j, value / scale_factor)
      glEnd()
    for i in range(40) :
      glBegin(GL_LINES)
      for j, row in enumerate(values) :
        if (j > 0) :
          glVertex3f(j-1, i, values[j-1][i] / scale_factor)
          glVertex3f(j, i, values[j][i] / scale_factor)
      glEnd()

class pixel_viewer_3d_frame (wx.MiniFrame) :
  def __init__ (self, *args, **kwds) :
    wx.MiniFrame.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self._viewer = pixel_viewer_3d(parent=self)
    szr.Add(self._viewer, 1, wx.EXPAND)
    szr.Fit(self._viewer)
    self.Fit()

  def __getattr__ (self, name) :
    return getattr(self._viewer, name)
