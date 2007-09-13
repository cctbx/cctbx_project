from __future__ import division
from gltbx import wx_viewer
import wx
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
from scitbx import iso_surface
from cctbx import maptbx
from libtbx import easy_pickle
import itertools

class map_view(wx_viewer.wxGLWindow):

  def __init__(self, map, iso_level, *args, **kwds):
    super(map_view, self).__init__(*args, **kwds)
    self.buffer_factor = 2.0
    self.back_colour = (0,)*4
    self._gl_has_been_initialised = False

    rho = map.real_map()
    self.density_stats = maptbx.statistics(rho)
    if iso_level is None: iso_level = self.density_stats.mean()

    unit_cell = map.unit_cell()
    o = unit_cell.orthogonalization_matrix()
    self.orthogonaliser = (  o[0:3] + (0,)
                           + o[3:6] + (0,)
                           + o[6:9] + (0,)
                           + (0,0,0,1) )
    a,b,c = unit_cell.parameters()[0:3]
    na, nb, nc = map.n_real()
    grid_size = (1/na, 1/nb, 1/nc)
    def f(iso_level):
      return iso_surface.triangulation(map.real_map(), iso_level, grid_size)
    self._compute_triangulation = f

    self._iso_level = None
    self.iso_level = iso_level

    p = (0,0,0)
    q = unit_cell.orthogonalize((1,1,1))
    r = unit_cell.orthogonalize((1,0,0))
    s = unit_cell.orthogonalize((0,1,1))
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([p,q,r,s]))

  def iso_level(self):
    return self._iso_level

  def set_iso_level(self, x):
    if self._iso_level == x: return
    self._iso_level = x
    self.triangulation = self._compute_triangulation(x)
    if self._gl_has_been_initialised:
      self.OnRedraw()

  iso_level = property(iso_level, set_iso_level)

  def InitGL(self):
    gltbx.util.handle_error()

    glClearColor(*self.back_colour)
    self.initialize_modelview()
    gltbx.util.rescale_normals(fallback_to_normalize=True).enable()
    glEnable(GL_CULL_FACE)
    glCullFace(GL_BACK)

    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
    glLightfv(GL_LIGHT0, GL_POSITION, [1, 1, 1, 0])

    gltbx.util.handle_error()
    self._gl_has_been_initialised = True

  def DrawGL(self):

    glShadeModel(GL_SMOOTH)
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                 [0.1, 0.5, 0.8, 1.])
    #glMaterialfv(GL_FRONT, GL_SPECULAR, [1., 1., 1., 1.])
    #glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.)

    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)

    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glMultTransposeMatrixd(self.orthogonaliser)

    self.draw_unit_cell()
    self.draw_triangulation()
    glPopMatrix()

  def draw_unit_cell(self):
    lw = [0.]
    glGetFloatv(GL_LINE_WIDTH, lw)
    glLineWidth(2)

    glBegin(GL_LINE_LOOP)
    glVertex3f(0,0,0)
    glVertex3f(1,0,0)
    glVertex3f(1,1,0)
    glVertex3f(0,1,0)
    glEnd()
    glBegin(GL_LINE_LOOP)
    glVertex3f(0,0,1)
    glVertex3f(1,0,1)
    glVertex3f(1,1,1)
    glVertex3f(0,1,1)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0,0,0)
    glVertex3f(0,0,1)
    glVertex3f(1,0,0)
    glVertex3f(1,0,1)
    glVertex3f(1,1,0)
    glVertex3f(1,1,1)
    glVertex3f(0,1,0)
    glVertex3f(0,1,1)
    glEnd()
    glLineWidth(lw[0])

  def draw_triangulation(self):
    #glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
    va = gltbx.util.vertex_array(self.triangulation.vertices,
                                 self.triangulation.normals)
    va.draw_triangles(self.triangulation.triangles)


class map_viewer(wx_viewer.App):
  """ Proof of concept for an electron density viewer with a slider to change
  the value of the iso-level.
  """

  def __init__(self, map, iso_level, *args, **kwds):
    self.map = map
    self.iso_level = iso_level
    super(map_viewer, self).__init__(*args, **kwds)

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = map_view(self.map, self.iso_level,
                                 self.frame, size=(600,600))
    self.iso_level = self.view_objects.iso_level
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.slider = wx.Slider(self.frame,
                            minValue=self.view_objects.density_stats.min(),
                            maxValue=self.view_objects.density_stats.max(),
                            value=self.iso_level,
                            style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    box.Add(self.slider, 0, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    self.slider.Bind(wx.EVT_SCROLL, self.OnSliderMoved)

  def OnSliderMoved(self, event):
    self.iso_level = self.slider.Value
    self.view_objects.iso_level = self.iso_level

if __name__ == '__main__':
  """ Loads the file map_coeff.pickle (see cctbx/examples/random_f_calc.py)
  and displays the FFT map based on these coefficients """
  import sys
  if sys.argv[1] == "--debug":
    iso_level = None
  else:
    iso_level = float(sys.argv[1])
  map_coeff = easy_pickle.load("map_coeff.pickle")
  map_coeff.show_summary()
  fft_map = map_coeff.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry)
  a = map_viewer(fft_map, iso_level,
                 title="Electron Density Viewer")
  a.MainLoop()
