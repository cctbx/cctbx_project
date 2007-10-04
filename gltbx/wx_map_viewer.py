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
import math

class map_view(wx_viewer.wxGLWindow):

  def __init__(self,
               fft_map=None,
               unit_cell=None, raw_map=None,
               frame=None, **kwds):
    assert fft_map is None or (unit_cell is None and raw_map is None)
    import time
    super(map_view, self).__init__(frame, **kwds)
    self.buffer_factor = 2.0
    self.back_colour = (0,)*4
    self._gl_has_been_initialised = False

    if fft_map is not None:
      unit_cell = fft_map.unit_cell()
    o = unit_cell.orthogonalization_matrix()
    self.orthogonaliser = (  o[0:3] + (0,)
                           + o[3:6] + (0,)
                           + o[6:9] + (0,)
                           + (0,0,0,1) )
    a,b,c = unit_cell.parameters()[0:3]
    if fft_map is not None:
      na, nb, nc = fft_map.n_real()
      rho = fft_map.real_map()
    else:
      na, nb, nc = raw_map.accessor().all()
      rho = raw_map
    print "Gridding: %i x %i x %i" % (na,nb,nc)
    def f(iso_level):
      return iso_surface.triangulation(rho, iso_level,
                                       map_extent=(1,1,1))
    self._compute_triangulation = f

    density_stats = maptbx.statistics(rho)
    self.min_density = density_stats.min()
    self.max_density = density_stats.max()
    self.iso_level = density_stats.sigma()
    print "Statistics:"
    print "min: %.3g" % self.min_density
    print "max: %.3g" % self.max_density
    print "sigma: %.3g" % self.iso_level

    p = (0,0,0)
    q = unit_cell.orthogonalize((1,1,1))
    r = unit_cell.orthogonalize((1,0,0))
    s = unit_cell.orthogonalize((0,1,1))
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([p,q,r,s]))

  def iso_level(self):
    return self._iso_level

  def set_iso_level(self, x):
    try:
      if self._iso_level == x: return
    except AttributeError:
      pass
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

  def __init__(self, fft_map=None, unit_cell=None, raw_map=None,
               **kwds):
    self.map = fft_map
    self.unit_cell = unit_cell
    self.raw_map = raw_map
    super(map_viewer, self).__init__(**kwds)

  def init_view_objects(self):
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = map_view(fft_map=self.map,
                                 unit_cell=self.unit_cell,
                                 raw_map=self.raw_map,
                                 frame=self.frame, size=(1200,800))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)

    range = self.view_objects.max_density - self.view_objects.min_density
    n = int(math.log10(range))-2
    p = 5
    self.amplitude = p*10**n
    self.slider = wx.Slider(self.frame,
                            minValue=self.view_objects.min_density
                                       /self.amplitude,
                            maxValue=self.view_objects.max_density
                                       /self.amplitude,
                            value=self.view_objects.iso_level
                                       /self.amplitude,
                            style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    box.Add(self.slider, 0, wx.EXPAND)
    self.multiplier = wx.StaticText(self.frame,
                                   style=wx.ALIGN_CENTER_HORIZONTAL,
                                   label="%i x 10^%i" % (p,n))
    box.Add(self.multiplier, 0, wx.EXPAND|wx.TOP, 10)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    self.slider.Bind(wx.EVT_SCROLL, self.OnSliderMoved)

  def OnSliderMoved(self, event):
    self.iso_level = self.slider.Value * self.amplitude
    self.view_objects.iso_level = self.iso_level



def show_help(): print """\
1. wx_map_viewer name.pickle
  1.a. That file contains an instance of miller.array
  1.b. That file contains a tuple (unit_cell, raw_map)
       where unit_cell is an instance of uctbx.unit_cell and
       raw_map is a 3D flex.double
2. wx_map_viewer name.mtz label

For 1.a and 2, the fft map of the structure factors is iso-surfaced.
For 1.b, the raw_map is displayed
In both cases, within the given unit cell
"""

def run():
  from iotbx import reflection_file_utils
  from iotbx import mtz
  from libtbx.utils import Sorry
  from cctbx import miller
  from cctbx import uctbx
  from scitbx.array_family import flex
  import sys
  import os.path

  file_name = os.path.expanduser(sys.argv[1])
  map_coeffs = None
  if file_name.endswith('.pickle'):
    something = easy_pickle.load(file_name)
    if type(something) is miller.array:
      map_coeffs = something
    else:
      if type(something) is not tuple:
        show_help()
        sys.exit(1)
      unit_cell, raw_map = something
      if (type(unit_cell) is not uctbx.unit_cell
          and raw_map.accessor().nd() != 3):
        show_help()
        sys.exit(1)
  elif file_name.endswith('.mtz'):
    if len(sys.argv[1:]) != 2:
      show_help()
      sys.exit(1)
    label = sys.argv[2]
    all_miller_arrays = mtz.object(file_name).as_miller_arrays()
    labels = reflection_file_utils.label_table(all_miller_arrays)
    map_coeffs = labels.match_data_label(label=label,
                                         command_line_switch="label")
  title = "Electron Density Viewer"
  if map_coeffs is not None:
    map_coeffs.show_summary()
    fft_map = map_coeffs.fft_map(
      resolution_factor=1/2,
      symmetry_flags=maptbx.use_space_group_symmetry)
    a = map_viewer(fft_map=fft_map, title=title)
  else:
    a = map_viewer(unit_cell=unit_cell, raw_map=raw_map, title=title)
  a.MainLoop()

if __name__ == '__main__':
  import sys
  if '--profile' in sys.argv:
    import profile
    import pstats
    sys.argv.remove('--profile')
    profile.run('run()', 'wx_map_viewer.prof')
    p = pstats.Stats('wx_map_viewer.prof')
    p.strip_dirs().sort_stats('time').print_stats(10)
  else:
    run()
