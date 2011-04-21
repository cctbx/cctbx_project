
from gltbx.wx_viewer import wxGLWindow
import gltbx.gl_managed
import gltbx.util
from gltbx.glu import *
from gltbx.gl import *
import wx
import sys

class hklview (wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wxGLWindow.__init__(self, *args, **kwds)
    # FIXME orthographic is definitely best for this application, but it isn't
    # working properly right now
    #self.orthographic = True
    self.buffer_factor = 2.0
    self.min_slab = 4
    self.min_viewport_use_fraction = 0.1
    self.min_dist = 4.0
    self.flag_show_fog = True
    self.flag_use_lights = True
    self.minimum_covering_sphere = None
    self.spheres_display_list = None
    self.miller_array = None
    self._points = None
    self._radii = None
    self._colors = None

  def InitGL(self):
    gltbx.util.handle_error()
    glClearColor(0.,0.,0.,0.)
    self.minimum_covering_sphere_display_list = None
    glDepthFunc(GL_LESS)
    glEnable(GL_ALPHA_TEST)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    glEnable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    self.initialize_modelview()
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    gltbx.util.handle_error()

  def OnRedrawGL (self, event=None) :
    if self.minimum_covering_sphere is None :
      gltbx.util.handle_error()
      glClear(GL_COLOR_BUFFER_BIT)
      glClear(GL_DEPTH_BUFFER_BIT)
      glFlush()
      self.SwapBuffers()
      gltbx.util.handle_error()
    else :
      wxGLWindow.OnRedrawGL(self, event)

  def initialize_modelview (self) :
    if self.minimum_covering_sphere is not None :
      wxGLWindow.initialize_modelview(self)
    else :
      self.setup_lighting()

  def DrawGL(self):
    if (self.GL_uninitialised) or (self.miller_array is None) :
      return
    self.draw_spheres()

  def set_data (self, miller_array) :
    from gltbx import viewer_utils
    from scitbx.array_family import flex
    self.miller_array = miller_array
    indices = miller_array.indices()
    data = miller_array.data()
    assert isinstance(data, flex.double)
    self._colors = viewer_utils.color_by_property(
      atom_properties=data,
      atoms_visible=flex.bool(data.size(), True),
      color_invisible_atoms=False,
      use_rb_color_gradient=False)
    uc = miller_array.unit_cell()
    self._points = uc.reciprocal_space_vector(indices) * 100.
    abc = uc.parameters()[0:3]
    min_radius = 0.02 / max(abc)
    max_radius = 40 / max(abc)
    scale = max_radius / flex.max(data)
    radii = data * scale
    for i in range(radii.size()) :
      if (radii[i] < min_radius) :
        radii[i] = min_radius
    self._radii = radii
    from scitbx.math import minimum_covering_sphere
    mcs = minimum_covering_sphere(points=self._points,
                                  epsilon=0.1)
    self.minimum_covering_sphere = mcs
    #if recenter_and_zoom :
    #  self.move_rotation_center_to_mcs_center()
    #  self.fit_into_viewport()

  def draw_spheres (self) :
    glMatrixMode(GL_MODELVIEW)
    glShadeModel(GL_SMOOTH)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_NORMALIZE)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    if (self.spheres_display_list is None) :
      self.spheres_display_list = gltbx.gl_managed.display_list()
      self.spheres_display_list.compile()
      colors = self._colors
      radii = self._radii
      for i, hkl in enumerate(self._points) :
        col = list(colors[i]) + [1.0]
        #glColor3f(*colors[i])
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
        glPushMatrix()
        glTranslated(*hkl)
        gltbx.util.SolidSphere(radius=radii[i], slices=20, stacks=20)
        glPopMatrix()
      self.spheres_display_list.end()
    self.spheres_display_list.call()

class App (wx.App) :
  def OnInit (self) :
    self.frame = wx.Frame(None, -1, "HKL viewer", pos=wx.DefaultPosition,
      size=(1024,768))
    self.frame.CreateStatusBar()
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = hklview(self.frame, size=(800,600))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    return True

if (__name__ == "__main__") :
  from iotbx import file_reader
  f = file_reader.any_file(sys.argv[-1])
  ma = None
  for array in f.file_server.miller_arrays :
    if array.is_xray_amplitude_array() or array.is_xray_intensity_array() :
      ma = array
      break
  a = App(0)
  a.view_objects.set_data(ma)
  a.frame.Show()
  a.MainLoop()
