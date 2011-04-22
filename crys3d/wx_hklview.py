
from gltbx.wx_viewer import wxGLWindow
import gltbx.gl_managed
import gltbx.util
from gltbx.glu import *
from gltbx.gl import *
import wxtbx.utils
import wx

class hklview (wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wxGLWindow.__init__(self, *args, **kwds)
    # FIXME orthographic is definitely best for this application, but it isn't
    # working properly right now
    #self.orthographic = True
    self.settings = settings()
    self.settings_window = None
    self.buffer_factor = 2.0
    self.min_slab = 4
    self.min_viewport_use_fraction = 0.1
    self.min_dist = 4.0
    self.flag_show_fog = True
    self.flag_use_lights = True
    self.minimum_covering_sphere = None
    self.spheres_display_list = None
    self.points_display_list = None
    self.miller_array = None
    self.d_min = None
    self._axis_lengths = [1,1,1]
    self._slice_axis = None
    self._slice_index = None
    self._points = None
    self._radii = None
    self._colors = None

  def set_miller_array (self, miller_array) :
    self.miller_array = miller_array
    self.d_min = miller_array.d_min()
    self.construct_reciprocal_space()

  def construct_reciprocal_space (self) :
    from gltbx import viewer_utils
    from scitbx.array_family import flex
    array = self.miller_array.map_to_asu()
    uc = array.unit_cell()
    if (self.settings.expand_data) :
      array = array.expand_to_p1().generate_bijvoet_mates()
    index_span = array.index_span()
    self.hkl_range = index_span.abs_range()
    axes = uc.reciprocal_space_vector(self.hkl_range)
    self._axis_lengths = (axes[0]*100., axes[1]*100., axes[2]*100.)
    indices = array.indices()
    data = array.data()
    assert isinstance(data, flex.double)
    if (self.settings.sqrt_scale_colors) :
      data_for_colors = flex.sqrt(data)
    else :
      data_for_colors = data
    self._colors = viewer_utils.color_by_property(
      atom_properties=data_for_colors,
      atoms_visible=flex.bool(data.size(), True),
      color_invisible_atoms=False,
      use_rb_color_gradient=False)
    if (self.settings.sqrt_scale_radii) :
      data = flex.sqrt(data)
    self._points = uc.reciprocal_space_vector(indices) * 100.
    abc = uc.parameters()[0:3]
    min_radius = 0.02 / max(abc)
    max_radius = 40 / max(abc)
    scale = max_radius / flex.max(data)
    radii = data * scale
    too_small = radii < min_radius
    radii.set_selected(too_small, flex.double(radii.size(), min_radius))
    self._radii = radii
    if (self.settings.show_missing_reflections) :
      missing = array.complete_set().lone_set(array).indices()
      n_missing = missing.size()
      if (n_missing > 0) :
        points_missing = uc.reciprocal_space_vector(missing) * 100.
        self._points.extend(points_missing)
        self._colors.extend(flex.vec3_double(n_missing, (1.,1.,1.)))
        self._radii.extend(flex.double(n_missing, max_radius / 2))
    from scitbx.math import minimum_covering_sphere
    mcs = minimum_covering_sphere(points=self._points,
                                  epsilon=0.1)
    self.minimum_covering_sphere = mcs
    self.spheres_display_list = None
    self.points_display_list = None

  #--- OpenGL methods
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
    if (self.settings.show_axes) :
      self.draw_axes()
    if (self.settings.display_as_spheres) :
      self.draw_spheres()
    else :
      self.draw_points()

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
      assert (colors.size() == radii.size() == self._points.size())
      for i, hkl in enumerate(self._points) :
        col = list(colors[i]) + [1.0]
        #glColor3f(*colors[i])
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
        glPushMatrix()
        glTranslated(*hkl)
        gltbx.util.SolidSphere(radius=radii[i],
          slices=self.settings.sphere_detail,
          stacks=self.settings.sphere_detail)
        glPopMatrix()
      self.spheres_display_list.end()
    self.spheres_display_list.call()

  def draw_points (self) :
    glMatrixMode(GL_MODELVIEW)
    glDisable(GL_LIGHTING)
    if (self.points_display_list is None) :
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      colors = self._colors
      radii = self._radii
      for i, hkl in enumerate(self._points) :
        glColor3f(*colors[i])
        glBegin(GL_LINES)
        glVertex3f(hkl[0] - radii[i], hkl[1], hkl[2])
        glVertex3f(hkl[0] + radii[i], hkl[1], hkl[2])
        glEnd()
        glBegin(GL_LINES)
        glVertex3f(hkl[0], hkl[1] - radii[i], hkl[2])
        glVertex3f(hkl[0], hkl[1] + radii[i], hkl[2])
        glEnd()
        glBegin(GL_LINES)
        glVertex3f(hkl[0], hkl[1], hkl[2] - radii[i])
        glVertex3f(hkl[0], hkl[1], hkl[2] + radii[i])
        glEnd()
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_axes (self) :
    glDisable(GL_LIGHTING)
    glColor3f(1.0, 1.0, 1.0)
    glLineWidth(2)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(self._axis_lengths[0], 0., 0.)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(0.,self._axis_lengths[0], 0.)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(0., 0., self._axis_lengths[0])
    glEnd()
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-self._axis_lengths[0], 0., 0.)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(0.,-self._axis_lengths[0], 0.)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(0., 0., -self._axis_lengths[0])
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  #--- user input and settings
  def process_key_stroke (self, key) :
    if (key == ord("`")) :
      self.edit_settings()

  def update_settings (self) :
    self.construct_reciprocal_space()
    self.Refresh()

  def edit_settings (self) :
    if (self.settings_window is None) :
      self.settings_window = settings_window(self, -1, "Settings")
      self.settings_window.Show()
    self.settings_window.Raise()

class settings (object) :
  def __init__ (self) :
    self.show_axes = True
    self.sqrt_scale_radii = True
    self.sqrt_scale_colors = False
    self.slice_mode = False
    self.expand_data = False
    self.display_as_spheres = True
    self.show_missing_reflections = False
    self.sphere_detail = 20

class settings_window (wxtbx.utils.SettingsToolBase) :
  def add_controls (self) :
    ctrls = self.create_controls(
      setting="show_axes",
      label="Show h,k,l axes")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_radii",
      label="Scale radii to sqrt(I)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_colors",
      label="Scale colors to sqrt(I)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="expand_data",
      label="Expand data to Friedel pairs in P1")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="display_as_spheres",
      label="Display reflections as spheres")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_missing_reflections",
      label="Show missing reflections")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sphere_detail",
      label="Sphere detail level",
      min=4,
      max=20)
    box = wx.BoxSizer(wx.HORIZONTAL)
    box.Add(ctrls[0], 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
    box.Add(ctrls[1], 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.panel_sizer.Add(box)
    return
    ctrls = self.create_controls(
      setting="slice_mode",
      label="Show only a slice through reciprocal space")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    box2.Add(wx.StaticText(self.panel, -1, "View slice:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.hkl_choice = wx.Choice(self.panel, -1, choices=["h","k","l"])
    box2.Add(self.hkl_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(wx.StaticText(self.panel, -1, "="), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.slice_index = wx.SpinCtrl(self.panel, -1)
    self.slice_index.SetValue(0)
    box2.Add(self.slice_index)
    self.panel_sizer.Add(box2)
    self.Bind(wx.EVT_CHOICE, self.OnSetSlice, self.hkl_choice)
    self.Bind(wx.EVT_SPINCTRL, self.OnSetSlice, self.slice_index)

  def OnSetSlice (self, event) :
    axis = self.hkl_choice.GetSelection()
    index = self.slice_index.GetValue()
    self.parent.set_slice(axis, index)

class HKLViewFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.statusbar = self.CreateStatusBar()
    self.toolbar = self.CreateToolBar()
    self.sizer = wx.BoxSizer(wx.VERTICAL)
    self.glwindow = hklview(self, size=(800,600))
    self.sizer.Add(self.glwindow, 1, wx.EXPAND)
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)

  def set_miller_array (self, array) :
    labels = array.info().label_string()
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    self.statusbar.SetStatusText("Data: %s  (Space group: %s  Unit Cell: %s)" %
      (labels, sg, uc))
    self.glwindow.set_miller_array(array)
