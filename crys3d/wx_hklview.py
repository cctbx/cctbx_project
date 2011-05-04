
# TODO:
#  - cached scenes

from gltbx.wx_viewer import wxGLWindow
import gltbx.gl_managed
import gltbx.util
import gltbx.fonts
from gltbx.glu import *
from gltbx.gl import *
from wxtbx import icons
import wxtbx.utils
import wx
from libtbx.utils import Sorry

class scene (object) :
  def __init__ (self, miller_array, settings) :
    self.miller_array = miller_array
    from cctbx import miller
    from scitbx import graphics_utils
    from scitbx.array_family import flex
    array = self.miller_array.map_to_asu()
    if (array.is_xray_intensity_array()) :
      data = array.data()
      data.set_selected(data < 0, flex.double(data.size(), 0.))
      array = array.customized_copy(data=data)
    if (settings.show_data_over_sigma) :
      if (array.sigmas() is None) :
        raise Sorry("sigmas not defined.")
      sigmas = array.sigmas()
      data = array.data()
      array = array.select(sigmas != 0).customized_copy(data=data/sigmas)
    uc = array.unit_cell()
    if (settings.expand_to_p1) :
      array = array.expand_to_p1()
    if (settings.expand_anomalous) :
      array = array.generate_bijvoet_mates()
    slice_selection = None
    if (settings.slice_mode) :
      slice_selection = miller.simple_slice(
        indices=array.indices(),
        slice_axis=settings.slice_axis,
        slice_index=settings.slice_index)
    index_span = array.index_span()
    self.hkl_range = index_span.abs_range()
    self.axes = [ uc.reciprocal_space_vector((self.hkl_range[0],0,0)),
                  uc.reciprocal_space_vector((0,self.hkl_range[1],0)),
                  uc.reciprocal_space_vector((0,0,self.hkl_range[2])) ]
    indices = array.indices()
    data = array.data()
    assert isinstance(data, flex.double) or isinstance(data, flex.bool)
    if isinstance(data, flex.bool) :
      data = flex.double(data.size(), 1.0)
    if (settings.sqrt_scale_colors) :
      data_for_colors = flex.sqrt(data)
    else :
      data_for_colors = data
    if (settings.sqrt_scale_radii) :
      data_for_radii = flex.sqrt(data)
    else :
      data_for_radii = data
    colors = graphics_utils.color_by_property(
      properties=data_for_colors,
      selection=flex.bool(data.size(), True),
      color_all=False,
      use_rb_color_gradient=False)
    if (slice_selection is not None) :
      data = data.select(slice_selection)
      if (data.size() == 0) :
        raise Sorry("No data selected!")
      indices = indices.select(slice_selection)
      if (settings.keep_constant_scale) :
        colors = colors.select(slice_selection)
      else :
        colors = graphics_utils.color_by_property(
          properties=data_for_colors.select(slice_selection),
          selection=flex.bool(data.size(), True),
          color_all=False,
          use_rb_color_gradient=False)
    self.colors = colors
    self.points = uc.reciprocal_space_vector(indices) * 100.
    abc = uc.parameters()[0:3]
    min_radius = 0.20 / max(abc)
    max_radius = 50 / max(abc)
    scale = max_radius / flex.max(data_for_radii)
    if (settings.sqrt_scale_radii) :
      data = flex.sqrt(data)
    radii = data * scale
    too_small = radii < min_radius
    radii.set_selected(too_small, flex.double(radii.size(), min_radius))
    self.radii = radii
    self.max_radius = flex.max(radii)
    self.missing = flex.bool(self.radii.size(), False)
    if (settings.show_missing_reflections) :
      if (settings.show_only_missing) :
        self.colors = flex.vec3_double()
        self.points = flex.vec3_double()
        self.radii = flex.double()
        self.missing = flex.bool()
      complete_set = array.complete_set()
      if (settings.slice_mode) :
        slice_selection = miller.simple_slice(
          indices=complete_set.indices(),
          slice_axis=settings.slice_axis,
          slice_index=settings.slice_index)
        missing = complete_set.select(slice_selection).lone_set(array).indices()
      else :
        missing = complete_set.lone_set(array).indices()
      n_missing = missing.size()
      if (n_missing > 0) :
        points_missing = uc.reciprocal_space_vector(missing) * 100.
        self.points.extend(points_missing)
        self.colors.extend(flex.vec3_double(n_missing, (1.,1.,1.)))
        self.radii.extend(flex.double(n_missing, max_radius / 2))
        self.missing.extend(flex.bool(n_missing, True))

class hklview (wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wxGLWindow.__init__(self, *args, **kwds)
    # FIXME orthographic is definitely best for this application, but it isn't
    # working properly right now
    #self.orthographic = True
    self.settings = self.GetParent().settings
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
    self.scene = None
    self.slice_axis = None
    self.slice_index = None

  def set_miller_array (self, miller_array) :
    self.miller_array = miller_array
    self.d_min = miller_array.d_min()
    self.construct_reciprocal_space()

  def construct_reciprocal_space (self) :
    self.scene = scene(miller_array=self.miller_array,
      settings=self.settings)
    from scitbx.math import minimum_covering_sphere
    mcs = minimum_covering_sphere(points=self.scene.points,
                                  epsilon=0.1)
    self.minimum_covering_sphere = mcs
    self.spheres_display_list = None
    self.points_display_list = None

  #--- OpenGL methods
  def InitGL(self):
    gltbx.util.handle_error()
    if (self.settings.black_background) :
      glClearColor(0.,0.,0.,0.)
    else :
      glClearColor(0.95,0.95,0.95,0.)
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
      colors = self.scene.colors
      radii = self.scene.radii
      max_radius = self.scene.max_radius
      assert (colors.size() == radii.size() == self.scene.points.size())
      for i, hkl in enumerate(self.scene.points) :
        col = list(colors[i]) + [1.0]
        detail = max(4, int(self.settings.sphere_detail*radii[i]/max_radius))
        #glColor3f(*colors[i])
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
        glPushMatrix()
        glTranslated(*hkl)
        gltbx.util.SolidSphere(radius=radii[i],
          slices=detail,
          stacks=detail)
        glPopMatrix()
      self.spheres_display_list.end()
    self.spheres_display_list.call()

  def draw_points (self) :
    glMatrixMode(GL_MODELVIEW)
    glDisable(GL_LIGHTING)
    if (self.points_display_list is None) :
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      colors = self.scene.colors
      radii = self.scene.radii
      assert (colors.size() == radii.size() == self.scene.points.size())
      for i, hkl in enumerate(self.scene.points) :
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
    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]
    gltbx.fonts.ucs_bitmap_8x13.setup_call_lists()
    glDisable(GL_LIGHTING)
    if (self.settings.black_background) :
      glColor3f(1.0, 1.0, 1.0)
    else :
      glColor3f(0.,0.,0.)
    glLineWidth(1.0)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(h_axis[0]*100, h_axis[1]*100, h_axis[2]*100)
    glEnd()
    glRasterPos3f(0.5+h_axis[0]*100, 0.2+h_axis[1]*100, 0.2+h_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("h")
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(k_axis[0]*100, k_axis[1]*100, k_axis[2]*100)
    glEnd()
    glRasterPos3f(0.2+k_axis[0]*100, 0.5+k_axis[1]*100, 0.2+k_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("k")
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(l_axis[0]*100, l_axis[1]*100, l_axis[2]*100)
    glEnd()
    glRasterPos3f(0.2+l_axis[0]*100, 0.5+l_axis[1]*100, 0.2+l_axis[2]*100)
    gltbx.fonts.ucs_bitmap_8x13.render_string("l")
    glEnable(GL_LINE_STIPPLE)
    glLineStipple(4, 0xAAAA)
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.,0.,0.)
    glVertex3f(-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100)
    glEnd()
    glDisable(GL_LINE_STIPPLE)

  #--- user input and settings
  def update_settings (self) :
    self.construct_reciprocal_space()
    if (self.settings.black_background) :
      glClearColor(0.,0.,0.,0.)
    else :
      glClearColor(0.95,0.95,0.95,0.)
    self.Refresh()

  def process_pick_points (self) :
    pass

class hklview_2d (wx.PyPanel) :
  def __init__ (self, *args, **kwds) :
    wx.PyPanel.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_PAINT, self.OnPaint)
    self.scene = None
    self.miller_array = None
    self.parent = self.GetParent()
    self.settings = self.parent.settings

  def set_miller_array (self, array) :
    self.miller_array = array
    if (array is not None) :
      self.construct_reciprocal_space()

  def construct_reciprocal_space (self) :
    self.scene = scene(miller_array=self.miller_array,
      settings=self.settings)

  def update_settings (self) :
    self.construct_reciprocal_space()
    self.Refresh()

  def OnPaint (self, event) :
    if (self.scene is None) :
      return
    if (self.settings.black_background) :
      self.SetBackgroundColour((0,0,0))
    else :
      self.SetBackgroundColour((255,255,255))
    dc = wx.AutoBufferedPaintDCFactory(self)
    dc.Clear()
    gc = wx.GraphicsContext.Create(dc)
    self.paint(gc)

  def paint (self, gc) :
    assert (self.settings.slice_mode)
    if (self.settings.slice_axis == 0) :
      i_x, i_y = 1, 2
      axes = ("k", "l")
    elif (self.settings.slice_axis == 1) :
      i_x, i_y = 0, 2
      axes = ("h", "l")
    else :
      i_x, i_y = 0, 1
      axes = ("h", "k")
    x_max = self.scene.axes[i_x][i_x] * 100.
    y_max = self.scene.axes[i_y][i_y] * 100.
    r = 300
    if (self.settings.show_axes) :
      # FIXME dimensions not right?
      x_end = self.scene.axes[i_x][i_x] * 300., self.scene.axes[i_x][i_y] * 300.
      x_axis = gc.CreatePath()
      x_axis.MoveToPoint(320, 320)
      x_axis.AddLineToPoint(320 + x_end[0], 320 - x_end[1])
      x_axis.CloseSubpath()
      if (self.settings.black_background) :
        gc.SetPen(wx.Pen('white'))
      else :
        gc.SetPen(wx.Pen('black'))
      gc.PushState()
      gc.StrokePath(x_axis)
      gc.PopState()
      y_end = self.scene.axes[i_y][i_x] * 300., self.scene.axes[i_y][i_y] * 300.
      y_axis = gc.CreatePath()
      y_axis.MoveToPoint(320, 320)
      y_axis.AddLineToPoint(320 + y_end[0], 320 - y_end[1])
      y_axis.CloseSubpath()
      gc.PushState()
      gc.StrokePath(y_axis)
      gc.PopState()
    gc.SetPen(wx.TRANSPARENT_PEN)
    if (self.settings.black_background) :
      main_pen = wx.Pen('white')
      main_brush = wx.Brush('white')
      if (self.settings.monochrome) :
        missing_pen = wx.Pen('red')
      else :
        missing_pen = wx.Pen('white')
    else :
      main_pen = wx.Pen('black')
      main_brush = wx.Brush('black')
      if (self.settings.monochrome) :
        missing_pen = wx.Pen('red')
      else :
        missing_pen = wx.Pen('black')
    for k, hkl in enumerate(self.scene.points) :
      x_, y_ = hkl[i_x], hkl[i_y]
      x = 320 + 300 * x_ / x_max
      y = 320 - 300 * y_ / y_max
      r = self.scene.radii[k] * 300 / max(x_max, y_max)
      path = gc.CreatePath()
      path.AddCircle(0, 0, r)
      path.CloseSubpath()
      gc.PushState()
      gc.Translate(x,y)
      if (self.scene.missing[k]) :
        gc.SetBrush(wx.TRANSPARENT_BRUSH)
        gc.SetPen(missing_pen)
        gc.StrokePath(path)
      else :
        if (self.settings.monochrome) :
          gc.SetBrush(main_brush)
          gc.SetPen(main_pen)
        else :
          c = self.scene.colors[k]
          gc.SetBrush(wx.Brush((c[0]*255,c[1]*255,c[2]*255)))
        gc.FillPath(path)
      gc.PopState()

  def save_screen_shot (self, **kwds) :
    pass

########################################################################
class settings (object) :
  def __init__ (self) :
    self.black_background = True
    self.show_axes = True
    self.show_data_over_sigma = False
    self.sqrt_scale_radii = True
    self.sqrt_scale_colors = False
    self.expand_to_p1 = False
    self.expand_anomalous = False
    self.display_as_spheres = True
    self.show_missing_reflections = False
    self.show_only_missing = False
    self.sphere_detail = 20
    self.slice_mode = False
    self.keep_constant_scale = True
    self.slice_axis = 0
    self.slice_index = 0
    self.monochrome = False

class settings_window (wxtbx.utils.SettingsPanel) :
  is_3d_view = True
  def add_controls (self) :
    ctrls = self.create_controls(
      setting="black_background",
      label="Black background")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_axes",
      label="Show h,k,l axes")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_data_over_sigma",
      label="Use I or F over sigma")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_radii",
      label="Scale radii to sqrt(I)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="sqrt_scale_colors",
      label="Scale colors to sqrt(I)")
    self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    if (self.is_3d_view) :
      ctrls = self.create_controls(
        setting="expand_to_p1",
        label="Expand data to P1")
      ctrls2 = self.create_controls(
        setting="expand_anomalous",
        label="show Friedel pairs")
      box = wx.BoxSizer(wx.HORIZONTAL)
      self.panel_sizer.Add(box)
      box.Add(ctrls[0], 0, wx.ALL, 5)
      box.Add(ctrls2[0], 0, wx.ALL, 5)
      ctrls = self.create_controls(
        setting="display_as_spheres",
        label="Display reflections as spheres")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    else :
      ctrls = self.create_controls(
        setting="monochrome",
        label="Monochrome display")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    ctrls = self.create_controls(
      setting="show_missing_reflections",
      label="Show missing reflections")
    ctrls2 = self.create_controls(
      setting="show_only_missing",
      label="only")
    box = wx.BoxSizer(wx.HORIZONTAL)
    self.panel_sizer.Add(box)
    box.Add(ctrls[0], 0, wx.ALL, 5)
    box.Add(ctrls2[0], 0, wx.ALL, 5)
    if (self.is_3d_view) :
      ctrls = self.create_controls(
        setting="sphere_detail",
        label="Sphere detail level",
        min=4,
        max=20)
      box = wx.BoxSizer(wx.HORIZONTAL)
      box.Add(ctrls[0], 0, wx.TOP|wx.BOTTOM|wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
      box.Add(ctrls[1], 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      self.panel_sizer.Add(box)
      ctrls = self.create_controls(
        setting="slice_mode",
        label="Show only a slice through reciprocal space")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
      self.slice_ctrl = ctrls[0]
      ctrls = self.create_controls(
        setting="keep_constant_scale",
        label="Keep scale constant across all slices")
      self.panel_sizer.Add(ctrls[0], 0, wx.ALL, 5)
    box2 = wx.BoxSizer(wx.HORIZONTAL)
    box2.Add(wx.StaticText(self.panel, -1, "View slice:"), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    ctrls = self.create_controls
    self.hkl_choice = wx.Choice(self.panel, -1, choices=["h","k","l"])
    self.hkl_choice.SetSelection(self.settings.slice_axis)
    box2.Add(self.hkl_choice, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    box2.Add(wx.StaticText(self.panel, -1, "="), 0,
      wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.slice_index = wx.SpinCtrl(self.panel, -1)
    self.slice_index.SetValue(self.settings.slice_index)
    box2.Add(self.slice_index)
    self.panel_sizer.Add(box2)
    self.Bind(wx.EVT_CHOICE, self.OnSetSlice, self.hkl_choice)
    self.Bind(wx.EVT_SPINCTRL, self.OnSetSlice, self.slice_index)

  def OnSetSlice (self, event) :
    self.settings.slice_axis = self.hkl_choice.GetSelection()
    self.settings.slice_index = self.slice_index.GetValue()
    if (not self.is_3d_view) or (self.slice_ctrl.GetValue()) :
      self.parent.update_settings()

class HKLViewFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.view_2d = None
    self.statusbar = self.CreateStatusBar()
    self.toolbar = self.CreateToolBar(style=wx.TB_3DBUTTONS|wx.TB_TEXT)
    self.toolbar.SetToolBitmapSize((32,32))
    self.sizer = wx.BoxSizer(wx.HORIZONTAL)
    self.settings = settings()
    self.create_settings_panel()
    self.sizer.Add(self.settings_panel, 0, wx.EXPAND)
    self.create_viewer_panel()
    self.sizer.Add(self.viewer, 1, wx.EXPAND)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Load file",
      bitmap=icons.hkl_file.GetBitmap(),
      shortHelp="Load file",
      kind=wx.ITEM_NORMAL)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Save image",
      bitmap=icons.save_all.GetBitmap(),
      shortHelp="Save image",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnSave, btn)
    menubar = wx.MenuBar(-1)
    self.file_menu = wx.Menu()
    menubar.Append(self.file_menu, "File")
    item = wx.MenuItem(self.file_menu, -1, "Load data...\tCtrl-O")
    self.file_menu.AppendItem(item)
    self.add_view_specific_functions()
    self.SetMenuBar(menubar)
    self.toolbar.Realize()
    self.SetSizer(self.sizer)
    self.sizer.SetSizeHints(self)

  def create_viewer_panel (self) :
    self.viewer = hklview(self, size=(800,600))

  def create_settings_panel (self) :
    self.settings_panel = settings_window(self, -1, style=wx.RAISED_BORDER)

  def add_view_specific_functions (self) :
    item = wx.MenuItem(self.file_menu, -1, "Show 2D view")
    self.file_menu.AppendItem(item)
    self.Bind(wx.EVT_MENU, self.OnShow2D, item)
    btn = self.toolbar.AddLabelTool(id=-1,
      label="Show 2D view",
      bitmap=icons.hklview_2d.GetBitmap(),
      shortHelp="Show 2D view",
      kind=wx.ITEM_NORMAL)
    self.Bind(wx.EVT_MENU, self.OnShow2D, btn)

  def set_miller_array (self, array) :
    if array.is_complex_array() or array.is_hendrickson_lattman_array() :
      raise Sorry("Complex (map) data and HL coefficients not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    self.statusbar.SetStatusText("Data: %s  (Space group: %s  Unit Cell: %s)" %
      (labels, sg, uc))
    self.viewer.set_miller_array(array)
    if (self.view_2d is not None) :
      self.view_2d.set_miller_array(array)

  def update_settings (self, *args, **kwds) :
    self.viewer.update_settings(*args, **kwds)

  def OnLoadFile (self, evt) :
    file_name = wx.FileSelector("Reflections file",
      wildcard="Reflection files (*.mtz, *.sca, *.hkl)|*.mtz;*.sca;*.hkl",
      default_path="",
      flags=wx.OPEN)
    if (file_name != "") :
      from iotbx import file_reader
      from iotbx.gui_tools.reflections import get_array_description
      f = file_reader.any_file(file_name, force_type="hkl")
      f.assert_file_type("hkl")
      arrays = f.file_server.miller_arrays
      valid_arrays = []
      array_info = []
      for array in arrays :
        if array.is_complex_array() or array.is_hendrickson_lattman_array() :
          continue
        labels = array.info().label_string()
        desc = get_array_description(array)
        array_info.append("%s (%s)" % (labels, desc))
        valid_arrays.append(array)
      if (len(valid_arrays) == 0) :
        raise Sorry("No arrays of the supported types in this file.")
      elif (len(valid_arrays) == 1) :
        self.set_miller_array(valid_arrays[0])
      else :
        #dlg = SelectArrayDialog(self, -1, "Select data")
        dlg = wx.SingleChoiceDialog(parent=self,
          message="Please select the data you wish to view:",
          caption="Select data",
          choices=array_info)
        if (dlg.ShowModal() == wx.ID_OK) :
          sel = dlg.GetSelection()
          self.set_miller_array(valid_arrays[sel])
        wx.CallAfter(dlg.Destroy)

  def OnSave (self, evt) :
    self.viewer.save_screen_shot(file_name="hklview",
      extensions=["png"])

  def OnShow2D (self, evt) :
    if (self.view_2d is None) :
      self.view_2d = HKLViewFrame2D(self, -1, "2D HKLview")
      self.view_2d.set_miller_array(self.viewer.miller_array)
      self.view_2d.Show()
    self.view_2d.Raise()

class settings_window_2d (settings_window) :
  is_3d_view = False

class HKLViewFrame2D (HKLViewFrame) :
  def create_viewer_panel (self) :
    self.viewer = hklview_2d(self, -1, size=(640,640))

  def create_settings_panel (self) :
    self.settings.expand_to_p1 = True
    self.settings.expand_anomalous = True
    self.settings.slice_mode = True
    self.settings.black_background = False
    self.settings_panel = settings_window_2d(self, -1, style=wx.RAISED_BORDER)

  def add_view_specific_functions (self) :
    pass

  def OnClose (self, evt) :
    self.Destroy()

  def OnDestroy (self, evt) :
    self.parent.view_2d = None
