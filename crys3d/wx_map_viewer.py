from __future__ import division

import wx

from gltbx import wx_viewer
import gltbx.util
from gltbx import gl_managed
from gltbx.gl import *
from gltbx.glu import *
from gltbx import wx_controllers
from gltbx import fonts

from crys3d import wx_extra
import crys3d.images

from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
from scitbx import iso_surface
from cctbx import maptbx
from libtbx import easy_pickle
from libtbx import object_oriented_patterns as oop

import itertools
import math
import unicodedata
import sys

# this class isn't very well suited to being a mix-in like the classes in
# wx_model_viewer; however, it's a nicely contained stand-alone program,
# so I added the class map_viewer_mixin further down for use in combination
# with the model/selection viewer mixins.  when this file is run directly
# it will still use the map_view class as before. --nat 2009-02-01

class map_view(wx_viewer.wxGLWindow):

  def barycentre_of_min_max(cls, x=0.8):
    return lambda stats: (1-x)*stats.min() + x*stats.max()
  barycentre_of_min_max = classmethod(barycentre_of_min_max)

  def __init__(self,
               fft_map=None,
               unit_cell=None, raw_map=None,
               from_here=None,
               to_there=None,
               periodic=False,
               iso_level=None,
               frame=None,
               wires=True,
               **kwds):
    assert fft_map is None or (unit_cell is None and raw_map is None)
    import time
    super(map_view, self).__init__(frame,
                                   animation_time=0.3,#second
                                   **kwds)
    if iso_level is None: iso_level = self.barycentre_of_min_max()
    self._gl_has_been_initialised = False
    self.buffer_factor = 2.0
    self.background_colour = (0,)*4
    self.wires = wires
    self.material = gl_managed.material_model()

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
    def f(iso_level):
      return iso_surface.triangulation(rho, iso_level,
                                       map_extent=(1,1,1),
                                       from_here=from_here,
                                       to_there=to_there,
                                       periodic=periodic,
                                       ascending_normal_direction=False
                                     )
    self._compute_triangulation = f

    density_stats = maptbx.statistics(rho)
    self.min_density = density_stats.min()
    self.max_density = density_stats.max()
    self.iso_level = iso_level(density_stats)

    if from_here is None: from_here = (0,0,0)
    if to_there is None: to_there = (1,1,1)
    p = unit_cell.orthogonalize(from_here)
    q = unit_cell.orthogonalize(to_there)
    r = unit_cell.orthogonalize((to_there[0], from_here[1], from_here[2]))
    s = unit_cell.orthogonalize((from_here[0], to_there[1], to_there[2]))
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

  def wires(self):
    return self._wires
  def set_wires(self, x):
    try:
      if self._wires == x: return
    except AttributeError:
      pass
    self._wires = x
    if self._gl_has_been_initialised:
      self.OnRedraw()
  wires = property(wires, set_wires)

  def on_lighting_or_material_change(self):
    if self._gl_has_been_initialised:
      self.OnRedraw()

  def InitGL(self):
    gltbx.util.handle_error()

    glClearColor(*self.background_colour)
    self.initialize_modelview()
    gltbx.util.rescale_normals(fallback_to_normalize=True).enable()

    glEnable(GL_DEPTH_TEST)
    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glShadeModel(GL_SMOOTH)
    glLightfv(GL_LIGHT0, GL_POSITION, [0, 0, 1, 0])

    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)

    self.unit_cell_label_fonts = fonts.ucs_bitmap("10x20")
    self.unit_cell_label_shift_from_axes = -0.02
    self.unit_cell_label_fonts.setup_call_lists()

    gltbx.util.handle_error()
    self._gl_has_been_initialised = True

  def DrawGL(self):
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glMultTransposeMatrixd(self.orthogonaliser)
    self.draw_unit_cell()
    self.draw_triangulation()
    glPopMatrix()

  def draw_unit_cell(self):
    glLightfv(GL_LIGHT0, GL_AMBIENT, [0.5, 0.5, 0.5, 1.])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (0., 0., 0., 1.))
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (0, 0, 0, 1))
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 0)

    (x0, y0, z0), (x1, y1, z1) = self.triangulation.bounds

    r,g,b = 0.9, 0.4, 0.3
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (r, g, b, 1.))
    e = self.unit_cell_label_shift_from_axes
    for pos, label in zip([((x0 + x1)/2, y0 + e, z0 + e),
                           (x0 + e, (y0 + y1)/2, z0 + e),
                           (x0 + e, y0 + e, (z0 + z1)/2)],
                          ['a','b','c']):
      self.unit_cell_label_fonts.render_text(pos, label, use_3d_position=True)

    lw = [0.]
    glGetFloatv(GL_LINE_WIDTH, lw)
    glLineWidth(3)

    r,g,b = (0.6,)*3
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (r, g, b, 1.))
    glBegin(GL_LINE_LOOP)
    glVertex3f(x0,y0,z0)
    glVertex3f(x1,y0,z0)
    glVertex3f(x1,y1,z0)
    glVertex3f(x0,y1,z0)
    glEnd()
    glBegin(GL_LINE_LOOP)
    glVertex3f(x0,y0,z1)
    glVertex3f(x1,y0,z1)
    glVertex3f(x1,y1,z1)
    glVertex3f(x0,y1,z1)
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(x0,y0,z0)
    glVertex3f(x0,y0,z1)
    glVertex3f(x1,y0,z0)
    glVertex3f(x1,y0,z1)
    glVertex3f(x1,y1,z0)
    glVertex3f(x1,y1,z1)
    glVertex3f(x0,y1,z0)
    glVertex3f(x0,y1,z1)
    glEnd()
    glLineWidth(lw[0])

  def draw_triangulation(self):
    glLightfv(GL_LIGHT0, GL_AMBIENT, [0., 0., 0., 1.])
    self.material.execute(specular=not self.wires)
    if self.wires:
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
    else:
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
    va = gltbx.util.vertex_array(self.triangulation.vertices,
                                 self.triangulation.normals)
    va.draw_triangles(self.triangulation.triangles)

  def process_pick_points(self):
    pass

########################################################################
# MIXIN FOR MAPS
#
# this is combined with a subclass of wx_model_viewer.selection_viewer_mixin
# to display the refinement status in phenix.refine
# TODO: replace map list and associated lists with dict
class map_viewer_mixin (wx_viewer.wxGLWindow) :
  initialize_map_viewer_super = True
  def __init__ (self, *args, **kwds) :
    if self.initialize_map_viewer_super :
      wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    # various data objects
    self.unit_cell = None
    self.orthogonaliser = None
    self.maps       = []
    self.iso_levels = []
    self.map_colors = []
    self.materials  = []
    self.triangles  = []
    self._triangle_tmp = []
    # user settings
    self.mesh_line_width = 0.25 # very buggy on OS X + NVidia (and ???)
    self.map_radius      = 10.0
    self.buffer_factor   = 2 # what does this even do?
    # flags
    self.flag_show_maps = True
    self.flag_use_materials = False
    self.flag_show_rotation_center = True
    # maps look much better when OpenGL materials are used, but this screws
    # up the model rendering so it's off by default
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([[0,0,0],[100,100,100],[100,0,0],[0,100,100]]))

  def InitGL (self) :
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    gltbx.util.rescale_normals(fallback_to_normalize=True).enable()
    glEnable(GL_DEPTH_TEST)
    glShadeModel(GL_SMOOTH)
    vendor = glGetString(GL_VENDOR)
    if sys.platform == "darwin" and vendor.startswith("NVIDIA") :
      glDisable(GL_LINE_SMOOTH)
    else :
      glEnable(GL_LINE_SMOOTH)
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)
    self.initialize_modelview()

  def DrawGL (self) :
    if self.unit_cell is None or self.orthogonaliser is None :
      return
    if self.flag_show_maps :
      self.draw_maps()
    if self.flag_show_rotation_center :
      self.draw_rotation_center()

  def OnTranslate (self, event) :
    wx_viewer.wxGLWindow.OnTranslate(self, event)
    self.update_map_objects()

  def draw_rotation_center(self):
    font = gltbx.fonts.ucs_bitmap_10x20
    font.setup_call_lists()
    glColor3f(0, 1.0, 0)
    glRasterPos3f(*self.rotation_center)
    font.render_string("+")
    return
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    rc = self.rotation_center
    (x,y,z) = (rc[0], rc[1], rc[2])
    glTranslatef(x,y,z)
    #glScalef(a, b, c)
    glBegin(GL_LINES)
    f = 0.5
    glColor3f(0, 1.0, 0)
    glVertex3f(-f, 0, 0)
    glVertex3f(f, 0 ,0)
    glVertex3f(0, -f, 0)
    glVertex3f(0, f, 0)
    glVertex3f(0, 0, -f)
    glVertex3f(0, 0, f)
    glEnd()
    glPopMatrix()

  def setup_fog (self) :
    if self.flag_show_fog :
      near, far = self.get_clipping_distances()
      # TODO: this needs work.
      fog_start = 0.25*(far - near) + near + self.fog_start_offset
      fog_end = self.far - self.fog_end_offset
      glMatrixMode(GL_MODELVIEW)
      glEnable(GL_FOG)
      glFogi(GL_FOG_MODE, GL_LINEAR)
      glFogf(GL_FOG_START, fog_start)
      glFogf(GL_FOG_END, fog_end)
      glFogfv(GL_FOG_COLOR, [self.r_back, self.g_back, self.b_back, 1.0])
    else :
      glDisable(GL_FOG)

  def OnMouseWheel (self, event) :
    scale = event.GetWheelRotation()
    if event.ShiftDown() :
      self.fog_end_offset -= scale
    else :
      self.clip_far -= scale

  def process_key_stroke (self, key) :
    if key == wx.WXK_UP :
      self.fog_start_offset += 1
      self.OnRedrawGL()
    elif key == wx.WXK_DOWN :
      self.fog_start_offset -= 1
      self.OnRedrawGL()
    elif key == wx.WXK_LEFT :
      self.fog_end_offset -= 1
      self.OnRedrawGL()
    elif key == wx.WXK_RIGHT :
      self.fog_end_offset += 1
      self.OnRedrawGL()

  def initialize_unit_cell (self, map) :
    self.unit_cell = map.unit_cell()
    o = self.unit_cell.orthogonalization_matrix()
    self.orthogonaliser = (  o[0:3] + (0,)
                           + o[3:6] + (0,)
                           + o[6:9] + (0,)
                           + (0,0,0,1) )
    p = self.unit_cell.orthogonalize((0,0,0))
    q = self.unit_cell.orthogonalize((1,1,1))
    r = self.unit_cell.orthogonalize((1, 0, 0))
    s = self.unit_cell.orthogonalize((0, 1, 1))
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([p,q,r,s]))

  def update_map_objects (self, redraw=True) :
    self.compute_triangulation()
    # this is done to prevent thread clashes (+ ensuing crashes)
    self.triangles = self._triangle_tmp
    if redraw :
      self.OnRedraw()

  def set_iso_levels (self, levels) :
    assert len(levels) == len(self.maps)
    for i, map_contour_level in enumerate(levels) :
      self.iso_levels[i] = map_contour_level

  def set_map_colors (self, colors) :
    assert len(colors) == len(self.maps)
    for i, map_color in enumerate(colors) :
      self.map_colors[i] = map_color
      self.materials[i] = gltbx.gl_managed.material_model(
        front_colour=map_color,
        back_colour=map_color)

  def set_maps (self, maps, reset_unit_cell=False) :
    assert len(maps) != 0
    self.maps = maps
    if (len(self.iso_levels) != len(maps) or
        len(self.map_colors) != len(maps)) :
      self.iso_levels = [ 1.0 for m in maps ]
      self.map_colors = [ (1., 1., 1.) for m in maps ]
      self.materials = [ None for m in maps ]
    if self.unit_cell is None or reset_unit_cell :
      self.initialize_unit_cell(self.maps[0])

  def compute_triangulation (self) :
    if len(self.maps) == 0 :
      return
    self._triangle_tmp = []
    r = self.map_radius
    c = self.rotation_center # set in wxGLWindow (default is origin)
    min = [ c[x] - float(r) for x in [0, 1, 2] ]
    max = [ c[x] + float(r) for x in [0, 1, 2] ]
    map_boundaries_cart = flex.vec3_double([min,max])
    bounds = self.unit_cell.fractionalize(sites_cart=map_boundaries_cart)
    for i, raw_map in enumerate(self.maps) :
      rho = raw_map.real_map()
      iso_level = self.iso_levels[i]
      triangulation = iso_surface.triangulation(rho,
                                iso_level,
                                map_extent=(1,1,1),
                                from_here=bounds[0],
                                to_there=bounds[1],
                                periodic=True,
                                ascending_normal_direction=False
                              )
      self._triangle_tmp.append(triangulation)

  def draw_maps (self) :
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glMultTransposeMatrixd(self.orthogonaliser)
    gltbx.util.handle_error()
    if self.flag_use_materials :
      glLightfv(GL_LIGHT0, GL_AMBIENT, [0., 0., 0., 1.])
    for i, triangulation in enumerate(self.triangles) :
      glLineWidth(0.2)
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
      if self.flag_use_materials :
        self.materials[i].execute(specular=False)
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
      else :
        glColor3f(*self.map_colors[i])
        if self.flag_use_lights :
          glDisable(GL_LIGHTING)
          glDisable(GL_LIGHT0)
          glDisable(GL_BLEND)
      va = gltbx.util.vertex_array(triangulation.vertices,
                                   triangulation.normals)
      va.draw_triangles(triangulation.triangles)
    glPopMatrix()

########################################################################
# CLASSES AND METHODS FOR RUNNING map_view
#
class App(wx_viewer.App):

  def __init__(self,
               fft_map=None,
               unit_cell=None, raw_map=None,
               from_here=None,
               to_there=None,
               periodic=False,
               iso_level=None,
               wires=True,
               default_size=(600,600),
               **kwds):
    self._make_view_objects = lambda: map_view(
      fft_map=fft_map,
      unit_cell=unit_cell, raw_map=raw_map,
      from_here=from_here,
      to_there=to_there,
      periodic=periodic,
      wires=wires,
      iso_level=iso_level,
      frame=self.frame,
      size=self.default_size)
    super(App, self).__init__(default_size=default_size, **kwds)


  def init_view_objects(self):
    # create widgets
    view = self.view_objects = self._make_view_objects()
    range = view.max_density - view.min_density
    if range != 0:
      n = int(math.log10(range))-3
      p = 1
      self.amplitude = p*10**n
    else:
      p,n = 1,0
      self.amplitude = 1

    if wx.Platform == '__WXGTK__': slider_size=(120,-1)
    else: slider_size=(-1,-1)

    self.tools = wx_extra.InspectorToolFrame(self.frame)

    iso_level_inspector = wx_extra.Inspector(
      self.tools, label="Electron density Iso-levels")
    iso_level_pane = iso_level_inspector.GetPane()
    self.iso_level_slider = wx.Slider(
                            iso_level_pane,
                            size=slider_size,
                            minValue=view.min_density / self.amplitude,
                            maxValue=view.max_density / self.amplitude,
                            value=view.iso_level / self.amplitude,
                            style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    self.front_colour_picker = wx.lib.colourselect.ColourSelect(
      iso_level_pane)
    self.back_colour_picker = wx.lib.colourselect.ColourSelect(
      iso_level_pane)
    self.swap_colour_btn = wx.Button(iso_level_pane, label="< swap >")
    self.wires_btn = wx.CheckBox(iso_level_pane, label="Wires")
    if self.view_objects.wires:
      self.wires_btn.Set3StateValue(wx.CHK_CHECKED)
    else:
      self.wires_btn.Set3StateValue(wx.CHK_UNCHECKED)
    self.wires_btn.Bind(wx.EVT_CHECKBOX,
                        self.on_wires_changed)
    lbl_fmt = "Iso-level %s %%s" % unicodedata.lookup('MULTIPLICATION SIGN')
    if n == 1:
      multiplier = "%i" % p
    else:
      multiplier = "%ie%i" % (p,n)
    lbl = lbl_fmt % multiplier
    self.multiplier = wx.StaticText(iso_level_pane,
                                    style=wx.ALIGN_CENTER_HORIZONTAL,
                                    label="label")#lbl)
    self.iso_level_slider.Bind(wx.EVT_SCROLL, self.on_iso_level_changed)

    opengl_inspector = wx_extra.Inspector(self.tools,
                                          label="Lighting and Materials")
    opengl_pane = opengl_inspector.GetPane()


    self.ambient_slider = wx.Slider(opengl_pane, size=slider_size,
                                    style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    self.diffuse_slider = wx.Slider(opengl_pane, size=slider_size,
                                    style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    self.specular_slider = wx.Slider(opengl_pane, size=slider_size,
                                     style=wx.SL_AUTOTICKS|wx.SL_LABELS)
    self.specular_focus_slider = wx.Slider(opengl_pane, size=slider_size,
                                           style=wx.SL_AUTOTICKS)

    self.material_ctrl = wx_controllers.material(view.material,
                                                 self.front_colour_picker,
                                                 self.back_colour_picker,
                                                 self.swap_colour_btn,
                                                 self.ambient_slider,
                                                 self.diffuse_slider,
                                                 self.specular_slider,
                                                 self.specular_focus_slider,
                                                 on_change=view.OnRedrawGL)

    # lay out main window
    box = wx.BoxSizer(wx.VERTICAL)
    box.Add(view, 1, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)

    # lay out toolbox palette

    # Iso-levels
    box = wx.BoxSizer(wx.VERTICAL)

    box.Add(self.wires_btn, 0, wx.ALL, 5)
    box.AddSpacer(10)

    surface_box = wx.BoxSizer(wx.HORIZONTAL)

    iso_level_box = wx.BoxSizer(wx.VERTICAL)
    iso_level_box.Add(self.iso_level_slider, 0, wx.BOTTOM, 5)
    iso_level_box.Add(self.multiplier, 0, wx.CENTER)
    surface_box.Add(iso_level_box, flag=wx.ALL, border=5)

    surface_box.AddSpacer(10)

    colour_box = wx.BoxSizer(wx.HORIZONTAL)
    colour_box.Add(wx.StaticText(iso_level_pane, label="Front"), 0,
                   wx.ALL|wx.CENTER, 2)
    colour_box.Add(self.front_colour_picker, 0, wx.ALL|wx.CENTER, 2)

    colour_box.Add(self.swap_colour_btn, 0, wx.ALL, 5)

    colour_box.Add(self.back_colour_picker, 0, wx.ALL|wx.CENTER, 2)
    colour_box.Add(wx.StaticText(iso_level_pane, label="Back"), 0,
                   wx.ALL|wx.CENTER, 2)
    surface_box.Add(colour_box, flag=wx.ALL, border=5)

    box.Add(surface_box)

    iso_level_pane.SetSizer(box)

    # Lighting and materials pane
    box = wx.FlexGridSizer(rows=0, cols=2, vgap=10, hgap=10)

    cell_box = wx.BoxSizer(wx.VERTICAL)
    cell_box.Add(self.ambient_slider, 0, wx.EXPAND|wx.BOTTOM, border=5)
    cell_box.Add(wx.StaticText(opengl_pane, label='Ambient (%)'),
                 flag=wx.CENTER)
    box.Add(cell_box)

    cell_box = wx.BoxSizer(wx.VERTICAL)
    cell_box.Add(self.diffuse_slider, 0, wx.EXPAND|wx.BOTTOM, border=5)
    cell_box.Add(wx.StaticText(opengl_pane, label='Diffuse (%)'),
                 flag=wx.CENTER)
    box.Add(cell_box)

    cell_box = wx.BoxSizer(wx.VERTICAL)
    cell_box.Add(self.specular_slider, 0, wx.EXPAND|wx.BOTTOM, border=5)
    cell_box.Add(wx.StaticText(opengl_pane, label='Specular (%)'),
                 flag=wx.CENTER)
    box.Add(cell_box)

    cell_box = wx.BoxSizer(wx.VERTICAL)
    box_sf = wx.BoxSizer(wx.VERTICAL)
    box_sf.Add(self.specular_focus_slider, 0, wx.EXPAND|wx.BOTTOM, border=5)
    box_sf_1 = wx.BoxSizer(wx.HORIZONTAL)
    box_sf_1.Add(wx.StaticText(opengl_pane, label="diffuse"), 0, 0)
    box_sf_1.AddStretchSpacer(10)
    box_sf_1.Add(wx.StaticText(opengl_pane, label="sharp"), 0, 0)
    box_sf.Add(box_sf_1, 0, wx.EXPAND)
    cell_box.Add(box_sf, flag=wx.BOTTOM, border=5)
    cell_box.Add(wx.StaticText(opengl_pane, label='Specular focus'),
                 flag=wx.CENTER)
    box.Add(cell_box)

    opengl_pane.SetSizer(box)

    # Layout tool window
    self.tools.Layout()

    # create toolbar button to show/hide the toolbox palette
    tb = self.frame.GetToolBar()
    self.toogle_tools_id = wx.NewId()
    tb.AddSimpleTool(
      id=self.toogle_tools_id,
      bitmap=crys3d.images.inspector_img.as_wx_Bitmap(),
      shortHelpString="Show/Hide inspectors",)
    tb.Realize()

    # final touch
    self.tools.move_parent_out_of_the_way()
    self.tools.Show()

  def OnToolClick(self, event):
    id = event.GetId()
    if id == self.toogle_tools_id:
      if self.tools.Shown: self.tools.Hide()
      else: self.tools.Show()
    else:
      super(App, self).OnToolClick(event)

  def on_iso_level_changed(self, event):
    self.view_objects.iso_level = self.iso_level_slider.Value * self.amplitude

  def on_wires_changed(self, e):
    self.view_objects.wires = self.wires_btn.IsChecked()

def display(*args, **kwds):
  App(*args, **kwds).MainLoop()


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
  fft_map = None
  unit_cell = None
  raw_map = None
  if file_name.endswith('.pickle'):
    something = easy_pickle.load(file_name)
    if type(something) is miller.array:
      map_coeffs = something
    else:
      if type(something) is not tuple:
        if type(something) is miller.fft_map:
          unit_cell, raw_map = something.unit_cell(), something.real_map()
        else:
          show_help()
          sys.exit(1)
      else:
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
  options = {}
  if 1:
    options['wires'] = False
  if map_coeffs is not None:
    map_coeffs.show_summary()
    fft_map = map_coeffs.fft_map(
      resolution_factor=1/2,
      symmetry_flags=maptbx.use_space_group_symmetry)
  display(fft_map=fft_map,
          unit_cell=unit_cell, raw_map=raw_map,
          **options)

if __name__ == '__main__':
  run()
