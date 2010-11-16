from __future__ import division

import sys, os
from crys3d.wx_selection_editor import selection_editor_mixin
import iotbx.phil
from gltbx.wx_viewer import wxGLWindow
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
from scitbx import iso_surface
from libtbx import adopt_init_args
import wx

viewer_phil = iotbx.phil.parse("""
  include scope crys3d.wx_selection_editor.viewer_phil
""", process_includes=True)

class map_data (object) :
  def __init__ (self, map, is_difference_map=False, radius=10.0) :
    adopt_init_args(self, locals())
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
    self.iso_levels = [1]
    self.colors = [(1,1,1)]

  def set_iso_levels (self, levels) :
    self.iso_levels = levels
    if len(self.colors) != len(levels) :
      self.colors = [ (1.0,1.0,1.0) for x in levels ]

  def increment_iso_levels (self, inc) :
    for i, cutoff in enumerate(self.iso_levels) :
      if self.is_difference_map :
        if cutoff < 0 :
          self.iso_levels[i] -= inc
        elif cutoff > 0 :
          self.iso_levels[i] += inc
      else :
        self.iso_levels[i] += inc

  def set_colors (self, colors) :
    assert len(self.iso_levels) == len(colors)
    self.colors = colors

  def update_map_data (self, map) :
    self.map = map

  def get_scene_data (self, rotation_center) :
    triangles = []
    r = self.radius
    c = rotation_center
    min = [ c[x] - float(r) for x in [0, 1, 2] ]
    max = [ c[x] + float(r) for x in [0, 1, 2] ]
    map_boundaries_cart = flex.vec3_double([min,max])
    bounds = self.unit_cell.fractionalize(sites_cart=map_boundaries_cart)
    rho = self.map.real_map()
    for iso_level in self.iso_levels :
      triangulation = iso_surface.triangulation(rho,
                                iso_level,
                                map_extent=(1,1,1),
                                from_here=bounds[0],
                                to_there=bounds[1],
                                periodic=True,
                                ascending_normal_direction=False
                              )
      triangles.append(triangulation)
    return map_scene(triangles, self.orthogonaliser, self.colors)

class map_scene (object) :
  def __init__ (self, triangles, orthogonaliser, colors) :
    assert len(triangles) == len(colors)
    adopt_init_args(self, locals())
    self.flag_use_materials = False
    self.clear_lists()

  def clear_lists (self) :
    self.mesh_display_list = None

  def draw_mesh (self) :
    if self.mesh_display_list is None :
      self.mesh_display_list = gltbx.gl_managed.display_list()
      self.mesh_display_list.compile()
      glMatrixMode(GL_MODELVIEW)
      glLineWidth(0.2)
      gltbx.util.handle_error()
      glPushMatrix()
      try :
        glMultTransposeMatrixd(self.orthogonaliser)
        for i, triangulation in enumerate(self.triangles) :
          glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
          if self.flag_use_materials :
            pass
          #  self.materials[i].execute(specular=False)
          #  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE)
          else :
            glColor3f(*self.colors[i])
          va = gltbx.util.vertex_array(triangulation.vertices,
                                     triangulation.normals)
          va.draw_triangles(triangulation.triangles)
      finally :
        glPopMatrix()
      self.mesh_display_list.end()
    self.mesh_display_list.call()

#-----------------------------------------------------------------------
class map_viewer_mixin (wxGLWindow) :
  initialize_map_viewer_super = True
  def __init__ (self, *args, **kwds) :
    if self.initialize_map_viewer_super :
      wxGLWindow.__init__(self, *args, **kwds)
    # various data objects
    self.map_ids     = []
    self.map_objects = []
    self.map_scenes  = {}
    self.show_object = {}
    self.map_ctrls = {}
    # user settings
    self.mesh_line_width = 0.25 # very buggy on OS X + NVidia (and ???)
    self.selected_map_id = None
    self.update_maps = False
    self.flag_show_maps = True
    self.flag_smooth_lines = True
    self.flag_use_materials = False
    self.flag_show_rotation_center = True
    self.minimum_covering_sphere = minimum_covering_sphere(
      flex.vec3_double([[0,0,0],[100,100,100],[100,0,0],[0,100,100]]))

  def InitGL (self) :
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    gltbx.util.rescale_normals(fallback_to_normalize=True).enable()
    glEnable(GL_DEPTH_TEST)
    glShadeModel(GL_SMOOTH)
    vendor = glGetString(GL_VENDOR)
    if sys.platform == "darwin" and vendor.startswith("NVIDIA") :
      self.flag_smooth_lines = False
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)
    self.initialize_modelview()

  def OnRedrawGL (self, event=None) :
    self.check_and_update_map_scenes()
    wxGLWindow.OnRedrawGL(self, event)

  def check_and_update_map_scenes (self) :
    if self.update_maps :
      self.update_map_scenes()
      self.update_maps = False

  def DrawGL (self) :
    if len(self.map_scenes) == 0 :
      return
    if self.flag_show_maps :
      self.draw_maps()
    if self.flag_show_rotation_center :
      self.draw_rotation_center()

  def OnTranslate (self, event) :
    wxGLWindow.OnTranslate(self, event)
    self.update_map_scenes()

  def draw_rotation_center(self):
    font = gltbx.fonts.ucs_bitmap_10x20
    font.setup_call_lists()
    glColor3f(0, 1.0, 0)
    glRasterPos3f(*self.rotation_center)
    font.render_string("+")

  def process_key_stroke (self, key) :
    pass

  def add_map (self, map_id, map, is_difference_map=False) :
    if map_id in self.map_ids :
      self.delete_map(map_id)
    map_object = map_data(map,
      is_difference_map=is_difference_map,
      radius=self.settings.opengl.map_radius)
    self.map_ids.append(map_id)
    self.map_objects.append(map_object)
    self.show_object[map_id] = True
    self.update_maps = True

  def delete_map (self, map_id) :
    if map_id in self.map_ids :
      i = self.map_ids.index(map_id)
      self.map_ids.pop(i)
      self.map_objects.pop(i)
      self.show_object.pop(map_id)
      if map_id in self.map_scenes :
        self.map_scenes.pop(map_id)
    self.update_maps = True

  def update_map (self, map_id, map) :
    if not map_id in self.map_ids :
      self.add_map(map_id, map)
    else :
      map_object = self.get_map(map_id)
      map_object.update_map_data(map)
      self.update_maps = True

  def show_map_ctrls (self, map_id) :
    map_object = self.get_map(map_id)
    if map_object is None :
      print "oops"
      return
    elif self.map_ctrls.get(map_id) is not None :
      self.map_ctrls[map_id].Raise()
    else :
      frame = MapEditor(self, map_object, map_id)
      frame.Show()
      self.map_ctrls[map_id] = frame

  def update_map_from_miller_array (self, map_id, map_coeffs,
      resolution_factor=0.33) :
    assert map_coeffs.is_complex_array()
    fft_map = map_coeffs.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    self.update_map(map_id, fft_map)

  def hide_maps (self, object_id=None) :
    for map_id in self.map_ids :
      if object_id is None or map_id == object_id :
        self.show_object[map_id] = False

  def iter_maps (self) :
    for (map_id, map_object) in zip(self.map_ids, self.map_objects) :
      yield (map_id, map_object)

  def get_map (self, map_id) :
    for (object_id, map_object) in self.iter_maps() :
      if object_id == map_id :
        return map_object

  def set_map_levels_and_colors (self, map_id, levels, colors) :
    map = self.get_map(map_id)
    map.set_iso_levels(levels)
    map.set_colors(colors)
    self.update_maps = True

  def set_selected_map (self, map_id) :
    if map_id is None :
      self.selected_map_id = None
    else :
      self.selected_map_id = map_id
    for other_id, ctrl_frame in self.map_ctrls.iteritems() :
      if ctrl_frame is not None :
        ctrl_frame.refresh_attach_mouse()

  def get_selected_map (self) :
    return self.get_map(self.selected_map_id)

  def increment_map_iso_levels (self, inc) :
    if self.get_selected_map() is not None :
      self.get_selected_map().increment_iso_levels(inc)
      if self.map_ctrls.get(self.selected_map_id) is not None :
        self.map_ctrls[self.selected_map_id].refresh_iso_levels()
      self.update_maps = True

  def update_map_scenes (self) :
    for object_id, map in self.iter_maps() :
      map.radius = self.settings.opengl.map_radius
      scene = map.get_scene_data(self.rotation_center)
      self.map_scenes[object_id] = scene

  def draw_maps (self) :
    gltbx.util.handle_error()
    if self.flag_use_materials :
      glLightfv(GL_LIGHT0, GL_AMBIENT, [0., 0., 0., 1.])
    if self.flag_use_lights :
      glDisable(GL_LIGHTING)
      glDisable(GL_LIGHT0)
      glDisable(GL_BLEND)
    if not self.flag_smooth_lines :
      glDisable(GL_LINE_SMOOTH)
    for map_id, scene in self.map_scenes.iteritems() :
      if self.show_object[map_id] :
        scene.draw_mesh()

#-----------------------------------------------------------------------
class model_and_map_viewer (selection_editor_mixin, map_viewer_mixin) :
  initialize_map_viewer_super = False
  def __init__ (self, *args, **kwds) :
    selection_editor_mixin.__init__(self, *args, **kwds)
    map_viewer_mixin.__init__(self, *args, **kwds)
    self.buffer_factor = 1
    self._debug_mode = ("GLTBX_DEBUG_MODE" in os.environ)

  def OnMouseWheel (self, event) :
    scale = event.GetWheelRotation()
    if self.selected_map_id is not None and not event.AltDown() :
      self.increment_map_iso_levels(scale * 0.1)
    elif event.ShiftDown() :
      self.fog_end_offset -= scale
    else :
      self.slab_scale += 0.01 * scale
      if self.slab_scale > 1.0 :
        self.slab_scale = 1.0
      elif self.slab_scale < 0.01 :
        self.slab_scale = 0.01
    self.OnRedrawGL()

  def InitGL (self) :
    selection_editor_mixin.InitGL(self)
    gltbx.util.rescale_normals(fallback_to_normalize=True).enable()
    glShadeModel(GL_SMOOTH)
    vendor = glGetString(GL_VENDOR)
    if sys.platform == "darwin" and vendor.startswith("NVIDIA") :
      self.flag_smooth_lines = False
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)

  def OnRedrawGL (self, event=None) :
    if self._debug_mode :
      self.show_stack_sizes()
    self.check_and_update_map_scenes()
    selection_editor_mixin.OnRedrawGL(self, event)

  def DrawGL (self) :
    selection_editor_mixin.DrawGL(self)
    map_viewer_mixin.DrawGL(self)

  def process_key_stroke (self, key) :
    selection_editor_mixin.process_key_stroke(self, key)
    map_viewer_mixin.process_key_stroke(self, key)
    self.OnRedrawGL()

  def update_mcs (self, *args, **kwds) :
    self.update_maps = True
    selection_editor_mixin.update_mcs(self, *args, **kwds)

  def recenter_on_atom (self, *args, **kwds) :
    self.update_maps = True
    selection_editor_mixin.recenter_on_atom(self, *args, **kwds)

  def hide_all (self, *args, **kwds) :
    self.hide_models(*args, **kwds)
    self.hide_maps(*args, **kwds)

  def hide_others (self, object_id=None) :
    if object_id is None :
      self.hide_models()
      self.hide_maps()
    else :
      for current_object_id in self.model_ids+self.map_ids :
        if current_object_id != object_id :
          self.show_object[current_object_id] = False

  def toggle_visibility (self, show_object, object_id=None) :
    for current_object_id in self.model_ids+self.map_ids :
      if current_object_id == object_id :
        self.show_object[current_object_id] = show_object

  def update_all_settings (self, params, redraw=False) :
    selection_editor_mixin.update_settings(self, params, redraw)
    self.update_maps = True

########################################################################
class MapEditor (wx.MiniFrame) :
  def __init__ (self, parent, map_object, map_id) :
    import wx.lib.colourselect
    adopt_init_args(self, locals())
    wx.MiniFrame.__init__(self, parent, -1, map_id,
      style=wx.CAPTION|wx.CLOSE_BOX|wx.RAISED_BORDER)
    panel = wx.Panel(self, -1)
    main_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(main_sizer)
    self.ctrls = []
    for iso_level, color in zip(map_object.iso_levels, map_object.colors) :
      level_sizer = wx.BoxSizer(wx.HORIZONTAL)
      slider = wx.Slider(panel,
                         size=(160,-1),
                         minValue=-100,
                         maxValue=100,
                         value=int(iso_level * 10),
                         style=wx.SL_AUTOTICKS)
      level_txt = wx.TextCtrl(panel,
                              size=(64,-1),
                              style=wx.TE_READONLY|wx.TE_RIGHT)
      level_txt.SetValue("%.2f" % iso_level)
      initial_color = [ int(x*255) for x in color ]
      color_ctrl = wx.lib.colourselect.ColourSelect(panel,
                                                    colour=initial_color)
      level_sizer.Add(slider, 0, wx.ALL, 5)
      level_sizer.Add(level_txt, 0, wx.ALL, 5)
      level_sizer.Add(color_ctrl, 0, wx.ALL, 5)
      main_sizer.Add(level_sizer, 1, wx.ALL, 5)
      self.ctrls.append((slider, level_txt, color_ctrl))
    self.attach_mouse_wheel = wx.CheckBox(panel,
      label="Attach mouse wheel to this map")
    self.attach_mouse_wheel.SetValue(map_id == parent.selected_map_id)
    main_sizer.Add(self.attach_mouse_wheel, 0, wx.ALL, 5)
    main_sizer.Fit(panel)
    self.Fit()
    self.Bind(wx.EVT_CHECKBOX, self.OnAttachMouse, self.attach_mouse_wheel)
    self.Bind(wx.EVT_SLIDER, self.OnUpdate)
    self.Bind(wx.lib.colourselect.EVT_COLOURSELECT, self.OnUpdate)
    self.Bind(wx.EVT_CLOSE, self.OnClose)
    self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

  def OnUpdate (self, event) :
    source_ctrl = event.GetEventObject()
    for i, (s, t, c) in enumerate(self.ctrls) :
      if source_ctrl is s :
        new_level = s.GetValue() / 10.0
        self.map_object.iso_levels[i] = new_level
        t.SetValue("%.2f" % new_level)
        self.parent.update_maps = True
      elif source_ctrl is c :
        new_color = c.GetValue()
        new_color = [ x / 255.0 for x in new_color ]
        self.map_object.colors[i] = new_color
        self.parent.update_maps = True
    self.parent.OnRedrawGL()

  def OnClose (self, event) :
    self.Destroy()

  def OnDestroy (self, event) :
    self.parent.map_ctrls.pop(self.map_id)

  def refresh_iso_levels (self) :
    levels = self.map_object.iso_levels
    assert len(levels) == len(self.ctrls)
    for i, (s, t, c) in enumerate(self.ctrls) :
      s.SetValue(int(10 * levels[i]))
      t.SetValue("%.2f" % levels[i])

  def refresh_attach_mouse (self) :
    self.attach_mouse_wheel.SetValue(self.map_id==self.parent.selected_map_id)

  def OnAttachMouse (self, event) :
    if self.attach_mouse_wheel.GetValue() == True :
      self.parent.set_selected_map(self.map_id)
    else :
      self.parent.set_selected_map(None)

#---end
