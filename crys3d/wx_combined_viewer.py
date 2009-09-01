from __future__ import division

import sys, os
from crys3d import wx_viewer_zoom
from crys3d.wx_selection_editor import selection_editor_mixin
import iotbx.phil
import gltbx.util
from gltbx.gl import *
from gltbx.glu import *
from cctbx import maptbx
from scitbx.math import minimum_covering_sphere
from scitbx.array_family import flex
from scitbx import iso_surface
from libtbx import adopt_init_args
import wx

viewer_phil = iotbx.phil.parse("""
  include scope crys3d.wx_selection_editor.viewer_phil
""", process_includes=True)

class map_data (object) :
  def __init__ (self, map) :
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
    self.radius = 10.0

  def set_iso_levels (self, levels) :
    self.iso_levels = levels
    if len(self.colors) != len(levels) :
      self.colors = [ (1.0,1.0,1.0) for x in levels ]

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

  def draw_mesh (self) :
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glMultTransposeMatrixd(self.orthogonaliser)
    gltbx.util.handle_error()
    for i, triangulation in enumerate(self.triangles) :
      glLineWidth(0.2)
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
    glPopMatrix()


class map_viewer_mixin (wx_viewer_zoom.viewer_with_automatic_zoom) :
  initialize_map_viewer_super = True
  def __init__ (self, *args, **kwds) :
    if self.initialize_map_viewer_super :
      wx_viewer_zoom.viewer_with_automatic_zoom.__init__(self, *args, **kwds)
    # various data objects
    self.map_ids     = []
    self.map_objects = []
    self.map_scenes  = {}
    self.show_object = {}
    # user settings
    self.mesh_line_width = 0.25 # very buggy on OS X + NVidia (and ???)
    self.buffer_factor = 2
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
    wx_viewer_zoom.viewer_with_automatic_zoom.OnRedrawGL(self, event)

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
    wx_viewer_zoom.viewer_with_automatic_zoom.OnTranslate(self, event)
    self.update_map_scenes()

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

  def OnMouseWheel (self, event) :
    scale = event.GetWheelRotation()
    if event.ShiftDown() :
      self.fog_end_offset -= scale
    else :
      self.clip_far -= scale
    self.OnRedrawGL()

  def process_key_stroke (self, key) :
    pass

  def add_map (self, map_id, map) :
    assert not map_id in self.map_ids
    map_object = map_data(map)
    self.map_ids.append(map_id)
    self.map_objects.append(map_object)
    self.show_object[map_id] = True
    self.update_maps = True

  def delete_map (self, map_id) :
    if map_id in self.map_ids :
      i = self.map_ids.index(map_id)
      self.map_ids.pop(i)
      self.map_objects.pop(i)
      if map_id in self.scene_objects :
        self.scene_objects.pop(map_id)
    self.update_maps = True

  def update_map (self, map_id, map) :
    if not map_id in self.map_ids :
      self.add_map(map_id, map)
    else :
      map_object = self.get_map(map_id)
      map_object.update_map_data(map)
      self.update_maps = True

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

  def update_map_scenes (self) :
    for object_id, map in self.iter_maps() :
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

class model_and_map_viewer (selection_editor_mixin, map_viewer_mixin) :
  initialize_map_viewer_super = False
  def __init__ (self, *args, **kwds) :
    selection_editor_mixin.__init__(self, *args, **kwds)
    map_viewer_mixin.__init__(self, *args, **kwds)

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
    self.check_and_update_map_scenes()
    selection_editor_mixin.OnRedrawGL(self, event)

  def DrawGL (self) :
    selection_editor_mixin.DrawGL(self)
    map_viewer_mixin.DrawGL(self)

  def process_key_stroke (self, key) :
    selection_editor_mixin.process_key_stroke(self, key)
    map_viewer_mixin.process_key_stroke(self, key)
    if key == wx.WXK_UP :
      self.fog_start_offset += 1
    elif key == wx.WXK_DOWN :
      self.fog_start_offset -= 1
    elif key == wx.WXK_LEFT :
      self.clip_near -= 1
    elif key == wx.WXK_RIGHT :
      self.clip_near += 1
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

#---end
