from __future__ import division

# XXX: To keep these classes as clean as possible, selections are handled
# entirely in wx_selection_editor.py.
# TODO: hide nonbonded point for any atom that has an ellipsoid drawn
# TODO: clean up handling of changes in atom count

import iotbx.phil
from cctbx import uctbx
from gltbx.wx_viewer import wxGLWindow
import gltbx.util
from gltbx import viewer_utils, quadrics
from gltbx.gl import *
from gltbx.glu import *
import gltbx
from scitbx.array_family import flex, shared
from scitbx.math import minimum_covering_sphere, sphere_3d
from libtbx.introspection import method_debug_log
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import wx
import sys, os

debug = method_debug_log()

opengl_phil = iotbx.phil.parse("""
opengl {
  line_width = 2
    .type = int
    .style = spinner min:1 max:10
  nonbonded_line_width = 2
    .type = int
    .style = spinner min:1 max:10
  map_radius = 10
    .type = int
    .style = spinner min:1 max:40
  base_atom_color = 1.0 1.0 0.0
    .type = floats
    .style = color
  background_color = 0.0 0.0 0.0
    .type = floats
    .style = color
  default_coloring = *rainbow element b chain
    .type = choice
  default_representation = *all_atoms trace bonded_only
    .type  = choice
  show_hydrogens = False
    .type = bool
  label_clicked_atom = True
    .type = bool
  use_atom_color_for_labels = True
    .type = bool
  orthographic = False
    .type = bool
    .style = noauto
  animate_zoom = False
    .type = bool
}
""")

draw_modes = [ ("trace", "Show trace"),
               ("all_atoms", "Show all atoms"),
               ("bonded_only", "Show bonded atoms"), ]
draw_flags = [ ("flag_show_hydrogens", "Show hydrogens"),
               ("flag_show_ellipsoids", "Show B-factor ellipsoids"),
               ("flag_show_labels", "Show labels"),
               ("flag_show_noncovalent_bonds", "Show non-covalent bonds"), ]
color_modes = [ ("rainbow", "Color rainbow"),
                ("b", "Color by B-factor"),
                ("element", "Color by element"),
                ("chain", "Color by chain"), ]

#-----------------------------------------------------------------------
# XXX: None of the data in this class is used directly in OpenGL calls;
# instead, the model_scene object contains the subset of information
# required for immediate display.
class model_data (object) :
  def __init__ (self, object_id, pdb_hierarchy, atomic_bonds,
      base_color=(0.0,1.0,1.0)) :
    self.object_id = object_id
    self.base_color = base_color
    self.draw_mode = None
    self.current_bonds = None
    self.noncovalent_bonds = None
    self.color_mode = None #"rainbow"
    self.flag_object_visible = True
    self._color_cache = {}
    self.flag_show_hydrogens = False
    self.flag_show_lines = True
    self.flag_show_labels = True
    self.flag_show_points = True
    self.flag_show_spheres = False
    self.flag_show_ellipsoids = False
    self.flag_show_noncovalent_bonds = False
    self.update_structure(pdb_hierarchy, atomic_bonds)
    self.use_u_aniso = flex.bool(self.atoms.size())
    #self.recalculate_visibility()

  def reset (self) :
    self.is_changed = False

  @debug
  def get_scene_data (self) :
    if self.atoms.size() != self.visibility.atoms_visible.size() :
      self.recalculate_visibility()
    return model_scene(bonds=self.current_bonds,
      points=self.atoms.extract_xyz(),
      b_iso=self.atoms.extract_b(),
      b_aniso=self.atoms.extract_uij(),
      atom_colors=self.atom_colors,
      atom_labels=self.atom_labels,
      atom_radii=self.atom_radii,
      visibility=self.visibility,
      noncovalent_bonds=self.noncovalent_bonds)

  def set_noncovalent_bonds (self, bonded_atoms) :
    self.noncovalent_bonds = bonded_atoms

  def update_scene_data (self, scene) :
    scene.update_bonds(self.current_bonds)
    scene.update_colors(self.atom_colors)
    scene.update_visibility(self.visibility)
    scene.clear_lists()

  @debug
  def update_xyz (self, xyz) :
    assert xyz.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.xyz = xyz[i_seq]
    self.is_changed = True

  @debug
  def update_u_iso (self, u_iso) :
    assert u_iso.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.b = adptbx.u_as_b(u_iso[i_seq])
    self.is_changed = True

  @debug
  def update_u_aniso (self, u_aniso, aniso_flag=None) :
    assert u_aniso.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.uij = u_aniso[i_seq]
    self.is_changed = True

  @debug
  def update_from_xray_structure (self, xray_structure) :
    sites_cart = xray_structure.sites_cart()
    u_iso = xray_structure.extract_u_iso_or_u_equiv()
    u_aniso = xray_structure.extract_u_cart_plus_u_iso()
    occ = xray_structure.scatterers().extract_occupancies()
    assert sites_cart.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.xyz = sites_cart[i_seq]
      atom.occ = occ[i_seq]
      atom.b = adptbx.u_as_b(u_iso[i_seq])
      atom.uij = u_aniso[i_seq]
    self.use_u_aniso = xray_structure.use_u_aniso()
    self._color_cache["b"] = None
    self.is_changed = True

  @debug
  def update_structure (self, pdb_hierarchy, atomic_bonds) :
    self.pdb_hierarchy = pdb_hierarchy
    self.atoms = pdb_hierarchy.atoms()
    self.atom_count = self.atoms.size()
    if atomic_bonds is None :
      atomic_bonds = flex.stl_set_unsigned(self.atom_count)
    self.atomic_bonds = atomic_bonds
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    #self.index_atoms()
    atom_index = []
    atom_labels = flex.std_string()
    for atom in self.pdb_hierarchy.atoms_with_labels() :
      atom_index.append(atom)
      atom_labels.append(format_atom_label(atom))
    self.atom_index = atom_index
    self.atom_labels = atom_labels
    self.trace_bonds = extract_trace(pdb_hierarchy) #, self.selection_cache)
    if self.draw_mode is None or self.draw_mode.startswith("trace") :
      self.current_bonds = self.trace_bonds
    else :
      self.current_bonds = self.atomic_bonds
    atom_radii = flex.double(self.atoms.size(), 1.5)
    hydrogen_flag = flex.bool(self.atoms.size(), False)
    for i_seq, atom in enumerate(self.atom_index) :
      if atom.element.strip() in ["H", "D"] :
        atom_radii[i_seq] = 0.75
        hydrogen_flag[i_seq] = True
    self.atom_radii = atom_radii
    self.hydrogen_flag = hydrogen_flag
    self._color_cache = {}
    self.is_changed = True

  @debug
  def recalculate_visibility (self) :
    c = 0
    if self.draw_mode == "spheres" :
      show_points = True
    else :
      show_points = self.flag_show_points
    if self.flag_show_hydrogens :
      atoms_drawable = flex.bool(self.atom_count, True)
    else :
      atoms_drawable = self.hydrogen_flag.__invert__()
      #atoms_drawable = flex.bool([ (atom.element != ' H') for atom in atoms ])
    self.visibility = viewer_utils.atom_visibility(
      bonds             = self.current_bonds,
      atoms_drawable    = atoms_drawable,
      flag_show_points  = show_points
    )
    self.visible_atom_count = self.visibility.visible_atoms_count

  def refresh (self) :
    self.is_changed = True
    self._color_cache = {}
    self.set_draw_mode(self.draw_mode)
    self.is_changed = False

  def toggle_hydrogens (self, show_hydrogens) :
    self.flag_show_hydrogens = show_hydrogens
    self.refresh()

  def toggle_ellipsoids (self, show_ellipsoids) :
    self.flag_show_ellipsoids = show_ellipsoids

  def set_draw_mode (self, draw_mode, color_mode=None) :
    if draw_mode == self.draw_mode and not self.is_changed :
      pass
    else :
      self.draw_mode = draw_mode
      show_points = True
      if draw_mode == "spheres" :
        self.flag_show_spheres = True
      else :
        if draw_mode in ["trace", "trace_and_nb"] :
          self.current_bonds = self.trace_bonds
        else :
          self.current_bonds = self.atomic_bonds
        if draw_mode in ["trace", "bonded_only"] :
          self.flag_show_points = False
        else :
          self.flag_show_points = True
      self.recalculate_visibility()
      if color_mode is not None :
        self.color_mode = color_mode
      self.set_color_mode(self.color_mode) # force re-coloring

  #---------------------------------------------------------------------
  # XXX: COLORING
  #
  def set_base_color (self, color) :
    self.base_color = color

  def set_color_mode (self, color_mode) :
    if color_mode == self.color_mode and not self.is_changed :
      pass
    else :
      self.color_mode = color_mode
      if color_mode == "mono" :
        self.color_mono()
      elif color_mode == "rainbow" :
        self.color_rainbow()
      elif color_mode == "b" :
        self.color_b()
      elif color_mode == "chain" :
        self.color_by_chain()
      elif color_mode == "element" :
        self.color_by_element()

  @debug
  def color_mono (self) :
    cached = self._color_cache.get("mono")
    if cached is not None :
      self.atom_colors = cached
    else :
      self.atom_colors = flex.vec3_double(
        [ self.base_color for i in xrange(0, self.atoms.size()) ]
      )
      self._color_cache["mono"] = self.atom_colors

  @debug
  def color_rainbow (self) :
    cached = self._color_cache.get("rainbow")
    if cached is not None :
      self.atom_colors = cached
    else :
      self.atom_colors = viewer_utils.color_rainbow(
        atoms_visible = self.visibility.atoms_visible,
        visible_atom_count = self.visible_atom_count
      )
      self._color_cache["rainbow"] = self.atom_colors

  @debug
  def color_b (self) :
    cached = self._color_cache.get("b")
    if cached is not None :
      self.atom_colors = cached
    else :
      self.atom_colors = viewer_utils.color_by_property(
        atom_properties       = self.atoms.extract_b(),
        atoms_visible         = self.visibility.atoms_visible,
        color_invisible_atoms = False,
        use_rb_color_gradient = False
      )
      self._color_cache["b"] = self.atom_colors

  @debug
  def color_by_chain (self) :
    cached = self._color_cache.get("chain")
    if cached is not None :
      self.atom_colors = cached
    else :
      c = 0
      for chain in self.pdb_hierarchy.chains() :
        c += 1
      rainbow = viewer_utils.make_rainbow_gradient(c)
      j = 0
      chain_shades = {}
      for chain in self.pdb_hierarchy.chains() :
        chain_shades[chain.id] = rainbow[j]
      atom_colors = flex.vec3_double()
      for atom in self.pdb_hierarchy.atoms_with_labels() :
        atom_colors.append(chain_shades[atom.chain_id])
      self.atom_colors = atom_colors
      self._color_cache["chain"] = atom_colors

  @debug
  def color_by_element (self) :
    cached = self._color_cache.get("element")
    if cached is not None :
      self.atom_colors = cached
    else :
      # these are approximations based on my (probably faulty) memory.
      # feel free to change to something more reasonable.
      element_shades = {' C' : self.base_color, # usually yellow or grey
                        ' H' : (0.95, 0.95, 0.95), # very light grey
                        ' N' : (0.0, 0.0, 1.0),    # blue
                        ' O' : (1.0, 0.0, 0.0),    # red
                        ' S' : (1.0, 0.5, 0.0),    # orange
                        ' P' : (1.0, 1.0, 0.0),    # yellow
                        'Se' : (0.0, 1.0, 0.0),    # green
                        'Mg' : (0.7, 0.7, 0.9),    # very pale blue
                        'Fe' : (0.8, 0.2, 0.0),    # rust
                        'Cl' : (0.8, 1.0, 0.2),    # yellow-green
                        'Na' : (0.7, 0.7, 0.7),    # light grey
                        'Ca' : (1.0, 1.0, 1.0),    # white
                        'Mn' : (1.0, 0.6, 0.8),    # lavender
                        'Zn' : (0.8, 0.9, 1.0),    # very pale cyan
                        'Ni' : (0.0, 0.8, 0.4),    # teal
                        'Cu' : (0.0, 0.8, 0.7),    # blue-green
                        'Co' : (0.0, 0.5, 0.6) }   # marine
      atom_colors = flex.vec3_double()
      for atom in self.pdb_hierarchy.atoms_with_labels() :
        element = atom.element
        if element == "  " :
          first_char = atom.name.strip()[0]
          element = " " + first_char
        color = element_shades.get(element)
        if color is None :
          color = (0.5, 0.5, 0.5)
        atom_colors.append(color)
      self.atom_colors = atom_colors
      self._color_cache["element"] = cached

#-----------------------------------------------------------------------
# XXX: this class contains only the information needed for OpenGL commands,
# which are also implemented as methods here.
class model_scene (object) :
  def __init__ (self, bonds, points, b_iso, b_aniso, atom_colors, atom_labels,
      atom_radii, visibility, noncovalent_bonds) :
    adopt_init_args(self, locals())
    self.clear_lists()
    self.clear_labels()
    self.update_visibility(visibility)

  @debug
  def clear_lists (self) :
    self.points_display_list = None
    self.lines_display_list = None
    self.spheres_display_list = None
    self.ellipsoid_display_list = None
    self.selection_display_list = None
    self.labels_display_list = None
    self.nc_display_list = None

  @debug
  def clear_labels (self) :
    self.show_labels = flex.bool(self.points.size(), False)
    self.labels_display_list = None

  @debug
  def add_label (self, i_seq) :
    self.show_labels[i_seq] = True
    self.labels_display_list = None

  @debug
  def update_colors (self, atom_colors) :
    assert atom_colors.size() == self.points.size()
    self.atom_colors = atom_colors

  @debug
  def update_bonds (self, bonds) :
    assert bonds.size() == self.points.size()
    self.bonds = bonds

  @debug
  def update_visibility (self, visibility) :
    assert visibility.atoms_visible.size() == self.points.size()
    self.atoms_visible = visibility.atoms_visible
    self.bonds_visible = visibility.bonds_visible
    self.points_visible = visibility.points_visible
    self.spheres_visible = self.atoms_visible # XXX: what to do here?
    self.visible_atom_count = visibility.visible_atoms_count
    self.clear_lists()

  def draw_points (self) :
    if self.points_display_list is None :
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      viewer_utils.draw_points(
        points = self.points,
        atom_colors = self.atom_colors,
        points_visible = self.points_visible)
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_lines (self) :
    if self.lines_display_list is None :
      self.lines_display_list = gltbx.gl_managed.display_list()
      self.lines_display_list.compile()
      viewer_utils.draw_bonds(
        points = self.points,
        bonds  = self.bonds,
        atom_colors = self.atom_colors,
        bonds_visible = self.bonds_visible)
      self.lines_display_list.end()
    self.lines_display_list.call()

  def draw_spheres (self, scale_factor=1.0) :
    if self.spheres_display_list is None :
      self.spheres_display_list = gltbx.gl_managed.display_list()
      self.spheres_display_list.compile()
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
      atom_radii = self.atom_radii
      atom_colors = self.atom_colors
      spheres_visible = self.spheres_visible
      for i_seq, point in enumerate(self.points) :
        if spheres_visible[i_seq] :
          glColor3f(*atom_colors[i_seq])
          glPushMatrix()
          glTranslated(*point)
          gltbx.util.SolidSphere(radius=atom_radii[i_seq] * scale_factor,
                                 slices=50, stacks=50)
          glPopMatrix()
      self.spheres_display_list.end()
    self.spheres_display_list.call()

  def draw_ellipsoids (self, proto_ellipsoid) :
    if self.ellipsoid_display_list is None :
      self.ellipsoid_display_list = gltbx.gl_managed.display_list()
      self.ellipsoid_display_list.compile()
      points = self.points
      atoms_visible = self.atoms_visible
      atom_colors = self.atom_colors
      for i_seq, uij in enumerate(self.b_aniso) :
        if atoms_visible[i_seq] and uij[0] != -1 :
          col = list(atom_colors[i_seq]) + [1.0]
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col)
          glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
          glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
          proto_ellipsoid.draw(points[i_seq], uij)
      self.ellipsoid_display_list.end()
    self.ellipsoid_display_list.call()

  def draw_labels (self, font, use_atom_color=False) :
    glDisable(GL_LIGHTING)
    if (self.labels_display_list is None) :
      self.labels_display_list = gltbx.gl_managed.display_list()
      self.labels_display_list.compile()
      points = self.points
      atoms_visible = self.atoms_visible
      atom_colors = self.atom_colors
      atom_labels = self.atom_labels
      for i_seq, show_label in enumerate(self.show_labels) :
        if atoms_visible[i_seq] and show_label :
          if use_atom_color :
            glColor3f(*atom_colors[i_seq])
          glRasterPos3f(*points[i_seq])
          font.render_string(atom_labels[i_seq])
      self.labels_display_list.end()
    self.labels_display_list.call()

  def draw_noncovalent_bonds (self) :
    if self.noncovalent_bonds is None :
      return
    if self.nc_display_list is None :
      self.nc_display_list = gltbx.gl_managed.display_list()
      self.nc_display_list.compile()
      points = self.points
      bonded_atoms = self.noncovalent_bonds
      for i_seq, j_seq in bonded_atoms :
        glBegin(GL_LINES)
        glVertex3f(*points[i_seq])
        glVertex3f(*points[j_seq])
        glEnd()
      self.nc_display_list.end()
    self.nc_display_list.call()

########################################################################
# VIEWER CLASS
#
UPDATE_MODEL_ID = wx.NewId()
ADD_MODEL_ID = wx.NewId()
class AddModelEvent (wx.PyEvent) :
  event_id = ADD_MODEL_ID
  recenter = True
  def __init__ (self, model_id, pdb_hierarchy, atomic_bonds) :
    wx.PyEvent.__init__(self)
    self.data = (model_id, pdb_hierarchy, atomic_bonds)
    self.SetEventType(self.event_id)

class UpdateModelEvent (AddModelEvent) :
  event_id = UPDATE_MODEL_ID
  recenter = False

class model_viewer_mixin (wxGLWindow) :
  def __init__ (self, *args, **kwds) :
    wxGLWindow.__init__(self, *args, **kwds)
    self.Bind(wx.EVT_LEFT_DCLICK, self.OnDoubleClick)
    self.Connect(-1, -1, UPDATE_MODEL_ID, self.OnUpdateModel)
    self.Connect(-1, -1, ADD_MODEL_ID, self.OnAddModel)
    self.minimum_covering_sphere = None
    self.show_object = {}
    self.pick_object = {}
    self.model_objects = []
    self.model_ids = []
    self.scene_objects = {}
    self.model_colors = {}
    self.model_reps = {}
    self.update_scene = False
    self.buffer_factor = 2 # see gltbx.wx_viewer
    self.min_slab = 4
    self.min_viewport_use_fraction = 0.1
    self.min_dist = 4.0
    self.sphere_scale_factor = 1.0
    self.update_settings(opengl_phil.extract())
    self.closest_point_i_seq     = None
    self.closest_point_model_id  = None
    # toggles for viewable objects
    self.flag_show_fog                     = True
    self.flag_show_lines                   = True
    self.flag_show_points                  = True
    self.flag_show_spheres                 = True
    self.flag_use_lights                   = True
    self.flag_show_labels                  = True
    self.flag_show_trace                   = False
    self.flag_show_noncovalent_bonds            = False
    self.flag_show_hydrogens               = False
    self.flag_show_ellipsoids              = True
    self.flag_smooth_lines                 = True
    self.flag_recenter_on_click            = False

  @debug
  def InitGL(self):
    gltbx.util.handle_error()
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    self.minimum_covering_sphere_display_list = None
    glDepthFunc(GL_LESS)
    glEnable(GL_ALPHA_TEST)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    # XXX: line smoothing is pretty essential for wireframe representation;
    # the problem with nvidia cards is really only a problem for the isomesh
    # in wx_map_viewer.py.
    glEnable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    self.initialize_modelview()
    if self.flag_use_lights :
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (0.2, 0.2, 0.2, 1.))
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1, 1, 1, 1.))
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    self.proto_ellipsoid = quadrics.proto_ellipsoid(slices=32, stacks=32)
    gltbx.util.handle_error()

  @debug
  def initialize_modelview (self) :
    if self.minimum_covering_sphere is not None :
      wxGLWindow.initialize_modelview(self)
    else :
      self.setup_lighting()

  def OnRedrawGL (self, event=None) :
    self.check_and_update_model_scenes()
    if self.minimum_covering_sphere is None :
      gltbx.util.handle_error()
      glClear(GL_COLOR_BUFFER_BIT)
      glClear(GL_DEPTH_BUFFER_BIT)
      glFlush()
      self.SwapBuffers()
      gltbx.util.handle_error()
    else :
      wxGLWindow.OnRedrawGL(self, event)

  def get_clipping_distances (self) :
    slab = self.far - self.near
    clip = (1.0 - self.slab_scale) * (slab / 2.0)
    near = self.near + clip
    far = self.far - clip
    if near < self.min_near :
      near = self.min_near
    if near > far or far < (near + self.min_slab) :
      far = near + self.min_slab
    return (near, far)

  def check_and_update_model_scenes (self) :
    if self.update_scene :
      self.update_scene_objects()
      self.update_scene = False

  def DrawGL(self):
    if self.GL_uninitialised or len(self.scene_objects) == 0 :
      return
    if self.flag_show_points :
      self.draw_points()
    if self.flag_show_lines :
      self.draw_lines()
    if self.flag_show_spheres :
      self.draw_spheres()
    if self.flag_show_ellipsoids :
      self.draw_ellipsoids()
    if self.flag_show_labels :
      self.draw_labels()
    if self.flag_show_noncovalent_bonds :
      self.draw_noncovalent_bonds()

  def draw_points (self) :
    glDisable(GL_LIGHTING)
    glLineWidth(self.settings.opengl.nonbonded_line_width)
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_points :
        self.scene_objects[model_id].draw_points()

  def draw_spheres (self) :
    glMatrixMode(GL_MODELVIEW)
    if self.flag_use_lights :
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glEnable(GL_LIGHT1)
      glEnable(GL_NORMALIZE)
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [1.0,1.0,1.0,1.0])
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0.1, 0.1, 0.1, 1.0])
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_spheres :
        self.scene_objects[model_id].draw_spheres(self.sphere_scale_factor)

  def draw_lines (self) :
    glEnable(GL_LINE_SMOOTH)
    #glDisable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    glDisable(GL_LIGHTING)
    glLineWidth(self.settings.opengl.line_width)
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_lines :
        self.scene_objects[model_id].draw_lines()
    if not self.flag_smooth_lines :
      glDisable(GL_LINE_SMOOTH)

  def draw_ellipsoids (self) :
    glMatrixMode(GL_MODELVIEW)
    if self.flag_use_lights :
      glShadeModel(GL_SMOOTH)
      glEnable(GL_DEPTH_TEST)
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glEnable(GL_LIGHT1)
      glEnable(GL_NORMALIZE)
      glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
      glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.5, 0.5, 0.5, 1.0])
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
    proto_ellipsoid = self.proto_ellipsoid
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_ellipsoids :
        self.scene_objects[model_id].draw_ellipsoids(proto_ellipsoid)

  def draw_labels (self) :
    glDisable(GL_LIGHTING)
    use_atom_color = self.settings.opengl.use_atom_color_for_labels
    if not use_atom_color :
      glColor3f(1.0, 1.0, 1.0)
    font = gltbx.fonts.ucs_bitmap_8x13
    font.setup_call_lists()
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_labels :
        self.scene_objects[model_id].draw_labels(font, use_atom_color)

  def draw_noncovalent_bonds (self) :
    glDisable(GL_LIGHTING)
    glColor3f(0.0, 1.0, 0.0)
    glLineStipple(4, 0xAAAA)
    glEnable(GL_LINE_STIPPLE)
    for model_id, model in self.iter_models() :
      if self.show_object[model_id] and model.flag_show_noncovalent_bonds :
        self.scene_objects[model_id].draw_noncovalent_bonds()
    glDisable(GL_LINE_STIPPLE)

  @debug
  def refresh_bg_color (self) :
    (r, g, b) = tuple(self.settings.opengl.background_color)
    glClearColor(r, g, b, 0.0)

  #---------------------------------------------------------------------
  # Non-OpenGL code below here
  #
  def update_settings (self, params, redraw=False) :
    assert (len(params.opengl.base_atom_color) ==
            len(params.opengl.background_color) == 3)
    self.settings = params
    if params.opengl.animate_zoom :
      self.animation_time = 1
    else :
      self.animation_time = 0
    self.toggle_hydrogens(params.opengl.show_hydrogens)
    if redraw :
      self.update_scene = True

  def iter_models (self) :
    for model_id, model_object in zip(self.model_ids, self.model_objects) :
      yield (model_id, model_object)

  def get_model (self, object_id) :
    for model in self.model_objects :
      if model.object_id == object_id :
        return model
    return None

  def show_model_controls_menu (self, object_id, source_widget) :
    model = self.get_model(object_id)
    if model is None :
      return
    menu = wx.Menu(title=object_id)
    for mode_name, mode_label in draw_modes :
      item = menu.AppendCheckItem(-1, mode_label)
      if model.draw_mode == mode_name :
        item.Check(True)
      source_widget.Bind(wx.EVT_MENU, self.OnModelMenu, item)
    menu.AppendSeparator()
    for flag_name, flag_label in draw_flags :
      item = menu.AppendCheckItem(-1, flag_label)
      if getattr(model, flag_name, False) :
        item.Check(True)
      source_widget.Bind(wx.EVT_MENU, self.OnModelMenu, item)
    menu.AppendSeparator()
    for color_name, color_label in color_modes :
      item = menu.AppendCheckItem(-1, color_label)
      if model.color_mode == color_name :
        item.Check(True)
      source_widget.Bind(wx.EVT_MENU, self.OnModelMenu, item)
    menu.AppendSeparator()
    item = menu.Append(-1, "Change base color. . .")
    source_widget.Bind(wx.EVT_MENU, self.OnChangeModelColor, item)
    source_widget.PopupMenu(menu)
    menu.Destroy()
    self.OnRedrawGL()

  def process_model_menu_event (self, object_id, menu_item) :
    model = self.get_model(object_id)
    if model is None :
      return
    item_label = menu_item.GetText()
    is_checked = menu_item.IsChecked()
    for mode, label in draw_modes :
      if label == item_label and is_checked :
        model.set_draw_mode(mode)
        self.update_scene = True
        return True
    for mode, label in color_modes :
      if label == item_label and is_checked :
        model.set_color_mode(mode)
        self.update_scene = True
        return True
    for flag, label in draw_flags :
      if label == item_label :
        setattr(model, flag, is_checked)
        model.refresh()
        self.update_scene = True
        return True
    return False

  def set_model_state (self, object_id, model_state) :
    model = self.get_model(object_id)
    for name in model_state :
      setattr(model, name, model_state[name])
    model.set_draw_mode(model_state["draw_mode"])
    #model.recalculate_visibility()
    self.update_scene = True

  @debug
  def add_model (self, model_id, pdb_hierarchy, atomic_bonds=None) :
    assert isinstance(model_id, str)
    model = model_data(model_id, pdb_hierarchy, atomic_bonds,
      base_color=self.settings.opengl.base_atom_color)
    self._add_model(model_id, model)

  def delete_model (self, model_id) :
    if model_id in self.model_ids :
      i = self.model_ids.index(model_id)
      self.model_ids.pop(i)
      self.model_objects.pop(i)
      if model_id in self.scene_objects :
        self.scene_objects.pop(model_id)
    self.update_scene = True

  def _add_model (self, model_id, model) :
    model.set_draw_mode(draw_mode=self.settings.opengl.default_representation,
      color_mode=self.settings.opengl.default_coloring)
    self.model_ids.append(model_id)
    self.model_objects.append(model)
    self.show_object[model_id] = True
    self.pick_object[model_id] = True
    self.update_scene = True

  def update_model (self, model_id, pdb_hierarchy, atomic_bonds) :
    model = self.get_model(model_id)
    if model is not None :
      model.update_structure(pdb_hierarchy, atomic_bonds)
      model.set_draw_mode(model.draw_mode)
      if model_id in self.scene_objects :
        self.scene_objects.pop(model_id)
      self.update_scene = True
    else :
      self.add_model(model_id, pdb_hierarchy, atomic_bonds)

  def update_model_from_xray_structure (self, model_id, xray_structure) :
    model = self.get_model(model_id)
    if model is not None :
      model.update_from_xray_structure(xray_structure)
      self.update_scene = True

  def set_noncovalent_bonds (self, model_id, bonded_atoms) :
    model = self.get_model(model_id)
    if model is not None :
      model.set_noncovalent_bonds(bonded_atoms)
      #model.flag_show_noncovalent_bonds = True

  def update_mcs (self, points, recenter_and_zoom=True, buffer=0) :
    mcs = minimum_covering_sphere(points=points,
                                  epsilon=0.1)
    if buffer > 0 :
      self.minimum_covering_sphere = sphere_3d(
        center=mcs.center(),
        radius=mcs.radius() + buffer)
    else :
      self.minimum_covering_sphere = mcs
    if recenter_and_zoom :
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  def zoom_object (self, object_id) :
    self.update_scene_objects()
    assert object_id in self.scene_objects
    self.update_mcs(self.scene_objects[object_id].points)

  @debug
  def unzoom (self, event=None) :
    self.update_scene_objects()
    if len(self.scene_objects) > 0 :
      points = flex.vec3_double()
      for object_id, scene in self.scene_objects.iteritems() :
        points.extend(scene.points)
      self.update_mcs(points)

  @debug
  def recenter_on_atom (self, object_id, i_seq) :
    assert object_id is not None and i_seq >= 0
    scene = self.scene_objects.get(object_id)
    if scene is not None and i_seq < scene.points.size() :
      self.rotation_center = scene.points[i_seq]
      self.move_to_center_of_viewport(self.rotation_center)

  @debug
  def process_pick_points (self) :
    self.closest_point_object_id = None
    self.closest_point_i_seq = None
    if self.pick_points is not None and len(self.scene_objects) > 0 :
      for object_id in self.model_ids :
        if self.show_object[object_id] and self.pick_object[object_id] :
          scene = self.scene_objects.get(object_id)
          if scene is None :
            continue
          closest_point_i_seq = viewer_utils.closest_visible_point(
            points = scene.points,
            atoms_visible = scene.atoms_visible,
            point0 = self.pick_points[0],
            point1 = self.pick_points[1]
          )
          if closest_point_i_seq is not None :
            self.closest_point_i_seq = closest_point_i_seq
            self.closest_point_object_id = object_id
            break
    if self.closest_point_i_seq is not None :
      clicked_scene = self.scene_objects.get(self.closest_point_object_id)
      if self.settings.opengl.label_clicked_atom :
        clicked_scene.add_label(self.closest_point_i_seq)
      if self.flag_recenter_on_click :
        self.recenter_on_atom(self.closest_point_object_id,
          self.closest_point_i_seq)

  @debug
  def update_scene_objects (self) :
    points = flex.vec3_double()
    for object_id, model in self.iter_models() :
      current_scene = self.scene_objects.get(object_id)
      if current_scene is None or model.is_changed :
        current_scene = model.get_scene_data()
        self.scene_objects[object_id] = current_scene
        model.reset()
      else :
        model.update_scene_data(current_scene)
      points.extend(current_scene.points)
    if points.size() == 0 :
      points.append((0,0,0))
    self.update_mcs(points, recenter_and_zoom=False)

  @debug
  def process_key_stroke (self, key) :
    if key == ord('u') :
      self.unzoom()
    elif key == ord('h') :
      self.flag_show_hydrogens = not self.flag_show_hydrogens
      self.toggle_hydrogens(self.flag_show_hydrogens)
      self.update_scene = True
    elif key == ord('e') :
      self.flag_show_ellipsoids = not self.flag_show_ellipsoids
      self.toggle_ellipsoids(self.flag_show_ellipsoids)
      self.update_scene = True
    elif key == ord('l') :
      self.flag_show_labels = not self.flag_show_labels
      self.update_scene = True
    elif key == ord('t') :
      self.flag_show_trace = not self.flag_show_trace
      if self.flag_show_trace :
        self.set_draw_mode("trace")
      else :
        self.set_draw_mode("all_atoms")
      self.update_scene = True
    elif key == ord('b') :
      self.set_color_mode("b")
    elif key == ord('r') :
      self.set_color_mode("rainbow")
    elif key == ord('y') :
      self.set_color_mode("element")
    elif key == 8 : # delete, at least on Mac
      self.clear_all_labels()
      self.update_scene = True
    elif key == ord('q') :
      app = wx.GetApp()
      app.Exit()
    else :
      pass
      #print key
    if self.update_scene :
      self.OnRedrawGL()

  def toggle_visibility (self, show_object, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        self.show_object[model_id] = show_object
    self.update_scene = True

  def hide_models (self, object_id=None) :
    for model_id in self.model_ids :
      if object_id is None or model_id == object_id :
        self.show_object[model_id] = False

  def show_all (self) :
    for model_id in self.model_ids :
      self.show_object[model_id] = True

  @debug
  def set_draw_mode (self, draw_mode, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        model.set_draw_mode(draw_mode)
    self.update_scene = True

  @debug
  def set_color_mode (self, color_mode, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        model.set_color_mode(color_mode)
    self.update_scene = True

  def set_model_base_color (self, color, object_id=None) :
    assert len(color) == 3
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        model.set_base_color(color)
    self.update_scene = True

  def toggle_ellipsoids (self, show_ellipsoids, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        model.toggle_ellipsoids(show_ellipsoids)

  def toggle_hydrogens (self, show_hydrogens, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        model.toggle_hydrogens(show_hydrogens)

  # TODO: something smarter - temporary toggle for draw_mode?
  def toggle_trace (self, show_trace, object_id=None) :
    for model_id, model in self.iter_models() :
      if object_id is None or object_id == model_id :
        if show_trace :
          model.set_draw_mode("trace")
        else :
          model.set_draw_mode("all_atoms")
    self.update_scene = True

  def toggle_labels (self, show_labels) :
    self.flag_show_labels = show_labels

  def clear_all_labels (self) :
    for object_id, scene in self.scene_objects.iteritems() :
      scene.clear_labels()

  #---------------------------------------------------------------------
  # EVENTS
  @debug
  def OnUpdate (self, event) :
    self.update_scene_objects()
    if getattr(event, "recenter", False) :
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  def OnUpdateModel (self, event) :
    (model_id, pdb_hierarchy, atomic_bonds) = event.data
    self.update_model(model_id, pdb_hierarchy, atomic_bonds)
    self.OnRedrawGL()

  def OnAddModel (self, event) :
    (model_id, pdb_hierarchy, atomic_bonds) = event.data
    self.add_model(model_id, pdb_hierarchy, atomic_bonds)
    self.OnRedrawGL()

  def OnDoubleClick (self, event) :
    self.get_pick_points((event.GetX(), event.GetY()))
    self.process_pick_points()
    if self.closest_point_i_seq is not None :
      self.recenter_on_atom(self.closest_point_object_id,
        self.closest_point_i_seq)

  def OnMouseWheel (self, event) :
    scale = event.GetWheelRotation()
    if event.ShiftDown() :
      self.fog_end_offset -= scale
    else :
      self.slab_scale += 0.01 * scale
      if self.slab_scale > 1.0 :
        self.slab_scale = 1.0
      elif self.slab_scale < 0.01 :
        self.slab_scale = 0.01
    self.OnRedrawGL()

  def OnModelMenu (self, event) :
    menu = event.GetEventObject()
    item = menu.FindItemById(event.GetId())
    model_id = menu.GetTitle()
    self.process_model_menu_event(model_id, item)

  def OnChangeModelColor (self, event) :
    menu = event.GetEventObject()
    model_id = menu.GetTitle()
    model = self.get_model(model_id)
    if model is not None :
      base_color = [ int(x*255) for x in model.base_color ]
      new_color = wx.GetColourFromUser(self, base_color)
      new_color = [ x / 255 for x in new_color ]
      model.set_base_color(new_color)
      model.refresh()
      self.update_scene = True

  def thread_safe_add_model (self, model_id, pdb_hierarchy, atomic_bonds) :
    wx.PostEvent(self, AddModelEvent(model_id, pdb_hierarchy, atomic_bonds))

  def thead_safe_update_model (self, model_id, pdb_hierarchy, atomic_bonds) :
    wx.PostEvent(self, UpdateModelEvent(model_id, pdb_hierarchy, atomic_bonds))

#-----------------------------------------------------------------------
# Utility functions
def extract_trace (pdb_hierarchy, selection_cache=None) :
  if selection_cache is None :
    selection_cache = pdb_hierarchy.atom_selection_cache()
  last_atom     = None
  isel = selection_cache.iselection
  selection_i_seqs = list(
    isel("(name ' CA ' or name ' P  ') and (altloc 'A' or altloc ' ')"))
  last_atom = None
  bonds = shared.stl_set_unsigned(pdb_hierarchy.atoms().size())
  for atom in pdb_hierarchy.atoms_with_labels() :
    if atom.i_seq in selection_i_seqs :
      if last_atom is not None :
        if (atom.chain_id        == last_atom.chain_id and
            atom.model_id        == last_atom.model_id and
            atom.resseq_as_int() == (last_atom.resseq_as_int() + 1) and
            ((atom.altloc == last_atom.altloc) or
             (atom.altloc == "A" and last_atom.altloc == "") or
             (atom.altloc == ""  and last_atom.altloc == "A"))) :
          bonds[last_atom.i_seq].append(atom.i_seq)
          bonds[atom.i_seq].append(last_atom.i_seq)
      last_atom = atom
  return bonds

def format_atom_label (atom_info) :
  return ("%s %s%s %s %s" % (atom_info.name, atom_info.altloc,
        atom_info.resname, atom_info.chain_id, atom_info.resid())).strip()
