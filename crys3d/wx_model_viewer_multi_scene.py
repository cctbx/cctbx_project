
from string import strip
from libtbx.utils import Sorry
import gltbx.util
from gltbx import wx_viewer, viewer_utils, quadrics
from gltbx.gl import *
from gltbx.glu import *
from scitbx.array_family import flex, shared
from scitbx.math import minimum_covering_sphere
from mmtbx.monomer_library import pdb_interpretation
from cctbx import uctbx
from libtbx.introspection import method_debug_log
from libtbx import adopt_init_args
import wx
import sys

debug = method_debug_log()

#-----------------------------------------------------------------------
# XXX: None of the data in this class is used directly in OpenGL calls;
# instead, the model_scene object contains the subset of information
# required for immediate display.
class model_data (object) :
  def __init__ (self, object_id, pdb_hierarchy, atomic_bonds,
      initial_color=(0.0,1.0,1.0)) :
    self.object_id = object_id
    self.update_structure(pdb_hierarchy, atomic_bonds)
    self.flag_object_visible = True
    self.flag_allow_selection = True
    self.use_u_aniso = flex.bool(self.atoms.size())
    self.trace_bonds = shared.stl_set_unsigned()
    self.current_bonds = self.atomic_bonds
    self.atom_colors = flex.vec3_double(self.atoms.size(), initial_color)
    self._color_cache = {}
    self.recalculate_visibility(False, True)

    # XXX: selections
    self.selection_string = "None"
    self.selection_covering_sphere = None
    self.selection_i_seqs = None
    self.selection_cache = None
    self.selected_points = []
    self.selection_display_list = None
    self.selection_colors = flex.vec3_double()
    self.atoms_selected = None

  def reset (self) :
    self.is_changed = False

  @debug
  def get_scene_data (self) :
    return model_scene(self.current_bonds, self.atoms.extract_xyz(),
      self.atoms.extract_b(), self.atoms.extract_uij(), self.atom_colors,
      self.visibility)

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
    occ = xray_structure.scatterrers().extract_occupancies()
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
    self.atomic_bonds = atomic_bonds
    self.atoms = pdb_hierarchy.atoms()
    self.atom_count = self.atoms.size()
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    atom_index = []
    for atom in self.pdb_hierarchy.atoms_with_labels() :
      atom_index.append(atom)
    self.atom_index = atom_index
    self.extract_trace()
    assert len(self.atom_index) == self.atom_count
    atom_radii = flex.double(self.atoms.size(), 1.5)
    hydrogen_flag = flex.bool(self.atoms.size(), False)
    for i_seq, atom in enumerate(self.atom_index) :
      if atom.element == ' H' :
        atom_radii[i_seq] = 0.75
        hydrogen_flag[i_seq] = True
    self.atom_radii = atom_radii
    self.hydrogen_flag = hydrogen_flag
    self.minimum_covering_sphere = minimum_covering_sphere(
      points=self.atoms.extract_xyz(),
      epsilon=0.1)
    self._color_cache = {}
    self.is_changed = True

  @debug
  def extract_trace (self) :
    if self.selection_cache is None : return
    last_atom     = None
    isel = self.selection_cache.iselection
    selection_i_seqs = list(
      isel("(name ' CA ' or name ' P  ') and (altloc 'A' or altloc ' ')"))
    last_atom = None
    bonds = shared.stl_set_unsigned(self.atoms.size())
    for atom in self.atom_index :
      if atom.i_seq in selection_i_seqs :
        if last_atom is not None :
          if atom.chain_id        == last_atom.chain_id and \
             atom.model_id        == last_atom.model_id and \
             atom.resseq_as_int() == (last_atom.resseq_as_int() + 1) and \
             compare_conformers(atom.altloc, last_atom.altloc) == True :
            bonds[last_atom.i_seq].append(atom.i_seq)
            bonds[atom.i_seq].append(last_atom.i_seq)
        last_atom = atom
    self.trace_bonds = bonds

  @debug
  def toggle_trace_mode (self, show_trace=True) :
    if show_trace :
      self.current_bonds = self.atomic_bonds
    else :
      self.current_bonds = self.trace_bonds

  @debug
  def recalculate_visibility (self, show_hydrogens=True, show_points=False) :
    c = 0
    atoms = self.atom_index
    if show_hydrogens :
      atoms_drawable = flex.bool(self.atom_count, True)
    else :
      atoms_drawable = flex.bool([ (atom.element != ' H') for atom in atoms ])
    self.visibility = viewer_utils.atom_visibility(
      bonds             = self.current_bonds,
      atoms_drawable    = atoms_drawable,
      flag_show_points  = show_points
    )
    self.visible_atom_count = self.visibility.visible_atoms_count

  #---------------------------------------------------------------------
  # XXX: COLORING
  #
  @debug
  def color_mono (self, base_atom_color=(0.0,1.0,1.0)) :
    cached = self._color_cache.get("mono")
    if cached is not None :
      self.atom_colors = cached
    else :
      self.atom_colors = flex.vec3_double(
        [ self.base_atom_color for i in xrange(0, self.points.size()) ]
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
  def color_b (self, scale_b_to_visible=True) :
    cached = self._color_cache.get("b")
    if cached is not None :
      self.atom_colors = cached
    else :
      self.atom_colors = viewer_utils.color_by_property(
        atom_properties       = self.atoms.extract_b(),
        atoms_visible         = self.visibility.atoms_visible,
        color_invisible_atoms = not scale_b_to_visible,
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
      for i_seq, atom_object in enumerate(self.atom_index) :
        atom_colors.append(chain_shades[atom_object.chain_id])
      self.atom_colors = atom_colors
      self._color_cache["chain"] = atom_colors

  @debug
  def color_by_element (self, carbon_atom_color=(1.0,1.0,0.0)) :
    cached = self._color_cache.get("element")
    if cached is not None :
      self.atom_colors = cached
    else :
      # these are approximations based on my (probably faulty) memory.
      # feel free to change to something more reasonable.
      element_shades = {' C' : carbon_atom_color, # usually yellow or grey
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
      for i_seq, atom_object in enumerate(self.atom_index) :
        element = atom_object.element
        color = element_shades.get(element)
        if color is None :
          color = (0.0, 1.0, 1.0)
        atom_colors.append(color)
      self.atom_colors = atom_colors
      self._color_cache["element"] = cached

#-----------------------------------------------------------------------
# XXX: this class contains only the information needed for OpenGL commands,
# which are also implemented as methods here.
class model_scene (object) :
  def __init__ (self, bonds, points, b_iso, b_aniso, atom_colors, visibility) :
    adopt_init_args(self, locals())
    self.clear_lists()
    self.update_visibility(visibility)
    self.minimum_covering_sphere = minimum_covering_sphere(points=points,
                                                           epsilon=0.1)
  @debug
  def clear_lists (self) :
    self.points_display_list = None
    self.lines_display_list = None
    self.spheres_display_list = None
    self.ellipsoid_display_list = None
    self.selection_display_list = None

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

  def draw_spheres (self) :
    if self.spheres_display_list is None :
      self.spheres_display_list = gltbx.gl_managed.display_list()
      self.spheres_display_list.compile()
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
      atom_radii = self.atom_radii
      spheres_visible = self.spheres_visible
      for i_seq, point in enumerate(points) :
        if self.spheres_visible[i_seq] :
          #glColor3f(*atom_colors[i_seq])
          glPushMatrix()
          glTranslated(*point)
          gltbx.util.SolidSphere(radius=atom_radii[i_seq],
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

########################################################################
# BASE CLASS FOR DISPLAYING STRUCTURES
#
class model_viewer_mixin (wx_viewer.wxGLWindow) :
  initialize_model_viewer_super = True
  def __init__ (self, *args, **kwds) :
    if self.initialize_model_viewer_super :
      wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    self.minimum_covering_sphere = None
    self.show_object = {}
    self.model_objects = []
    self.model_ids = []
    self.scene_objects = {}
    self.update_scene = False
    # various settings; override these in subclasses and/or provide
    # a mechanism for changing their values
    self.buffer_factor           = 2
    self.line_width              = 2
    self.nonbonded_line_width    = 2
    self.base_atom_color         = (0.6, 0.6, 0.6) # grey
    self.orthographic            = False
    self.draw_mode               = "all_atoms"
    self.color_mode              = "b_factors"
    self.recolor                 = self.color_b
    self.scale_b_to_visible      = True
    self.carbon_atom_color       = (1.0, 1.0, 0.0) # yellow
    self.closest_point_i_seq     = None
    self.closest_point_model_id  = None
    # toggles for viewable objects
    self.flag_show_fog                     = True
    self.flag_show_lines                   = True
    self.flag_show_points                  = True
    self.flag_show_spheres                 = False
    self.flag_show_ellipsoids              = False
    self.flag_use_lights                   = True
    self.flag_show_trace                   = False
    self.flag_show_hydrogens               = False
    self.flag_show_ellipsoids              = False

  def iter_models (self) :
    for model_id, model_object in zip(self.model_ids, self.model_objects) :
      yield (model_id, model_object)

  @debug
  def InitGL(self):
    gltbx.util.handle_error()
    glClearColor(self.r_back, self.g_back, self.b_back, 0.0)
    self.minimum_covering_sphere_display_list = None
    glDepthFunc(GL_LESS)
    glEnable(GL_ALPHA_TEST)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_BLEND)
    vendor = glGetString(GL_VENDOR)
    if sys.platform == "darwin" and vendor.startswith("NVIDIA") :
      glDisable(GL_LINE_SMOOTH)
    else :
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
      wx_viewer.wxGLWindow.initialize_modelview(self)
    else :
      self.setup_lighting()

  def DrawGL(self):
    if self.GL_uninitialised or len(self.scene_objects) == 0 :
      return
    if (self.flag_show_points):
      self.draw_points()
    if (self.flag_show_lines):
      self.draw_lines()
    if (self.flag_show_spheres):
      self.draw_spheres()
    if self.flag_show_ellipsoids :
      self.draw_ellipsoids()

  def OnRedrawGL (self, event=None) :
    if self.update_scene :
      self.update_scene_objects()
      self.update_scene = False
    if self.minimum_covering_sphere is None :
      gltbx.util.handle_error()
      glClear(GL_COLOR_BUFFER_BIT)
      glClear(GL_DEPTH_BUFFER_BIT)
      glFlush()
      self.SwapBuffers()
      gltbx.util.handle_error()
    else :
      wx_viewer.wxGLWindow.OnRedrawGL(self, event)

  @debug
  def add_model (self, model_id, pdb_hierarchy, atomic_bonds=None) :
    model = model_data(model_id, pdb_hierarchy, atomic_bonds)
    self.model_ids.append(model_id)
    self.model_objects.append(model)
    self.show_object[model_id] = True
    self.recolor()
    self.update_scene = True

  @debug
  def unzoom (self, event=None) :
    if len(self.scene_objects) > 0 :
      points = flex.vec3_double()
      for scene in self.scene_objects :
        points.extend(scene.points)
      self.minimum_covering_sphere = minimum_covering_sphere(
                                       points=points,
                                       epsilon=0.1)
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  @debug
  def recenter_on_atom (self, object_id, i_seq) :
    assert object_id is not None and i_seq >= 0
    scene = self.scene_objects.get(object_id)
    if scene is not None and i_seq < scene.points.size() :
      self.rotation_center = scene.points[i_seq]
      self.move_to_center_of_viewport(self.rotation_center)

  @debug
  def process_key_stroke (self, key) :
    if key == ord('u') :
      self.unzoom()
    elif key == ord('h') :
      self.flag_show_hydrogens = not self.flag_show_hydrogens
      self.update_scene = True
    elif key == ord('e') :
      self.flag_show_ellipsoids = not self.flag_show_ellipsoids
      self.update_scene = True
    elif key == ord('p') :
      self.flag_show_points = not self.flag_show_points
      self.update_scene = True
    elif key == ord('b') :
      if self.color_mode == "b_factors" :
        self.set_color_mode("element")
      else :
        self.set_color_mode("b_factors")
      self.update_scene = True
    if self.update_scene :
      self.OnRedrawGL()

  @debug
  def process_pick_points (self) :
    if self.pick_points is not None and len(self.scene_objects) > 0 :
      for object_id in self.model_ids :
        if self.show_object[object_id] :
          scene = self.scene_objects.get(object_id)
          if scene is None :
            continue
          self.closest_point_object_id = object_id
          self.closest_point_i_seq = viewer_utils.closest_visible_point(
            points = scene.points,
            atoms_visible = scene.atoms_visible,
            point0 = self.pick_points[0],
            point1 = self.pick_points[1]
          )
          break

  @debug
  def OnUpdate (self, event) :
    self.update_scene_objects()
    if (event is not None and hasattr(event, "recenter") and
        event.recenter == True) or recenter == True :
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

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
    self.minimum_covering_sphere = minimum_covering_sphere(points=points,
                                                           epsilon=0.1)

  @debug
  def update_model_objects (self) :
    for object_id, model in self.iter_models() :
      model.recalculate_visibility(self.flag_show_hydrogens,
        self.flag_show_points)
    self.recolor()

  @debug
  def set_draw_mode (self, draw_mode, redraw=False) :
    if not self.minimum_covering_sphere :
      return
    self.draw_mode = draw_mode
    if draw_mode == "trace" :
      self.flag_show_points = False
      for object in self.model_objects :
        object.toggle_trace_mode(True)
    elif self.draw_mode == "bonded_only" :
      self.flag_show_points = False
      for object in self.model_objects :
        object.toggle_trace_mode(False)
    else :
      self.flag_show_points = True
      for object in self.model_objects :
        object.toggle_trace_mode(False)
    self.update_model_objects()
    self.set_color_mode(self.color_mode, redraw=redraw)

  #---------------------------------------------------------------------
  # coloring
  @debug
  def set_bg_color (self) :
    (r,g,b) = self.bg_color
    glClearColor(r, g, b, 0.0)

  @debug
  def set_color_mode (self, color_mode) :
    self.color_mode = color_mode
    if color_mode == "rainbow" :
      self.recolor = self.color_rainbow
    elif color_mode == "element" :
      self.recolor = self.color_by_element
    elif color_mode == "b_factors" :
      self.recolor = self.color_b
    elif color_mode == "chain" :
      self.recolor = self.color_by_chain
    elif color_mode == "single_color" :
      self.recolor = self.color_mono
    else :
      pass

  def color_mono (self) :
    for object in self.model_objects :
      object.color_mono(self.base_atom_color)

  def color_rainbow (self) :
    for object in self.model_objects :
      object.color_rainbow()

  def color_b (self) :
    for object in self.model_objects :
      object.color_b(self.scale_b_to_visible)

  def color_by_chain (self) :
    for object in self.model_objects :
      object.color_by_chain()

  def color_by_element (self) :
    for object in self.model_objects :
      object.color_by_element(self.carbon_atom_color)

  #---------------------------------------------------------------------
  # DRAWING ROUTINES
  def draw_points (self) :
    glDisable(GL_LIGHTING)
    glLineWidth(self.nonbonded_line_width)
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_points()

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
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_spheres()

  def draw_lines (self) :
    glDisable(GL_LIGHTING)
    glLineWidth(self.line_width)
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_lines()

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
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_ellipsoids(proto_ellipsoid)

########################################################################
# UTILITY FUNCTIONS
#
def compare_conformers (altloc1, altloc2) :
  if altloc1 == altloc2 :
    return True
  elif altloc1 == "A" and altloc2 == "" :
    return True
  elif altloc2 == "A" and altloc1 == "" :
    return True
  else :
    return False

########################################################################
# CLASSES AND METHODS FOR STANDALONE VIEWER
#
class App (wx.App) :
  def __init__ (self, title="crys3d.wx_model_viewer", default_size=(800,600)) :
    self.title = title
    self.default_size = default_size
    wx.App.__init__(self, 0)

  def OnInit (self) :
    self.frame = wx.Frame(None, -1, self.title, pos=wx.DefaultPosition,
      size=self.default_size)
    self.frame.CreateStatusBar()
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = model_viewer_mixin(self.frame, size=(800,600))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    return True

def run (args) :
  if len(args) == 0 :
    print "Please specify a PDB file (and optional CIFs) on the command line."
    return
  a = App()
  processed_pdb_file = pdb_interpretation.run(args=args)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  pdb_hierarchy.atoms().reset_i_seq()
  grm = processed_pdb_file.geometry_restraints_manager()
  if grm is None or grm.shell_sym_tables is None :
    raise Sorry("Atomic bonds could not be calculated for this model. "+
      "This is probably due to a missing CRYST1 record in the PDB file.")
  atomic_bonds = grm.shell_sym_tables[0].full_simple_connectivity()
  #a.frame.Show()
  a.view_objects.add_model("1", pdb_hierarchy, atomic_bonds)
  a.frame.Show()
  a.view_objects.force_update(recenter=True)
  #a.view_objects.OnUpdate()
  a.MainLoop()

if __name__ == "__main__" :
  import sys
  run(sys.argv[1:])
