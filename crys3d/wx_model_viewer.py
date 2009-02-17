
from string import strip
from libtbx.utils import Sorry
import gltbx.util
from gltbx import wx_viewer, viewer_utils
from gltbx.gl import *
from gltbx.glu import *
from scitbx.array_family import flex, shared
from scitbx.math import minimum_covering_sphere
from mmtbx.monomer_library import pdb_interpretation
import wx
import sys

########################################################################
# BASE CLASS FOR DISPLAYING STRUCTURES
#
class model_viewer_mixin (wx_viewer.wxGLWindow) :
  initialize_model_viewer_super = True
  def __init__ (self, *args, **kwds) :
    if self.initialize_model_viewer_super :
      wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    self.pdb_hierarchy           = None
    self.atomic_bonds            = None # from geometry restraints manager
    self.selection_cache         = None # this is used by extract_trace
    self.atoms                   = None
    self.atom_count              = 0
    self.points                  = None # basic 3d data - atoms.extract_xyz()
    self.b_cache                 = None # atoms.extract_b()
    self.atom_index              = []   # stores atoms_with_labels data
    self.current_atom_i_seq      = None
    self.closest_point_i_seq     = None # usually set by mouse clicks
    self.minimum_covering_sphere = None
    # display lists - resetting these to None will force a redraw
    self.points_display_list     = None
    self.lines_display_list      = None
    self.spheres_display_list    = None
    # various settings; override these in subclasses and/or provide
    # a mechanism for changing their values
    self.line_width              = 1
    self.buffer_factor           = 2.0
    self.nonbonded_line_width    = 1
    self._fog_start              = 50
    self._fog_end                = 200
    self.draw_mode               = "all_atoms"
    self.color_mode              = "rainbow"
    self.recolor                 = self.color_rainbow
    self.base_atom_color         = (0.6, 0.6, 0.6) # grey
    self.carbon_atom_color       = (1.0, 1.0, 0.0) # yellow
    self.orthographic            = False
    # toggles for viewable objects
    self.flag_show_fog                     = True
    self.flag_show_lines                   = True
    self.flag_show_points                  = True
    self.flag_show_spheres                 = False
    self.flag_show_trace                   = False
    self.flag_show_hydrogens               = True
    self.flag_show_minimum_covering_sphere = False
    self.flag_show_rotation_center         = False

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
    gltbx.util.handle_error()

  def initialize_modelview (self) :
    if self.minimum_covering_sphere is not None :
      wx_viewer.wxGLWindow.initialize_modelview(self)

  def DrawGL(self):
    if (self.flag_show_points):
      self.draw_points()
    if (self.flag_show_lines):
      self.draw_lines()
    if (self.flag_show_spheres):
      self.draw_spheres()

  def OnRedrawGL (self, event=None) :
    if self.minimum_covering_sphere is None :
      gltbx.util.handle_error()
      glClear(GL_COLOR_BUFFER_BIT)
      glClear(GL_DEPTH_BUFFER_BIT)
      glFlush()
      self.SwapBuffers()
      gltbx.util.handle_error()
    else :
      wx_viewer.wxGLWindow.OnRedrawGL(self, event)

  def unzoom (self, event=None) :
    if self.atoms is not None :
      self.minimum_covering_sphere = minimum_covering_sphere(
                                       points=self.atoms.extract_xyz(),
                                       epsilon=0.1)
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  def recenter_on_atom (self, i_seq) :
    if self.points is not None and i_seq < self.points.size() :
      self.rotation_center = self.points[i_seq]

  def process_key_stroke (self, key) :
    if key == ord('u') :
      self.unzoom()

  def process_pick_points(self):
    if self.pick_points is not None :
      self.closest_point_i_seq = viewer_utils.closest_visible_point(
        points = self.points,
        atoms_visible = self.atoms_visible,
        point0 = self.pick_points[0],
        point1 = self.pick_points[1]
      )

  def update_view (self, redraw_points=True, redraw_lines=False) :
    self.get_drawn_atom_count()
    self.recolor()
    if redraw_lines :
      self.lines_display_list = None
    if redraw_points :
      self.points_display_list = None
    self.OnRedraw()

  def OnUpdate (self, event=None, recenter=False) :
    if self._structure_was_updated :
      self.extract_atom_data()
      self.extract_trace()
    else :
      self.update_coords()
    self._structure_was_updated = False
    self.update_view(True, True)
    if (event is not None and event.recenter == True) or recenter == True :
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  # Hopefully none of the remaining functions need to be overridden...

  #---------------------------------------------------------------------
  # coordinates, bonds, etc.
  def update_structure (self, pdb_hierarchy, atomic_bonds, redraw=False) :
    self.pdb_hierarchy = pdb_hierarchy
    self.atomic_bonds = atomic_bonds
    self.extract_atom_data()
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.set_draw_mode(self.draw_mode, redraw=False)
    self._structure_was_updated = True
    if redraw :
      self.OnUpdate(recenter=True)

  def extract_atom_data (self) :
    self.update_atoms(self.pdb_hierarchy.atoms(), redraw=False)
    self.atom_count = self.atoms.size()
    self.update_coords()
    atom_index = []
    for atom in self.pdb_hierarchy.atoms_with_labels() :
      atom_index.append(atom)
    self.atom_index = atom_index
    assert len(self.atom_index) == self.atom_count
    atom_radii = flex.double(self.points.size(), 1.5)
    for i_seq, atom in enumerate(self.atom_index) :
      if atom.element == ' H' :
        atom_radii[i_seq] = 0.75
    self.atom_radii = atom_radii

  def OnUpdate (self, event=None, recenter=False) :
    if self._structure_was_updated :
      self.extract_atom_data()
      self.extract_trace()
    else :
      self.update_coords()
    self._structure_was_updated = False
    self.update_view(True, True)
    if (recenter or
        (event is not None and getattr(event, "recenter", None) == True)) :
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()

  def update_coords (self) :
    self.points = self.atoms.extract_xyz()
    self.b_cache = self.atoms.extract_b()
    self.lines_display_list = None
    self.points_display_list = None
    self.selection_display_list = None

  def update_atoms (self, atoms, redraw=False) :
    self.atoms = atoms
    self.minimum_covering_sphere = minimum_covering_sphere(atoms.extract_xyz(),
      epsilon=1.e-1)
    if redraw :
      self.OnUpdate()

  def extract_trace (self) :
    last_atom     = None
    isel = self.selection_cache.iselection
    selection_i_seqs = list(
      isel("(name ' CA ' or name ' P  ') and (altloc 'A' or altloc ' ')"))
    last_atom = None
    bonds = shared.stl_set_unsigned(self.atoms.size())
    atoms = self.pdb_hierarchy.atoms_with_labels()
    for atom in atoms :
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

  def set_draw_mode (self, draw_mode, redraw=False) :
    if not self.minimum_covering_sphere :
      return
    self.draw_mode = draw_mode
    if draw_mode == "trace" :
      self.flag_show_points = False
      self.bonds = self.trace_bonds
    elif self.draw_mode == "bonded_only" :
      self.flag_show_points = False
      self.bonds = self.atomic_bonds
    else :
      self.flag_show_points = True
      self.bonds = self.atomic_bonds
    self.get_drawn_atom_count()
    self.set_color_mode(self.color_mode, redraw=redraw)

  def get_drawn_atom_count (self) :
    c = 0
    atoms = self.pdb_hierarchy.atoms_with_labels()
    if self.flag_show_hydrogens :
      atoms_drawable = flex.bool(self.atom_count, True)
    else :
      atoms_drawable = flex.bool([ (atom.element != ' H') for atom in atoms ])
    self.visibility = viewer_utils.atom_visibility(
      bonds             = self.bonds,
      atoms_drawable    = atoms_drawable,
      flag_show_points  = self.flag_show_points
    )
    self.atoms_visible = self.visibility.atoms_visible
    self.bonds_visible = self.visibility.bonds_visible
    self.points_visible = self.visibility.points_visible
    self.visible_atom_count = self.visibility.visible_atoms_count

  #---------------------------------------------------------------------
  # coloring
  def set_bg_color (self) :
    (r,g,b) = self.bg_color
    glClearColor(r, g, b, 0.0)

  def set_color_mode (self, color_mode, redraw=False) :
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
    if redraw :
      self.update_view(True, True)

  def color_mono (self) :
    self.atom_colors = flex.vec3_double(
      [ self.base_atom_color for i in xrange(0, self.points.size()) ]
    )

  def color_rainbow (self) :
    self.atom_colors = viewer_utils.color_rainbow(
      atoms_visible = self.atoms_visible,
      visible_atom_count = self.visible_atom_count
    )

  def color_b (self) :
    self.atom_colors = viewer_utils.color_by_property(
      atom_properties       = self.b_cache,
      atoms_visible         = self.atoms_visible,
      color_invisible_atoms = not self.scale_b_to_visible,
      use_rb_color_gradient = False
    )

  def color_by_chain (self) :
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

  def color_by_element (self) :
    # these are approximations based on my (probably faulty) memory.
    # feel free to change to something more reasonable.
    element_shades = {' C' : self.carbon_atom_color, # usually yellow or grey
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

  #---------------------------------------------------------------------
  # drawing functions
  def draw_lines (self):
    if self.lines_display_list is None :
      self.lines_display_list = gltbx.gl_managed.display_list()
      self.lines_display_list.compile()
      glLineWidth(self.line_width)
      viewer_utils.draw_bonds(
        points = self.points,
        bonds  = self.bonds,
        atom_colors = self.atom_colors,
        bonds_visible = self.bonds_visible)
      self.lines_display_list.end()
    self.lines_display_list.call()

  def draw_points (self) :
    if self.points_display_list is None :
      self.points_display_list = gltbx.gl_managed.display_list()
      self.points_display_list.compile()
      glLineWidth(self.nonbonded_line_width)
      viewer_utils.draw_points(
        points = self.points,
        atom_colors = self.atom_colors,
        points_visible = self.points_visible
      )
      self.points_display_list.end()
    self.points_display_list.call()

  def draw_spheres (self) :
    pass

  def _draw_spheres(self, spheres_visible, atom_colors, atom_radii,
      solid=False):
    glMatrixMode(GL_MODELVIEW)
    if (solid):
      glEnable(GL_LIGHTING)
      glEnable(GL_LIGHT0)
      glLightfv(GL_LIGHT0, GL_AMBIENT, [1, 1, 1, 1])
      glLightfv(GL_LIGHT0, GL_POSITION, [0, 0, 1, 0])
      glEnable(GL_BLEND)
      #glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
      #glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
      glMaterialfv(GL_FRONT, GL_DIFFUSE, [1,1,1,1.0])
      glMaterialfv(GL_FRONT, GL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
      sphere = gltbx.util.SolidSphere
      grid = 50
    else:
      sphere = gltbx.util.WireSphere
      grid = 20
    for i_seq, point in enumerate(self.points) :
      if spheres_visible[i_seq] :
        glColor3f(*atom_colors[i_seq])
        glPushMatrix()
        glTranslated(*point)
        sphere(radius=atom_radii[i_seq], slices=grid, stacks=grid)
        glPopMatrix()
    if (solid):
      glDisable(GL_LIGHTING)
      glDisable(GL_LIGHT0)
      glDisable(GL_BLEND)

########################################################################
# ATOM SELECTION VIEWER
# this will handle any valid atom selection recognized by cctbx.
#
class selection_viewer_mixin (model_viewer_mixin) :
  initialize_model_viewer_super = True
  def __init__ (self, *args, **kwds) :
    model_viewer_mixin.__init__(self, *args, **kwds)
    # objects used internally
    self.selection_covering_sphere = None
    self.selection_i_seqs = None
    self.selection_cache = None
    self.selected_points = []
    self.selection_display_list = None
    self.selection_colors = None
    self.atoms_selected = None
    # various settings
    self.selection_string = "None"
    self.selection_color = (1.0, 1.0, 1.0)
    self.selection_draw_mode = "bonds_and_atoms"
    # flags
    self.flag_show_selection = True
    self.flag_recolor_selection = True

  def DrawGL (self) :
    model_viewer_mixin.DrawGL(self)
    if self.flag_show_selection :
      self.draw_selection()

  def process_key_stroke (self, key) :
    model_viewer_mixin.process_key_stroke(self, key)
    if key == ord('z') :
      self.zoom_selection()

  def update_view (self, redraw_points=True, redraw_lines=False) :
    self.selection_display_list = None
    model_viewer_mixin.update_view(self, redraw_points, redraw_lines)

  def zoom_selection (self, event=None) :
    self.set_selection_sphere()
    if self.selection_covering_sphere is not None :
      self.minimum_covering_sphere = self.selection_covering_sphere
      self.move_rotation_center_to_mcs_center()
      self.fit_into_viewport()
    else :
      self.unzoom()

  def update_selection (self, selection_string=None) :
    self.selection_string = str(selection_string)
    self.apply_selection(self.selection_string)

  def apply_selection (self, selection_string) :
    if self.selection_cache is None :
      return
    sel = self.selection_cache.selection
    try :
      self.atoms_selected = sel(selection_string)
    except :
      self.atoms_selected = sel("none")
      raise Sorry("The string '%s' is not a valid selection."%selection_string)
    self.selection_i_seqs = self.atoms_selected.iselection()
    self.update_view()

  def get_selected_atom_count (self) :
    return self.selection_i_seqs.size()

  def set_selection_sphere (self) :
    if self.selection_i_seqs is None or self.selection_i_seqs.size() == 0 :
      self.selection_covering_sphere = None
      return
    if self.selection_i_seqs.size() == 0 :
      self.selection_covering_sphere = minimum_covering_sphere(
                                         points=self.atoms.extract_xyz(),
                                         epsilon=0.1)
    else :
      selected_points = flex.vec3_double()
      points = self.points
      for i_seq in self.selection_i_seqs :
        selected_points.append(points[i_seq])
      self.selection_covering_sphere = minimum_covering_sphere(
                                         points=selected_points,
                                         epsilon=0.1)

  def set_color_mode (self, color_mode, redraw=False) :
    self.set_sel_color(color_tuple=self.selection_color)
    model_viewer_mixin.set_color_mode(self, color_mode, redraw)

  def set_sel_color (self, color_tuple=(1.0,1.0,1.0)) :
    self.selection_color = color_tuple
    sele_color_list = flex.vec3_double()
    # does this need to be done in c++?
    for i_seq in xrange(0, self.atom_count) :
      sele_color_list.append(color_tuple)
    self.selection_colors = sele_color_list

  def set_selected_atom (self, closest_point_i_seq) :
    self.current_atom_i_seq = closest_point_i_seq

  #---------------------------------------------------------------------
  # DRAWING
  def draw_selection (self) :
    selection_i_seqs = self.selection_i_seqs
    if selection_i_seqs is None or selection_i_seqs.size() == 0 :
      return
    selection_colors = self.atom_colors
    if self.flag_recolor_selection and self.selection_colors is not None :
      selection_colors = self.selection_colors
    if self.selection_display_list is None :
      self.selection_display_list = gltbx.gl_managed.display_list()
      self.selection_display_list.compile()
      draw_mode = self.selection_draw_mode
      if draw_mode is None :
        draw_mode = "bonds_and_atoms"
      if draw_mode == "bonds_and_atoms" :
        self.visibility.get_selection_visibility(
          bonds          = self.bonds,
          atoms_selected = self.atoms_selected
        )
        glLineWidth(self.line_width + 2)
        viewer_utils.draw_points(
          points         = self.points,
          atom_colors    = selection_colors,
          points_visible = self.visibility.selected_points_visible,
          cross_radius   = 0.4
        )
        glLineWidth(self.line_width + 3)
        viewer_utils.draw_bonds(
          points        = self.points,
          bonds         = self.bonds,
          atom_colors   = selection_colors,
          bonds_visible = self.visibility.selected_bonds_visible
        )
      elif draw_mode == "points" :
        glLineWidth(self.line_width + 2)
        viewer_utils.draw_points(
          points         = self.points,
          atom_colors    = selection_colors,
          points_visible = self.atoms_selected,
          cross_radius   = 0.4
        )
      elif draw_mode == "spheres" :
        self._draw_spheres(
          spheres_visible = self.atoms_selected,
          atom_colors     = selection_colors,
          atom_radii      = self.atom_radii,
          solid           = True
        )
      self.selection_display_list.end()
    self.selection_display_list.call()

########################################################################
# MIXIN FOR DRAWING ATOM LABELS
#
# This is meant to be combined with other viewer mixins, and is much
# simpler as a result.
class atom_label_mixin (wx_viewer.wxGLWindow) :
  initialize_atom_label_super = False
  def __init__ (self, *args, **kwds) :
    if self.initialize_atom_label_super :
      wx_viewer.wxGLWindow.__init__(self, *args, **kwds)
    self.label_xyz = []
    self.label_text = []
    self.label_color = (1.0, 1.0, 1.0)
    self.labels_display_list = None
    self.flag_show_labels = True

  # No InitGL method here.

  def DrawGL (self) :
    if self.flag_show_labels :
      self.draw_labels()

  def clear_labels (self, event=None) :
    self.label_xyz = []
    self.label_text = []
    self.labels_display_list = None

  def draw_labels (self) :
    if (self.labels_display_list is None) :
      font = gltbx.fonts.ucs_bitmap_8x13
      font.setup_call_lists()
      self.labels_display_list = gltbx.gl_managed.display_list()
      self.labels_display_list.compile()
      glColor3f(*self.label_color)
      for label,point in zip(self.label_text, self.label_xyz):
        glRasterPos3f(*point)
        font.render_string(label)
      self.labels_display_list.end()
    self.labels_display_list.call()

  def show_atom_label (self, i_seq) :
    if self.points is None or i_seq >= self.points.size() :
      return
    self.label_xyz.append(self.points[i_seq])
    a = self.atom_index[i_seq]
    atom_str = "%s %s%s %s" % (strip(a.name), a.chain_id, strip(a.resseq),
      strip(a.resname))
    self.label_text.append(atom_str)
    self.OnRedrawGL()

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
  a.view_objects.update_structure(pdb_hierarchy, atomic_bonds, redraw=True)
  a.frame.Show()
  a.MainLoop()

if __name__ == "__main__" :
  import sys
  run(sys.argv[1:])
