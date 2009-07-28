
from crys3d.wx_model_viewer_multi_scene import model_data, model_scene, \
    model_viewer_mixin
import gltbx.gl_managed
from gltbx.gl import *
from gltbx.glu import *
from gltbx import viewer_utils
from scitbx.array_family import flex
import iotbx.phil
from libtbx.utils import Sorry
import wx

viewer_phil = iotbx.phil.parse("""
include scope crys3d.wx_model_viewer_multi_scene.opengl_phil
selections {
  selection_color = 1.0 1.0 1.0
    .type = floats
    .style = color
  use_global_selection_color = True
    .type = bool
  line_padding = 4
    .type = int
}
""", process_includes=True)

class chain_selection_info (object) :
  def __init__ (self, chain_id) :
    adopt_init_args(self, locals())

  def contains (self, selection_info) :
    if selection_info.chain_id == self.chain_id :
      return True
    return False

  def __str__ (self) :
    return "chain '%s'" % self.chain_id

class residue_selection_info (object) :
  def __init__ (self, chain_id, resid, altloc=None) :
    adopt_init_args(self, locals())

  def contains (self, selection_info) :
    if (selection_info.chain_id == self.chain_id and
        selection_info.resid == self.resid) :
      return True
    return False

  def __str__ (self) :
    if self.altloc is None :
      return "(chain '%s' and resid '%s')" % (self.chain_id, self.resid)
    else :
      return "(chain '%s' and resid '%s' and altloc '%s')" % (self.chain_id,
        self.resid, self.altloc)

class atom_selection_info (object) :
  def __init__ (self, chain_id, resid, atom_name, altloc) :
    adopt_init_args(self, locals())

  def __str__ (self) :
    return "(chain '%s' and resid '%s' and name '%s' and altloc '%s')"

class model_data_with_selection (model_data) :
  def __init__ (self, *args, **kwds) :
    self.mmtbx_handler = None
    self.flag_show_selection = True
    self.flag_allow_selection = True
    self.selection_string = "none"
    self.saved_selection = "none" # for recovery later
    self.selection_i_seqs = None
    self.selection_callback = None
    self.start_i_seq = None
    self.end_i_seq = None
    model_data.__init__(self, *args, **kwds)

  def set_mmtbx_selection_handler (self, mmtbx_handler) :
    self.mmtbx_handler = mmtbx_handler

  def set_selection_callback (self, callback) :
    self.selection_callback = callback

  def get_scene_data (self) :
    scene = model_scene_with_selection(bonds=self.current_bonds,
      points=self.atoms.extract_xyz(),
      b_iso=self.atoms.extract_b(),
      b_aniso=self.atoms.extract_uij(),
      atom_colors=self.atom_colors,
      atom_labels=self.atom_labels,
      atom_radii=self.atom_radii,
      visibility=self.visibility)
    self.update_scene_data(scene)
    return scene

  def update_scene_data (self, scene) :
    scene.update_selection(self.atoms_selected)
    model_data.update_scene_data(self, scene)

  def update_structure (self, pdb_hierarchy, atomic_bonds,
      mmtbx_handler=None) :
    model_data.update_structure(self, pdb_hierarchy, atomic_bonds)
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.mmtbx_handler = mmtbx_handler
    self.clear_selection()
    #self.apply_selection(self.saved_selection)

  def recalculate_visibility (self) :
    model_data.recalculate_visibility(self)
    self.visibility.get_selection_visibility(
      bonds          = self.current_bonds,
      atoms_selected = self.atoms_selected)

  #---------------------------------------------------------------------
  #
  def clear_selection (self) :
    self.selected_chains = {}
    self.deselected_chains = {}
    self.selected_atoms = []
    self.deselected_atoms = []
    self.selected_residues = []
    self.deselected_residues = []
    self._cached_selections = {}
    self._cached_deselections = {}
    self.apply_selection(self.saved_selection)

  def selection_size (self) :
    return self.selection_i_seqs.size()

  def apply_selection (self, selection_string) :
    atoms_selected = self.get_atom_selection(selection_string)
    if atoms_selected is None :
      self.selection_string = "none"
      atoms_selected = self.selection_cache.selection("none")
      raise Sorry("Invalid selection '%s'."%selection_string)
    else :
      self.selection_string = selection_string
    self.atoms_selected = atoms_selected
    self.selection_i_seqs = atoms_selected.iselection()
    self.recalculate_visibility()
    if self.selection_callback is not None :
      self.selection_callback(self.selection_string, self.atoms_selected)

  def get_atom_selection (self, selection_string) :
    try :
      if self.mmtbx_handler is not None :
        atoms_selected = self.mmtbx_handler.selection(
          string=selection_string,
          cache=self.selection_cache)
      else :
        atoms_selected = self.selection_cache.selection(selection_string)
    except KeyboardInterrupt :
      raise
    except Exception, e :
      atoms_selected =None
    return atoms_selected

  def revert_selection (self) :
    self.apply_selection(self.saved_selection)

  def set_initial_selection (self, selection_string) :
    self.saved_selection = selection_string
    self.apply_selection(selection_string)

  def start_selection (self, i_seq) :
    self.start_i_seq = i_seq

  def end_range_selection (self, i_seq, deselect=False) :
    pass

  def toggle_chain_selection (self, i_seq) :
    atom = self.atom_index[i_seq]
    chain_id = atom.chain_id
    if chain_id in self.selected_chains :
      self.selected_chains.pop(chain_id)
    else :
      self.selected_chains[chain_id] = True
    self.construct_selection()

  def toggle_residue_selection (self, i_seq) :
    atom = self.atom_index[i_seq]
    resi_sel_str = "chain '%s' and resid '%s'" % (atom.chain_id, atom.resid())
    resi_sel = self.selection_cache.selection(resi_sel_str)
    if self.atoms_selected.is_super_set(resi_sel) :
      self.deselected_residues.append(resi_sel_str)
    else :
      self.selected_residues.append(resi_sel_str)
    self.construct_selection()

  def toggle_atom_selection (self, i_seq, deselect=False) :
    pass

  def construct_selection (self) :
    final_selection = ""
    # Part 1: stuff we want
    positive_selection = ""
    clauses = []
    for chain in sorted(self.selected_chains.keys()) :
      chain_selection = "chain '%s'" % chain
      clauses.append(chain_selection)
    for residue in self.selected_residues :
      clauses.append(residue)
    for atom in self.selected_atoms :
      clauses.append(atom)
    if len(clauses) == 1 : # not necessary, but looks nicer
      positive_selection = clauses[0]
    elif len(clauses) > 1 :
      positive_selection = "(" + ") or (".join(clauses) + ")"
    # Part 2: stuff we don't want
    negative_selection = ""
    clauses = []
    for residue in self.deselected_residues :
      clauses.append(residue)
    for atom in self.deselected_atoms :
      clauses.append(atom)
    if len(clauses) == 1 :
      negative_selection = clauses[0]
    elif len(clauses) > 1 :
      negative_selection = "(" + ") or (".join(clauses) + ")"
    # assemble final selection
    if positive_selection != "" :
      if negative_selection != "" :
        final_selection = "(%s) and not (%s)" % (positive_selection,
                                                 negative_selection)
      else :
        final_selection = positive_selection
    #elif negative_selection != "" :
    #  final_selection = "not (%s)" % negative_selection
    else :
      final_selection = "none"
    self.apply_selection(final_selection)

#-----------------------------------------------------------------------
class model_scene_with_selection (model_scene) :
  def __init__ (self, *args, **kwds) :
    model_scene.__init__(self, *args, **kwds)
    self.selection_draw_mode = "bonds_and_atoms"
    self.update_selection(flex.bool(self.points.size(), False))

  def clear_lists (self) :
    self.selection_display_list = None
    model_scene.clear_lists(self)

  def update_selection (self, atoms_selected) :
    self.atoms_selected = atoms_selected
    self.selected_atom_count = atoms_selected.iselection().size()
    self.selection_display_list = None

  def update_visibility (self, visibility) :
    model_scene.update_visibility(self, visibility)
    self.selected_points_visible = visibility.selected_points_visible
    self.selected_bonds_visible = visibility.selected_bonds_visible

  def get_selected_xyz (self) :
    points = self.points
    for i_seq in self.selection_i_seqs :
      yield points[i_seq]

  # XXX: this is still gross.
  def draw_selection (self, color, use_global_color=False) :
    if self.selected_atom_count == 0 :
      return
    if self.selection_display_list is None :
      if use_global_color :
        selection_colors = flex.vec3_double(self.points.size(), color)
      else :
        selection_colors = self.atom_colors
      self.selection_display_list = gltbx.gl_managed.display_list()
      self.selection_display_list.compile()
      draw_mode = self.selection_draw_mode
      if draw_mode is None :
        draw_mode = "bonds_and_atoms"
      if draw_mode == "bonds_and_atoms" :
        viewer_utils.draw_points(
          points         = self.points,
          atom_colors    = selection_colors,
          points_visible = self.selected_points_visible,
          cross_radius   = 0.4)
        viewer_utils.draw_bonds(
          points        = self.points,
          bonds         = self.bonds,
          atom_colors   = selection_colors,
          bonds_visible = self.selected_bonds_visible)
      elif draw_mode == "points" :
        viewer_utils.draw_points(
          points         = self.points,
          atom_colors    = selection_colors,
          points_visible = self.atoms_selected,
          cross_radius   = 0.4)
      elif draw_mode == "spheres" :
        atoms_selected = self.atoms_selected
        atom_radii = self.atom_radii
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        for i_seq, point in enumerate(self.points) :
          if atoms_selected[i_seq] :
            glColor3f(*selection_colors[i_seq])
            glPushMatrix()
            glTranslated(*point)
            gltbx.util.SolidSphere(radius=atom_radii[i_seq],
                                   slices=50, stacks=50)
            glPopMatrix()
      self.selection_display_list.end()
    self.selection_display_list.call()

########################################################################
# VIEWER CLASS
mouse_modes = ["Rotate view", "Show selection menu", "Toggle chain selection",
               "Select residue", "Select atom", "Select range",
                "Deselect range"]
class selection_editor_mixin (model_viewer_mixin) :
  def __init__ (self, *args, **kwds) :
    self.left_button_mode = 0
    self.flag_select_all_conformers = True
    self.flag_show_selections = True
    self.flag_enable_mouse_selection = True
    self.current_atom_i_seq = None
    self.current_object_id = None
    model_viewer_mixin.__init__(self, *args,**kwds)
    self.settings = viewer_phil.extract()

  def DrawGL (self) :
    model_viewer_mixin.DrawGL(self)
    if self.flag_show_selections :
      self.draw_selections()

  def draw_selections (self) :
    line_width = (self.settings.opengl.line_width +
                  self.settings.selections.line_padding)
    glLineWidth(line_width)
    glEnable(GL_LINE_SMOOTH)
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_selection(color=self.settings.selections.selection_color,
          use_global_color=self.settings.selections.use_global_selection_color)

  #---------------------------------------------------------------------
  def zoom_selections (self) :
    points = flex.vec3_double()
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        for point in scene.get_selected_xyz() :
          points.append(point)
    if points.size() != 0 :
      self.update_mcs(points)

  def set_selection (self, object_id, selection_string) :
    for model_id, model in self.iter_models() :
      if model_id == object_id :
        model.apply_selection(selection_string)
        self.update_scene = True

  def clear_selections (self) :
    for model_id, model in self.iter_models() :
      model.clear_selection()

  def add_model (self, model_id, pdb_hierarchy, atomic_bonds,
      mmtbx_handler=None) :
    assert isinstance(model_id, str)
    model = model_data_with_selection(model_id, pdb_hierarchy, atomic_bonds,
      base_color=self.settings.opengl.base_atom_color)
    model.set_mmtbx_selection_handler(mmtbx_handler)
    model.set_selection_callback(print_cb)
    self._add_model(model_id, model)

  def process_key_stroke (self, key) :
    if key == ord('S') :
      self.flag_show_selections = not self.flag_show_selections
      self.update_scene = True
    elif key == ord('z') :
      self.zoom_selections()
    elif key >= 49 and key <= 54 :
      self.left_button_mode = key - 49
    elif key == 27 : # escape
      self.clear_selections()
      self.update_scene = True
    model_viewer_mixin.process_key_stroke(self,key)

  #---------------------------------------------------------------------
  def pick_selection_object (self, object_id) :
    for model_id in self.pick_object :
      if model_id == object_id :
        self.pick_object[model_id] = True
        self.show_object[model_id] = True
      else :
        self.pick_object[model_id] = False

  def save_selected_atom (self) : #, object_id, i_seq) :
    self.current_object_id = self.closest_point_object_id
    self.current_atom_i_seq = self.closest_point_i_seq

  def toggle_chain_selection (self) :
    model = self.get_model(self.current_object_id)
    if model is not None :
      model.toggle_chain_selection(self.current_atom_i_seq)
      self.update_scene = True

  def toggle_residue_selection (self) :
    model = self.get_model(self.current_object_id)
    if model is not None :
      model.toggle_residue_selection(self.current_atom_i_seq)
      self.update_scene = True

  def process_range_selection (self) :
    pass

  def context_selection_menu (self) :
    menu = wx.Menu()
    toggle_chain = menu.Append(-1, "Toggle chain selection")
    self.Bind(wx.EVT_MENU, self.OnToggleChain, toggle_chain)
    self.PopupMenu(menu)
    menu.Destroy()

  def OnLeftClick (self, event) :
    model_viewer_mixin.OnLeftClick(self, event)
    if self.left_button_mode != 0 and not self.was_dragged :
      self.get_pick_points((self.xmouse, self.ymouse))
      self.process_pick_points()
      if (self.closest_point_i_seq is not None and
          self.flag_enable_mouse_selection) :
        self.save_selected_atom()
        if self.left_button_mode == 1 :   # Selection menu
          self.context_selection_menu()
        elif self.left_button_mode == 2 : # (de)select chain
          self.toggle_chain_selection()
        elif self.left_button_mode == 3 : # (de)select residue
          self.toggle_residue_selection()
        elif self.left_button_mode == 4 : # (de)select atom
          self.toggle_atom_selection()
        elif self.left_button_mode == 5 : # select range
          self.process_range_selection()
        else :                            # deselect range
          self.process_range_deselection()

  def OnToggleChain (self, event) :
    pass

  def OnToggleResidue (self, event) :
    pass

def print_cb (selection_string, atoms_selected) :
  print "%s (%s)" % (selection_string, atoms_selected.iselection().size())

def is_subset (selection1, selection2) :
  union = selection1 & selection2
  return union.count(True) == selection2.count(True)

#---end
