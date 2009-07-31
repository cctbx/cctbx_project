
# TODO: move selection logic to separate module

from crys3d.wx_model_viewer_multi_scene import model_data, model_scene, \
    model_viewer_mixin
from crys3d.reverse_selection import mouse_selection_manager
import gltbx.gl_managed
from gltbx.gl import *
from gltbx.glu import *
from gltbx import viewer_utils
from scitbx.array_family import flex
import iotbx.phil
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import wx

viewer_phil = iotbx.phil.parse("""
include scope crys3d.wx_model_viewer_multi_scene.opengl_phil
selections {
  selection_color = 1.0 1.0 1.0
    .type = floats
    .style = color
  use_global_selection_color = True
    .type = bool
  line_padding = 2
    .type = int
}
""", process_includes=True)

class model_data_with_selection (model_data, mouse_selection_manager) :
  def __init__ (self, *args, **kwds) :
    self.flag_show_selection = True
    self.flag_allow_selection = True
    mouse_selection_manager.__init__(self)
    model_data.__init__(self, *args, **kwds)

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
    scene.update_selection(self.atom_selection)
    model_data.update_scene_data(self, scene)

  def update_structure (self, pdb_hierarchy, atomic_bonds,
      mmtbx_selection_function=None) :
    model_data.update_structure(self, pdb_hierarchy, atomic_bonds)
    mouse_selection_manager.update_selection_handlers(self, pdb_hierarchy,
      mmtbx_selection_function)

  def recalculate_visibility (self) :
    model_data.recalculate_visibility(self)
    self.visibility.get_selection_visibility(
      bonds          = self.current_bonds,
      atoms_selected = self.atom_selection)

  def selection_callback (self, *args, **kwds) :
    self.recalculate_visibility()
    mouse_selection_manager.selection_callback(self, *args, **kwds)

#-----------------------------------------------------------------------
class model_scene_with_selection (model_scene) :
  def __init__ (self, *args, **kwds) :
    model_scene.__init__(self, *args, **kwds)
    self.selection_draw_mode = "bonds_and_atoms"
    self.update_selection(flex.bool(self.points.size(), False))

  def clear_lists (self) :
    self.selection_display_list = None
    model_scene.clear_lists(self)

  def update_selection (self, atom_selection) :
    self.atom_selection = atom_selection
    self.selected_atom_count = atom_selection.iselection().size()
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
          points_visible = self.atom_selection,
          cross_radius   = 0.4)
      elif draw_mode == "spheres" :
        atom_selection = self.atom_selection
        atom_radii = self.atom_radii
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        for i_seq, point in enumerate(self.points) :
          if atom_selection[i_seq] :
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
mouse_modes = ["Rotate view", "Toggle chain", "Toggle residue", "Toggle atom",
  "Select range", "Deselect range", "Show selection menu"]
class selection_editor_mixin (model_viewer_mixin) :
  def __init__ (self, *args, **kwds) :
    self.left_button_mode = 0
    self.flag_select_all_conformers = True
    self.flag_show_selections = True
    self.flag_enable_mouse_selection = True
    self.current_atom_i_seq = None
    self.current_object_id = None
    self._in_range_selection = False
    model_viewer_mixin.__init__(self, *args,**kwds)
    self.settings = viewer_phil.extract()

  def DrawGL (self) :
    model_viewer_mixin.DrawGL(self)
    if self.flag_show_selections :
      self.draw_selections()

  def draw_selections (self) :
    line_width = (self.settings.opengl.line_width +
                  (2 * self.settings.selections.line_padding))
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
      mmtbx_selection_function=None) :
    assert isinstance(model_id, str)
    model = model_data_with_selection(model_id, pdb_hierarchy, atomic_bonds,
      base_color=self.settings.opengl.base_atom_color)
    model.set_mmtbx_selection_function(mmtbx_selection_function)
    model.set_selection_callback(print_cb)
    self._add_model(model_id, model)

  def process_key_stroke (self, key) :
    if key == ord('S') :
      self.flag_show_selections = not self.flag_show_selections
      self.update_scene = True
    elif key == ord('z') :
      self.zoom_selections()
    elif key >= 49 and key <= 55 : # 1-7
      self.set_left_button_mode(key - 49)
    elif key == 27 : # escape
      self.clear_selections()
      self.update_scene = True
    model_viewer_mixin.process_key_stroke(self,key)

  def set_left_button_mode (self, mode) :
    self.left_button_mode = mode
    print "left button mode is %s" % mouse_modes[mode]
    self._in_range_selection = False

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
      model.toggle_residue_selection(self.current_atom_i_seq,
        ignore_altloc=self.flag_select_all_conformers)
      self.update_scene = True

  def toggle_atom_selection (self) :
    model = self.get_model(self.current_object_id)
    if model is not None :
      model.toggle_atom_selection(self.current_atom_i_seq)
      self.update_scene = True

  def process_range_selection (self, deselect=False) :
    model = self.get_model(self.current_object_id)
    if model is not None :
      if self._in_range_selection :
        try :
          model.end_range_selection(self.current_atom_i_seq, deselect,
            ignore_altloc=self.flag_select_all_conformers)
        finally :
          self._in_range_selection = False
      else :
        model.start_range_selection(self.current_atom_i_seq)
        self._in_range_selection = True
      self.update_scene = True

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
        if self.left_button_mode == 1 :   # (de)select chain
          self.toggle_chain_selection()
        elif self.left_button_mode == 2 : # (de)select residue
          self.toggle_residue_selection()
        elif self.left_button_mode == 3 : # (de)select atom
          self.toggle_atom_selection()
        elif self.left_button_mode == 4 : # select range
          self.process_range_selection()
        elif self.left_button_mode == 5 : # deselect range
          self.process_range_selection(deselect=True)
        else:                             # Selection menu
          self.context_selection_menu()

  def OnToggleChain (self, event) :
    pass

  def OnToggleResidue (self, event) :
    pass

def print_cb (selection_string, atom_selection) :
  print "%s (%s)" % (selection_string, atom_selection.iselection().size())

#---end
