
from crys3d.wx_model_viewer_multi_scene import model_data, model_scene, \
    model_viewer_mixin
from gltbx.gl import *
from gltbx.glu import *
from gltbx import viewer_utils
from scitbx.array_family import flex
import iotbx.phil
from libtbx.utils import Sorry

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

class model_data_with_selection (model_data) :
  def __init__ (self, *args, **kwds) :
    self.mmtbx_handler = None
    self.flag_show_selection = False
    self.flag_allow_selection = True
    self.selection_string = "none"
    self.saved_selection = "none" # for recovery later
    self.selection_i_seqs = None
    model_data.__init__(self, *args, **kwds)
    self.selection_colors = flex.vec3_double(self.atoms.size(), (1,1,1))
    self.atoms_selected = flex.bool(self.atoms.size(), False)
    self.selection_i_seqs = self.atoms_selected.iselection()

  def set_mmtbx_selection_handler (self, mmtbx_handler) :
    self.mmtbx_handler = mmtbx_handler

  def get_scene_data (self) :
    return model_scene_with_selection(bonds=self.current_bonds,
      points=self.atoms.extract_xyz(),
      b_iso=self.atoms.extract_b(),
      b_aniso=self.atoms.extract_uij(),
      atom_colors=self.atom_colors,
      atom_labels=self.atom_labels,
      atom_radii=self.atom_radii,
      visibility=self.visibility)

  def update_scene_data (self, scene) :
    scene.update_selection(self.atoms_selected)
    model_data.update_scene_data(self, scene)

  def update_structure (self, pdb_hierarchy, atomic_bonds,
      mmtbx_handler=None) :
    model_data.update_structure(self, pdb_hierarchy, atomic_bonds)
    self.selection_cache = pdb_hierarchy.atom_selection_cache()
    self.mmtbx_handler = mmtbx_handler
    self.apply_selection(self.saved_selection)

  def selection_size (self) :
    return self.selection_i_seqs.size()

  def apply_selection (self, selection_string) :
    try :
      if self.mmtbx_handler is not None :
        atoms_selected = self.mmtbx_handler.selection(
          string=selection_string,
          cache=self.selection_cache)
      else :
        atoms_selected = self.selection_cache.selection(selection_string)
    except KeyboardInterrupt :
      raise
    except Exception :
      self.selection_string = "none"
      atoms_selected = self.selection_cache.selection("none")
      raise Sorry("Invalid selection '%s'."%selection_string)
    else :
      self.selection_string = selection_string
    finally :
      self.atoms_selected = atoms_selected
      self.selection_i_seqs = atoms_selected.iselection()

  def revert_selection (self) :
    self.apply_selection(self.saved_selection)

  def set_initial_selection (self, selection_string) :
    self.saved_selection = selection_string
    self.apply_selection(selection_string)

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
    self.selection_i_seqs = atoms_selected.iselection()
    self.selection_display_list = None

  def get_selected_xyz (self) :
    points = self.points
    for i_seq in self.selection_i_seqs :
      yield points[i_seq]

  # XXX: this is still gross.
  def draw_selection (self, color, use_atom_colors=False) :
    selection_i_seqs = self.selection_i_seqs
    if selection_i_seqs is None or selection_i_seqs.size() == 0 :
      return
    if self.selection_display_list is None :
      if use_atom_colors :
        selection_colors = self.atom_colors
      else :
        selection_colors = flex.vec3_double(self.points.size(), color)
      self.selection_display_list = gltbx.gl_managed.display_list()
      self.selection_display_list.compile()
      draw_mode = self.selection_draw_mode
      if draw_mode is None :
        draw_mode = "bonds_and_atoms"
      if draw_mode == "bonds_and_atoms" :
        self.visibility.get_selection_visibility(
          bonds          = self.bonds,
          atoms_selected = self.atoms_selected)
        viewer_utils.draw_points(
          points         = self.points,
          atom_colors    = selection_colors,
          points_visible = self.visibility.selected_points_visible,
          cross_radius   = 0.4)
        viewer_utils.draw_bonds(
          points        = self.points,
          bonds         = self.bonds,
          atom_colors   = selection_colors,
          bonds_visible = self.visibility.selected_bonds_visible)
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
# mouse controls are implemented separately in selection_editor_mixin.
class selection_viewer_mixin (model_viewer_mixin) :
  def __init__ (self, *args, **kwds) :
    model_viewer_mixin.__init__(self, *args,**kwds)
    self.selection_sphere = None
    self.flag_show_selections = True
    self.settings = viewer_phil.extract()

  def DrawGL (self) :
    model_viewer_mixin.DrawGL(self)
    if self.flag_show_selections :
      self.draw_selections()

  def draw_selections (self) :
    line_width = (self.settings.opengl.line_width +
                  self.settings.selections.line_padding)
    glLineWidth(line_width)
    for object_id, scene in self.scene_objects.iteritems() :
      if self.show_object[object_id] :
        scene.draw_selection(color=self.settings.selections.selection_color,
          use_atom_colors=self.settings.selections.use_global_selection_color)

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

  def add_model (self, model_id, pdb_hierarchy, atomic_bonds,
      mmtbx_handler=None) :
    assert isinstance(model_id, str)
    model = model_data_with_selection(model_id, pdb_hierarchy, atomic_bonds,
      base_color=self.settings.opengl.base_atom_color)
    self.model_ids.append(model_id)
    self.model_objects.append(model)
    self.show_object[model_id] = True
    model.set_mmtbx_selection_handler(mmtbx_handler)
    self.update_scene = True

  def process_key_stroke (self, key) :
    if key == ord('s') :
      self.flag_show_selections = not self.flag_show_selections
      self.update_scene = True
    elif key == ord('z') :
      self.zoom_selections()
    model_viewer_mixin.process_key_stroke(self,key)

  def set_selected_atom (self, i_seq) :
    self.current_atom_i_seq = i_seq

#-----------------------------------------------------------------------
class selection_editor_mixin (selection_viewer_mixin) :
  pass

#---end
