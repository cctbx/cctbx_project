from __future__ import absolute_import, division, print_function
from six.moves import range

# TODO: clean up handling of changes in atom count

# these are approximations based on my (probably faulty) memory.
# feel free to change to something more reasonable.
# note: carbon is assigned the base color.
element_shades = {'H'  : (0.95, 0.95, 0.95), # very light grey
                  #'C'  : (0.8, 0.8, 0.8),    # light grey
                  'N'  : (0.0, 0.0, 1.0),    # blue
                  'O'  : (1.0, 0.0, 0.0),    # red
                  'S'  : (1.0, 0.5, 0.0),    # orange
                  'P'  : (1.0, 1.0, 0.0),    # yellow
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

class model_data (object) :
  def __init__ (self, object_id, pdb_hierarchy, atomic_bonds,
      special_position_settings=None,
      base_color=(0.0,1.0,1.0)) :
    self.object_id = object_id
    self.base_color = base_color
    self.draw_mode = None
    self.current_bonds = None
    self.noncovalent_bonds = None
    self.ribbon = None
    self.color_mode = None #"rainbow"
    self.flag_object_visible = True
    self._color_cache = {}
    self.flag_show_hydrogens = False
    self.flag_show_lines = True
    self.flag_show_labels = True
    self.flag_show_points = True
    self.flag_show_spheres = False
    self.flag_show_ribbon = False
    self.flag_show_ellipsoids = False
    self.flag_show_noncovalent_bonds = False
    self.update_structure(pdb_hierarchy=pdb_hierarchy,
      atomic_bonds=atomic_bonds,
      special_position_settings=special_position_settings)
    from scitbx.array_family import flex
    self.use_u_aniso = flex.bool(self.atoms.size())
    #self.recalculate_visibility()

  def reset (self) :
    self.is_changed = False

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
      noncovalent_bonds=self.noncovalent_bonds,
      atomic_bonds=self.atomic_bonds,
      ribbon=self.ribbon)

  def set_noncovalent_bonds (self, bonded_atoms) :
    self.noncovalent_bonds = bonded_atoms
    self.is_changed = True

  def update_scene_data (self, scene) :
    scene.update_bonds(self.current_bonds)
    scene.update_colors(self.atom_colors)
    scene.update_visibility(self.visibility)
    scene.clear_lists()

  def update_xyz (self, xyz) :
    assert xyz.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.xyz = xyz[i_seq]
    self.is_changed = True

  def update_u_iso (self, u_iso) :
    assert u_iso.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.b = adptbx.u_as_b(u_iso[i_seq])
    self.is_changed = True

  def update_u_aniso (self, u_aniso, aniso_flag=None) :
    assert u_aniso.size() == self.atoms.size()
    for i_seq, atom in enumerate(self.atoms) :
      atom.uij = u_aniso[i_seq]
    self.is_changed = True

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

  def update_structure (self, pdb_hierarchy, atomic_bonds,
      special_position_settings=None) :
    from scitbx.array_family import flex
    self.pdb_hierarchy = pdb_hierarchy
    self.atoms = pdb_hierarchy.atoms()
    self.atom_count = self.atoms.size()
    assert (self.atom_count == len(atomic_bonds))
    if atomic_bonds is None :
      atomic_bonds = flex.stl_set_unsigned(self.atom_count)
    self.atomic_bonds = atomic_bonds
    self.selection_cache = pdb_hierarchy.atom_selection_cache(
      special_position_settings=special_position_settings)
    #self.index_atoms()
    atom_index = []
    atom_labels = flex.std_string()
    for atom in self.pdb_hierarchy.atoms_with_labels() :
      atom_index.append(atom)
      atom_labels.append(format_atom_label(atom))
    self.atom_index = atom_index
    self.atom_labels = atom_labels
    self.trace_bonds = extract_trace(pdb_hierarchy, self.selection_cache)
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

  def recalculate_visibility (self) :
    from scitbx.array_family import flex
    from gltbx import viewer_utils
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

  def initialize_cartoon (self, sec_str=None) :
    if (sec_str is None) :
      from mmtbx import secondary_structure
      manager = secondary_structure.manager(
        pdb_hierarchy=self.pdb_hierarchy,
        xray_structure=None)
      sec_str = manager.selections_as_ints()
    from crys3d import ribbon
    self.ribbon = ribbon.cartoon(pdb_hierarchy=self.pdb_hierarchy,
      sec_str=sec_str)
    self.ribbon.construct_geometry()

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
      elif draw_mode == "ribbon" :
        self.flag_show_ribbon = True
        self.flag_show_lines = False
        self.flag_show_points = False
        if (self.ribbon is None) :
          self.initialize_cartoon()
      else :
        self.flag_show_ribbon = False
        self.flag_show_lines = True
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
      n_visible = self.visibility.atoms_visible.count(True)
      if (n_visible == 0) :
        return
      elif (n_visible == 1) :
        color_mode = "mono"
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

  def color_mono (self) :
    cached = self._color_cache.get("mono")
    if cached is not None :
      self.atom_colors = cached
    else :
      from scitbx.array_family import flex
      self.atom_colors = flex.vec3_double(
        [ self.base_color for i in range(0, self.atoms.size()) ]
      )
      self._color_cache["mono"] = self.atom_colors

  def color_rainbow (self) :
    cached = self._color_cache.get("rainbow")
    if cached is not None :
      self.atom_colors = cached
    else :
      from scitbx import graphics_utils
      if (self.visibility.atoms_visible.count(True) == 1) :
        from scitbx.array_family import flex
        self.atom_colors = flex.vec3_double([(0,0,1)])
        return
      self.atom_colors = graphics_utils.color_rainbow(
        selection=self.visibility.atoms_visible)
      self._color_cache["rainbow"] = self.atom_colors

  def color_b (self) :
    cached = self._color_cache.get("b")
    if cached is not None :
      self.atom_colors = cached
    else :
      from scitbx import graphics_utils
      self.atom_colors = graphics_utils.color_by_property(
        properties=self.atoms.extract_b(),
        selection=self.visibility.atoms_visible,
        color_all=False,
        gradient_type="rainbow")
      self._color_cache["b"] = self.atom_colors

  def color_by_chain (self) :
    cached = self._color_cache.get("chain")
    if cached is not None :
      self.atom_colors = cached
    else :
      from scitbx.array_family import flex
      from scitbx import graphics_utils
      c = 0
      for chain in self.pdb_hierarchy.chains() :
        c += 1
      rainbow = graphics_utils.make_rainbow_gradient(c)
      j = 0
      chain_shades = {}
      for chain in self.pdb_hierarchy.chains() :
        chain_shades[chain.id] = rainbow[j]
      atom_colors = flex.vec3_double()
      for atom in self.pdb_hierarchy.atoms_with_labels() :
        atom_colors.append(chain_shades[atom.chain_id])
      self.atom_colors = atom_colors
      self._color_cache["chain"] = atom_colors

  def color_by_element (self) :
    cached = self._color_cache.get("element")
    if cached is not None :
      self.atom_colors = cached
    else :
      from scitbx.array_family import flex
      atom_colors = flex.vec3_double()
      for atom in self.pdb_hierarchy.atoms_with_labels() :
        element = atom.element.strip()
        color = element_shades.get(element, self.base_color)
        atom_colors.append(color)
      self.atom_colors = atom_colors
      self._color_cache["element"] = cached

#-----------------------------------------------------------------------
# Utility functions
def extract_trace (pdb_hierarchy, selection_cache=None) :
  from scitbx.array_family import shared
  if selection_cache is None :
    selection_cache = pdb_hierarchy.atom_selection_cache()
  last_atom     = None
  vertices = selection_cache.selection(
    "(name ' CA ' or name ' P  ') and (altloc 'A' or altloc ' ')")
  last_i_seq = None
  last_labels = None
  atoms = pdb_hierarchy.atoms()
  bonds = shared.stl_set_unsigned(atoms.size())
  for i_seq, atom in enumerate(atoms) :
    labels = atom.fetch_labels()
    if vertices[i_seq] :
      if last_i_seq is not None :
        if (labels.chain_id        == last_labels.chain_id and
            labels.model_id        == last_labels.model_id and
            labels.resseq_as_int() == (last_labels.resseq_as_int() + 1) and
            ((labels.altloc == last_labels.altloc) or
             (labels.altloc == "A" and last_labels.altloc == "") or
             (labels.altloc == ""  and last_labels.altloc == "A"))) :
          bonds[last_i_seq].append(i_seq)
          bonds[i_seq].append(last_i_seq)
      last_i_seq = i_seq
      last_labels = labels
  return bonds

def format_atom_label (atom_info) :
  return ("%s %s%s %s %s" % (atom_info.name, atom_info.altloc,
        atom_info.resname, atom_info.chain_id, atom_info.resid())).strip()
