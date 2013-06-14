
from __future__ import division
from libtbx import slots_getstate_setstate

class water_properties (slots_getstate_setstate) :
  __slots__ = [
    "i_seq",
    "id_str",
    "site_cart",
    "b_iso",
    "occ",
    "nearest_contact",
    "map_cc",
    "two_fofc_map",
    "fofc_map",
    "anom_map",
    "n_hbonds",
  ]
  def __init__ (self, **kwds) :
    for name in self.__slots__ :
      setattr(self, name, kwds.get(name, None))

  def is_bad_water (self) :
    return ((self.map_cc < 0.8) or (self.fofc_map < -3.0) or (self.occ < 0.5)
            or (self.two_fofc_map < 1.0))

  def is_heavy_atom (self) :
    return (self.fofc_map > 3.0) or (self.anom_map > 3.0)

def analyze_waters (pdb_hierarchy, xray_structure, fmodel,
    distance_cutoff=4.0) :
  """
  Assess the properties of solvent atoms and return a list of objects
  containing their attributes, including local environment and electron
  density.
  """
  from mmtbx.real_space_correlation import extract_map_stats_for_single_atoms
  from cctbx import adptbx
  pdb_atoms = pdb_hierarchy.atoms()
  assert (not pdb_atoms.extract_i_seq().all_eq(0))
  unit_cell = xray_structure.unit_cell()
  pair_asu_table = xray_structure.pair_asu_table(
    distance_cutoff = distance_cutoff)
  asu_mappings = pair_asu_table.asu_mappings()
  asu_table = pair_asu_table.table()
  u_isos = xray_structure.extract_u_iso_or_u_equiv()
  occupancies = xray_structure.scatterers().extract_occupancies()
  sites_cart = xray_structure.sites_cart()
  sel_cache = pdb_hierarchy.atom_selection_cache()
  water_sel = sel_cache.selection("resname HOH and name O")
  map_stats = extract_map_stats_for_single_atoms(
    pdb_atoms=pdb_atoms,
    xray_structure=xray_structure,
    fmodel=fmodel,
    selection=water_sel)
  waters = []
  for i_seq, atom in enumerate(pdb_atoms) :
    if (water_sel[i_seq]) :
      asu_dict = asu_table[i_seq]
      site_cart = sites_cart[i_seq]
      water = water_properties(
        i_seq=i_seq,
        id_str=atom.id_str(),
        site_cart=site_cart,
        b_iso=adptbx.u_as_b(u_isos[i_seq]),
        occ=occupancies[i_seq],
        nearest_contact=None, #nearest_contact,
        map_cc=map_stats.two_fofc_ccs[i_seq],
        two_fofc_map=map_stats.two_fofc_values[i_seq],
        fofc_map=map_stats.fofc_values[i_seq],
        anom_map=map_stats.anom_values[i_seq],
        n_hbonds=None) # TODO
      waters.append(water)
  return waters
