
from __future__ import absolute_import, division, print_function
from libtbx import slots_getstate_setstate, \
    slots_getstate_setstate_default_initializer
import sys

class conformer(slots_getstate_setstate):
  __slots__ = ["atom_groups", "label", "selection"]
  def __init__(self, atom_group):
    from scitbx.array_family import flex
    assert (atom_group.altloc != '')
    self.atom_groups = [ atom_group ]
    self.selection = flex.size_t()
    self.label = atom_group.altloc

  def append_atom_group(self, atom_group):
    assert (atom_group.altloc == self.labl)
    self.atom_groups.append(atom_group)
    self.selection.extend(atom_group.atoms().extract_i_seq())

class conformer_group(slots_getstate_setstate):
  __slots__ = ["residue_groups"]
  def __init__(self, residue_groups):
    self.residue_groups = residue_groups
    self.conformers = []
    for i_res, residue_group in enumerate(residue_groups):
      for i_group, atom_group in enumerate(residue_group.atom_groups()):
        if (i_res == 0):
          conf = conformer(atom_group)
          self.conformers.append(conf)
        else :
          self.conformers[i_group].append_atom_group(atom_group)

class nb_clash(slots_getstate_setstate_default_initializer):
  __slots__ = ["i_seq", "j_seq", "id_str_i", "id_str_j", "overlap"]

  def show(self, out=sys.stdout, prefix=""):
    print(prefix+"%s  %s  : %.3f" % (self.id_str_i, self.id_str_j,
      self.overlap), file=out)

def get_sorted_clashes(
    pdb_atoms,
    xray_structure,
    selection,
    geometry_restraints_manager,
    clash_min=0.2,
    out=sys.stdout):
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  pair_proxies = geometry_restraints_manager.pair_proxies(
    sites_cart=sites_cart,
    site_labels=site_labels)
  asu_proxies = pair_proxies.nonbonded_proxies.asu
  deltas = pair_proxies.nonbonded_proxies.deltas(sites_cart=sites_cart)
  clashes = []
  for proxy, delta in zip(asu_proxies, deltas):
    if (not selection[proxy.i_seq]) or (not selection[proxy.j_seq]):
      continue
    overlap = proxy.vdw_distance - delta
    if (overlap >= clash_min):
      #print proxy.vdw_distance, delta, clash_min
      clash = nb_clash(
        i_seq=proxy.i_seq,
        j_seq=proxy.j_seq,
        id_str_i=site_labels[proxy.i_seq], #pdb_atoms[proxy.i_seq].id_str(),
        id_str_j=site_labels[proxy.j_seq], #pdb_atoms[proxy.j_seq].id_str(),
        overlap=overlap)
      clashes.append(clash)
  clashes.sort(lambda a,b: cmp(b.overlap, a.overlap))
  return clashes

def show_altloc_clashes(
    pdb_hierarchy,
    xray_structure,
    geometry_restraints_manager,
    clash_min=0.2,
    out=sys.stdout):
  sel_cache = pdb_hierarchy.atom_selection_cache()
  alt_conf_sel = sel_cache.selection("not altloc ' '")
  clashes = get_sorted_clashes(
    pdb_atoms=pdb_hierarchy.atoms(),
    xray_structure=xray_structure,
    selection=alt_conf_sel,
    geometry_restraints_manager=geometry_restraints_manager,
    clash_min=clash_min,
    out=out)
  if (len(clashes) > 0):
    print("Clashing atoms in alternate conformations:", file=out)
    for clash in clashes :
      clash.show(out=out, prefix="  ")
  return clashes
