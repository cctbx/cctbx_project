# -*- coding: utf-8; py-indent-offset: 2 -*-
from __future__ import division
from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
from libtbx.utils import null_out
from math import sqrt
import sys

def anonymize_ions (hierarchy, log=sys.stdout) :
  """
  Convert any elemental ions in the PDB hierarchy to water, resetting the
  occupancy and scaling the B-factor.  The atom segids will be set to the old
  resname.  NOTE: this does not change the corresponding scatterer in the xray
  structure, but a new xray structure can be obtained by calling
  hierarchy.extract_xray_structure(crystal_symmetry).
  """
  from cctbx.eltbx import sasaki, chemical_elements
  import mmtbx.ions
  ion_resnames = set(chemical_elements.proper_upper_list())
  for resname in mmtbx.ions.server.params["_lib_charge.resname"]:
    if resname not in WATER_RES_NAMES:
      ion_resnames.add(resname)
  n_converted = 0
  hierarchy = hierarchy.deep_copy()
  for model in hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        for atom_group in residue_group.atom_groups() :
          if (atom_group.resname.strip() in ion_resnames) :
            atoms = atom_group.atoms()
            id_strs = []
            for atom in atoms :
              elem = atom.element.strip()
              atomic_number = sasaki.table(elem).atomic_number()
              id_strs.append(atom.id_str())
              atom.segid = atom_group.resname
              atom.name = " O  "
              atom.element = "O"
              atom.charge = ""
              atom.occupancy = 1.0
              atom.b = atom.b * (10 / atomic_number)
            atom_group.resname = "HOH"
            for atom, id_str in zip(atoms, id_strs) :
              print >> log, "%s --> %s, B-iso = %.2f" % (id_str, atom.id_str(),
                atom.b)
              n_converted += 1
  return hierarchy, n_converted

def sort_atoms_permutation (pdb_atoms, xray_structure) :
  assert (pdb_atoms.size() == xray_structure.scatterers().size())
  from scitbx.array_family import flex
  pdb_atoms.reset_i_seq()
  atoms_sorted = sorted(pdb_atoms, cmp_atom)
  sele = flex.size_t([ atom.i_seq for atom in atoms_sorted ])
  return sele

# compare mass, then occupancy, then B_iso
def cmp_atom (a, b) :
  from cctbx.eltbx import sasaki
  mass_a = sasaki.table(a.element.strip().upper()).atomic_number()
  mass_b = sasaki.table(b.element.strip().upper()).atomic_number()
  if (mass_a == mass_b) :
    if (a.occ == b.occ) :
      return cmp(b.b, a.b)
    else :
      return cmp(b.occ, a.occ)
  else :
    return cmp(mass_b, mass_a)

def collect_ions (pdb_hierarchy) :
  from cctbx.eltbx import chemical_elements
  elements = chemical_elements.proper_upper_list()
  ions = []
  for model in pdb_hierarchy.models() :
    for chain in model.chains() :
      for residue_group in chain.residue_groups() :
        for atom_group in residue_group.atom_groups() :
          if (atom_group.resname.strip() in elements) :
            atoms = atom_group.atoms()
            assert (len(atoms) == 1)
            for atom in atoms:
              ions.append(atom)
  return ions

# TODO add test
def compare_ions (hierarchy, reference_hierarchy, reference_xrs,
    distance_cutoff=2.0, log=None, ignore_elements=(), only_elements=(),
    sel_str_base="segid ION") :
  if (log is None) : log = null_out()
  sel_cache = hierarchy.atom_selection_cache()
  sel_str = sel_str_base
  if (len(only_elements) > 0) :
    sel_str += " and (%s)" % " or ".join(
      [ "element %s" % e for e in only_elements ])
  elif (len(ignore_elements) > 0) :
    sel_str += " and not (%s)" % " or ".join(
      [ "element %s" % e for e in ignore_elements ])
  ion_isel = sel_cache.selection(sel_str).iselection()
  if (len(ion_isel) == 0) :
    return [], []
  from cctbx import crystal
  pdb_atoms = hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  ions = pdb_atoms.select(ion_isel)
  asu_mappings = reference_xrs.asu_mappings(
    buffer_thickness=distance_cutoff+0.1)
  unit_cell = reference_xrs.unit_cell()
  sites_cart = ions.extract_xyz()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  asu_mappings.process_sites_frac(sites_frac,
    min_distance_sym_equiv = reference_xrs.min_distance_sym_equiv())
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings = asu_mappings,
    distance_cutoff = distance_cutoff)
  reference_atoms = reference_hierarchy.atoms()
  n_xray = reference_xrs.scatterers().size()
  ion_ref_i_seqs = []
  for k in range(len(ions)) :
    ion_ref_i_seqs.append([])
  for pair in pair_generator:
    if ((pair.i_seq < n_xray and pair.j_seq < n_xray) or
        (pair.i_seq >= n_xray and pair.j_seq >= n_xray)) :
      continue
    if (pair.i_seq < n_xray) :
      ion_seq, ref_seq = pair.j_seq, pair.i_seq
    else :
      ion_seq, ref_seq = pair.i_seq, pair.j_seq
    site_frac = sites_frac[ion_seq - n_xray]
    dxyz = sqrt(pair.dist_sq)
    j_seq = ion_seq - n_xray
    ion_ref_i_seqs[j_seq].append(ref_seq)
  # FIXME better filtering needed - right now we risk double-counting ions in
  # the reference model, although I haven't found a case of this yet
  matched = []
  missing = []
  for i_seq, ref_i_seqs in enumerate(ion_ref_i_seqs) :
    ion = ions[i_seq]
    if (len(ref_i_seqs) == 0) :
      print >> log, "No match for %s" % ion.id_str()
      missing.append(ion.id_str())
    else :
      ref_ions = []
      for i_seq in ref_i_seqs :
        ref_atom = reference_atoms[i_seq]
        if (ion.element.upper() == ref_atom.element.upper()) :
          ref_ions.append(ref_atom.id_str())
      if (len(ref_ions) >= 1) :
        matched.append(ion.id_str())
        if (len(ref_ions) > 1) :
          print >> log, "Multiple matches for %s:" % ion.id_str()
          for ref_ion in ref_ions :
            print >> log, "  %s" % ref_ion
        else :
          print >> log, "Ion %s matches %s" % (ion.id_str(),
            ref_ions[0])
      else :
        print >> log, "No match for %s" % ion.id_str()
        missing.append(ion.id_str())
  return matched, missing
