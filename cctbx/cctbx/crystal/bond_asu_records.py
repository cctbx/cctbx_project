from iotbx.kriber import strudat
import iotbx.pdb
from iotbx.option_parser import iotbx_option_parser
from cctbx import crystal
import libtbx.itertbx
import math
import sys, os

def start_bond_table(structure, distance_cutoff, verbose=0):
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
  bond_table = [{} for i_seq in xrange(structure.scatterers().size())]
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    minimal=0001)
  for pair in pair_generator:
    rt_mx_i = asu_mappings.get_rt_mx(i_seq=pair.i_seq, i_sym=0)
    rt_mx_j = asu_mappings.get_rt_mx(i_seq=pair.j_seq, i_sym=pair.j_sym)
    rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
    is_new = process_bond(
      asu_mappings=asu_mappings,
      bond_table=bond_table,
      i_seq=pair.i_seq,
      j_seq=pair.j_seq,
      rt_mx_ji=rt_mx_ji,
      verbose=verbose)
    if (is_new and pair.i_seq != pair.j_seq):
      is_new = process_bond(
        asu_mappings=asu_mappings,
        bond_table=bond_table,
        i_seq=pair.j_seq,
        j_seq=pair.i_seq,
        rt_mx_ji=rt_mx_ji.inverse_cancel(),
        verbose=verbose)
      assert is_new
  return asu_mappings, bond_table

def process_bond(asu_mappings, bond_table,
                 i_seq, j_seq, rt_mx_ji, verbose=0):
  rt_mx_i = asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0)
  rt_mx_j = rt_mx_i.multiply(rt_mx_ji)
  j_sym = asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j)
  assert j_sym >= 0
  j_sym_groups = bond_table[i_seq].setdefault(j_seq, [])
  for j_sym_group in j_sym_groups:
    if (j_sym in j_sym_group):
      return 00000
  j_syms = [j_sym]
  j_sym_groups.append(j_syms)
  if (0 or verbose):
    print "primary:     i_seq, j_seq, j_sym", i_seq, j_seq, j_sym
  for mi in asu_mappings.site_symmetry_table().get(i_seq).matrices():
    if (i_seq == j_seq):
      rt_mx_j_eq = rt_mx_i.multiply(rt_mx_ji.multiply(mi).inverse_cancel())
      j_sym_eq = asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j_eq)
      assert j_sym_eq >= 0
      if (not j_sym_eq in j_syms):
        j_syms.append(j_sym_eq)
        if (0 or verbose):
          print "    equiv a: i_seq, j_seq, j_sym", i_seq, j_seq, j_sym_eq
    rt_mx_j_eq = rt_mx_i.multiply(mi.multiply(rt_mx_ji))
    j_sym_eq = asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j_eq)
    assert j_sym_eq >= 0
    if (not j_sym_eq in j_syms):
      j_syms.append(j_sym_eq)
      if (0 or verbose):
        print "    equiv b: i_seq, j_seq, j_sym", i_seq, j_seq, j_sym_eq
  return 0001

def is_sym_equiv_interaction_simple(unit_cell,
                                    i_seq,
                                    site_frac_i,
                                    j_seq,
                                    site_frac_j,
                                    special_op_j,
                                    rt_mx_ji_1,
                                    rt_mx_ji_2):
  f = unit_cell.shortest_vector_sq()**.5*.1
  trial_shifts = [f*x for x in [math.sqrt(2),math.sqrt(3),math.sqrt(5)]]
  frac = unit_cell.fractionalize
  orth = unit_cell.orthogonalize
  dist = unit_cell.distance
  for shifts in [[0,0,0], trial_shifts]:
    site_j_mod = special_op_j * frac([x+s
      for x,s in zip(orth(site_frac_j),shifts)])
    if (shifts == [0,0,0] or j_seq != i_seq):
      site_i_mod = site_frac_i
    else:
      site_i_mod = site_j_mod
    d1 = dist(rt_mx_ji_1 * site_j_mod, site_i_mod)
    d2 = dist(rt_mx_ji_2 * site_j_mod, site_i_mod)
    if (shifts == [0,0,0]):
      if (abs(d1-d2) >= 1.e-3):
        return 00000
  return abs(d1-d2) < 1.e-3

def check_sym_equiv(structure, asu_mappings, bond_table, weak=00000):
  unit_cell = structure.unit_cell()
  sites_frac = structure.scatterers().extract_sites()
  for i_seq,records in enumerate(bond_table):
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
    for j_seq,j_sym_groups in records.items():
      i_group_rt_mx_jis = []
      for i_group,j_sym_group in enumerate(j_sym_groups):
        for j_sym in j_sym_group:
          rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
          i_group_rt_mx_jis.append((i_group,rt_mx_ji))
      for gi,ri in i_group_rt_mx_jis:
        for gj,rj in i_group_rt_mx_jis:
          is_sym_equiv = is_sym_equiv_interaction_simple(
            unit_cell=unit_cell,
            i_seq=i_seq,
            site_frac_i=sites_frac[i_seq],
            j_seq=j_seq,
            site_frac_j=sites_frac[j_seq],
            special_op_j=asu_mappings.special_op(j_seq),
            rt_mx_ji_1=ri,
            rt_mx_ji_2=rj)
          if (is_sym_equiv):
            if (not weak): assert gi == gj
          else:
            assert gi != gj

def check_connectivities(bond_table, connectivities, verbose=0):
  n_mismatches = 0
  for records,connectivity in zip(bond_table, connectivities):
    n = 0
    for j_seq,j_sym_groups in records.items():
      for j_sym_group in j_sym_groups:
        n += len(j_sym_group)
    if (0 or verbose):
      print "n, connectivity:", n, connectivity
    assert n == connectivity

def exercise(
      structure,
      distance_cutoff,
      connectivities=None,
      weak_check_sym_equiv=00000,
      verbose=0):
  if (0 or verbose):
    print "distance_cutoff:", distance_cutoff
  asu_mappings, bond_table = start_bond_table(
    structure=structure,
    distance_cutoff=distance_cutoff,
    verbose=verbose)
  if (connectivities is not None):
    check_connectivities(bond_table, connectivities, verbose)
  check_sym_equiv(
    structure=structure,
    asu_mappings=asu_mappings,
    bond_table=bond_table,
    weak=weak_check_sym_equiv)

def run():
  command_line = (iotbx_option_parser(
    usage="python neighbors.py [options] studat_file [...]",
    description="Example: python neighbors.py strudat --tag=SOD")
    .enable_symmetry_comprehensive()
    .option(None, "--tag",
      action="store",
      type="string",
      dest="tag",
      help="tag as it appears in the strudat file")
    .option(None, "--distance_cutoff",
      action="store",
      type="float",
      dest="distance_cutoff",
      metavar="FLOAT")
    .option(None, "--verbose",
      action="store_true",
      dest="verbose",
      help="produce output")
  ).process()
  default_distance_cutoff=3.5
  file_names = command_line.args
  if (len(file_names) == 0):
    regression_misc = os.path.join(
      os.environ["LIBTBX_DIST_ROOT"], "regression", "misc")
    for file_name in ["strudat_zeolite_atlas", "strudat_special_bonds"]:
      path = os.path.join(regression_misc, file_name)
      if (os.path.isfile(path)):
        file_names.append(path)
  for file_name in file_names:
    try:
      strudat_entries = strudat.read_all_entries(open(file_name))
    except:
      strudat_entries = None
    if (strudat_entries is not None and len(strudat_entries.entries) > 0):
      for entry in strudat_entries.entries:
        if (    command_line.options.tag is not None
            and command_line.options.tag != entry.tag):
          continue
        if (0 or command_line.options.verbose):
          print "strudat tag:", entry.tag
        structure = entry.as_xray_structure()
        if (0 or command_line.options.verbose):
          structure.show_summary().show_scatterers()
        if (command_line.options.distance_cutoff is not None):
          distance_cutoff = command_line.options.distance_cutoff
        elif (entry.title.startswith("cutoff")):
          distance_cutoff = float(entry.title.split()[1])
        else:
          distance_cutoff = default_distance_cutoff
        connectivities = []
        weak_check_sym_equiv = (
          entry.reference.find("weak_check_sym_equiv") >= 0)
        for atom in entry.atoms:
          if (atom.connectivity is None):
            if (len(connectivities) != 0):
              raise AssertionError(
                "Tag %s: some atoms are missing the bond count." % entry.tag)
            connectivities = None
            break
          connectivities.append(atom.connectivity)
        exercise(
          structure=structure,
          distance_cutoff=distance_cutoff,
          connectivities=connectivities,
          weak_check_sym_equiv=weak_check_sym_equiv,
          verbose=command_line.options.verbose)
        if (0 or command_line.options.verbose):
          print
    else:
      try:
        structure = iotbx.pdb.as_xray_structure(
          file_name=file_name,
          crystal_symmetry=command_line.symmetry)
      except:
        raise RuntimeError("Coordinate file %s: unknown format." % file_name)
      else:
        if (0 or command_line.options.verbose):
          structure.show_summary().show_scatterers()
        exercise(
          structure=structure,
          distance_cutoff=distance_cutoff,
          verbose=command_line.options.verbose)
  print "OK"

if (__name__ == "__main__"):
  run()
