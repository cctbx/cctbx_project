from iotbx.kriber import strudat
from cctbx import crystal
from scitbx.python_utils.misc import adopt_init_args
import math
import sys, os

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

def check_sym_equiv(structure, bond_asu_table, weak=00000):
  unit_cell = structure.unit_cell()
  asu_mappings = bond_asu_table.asu_mappings()
  sites_frac = structure.scatterers().extract_sites()
  for i_seq,records in enumerate(bond_asu_table.table()):
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

def check_connectivities(bond_asu_table, connectivities, verbose=0):
  n_mismatches = 0
  for records,connectivity in zip(bond_asu_table.table(), connectivities):
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
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff)
  for i_pass in xrange(2):
    if (i_pass == 0):
      bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_all_pairs(
        distance_cutoff=distance_cutoff)
    else:
      bond_sym_table = bond_asu_table.extract_pair_sym_table()
      bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      bond_asu_table.add_pair_sym_table(
        sym_table=bond_sym_table)
    if (connectivities is not None):
      check_connectivities(bond_asu_table, connectivities, verbose)
    check_sym_equiv(
      structure=structure,
      bond_asu_table=bond_asu_table,
      weak=weak_check_sym_equiv)

def run():
  verbose = "--Verbose" in sys.argv[1:]
  default_distance_cutoff=3.5
  regression_misc = os.path.join(
    os.environ["LIBTBX_DIST_ROOT"], "regression", "misc")
  file_names = []
  for file_name in ["strudat_zeolite_atlas", "strudat_special_bonds"]:
    path = os.path.join(regression_misc, file_name)
    if (os.path.isfile(path)):
      file_names.append(path)
  for file_name in file_names:
    strudat_entries = strudat.read_all_entries(open(file_name))
    for entry in strudat_entries.entries:
      if (0 or verbose):
        print "strudat tag:", entry.tag
      structure = entry.as_xray_structure()
      if (0 or verbose):
        structure.show_summary().show_scatterers()
      if (entry.title.startswith("cutoff")):
        distance_cutoff = float(entry.title.split()[1])
      else:
        distance_cutoff = default_distance_cutoff
      weak_check_sym_equiv = (
        entry.reference.find("weak_check_sym_equiv") >= 0)
      connectivities = entry.connectivities(all_or_nothing=0001)
      exercise(
        structure=structure,
        distance_cutoff=distance_cutoff,
        connectivities=connectivities,
        weak_check_sym_equiv=weak_check_sym_equiv,
        verbose=verbose)
      if (0 or verbose):
        print
  print "OK"

if (__name__ == "__main__"):
  run()
