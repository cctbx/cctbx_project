from iotbx.kriber import strudat
import iotbx.pdb
from iotbx.option_parser import iotbx_option_parser
from cctbx import restraints
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
import libtbx.itertbx
import math
import sys, os

class j_sym_dist:

  def __init__(self, j_sym, dist):
    adopt_init_args(self, locals())

def get_bond_asu_proxies(asu_mappings, distance_cutoff):
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=asu_mappings,
    distance_cutoff=distance_cutoff,
    minimal=0001)
  bond_dicts = [{} for i in xrange(asu_mappings.mappings().size())]
  for pair in pair_generator:
    bond_dicts[pair.i_seq].setdefault(pair.j_seq, []).append(
      j_sym_dist(pair.j_sym, pair.dist_sq**.5))
  return bond_dicts

class bond_sym_proxy:

  def __init__(self, i_seqs, rt_mx, distance_ideal=0, weight=0):
    adopt_init_args(self, locals())

def setup_bond_sym_proxies(structure, distance_cutoff, tolerance=1.e-6,
                           short_cut=00000, verbose=0):
  if (0 or verbose):
    print "distance_cutoff:", distance_cutoff
  distance_cutoff_plus = distance_cutoff * (1 + tolerance)
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff_plus)
  bond_dicts = get_bond_asu_proxies(
    asu_mappings=asu_mappings, distance_cutoff=distance_cutoff_plus)
  bond_sym_proxies = []
  for i_seq,bond_dict in enumerate(bond_dicts):
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0).inverse()
    so_i = asu_mappings.special_op(i_seq)
    ss_i = structure.site_symmetry(structure.scatterers()[i_seq].site)
    assert ss_i.special_op() == so_i
    for j_seq,j_sym_dists in bond_dict.items():
      if (0 or verbose):
        print "i_seq:", i_seq
      rt_mx_jis = [rt_mx_i_inv.multiply(
        asu_mappings.get_rt_mx(i_seq=j_seq, i_sym=jd.j_sym))
          for jd in j_sym_dists]
      so_j = asu_mappings.special_op(j_seq)
      is_processed = [-1] * len(j_sym_dists)
      for i in xrange(len(j_sym_dists)):
        jd_i = j_sym_dists[i]
        if (is_processed[i] >= 0): continue
        if (0 or verbose):
          print "  indep", i
        bond_sym_proxies.append(bond_sym_proxy(
          i_seqs=(i_seq, j_seq),
          rt_mx=rt_mx_jis[i]))
        for j in xrange(i+1,len(j_sym_dists)):
          jd_j = j_sym_dists[j]
          if (short_cut and abs(jd_i.dist - jd_j.dist) > 1.e-6): continue
          # the case j_seq == i_seq is not handled correctly here
          # see process_bond for correct treatment
          if (j_seq == i_seq and rt_mx_jis[i].inverse_cancel()==rt_mx_jis[j]):
            is_processed[j] = i
            if (0 or verbose):
              print "    symeq a:", j
            assert abs(jd_i.dist - jd_j.dist) <= 1.e-6
          else:
            ri = rt_mx_jis[i]
            rj = rt_mx_jis[j]
            for mi in ss_i.matrices():
              if (mi.multiply(rj).multiply(so_j) == ri.multiply(so_j)):
                is_processed[j] = i
                if (0 or verbose):
                  print "    symeq b:", j
                assert abs(jd_i.dist - jd_j.dist) <= 1.e-6
                break
      if (0 or verbose):
        print "  is_processed:", is_processed
  return bond_sym_proxies

def get_max_bond_length(structure, bond_sym_proxies):
  max_distance = 0
  sites_frac = structure.scatterers().extract_sites()
  unit_cell = structure.unit_cell()
  for proxy in bond_sym_proxies:
    max_distance = max(max_distance, unit_cell.distance(
      sites_frac[proxy.i_seqs[0]], proxy.rt_mx * sites_frac[proxy.i_seqs[1]]))
  return max_distance

def expand_bond_sym_proxies(structure, bond_sym_proxies, distance_cutoff=0,
                            tolerance=1.e-6, verbose=0):
  bond_asu_records = [{} for i_seq in xrange(structure.scatterers().size())]
  if (len(bond_sym_proxies) == 0):
    return None, bond_asu_records
  max_bond_length = get_max_bond_length(
    structure=structure, bond_sym_proxies=bond_sym_proxies)
  if (0 or verbose):
    print "max_bond_length:", max_bond_length
  assert max_bond_length > 0
  distance_cutoff = max(distance_cutoff, max_bond_length)
  if (0 or verbose):
    print "distance_cutoff:", distance_cutoff
  distance_cutoff_plus = distance_cutoff * (1 + tolerance)
  asu_mappings = structure.asu_mappings(buffer_thickness=distance_cutoff_plus)
  for sym_proxy in bond_sym_proxies:
    i_seq, j_seq = sym_proxy.i_seqs
    is_new = process_bond(
      asu_mappings=asu_mappings,
      site_symmetry_i=structure.site_symmetry(
        structure.scatterers()[i_seq].site),
      bond_asu_records=bond_asu_records,
      i_seq=i_seq,
      j_seq=j_seq,
      rt_mx_ji=sym_proxy.rt_mx,
      verbose=verbose)
    if (is_new and i_seq != j_seq):
      is_new = process_bond(
        asu_mappings=asu_mappings,
        site_symmetry_i=structure.site_symmetry(
          structure.scatterers()[j_seq].site),
        bond_asu_records=bond_asu_records,
        i_seq=j_seq,
        j_seq=i_seq,
        rt_mx_ji=sym_proxy.rt_mx.inverse_cancel(),
        verbose=verbose)
      assert is_new
  return asu_mappings, bond_asu_records

def process_bond(asu_mappings, site_symmetry_i, bond_asu_records,
                 i_seq, j_seq, rt_mx_ji, verbose=0):
  rt_mx_i = asu_mappings.get_rt_mx(i_seq=i_seq, i_sym=0)
  rt_mx_j = rt_mx_i.multiply(rt_mx_ji)
  j_sym = asu_mappings.find_i_sym(i_seq=j_seq, rt_mx=rt_mx_j)
  assert j_sym >= 0
  j_sym_groups = bond_asu_records[i_seq].setdefault(j_seq, [])
  for j_sym_group in j_sym_groups:
    if (j_sym in j_sym_group):
      return 00000
  j_syms = [j_sym]
  j_sym_groups.append(j_syms)
  if (0 or verbose):
    print "primary:     i_seq, j_seq, j_sym", i_seq, j_seq, j_sym
  for mi in site_symmetry_i.matrices():
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

def check_sym_equiv(structure, asu_mappings, bond_asu_records, weak=00000):
  unit_cell = structure.unit_cell()
  sites_frac = structure.scatterers().extract_sites()
  for i_seq,records in enumerate(bond_asu_records):
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

def check_connectivities(bond_asu_records, connectivities, verbose=0):
  n_mismatches = 0
  for records,connectivity in zip(bond_asu_records, connectivities):
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
  bond_sym_proxies = setup_bond_sym_proxies(
    structure=structure,
    distance_cutoff=distance_cutoff,
    verbose=verbose)
  if (0 or verbose):
    print "bond_sym_proxies.size():", len(bond_sym_proxies)
  for dist_cutoff in [distance_cutoff, 0]:
    asu_mappings, bond_asu_records = expand_bond_sym_proxies(
      structure=structure,
      bond_sym_proxies=bond_sym_proxies,
      distance_cutoff=dist_cutoff,
      verbose=verbose)
    if (connectivities is not None):
      check_connectivities(bond_asu_records, connectivities, verbose)
    if (asu_mappings is not None):
      check_sym_equiv(
        structure=structure,
        asu_mappings=asu_mappings,
        bond_asu_records=bond_asu_records,
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
