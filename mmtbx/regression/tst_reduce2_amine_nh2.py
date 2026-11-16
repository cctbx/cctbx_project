"""reduce2: amine NH2 (A17). The riding H code places a primary amine whose
dictionary has a neutral NH2 (N with 2 H and one sp3 C, N or O neighbour, no
plane restraint) pyramidal at its dictionary angles. Before: a planar pair
(H-N-H ~141), or no H at all when the dictionary H-N-H is 107-111 deg.
Fixtures: CCD ideal heavy atoms in a P1 30 A box. GeoStd PMP and the CCD
fallback of A71 (no GeoStd file) get both H; the NH2 torsions of HSM (one
restraint on each H, stored atom-swapped by the proxy registry) are met, and
a fresh riding idealization (refinement's first step) moves nothing; the
result does not depend on a rigid-body move of the input (EDN); H and D of
the neutron options are placed alike (PUT); NME and HOA, with no anchor for
riding H, get the dictionary H-N-H too. The H follow the dictionary's
protonation state.
"""
from __future__ import absolute_import, division, print_function
import math
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx import matrix
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from cctbx import geometry_restraints
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen

CRYST1 = "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"

def pdb(resname, atoms):
  lines = [CRYST1.strip()]
  for i, (name, x, y, z) in enumerate(atoms):
    element = 'PT' if name == 'PT' else name[0]
    lines.append("HETATM%5d %-4s %3s A   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "          %2s" % (i + 1, (name if len(name) == 4 or len(element) == 2
      else ' ' + name), resname, x, y, z, element))
  return "\n".join(lines) + "\nEND\n"

pdb_pmp = pdb('PMP', [
  ('N1', 16.962, 14.767, 16.724), ('C2', 16.168, 14.918, 17.766),
  ('C2A', 16.739, 14.829, 19.158), ('C3', 14.813, 15.149, 17.594),
  ('O3', 13.998, 15.303, 18.670), ('C4', 14.300, 15.232, 16.303),
  ('C4A', 12.834, 15.489, 16.070), ('N4A', 12.115, 14.210, 16.008),
  ('C5', 15.169, 15.068, 15.237), ('C6', 16.508, 14.834, 15.489),
  ('C5A', 14.662, 15.143, 13.821), ('O4P', 15.747, 14.943, 12.914),
  ('P', 15.136, 15.034, 11.428), ('O1P', 14.527, 16.367, 11.228),
  ('O2P', 16.311, 14.816, 10.350), ('O3P', 14.014, 13.898, 11.237)])

pdb_a71 = pdb('A71', [
  ('N1', 11.062, 15.675, 14.364), ('C2', 12.133, 15.638, 13.357),
  ('C3', 13.299, 14.734, 13.715), ('C4', 14.433, 14.718, 12.685),
  ('C5', 15.804, 14.264, 13.194), ('C6', 16.755, 15.368, 13.675),
  ('C7', 17.781, 14.933, 14.712), ('N8', 17.554, 15.439, 16.078),
  ('N9', 16.512, 13.116, 17.272), ('N10', 14.978, 16.785, 16.505),
  ('PT', 15.740, 14.936, 16.837), ('N11', 13.944, 14.399, 17.604)])

pdb_hsm = pdb('HSM', [
  ('N', 14.997, 15.142, 18.490), ('CA', 14.719, 15.442, 17.079),
  ('CB', 15.476, 14.456, 16.188), ('CG', 15.190, 14.764, 14.742),
  ('ND1', 13.985, 14.745, 14.156), ('CD2', 16.102, 15.099, 13.809),
  ('CE1', 14.113, 15.063, 12.897), ('NE2', 15.419, 15.290, 12.639)])

pdb_edn = pdb('EDN', [
  ('NAA', 16.861, 14.833, 15.000), ('CAB', 15.562, 15.519, 15.000),
  ('CAC', 14.438, 14.481, 15.000), ('NAD', 13.139, 15.167, 15.000)])

pdb_put = pdb('PUT', [
  ('N2', 14.750, 15.000, 18.089), ('C4', 15.535, 15.000, 16.849),
  ('C3', 14.591, 15.000, 15.646), ('C2', 15.409, 15.000, 14.354),
  ('C1', 14.465, 15.000, 13.151), ('N1', 15.250, 15.000, 11.911)])

# NME (methylamine) and HOA (hydroxylamine): two heavy atoms, no anchor for
# riding H (placement's place_anchorless_h())
pdb_nme = pdb('NME', [('N', 20.712, 20.012, 19.934), ('C', 19.258, 19.979,
  20.064)])
pdb_hoa = pdb('HOA', [('N', 20.740, 20.031, 19.924), ('O', 19.284, 19.933,
  20.030)])

def run(prefix, pdb_str, args=()):
  fn, out = prefix + '_in.pdb', prefix + '_out.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  if os.path.exists(out): os.remove(out)
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[fn,
    'output.filename=' + out, 'output.description_file_name=%s.txt' % prefix,
    'output.overwrite=True'] + list(args))
  return out

def sites(h, altloc=''):
  return dict((a.name.strip(), matrix.col(a.xyz)) for a in h.atoms()
    if a.parent().altloc in ('', altloc))

def angle(a, b, c):
  return math.degrees((a - b).angle(c - b))

def check_nh2(s, n, x, hs, hnh, tol=0.2):
  """two H, the dictionary H-N-H, N out of the X-H-H plane"""
  assert all(h in s for h in hs), (n, hs, sorted(s))
  h1, h2 = [s[h] for h in hs]
  value = angle(h1, s[n], h2)
  assert abs(value - hnh) < tol, (n, value, hnh)
  normal = (h1 - s[x]).cross(h2 - s[x]).normalize()
  height = abs((s[n] - s[x]).dot(normal))
  assert height > 0.3, (n, height)

def fresh(out):
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(out), log=null_out())
  reduce_hydrogen.add_ccd_restraint_objects(m)
  m.process(make_restraints=True)
  return m

def test_refused_dictionaries():
  # GeoStd PMP N4A (H-N-H 107.2) and the CCD fallback of A71 N1 (107.5):
  # inside the propeller window, so no H before
  for code, pdb_str, n, x, hs, hnh in (
      ('pmp', pdb_pmp, 'N4A', 'C4A', ('HNA1', 'HNA2'), 107.2),
      ('a71', pdb_a71, 'N1', 'C2', ('HN11', 'HN12'), 107.5)):
    out = run('tst_reduce2_amine_nh2_' + code, pdb_str)
    check_nh2(sites(iotbx.pdb.input(out).construct_hierarchy()), n, x, hs,
      hnh)

def test_torsions_and_rerun():
  # HSM: each NH2 H meets its own torsion restraint; idealizing a fresh
  # model of the output (refinement) moves nothing (before: 1.45 A, and the
  # torsion of H1 55 deg off)
  out = run('tst_reduce2_amine_nh2_hsm', pdb_hsm)
  m = fresh(out)
  atoms = m.get_hierarchy().atoms()
  sc = m.get_sites_cart()
  check_nh2(sites(m.get_hierarchy()), 'N', 'CA', ('H1', 'H2'), 106.1)
  nh2_h = set(a.i_seq for a in atoms if a.name.strip() in ('H1', 'H2'))
  n_torsions = 0
  for p in m.get_restraints_manager().geometry.dihedral_proxies:
    if nh2_h.intersection(p.i_seqs):
      delta = geometry_restraints.dihedral(sites_cart=sc, proxy=p).delta
      assert abs(delta) < 1.0, ([atoms[i].name for i in p.i_seqs], delta)
      n_torsions += 1
  assert n_torsions == 2, n_torsions
  m.setup_riding_h_manager(idealize=False)
  sc2 = sc.deep_copy()
  m.riding_h_manager.idealize_riding_h_positions(sites_cart=sc2)
  d = flex.sqrt((sc2 - sc).dot())
  assert flex.max(d) < 0.01, flex.max(d)

def test_rigid_body_move():
  # EDN rotated and shifted: the same H by name (before: 1.45 A apart)
  out = run('tst_reduce2_amine_nh2_edn', pdb_edn)
  ref = sites(iotbx.pdb.input(out).construct_hierarchy())
  check_nh2(ref, 'NAA', 'CAB', ('HAA1', 'HAA2'), 105.5)
  check_nh2(ref, 'NAD', 'CAC', ('HAD1', 'HAD2'), 105.5)
  rot = matrix.col((1, -2, 0.5)).normalize(
    ).axis_and_angle_as_r3_rotation_matrix(125, deg=True)
  center, shift = matrix.col((15, 15, 15)), matrix.col((1.3, -0.7, 2.1))
  h = iotbx.pdb.input(source_info=None, lines=pdb_edn).construct_hierarchy()
  for a in h.atoms():
    a.xyz = tuple(rot * (matrix.col(a.xyz) - center) + center + shift)
  out = run('tst_reduce2_amine_nh2_edn_moved', CRYST1 + h.as_pdb_string())
  moved = sites(iotbx.pdb.input(out).construct_hierarchy())
  for name in ('HAA1', 'HAA2', 'HAD1', 'HAD2'):
    back = rot.transpose() * (moved[name] - center - shift) + center
    # PDB coordinates: 0.001 A rounding of the moved input
    assert (back - ref[name]).length() < 0.005, (name,
      (back - ref[name]).length())

def test_neutron_options():
  # hd_and_h: H (conformer A) and D (B) at the same sites; all_d: D only
  out = run('tst_reduce2_amine_nh2_put_hd', pdb_put,
    ['neutron_option=hd_and_h'])
  h = iotbx.pdb.input(out).construct_hierarchy()
  a, b = sites(h, 'A'), sites(h, 'B')
  for n, x, hs in (('N1', 'C1', ('HN11', 'HN12')),
                   ('N2', 'C4', ('HN21', 'HN22'))):
    ds = tuple('D' + name[1:] for name in hs)
    check_nh2(a, n, x, hs, 105.45)
    check_nh2(b, n, x, ds, 105.45)
    for hn, dn in zip(hs, ds):
      assert (a[hn] - b[dn]).length() < 1.e-3, (hn, dn)
  out = run('tst_reduce2_amine_nh2_put_d', pdb_put, ['neutron_option=all_d'])
  s = sites(iotbx.pdb.input(out).construct_hierarchy())
  check_nh2(s, 'N1', 'C1', ('DN11', 'DN12'), 105.45)
  check_nh2(s, 'N2', 'C4', ('DN21', 'DN22'), 105.45)

def test_unanchored():
  # no heavy atom besides N and X: pyramidal at the dictionary H-N-H, the
  # orientation a fixed convention (before: 120 deg apart about X-N, H-N-H
  # 109.6 for NME, 113.5 for HOA)
  for code, pdb_str, x, hnh in (('nme', pdb_nme, 'C', 105.4),
                                ('hoa', pdb_hoa, 'O', 105.3)):
    out = run('tst_reduce2_amine_nh2_' + code, pdb_str)
    check_nh2(sites(iotbx.pdb.input(out).construct_hierarchy()), 'N', x,
      ('HN1', 'HN2'), hnh)

if __name__ == '__main__':
  test_refused_dictionaries()
  test_torsions_and_rerun()
  test_rigid_body_move()
  test_neutron_options()
  test_unanchored()
  print('OK')
