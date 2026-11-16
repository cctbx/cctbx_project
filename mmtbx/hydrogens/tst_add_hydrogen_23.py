"""
add_link_h_restraints: the amide H of a peptide-linked N without CA. 9eor HT7
is a beta-amino acid, N bonded to CB: riding needs a C(i-1)-N-X angle, and the
H lies in the C, N, X plane at the peptide N-H length, which riding uses as
its target on a first run and on a keep_existing_H rerun. Fragment: 9eor G 2-3.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.cif
import iotbx.pdb
import libtbx.load_env
import mmtbx.model
from cctbx import geometry_restraints
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx import matrix
from scitbx.array_family import flex
from mmtbx.hydrogens import reduce_hydrogen

pdb_str = """\
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
HETATM14142  N   DNE G   2     -19.536 -15.401  62.306  1.00 19.63           N
HETATM14143  CA  DNE G   2     -20.187 -14.689  63.387  1.00 20.76           C
HETATM14144  C   DNE G   2     -21.461 -15.483  63.644  1.00 19.22           C
HETATM14145  O   DNE G   2     -21.379 -16.625  64.070  1.00 20.98           O
HETATM14146  CB  DNE G   2     -19.292 -14.498  64.611  1.00 21.84           C
HETATM14147  CG  DNE G   2     -17.813 -14.327  64.273  1.00 23.05           C
HETATM14148  CD  DNE G   2     -17.239 -12.964  64.370  1.00 23.31           C
HETATM14149  CE  DNE G   2     -16.122 -12.745  63.372  1.00 25.60           C
HETATM14150  N   HT7 G   3     -22.623 -14.925  63.307  1.00 17.64           N
HETATM14151  CA  HT7 G   3     -24.843 -15.694  62.547  1.00 21.80           C
HETATM14152  C   HT7 G   3     -24.240 -16.569  61.471  1.00 24.00           C
HETATM14153  O   HT7 G   3     -23.164 -17.136  61.638  1.00 27.40           O
HETATM14154  CB  HT7 G   3     -23.893 -15.480  63.744  1.00 18.87           C
HETATM14155  CG  HT7 G   3     -24.479 -14.557  64.804  1.00 16.66           C
HETATM14156  CD  HT7 G   3     -23.598 -14.420  66.002  1.00 18.30           C
HETATM14157  CE1 HT7 G   3     -23.619 -15.222  67.111  1.00 19.71           C
HETATM14158  CE2 HT7 G   3     -22.550 -13.456  66.245  1.00 18.01           C
HETATM14159  NZ1 HT7 G   3     -22.684 -14.812  68.030  1.00 19.29           N
HETATM14160  CZ2 HT7 G   3     -22.029 -13.716  67.535  1.00 18.92           C
HETATM14161  CZ3 HT7 G   3     -22.008 -12.396  65.526  1.00 18.29           C
HETATM14162  CH2 HT7 G   3     -20.986 -12.984  68.092  1.00 18.64           C
HETATM14163  CH3 HT7 G   3     -20.988 -11.657  66.089  1.00 19.09           C
HETATM14164  CT2 HT7 G   3     -20.469 -11.960  67.343  1.00 18.11           C
END
"""

def ht7_cif(extra_bond=None, resname='HT7'):
  '''The geostd HT7 dictionary as a user CIF, called resname; extra_bond: a
  line to add.'''
  fn = libtbx.env.find_in_repositories(
    relative_path='chem_data/geostd/h/data_HT7.cif', test=os.path.isfile)
  with open(fn) as f:
    text = f.read().replace('HT7', resname)
  if extra_bond is not None:
    anchor = ' %s  N    CB   single' % resname
    assert anchor in text
    text = text.replace(anchor, extra_bond + '\n' + anchor, 1)
  return [(resname + '.cif', iotbx.cif.reader(input_string=text).model())]

def place(pdb_string=pdb_str, restraint_objects=None, cls=None,
          generated=None, **kw):
  '''generated: restraint objects the caller made from the CCD (reduce2's
  generated_restraint_objects).'''
  inp = iotbx.pdb.input(lines=pdb_string.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp,
    restraint_objects=restraint_objects, log=null_out())
  if cls is None: cls = reduce_hydrogen.place_hydrogens
  o = cls(model=model, **kw)
  if generated is not None: o.generated_restraint_objects = generated
  o.run()
  return o

def ht7_atoms(o):
  r = {}
  for a in o.get_model().get_hierarchy().atoms():
    if a.parent().resname.strip() == 'HT7': r[a.name.strip()] = a
    elif a.name.strip() == 'C': r['C_prev'] = a
  return r

def angle(a, b, c):
  v = matrix.col(b.xyz)
  return (matrix.col(a.xyz)-v).angle(matrix.col(c.xyz)-v, deg=True)

def assert_planar_peptide_h(o, d_nh):
  a = ht7_atoms(o)
  assert 'HA2' not in a
  h, n, c, cb = a['HN1'], a['N'], a['C_prev'], a['CB']
  assert approx_equal(h.distance(n), d_nh, eps=0.005), h.distance(n)
  total = angle(h, n, c) + angle(h, n, cb) + angle(c, n, cb)
  assert approx_equal(total, 360, eps=1), total
  # the target riding H uses, not only the written length
  rc = o.get_model().riding_h_manager.h_parameterization[h.i_seq]
  assert approx_equal(rc.disth, d_nh, eps=0.005), rc.disth

def assert_link_restraints(o):
  '''One C-N-H, C-N-X, X-N-H and plane C, N, X, H each, at the link values
  (X-N-H = 360 - 124.3 - 121.7): riding puts the H on the bisector whatever
  the angles are, so the positions alone would not show them.'''
  a = ht7_atoms(o)
  h, n, c, cb = [a[k].i_seq for k in ('HN1', 'N', 'C_prev', 'CB')]
  grm = o.get_model().get_restraints_manager().geometry
  for i_seqs, ideal in (((c, n, h), 124.3), ((c, n, cb), 121.7),
                        ((cb, n, h), 114.0)):
    found = [p.angle_ideal for p in grm.angle_proxies
             if sorted(p.i_seqs) == sorted(i_seqs)]
    assert len(found) == 1, (i_seqs, found)
    assert approx_equal(found[0], ideal, eps=0.01), (i_seqs, found)
  planes = [p for p in grm.planarity_proxies
            if set(p.i_seqs) == set((c, n, cb, h))]
  assert len(planes) == 1, len(planes)

def rerun(o, **kw):
  return place(o.get_model().model_as_pdb(), keep_existing_H=True, **kw)

def test_beta_amino_acid():
  o = place()
  o2 = rerun(o)
  for x in (o, o2):
    assert_planar_peptide_h(x, 0.86)
    assert_link_restraints(x)

def test_neutron_length():
  o = place(use_neutron_distances=True)
  o2 = rerun(o, use_neutron_distances=True)
  for x in (o, o2):
    assert_planar_peptide_h(x, 1.02)
    assert_link_restraints(x)

class supplied_c_n_h(reduce_hydrogen.place_hydrogens):
  '''A C(i-1)-N-H angle that exists before add_link_h_restraints.'''
  def add_link_h_restraints(self):
    grm = self.model.get_restraints_manager().geometry
    name = {}
    for atom in self.model.get_hierarchy().atoms():
      name[atom.name.strip(), atom.parent().resname.strip()] = atom.i_seq
    grm.add_angles_in_place([geometry_restraints.angle_proxy(
      i_seqs      = [name['C', 'DNE'], name['N', 'HT7'], name['HN1', 'HT7']],
      angle_ideal = 124.3,
      weight      = 1./3.**2)])
    reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)

def test_supplied_angle():
  o = place(cls=supplied_c_n_h)
  assert_planar_peptide_h(o, 0.86)
  assert_link_restraints(o)

class supplied_plane(reduce_hydrogen.place_hydrogens):
  '''A plane C(i-1), N, CB, H that exists before add_link_h_restraints.'''
  def add_link_h_restraints(self):
    grm = self.model.get_restraints_manager().geometry
    name = {}
    for atom in self.model.get_hierarchy().atoms():
      name[atom.name.strip(), atom.parent().resname.strip()] = atom.i_seq
    i_seqs = [name['C', 'DNE'], name['N', 'HT7'], name['CB', 'HT7'],
              name['HN1', 'HT7']]
    grm.add_planarities_in_place([geometry_restraints.planarity_proxy(
      i_seqs  = flex.size_t(i_seqs),
      weights = flex.double([1./0.02**2]*4))])
    reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)

def test_supplied_plane():
  o = place(cls=supplied_plane)
  assert_planar_peptide_h(o, 0.86)
  assert_link_restraints(o)

def test_user_dictionary_keeps_n_h():
  o = place(restraint_objects=ht7_cif())
  assert_planar_peptide_h(o, 0.89)
  assert_planar_peptide_h(rerun(o, restraint_objects=ht7_cif()), 0.89)

def test_generated_dictionary():
  # a dictionary generated from the CCD (HT7's as a stand-in) is not the
  # user's: peptide N-H
  for neutron, d_nh in ((False, 0.86), (True, 1.02)):
    cif = ht7_cif()[0][1]
    ro = reduce_hydrogen.with_auto_restraint_object([], 'HT7', cif, [])
    o = place(restraint_objects=ro, generated=[cif],
              use_neutron_distances=neutron)
    assert_planar_peptide_h(o, d_nh)
    assert_link_restraints(o)

def test_user_file_named_like_generated():
  # the name does not decide: a user file called like a generated object
  # keeps its own N-H
  for name in ('auto_HT7.cif', 'auto_HT7'):
    cif = ht7_cif()[0][1]
    o = place(restraint_objects=[(name, cif)])
    assert_planar_peptide_h(o, 0.89)
    assert 'HT7' in o._user_defined_resnames()

def test_reduce2_generated_objects():
  # reduce2 hands its CCD objects (BGM, CCD-only) to placement; ZQ4 (code 04)
  # comes from a user CIF
  from mmtbx.regression import reduce2_neutron_utils as u
  u.set_prefix('tst_add_hydrogen_23')
  zq4 = u.heavy_pdb('04', u.zq4_atoms, offset=(10., 10., 2.), cryst1=False)
  pdb = (u.heavy_pdb('BGM', u.bgm_atoms(), offset=(10., 10., 15.)) +
         zq4.replace(' A   1', ' B   1'))
  seen = []
  run = reduce_hydrogen.place_hydrogens.run
  def record(self):
    seen.append(self)
    return run(self)
  reduce_hydrogen.place_hydrogens.run = record
  try:
    u.run_prog('tst_add_hydrogen_23_bgm', pdb, u.zq4_cif)
  finally:
    reduce_hydrogen.place_hydrogens.run = run
  assert len(seen) == 1, len(seen)
  names = seen[0]._user_defined_resnames()
  assert names == set(['04']), names
  generated = [o for n, o in seen[0].get_model().get_restraint_objects()
               if n == 'auto_BGM']
  assert len(generated) == 1
  assert [g for g in seen[0].generated_restraint_objects
          if g is generated[0]]

def test_two_candidates():
  # N with CB and CG in the residue: which one is X is not known. Class other,
  # so that h_on_tertiary_amide_n does not remove the H first.
  o = place(pdb_str.replace('HT7', 'X7X'), restraint_objects=ht7_cif(
    ' X7X  N    CG   single    1.456  0.020  1.456', 'X7X'))
  names = [a.name.strip() for a in o.get_model().get_hierarchy().atoms()
           if a.parent().resname.strip() == 'X7X']
  assert 'HN1' not in names and 'HA2' not in names, names
  assert [l for l in o.site_labels_no_para if 'HN1' in l], o.site_labels_no_para

if (__name__ == "__main__"):
  test_beta_amino_acid()
  test_neutron_length()
  test_supplied_angle()
  test_supplied_plane()
  test_user_dictionary_keeps_n_h()
  test_generated_dictionary()
  test_user_file_named_like_generated()
  test_reduce2_generated_objects()
  test_two_candidates()
  print("OK")
