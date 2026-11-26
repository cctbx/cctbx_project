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

# GeoStd HT7 as it was until 2026-11-25 (GeoStd 55c1ed3ca3): the amide H is
# HN1, with HA2 on N, HA on CB and HB1, HC2 on CA. GeoStd now has H, HA2 on
# N, HB on CB and HA, HC2 on CA; the TRANS link finds an H named H itself
cif_ht7_old = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 HT7  HT7  '(3S)-3-amino-4-(1H-indol-3-yl)butanoic acid'  L-peptide  29  16  .  2022-03-18  2023-10-30
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as GRAND
Added dihedrals for adding hydrogens : 2023-10-30
;

data_comp_HT7
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 HT7  N    N  NH2    0  -0.790   1.3266   0.9214   0.5800
 HT7  CB   C  CH1    0  -0.077   1.6663  -0.3222  -0.0963
 HT7  CG   C  CH2    0  -0.528   0.5913  -0.8239  -1.0805
 HT7  CD   C  CR5    0   0.013  -0.7566  -0.9900  -0.4580
 HT7  CE2  C  CR56   0   0.013  -1.7868   0.0015  -0.3406
 HT7  CZ3  C  CR16   0  -0.374  -1.9069   1.3293  -0.7697
 HT7  CH3  C  CR16   0  -0.269  -3.0675   2.0233  -0.4901
 HT7  CT2  C  CR16   0  -0.286  -4.1205   1.4184   0.2177
 HT7  CH2  C  CR16   0  -0.306  -4.0338   0.1099   0.6520
 HT7  CZ2  C  CR56   0   0.351  -2.8607  -0.5879   0.3648
 HT7  NZ1  N  NR15   0  -0.686  -2.4965  -1.8677   0.6586
 HT7  CE1  C  CR15   0  -0.084  -1.2423  -2.1054   0.1668
 HT7  CA   C  CH2    0  -0.556   3.0038  -0.1565  -0.8050
 HT7  C    C  C      0   0.409   4.1609   0.1765   0.1502
 HT7  O    O  O      0  -0.624   4.3071  -0.5823   1.1315
 HT7  OXT  O  OC    -1  -0.623   4.8798   1.1561  -0.1349
 HT7  HN1  H  HNH2   0   0.328   1.2310   1.6635  -0.1053
 HT7  HA2  H  HNH2   0   0.345   0.4183   0.8342   1.0209
 HT7  HA   H  HCH1   0   0.208   1.7843  -1.0882   0.6761
 HT7  HA1  H  HCH2   0   0.221   0.5293  -0.1233  -1.9216
 HT7  HB2  H  HCH2   0   0.231   0.9235  -1.7775  -1.5012
 HT7  HB1  H  HCH2   0   0.221   3.2403  -1.0909  -1.3269
 HT7  HC2  H  HCH2   0   0.196   2.9076   0.6130  -1.5783
 HT7  HE1  H  HCR5   0   0.277  -0.7800  -3.0719   0.2969
 HT7  HZ3  H  HCR6   0   0.237  -1.1023   1.8080  -1.3148
 HT7  HH3  H  HCR6   0   0.244  -3.1725   3.0498  -0.8173
 HT7  HT2  H  HCR6   0   0.247  -5.0163   1.9891   0.4266
 HT7  HH2  H  HCR6   0   0.257  -4.8465  -0.3551   1.1958
 HT7  HZ1  H  HNR5   0   0.406  -3.0576  -2.5392   1.1545

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 HT7  N    CB   single    1.456  0.020  1.456
 HT7  CB   CG   single    1.541  0.020  1.541
 HT7  CB   CA   single    1.530  0.020  1.530
 HT7  CG   CD   single    1.494  0.020  1.494
 HT7  CD   CE2  aromatic  1.435  0.020  1.435
 HT7  CD   CE1  aromatic  1.368  0.020  1.368
 HT7  CE2  CZ3  aromatic  1.401  0.020  1.401
 HT7  CE2  CZ2  aromatic  1.414  0.020  1.414
 HT7  CZ3  CH3  aromatic  1.381  0.020  1.381
 HT7  CH3  CT2  aromatic  1.406  0.020  1.406
 HT7  CT2  CH2  aromatic  1.382  0.020  1.382
 HT7  CH2  CZ2  aromatic  1.395  0.020  1.395
 HT7  CZ2  NZ1  aromatic  1.363  0.020  1.363
 HT7  NZ1  CE1  aromatic  1.368  0.020  1.368
 HT7  CA   C    single    1.525  0.021  1.525
 HT7  C    O    deloc     1.231  0.020  1.231
 HT7  C    OXT  deloc     1.248  0.020  1.248
 HT7  N    HN1  single    0.890  0.020  1.040
 HT7  N    HA2  single    0.890  0.020  1.040
 HT7  CB   HA   single    0.970  0.020  1.090
 HT7  CG   HA1  single    0.970  0.020  1.090
 HT7  CG   HB2  single    0.970  0.020  1.090
 HT7  CA   HB1  single    0.970  0.020  1.090
 HT7  CA   HC2  single    0.970  0.020  1.090
 HT7  CE1  HE1  single    0.930  0.020  1.080
 HT7  CZ3  HZ3  single    0.930  0.020  1.080
 HT7  CH3  HH3  single    0.930  0.020  1.080
 HT7  CT2  HT2  single    0.930  0.020  1.080
 HT7  CH2  HH2  single    0.930  0.020  1.080
 HT7  NZ1  HZ1  single    0.860  0.020  1.020

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 HT7  HA2  N    HN1  105.81  3.000
 HT7  HA2  N    CB   109.74  3.000
 HT7  HN1  N    CB   109.45  3.000
 HT7  HA   CB   CA   109.00  3.000
 HT7  HA   CB   CG   107.34  3.000
 HT7  CA   CB   CG   110.54  3.000
 HT7  HA   CB   N    107.17  3.000
 HT7  CA   CB   N    110.50  1.700
 HT7  CG   CB   N    114.32  3.000
 HT7  HB2  CG   HA1  106.21  3.000
 HT7  HB2  CG   CD   109.71  3.000
 HT7  HA1  CG   CD   109.86  3.000
 HT7  HB2  CG   CB   108.51  3.000
 HT7  HA1  CG   CB   108.74  3.000
 HT7  CD   CG   CB   113.54  3.000
 HT7  CE1  CD   CE2  105.74  3.000
 HT7  CE1  CD   CG   126.98  3.000
 HT7  CE2  CD   CG   127.24  3.000
 HT7  CZ2  CE2  CZ3  118.88  3.000
 HT7  CZ2  CE2  CD   107.35  3.000
 HT7  CZ3  CE2  CD   133.77  3.000
 HT7  HZ3  CZ3  CH3  120.27  3.000
 HT7  HZ3  CZ3  CE2  120.62  3.000
 HT7  CH3  CZ3  CE2  119.11  3.000
 HT7  HH3  CH3  CT2  119.19  3.000
 HT7  HH3  CH3  CZ3  119.79  3.000
 HT7  CT2  CH3  CZ3  121.02  3.000
 HT7  HT2  CT2  CH2  119.38  3.000
 HT7  HT2  CT2  CH3  119.36  3.000
 HT7  CH2  CT2  CH3  121.26  3.000
 HT7  HH2  CH2  CZ2  121.32  3.000
 HT7  HH2  CH2  CT2  121.17  3.000
 HT7  CZ2  CH2  CT2  117.51  3.000
 HT7  NZ1  CZ2  CH2  130.56  3.000
 HT7  NZ1  CZ2  CE2  107.22  3.000
 HT7  CH2  CZ2  CE2  122.21  3.000
 HT7  HZ1  NZ1  CE1  124.94  3.000
 HT7  HZ1  NZ1  CZ2  125.75  3.000
 HT7  CE1  NZ1  CZ2  109.31  3.000
 HT7  HE1  CE1  NZ1  120.33  3.000
 HT7  HE1  CE1  CD   129.30  3.000
 HT7  NZ1  CE1  CD   110.37  3.000
 HT7  HC2  CA   HB1  106.37  3.000
 HT7  HC2  CA   C    110.66  3.000
 HT7  HB1  CA   C    108.55  3.000
 HT7  HC2  CA   CB   109.15  3.000
 HT7  HB1  CA   CB   108.56  3.000
 HT7  C    CA   CB   110.10  1.900
 HT7  OXT  C    O    126.07  3.000
 HT7  OXT  C    CA   117.50  3.000
 HT7  O    C    CA   120.80  1.700

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 HT7  CONST_01  CH3  CZ3  CE2  CD    180.00   0.0  0
 HT7  CONST_02  CH2  CZ2  CE2  CD    180.00   0.0  0
 HT7  CONST_03  NZ1  CZ2  CE2  CD      0.00   0.0  0
 HT7  CONST_04  CZ2  NZ1  CE1  CD      0.00   0.0  0
 HT7  CONST_05  NZ1  CE1  CD   CE2     0.00   0.0  0
 HT7  CONST_06  CT2  CH3  CZ3  CE2     0.00   0.0  0
 HT7  CONST_07  CT2  CH2  CZ2  CE2     0.00   0.0  0
 HT7  CONST_08  CE1  NZ1  CZ2  CE2     0.00   0.0  0
 HT7  CONST_09  CE1  CD   CE2  CZ3   180.00   0.0  0
 HT7  CONST_10  CH2  CZ2  CE2  CZ3     0.00   0.0  0
 HT7  CONST_11  NZ1  CZ2  CE2  CZ3   180.00   0.0  0
 HT7  CONST_12  CH2  CT2  CH3  CZ3     0.00   0.0  0
 HT7  CONST_13  CZ2  CE2  CZ3  CH3     0.00   0.0  0
 HT7  CONST_14  CZ2  CH2  CT2  CH3     0.00   0.0  0
 HT7  CONST_15  NZ1  CZ2  CH2  CT2   180.00   0.0  0
 HT7  CONST_16  CE1  NZ1  CZ2  CH2   180.00   0.0  0
 HT7  CONST_17  CE1  CD   CE2  CZ2     0.00   0.0  0
 HT7  CONST_18  CZ3  CE2  CD   CG      0.00   0.0  0
 HT7  CONST_19  CZ2  CE2  CD   CG    180.00   0.0  0
 HT7  CONST_20  NZ1  CE1  CD   CG    180.00   0.0  0
 HT7  CONST_21  HZ3  CZ3  CE2  CD      0.00   0.0  0
 HT7  CONST_22  HZ1  NZ1  CE1  CD    180.00   0.0  0
 HT7  CONST_23  HE1  CE1  CD   CE2   180.00   0.0  0
 HT7  CONST_24  HH3  CH3  CZ3  CE2   180.00   0.0  0
 HT7  CONST_25  HH2  CH2  CZ2  CE2   180.00   0.0  0
 HT7  CONST_26  HT2  CT2  CH3  CZ3   180.00   0.0  0
 HT7  Var_01    CE2  CD   CG   CB     88.59  30.0  2
 HT7  Var_02    O    C    CA   CB     52.13  30.0  3
 HT7  Var_03    HA   CB   CG   CD     63.26  30.0  1
 HT7  Var_04    HB1  CA   CB   N    -178.89  30.0  3
 HT7  Var_05    HC2  CA   CB   N     -63.34  30.0  3
 HT7  Var_06    HN1  N    CB   CG    -62.55  30.0  3
 HT7  Var_07    HA2  N    CB   CG     53.15  30.0  3

loop_
_chem_comp_chir.comp_id
_chem_comp_chir.id
_chem_comp_chir.atom_id_centre
_chem_comp_chir.atom_id_1
_chem_comp_chir.atom_id_2
_chem_comp_chir.atom_id_3
_chem_comp_chir.volume_sign
 HT7  chir_01  CB  N  CG  CA  negativ

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 HT7  plan-1  CG   0.020
 HT7  plan-1  CD   0.020
 HT7  plan-1  CE2  0.020
 HT7  plan-1  CZ3  0.020
 HT7  plan-1  CH3  0.020
 HT7  plan-1  CT2  0.020
 HT7  plan-1  CH2  0.020
 HT7  plan-1  CZ2  0.020
 HT7  plan-1  NZ1  0.020
 HT7  plan-1  CE1  0.020
 HT7  plan-1  HE1  0.020
 HT7  plan-1  HZ3  0.020
 HT7  plan-1  HH3  0.020
 HT7  plan-1  HT2  0.020
 HT7  plan-1  HH2  0.020
 HT7  plan-1  HZ1  0.020
 HT7  plan-2  CA   0.020
 HT7  plan-2  C    0.020
 HT7  plan-2  O    0.020
 HT7  plan-2  OXT  0.020
"""

def old_ht7_library(f):
  '''f with cif_ht7_old in the place of the library's HT7: the monomer
  library server gets it where it reads geostd/h/data_HT7.cif, and
  reduce_hydrogen._dictionary_sites() (which reads that file again for the
  ideal sites) a copy of it, so HT7 is a library entry as before (no
  restraint file, not user-defined).'''
  def wrapper(*args, **kw):
    from mmtbx.monomer_library import server
    read_cif = server.read_cif
    sites = reduce_hydrogen._dictionary_sites
    def read_old(file_name):
      if os.path.basename(file_name) == 'data_HT7.cif':
        return iotbx.cif.reader(input_string=cif_ht7_old, strict=False).model()
      return read_cif(file_name=file_name)
    def sites_old(source_info, resname):
      if (source_info and source_info.startswith('file:') and
          os.path.basename(source_info[5:].strip()) == 'data_HT7.cif'):
        fn = os.path.abspath('tst_add_hydrogen_23_old_data_HT7.cif')
        with open(fn, 'w') as out:
          out.write(cif_ht7_old)
        source_info = 'file: ' + fn
      return sites(source_info, resname)
    server.read_cif = read_old
    reduce_hydrogen._dictionary_sites = sites_old
    try:
      return f(*args, **kw)
    finally:
      server.read_cif = read_cif
      reduce_hydrogen._dictionary_sites = sites
  wrapper.__name__ = f.__name__
  return wrapper

def ht7_cif(extra_bond=None, resname='HT7'):
  '''The old geostd HT7 dictionary (cif_ht7_old) as a user CIF, called
  resname; extra_bond: a line to add.'''
  text = cif_ht7_old.replace('HT7', resname)
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

@old_ht7_library
def test_beta_amino_acid():
  o = place()
  o2 = rerun(o)
  for x in (o, o2):
    assert_planar_peptide_h(x, 0.86)
    assert_link_restraints(x)

@old_ht7_library
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

@old_ht7_library
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

@old_ht7_library
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
