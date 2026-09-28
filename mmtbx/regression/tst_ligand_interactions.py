from __future__ import absolute_import, division, print_function
'''
Tests for mmtbx.validation.ligand_interactions.

Model: ligand EDO A 1 and other EDO copies placed around it (ideal geostd
geometry, H at electron-cloud X-H lengths after processing):
  EDO A 2   O1 accepts an H-bond from the ligand's O1-HO1 (H...A 2.01 A, 180 deg)
  EDO A 3   H12 clashes with the ligand's H11 (1.75 A, overlap 0.69 A)
  EDO A 4   C1 and H12 in van der Waals contact with the ligand's O2 and HO2
  EDO A 0   15 A away, first in the file: outside validate_ligands' 3 A region,
            so its records' i_seqs must be mapped to the full model
The other copies are environment: probe2 is run by the ligand's selection, not
its resname.
'''
import json
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from scitbx.array_family import flex
from cctbx import sgtbx
from mmtbx.validation import ligand_interactions as LI

LIG_SEL = 'chain A and resseq 1 and resname EDO'

model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   0      30.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   0      31.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   0      28.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   0      28.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   0      30.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   0      30.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   0      31.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   0      28.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   0      28.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   0      27.463  10.001  13.644  1.00 20.00           H
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
HETATM   11  C1  EDO A   2      19.869   9.698  14.350  1.00 20.00           C
HETATM   12  O1  EDO A   2      18.690   9.535  15.106  1.00 20.00           O
HETATM   13  C2  EDO A   2      19.563  10.454  13.082  1.00 20.00           C
HETATM   14  O2  EDO A   2      19.185  11.776  13.394  1.00 20.00           O
HETATM   15  H11 EDO A   2      20.643  10.236  14.914  1.00 20.00           H
HETATM   16  H12 EDO A   2      20.303   8.733  14.059  1.00 20.00           H
HETATM   17  HO1 EDO A   2      18.917   9.078  15.920  1.00 20.00           H
HETATM   18  H21 EDO A   2      20.466  10.425  12.458  1.00 20.00           H
HETATM   19  H22 EDO A   2      18.778   9.926  12.526  1.00 20.00           H
HETATM   20  HO2 EDO A   2      18.957  12.223  12.574  1.00 20.00           H
HETATM   21  C1  EDO A   3      15.673   6.534  16.101  1.00 20.00           C
HETATM   22  O1  EDO A   3      16.471   6.047  15.046  1.00 20.00           O
HETATM   23  C2  EDO A   3      14.542   5.574  16.373  1.00 20.00           C
HETATM   24  O2  EDO A   3      15.055   4.361  16.877  1.00 20.00           O
HETATM   25  H11 EDO A   3      16.256   6.673  17.022  1.00 20.00           H
HETATM   26  H12 EDO A   3      15.226   7.507  15.860  1.00 20.00           H
HETATM   27  HO1 EDO A   3      17.199   6.660  14.913  1.00 20.00           H
HETATM   28  H21 EDO A   3      13.865   6.060  17.088  1.00 20.00           H
HETATM   29  H22 EDO A   3      13.968   5.421  15.449  1.00 20.00           H
HETATM   30  HO2 EDO A   3      14.320   3.755  17.002  1.00 20.00           H
HETATM   31  C1  EDO A   4      10.915   8.639  11.757  1.00 20.00           C
HETATM   32  O1  EDO A   4      11.786   9.042  10.725  1.00 20.00           O
HETATM   33  C2  EDO A   4       9.559   8.312  11.184  1.00 20.00           C
HETATM   34  O2  EDO A   4       9.650   7.170  10.362  1.00 20.00           O
HETATM   35  H11 EDO A   4      11.299   7.763  12.296  1.00 20.00           H
HETATM   36  H12 EDO A   4      10.772   9.431  12.504  1.00 20.00           H
HETATM   37  HO1 EDO A   4      12.649   9.215  11.112  1.00 20.00           H
HETATM   38  H21 EDO A   4       8.876   8.153  12.028  1.00 20.00           H
HETATM   39  H22 EDO A   4       9.181   9.182  10.631  1.00 20.00           H
HETATM   40  HO2 EDO A   4       8.785   7.008   9.978  1.00 20.00           H
END
'''


# The ligand alone in a 6 A cell along x: its x+1 and x-1 copies make O-H...O
# H-bonds (pnp, symmetry); EDO A 2 is 15 A away along y as probe2's target.
sym_model_str = '''
CRYST1    6.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
HETATM   11  C1  EDO A   2      15.177  25.059  15.375  1.00 20.00           C
HETATM   12  O1  EDO A   2      16.075  25.413  14.347  1.00 20.00           O
HETATM   13  C2  EDO A   2      13.776  25.464  14.992  1.00 20.00           C
HETATM   14  O2  EDO A   2      13.340  24.696  13.892  1.00 20.00           O
HETATM   15  H11 EDO A   2      15.198  23.980  15.580  1.00 20.00           H
HETATM   16  H12 EDO A   2      15.418  25.564  16.319  1.00 20.00           H
HETATM   17  HO1 EDO A   2      16.953  25.118  14.602  1.00 20.00           H
HETATM   18  H21 EDO A   2      13.139  25.312  15.873  1.00 20.00           H
HETATM   19  H22 EDO A   2      13.758  26.539  14.771  1.00 20.00           H
HETATM   20  HO2 EDO A   2      12.463  25.001  13.644  1.00 20.00           H
END
'''

# EDO A 1 folded (O1-C1-C2-O2 52 deg, both hydroxyl H turned inwards): HO1...HO2
# 1.61 A, an intramolecular clash (pnp overlap 0.49 A); EDO A 2, a copy 15 A away
# along y, is probe2's target.
fold_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.467   9.962  13.710  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      15.708  10.111  13.512  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.099  10.066  15.759  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.692  11.558  15.033  1.00 20.00           H
HETATM   10  HO2 EDO A   1      14.288   9.874  13.218  1.00 20.00           H
HETATM   11  C1  EDO A   2      15.177  25.059  15.375  1.00 20.00           C
HETATM   12  O1  EDO A   2      16.075  25.413  14.347  1.00 20.00           O
HETATM   13  C2  EDO A   2      13.776  25.464  14.992  1.00 20.00           C
HETATM   14  O2  EDO A   2      13.467  24.962  13.710  1.00 20.00           O
HETATM   15  H11 EDO A   2      15.198  23.980  15.580  1.00 20.00           H
HETATM   16  H12 EDO A   2      15.418  25.564  16.319  1.00 20.00           H
HETATM   17  HO1 EDO A   2      15.708  25.111  13.512  1.00 20.00           H
HETATM   18  H21 EDO A   2      13.099  25.066  15.759  1.00 20.00           H
HETATM   19  H22 EDO A   2      13.692  26.558  15.033  1.00 20.00           H
HETATM   20  HO2 EDO A   2      14.288  24.874  13.218  1.00 20.00           H
END
'''

# A probe2-only H-bond (D-H...A 113 deg, below pnp's 120): EDO A 1 O1-HO1 donates
# to EDO A 2 O1. EDO A 1 has conformers A and B, B shifted 0.45 A along O1->HO1,
# so conformer B's O1 is the heavy atom nearest to conformer A's HO1.
alt_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1 AEDO A   1      15.177  10.059  15.375  0.50 20.00           C
HETATM    2  O1 AEDO A   1      16.075  10.413  14.347  0.50 20.00           O
HETATM    3  C2 AEDO A   1      13.776  10.464  14.992  0.50 20.00           C
HETATM    4  O2 AEDO A   1      13.340   9.696  13.892  0.50 20.00           O
HETATM    5  H11AEDO A   1      15.198   8.980  15.580  0.50 20.00           H
HETATM    6  H12AEDO A   1      15.418  10.564  16.319  0.50 20.00           H
HETATM    7  HO1AEDO A   1      16.953  10.118  14.602  0.50 20.00           H
HETATM    8  H21AEDO A   1      13.139  10.312  15.873  0.50 20.00           H
HETATM    9  H22AEDO A   1      13.758  11.539  14.771  0.50 20.00           H
HETATM   10  HO2AEDO A   1      12.463  10.001  13.644  0.50 20.00           H
HETATM   11  C1 BEDO A   1      15.588   9.921  15.494  0.50 20.00           C
HETATM   12  O1 BEDO A   1      16.486  10.275  14.466  0.50 20.00           O
HETATM   13  C2 BEDO A   1      14.187  10.326  15.111  0.50 20.00           C
HETATM   14  O2 BEDO A   1      13.751   9.558  14.011  0.50 20.00           O
HETATM   15  H11BEDO A   1      15.609   8.842  15.699  0.50 20.00           H
HETATM   16  H12BEDO A   1      15.829  10.426  16.438  0.50 20.00           H
HETATM   17  HO1BEDO A   1      17.364   9.980  14.721  0.50 20.00           H
HETATM   18  H21BEDO A   1      13.550  10.174  15.992  0.50 20.00           H
HETATM   19  H22BEDO A   1      14.169  11.401  14.890  0.50 20.00           H
HETATM   20  HO2BEDO A   1      12.874   9.863  13.763  0.50 20.00           H
HETATM   21  C1  EDO A   2      19.271  12.379  15.630  1.00 20.00           C
HETATM   22  O1  EDO A   2      18.382  11.848  14.673  1.00 20.00           O
HETATM   23  C2  EDO A   2      20.053  11.264  16.277  1.00 20.00           C
HETATM   24  O2  EDO A   2      20.908  10.666  15.327  1.00 20.00           O
HETATM   25  H11 EDO A   2      19.970  13.098  15.181  1.00 20.00           H
HETATM   26  H12 EDO A   2      18.740  12.913  16.428  1.00 20.00           H
HETATM   27  HO1 EDO A   2      17.920  12.579  14.256  1.00 20.00           H
HETATM   28  H21 EDO A   2      20.615  11.700  17.113  1.00 20.00           H
HETATM   29  H22 EDO A   2      19.355  10.537  16.711  1.00 20.00           H
HETATM   30  HO2 EDO A   2      21.363   9.935  15.753  1.00 20.00           H
END
'''

# Names probe2's raw atom field runs together: a 4-character atom name with an
# altloc, a 4-digit residue number with a one-letter chain, an insertion code.
names_str = '''
ATOM      1 HD21AASN A1001B      1.000   2.000   3.000  1.00 20.00           H
ATOM      2 HD21BASN A1001B      1.300   2.000   3.000  1.00 20.00           H
ATOM      3  ND2 ASN A1001B      1.500   2.500   3.000  1.00 20.00           N
ATOM      4  CA  GLY B  52A      5.000   2.000   3.000  1.00 20.00           C
END
'''

# ------------------------------------------------------------------------------

def get_model(lines=None):
  if lines is None:
    lines = model_str.split('\n')
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(lines=lines, source_info=None), log=null_out())
  model.process(make_restraints=True)
  model.set_hydrogen_bond_length(use_neutron_distances=False)
  return model

def get_manager(model, sel=LIG_SEL, **kwargs):
  params = LI.master_params().extract().ligand_interactions
  for k, v in kwargs.items():
    setattr(params, k, v)
  isel = model.selection(sel).iselection()
  return LI.manager(model, isel, sel, params=params).run()

def short(label):
  '''"A EDO 3 H12" -> "EDO 3 H12"'''
  return None if label is None else label.split(None, 1)[1]

def entries_by_type(m):
  result = {}
  for e in m.entries:
    result.setdefault(e['type'], []).append(e)
  return result

# ------------------------------------------------------------------------------

def exercise_entries_and_counts(model):
  '''
  One H-bond and one clash (pnp and probe2 agree), vdW contacts with probe's
  class as subtype; counts; criteria recorded; as_dict is JSON-able; show() runs.
  '''
  m = get_manager(model)
  by = entries_by_type(m)
  assert len(by['hbond']) == 1, by['hbond']
  hb = by['hbond'][0]
  assert [short(l) for l in hb['labels']] == ['EDO 1 O1', 'EDO 1 HO1', 'EDO 2 O1']
  assert hb['sources'] == ['pnp', 'probe2']
  assert hb['cross_check'] == 'pnp and probe2'
  assert hb['residue'] == 'A EDO 2' and hb['symop'] is None
  assert approx_equal(hb['geometry']['pnp']['d_HA'], 2.011, eps=0.01)
  assert hb['geometry']['probe2']['dots']['hb'] > 0
  assert len(by['clash']) == 1, by['clash']
  cl = by['clash'][0]
  assert [short(l) for l in cl['labels']] == ['EDO 1 H11', 'EDO 3 H12']
  assert cl['sources'] == ['pnp', 'probe2']
  assert approx_equal(cl['geometry']['pnp']['overlap'], 0.69, eps=0.01)
  assert approx_equal(cl['geometry']['probe2']['min_gap'], -0.69, eps=0.01)
  vdw = dict([(tuple([short(l) for l in e['labels']]), e) for e in by['vdw']])
  assert vdw[('EDO 1 HO2', 'EDO 4 H12')]['subtype'] == 'so'
  assert vdw[('EDO 1 HO2', 'EDO 4 C1')]['subtype'] == 'cc'
  assert vdw[('EDO 1 O2', 'EDO 4 C1')]['subtype'] == 'wc'
  for e in by['vdw']:
    assert e['sources'] == ['probe2'] and e['cross_check'] is None
    assert e['geometry']['probe2']['min_gap'] is not None
  for e in m.entries:
    assert e['model_support'] is None
  assert m.disagreements == [] and m.internal == []
  assert m.probe_unmapped == set()
  c = m.counts()
  assert c.per_type['hbond'] == 1 and c.per_type['clash'] == 1
  assert c.per_type['vdw:so'] == 1 and c.per_type['vdw:cc'] == 1
  assert sum(c.per_type.values()) == len(m.entries)
  assert c.per_residue['A EDO 3'] == {'clash': 1}
  assert c.per_ligand_atom['A EDO 1 H11'] == {'clash': 1}
  d = json.loads(json.dumps(m.as_dict()))
  assert d['pair_class_order'] == ['bo', 'hb', 'so', 'cc', 'wc']
  assert d['probe']['density'] == 16
  assert d['hbond_criteria']['d_HA_cutoff'] == [1.4, 2.8]
  assert d['hbond_criteria']['d_DA_cutoff'] == [2.4, 4.1]
  assert d['hbond_criteria']['a_DHA_cutoff'] == 120
  assert d['hbond_criteria']['min_bonds_H_A'] == 5
  assert d['clash_criteria']['min_overlap'] == 0.4
  assert len(d['entries']) == len(m.entries)
  log = StringIO()
  m.show(log=log)
  assert 'contact patches' in log.getvalue()
  assert 'd_HA_cutoff=[1.4, 2.8]' in log.getvalue()
  return m

def exercise_pair_class_order(model):
  '''
  Plumbing: pair_class_order changes the label (the clash pair has so and bo
  dots; with so first it becomes a vdW contact and pnp's clash a disagreement)
  and invalid orders are refused. Not a recommended setting.
  '''
  m = get_manager(model, pair_class_order=['so', 'bo', 'hb', 'cc', 'wc'])
  by = entries_by_type(m)
  cl = by['clash'][0]
  assert cl['sources'] == ['pnp'] and cl['cross_check'] == 'pnp only', cl
  pair = [e for e in by['vdw'] if [short(l) for l in e['labels']] ==
    ['EDO 1 H11', 'EDO 3 H12']]
  assert len(pair) == 1 and pair[0]['subtype'] == 'so'
  assert pair[0]['geometry']['probe2']['dots']['bo'] > 0
  assert len(m.disagreements) == 1
  d = m.disagreements[0]
  assert d['type'] == 'clash' and d['missing'] == ['probe2']
  assert d['probe_class'] == 'so'
  # each class once
  for bad in (['bo', 'hb', 'so', 'cc'], ['bo', 'hb', 'so', 'cc', 'wc', 'wc'],
              ['bo', 'hb', 'so', 'cc', 'xx']):
    try:
      get_manager(model, pair_class_order=bad)
    except Sorry:
      pass
    else:
      raise AssertionError('pair_class_order %s accepted' % bad)

def exercise_internal():
  '''
  An intramolecular clash (folded EDO, HO1...HO2): listed as internal, not an
  entry, not a disagreement.
  '''
  m = get_manager(get_model(fold_model_str.split('\n')), sel='chain A and resseq 1')
  assert len(m.overlaps.clash_records) == 1
  assert [(d['type'], [short(l) for l in d['labels']]) for d in m.internal] == \
    [('clash', ['EDO 1 HO1', 'EDO 1 HO2'])], m.internal
  assert approx_equal(m.internal[0]['geometry']['pnp']['overlap'], 0.49, eps=0.01)
  assert m.entries == [] and m.disagreements == []
  assert len(json.loads(json.dumps(m.as_dict()))['internal']) == 1

def exercise_symmetry():
  '''
  H-bonds between the ligand and its symmetry-related copies: kept with the
  operator, "symmetry, probe2 not applicable", not disagreements, not internal.
  '''
  model = get_model(sym_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 1')
  by = entries_by_type(m)
  assert len(by['hbond']) == 2, m.entries
  # pnp's operator belongs to its pair: x+1,y,z for both records
  for r in m.overlaps.hbond_records:
    assert r['symop'] == 'x+1,y,z'
  uc = model.crystal_symmetry().unit_cell()
  atoms = model.get_hierarchy().atoms()
  def site(i, op):
    return uc.orthogonalize(sgtbx.rt_mx(op) * uc.fractionalize(atoms[i].xyz))
  ops = []
  for e in by['hbond']:
    assert e['cross_check'] == 'symmetry, probe2 not applicable'
    assert e['sources'] == ['pnp']
    d, h, a = e['atoms']
    # the ligand donor stays, the acceptor (its own copy) moves
    assert e['operators'][:2] == ['x,y,z', 'x,y,z'] and e['operators'][2] == e['symop']
    d_HA = e['geometry']['pnp']['d_HA']
    assert approx_equal(d_HA, 2.62, eps=0.01)
    assert approx_equal(atoms[h].distance(site(a, e['operators'][2])), d_HA, eps=1.e-6)
    # pnp's own operator on the acceptor gives the other copy
    assert approx_equal(atoms[h].distance(site(a, 'x+1,y,z')) if e['symop'] != 'x+1,y,z'
      else atoms[h].distance(site(a, 'x-1,y,z')), 9.55, eps=0.01)
    assert e['residue'] == 'A EDO 1 (%s)' % e['symop']
    assert e['labels'][2].endswith(' (%s)' % e['symop'])
    assert e['ligand_atoms'] == [e['labels'][0], e['labels'][1]]
    ops.append(e['symop'])
  assert sorted(ops) == ['x+1,y,z', 'x-1,y,z'], ops
  assert m.disagreements == [] and m.internal == []
  per_residue = m.counts().per_residue
  assert sorted(per_residue) == ['A EDO 1 (x+1,y,z)', 'A EDO 1 (x-1,y,z)'], per_residue
  assert '-' not in per_residue and None not in per_residue

def exercise_donor_conformers():
  '''
  A probe2-only H-bond from a donor with alternate conformers. probe2's hb class
  is overlap-based and ignores the angle, which is why only probe2 reports this
  H-bond (D-H...A 113 deg, below pnp's 120; recorded from the model). The donor
  is the H's bonded O1 of the same conformer (restraints' connectivity), not the
  nearest heavy atom (conformer B's O1); the partner residue is right with the
  donor, the acceptor or one donor conformer selected as the ligand.
  '''
  model = get_model(alt_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  by_name = dict([((a.parent().parent().resseq.strip(), a.name.strip(),
    a.parent().altloc), a) for a in atoms])
  h_a = by_name[('1', 'HO1', 'A')]
  o_a, o_b = by_name[('1', 'O1', 'A')], by_name[('1', 'O1', 'B')]
  assert h_a.distance(o_b) < h_a.distance(o_a) < 1.3  # the old lookup's choice
  for sel, residue in (('chain A and resseq 2', 'A EDO 1'),
                       ('chain A and resseq 1', 'A EDO 2'),
                       ('chain A and resseq 1 and altloc A', 'A EDO 2')):
    m = get_manager(model, sel=sel)
    assert m.overlaps.hbond_records == [] and m.unresolved_donors == []
    hb = [e for e in m.entries if e['type'] == 'hbond' and e['atoms'][1] == h_a.i_seq]
    assert len(hb) == 1, (sel, m.entries)
    e = hb[0]
    assert e['atoms'][0] == o_a.i_seq, (sel, e['labels'])
    assert e['labels'] == ['A EDO 1 O1 alt A', 'A EDO 1 HO1 alt A', 'A EDO 2 O1'], e
    assert e['residue'] == residue, (sel, e['residue'])
    assert e['cross_check'] == 'probe2 only'
    assert e['operators'] == ['x,y,z'] * 3 and e['symop'] is None
    g = e['geometry']['probe2']
    assert approx_equal(g['a_DHA'], h_a.angle(by_name[('2', 'O1', '')], o_a, deg=True),
      eps=1.e-6)
    assert 110 < g['a_DHA'] < 120, g['a_DHA']
    assert approx_equal(g['d_HA'], 2.28, eps=0.02), g['d_HA']
    n_b = len([x for x in m.entries if x['type'] == 'hbond'])
    assert n_b == (1 if 'altloc A' in sel else 2), n_b
  # an H without a bonded heavy atom: reported, donor unset
  m._fsc0 = [[] for a in atoms]
  assert m._donor_of(h_a.i_seq) is None
  assert m.unresolved_donors == [dict(h='A EDO 1 HO1 alt A', bonded_heavy=[])]

def exercise_probe_names():
  '''
  probe2's atom field (its own format) for names that run together; the key
  parsed by columns equals atom_key, the table maps each field to its atom.
  '''
  h = iotbx.pdb.input(lines=names_str.split('\n'), source_info=None).construct_hierarchy()
  atoms = h.atoms()
  fields = [LI.probe_atom_field(a) for a in atoms]
  by_name = dict([((a.name.strip(), a.parent().altloc), f) for a, f in zip(atoms, fields)])
  assert by_name[('HD21', 'A')] == ' A1001BASN HD21A', by_name
  assert by_name[('HD21', 'B')] == ' A1001BASN HD21B'
  assert by_name[('CA', '')] == ' B  52AGLY  CA  '
  assert by_name[('HD21', 'A')].split() == ['A1001BASN', 'HD21A']  # why not split
  for a, f in zip(atoms, fields):
    assert LI.probe_key(f) == LI.atom_key(a), (f, LI.probe_key(f), LI.atom_key(a))
  assert LI.probe_key(by_name[('HD21', 'A')]) == ('A', '1001B', 'ASN', 'HD21', 'A')
  assert LI.probe_key(by_name[('CA', '')]) == ('B', '52A', 'GLY', 'CA', '')
  table = LI.probe_atom_table(atoms)
  assert [table[f] for f in fields] == [a.i_seq for a in atoms]

def exercise_probe_names_in_run():
  '''
  The same on a probe2 run: ligand A1001, H-bond partner 52A (insertion code),
  vdW partner EDO 4 in two altlocs; every field maps to the atom it names.
  '''
  lines = []
  for l in model_str.split('\n'):
    if l[17:26] == 'EDO A   1':
      l = l[:22] + '1001' + l[26:]
    elif l[17:26] == 'EDO A   2':
      l = l[:22] + '  52A' + l[27:]
    elif l[17:26] == 'EDO A   4':
      lb = l[:16] + 'B' + l[17:46] + '%8.3f' % (float(l[46:54]) - 0.3) + l[54:]
      l = l[:16] + 'A' + l[17:54] + '  0.50' + l[60:]
      lb = lb[:54] + '  0.50' + lb[60:]
      lines.append(l)
      l = lb
    lines.append(l)
  model = get_model(lines)
  m = get_manager(model, sel='chain A and resseq 1001')
  atoms = model.get_hierarchy().atoms()
  assert m.probe_unmapped == set()
  table = LI.probe_atom_table(atoms)
  fields = set()
  for line in m.probe_output.splitlines():
    f = line.split(':')
    if len(f) > 4 and f[2] in LI.probe_classes:
      fields.update([f[3], f[4]])
  for f in fields:
    assert LI.probe_key(f) == LI.atom_key(atoms[table[f]]), f
  assert any([f.startswith(' A1001 EDO') for f in fields])
  assert any([f.startswith(' A  52AEDO') for f in fields])
  alts = set([f[-1] for f in fields if f.startswith(' A   4 EDO')])
  assert alts == set(['A', 'B']), alts
  hb = entries_by_type(m)['hbond']
  assert [short(l) for l in hb[0]['labels']] == ['EDO 1001 O1', 'EDO 1001 HO1',
    'EDO 52A O1'], hb
  assert hb[0]['cross_check'] == 'pnp and probe2'
  vdw_alts = set([l.split()[-1] for e in m.entries if e['type'] == 'vdw'
    for l in e['labels'] if ' EDO 4 ' in l])
  assert vdw_alts == set(['A', 'B']), vdw_alts

def exercise_dot_patches():
  '''Two clusters 2 A apart; a chain of points 0.4 A apart links them at 0.5 A only.'''
  a = [(x * 0.25, y * 0.25, 0) for x in range(4) for y in range(4)]
  b = [(p[0] + 2.75, p[1], p[2]) for p in a]
  groups = LI.dot_patches(flex.vec3_double(a + b), 0.5)
  assert [len(g) for g in groups] == [16, 16]
  assert sorted(groups[0] + groups[1]) == list(range(32))
  chain = [(0.75 + 0.4 * k, 0, 0) for k in range(1, 5)]
  groups = LI.dot_patches(flex.vec3_double(a + b + chain), 0.5)
  assert [len(g) for g in groups] == [36]
  groups = LI.dot_patches(flex.vec3_double(a + b + chain), 0.35)
  assert [len(g) for g in groups] == [16, 16, 1, 1, 1, 1]
  assert LI.dot_patches(flex.vec3_double(), 0.5) == []

def exercise_patches(m):
  '''
  Patches on the model: each patch sees one environment copy here. Dots are
  at their location on the ligand atom's surface (vdW radius from the atom), not
  at the spike end, which overlap dots push into the neighbour.
  '''
  lig = set(m.ligand_isel)
  out = [d for d in m.probe_dots if d['source'] in lig and d['target'] not in lig]
  assert sum([p['n_dots'] for p in m.patches]) == len(out)
  residues = [p['residues'] for p in m.patches]
  assert ['A EDO 3'] in residues and ['A EDO 4'] in residues
  big = [p for p in m.patches if p['residues'] == ['A EDO 3']][0]
  assert big['ligand_atoms'] == ['A EDO 1 H11']
  assert big['dots']['bo'] > 0 and big['dots']['so'] > 0
  assert sum(big['dots'].values()) == big['n_dots']
  assert sum([q['dots'] for q in big['pairs']]) == big['n_dots']
  atoms = m.model.get_hierarchy().atoms()
  h11 = [d for d in out if atoms[d['source']].name.strip() == 'H11']
  r_loc = [atoms[d['source']].distance(d['loc']) for d in h11]
  r_spike = [atoms[d['source']].distance(d['spike']) for d in h11]
  assert max(r_loc) - min(r_loc) < 0.01, (min(r_loc), max(r_loc))
  assert max(r_spike) - min(r_spike) > 0.1
  xyz = flex.vec3_double([d['loc'] for d in out])
  groups = LI.dot_patches(xyz, 0.5)
  assert [len(g) for g in groups] == [p['n_dots'] for p in m.patches]
  for g, p in zip(groups, m.patches):
    assert approx_equal(xyz.select(flex.size_t(g)).mean(), p['center'], eps=1.e-9)

def exercise_library_vs_command_line(model):
  '''
  probe2 as a library call gives the same dots as mmtbx.probe2 on the file (both
  from the file: the PDB format rounds coordinates to 0.001 A).
  '''
  src, tgt = '(%s)' % LIG_SEL, 'not (%s)' % LIG_SEL
  fn = 'tst_ligand_interactions.pdb'
  with open(fn, 'w') as f:
    f.write(model.model_as_pdb())
  from_file = mmtbx.model.manager(model_input=iotbx.pdb.input(file_name=fn),
    log=null_out())
  from_file.process(make_restraints=True)
  text = LI.run_probe2(from_file, src, tgt)
  out = 'tst_ligand_interactions_probe2.txt'
  cmd = ('mmtbx.probe2 %s approach=both "source_selection=%s" "target_selection=%s"'
    ' output.format=raw output.filename=%s --overwrite' % (fn, src, tgt, out))
  r = easy_run.fully_buffered(cmd)
  assert r.return_code == 0, '\n'.join(r.stderr_lines)
  with open(out) as f:
    cli = f.read()
  lib_lines = [l for l in text.splitlines() if l.strip()]
  cli_lines = [l for l in cli.splitlines() if l.strip()]
  assert len(lib_lines) > 100
  assert lib_lines == cli_lines

def exercise_ligand_overlaps(model):
  '''
  ligand_overlaps is what validate_ligands' get_overlaps reports; its records
  point to the right atoms of the full model. EDO 0 comes first in the file but
  lies outside the 3 A region, so model_within's atom numbers differ from the
  full model's by 10; the assertions fail if the mapping back to the full model
  is missing.
  '''
  from mmtbx.validation import validate_ligands
  isel = model.selection(LIG_SEL).iselection()
  ov = LI.ligand_overlaps(model, LIG_SEL)
  params = validate_ligands.master_params().extract().validate_ligands
  lr = validate_ligands.ligand_result(model=model, fmodel=None, map_manager=None,
    ligand_isel=isel, sel_str=LIG_SEL, params=params)
  vl = lr.get_overlaps()
  for k in ('n_clashes', 'clashscore', 'n_clashes_sym', 'clashes_str', 'n_hbonds',
            'clash_records', 'hbond_records'):
    assert getattr(vl, k) == getattr(ov, k), k
  atoms = model.get_hierarchy().atoms()
  assert ov.n_clashes == len(ov.clash_records) == 1
  assert ov.n_hbonds == len(ov.hbond_records) == 1
  r = ov.clash_records[0]
  assert r['symop'] == ''
  names = [atoms[r['i_seq']].id_str(), atoms[r['j_seq']].id_str()]
  assert 'H11 EDO A   1' in names[0] and 'H12 EDO A   3' in names[1], names
  assert approx_equal(atoms[r['i_seq']].distance(atoms[r['j_seq']]), r['distance'],
    eps=1.e-6)
  assert '%.2f' % r['distance'] in ov.clashes_str
  assert '%.2f' % r['overlap'] in ov.clashes_str
  for field in ('H11 EDO A   1', 'H12 EDO A   3'):
    assert field in ov.clashes_str
  h = ov.hbond_records[0]
  assert [atoms[h[k]].name.strip() for k in ('d_seq', 'h_seq', 'a_seq')] == \
    ['O1', 'HO1', 'O1']
  assert [atoms[h[k]].parent().parent().resseq.strip() for k in
    ('d_seq', 'h_seq', 'a_seq')] == ['1', '1', '2']
  assert approx_equal(atoms[h['h_seq']].distance(atoms[h['a_seq']]), h['d_HA'],
    eps=1.e-6)

def exercise_imports():
  '''
  A fresh process that imports ligand_interactions and runs it loads neither
  validate_ligands nor mmtbx.nci.hbond through the module itself (pnp imports
  mmtbx.nci.hbond), nor anything from nci_analysis.
  '''
  code = '\n'.join([
    'from __future__ import print_function',
    'import sys',
    'import mmtbx.validation.ligand_interactions',
    'from mmtbx.regression import tst_ligand_interactions as T',
    'T.get_manager(T.get_model())',
    'print(" ".join(sorted([k for k in sys.modules if',
    '  k == "mmtbx.validation.validate_ligands" or',
    '  k.split(".")[0] == "nci_analysis"])))'])
  fn = 'tst_ligand_interactions_imports.py'
  with open(fn, 'w') as f:
    f.write(code + '\n')
  r = easy_run.fully_buffered('libtbx.python %s' % fn)
  assert r.return_code == 0, '\n'.join(r.stderr_lines)
  assert [l for l in r.stdout_lines if l.strip()] == [], r.stdout_lines
  src = open(LI.__file__.replace('.pyc', '.py')).read()
  for name in ('validate_ligands', 'mmtbx.nci', 'nci_analysis'):
    assert ('import %s' % name) not in src and ('from %s' % name) not in src, name

def exercise_errors(model):
  isel = model.selection(LIG_SEL).iselection()
  try:
    LI.manager(model, isel, 'chain A and resseq 2').run()
  except Sorry:
    pass
  else:
    raise AssertionError('mismatched sel_str accepted')
  no_h = model.select(~model.get_hierarchy().atom_selection_cache().selection(
    'element H'))
  try:
    LI.manager(no_h, no_h.selection(LIG_SEL).iselection(), LIG_SEL)
  except Sorry:
    pass
  else:
    raise AssertionError('model without H accepted')
  # one molecule: two separate EDO residues refused
  sel = 'chain A and resname EDO and (resseq 1 or resseq 3)'
  try:
    LI.manager(model, model.selection(sel).iselection(), sel).run()
  except Sorry as e:
    assert 'not one molecule' in str(e) and '2 fragments' in str(e), str(e)
  else:
    raise AssertionError('two molecules accepted')

# ------------------------------------------------------------------------------

def run():
  model = get_model()
  exercise_dot_patches()
  exercise_probe_names()
  m = exercise_entries_and_counts(model)
  exercise_patches(m)
  exercise_pair_class_order(model)
  exercise_internal()
  exercise_symmetry()
  exercise_donor_conformers()
  exercise_probe_names_in_run()
  exercise_library_vs_command_line(model)
  exercise_ligand_overlaps(model)
  exercise_imports()
  exercise_errors(model)
  print('OK')

if __name__ == '__main__':
  run()
