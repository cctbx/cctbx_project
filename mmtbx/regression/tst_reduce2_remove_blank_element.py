"""
approach=remove selects H and D by element. Where element columns are blank
it takes the elements from an interpretation of the run's model, with and
without a restraint file: H and D are removed, a heavy atom keeps the
element of its dictionary, or of its name where there is none. The caller's
model is not changed. An input that cannot be loaded (a blank element and a
name without an element guess) stops as before.
"""
import os
import iotbx.pdb
from libtbx.utils import Sorry, null_out
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import CCTBXParser
from mmtbx.programs import reduce2
from mmtbx.regression import tst_monomer_library_terminal_nh3 as tml

CRYST1 = 'CRYST1   90.000   90.000   90.000  90.00  90.00  90.00 P 1\n'
MERCURY = ('HETATM  999 HG    HG C   1      40.000  40.000  40.000  1.00 '
  '20.00          HG\n')
LINE = ('Blank element columns in %d atom(s): the model was interpreted for '
  'their elements; %d still blank.')

# 1pq7 ALA A 1 and SER A 2 (conformer A) with reduce2's H (place)
ala_ser = """\
ATOM      1  N   ALA A   1      20.064 -26.291   0.046  1.00  6.09           N
ATOM      2  CA  ALA A   1      21.152 -26.951   0.767  1.00  6.82           C
ATOM      3  C   ALA A   1      21.867 -26.032   1.762  1.00  6.66           C
ATOM      4  O   ALA A   1      22.541 -26.551   2.653  1.00  8.52           O
ATOM      5  CB  ALA A   1      22.184 -27.500  -0.229  1.00  8.10           C
ATOM      6  H1  ALA A   1      19.667 -26.877  -0.494  1.00  6.09           H
ATOM      7  H2  ALA A   1      19.470 -25.976   0.629  1.00  6.09           H
ATOM      8  H3  ALA A   1      20.392 -25.619  -0.436  1.00  6.09           H
ATOM      9  HA  ALA A   1      20.756 -27.676   1.276  1.00  6.82           H
ATOM     10  HB1 ALA A   1      22.541 -26.765  -0.751  1.00  8.10           H
ATOM     11  HB2 ALA A   1      22.898 -27.934   0.263  1.00  8.10           H
ATOM     12  HB3 ALA A   1      21.749 -28.140  -0.814  1.00  8.10           H
ATOM     13  N   SER A   2      21.784 -24.715   1.626  1.00  6.65           N
ATOM     14  CA  SER A   2      22.466 -23.794   2.534  1.00  6.81           C
ATOM     15  C   SER A   2      21.522 -22.817   3.228  1.00  7.07           C
ATOM     16  O   SER A   2      21.992 -21.784   3.761  1.00  8.65           O
ATOM     17  CB  SER A   2      23.532 -22.969   1.784  1.00  9.30           C
ATOM     18  OG  SER A   2      24.587 -23.797   1.419  1.00  9.73           O
ATOM     19  H   SER A   2      21.333 -24.321   1.008  1.00  6.65           H
ATOM     20  HA  SER A   2      22.885 -24.344   3.215  1.00  6.81           H
ATOM     21  HB2 SER A   2      23.137 -22.582   0.987  1.00  9.30           H
ATOM     22  HB3 SER A   2      23.862 -22.266   2.365  1.00  9.30           H
ATOM     23  HG  SER A   2      24.934 -24.139   2.103  1.00  9.73           H
"""

# 1pq7 VAL A 17, GLN A 67 and ARG A 35 with reduce2's H (add)
val = """\
ATOM     20  N   VAL A   1      -2.773   8.494 -10.089  1.00  5.27           N
ATOM     21  CA  VAL A   1      -3.068   9.061 -11.398  1.00  5.73           C
ATOM     22  C   VAL A   1      -4.514   8.748 -11.734  1.00  5.60           C
ATOM     23  O   VAL A   1      -4.966   7.602 -11.604  1.00  5.67           O
ATOM     24  CB  VAL A   1      -2.140   8.466 -12.459  1.00  6.26           C
ATOM     25  CG1 VAL A   1      -2.427   9.068 -13.828  1.00  7.88           C
ATOM     26  CG2 VAL A   1      -0.688   8.663 -12.119  1.00  8.69           C
ATOM     27  H   VAL A   1      -2.622   7.648 -10.122  1.00  5.27           H
ATOM     28  HA  VAL A   1      -2.934  10.022 -11.388  1.00  5.73           H
ATOM     29  HB  VAL A   1      -2.313   7.512 -12.485  1.00  6.26           H
ATOM     30 HG11 VAL A   1      -2.287  10.027 -13.787  1.00  7.88           H
ATOM     31 HG12 VAL A   1      -1.826   8.673 -14.479  1.00  7.88           H
ATOM     32 HG13 VAL A   1      -3.347   8.878 -14.070  1.00  7.88           H
ATOM     33 HG21 VAL A   1      -0.500   8.230 -11.272  1.00  8.69           H
ATOM     34 HG22 VAL A   1      -0.144   8.270 -12.819  1.00  8.69           H
ATOM     35 HG23 VAL A   1      -0.505   9.613 -12.052  1.00  8.69           H
"""

gln = """\
ATOM    746  N   GLN A   2       7.441   9.043  10.223  1.00  5.63           N
ATOM    747  CA  GLN A   2       6.326   9.435  11.074  1.00  5.80           C
ATOM    748  C   GLN A   2       5.034   8.999  10.438  1.00  5.03           C
ATOM    749  O   GLN A   2       4.748   9.347   9.270  1.00  5.75           O
ATOM    750  CB  GLN A   2       6.294  10.934  11.362  1.00  6.79           C
ATOM    751  CG  GLN A   2       5.234  11.296  12.397  1.00  8.45           C
ATOM    752  CD  GLN A   2       5.197  12.753  12.712  1.00  9.35           C
ATOM    753  OE1 GLN A   2       4.945  13.549  11.785  1.00 11.19           O
ATOM    754  NE2 GLN A   2       5.448  13.143  13.863  1.00 15.18           N
ATOM    755  H   GLN A   2       7.323   9.218   9.389  1.00  5.63           H
ATOM    756  HA  GLN A   2       6.432   9.005  11.937  1.00  5.80           H
ATOM    757  HB2 GLN A   2       7.159  11.211  11.704  1.00  6.79           H
ATOM    758  HB3 GLN A   2       6.093  11.411  10.542  1.00  6.79           H
ATOM    759  HG2 GLN A   2       4.362  11.042  12.056  1.00  8.45           H
ATOM    760  HG3 GLN A   2       5.421  10.817  13.220  1.00  8.45           H
ATOM    761 HE21 GLN A   2       5.643  12.573  14.477  1.00 15.18           H
ATOM    762 HE22 GLN A   2       5.427  13.983  14.047  1.00 15.18           H
"""

arg = """\
ATOM    260  N   ARG A   3       8.711  10.689   5.662  1.00  5.23           N
ATOM    261  CA  ARG A   3      10.118  10.342   5.920  1.00  5.31           C
ATOM    262  C   ARG A   3      10.893  11.642   5.842  1.00  5.81           C
ATOM    263  O   ARG A   3      10.772  12.404   4.877  1.00  6.16           O
ATOM    264  CB  ARG A   3      10.652   9.339   4.912  1.00  5.42           C
ATOM    265  CG  ARG A   3      12.066   8.888   5.233  1.00  5.75           C
ATOM    266  CD  ARG A   3      12.602   7.946   4.160  1.00  5.92           C
ATOM    267  NE  ARG A   3      14.003   7.702   4.417  1.00  5.76           N
ATOM    268  CZ  ARG A   3      14.894   7.491   3.455  1.00  5.92           C
ATOM    269  NH1 ARG A   3      14.526   7.282   2.225  1.00  7.29           N
ATOM    270  NH2 ARG A   3      16.197   7.472   3.775  1.00  6.96           N
ATOM    271  H   ARG A   3       8.590  11.244   5.016  1.00  5.23           H
ATOM    272  HA  ARG A   3      10.225   9.914   6.784  1.00  5.31           H
ATOM    273  HB2 ARG A   3      10.079   8.556   4.910  1.00  5.42           H
ATOM    274  HB3 ARG A   3      10.658   9.747   4.032  1.00  5.42           H
ATOM    275  HG2 ARG A   3      12.649   9.662   5.279  1.00  5.75           H
ATOM    276  HG3 ARG A   3      12.070   8.419   6.082  1.00  5.75           H
ATOM    277  HD2 ARG A   3      12.123   7.103   4.189  1.00  5.92           H
ATOM    278  HD3 ARG A   3      12.506   8.352   3.284  1.00  5.92           H
ATOM    279  HE  ARG A   3      14.273   7.693   5.233  1.00  5.76           H
ATOM    280 HH11 ARG A   3      13.691   7.279   2.021  1.00  7.29           H
ATOM    281 HH12 ARG A   3      15.118   7.147   1.616  1.00  7.29           H
ATOM    282 HH21 ARG A   3      16.441   7.595   4.590  1.00  6.96           H
ATOM    283 HH22 ARG A   3      16.788   7.337   3.165  1.00  6.96           H
"""

def fixture(residues, isotope='H', blank=True):
  """residues and the mercury, the H as isotope (H or D) with blank or given
  element columns"""
  result = []
  for l in residues.splitlines():
    if l[76:78].strip() == 'H':
      name = l[12:16]
      if isotope == 'D':
        i = name.index('H')
        name = name[:i] + 'D' + name[i + 1:]
      l = l[:12] + name + l[16:76] + ('  ' if blank else ' ' + isotope)
    result.append(l)
  return CRYST1 + '\n'.join(result) + '\n' + MERCURY + 'END\n'

def atoms(pdb_text):
  """(chain, resseq, residue name, atom name, element) of every atom"""
  h = iotbx.pdb.input(source_info=None,
    lines=pdb_text.splitlines()).construct_hierarchy()
  return sorted((a.parent().parent().parent().id, a.parent().parent().resseq,
    a.parent().resname, a.name, a.element.strip().upper())
    for a in h.atoms())

def program(prefix, pdb_text, cif_text=None, args=()):
  """the parser and the output names of a reduce2 approach=remove run on
  pdb_text, with cif_text as restraint file"""
  prefix = 'tst_reduce2_remove_blank_element_' + prefix
  fn, out, txt = prefix + '.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_text)
  for p in (out, txt):
    if os.path.exists(p): os.remove(p)
  files = [fn]
  if cif_text is not None:
    files.append(prefix + '.cif')
    with open(files[-1], 'w') as f:
      f.write(cif_text)
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args(files + ['approach=remove', 'output.filename=' + out,
    'output.description_file_name=' + txt, 'overwrite=True'] + list(args))
  return parser, out, txt

def execute(parser, out, txt):
  """one run: the written model, the description, the log"""
  log = StringIO()
  task = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    master_phil=parser.master_phil, logger=log)
  task.validate()
  task.run()
  task.get_results()
  with open(out) as f:
    written = f.read()
  with open(txt) as f:
    description = f.read()
  return written, description, log.getvalue()

def run(prefix, pdb_text, cif_text=None, args=()):
  return execute(*program(prefix, pdb_text, cif_text, args))

def test_blank_element_removed():
  """ALA-SER, VAL (HG21: Hg by its name) and GLN (HE21: He) with blank H
  elements, ALA-SER with blank D elements, with and without an unrelated
  restraint file: exactly the non-H atoms are written, by name and element,
  the mercury included"""
  for label, residues, isotope in (('ala_ser', ala_ser, 'H'),
      ('val', val, 'H'), ('gln', gln, 'H'), ('ala_ser', ala_ser, 'D')):
    given = atoms(fixture(residues, isotope, blank=False))
    heavy = [a for a in given if a[4] not in ('H', 'D')]
    n_h = len(given) - len(heavy)
    assert ('C', '   1', ' HG', 'HG  ', 'HG') in heavy and n_h > 7
    for cif_text in (None, tml.cif_0a0):
      prefix = '%s_%s_%s' % (label, isotope, 'cif' if cif_text else 'nocif')
      written, description, log = run(prefix, fixture(residues, isotope),
        cif_text)
      assert atoms(written) == heavy, (prefix, atoms(written))
      if cif_text is None:
        assert LINE % (n_h, 0) in log, (prefix, log)

# a ligand without a library dictionary: two C, an H with its element, and
# an atom with a blank element whose name reads as another element than the
# restraint file gives it
def j9j(name, element=''):
  return CRYST1 + """\
HETATM    1  C1  J9J A   1      10.000  10.000  10.000  1.00 20.00           C
HETATM    2  C2  J9J A   1      11.500  10.000  10.000  1.00 20.00           C
HETATM    3  H1  J9J A   1       9.600  11.000  10.000  1.00 20.00           H
HETATM    4 %-4s J9J A   1      13.000  11.500  10.000  1.00 20.00          %2s
END
""" % (name, element)

def j9j_cif(name, symbol, energy):
  return """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 J9J J9J 'test ligand' non-polymer 4 3 .
data_comp_J9J
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
 J9J C1 C CH3 0
 J9J C2 C CH2 0
 J9J H1 H HCH3 0
 J9J %s %s %s 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 J9J C1 C2 single 1.500 0.020
 J9J C1 H1 single 0.970 0.020
 J9J C2 %s single 2.100 0.020
""" % (name, symbol, energy, name)

def j9j_atoms(written):
  return [(a[3], a[4]) for a in atoms(written)]

def test_heavy_atom_by_restraint_file():
  """HO1 with a blank element, mercury in the restraint file: written as
  Hg (its name reads Ho), as with its element given"""
  cif = j9j_cif('HO1', 'HG', 'HG')
  written, description, log = run('ho1_blank', j9j('HO1'), cif)
  assert j9j_atoms(written) == [(' C1 ', 'C'), (' C2 ', 'C'),
    ('HO1 ', 'HG')], j9j_atoms(written)
  given, description, log = run('ho1_given', j9j('HO1', 'HG'), cif)
  assert j9j_atoms(given) == j9j_atoms(written), j9j_atoms(given)

def test_four_character_name_kept():
  """HG12 with a blank element in a residue without a dictionary: not an H
  by its name alone; kept, as Hg"""
  written, description, log = run('hg12_blank', j9j('HG12'))
  assert j9j_atoms(written) == [(' C1 ', 'C'), (' C2 ', 'C'),
    ('HG12', 'HG')], j9j_atoms(written)

def test_callers_model():
  """the caller's model keeps its blank element fields; later runs on it,
  without the first run's restraint file (the caller's monomer library
  server in use before the first run) and with one that gives HG12 another
  element, equal the runs on a model read again"""
  as_h, as_hg = j9j_cif('HG12', 'H', 'HCH2'), j9j_cif('HG12', 'HG', 'HG')
  parser, out, txt = program('reused', j9j('HG12'), as_h)
  model = parser.data_manager.get_model()
  model.get_mon_lib_srv()
  def elements():
    return [a.element for a in model.get_hierarchy().atoms()]
  before = elements()
  assert before.count('  ') == 1, before
  written, description, log = execute(parser, out, txt)
  assert j9j_atoms(written) == [(' C1 ', 'C'), (' C2 ', 'C')], \
    j9j_atoms(written)
  assert elements() == before, elements()
  dm = parser.data_manager
  for name in dm.get_restraint_names():
    dm.remove_restraint(name)
  assert dm.get_model() is model
  withdrawn, description, log = execute(parser, out, txt)
  fresh, description, log = run('fresh_nocif', j9j('HG12'))
  assert j9j_atoms(withdrawn) == j9j_atoms(fresh) == [(' C1 ', 'C'),
    (' C2 ', 'C'), ('HG12', 'HG')], (j9j_atoms(withdrawn), j9j_atoms(fresh))
  assert elements() == before, elements()
  fn = 'tst_reduce2_remove_blank_element_reused_hg.cif'
  with open(fn, 'w') as f:
    f.write(as_hg)
  dm.process_restraint_file(fn)
  assert dm.get_model() is model
  second, description, log = execute(parser, out, txt)
  fresh, description, log = run('fresh', j9j('HG12'), as_hg)
  assert j9j_atoms(second) == j9j_atoms(fresh) == [(' C1 ', 'C'),
    (' C2 ', 'C'), ('HG12', 'HG')], (j9j_atoms(second), j9j_atoms(fresh))
  assert elements() == before, elements()

def test_model_id():
  """a blank element only in model 2: with model_id=1 nothing is said
  about it, with model_id=2 the one line counts it"""
  body = lambda element: ''.join(l + '\n' for l in j9j('HG12',
    element).splitlines()[1:-1])
  two = CRYST1 + 'MODEL        1\n' + body('HG') + 'ENDMDL\n' + \
    'MODEL        2\n' + body('') + 'ENDMDL\nEND\n'
  for model_id, said in ((1, False), (2, True)):
    written, description, log = run('model_%d' % model_id, two,
      args=['model_id=%d' % model_id])
    assert j9j_atoms(written) == [(' C1 ', 'C'), (' C2 ', 'C'),
      ('HG12', 'HG')], (model_id, j9j_atoms(written))
    assert ('Blank element' in log) == said, (model_id, log)
    assert (LINE % (1, 1) in log) == said, (model_id, log)
    assert 'Blank element' not in description, description

def test_not_loaded():
  """ARG with blank H elements: HH11 has no element guess, the model is
  not loaded; the stop and its message are unchanged, with and without a
  restraint file"""
  for cif_text in (None, tml.cif_0a0):
    try:
      run('arg_%s' % ('cif' if cif_text else 'nocif'), fixture(arg), cif_text)
    except Sorry as e:
      assert str(e).startswith('Unknown chemical element type:'), str(e)
      assert 'HH11 ARG A   3' in str(e), str(e)
    else:
      raise AssertionError('no stop')

if __name__ == '__main__':
  test_blank_element_removed()
  test_heavy_atom_by_restraint_file()
  test_four_character_name_kept()
  test_callers_model()
  test_model_id()
  test_not_loaded()
  print('OK')
