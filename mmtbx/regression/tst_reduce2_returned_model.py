"""reduce2's returned and written model.

The run ends with a fresh, unprocessed model built from the finished
hierarchy (no second full interpretation): atoms sorted, serials and i_seq
consecutive, elements upper case, no restraints manager. Selections that do
not depend on the interpretation work on the returned model, also after the
Optimizer deleted an H; interpretation-dependent ones (backbone, sidechain)
need process(make_restraints=False) first.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.regression.tst_reduce2_place import pdb_str as peptide
from mmtbx.regression.tst_reduce2_metal_thiol import m_zn_fixture_H as cys_h
from mmtbx.regression.reduce2_neutron_utils import his_no_h

# his_no_h with a Zn ion and ANISOU records, for the mmCIF input
his_zn_anisou = """\
CRYST1   12.670   14.862   15.186  90.00  90.00  90.00 P 1
ATOM    780  N   HIS A  45       5.209   7.219   5.800  1.00  9.38           N
ANISOU  780  N   HIS A  45     1105   1206   1252     12    -30     25       N
ATOM    781  CA  HIS A  45       6.450   6.586   6.229  1.00 10.30           C
ATOM    782  C   HIS A  45       6.629   5.221   5.572  1.00 11.06           C
ATOM    783  O   HIS A  45       7.669   4.940   4.976  1.00 12.76           O
ATOM    784  CB  HIS A  45       6.482   6.444   7.752  1.00 11.71           C
ATOM    785  CG  HIS A  45       6.466   7.752   8.482  1.00 17.20           C
ATOM    786  ND1 HIS A  45       6.374   7.837   9.854  1.00 28.53           N
ATOM    787  CD2 HIS A  45       6.532   9.027   8.030  1.00 22.78           C
ATOM    788  CE1 HIS A  45       6.383   9.107  10.217  1.00 26.38           C
ATOM    789  NE2 HIS A  45       6.478   9.850   9.129  1.00 20.23           N
HETATM  790 ZN    ZN A 101       1.500  12.000  12.500  1.00 15.00          Zn
ANISOU  790 ZN    ZN A 101     1900   2000   1800      0      0      0      Zn
END
"""

def run(prefix, pdb_str, args=(), cif_input=False):
  '''reduce2 with write_files=False unless an output file name is given.
  :return: The program's results.'''
  fn = prefix + ('_in.cif' if cif_input else '_in.pdb')
  with open(fn, 'w') as f:
    if cif_input:
      f.write(iotbx.pdb.input(source_info=None, lines=pdb_str)
        .construct_hierarchy().as_mmcif_string(
          crystal_symmetry=iotbx.pdb.input(source_info=None,
            lines=pdb_str).crystal_symmetry()))
    else:
      f.write(pdb_str)
  args = list(args)
  if not any(a.startswith('output.filename=') for a in args):
    args.append('output.write_files=False')
  return run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.overwrite=True'] + args)

def check_order(model):
  atoms = model.get_hierarchy().atoms()
  n = atoms.size()
  assert list(atoms.extract_i_seq()) == list(range(n)), list(atoms.extract_i_seq())
  assert list(atoms.extract_serial()) == ['%5d' % i for i in range(1, n+1)], \
    list(atoms.extract_serial())
  h = model.get_hierarchy().deep_copy()
  h.sort_atoms_in_place()
  assert [a.id_str() for a in h.atoms()] == [a.id_str() for a in atoms]

def test_peptide_selections():
  m = run('tst_reduce2_returned_model_peptide', peptide,
    ['approach=place']).model
  assert not m.processed()
  assert not m.restraints_manager_available()
  n = m.get_number_of_atoms()
  assert n > 21, n
  assert m.selection('all').count(True) == n
  assert m.apply_selection_string('not water').get_number_of_atoms() == n
  check_order(m)
  # interpretation-dependent selections need process() first
  try:
    m.selection('backbone')
  except Exception as e:
    assert 'backbone' in str(e), str(e)
  else:
    assert False, 'backbone selection on the unprocessed model'
  m.process(make_restraints=False)
  assert m.selection('backbone').count(True) > 0

def test_remove_has_hd():
  m = run('tst_reduce2_returned_model_cys', cys_h, ['approach=remove']).model
  assert m.get_number_of_atoms() == 7, m.get_number_of_atoms()
  assert not m.has_hd()
  check_order(m)

def test_his_deletion_selection():
  r = run('tst_reduce2_returned_model_his', his_no_h,
    ['approach=add', 'add_flip_movers=True'])
  m = r.model
  names = [a.name.strip() for a in m.get_hierarchy().atoms()]
  # one of HD1/HE2 deleted by the Optimizer (His flip Mover)
  assert ('HD1' in names) != ('HE2' in names), names
  n = m.get_number_of_atoms()
  assert n == 17, n
  assert m.apply_selection_string('not water').get_number_of_atoms() == n
  assert m.apply_selection_string('element H').get_number_of_atoms() == 7
  check_order(m)

def test_written_files():
  prefix = 'tst_reduce2_returned_model_files'
  for ext in ('pdb', 'cif'):
    out, txt = '%s_out.%s' % (prefix, ext), '%s.txt' % prefix
    for f in (out, txt):
      if os.path.exists(f): os.remove(f)
    r = run(prefix, his_zn_anisou, ['approach=add', 'add_flip_movers=True',
      'output.filename=' + out, 'output.description_file_name=' + txt],
      cif_input=True)
    with open(prefix + '_in.cif') as f:
      assert ' Zn ' in f.read()
    with open(txt) as f:
      assert 'Deleting chain A HIS 45 H' in f.read()
    m = r.model
    assert not m.restraints_manager_available()
    check_order(m)
    with open(out) as f:
      text = f.read()
    assert ('ANISOU' if ext == 'pdb' else '_atom_site_anisotrop') in text
    h = iotbx.pdb.input(out).construct_hierarchy()
    atoms = h.atoms()
    assert atoms.size() == 18, atoms.size()
    assert [a.id_str() for a in atoms] == \
      [a.id_str() for a in m.get_hierarchy().atoms()]
    assert [a.serial.strip() for a in atoms] == \
      [str(i) for i in range(1, atoms.size()+1)]
    for a in atoms:
      assert a.element == a.element.upper(), a.id_str()
    zn = [a for a in atoms if a.name.strip() == 'ZN']
    assert len(zn) == 1 and zn[0].element.strip() == 'ZN'
    assert zn[0].uij_is_defined()

if __name__ == '__main__':
  test_peptide_selections()
  test_remove_has_hd()
  test_his_deletion_selection()
  test_written_files()
  print('OK')
