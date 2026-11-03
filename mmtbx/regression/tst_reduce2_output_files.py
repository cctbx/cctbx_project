"""reduce2 output files.

A .pdb request for a model that does not fit PDB format (chain ID longer than
2 characters) writes mmCIF with the same base name and says so; before, the
.pdb was empty. An existing mmCIF of that name is not overwritten when
overwrite is False.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from six.moves import cStringIO as StringIO
from libtbx.utils import Sorry, null_out
from iotbx.cli_parser import CCTBXParser
from mmtbx.programs import reduce2

pdb_str = """\
ATOM      1  N   ALA A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       6.458   5.000   5.000  1.00 10.00           C
ATOM      3  C   ALA A   1       7.009   6.420   5.000  1.00 10.00           C
ATOM      4  O   ALA A   1       6.252   7.390   5.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       6.986   4.235   6.210  1.00 10.00           C
ATOM      6  N   GLY A   2       8.330   6.540   5.000  1.00 10.00           N
ATOM      7  CA  GLY A   2       8.960   7.850   5.000  1.00 10.00           C
ATOM      8  C   GLY A   2      10.470   7.740   5.000  1.00 10.00           C
ATOM      9  O   GLY A   2      11.050   6.650   5.000  1.00 10.00           O
END
"""

def write_input(prefix, chain_id):
  '''mmCIF input without a cell; returns its name and the heavy atoms.'''
  h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  h.only_chain().id = chain_id
  fn = prefix + '_in.cif'
  with open(fn, 'w') as f:
    f.write(h.as_mmcif_string())
  return fn, heavy(h)

def heavy(h):
  return dict(((a.parent().parent().resseq.strip(), a.name.strip()), a.xyz)
    for a in h.atoms() if a.element.strip() not in ('H', 'D'))

def remove(*names):
  for n in names:
    if os.path.exists(n): os.remove(n)

def make_program(args, log):
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args(args)
  return reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    master_phil=parser.master_phil, logger=log)

def run(args):
  log = StringIO()
  task = make_program(args, log)
  task.validate()
  task.run()
  return task, log.getvalue()

def test_pdb_request_does_not_fit():
  prefix = 'tst_reduce2_output_files_fallback'
  fn_in, heavy_in = write_input(prefix, 'ABC')
  pdb, cif, txt = prefix + '_out.pdb', prefix + '_out.cif', prefix + '.txt'
  remove(pdb, cif, txt)
  task, log = run([fn_in, 'approach=add', 'output.filename=' + pdb,
    'output.description_file_name=' + txt])
  assert not os.path.exists(pdb), pdb
  assert os.path.exists(cif), cif
  assert ('Model cannot be represented in PDB format (chain ID longer than 2 '
    'characters: ABC); writing mmCIF to %s instead of %s.' % (cif, pdb)) in log, log
  assert task.params.output.filename == cif, task.params.output.filename
  assert 'Wrote %s and %s' % (cif, txt) in log, log
  inp = iotbx.pdb.input(cif)
  h = inp.construct_hierarchy()
  assert [c.id for c in h.chains()] == ['ABC'], [c.id for c in h.chains()]
  heavy_out = heavy(h)
  assert sorted(heavy_out) == sorted(heavy_in), (sorted(heavy_out), sorted(heavy_in))
  for k in heavy_in:
    assert max(abs(x - y) for x, y in zip(heavy_out[k], heavy_in[k])) < 1.e-3, \
      (k, heavy_out[k], heavy_in[k])
  assert sum(1 for a in h.atoms() if a.element.strip() == 'H') > 0
  # no cell in the input: the working box is not written
  assert inp.crystal_symmetry() is None or \
    inp.crystal_symmetry().unit_cell() is None, inp.crystal_symmetry()
  with open(cif) as f:
    assert '_cell.length_a' not in f.read()

  # The mmCIF already present, no .pdb, overwrite=False: Sorry at write time
  remove(txt)
  # reduce2 turns overwrite on when this is set (test servers set it)
  overwrite_all = os.environ.pop('PHENIX_OVERWRITE_ALL', None)
  try:
    try:
      run([fn_in, 'approach=add', 'output.filename=' + pdb,
        'output.description_file_name=' + txt])
    except Sorry as e:
      assert cif in str(e), str(e)
    else:
      assert False, 'no Sorry for existing ' + cif
  finally:
    if overwrite_all is not None:
      os.environ['PHENIX_OVERWRITE_ALL'] = overwrite_all
  assert not os.path.exists(pdb), pdb

if __name__ == '__main__':
  test_pdb_request_does_not_fit()
  print('OK')
