from __future__ import absolute_import, division, print_function
import sys
import mmtbx.model
from iotbx import pdb
from libtbx.utils import null_out

pdb_str = '''
CRYST1   74.901   73.950   81.107  90.00  90.00  90.00 P 1
SCALE1      0.013351  0.000000  0.000000        0.00000
SCALE2      0.000000  0.013523  0.000000        0.00000
SCALE3      0.000000  0.000000  0.012329        0.00000
ATOM      1  N   ALA A  10      20.479  70.180  20.146  1.00 35.17           N
ATOM      2  CA  ALA A  10      21.143  71.035  21.111  1.00 34.92           C
ATOM      3  C   ALA A  10      21.840  70.184  22.154  1.00 43.85           C
ATOM      4  O   ALA A  10      21.298  69.164  22.607  1.00 46.47           O
ATOM      5  CB  ALA A  10      20.146  71.953  21.787  1.00 44.05           C
ATOM      6  H   ALA A  10      19.794  69.554  20.569  0.00 35.17           H
ATOM      7  HA  ALA A  10      21.895  71.643  20.608  0.00 34.92           H
ATOM      8  HB1 ALA A  10      20.673  72.582  22.505  0.00 44.05           H
ATOM      9  HB2 ALA A  10      19.666  72.574  21.031  0.00 44.05           H
ATOM     10  HB3 ALA A  10      19.399  71.348  22.301  0.00 44.05           H
TER
ATOM         H2  ALA A  10      20.046  70.688  19.542  1.00 35.17           H
TER
ATOM         HC  ALA A  10      22.748  70.442  22.485  1.00 43.85           H
'''

def main(filename=None):
  if filename is None:
    filename = 'ala_term.pdb'
    with open(filename, 'w') as f:
      f.write(pdb_str)
  pdb_inp = pdb.input(filename)
  hierarchy = pdb_inp.construct_hierarchy()
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  params.pdb_interpretation.restraints_library.cdl = False
  print(len(hierarchy.atoms()))
  asc = hierarchy.atom_selection_cache()
  sel = asc.selection("element H or element D")
  print('sel len',len(sel))
  print(len(hierarchy.atoms()))
  hierarchy.show()

  model = mmtbx.model.manager(
    model_input               = None,
    pdb_hierarchy             = hierarchy,
    crystal_symmetry          = pdb_inp.crystal_symmetry(),
    log                       = null_out(),
    )
  print('m1',len(model.get_hierarchy().atoms()))
  for atom in model.get_hierarchy().atoms(): print(atom.format_atom_record())
  model.process(make_restraints=True,
                grm_normalization=True,
                pdb_interpretation_params = params)
  print('m2',len(model.get_hierarchy().atoms()))
  for atom in model.get_hierarchy().atoms(): print(atom.format_atom_record())
  model.idealize_h_riding()
  model.set_occupancies(0., selection=sel)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
