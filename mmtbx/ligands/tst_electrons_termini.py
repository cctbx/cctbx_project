from __future__ import absolute_import, division, print_function
'''
Charges at peptide chain ends. The restraint files describe in-chain residues, with
no OXT or terminal H, so what is at a chain end has to come from the model: OXT
closes the valence of C and extra H on N make a free N-terminus, neither needing a
cap for an absent neighbour residue. Hydrogen atoms from phenix.ready_set.
'''
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from iotbx.cli_parser import run_program
from libtbx.utils import null_out
from mmtbx.ligands import electrons
from mmtbx.ligands import rdkit_utils

ala_ala_ala = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  ALA A   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   ALA A   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   ALA A   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      5  CB  ALA A   1      29.170  -1.266   3.162  1.00  0.00           C
ATOM      6  H1  ALA A   1      27.235   0.288   2.456  1.00  0.00           H
ATOM      7  H2  ALA A   1      27.942   1.373   1.796  1.00  0.00           H
ATOM      8  H3  ALA A   1      27.977   0.034   1.233  1.00  0.00           H
ATOM      9  HA  ALA A   1      29.128   0.745   3.579  1.00  0.00           H
ATOM     10  HB1 ALA A   1      29.966  -1.454   3.683  1.00  0.00           H
ATOM     11  HB2 ALA A   1      28.383  -1.470   3.690  1.00  0.00           H
ATOM     12  HB3 ALA A   1      29.172  -1.813   2.361  1.00  0.00           H
ATOM     13  N   ALA A   2      31.123   1.587   2.474  1.00  0.00           N
ATOM     14  CA  ALA A   2      32.355   2.031   1.832  1.00  0.00           C
ATOM     15  C   ALA A   2      33.552   1.851   2.758  1.00  0.00           C
ATOM     16  O   ALA A   2      33.675   2.539   3.772  1.00  0.00           O
ATOM     17  CB  ALA A   2      32.231   3.487   1.408  1.00  0.00           C
ATOM     18  H   ALA A   2      30.901   2.042   3.170  1.00  0.00           H
ATOM     19  HA  ALA A   2      32.499   1.484   1.044  1.00  0.00           H
ATOM     20  HB1 ALA A   2      33.055   3.772   0.983  1.00  0.00           H
ATOM     21  HB2 ALA A   2      31.496   3.581   0.782  1.00  0.00           H
ATOM     22  HB3 ALA A   2      32.063   4.038   2.188  1.00  0.00           H
ATOM     23  N   ALA A   3      34.434   0.922   2.403  1.00  0.00           N
ATOM     24  CA  ALA A   3      35.624   0.650   3.202  1.00  0.00           C
ATOM     25  C   ALA A   3      36.893   0.975   2.422  1.00  0.00           C
ATOM     26  O   ALA A   3      37.227   0.298   1.450  1.00  0.00           O
ATOM     27  CB  ALA A   3      35.637  -0.806   3.646  1.00  0.00           C
ATOM     28  H   ALA A   3      34.361   0.435   1.698  1.00  0.00           H
ATOM     29  HA  ALA A   3      35.598   1.219   3.987  1.00  0.00           H
ATOM     30  HB1 ALA A   3      36.432  -0.976   4.176  1.00  0.00           H
ATOM     31  HB2 ALA A   3      34.848  -0.989   4.180  1.00  0.00           H
ATOM     32  HB3 ALA A   3      35.640  -1.382   2.866  1.00  0.00           H
ATOM    900  OXT ALA A   3      37.603   1.941   2.774  1.00  0.00           O
ATOM    901  HXT ALA A   3      38.429   2.308   2.422  1.00  0.00           H
'''

pro_ala_ala = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  PRO A   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   PRO A   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   PRO A   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      5  CB  PRO A   1      29.067  -1.304   2.986  1.00  0.00           C
ATOM      6  CG  PRO A   1      27.616  -1.613   2.881  1.00  0.00           C
ATOM      7  CD  PRO A   1      27.077  -0.688   1.829  1.00  0.00           C
ATOM      8  H2  PRO A   1      28.195   0.791   1.212  1.00  0.00           H
ATOM      9  H3  PRO A   1      27.520   1.137   2.368  1.00  0.00           H
ATOM     10  HA  PRO A   1      29.192   0.709   3.601  1.00  0.00           H
ATOM     11  HB2 PRO A   1      29.580  -1.785   2.318  1.00  0.00           H
ATOM     12  HB3 PRO A   1      29.421  -1.559   3.852  1.00  0.00           H
ATOM     13  HG2 PRO A   1      27.475  -2.540   2.635  1.00  0.00           H
ATOM     14  HG3 PRO A   1      27.168  -1.473   3.730  1.00  0.00           H
ATOM     15  HD2 PRO A   1      27.140  -1.072   0.940  1.00  0.00           H
ATOM     16  HD3 PRO A   1      26.143  -0.471   1.979  1.00  0.00           H
ATOM     17  N   ALA A   2      31.123   1.587   2.474  1.00  0.00           N
ATOM     18  CA  ALA A   2      32.355   2.031   1.832  1.00  0.00           C
ATOM     19  C   ALA A   2      33.552   1.851   2.758  1.00  0.00           C
ATOM     20  O   ALA A   2      33.675   2.539   3.772  1.00  0.00           O
ATOM     21  CB  ALA A   2      32.231   3.487   1.408  1.00  0.00           C
ATOM     22  H   ALA A   2      30.901   2.042   3.170  1.00  0.00           H
ATOM     23  HA  ALA A   2      32.499   1.484   1.044  1.00  0.00           H
ATOM     24  HB1 ALA A   2      33.055   3.772   0.983  1.00  0.00           H
ATOM     25  HB2 ALA A   2      31.496   3.581   0.782  1.00  0.00           H
ATOM     26  HB3 ALA A   2      32.063   4.038   2.188  1.00  0.00           H
ATOM     27  N   ALA A   3      34.434   0.922   2.403  1.00  0.00           N
ATOM     28  CA  ALA A   3      35.624   0.650   3.202  1.00  0.00           C
ATOM     29  C   ALA A   3      36.893   0.975   2.422  1.00  0.00           C
ATOM     30  O   ALA A   3      37.227   0.298   1.450  1.00  0.00           O
ATOM     31  CB  ALA A   3      35.637  -0.806   3.646  1.00  0.00           C
ATOM     32  H   ALA A   3      34.361   0.435   1.698  1.00  0.00           H
ATOM     33  HA  ALA A   3      35.598   1.219   3.987  1.00  0.00           H
ATOM     34  HB1 ALA A   3      36.432  -0.976   4.176  1.00  0.00           H
ATOM     35  HB2 ALA A   3      34.848  -0.989   4.180  1.00  0.00           H
ATOM     36  HB3 ALA A   3      35.640  -1.382   2.866  1.00  0.00           H
'''

# name, model, atoms removed, charge of each residue, residues with a missing cap
cases = [
  ('NH3+ COOH',    ala_ala_ala, [],                   [1, 0, 0],  []),
  ('NH3+ COO-',    ala_ala_ala, ['HXT'],              [1, 0, -1], []),
  ('NH2 COO-',     ala_ala_ala, ['H3', 'HXT'],        [0, 0, -1], []),
  # one H on N and no OXT: chain breaks, the absent residues are capped
  ('capped',       ala_ala_ala, ['H1', 'H3', 'OXT', 'HXT'], [0, 0, 0], [0, 2]),
  ('proline NH2+', pro_ala_ala, [],                   [1, 0, 0],  [2]),
  ('proline NH',   pro_ala_ala, ['H3'],               [0, 0, 0],  [2]),
]

def get_pdb_str(pdb_str, removed):
  return '\n'.join([line for line in pdb_str.split('\n')
    if line[12:16].strip() not in removed])

def exercise_residues(name, pdb_str, charges, capped):
  inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  model.process(make_restraints=True)
  for i, rg in enumerate(model.get_hierarchy().residue_groups()):
    r = rdkit_utils.residue_rigid_components(model=model, residue_group=rg).molecule
    assert r.ok, (name, i, r.reason)
    assert r.total_charge == charges[i], (name, i, r.total_charge)
    n_caps = len([c for c in r.caps if c['kind'] == 'missing'])
    assert n_caps == (i in capped), (name, i, r.caps)

def exercise_program(name, pdb_str, charges):
  model_filename = 'tst_electrons_termini.pdb'
  with open(model_filename, 'w') as f:
    f.write(pdb_str)
  log = StringIO()
  rc = run_program(program_class=electrons.Program,
                   args=[model_filename],
                   logger=log,
                  )
  assert rc.total_charge == sum(charges), (name, rc.total_charge)
  # the listing names the atoms of each residue, not those of the first residue
  listed = [line for line in log.getvalue().split('\n') if 'charge=' in line]
  assert len(listed) == len([q for q in charges if q]), (name, listed)
  for resseq, charge in enumerate(charges):
    if charge:
      found = [line for line in listed if ' A%4d ' % (resseq + 1) in line]
      assert len(found) == 1, (name, resseq, listed)
      assert found[0].endswith('charge=%2d' % charge), (name, found)

def main():
  for name, pdb_str, removed, charges, capped in cases:
    pdb_str = get_pdb_str(pdb_str, removed)
    exercise_residues(name, pdb_str, charges, capped)
    exercise_program(name, pdb_str, charges)
  print('OK')

if __name__ == '__main__':
  main()
