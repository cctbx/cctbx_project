from mmtbx import utils
import mmtbx.model
import mmtbx.restraints
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import format_cpu_times
from cStringIO import StringIO
import sys

input_model = """\
CRYST1   15.000   15.000   15.000  80.00  70.00 100.00 P 1

HETATM    1  O   HOH X   1      20.319  10.959   7.882  1.00  1.00           O
HETATM    2  H1  HOH X   2      20.158  13.476  10.359  1.00  2.00           H
HETATM    3  O   HOH X   2      19.563  13.548   9.584  1.00  1.00           O

HETATM    4  O   HOH     1      17.253  15.650  10.892  1.00  2.00           O
HETATM    5  H1  HOH     1      18.134  15.627  11.321  1.00  2.00           H
HETATM    6  H2  HOH     1      17.375  15.320   9.977  1.00  2.00           H
HETATM    7  O   HOH     2       5.150  11.586  12.474  1.00  3.00           O

ATOM      8  CB  PHE A   1      11.914  10.410  11.811  1.00  2.00           C
ATOM      9  CG  PHE A   1      11.204   9.472  12.746  1.00  2.00           C
ATOM     10  CD1 PHE A   1      10.636   8.301  12.273  1.00  2.00           C
ATOM     11  CD2 PHE A   1      11.105   9.762  14.096  1.00  2.00           C
ATOM     12  CE1 PHE A   1       9.982   7.436  13.131  1.00  2.00           C
ATOM     13  CE2 PHE A   1      10.452   8.901  14.958  1.00  2.00           C
ATOM     14  CZ  PHE A   1       9.890   7.737  14.475  1.00  2.00           C
ATOM     15  C   PHE A   1      11.828  12.443  10.351  1.00  2.00           C
ATOM     16  O   PHE A   1      11.808  12.365   9.123  1.00  2.00           O
ATOM     17  OXT PHE A   1      12.531  13.314  10.864  1.00  2.00           O
ATOM     18  N   PHE A   1       9.929  10.880  10.444  1.00  2.00           N
ATOM     19  CA  PHE A   1      11.008  11.488  11.213  1.00  2.00           C

ATOM     20  O   HOH Q   1      15.736  13.074   8.108  1.00  4.00           O
ATOM     21  D1  HOH Q   1      14.756  13.074   8.108  1.00  4.00           D

ATOM     22  O   HOH S   1      14.337  16.126   9.974  1.00  5.00           O
ATOM     23  H2  HOH S   1      15.083  16.615   9.567  1.00  5.00           H
ATOM     24  O   HOH S   2       9.491  12.823  15.494  1.00  6.00           O
END
"""

expected_result = """\
CRYST1   15.000   15.000   15.000  80.00  70.00 100.00 P 1
SCALE1      0.066667  0.011755 -0.028131        0.00000
SCALE2      0.000000  0.067695 -0.017615        0.00000
SCALE3      0.000000  0.000000  0.073308        0.00000
HETATM    1  O   HOH X   1      20.319  10.959   7.882  1.00  1.00           O
HETATM    2 H1   HOH X   1      19.339  10.959   7.882  0.00  1.00           H
HETATM    3 H2   HOH X   1      19.339  10.959   7.882  0.00  1.00           H
HETATM    4  H1  HOH X   2      20.024  13.469  10.445  1.00  2.00           H
HETATM    5  O   HOH X   2      19.563  13.548   9.584  1.00  1.00           O
HETATM    6 H2   HOH X   2      18.608  13.419   9.761  0.00  1.00           H
HETATM    7  O   HOH     1      17.253  15.650  10.892  1.00  2.00           O
HETATM    8  H1  HOH     1      18.146  15.616  11.293  1.00  2.00           H
HETATM    9  H2  HOH     1      17.357  15.408   9.948  1.00  2.00           H
HETATM   10  O   HOH     2       5.150  11.586  12.474  1.00  3.00           O
HETATM   11 H1   HOH     2       4.364  12.000  12.888  0.00  3.00           H
HETATM   12 H2   HOH     2       5.938  11.998  12.886  0.00  3.00           H
ATOM     13  CB  PHE A   1      11.914  10.410  11.811  1.00  2.00           C
ATOM     14  CG  PHE A   1      11.204   9.472  12.746  1.00  2.00           C
ATOM     15  CD1 PHE A   1      10.636   8.301  12.273  1.00  2.00           C
ATOM     16  CD2 PHE A   1      11.105   9.762  14.096  1.00  2.00           C
ATOM     17  CE1 PHE A   1       9.982   7.436  13.131  1.00  2.00           C
ATOM     18  CE2 PHE A   1      10.452   8.901  14.958  1.00  2.00           C
ATOM     19  CZ  PHE A   1       9.890   7.737  14.475  1.00  2.00           C
ATOM     20  C   PHE A   1      11.828  12.443  10.351  1.00  2.00           C
ATOM     21  O   PHE A   1      11.808  12.365   9.123  1.00  2.00           O
ATOM     22  OXT PHE A   1      12.531  13.314  10.864  1.00  2.00           O
ATOM     23  N   PHE A   1       9.929  10.880  10.444  1.00  2.00           N
ATOM     24  CA  PHE A   1      11.008  11.488  11.213  1.00  2.00           C
ATOM     25  O   HOH Q   1      15.736  13.074   8.108  1.00  4.00           O
ATOM     26 D2   HOH Q   1      16.514  12.892   7.541  0.00  4.00           D
ATOM     27  D1  HOH Q   1      14.941  12.899   7.563  1.00  4.00           D
ATOM     28  O   HOH S   1      14.337  16.126   9.974  1.00  5.00           O
ATOM     29 H1   HOH S   1      14.734  15.233  10.041  0.00  5.00           H
ATOM     30  H2  HOH S   1      14.968  16.680   9.469  1.00  5.00           H
ATOM     31  O   HOH S   2       9.491  12.823  15.494  1.00  6.00           O
ATOM     32 H1   HOH S   2       8.899  13.604  15.494  0.00  6.00           H
ATOM     33 H2   HOH S   2       8.914  12.031  15.494  0.00  6.00           H
END
"""

def run():
  pdb_file_name = "add_h_to_hoh.pdb"
  tmp_f = open(pdb_file_name, "w")
  tmp_f.write(input_model)
  tmp_f.close()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log = StringIO())
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file_name])
  xray_structure = processed_pdb_file.xray_structure()
  aal = processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  #
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = False, assume_hydrogens_all_missing = True)
  restraints_manager = mmtbx.restraints.manager(geometry = geometry)
  #
  model = mmtbx.model.manager(
    refinement_flags        = None,
    processed_pdb_files_srv = processed_pdb_files_srv,
    restraints_manager      = restraints_manager,
    xray_structure          = xray_structure,
    atom_attributes_list    = aal,
    log                     = None)
  ####
  model.add_hydrogens()
  result = model.write_pdb_file()
  ####
  result1 = []
  for r1 in result:
    if(r1.startswith("ATOM") or r1.startswith("HETATM")): result1.append(r1)
  result2 = []
  for r2 in expected_result.splitlines():
    if(r2.startswith("ATOM") or r2.startswith("HETATM")): result2.append(r2)
  assert len(result1) == len(result2)
  for r1, r2 in zip(result1, result2):
    r1 = r1[:30] + r1[55:]
    r2 = r2[:30] + r2[55:]
    assert not show_diff(r1, r2)
  ####
  cntr = 0
  xrs1 = iotbx.pdb.input(source_info = None, lines = flex.std_string(
    expected_result.splitlines())).xray_structure_simple()
  xrs2 = iotbx.pdb.input(source_info = None, lines = flex.std_string(result)
    ).xray_structure_simple()
  for s1, s2 in zip(xrs1.scatterers(), xrs2.scatterers()):
    if(s1.element_symbol().strip() not in ['H','D']):
      assert s1.element_symbol().strip() == s2.element_symbol().strip()
      assert approx_equal(s1.site, s2.site)
      cntr += 1
  assert cntr == 19
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
