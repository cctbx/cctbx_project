from __future__ import absolute_import, division, print_function
import iotbx.pdb
from cctbx.array_family import flex
from libtbx import easy_run
from cctbx import miller
from libtbx.test_utils import approx_equal

pdb_str = """
CRYST1   73.110   85.164   25.674  90.00  90.00  90.00 P 1
MODEL        1
ATOM      1  N   MET A   1       6.018   9.573   5.000  1.00 30.00           N
ATOM      2  CA  MET A   1       6.167   8.775   6.243  1.00 30.00           C
ATOM      3  C   MET A   1       5.758   9.566   7.448  1.00 30.00           C
ATOM      4  O   MET A   1       5.000   9.074   8.283  1.00 30.00           O
ATOM      5  CB  MET A   1       7.620   8.264   6.399  1.00 30.00           C
ATOM      6  CG  MET A   1       7.823   7.287   7.573  1.00 30.00           C
ATOM      7  SD  MET A   1       6.760   5.817   7.493  1.00 30.00           S
ATOM      8  CE  MET A   1       7.423   5.000   8.970  1.00 30.00           C
TER
ENDMDL
MODEL        2
ATOM      1  N   MET A   1      29.179  17.377   7.934  1.00 30.00           N
ATOM      2  CA  MET A   1      29.606  16.685   9.177  1.00 30.00           C
ATOM      3  C   MET A   1      28.940  17.277  10.381  1.00 30.00           C
ATOM      4  O   MET A   1      28.410  16.545  11.217  1.00 30.00           O
ATOM      5  CB  MET A   1      31.145  16.731   9.331  1.00 30.00           C
ATOM      6  CG  MET A   1      31.686  15.890  10.506  1.00 30.00           C
ATOM      7  SD  MET A   1      31.222  14.137  10.425  1.00 30.00           S
ATOM      8  CE  MET A   1      32.135  13.614  11.902  1.00 30.00           C
TER
ENDMDL
MODEL        3
ATOM      1  N   MET A   1      47.992  33.026  10.855  1.00 30.00           N
ATOM      2  CA  MET A   1      48.638  32.534  12.098  1.00 30.00           C
ATOM      3  C   MET A   1      47.804  32.846  13.303  1.00 30.00           C
ATOM      4  O   MET A   1      47.573  31.973  14.138  1.00 30.00           O
ATOM      5  CB  MET A   1      50.059  33.129  12.253  1.00 30.00           C
ATOM      6  CG  MET A   1      50.865  32.539  13.426  1.00 30.00           C
ATOM      7  SD  MET A   1      51.062  30.736  13.345  1.00 30.00           S
ATOM      8  CE  MET A   1      52.103  30.574  14.821  1.00 30.00           C
TER
ENDMDL
MODEL        4
ATOM      1  N   MET A   1      59.934  54.401  13.778  1.00 30.00           N
ATOM      2  CA  MET A   1      60.716  54.173  15.021  1.00 30.00           C
ATOM      3  C   MET A   1      59.826  54.166  16.226  1.00 30.00           C
ATOM      4  O   MET A   1      59.923  53.268  17.062  1.00 30.00           O
ATOM      5  CB  MET A   1      61.828  55.239  15.176  1.00 30.00           C
ATOM      6  CG  MET A   1      62.793  54.977  16.348  1.00 30.00           C
ATOM      7  SD  MET A   1      63.625  53.365  16.266  1.00 30.00           S
ATOM      8  CE  MET A   1      64.654  53.587  17.743  1.00 30.00           C
TER
ENDMDL
MODEL        5
ATOM      1  N   MET A   1      63.411  78.601  16.713  1.00 30.00           N
ATOM      2  CA  MET A   1      64.223  78.668  17.954  1.00 30.00           C
ATOM      3  C   MET A   1      63.396  78.342  19.161  1.00 30.00           C
ATOM      4  O   MET A   1      63.810  77.538  19.995  1.00 30.00           O
ATOM      5  CB  MET A   1      64.878  80.062  18.109  1.00 30.00           C
ATOM      6  CG  MET A   1      65.874  80.164  19.281  1.00 30.00           C
ATOM      7  SD  MET A   1      67.228  78.958  19.199  1.00 30.00           S
ATOM      8  CE  MET A   1      68.110  79.536  20.674  1.00 30.00           C
TER
ENDMDL
END
"""

def r(x,y):
  x = abs(x.deep_copy()).data()
  y = abs(y.deep_copy()).data()
  n = flex.sum(flex.abs(x-y))
  d = flex.sum(flex.abs(x+y))
  return n/d

def exercise(prefix="tst_models_to_from_chains"):
  # inputs
  expected_n = 5
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  input_file_name = "%s.pdb"%prefix
  of = open(input_file_name,"w")
  print(pdb_str, file=of)
  of.close()
  mi = flex.miller_index(((0,0,1), ))
  ms = miller.set(xrs.crystal_symmetry(), mi,  anomalous_flag=False)
  complete_set = ms.complete_set(d_min=3)
  fc = complete_set.structure_factors_from_scatterers(
    xray_structure=xrs).f_calc()
  # models -> chains
  assert not easy_run.call("phenix.models_as_chains %s"%input_file_name)
  pdb_inp = iotbx.pdb.input(file_name="chains_"+input_file_name)
  h = pdb_inp.construct_hierarchy()
  assert len(list(h.chains()))==expected_n
  assert len(list(h.models()))==1
  xrs_c = pdb_inp.xray_structure_simple()
  fc_c = complete_set.structure_factors_from_scatterers(
    xray_structure=xrs_c).f_calc()
  #
  assert approx_equal(0, r(fc, fc_c)  )

if (__name__ == "__main__"):
  exercise()
