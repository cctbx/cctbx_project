from __future__ import absolute_import, division, print_function

import mmtbx
from libtbx.utils import null_out
import mmtbx.monomer_library.pdb_interpretation
import mmtbx
import mmtbx.secondary_structure
import iotbx
import libtbx.load_env

pdb_str = """\
CRYST1   18.879   16.714   25.616  90.00  90.00  90.00 P 1
ATOM      1  N   ALA     1      13.515   7.809  20.095  1.00  0.00           N
ATOM      2  CA  ALA     1      13.087   6.532  19.536  1.00  0.00           C
ATOM      3  C   ALA     1      11.716   6.653  18.880  1.00  0.00           C
ATOM      4  O   ALA     1      11.425   5.972  17.896  1.00  0.00           O
ATOM      5  CB  ALA     1      13.065   5.461  20.616  1.00  0.00           C
ATOM      6  N   ALA     2      10.876   7.524  19.431  1.00  0.00           N
ATOM      7  CA  ALA     2       9.535   7.735  18.900  1.00  0.00           C
ATOM      8  C   ALA     2       9.565   8.647  17.678  1.00  0.00           C
ATOM      9  O   ALA     2       8.787   8.471  16.741  1.00  0.00           O
ATOM     10  CB  ALA     2       8.626   8.316  19.973  1.00  0.00           C
ATOM     11  N   ALA     3      10.469   9.622  17.697  1.00  0.00           N
ATOM     12  CA  ALA     3      10.606  10.565  16.593  1.00  0.00           C
ATOM     13  C   ALA     3      11.132   9.871  15.340  1.00  0.00           C
ATOM     14  O   ALA     3      10.687  10.157  14.228  1.00  0.00           O
ATOM     15  CB  ALA     3      11.520  11.714  16.987  1.00  0.00           C
ATOM     16  N   ALA     4      12.081   8.960  15.530  1.00  0.00           N
ATOM     17  CA  ALA     4      12.656   8.209  14.421  1.00  0.00           C
ATOM     18  C   ALA     4      11.624   7.266  13.812  1.00  0.00           C
ATOM     19  O   ALA     4      11.570   7.090  12.595  1.00  0.00           O
ATOM     20  CB  ALA     4      13.879   7.432  14.883  1.00  0.00           C
ATOM     21  N   ALA     5      10.805   6.663  14.669  1.00  0.00           N
ATOM     22  CA  ALA     5       9.753   5.761  14.218  1.00  0.00           C
ATOM     23  C   ALA     5       8.662   6.528  13.481  1.00  0.00           C
ATOM     24  O   ALA     5       8.107   6.045  12.494  1.00  0.00           O
ATOM     25  CB  ALA     5       9.165   5.000  15.396  1.00  0.00           C
ATOM     26  N   ALA     6       8.360   7.728  13.967  1.00  0.00           N
ATOM     27  CA  ALA     6       7.358   8.580  13.338  1.00  0.00           C
ATOM     28  C   ALA     6       7.860   9.106  11.998  1.00  0.00           C
ATOM     29  O   ALA     6       7.078   9.322  11.072  1.00  0.00           O
ATOM     30  CB  ALA     6       6.986   9.732  14.257  1.00  0.00           C
ATOM     31  N   ALA     7       9.169   9.311  11.903  1.00  0.00           N
ATOM     32  CA  ALA     7       9.781   9.787  10.668  1.00  0.00           C
ATOM     33  C   ALA     7       9.912   8.655   9.655  1.00  0.00           C
ATOM     34  O   ALA     7       9.905   8.886   8.446  1.00  0.00           O
ATOM     35  CB  ALA     7      11.141  10.405  10.952  1.00  0.00           C
ATOM     36  N   ALA     8      10.030   7.429  10.157  1.00  0.00           N
ATOM     37  CA  ALA     8      10.152   6.258   9.297  1.00  0.00           C
ATOM     38  C   ALA     8       8.788   5.809   8.786  1.00  0.00           C
ATOM     39  O   ALA     8       8.667   5.312   7.666  1.00  0.00           O
ATOM     40  CB  ALA     8      10.839   5.124  10.041  1.00  0.00           C
ATOM     41  N   ALA     9       7.762   5.988   9.613  1.00  0.00           N
ATOM     42  CA  ALA     9       6.405   5.603   9.243  1.00  0.00           C
ATOM     43  C   ALA     9       5.816   6.576   8.228  1.00  0.00           C
ATOM     44  O   ALA     9       5.000   6.195   7.389  1.00  0.00           O
ATOM     45  CB  ALA     9       5.521   5.526  10.478  1.00  0.00           C
ATOM     46  N   ALA    10       6.235   7.835   8.309  1.00  0.00           N
ATOM     47  CA  ALA    10       5.751   8.864   7.397  1.00  0.00           C
ATOM     48  C   ALA    10       6.434   8.760   6.038  1.00  0.00           C
ATOM     49  O   ALA    10       5.773   8.734   5.000  1.00  0.00           O
ATOM     50  CB  ALA    10       5.966  10.246   7.995  1.00  0.00           C
TER
END
"""

def exercise_geo_out():
  log = null_out()
  # import sys
  # log = sys.stdout
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  pdb_int_pars = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  pdb_int_pars.secondary_structure.enabled=True
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=pdb_int_pars,
    raw_records=pdb_str.split('\n'),
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  acp = processed_pdb_file.all_chain_proxies
  source = acp.pdb_inp.source_info()
  defpars = mmtbx.secondary_structure.sec_str_master_phil.fetch()
  custom_pars = defpars.fetch(source = iotbx.phil.parse(
      "h_bond_restraints.remove_outliers=False\n"))
  custom_pars_extr = custom_pars.extract()
  ss_manager = mmtbx.secondary_structure.manager(
      pdb_hierarchy=acp.pdb_hierarchy,
      params=custom_pars_extr.secondary_structure)
  grm = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      assume_hydrogens_all_missing = True,
      )
  n_hbp = grm.get_n_hbond_proxies()
  assert n_hbp == 4

if (__name__ == "__main__"):
  if libtbx.env.has_module("ksdssp"):
    # removed check for old h-bond proxies
    exercise_geo_out()
    print("OK")
  else :
    print("skipping test, ksdssp not available")
