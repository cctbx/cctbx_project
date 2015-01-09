from __future__ import division

import mmtbx
from libtbx.utils import null_out
import mmtbx.monomer_library.pdb_interpretation
import mmtbx
import mmtbx.secondary_structure
import iotbx
import libtbx.load_env

pdb_str = """\
ATOM      1  N   ALA     1       1.558  -2.549  -1.143  1.00  0.00           N
ATOM      2  CA  ALA     1       1.130  -3.826  -1.702  1.00  0.00           C
ATOM      3  C   ALA     1      -0.241  -3.705  -2.358  1.00  0.00           C
ATOM      4  O   ALA     1      -0.532  -4.386  -3.342  1.00  0.00           O
ATOM      5  CB  ALA     1       1.108  -4.897  -0.622  1.00  0.00           C
ATOM      6  N   ALA     2      -1.081  -2.834  -1.807  1.00  0.00           N
ATOM      7  CA  ALA     2      -2.422  -2.623  -2.338  1.00  0.00           C
ATOM      8  C   ALA     2      -2.392  -1.711  -3.560  1.00  0.00           C
ATOM      9  O   ALA     2      -3.170  -1.887  -4.497  1.00  0.00           O
ATOM     10  CB  ALA     2      -3.331  -2.042  -1.265  1.00  0.00           C
ATOM     11  N   ALA     3      -1.488  -0.736  -3.541  1.00  0.00           N
ATOM     12  CA  ALA     3      -1.351   0.207  -4.645  1.00  0.00           C
ATOM     13  C   ALA     3      -0.825  -0.487  -5.898  1.00  0.00           C
ATOM     14  O   ALA     3      -1.270  -0.201  -7.010  1.00  0.00           O
ATOM     15  CB  ALA     3      -0.437   1.356  -4.251  1.00  0.00           C
ATOM     16  N   ALA     4       0.124  -1.398  -5.708  1.00  0.00           N
ATOM     17  CA  ALA     4       0.699  -2.149  -6.817  1.00  0.00           C
ATOM     18  C   ALA     4      -0.333  -3.092  -7.426  1.00  0.00           C
ATOM     19  O   ALA     4      -0.387  -3.268  -8.643  1.00  0.00           O
ATOM     20  CB  ALA     4       1.922  -2.926  -6.355  1.00  0.00           C
ATOM     21  N   ALA     5      -1.152  -3.695  -6.569  1.00  0.00           N
ATOM     22  CA  ALA     5      -2.204  -4.597  -7.020  1.00  0.00           C
ATOM     23  C   ALA     5      -3.295  -3.830  -7.757  1.00  0.00           C
ATOM     24  O   ALA     5      -3.850  -4.313  -8.744  1.00  0.00           O
ATOM     25  CB  ALA     5      -2.792  -5.358  -5.842  1.00  0.00           C
ATOM     26  N   ALA     6      -3.597  -2.630  -7.271  1.00  0.00           N
ATOM     27  CA  ALA     6      -4.599  -1.778  -7.900  1.00  0.00           C
ATOM     28  C   ALA     6      -4.097  -1.252  -9.240  1.00  0.00           C
ATOM     29  O   ALA     6      -4.879  -1.036 -10.166  1.00  0.00           O
ATOM     30  CB  ALA     6      -4.971  -0.626  -6.981  1.00  0.00           C
ATOM     31  N   ALA     7      -2.788  -1.047  -9.335  1.00  0.00           N
ATOM     32  CA  ALA     7      -2.176  -0.571 -10.570  1.00  0.00           C
ATOM     33  C   ALA     7      -2.045  -1.703 -11.583  1.00  0.00           C
ATOM     34  O   ALA     7      -2.052  -1.472 -12.792  1.00  0.00           O
ATOM     35  CB  ALA     7      -0.816   0.047 -10.286  1.00  0.00           C
ATOM     36  N   ALA     8      -1.927  -2.929 -11.081  1.00  0.00           N
ATOM     37  CA  ALA     8      -1.805  -4.100 -11.941  1.00  0.00           C
ATOM     38  C   ALA     8      -3.169  -4.549 -12.452  1.00  0.00           C
ATOM     39  O   ALA     8      -3.290  -5.046 -13.572  1.00  0.00           O
ATOM     40  CB  ALA     8      -1.118  -5.234 -11.197  1.00  0.00           C
ATOM     41  N   ALA     9      -4.195  -4.370 -11.625  1.00  0.00           N
ATOM     42  CA  ALA     9      -5.552  -4.755 -11.995  1.00  0.00           C
ATOM     43  C   ALA     9      -6.141  -3.782 -13.010  1.00  0.00           C
ATOM     44  O   ALA     9      -6.957  -4.163 -13.849  1.00  0.00           O
ATOM     45  CB  ALA     9      -6.436  -4.832 -10.760  1.00  0.00           C
ATOM     46  N   ALA    10      -5.722  -2.523 -12.929  1.00  0.00           N
ATOM     47  CA  ALA    10      -6.206  -1.494 -13.841  1.00  0.00           C
ATOM     48  C   ALA    10      -5.523  -1.598 -15.200  1.00  0.00           C
ATOM     49  O   ALA    10      -6.184  -1.624 -16.238  1.00  0.00           O
ATOM     50  CB  ALA    10      -5.991  -0.112 -13.243  1.00  0.00           C
TER
"""

def exercise_geo_out():
  log = null_out()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=pdb_str.split('\n'),
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=log)
  acp = processed_pdb_file.all_chain_proxies
  source = acp.pdb_inp.source_info()

  defpars = mmtbx.secondary_structure.sec_str_master_phil.fetch()
  custom_pars = defpars.fetch(source = iotbx.phil.parse(
      "h_bond_restraints.remove_outliers=False\n"))
  ss_manager = mmtbx.secondary_structure.manager(
      pdb_hierarchy=acp.pdb_hierarchy, params=custom_pars.extract())
  ss_manager.find_automatically(log=log)
  ss_manager.initialize(log=log)
  proxies_for_grm = ss_manager.create_hbond_proxies(
      log          = log,
      as_python_objects = False)
  custom_nb_excl = proxies_for_grm.exclude_nb_list
  grm = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      hydrogen_bond_proxies=proxies_for_grm.proxies,
      custom_nonbonded_exclusions = custom_nb_excl,
      assume_hydrogens_all_missing = True)
  from mmtbx.geometry_restraints import ramachandran
  params = ramachandran.master_phil.fetch().extract()
  params.rama_potential = "emsley"
  #params.rama_weight = sigma_on_ramachandran
  proxies = ramachandran.extract_proxies(acp.pdb_hierarchy, log=log)
  rama_lookup = ramachandran.lookup_manager(params)
  restraints_helper = mmtbx.geometry_restraints.manager(
      ramachandran_proxies=proxies,
      ramachandran_lookup=rama_lookup,
      hydrogen_bond_proxies=proxies_for_grm.proxies,
      hydrogen_bond_params=None)
  grm.set_generic_restraints(restraints_helper)

  atoms = acp.pdb_hierarchy.atoms()
  assert grm.generic_restraints_manager.get_n_hbonds() == 4
  import StringIO
  str_result = StringIO.StringIO()
  grm.generic_restraints_manager.show_sorted_hbonds(
      by_value="residual",
      sites_cart = atoms.extract_xyz(),
      site_labels= [a.id_str() for a in atoms],
      f=str_result)
  assert str_result.getvalue() == """\
Hydrogen bond restraints: 4
Sorted by residual:
hbond pdb=" O   ALA     4 "
      pdb=" N   ALA     8 "
  ideal  model  delta     sigma   weight residual
  2.900  2.904 -0.004  5.00e-02 4.00e+02 4.93e-03
hbond pdb=" O   ALA     2 "
      pdb=" N   ALA     6 "
  ideal  model  delta     sigma   weight residual
  2.900  2.903 -0.003  5.00e-02 4.00e+02 4.49e-03
hbond pdb=" O   ALA     3 "
      pdb=" N   ALA     7 "
  ideal  model  delta     sigma   weight residual
  2.900  2.903 -0.003  5.00e-02 4.00e+02 2.92e-03
hbond pdb=" O   ALA     5 "
      pdb=" N   ALA     9 "
  ideal  model  delta     sigma   weight residual
  2.900  2.902 -0.002  5.00e-02 4.00e+02 1.84e-03

"""
  str_result = StringIO.StringIO()
  grm.generic_restraints_manager.show_sorted_ramachandran(
      by_value="residual",
      sites_cart = atoms.extract_xyz(),
      site_labels= [a.id_str() for a in atoms],
      f=str_result)
  assert str_result.getvalue() == """\
Ramachandran plot restraints: 8
Sorted by residual:
phi-psi angles formed by             residual
    pdb=" C   ALA     7 "            1.53e+01
    pdb=" N   ALA     8 "
    pdb=" CA  ALA     8 "
    pdb=" C   ALA     8 "
    pdb=" N   ALA     9 "
phi-psi angles formed by             residual
    pdb=" C   ALA     1 "            1.52e+01
    pdb=" N   ALA     2 "
    pdb=" CA  ALA     2 "
    pdb=" C   ALA     2 "
    pdb=" N   ALA     3 "
phi-psi angles formed by             residual
    pdb=" C   ALA     6 "            1.20e+01
    pdb=" N   ALA     7 "
    pdb=" CA  ALA     7 "
    pdb=" C   ALA     7 "
    pdb=" N   ALA     8 "
phi-psi angles formed by             residual
    pdb=" C   ALA     2 "            1.14e+01
    pdb=" N   ALA     3 "
    pdb=" CA  ALA     3 "
    pdb=" C   ALA     3 "
    pdb=" N   ALA     4 "
phi-psi angles formed by             residual
    pdb=" C   ALA     8 "            1.06e+01
    pdb=" N   ALA     9 "
    pdb=" CA  ALA     9 "
    pdb=" C   ALA     9 "
    pdb=" N   ALA    10 "
phi-psi angles formed by             residual
    pdb=" C   ALA     4 "            1.06e+01
    pdb=" N   ALA     5 "
    pdb=" CA  ALA     5 "
    pdb=" C   ALA     5 "
    pdb=" N   ALA     6 "
phi-psi angles formed by             residual
    pdb=" C   ALA     3 "            1.03e+01
    pdb=" N   ALA     4 "
    pdb=" CA  ALA     4 "
    pdb=" C   ALA     4 "
    pdb=" N   ALA     5 "
phi-psi angles formed by             residual
    pdb=" C   ALA     5 "            7.58e+00
    pdb=" N   ALA     6 "
    pdb=" CA  ALA     6 "
    pdb=" C   ALA     6 "
    pdb=" N   ALA     7 "

"""

if (__name__ == "__main__") :
  if libtbx.env.has_module("ksdssp") :
    exercise_geo_out()
    print "OK"
  else :
    print "skipping test, ksdssp not available"
