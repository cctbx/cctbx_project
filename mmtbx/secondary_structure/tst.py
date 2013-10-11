from __future__ import division

from mmtbx.secondary_structure import manager
from mmtbx.secondary_structure.base_pairing import pair_database
from mmtbx.geometry_restraints import hbond
from iotbx import file_reader
import iotbx.pdb.hierarchy
import libtbx.load_env
from libtbx.utils import null_out
import os

def exercise_protein () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_file_h = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print "Skipping exercise(): input file not available."
    return False
  if pdb_file_h is None :
    print "Skipping exercise(): input file not available."
    return False
  run_ksdssp = True
  if (not libtbx.env.has_module(name="ksdssp")):
    print "Skipping KSDSSP tests: ksdssp module not available."
    run_ksdssp = False
  log = null_out()
  expected_distances = [2.9, 1.975]
  for k, file_name in enumerate([pdb_file, pdb_file_h]) :
    pdb_in = file_reader.any_file(file_name, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    m.find_automatically(log=log)
    m.params.h_bond_restraints.remove_outliers = False
    hbond_params = hbond.master_phil.extract()
    build_proxies = m.create_hbond_proxies(hbond_params=hbond_params,
      log=log)
    proxies = build_proxies.proxies
    assert (len(proxies) == len(build_proxies.exclude_nb_list) == 109)
    assert (type(proxies[0]).__name__ == "h_bond_simple_proxy")
    assert (proxies[0].distance_ideal == expected_distances[k])
    (frac_alpha, frac_beta) = m.calculate_structure_content()
    assert ("%.3f" % frac_alpha == "0.643")
    assert ("%.3f" % frac_beta == "0.075")
    # Make sure the hydrogen auto-detection override is working
    m.params.h_bond_restraints.substitute_n_for_h = True
    build_proxies = m.create_hbond_proxies(
      log=log,
      as_python_objects=True)
    proxies = build_proxies.proxies
    atom_ids = []
    for i_seq in proxies[0].i_seqs :
      atom_ids.append(pdb_hierarchy.atoms()[i_seq].id_str())
    assert (atom_ids == ['pdb=" N   ARG A  41 "', 'pdb=" O   ASP A  37 "'])
    assert (proxies[0].distance_ideal == 2.9)
    if (run_ksdssp) :
      m = manager(pdb_hierarchy=pdb_hierarchy,
        sec_str_from_pdb_file=None)
      m.find_automatically(log=log)
      m.params.h_bond_restraints.remove_outliers = False
      hbond_params.restraint_type = "simple"
      build_proxies = m.create_hbond_proxies(hbond_params=hbond_params, log=log)
      assert (build_proxies.proxies.size() == 90)

def exercise_nucleic_acids () :
  pdb_file_rna = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1u8d.pdb",
    test=os.path.isfile)
  pdb_file_dna = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/dnatest.pdb",
    test=os.path.isfile)
  if (pdb_file_rna is None) :
    print "Skipping exercise_nucleic_acids(): input file not available."
    return False
  # Nucleic acids (requires REDUCE and PROBE)
  if (libtbx.env.has_module(name="reduce") and
      libtbx.env.has_module(name="probe")):
    pdb_in = file_reader.any_file(pdb_file_rna, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    log = null_out()
    m.find_automatically(log=log)
    build_proxies = m.create_hbond_proxies(log=log)
    assert (build_proxies.proxies.size() == 70)
    db = pair_database()
    assert (db.get_atoms("A-U", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
    assert (db.get_atoms("DA-DT", "WWT", False) == [('N6', 'O4'), ('N1', 'N3')])
    assert (db.get_atoms("DT-DA", "WWT", False) == [('O4', 'N6'), ('N3', 'N1')])
    assert (db.get_atoms("G-C", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
      ('N2', 'O2')])
    assert (db.get_atoms("C-G", "WWT", False) == [('N4', 'O6'), ('N3', 'N1'),
      ('O2', 'N2')])
    assert db.get_pair_type("A-U", [('H61', 'O4'), ('N1', 'H3')], True) == "XX"
    assert db.get_pair_type("A-U", [('N1', 'H3'), ('H61', 'O4')], True) == "XX"
    assert db.get_pair_type("C-G", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "XIX"
    # DNA
    pdb_in = file_reader.any_file(pdb_file_dna, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    log = null_out()
    m.find_automatically(log=log)
    assert (len(m.params.nucleic_acids.base_pair) == 4)
  else:
    print "Skipping base-pairing tests: reduce or probe module not available."
    pass
  print "OK"

def exercise_sheet_ends () :
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string="""\
ATOM      1  N   ALA A  37      27.738  -0.961  14.033  1.00 20.00           N
ATOM      2  CA  ALA A  37      27.940   0.419  14.337  1.00 20.00           C
ATOM      3  C   ALA A  37      26.683   0.986  14.948  1.00 20.00           C
ATOM      4  O   ALA A  37      26.735   1.984  15.528  1.00 20.00           O
ATOM      5  CB  ALA A  37      28.398   1.219  13.140  1.00 20.00           C
ATOM      6  N   ALA A  38      25.588   0.264  14.880  1.00 20.00           N
ATOM      7  CA  ALA A  38      24.384   0.731  15.506  1.00 20.00           C
ATOM      8  C   ALA A  38      23.201   0.011  14.926  1.00 20.00           C
ATOM      9  O   ALA A  38      23.327  -1.063  14.341  1.00 20.00           O
ATOM     10  N   ALA A  39      22.037   0.615  15.074  1.00 20.00           N
ATOM     11  CA  ALA A  39      20.822  -0.033  14.646  1.00 20.00           C
ATOM     12  C   ALA A  39      19.950   0.873  13.809  1.00 20.00           C
ATOM     13  O   ALA A  39      20.035   2.092  13.903  1.00 20.00           O
ATOM     14  CB  ALA A  39      20.054  -0.530  15.864  1.00 20.00           C
ATOM     15  N   ALA A  40      19.119   0.270  12.969  1.00 20.00           N
ATOM     16  CA  ALA A  40      18.164   1.039  12.204  1.00 20.00           C
ATOM     17  C   ALA A  40      16.761   0.571  12.489  1.00 20.00           C
ATOM     18  O   ALA A  40      16.444  -0.583  12.263  1.00 20.00           O
ATOM     19  CB  ALA A  40      18.436   0.930  10.717  1.00 20.00           C
ATOM     20  N   ALA A  41      15.938   1.493  12.973  1.00 20.00           N
ATOM     21  CA  ALA A  41      14.533   1.270  13.264  1.00 20.00           C
ATOM     22  C   ALA A  41      13.658   1.820  12.148  1.00 20.00           C
ATOM     23  O   ALA A  41      13.680   3.008  11.883  1.00 20.00           O
ATOM     24  CB  ALA A  41      14.135   1.952  14.582  1.00 20.00           C
ATOM     25  N   ALA A  42      12.897   0.961  11.489  1.00 20.00           N
ATOM     26  CA  ALA A  42      11.938   1.402  10.495  1.00 20.00           C
ATOM     27  C   ALA A  42      10.559   1.417  11.119  1.00 20.00           C
ATOM     28  O   ALA A  42      10.172   0.468  11.765  1.00 20.00           O
ATOM     29  N   ALA A  43       9.825   2.465  10.888  1.00 20.00           N
ATOM     30  CA  ALA A  43       8.546   2.641  11.492  1.00 20.00           C
ATOM     31  C   ALA A  43       7.468   3.106  10.585  1.00 20.00           C
ATOM     32  O   ALA A  43       7.710   3.854   9.738  1.00 20.00           O
ATOM     33  CB  ALA A  43       8.646   3.437  12.767  1.00 20.00           C
ATOM     34  N   ALA A  44       6.285   2.569  10.752  1.00 20.00           N
ATOM     35  CA  ALA A  44       5.209   2.771   9.823  1.00 20.00           C
ATOM     36  C   ALA A  44       3.822   2.663  10.467  1.00 20.00           C
ATOM     37  O   ALA A  44       3.610   1.972  11.406  1.00 20.00           O
ATOM     38  CB  ALA A  44       5.381   1.912   8.510  1.00 20.00           C
ATOM     39  N   ALA A  45       2.883   3.390   9.914  1.00 20.00           N
ATOM     40  CA  ALA A  45       1.539   3.355  10.375  1.00 20.00           C
ATOM     41  C   ALA A  45       0.845   2.276   9.614  1.00 20.00           C
ATOM     42  O   ALA A  45       1.284   1.970   8.558  1.00 20.00           O
ATOM     43  CB  ALA A  45       0.868   4.648  10.084  1.00 20.00           C
ATOM     44  N   ALA A  74      -2.478   6.652  17.490  1.00 20.00           N
ATOM     45  CA  ALA A  74      -2.170   5.880  16.289  1.00 20.00           C
ATOM     46  C   ALA A  74      -1.360   4.598  16.476  1.00 20.00           C
ATOM     47  O   ALA A  74      -0.341   4.581  17.153  1.00 20.00           O
ATOM     48  CB  ALA A  74      -1.405   6.748  15.290  1.00 20.00           C
ATOM     49  N   ALA A  75      -1.809   3.528  15.830  1.00 20.00           N
ATOM     50  CA  ALA A  75      -1.038   2.294  15.762  1.00 20.00           C
ATOM     51  C   ALA A  75       0.178   2.394  14.825  1.00 20.00           C
ATOM     52  O   ALA A  75       0.029   2.691  13.638  1.00 20.00           O
ATOM     53  CB  ALA A  75      -1.942   1.143  15.310  1.00 20.00           C
ATOM     54  N   ALA A  76       1.369   2.129  15.363  1.00 20.00           N
ATOM     55  CA  ALA A  76       2.612   2.062  14.578  1.00 20.00           C
ATOM     56  C   ALA A  76       3.222   0.655  14.618  1.00 20.00           C
ATOM     57  O   ALA A  76       3.065  -0.086  15.584  1.00 20.00           O
ATOM     58  CB  ALA A  76       3.655   3.064  15.089  1.00 20.00           C
ATOM     59  N   ALA A  77       3.933   0.308  13.558  1.00 20.00           N
ATOM     60  CA  ALA A  77       4.639  -0.950  13.487  1.00 20.00           C
ATOM     61  C   ALA A  77       6.146  -0.691  13.381  1.00 20.00           C
ATOM     62  O   ALA A  77       6.547   0.319  12.808  1.00 20.00           O
ATOM     63  CB  ALA A  77       4.125  -1.752  12.300  1.00 20.00           C
ATOM     64  N   ALA A  78       6.978  -1.586  13.921  1.00 20.00           N
ATOM     65  CA  ALA A  78       8.426  -1.349  13.964  1.00 20.00           C
ATOM     66  C   ALA A  78       9.300  -2.568  13.613  1.00 20.00           C
ATOM     67  O   ALA A  78       8.926  -3.706  13.886  1.00 20.00           O
ATOM     68  CB  ALA A  78       8.844  -0.866  15.364  1.00 20.00           C
ATOM     69  N   ALA A  79      10.475  -2.320  13.029  1.00 20.00           N
ATOM     70  CA  ALA A  79      11.544  -3.327  12.991  1.00 20.00           C
ATOM     71  C   ALA A  79      12.910  -2.716  13.258  1.00 20.00           C
ATOM     72  O   ALA A  79      13.172  -1.578  12.910  1.00 20.00           O
ATOM     73  CB  ALA A  79      11.648  -4.067  11.649  1.00 20.00           C
ATOM     74  N   ALA A  80      13.783  -3.507  13.858  1.00 20.00           N
ATOM     75  CA  ALA A  80      15.135  -3.088  14.163  1.00 20.00           C
ATOM     76  C   ALA A  80      16.153  -4.001  13.478  1.00 20.00           C
ATOM     77  O   ALA A  80      16.078  -5.224  13.575  1.00 20.00           O
ATOM     78  CB  ALA A  80      15.395  -3.117  15.671  1.00 20.00           C
ATOM     79  N   ALA A  81      17.119  -3.408  12.798  1.00 20.00           N
ATOM     80  CA  ALA A  81      18.232  -4.178  12.289  1.00 20.00           C
ATOM     81  C   ALA A  81      19.542  -3.491  12.560  1.00 20.00           C
ATOM     82  O   ALA A  81      19.643  -2.276  12.472  1.00 20.00           O
ATOM     83  CB  ALA A  81      18.067  -4.425  10.802  1.00 20.00           C
ATOM     84  N   ALA A  82      20.541  -4.285  12.908  1.00 20.00           N
ATOM     85  CA  ALA A  82      21.892  -3.793  13.044  1.00 20.00           C
ATOM     86  C   ALA A  82      22.410  -3.375  11.697  1.00 20.00           C
ATOM     87  O   ALA A  82      21.897  -3.798  10.685  1.00 20.00           O
ATOM     88  CB  ALA A  82      22.786  -4.846  13.647  1.00 20.00           C
ATOM     89  N   ALA A  83      23.426  -2.532  11.679  1.00 20.00           N
ATOM     90  CA  ALA A  83      24.083  -2.208  10.429  1.00 20.00           C
ATOM     91  C   ALA A  83      25.543  -2.018  10.735  1.00 20.00           C
ATOM     92  O   ALA A  83      25.905  -1.692  11.853  1.00 20.00           O
ATOM     93  CB  ALA A  83      23.488  -0.951   9.765  1.00 20.00           C
ATOM     94  N   ALA A  84      26.388  -2.269   9.753  1.00 20.00           N
ATOM     95  CA  ALA A  84      27.766  -1.859   9.852  1.00 20.00           C
ATOM     96  C   ALA A  84      28.199  -1.539   8.456  1.00 20.00           C
ATOM     97  O   ALA A  84      27.350  -1.415   7.580  1.00 20.00           O
ATOM     98  CB  ALA A  84      28.646  -2.929  10.507  1.00 20.00           C
ATOM     99  N   ALA A  94       8.539  -8.036  13.764  1.00 20.00           N
ATOM    100  CA  ALA A  94       7.784  -6.815  13.849  1.00 20.00           C
ATOM    101  C   ALA A  94       7.141  -6.661  15.220  1.00 20.00           C
ATOM    102  O   ALA A  94       6.628  -7.611  15.776  1.00 20.00           O
ATOM    103  CB  ALA A  94       6.745  -6.817  12.744  1.00 20.00           C
ATOM    104  N   ALA A  95       7.188  -5.460  15.782  1.00 20.00           N
ATOM    105  CA  ALA A  95       6.447  -5.181  17.013  1.00 20.00           C
ATOM    106  C   ALA A  95       5.605  -3.906  16.830  1.00 20.00           C
ATOM    107  O   ALA A  95       5.855  -3.110  15.933  1.00 20.00           O
ATOM    108  CB  ALA A  95       7.387  -5.044  18.220  1.00 20.00           C
ATOM    109  N   ALA A  96       4.581  -3.743  17.658  1.00 20.00           N
ATOM    110  CA  ALA A  96       3.609  -2.664  17.488  1.00 20.00           C
ATOM    111  C   ALA A  96       3.285  -1.911  18.775  1.00 20.00           C
ATOM    112  O   ALA A  96       3.129  -2.501  19.839  1.00 20.00           O
ATOM    113  CB  ALA A  96       2.316  -3.219  16.906  1.00 20.00           C
ATOM    114  N   ALA A  97       3.188  -0.599  18.670  1.00 20.00           N
ATOM    115  CA  ALA A  97       2.653   0.184  19.763  1.00 20.00           C
ATOM    116  C   ALA A  97       1.705   1.289  19.259  1.00 20.00           C
ATOM    117  O   ALA A  97       1.820   1.731  18.115  1.00 20.00           O
ATOM    118  CB  ALA A  97       3.779   0.765  20.569  1.00 20.00           C
ATOM    119  N   ALA A  98       0.739   1.686  20.047  1.00 20.00           N
ATOM    120  CA  ALA A  98      -0.189   2.723  19.701  1.00 20.00           C
ATOM    121  C   ALA A  98       0.076   3.866  20.555  1.00 20.00           C
ATOM    122  O   ALA A  98       0.025   3.751  21.692  1.00 20.00           O
ATOM    123  CB  ALA A  98      -1.627   2.343  20.021  1.00 20.00           C
""")
  pdb_hierarchy = pdb_in.hierarchy
  pdb_hierarchy.atoms().reset_i_seq()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    sec_str_from_pdb_file=None)
  m.params.input.use_ksdssp = False
  m.find_automatically(log=null_out())
  m.params.h_bond_restraints.remove_outliers = False
  hbond_params = hbond.master_phil.extract()
  build_proxies = m.create_hbond_proxies(hbond_params=hbond_params,
    log=null_out())
  proxies = build_proxies.proxies
  assert len(proxies) == 12

if __name__ == "__main__" :
  exercise_protein()
  exercise_sheet_ends()
  exercise_nucleic_acids()
