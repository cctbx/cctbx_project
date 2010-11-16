from mmtbx.hydrogens import find as find_hydrogens
import mmtbx.utils
import mmtbx.model
import mmtbx.restraints
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import iotbx.pdb
from cctbx import miller
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import format_cpu_times
from cStringIO import StringIO

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
HETATM    2 H1   HOH X   1      19.339  10.959   7.882  1.00  1.00           H
HETATM    3 H2   HOH X   1      19.339  10.959   7.882  1.00  1.00           H
HETATM    4  H1  HOH X   2      20.024  13.469  10.445  1.00  2.00           H
HETATM    5  O   HOH X   2      19.563  13.548   9.584  1.00  1.00           O
HETATM    6 H2   HOH X   2      18.608  13.419   9.761  1.00  1.00           H
HETATM    7  O   HOH     1      17.253  15.650  10.892  1.00  2.00           O
HETATM    8  H1  HOH     1      18.146  15.616  11.293  1.00  2.00           H
HETATM    9  H2  HOH     1      17.357  15.408   9.948  1.00  2.00           H
HETATM   10  O   HOH     2       5.150  11.586  12.474  1.00  3.00           O
HETATM   11 H1   HOH     2       4.364  12.000  12.888  1.00  3.00           H
HETATM   12 H2   HOH     2       5.938  11.998  12.886  1.00  3.00           H
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
ATOM     26  D1  HOH Q   1      14.941  12.899   7.563  1.00  4.00           D
ATOM     27 D2   HOH Q   1      16.514  12.892   7.541  1.00  4.00           D
ATOM     28  O   HOH S   1      14.337  16.126   9.974  1.00  5.00           O
ATOM     29  H2  HOH S   1      14.968  16.680   9.469  1.00  5.00           H
ATOM     30 H1   HOH S   1      14.734  15.233  10.041  1.00  5.00           H
ATOM     31  O   HOH S   2       9.491  12.823  15.494  1.00  6.00           O
ATOM     32 H1   HOH S   2       8.899  13.604  15.494  1.00  6.00           H
ATOM     33 H2   HOH S   2       8.914  12.031  15.494  1.00  6.00           H
END
"""

def exercise_01():
  pdb_file_name = "add_h_to_hoh.pdb"
  tmp_f = open(pdb_file_name, "w")
  tmp_f.write(input_model)
  tmp_f.close()
  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(log = StringIO())
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file_name])
  xray_structure = processed_pdb_file.xray_structure()
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
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log                     = None)
  ####
  model.add_hydrogens()
  result = StringIO()
  model.write_pdb_file(out=result)
  result = result.getvalue().splitlines()
  ####
  result1 = []
  for r1 in result:
    if(r1.startswith("ATOM") or r1.startswith("HETATM")): result1.append(r1)
  result2 = []
  for r2 in expected_result.splitlines():
    if(r2.startswith("ATOM") or r2.startswith("HETATM")): result2.append(r2)
  assert len(result1) == len(result2)
  for r1, r2 in zip(result1, result2):
    r1 = r1[:30] + r1[60:]
    r2 = r2[:30] + r2[60:]
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

model_good = """\
CRYST1   15.000   15.000   15.000  80.00  70.00 100.00 P 1
ATOM      1  CB  PHE A   1      12.073   9.948  11.650  1.00  5.00           C
ATOM      2  CG  PHE A   1      11.316   9.114  12.643  1.00  5.00           C
ATOM      3  CD1 PHE A   1      10.646   7.971  12.240  1.00  5.00           C
ATOM      4  CD2 PHE A   1      11.276   9.472  13.980  1.00  5.00           C
ATOM      5  CE1 PHE A   1       9.949   7.202  13.152  1.00  5.00           C
ATOM      6  CE2 PHE A   1      10.581   8.706  14.897  1.00  5.00           C
ATOM      7  CZ  PHE A   1       9.916   7.570  14.482  1.00  5.00           C
ATOM      8  C   PHE A   1      12.110  11.913  10.097  1.00  5.00           C
ATOM      9  O   PHE A   1      12.051  11.781   8.875  1.00  5.00           O
ATOM     10  OXT PHE A   1      12.892  12.748  10.551  1.00  5.00           O
ATOM     11  N   PHE A   1      10.096  10.513  10.307  1.00  5.00           N
ATOM     12  CA  PHE A   1      11.240  11.067  11.021  1.00  5.00           C
HETATM   13  O   HOH     1      13.866  16.009  12.098  1.00  3.00           O
HETATM   14  H1  HOH     1      13.327  16.140  12.905  1.00  3.00           H
HETATM   15  H2  HOH     1      14.117  16.901  11.777  1.00  3.00           H
HETATM   16  O   HOH     2      17.215  16.288  11.122  1.00  3.00           O
HETATM   17  H1  HOH     2      17.912  15.900  10.553  1.00  3.00           H
HETATM   18  H2  HOH     2      17.641  16.523  11.973  1.00  3.00           H
HETATM   19  O   HOH     3       8.927  12.312  13.459  1.00  3.00           O
HETATM   20  H1  HOH     3       8.879  11.933  14.362  1.00  3.00           H
HETATM   21  H2  HOH     3       8.721  11.584  12.835  1.00  3.00           H
HETATM   22  O   HOH     4      16.005  11.974   8.964  1.00  3.00           O
HETATM   23  H1  HOH     4      16.817  12.516   9.046  1.00  3.00           H
HETATM   24  H2  HOH     4      15.632  12.153   8.076  0.00  3.00           H
HETATM   25  O   HOH     5      13.626   9.207   8.521  1.00  3.00           O
HETATM   26  H1  HOH     5      13.711   9.868   7.803  1.00  3.00           H
HETATM   27  H2  HOH     5      14.464   9.228   9.028  0.00  3.00           H
HETATM   28  O   HOH     6       9.841  14.509  11.210  1.00  3.00           O
HETATM   29  H1  HOH     6      10.096  14.540  12.156  1.00  3.00           H
HETATM   30  H2  HOH     6       9.341  13.676  11.079  0.00  3.00           H
END
"""
model_bad = """\
CRYST1   15.000   15.000   15.000  80.00  70.00 100.00 P 1
ATOM      1  CB  PHE A   1      12.073   9.948  11.650  1.00  5.00           C
ATOM      2  CG  PHE A   1      11.316   9.114  12.643  1.00  5.00           C
ATOM      3  CD1 PHE A   1      10.646   7.971  12.240  1.00  5.00           C
ATOM      4  CD2 PHE A   1      11.276   9.472  13.980  1.00  5.00           C
ATOM      5  CE1 PHE A   1       9.949   7.202  13.152  1.00  5.00           C
ATOM      6  CE2 PHE A   1      10.581   8.706  14.897  1.00  5.00           C
ATOM      7  CZ  PHE A   1       9.916   7.570  14.482  1.00  5.00           C
ATOM      8  C   PHE A   1      12.110  11.913  10.097  1.00  5.00           C
ATOM      9  O   PHE A   1      12.051  11.781   8.875  1.00  5.00           O
ATOM     10  OXT PHE A   1      12.892  12.748  10.551  1.00  5.00           O
ATOM     11  N   PHE A   1      10.096  10.513  10.307  1.00  5.00           N
ATOM     12  CA  PHE A   1      11.240  11.067  11.021  1.00  5.00           C
HETATM   13  O   HOH     1      13.866  16.009  12.098  1.00  3.00           O
HETATM   14  H1  HOH     1      13.114  16.628  12.001  0.00  3.00           H
HETATM   15  H2  HOH     1      13.771  15.335  11.392  0.00  3.00           H
HETATM   16  O   HOH     2      17.215  16.288  11.122  1.00  3.00           O
HETATM   17  H1  HOH     2      16.474  16.378  10.487  0.00  3.00           H
HETATM   18  H2  HOH     2      17.045  15.476  11.644  0.00  3.00           H
HETATM   19  O   HOH     3       8.927  12.312  13.459  1.00  3.00           O
HETATM   20  H1  HOH     3       8.895  11.898  14.347  0.00  3.00           H
HETATM   21  H2  HOH     3       7.999  12.412  13.159  0.00  3.00           H
HETATM   22  O   HOH     4      16.005  11.974   8.964  1.00  3.00           O
HETATM   23  H1  HOH     4      16.427  11.730   9.814  0.00  3.00           H
HETATM   24  H2  HOH     4      16.642  11.748   8.255  0.00  3.00           H
HETATM   25  O   HOH     5      13.626   9.207   8.521  1.00  3.00           O
HETATM   26  H1  HOH     5      13.711   9.868   7.803  0.00  3.00           H
HETATM   27  H2  HOH     5      14.464   9.228   9.028  0.00  3.00           H
HETATM   28  O   HOH     6       9.839  14.506  11.213  1.00  3.00           O
HETATM   29  H1  HOH     6      10.483  13.917  10.766  0.00  3.00           H
HETATM   30  H2  HOH     6       9.067  14.586  10.614  0.00  3.00           H
END
"""
expected_result2 = """\

==================== Fit water hydrogens into residual map ====================


                    ----------find peak-candidates----------

Number of peaks found at mFobs-DFmodel map (map cutoff=6.50 sigma)= 9
Filter by distance & map next to the model:
mapped sites are within: 0.980 - 1.009
number of sites selected in [dist_min= 0.70, dist_max= 1.05]: 9 from: 9
mapped sites are within: 0.980 - 1.009

peak=   23.135 closest distance to pdb=" O   HOH     1 " =    0.990
peak=   20.851 closest distance to pdb=" O   HOH     1 " =    0.981
peak=   22.681 closest distance to pdb=" O   HOH     2 " =    0.999
peak=   21.340 closest distance to pdb=" O   HOH     2 " =    0.997
peak=   22.116 closest distance to pdb=" O   HOH     3 " =    0.990
peak=   21.072 closest distance to pdb=" O   HOH     3 " =    0.980
peak=   21.269 closest distance to pdb=" O   HOH     4 " =    0.999
peak=   22.049 closest distance to pdb=" O   HOH     5 " =    0.988
peak=   20.580 closest distance to pdb=" O   HOH     6 " =    1.009

----------6D rigid body fit of HOH----------

Fit quality:
0.028
0.041
0.022
0.019
0.008
0.028
"""

def exercise_02():
  for file_name, input_model in [("m_good.pdb",model_good), ("m_bad.pdb",model_bad)]:
    tmp_f = open(file_name, "w")
    tmp_f.write(input_model)
    tmp_f.close()
  xrs_exact = iotbx.pdb.input(file_name = "m_good.pdb").xray_structure_simple()
  xrs_part = iotbx.pdb.input(file_name = "m_bad.pdb").xray_structure_simple()
  miller_set = miller.build_set(
    crystal_symmetry = xrs_exact.crystal_symmetry(),
    anomalous_flag   = False,
    d_min            = 0.6)
  f_obs = abs(miller_set.structure_factors_from_scatterers(
    xray_structure = xrs_exact,
    algorithm      = "direct",
    cos_sin_table  = False).f_calc())
  sf_par = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sf_par.algorithm = "direct"
  sf_par.cos_sin_table = False
  flags = f_obs.array(data=flex.bool(f_obs.data().size(),False))
  fmodel = mmtbx.f_model.manager(
    xray_structure               = xrs_part,
    sf_and_grads_accuracy_params = sf_par,
    r_free_flags                 = flags,
    target_name                  = "ls_wunit_k1",
    f_obs                        = f_obs)
  #
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = ener_lib,
    file_name                = "m_bad.pdb")
  model = mmtbx.model.manager(
    xray_structure             = xrs_part,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy,
    log                        = None)
  #
  out = StringIO()
  params = find_hydrogens.all_master_params().extract()
  params.map_cutoff=6.5
  find_hydrogens.run(fmodel=fmodel, model=model, log=out, params=params)
  for a,b in zip(out.getvalue().splitlines(), expected_result2.splitlines()):
    a = a.strip()
    b = b.strip()
    try:
      a = float(a)
      b = float(b)
      assert a < 0.05
      assert b < 0.05
    except:
      assert not show_diff(a, b)

def exercise():
  exercise_01()
  exercise_02()
  print format_cpu_times()

if (__name__ == "__main__"):
  exercise()
