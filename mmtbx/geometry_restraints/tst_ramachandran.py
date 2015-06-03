from __future__ import division
from scitbx.array_family import flex
from mmtbx.geometry_restraints import ramachandran
import mmtbx.geometry_restraints
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.monomer_library import server, pdb_interpretation
from mmtbx.refinement import geometry_minimization
import iotbx.pdb
import cctbx.geometry_restraints
import scitbx.lbfgs
from scitbx.array_family import flex
import boost.python
from libtbx.test_utils import approx_equal
import libtbx.load_env
from libtbx import group_args
from cStringIO import StringIO
import sys
import os

def exercise_basic () :
  t = ramachandran.load_tables()["ala"]
  assert approx_equal(t.get_score(0.0,0.0), -26.16, eps=0.01)
  assert approx_equal(t.get_score(-60,120), 10.41, eps=0.01)
  assert approx_equal(t.get_score(90,90), -7.43, eps=0.01)
  assert approx_equal(t.get_energy(0.0,0.0), 53.81, eps=0.01)
  assert approx_equal(t.get_energy(-60,120), 17.24, eps=0.01)
  e1 = t.get_energy(-85.0, 85.0)
  e2 = t.get_energy(-85.000001, 85.000001)
  assert approx_equal(e1, e2, eps=0.000001)
  assert approx_equal(t.get_energy(-85.0, 86.0), 21.3345, eps=0.001)
  assert approx_equal(t.get_energy(-86.0, 85.0), 21.389, eps=0.001)
  ext = boost.python.import_ext("mmtbx_ramachandran_restraints_ext")
  proxies = ext.shared_phi_psi_proxy()
  proxies.append(
    ext.phi_psi_proxy(
      i_seqs=[0,1,2,3,4],
      residue_name="ALA",
      residue_type="general"))
  proxies.append(
    ext.phi_psi_proxy(
      i_seqs=[4,5,6,7,8],
      residue_name="ALA",
      residue_type="prepro"))
  proxies.append(
    ext.phi_psi_proxy(
      i_seqs=[8,9,10,11,12],
      residue_name="ALA",
      residue_type="general"))
  selected = proxies.proxy_select(n_seq=13,
    iselection=flex.size_t(range(9)))
  assert (selected.size() == 2)
  assert list(selected[0].get_i_seqs()) == [0,1,2,3,4]
  assert selected[0].residue_name == "ALA"
  assert selected[0].residue_type == "general"
  assert list(selected[1].get_i_seqs()) == [4,5,6,7,8]
  assert selected[1].residue_name == "ALA"
  assert selected[1].residue_type == "prepro"

pdb1 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      2  CA  ALA A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      3  C   ALA A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      4  O   ALA A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      5  CB  ALA A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM      6  H   ALA A   2      -8.111   3.328   2.362  1.00 15.02           H
ATOM      7  N   ALA A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM      8  CA  ALA A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM      9  C   ALA A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     10  O   ALA A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     11  CB  ALA A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     12  H   ALA A   3      -4.629   0.611   3.830  1.00 12.26           H
ATOM     13  N   ALA A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     14  CA  ALA A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     15  C   ALA A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     16  O   ALA A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     17  CB  ALA A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     18  H   ALA A   4      -1.172   3.214   3.626  1.00 10.29           H
END
"""
pdb2 = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   ALA A   2      -7.133   3.038   1.936  1.00 15.02           N
ATOM      2  CA  ALA A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      3  C   ALA A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      4  O   ALA A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      5  CB  ALA A   2      -6.260   0.736   2.116  1.00 15.38           C
ATOM      6  H   ALA A   2      -7.242   2.739   0.988  1.00 15.02           H
ATOM      7  N   ALA A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM      8  CA  ALA A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM      9  C   ALA A   3      -2.919   3.406   4.698  1.00 11.10           C
ATOM     10  O   ALA A   3      -2.901   4.120   3.684  1.00 10.42           O
ATOM     11  CB  ALA A   3      -2.017   1.217   3.856  1.00 12.15           C
ATOM     12  H   ALA A   3      -4.629   0.611   3.830  1.00 12.26           H
ATOM     13  N   ALA A   4      -2.702   3.855   5.940  1.00 10.29           N
ATOM     14  CA  ALA A   4      -2.136   5.187   6.272  1.00 10.53           C
ATOM     15  C   ALA A   4      -2.252   5.506   7.787  1.00 10.24           C
ATOM     16  O   ALA A   4      -2.330   4.614   8.648  1.00  8.86           O
ATOM     17  CB  ALA A   4      -2.740   6.317   5.426  1.00  9.80           C
ATOM     18  H   ALA A   4      -2.913   3.309   6.751  1.00 10.29           H
END
"""
pdb3 = """
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   TRP A   2      -7.117   3.021   1.937  1.00  0.00           N
ATOM      2  CA  TRP A   2      -6.531   2.047   2.853  1.00  0.00           C
ATOM      3  C   TRP A   2      -5.234   2.561   3.433  1.00  0.00           C
ATOM      4  O   TRP A   2      -4.978   3.742   3.426  1.00  0.00           O
ATOM      5  CB  TRP A   2      -6.274   0.719   2.087  1.00  0.00           C
ATOM      6  CG  TRP A   2      -7.540  -0.019   1.642  1.00  0.00           C
ATOM      7  CD1 TRP A   2      -8.090   0.002   0.344  1.00  0.00           C
ATOM      8  CD2 TRP A   2      -8.418  -0.752   2.413  1.00  0.00           C
ATOM      9  CE2 TRP A   2      -9.485  -1.169   1.578  1.00  0.00           C
ATOM     10  CE3 TRP A   2      -8.409  -1.083   3.792  1.00  0.00           C
ATOM     11  NE1 TRP A   2      -9.302  -0.712   0.282  1.00  0.00           N
ATOM     12  CZ2 TRP A   2     -10.545  -1.931   2.116  1.00  0.00           C
ATOM     13  CZ3 TRP A   2      -9.464  -1.842   4.298  1.00  0.00           C
ATOM     14  CH2 TRP A   2     -10.516  -2.263   3.473  1.00  0.00           C
ATOM     25  N   TRP A   3      -4.433   1.613   3.906  1.00  0.00           N
ATOM     26  CA  TRP A   3      -3.184   1.900   4.605  1.00  0.00           C
ATOM     27  C   TRP A   3      -2.939   3.389   4.685  1.00  0.00           C
ATOM     28  O   TRP A   3      -2.901   4.120   3.684  1.00  0.00           O
ATOM     29  CB  TRP A   3      -2.011   1.214   3.851  1.00  0.00           C
ATOM     30  CG  TRP A   3      -0.639   1.384   4.511  1.00  0.00           C
ATOM     31  CD1 TRP A   3       0.385   2.247   4.070  1.00  0.00           C
ATOM     32  CD2 TRP A   3      -0.111   0.684   5.577  1.00  0.00           C
ATOM     33  CE2 TRP A   3       1.220   1.131   5.773  1.00  0.00           C
ATOM     34  CE3 TRP A   3      -0.658  -0.336   6.396  1.00  0.00           C
ATOM     35  NE1 TRP A   3       1.552   2.107   4.845  1.00  0.00           N
ATOM     36  CZ2 TRP A   3       2.010   0.567   6.798  1.00  0.00           C
ATOM     37  CZ3 TRP A   3       0.140  -0.873   7.407  1.00  0.00           C
ATOM     38  CH2 TRP A   3       1.454  -0.427   7.608  1.00  0.00           C
ATOM     49  N   TRP A   4      -2.711   3.897   5.951  1.00  0.00           N
ATOM     50  CA  TRP A   4      -2.097   5.179   6.287  1.00  0.00           C
ATOM     51  C   TRP A   4      -2.267   5.492   7.755  1.00  0.00           C
ATOM     52  O   TRP A   4      -2.330   4.614   8.648  1.00  0.00           O
ATOM     53  CB  TRP A   4      -2.755   6.297   5.432  1.00  0.00           C
ATOM     54  CG  TRP A   4      -2.461   6.213   3.931  1.00  0.00           C
ATOM     55  CD1 TRP A   4      -3.330   5.690   2.951  1.00  0.00           C
ATOM     56  CD2 TRP A   4      -1.294   6.534   3.269  1.00  0.00           C
ATOM     57  CE2 TRP A   4      -1.467   6.208   1.901  1.00  0.00           C
ATOM     58  CE3 TRP A   4      -0.064   7.064   3.735  1.00  0.00           C
ATOM     59  NE1 TRP A   4      -2.730   5.681   1.677  1.00  0.00           N
ATOM     60  CZ2 TRP A   4      -0.413   6.421   0.986  1.00  0.00           C
ATOM     61  CZ3 TRP A   4       0.960   7.272   2.810  1.00  0.00           C
ATOM     62  CH2 TRP A   4       0.788   6.958   1.455  1.00  0.00           C
END
"""
def exercise_lbfgs_simple (mon_lib_srv, ener_lib, verbose=False) :
  # three peptides:
  #  1 = poly-ALA, favored
  #  2 = poly-ALA, outlier
  #  3 = poly-TRP, outlier
  #
  # Note that the ramalyze score for the first actually gets slightly worse,
  # but it's still good and we're starting from an excellent score anyway.
  #
  residuals = [0.00024512, 307.616444, 294.913714]
  for i, peptide in enumerate([pdb1, pdb2, pdb3]) :
    pdb_in = iotbx.pdb.input(source_info="peptide",
      lines=flex.split_lines(peptide))
    params = pdb_interpretation.master_params.extract()
    processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      pdb_inp=pdb_in,
      log=StringIO())
    log = StringIO()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    atoms = pdb_hierarchy.atoms()
    sites_cart_1 = atoms.extract_xyz().deep_copy()
    gradients_fd = flex.vec3_double(sites_cart_1.size(), (0,0,0))
    gradients_an = flex.vec3_double(sites_cart_1.size(), (0,0,0))
    params = ramachandran.master_phil.fetch().extract()
    rama_manager = ramachandran.ramachandran_manager(
        pdb_hierarchy, None, params, log)
    assert rama_manager.get_n_proxies() == 1
    residual_an = rama_manager.target_and_gradients(
      unit_cell=None,
      sites_cart=sites_cart_1,
      gradient_array=gradients_an)
    assert approx_equal(residual_an, residuals[i], eps=0.00001)
  if verbose :
    print ""
  for i, peptide in enumerate([pdb1, pdb2, pdb3]) :
    pdb_in = iotbx.pdb.input(source_info="peptide",
      lines=flex.split_lines(peptide))
    o = benchmark_structure(pdb_in, mon_lib_srv, ener_lib, verbose)
    phi0, psi0 = o.r0.results[0].phi, o.r0.results[0].psi
    phi1, psi1 = o.r1.results[0].phi, o.r1.results[0].psi
    phi2, psi2 = o.r2.results[0].phi, o.r2.results[0].psi
    r0 = o.r0.results[0].score
    r1 = o.r1.results[0].score
    r2 = o.r2.results[0].score
    if verbose :
      print "peptide %d" % (i+1)
      print " before: rmsd_bonds=%-6.4f rmsd_angles=%-6.3f" % (o.b0,o.a0)
      print "         phi=%-6.1f psi=%-6.1f score=%-.2f" % (phi0, psi0, r0)
      print " simple: rmsd_bonds=%-6.4f rmsd_angles=%-6.3f" % (o.b1,o.a1)
      print "         phi=%-6.1f psi=%-6.1f score=%-.2f" % (phi1, psi1, r1)
      print " + Rama: rmsd_bonds=%-6.4f rmsd_angles=%-6.3f" % (o.b2,o.a2)
      print "         phi=%-6.1f psi=%-6.1f score=%-.2f" % (phi2, psi2, r2)
      print ""

def exercise_lbfgs_big (verbose=False) :
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3mku.pdb",
    test=os.path.isfile)
  if (file_name is None) :
    print "Skipping big test."
    return
  pdb_in = iotbx.pdb.input(source_info="peptide",
    file_name=file_name)
  o = benchmark_structure(pdb_in, verbose, 1.0)
  if verbose :
    show_results(o, "3mhk")

def show_results (o, structure_name) :
  print structure_name
  print " before: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b0, o.a0, o.r0.percent_outliers, o.r0.percent_favored)
  print " simple: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b1, o.a1, o.r1.percent_outliers, o.r1.percent_favored)
  print " + Rama: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b2, o.a2, o.r2.percent_outliers, o.r2.percent_favored)
  print ""

def benchmark_structure (pdb_in, mon_lib_srv, ener_lib, verbose=False, w=1.0) :
  params = pdb_interpretation.master_params.extract()
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    pdb_inp=pdb_in,
    log=StringIO())
  log = StringIO()
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  r0 = ramalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  atoms = pdb_hierarchy.atoms()
  sites_cart_1 = atoms.extract_xyz().deep_copy()
  sites_cart_2 = sites_cart_1.deep_copy()
  grm = processed_pdb_file.geometry_restraints_manager()
  assert (grm is not None)
  e = grm.energies_sites(sites_cart=sites_cart_1)
  b0 = e.bond_deviations()[-1]
  a0 = e.angle_deviations()[-1]
  flags = cctbx.geometry_restraints.flags.flags(default=True)
  lbfgs = geometry_minimization.lbfgs(
    sites_cart=sites_cart_1,
    geometry_restraints_manager=grm,
    geometry_restraints_flags=flags,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=500))
  a1 = lbfgs.rmsd_angles
  b1 = lbfgs.rmsd_bonds
  atoms.set_xyz(sites_cart_1)
  r1 = ramalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  rama_params = ramachandran.master_phil.fetch().extract()
  rama_manager = ramachandran.ramachandran_manager(
      pdb_hierarchy, None, rama_params, log)
  grm.set_ramachandran_restraints(rama_manager)
  lbfgs = geometry_minimization.lbfgs(
    sites_cart=sites_cart_2,
    geometry_restraints_manager=grm,
    geometry_restraints_flags=flags,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=500))
  a2 = lbfgs.rmsd_angles
  b2 = lbfgs.rmsd_bonds
  atoms.set_xyz(sites_cart_2)
  r2 = ramalyze(pdb_hierarchy=pdb_hierarchy, outliers_only=False)
  return group_args(
    a0=a0,
    a1=a1,
    a2=a2,
    b0=b0,
    b1=b1,
    b2=b2,
    r0=r0,
    r1=r1,
    r2=r2)

def exercise_geo_output(mon_lib_srv, ener_lib):
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
  pdb_inp = iotbx.pdb.input(source_info="peptide",lines=flex.split_lines(pdb_str))
  hierarchy = pdb_inp.construct_hierarchy()
  atoms = hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  params = ramachandran.master_phil.fetch().extract()
  params.rama_potential = "emsley"
  rama_manager = ramachandran.ramachandran_manager(
      hierarchy, None, params, StringIO())
  out = StringIO()
  rama_manager.show_sorted(
      by_value="residual",
      sites_cart=sites_cart,
      site_labels=[a.id_str() for a in atoms],
      f=out)
  gv = out.getvalue()
  # print out.getvalue()
  assert gv == """\
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

  params.rama_potential = "oldfield"
  rama_manager = ramachandran.ramachandran_manager(
      hierarchy, None, params, StringIO())
  out = StringIO()
  rama_manager.show_sorted(
      by_value="residual",
      sites_cart=sites_cart,
      site_labels=[a.id_str() for a in atoms],
      f=out)
  gv = out.getvalue()
  # print out.getvalue()
  assert gv == """\
Ramachandran plot restraints: 8
Sorted by residual:
phi-psi angles formed by             residual
    pdb=" C   ALA     7 "            6.90e-01
    pdb=" N   ALA     8 "
    pdb=" CA  ALA     8 "
    pdb=" C   ALA     8 "
    pdb=" N   ALA     9 "
phi-psi angles formed by             residual
    pdb=" C   ALA     1 "            2.01e-01
    pdb=" N   ALA     2 "
    pdb=" CA  ALA     2 "
    pdb=" C   ALA     2 "
    pdb=" N   ALA     3 "
phi-psi angles formed by             residual
    pdb=" C   ALA     5 "            2.27e-02
    pdb=" N   ALA     6 "
    pdb=" CA  ALA     6 "
    pdb=" C   ALA     6 "
    pdb=" N   ALA     7 "
phi-psi angles formed by             residual
    pdb=" C   ALA     6 "            1.71e-02
    pdb=" N   ALA     7 "
    pdb=" CA  ALA     7 "
    pdb=" C   ALA     7 "
    pdb=" N   ALA     8 "
phi-psi angles formed by             residual
    pdb=" C   ALA     2 "            1.43e-02
    pdb=" N   ALA     3 "
    pdb=" CA  ALA     3 "
    pdb=" C   ALA     3 "
    pdb=" N   ALA     4 "
phi-psi angles formed by             residual
    pdb=" C   ALA     8 "            1.00e-02
    pdb=" N   ALA     9 "
    pdb=" CA  ALA     9 "
    pdb=" C   ALA     9 "
    pdb=" N   ALA    10 "
phi-psi angles formed by             residual
    pdb=" C   ALA     4 "            3.55e-03
    pdb=" N   ALA     5 "
    pdb=" CA  ALA     5 "
    pdb=" C   ALA     5 "
    pdb=" N   ALA     6 "
phi-psi angles formed by             residual
    pdb=" C   ALA     3 "            3.16e-03
    pdb=" N   ALA     4 "
    pdb=" CA  ALA     4 "
    pdb=" C   ALA     4 "
    pdb=" N   ALA     5 "

"""


def exercise_other (mon_lib_srv, ener_lib) :
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3mku.pdb",
    test=os.path.isfile)
  if (file_name is None) :
    print "Skipping test."
    return
  params = pdb_interpretation.master_params.fetch().extract()
  params.peptide_link.ramachandran_restraints = True
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    pdb_inp=pdb_in,
    log=StringIO())
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  params2 = ramachandran.refine_opt_params.fetch().extract()
  params2.exclude_secondary_structure = True
  ramachandran.process_refinement_settings(
    params=params2,
    pdb_hierarchy=pdb_hierarchy,
    secondary_structure_manager=ss_mgr,
    d_min=2.9,
    log=StringIO())

def exercise_ramachandran_selections(mon_lib_srv, ener_lib):
  # Just check overall rama proxies
  file_name = libtbx.env.find_in_repositories(
    # relative_path="phenix_regression/pdb/3mku.pdb",
    relative_path="phenix_regression/pdb/fab_a_cut.pdb",
    test=os.path.isfile)
  if (file_name is None) :
    print "Skipping test."
    return
  params = pdb_interpretation.master_params.fetch().extract()
  params.peptide_link.ramachandran_restraints = True
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    file_name=file_name,
    log=StringIO())
  grm = processed_pdb_file.geometry_restraints_manager()
  # print grm.ramachandran_manager.get_n_proxies() # 47 is wrong here

  # simple selection
  params = pdb_interpretation.master_params.fetch().extract()
  params.peptide_link.ramachandran_restraints = True
  params.peptide_link.rama_selection = "chain A and resid 1:7"
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    file_name=file_name,
    log=StringIO())
  grm = processed_pdb_file.geometry_restraints_manager()
  nprox = grm.ramachandran_manager.get_n_proxies()
  assert nprox == 5, ""+\
      "Want to get 5 rama proxies, got %d" % nprox

  # 7 residues: there are insertion codes
  params = pdb_interpretation.master_params.fetch().extract()
  params.peptide_link.ramachandran_restraints = True
  params.peptide_link.rama_selection = "chain A and resid 27:28"
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=params,
    file_name=file_name,
    log=StringIO())
  grm = processed_pdb_file.geometry_restraints_manager()
  # print grm.ramachandran_manager.get_n_proxies() 0 is wrong here


if __name__ == "__main__" :
  import time
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()

  t0 = time.time()
  exercise_basic()
  t1 = time.time()
  exercise_lbfgs_simple(mon_lib_srv, ener_lib,
      ("--verbose" in sys.argv) or ("-v" in sys.argv))
  t2 = time.time()
  if ("--full" in sys.argv):
    exercise_lbfgs_big(mon_lib_srv, ener_lib,
        ("--verbose" in sys.argv) or ("-v" in sys.argv))
  t3 = time.time()
  exercise_geo_output(mon_lib_srv, ener_lib)
  t4 = time.time()
  #
  # XXX FIXME!!! Ramachandran selections does not work properly with insertion
  # codes. Actually, those are not getting restrained due to flawed logic in
  # mmtbx/rotamer/__init__.py: extract_phi_psi function.
  #
  exercise_ramachandran_selections(mon_lib_srv, ener_lib)
  t5 = time.time()
  print "Times: %.3f, %.3f, %.3f, %.3f, %.3f. Total: %.3f" % (
      t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t5-t0)
  print "OK"
