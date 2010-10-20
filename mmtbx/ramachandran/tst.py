
from mmtbx import ramachandran
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.monomer_library import server, pdb_interpretation
from mmtbx.command_line import geometry_minimization
import iotbx.pdb
import cctbx.geometry_restraints
import scitbx.lbfgs
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from cStringIO import StringIO

def exercise_basic () :
  t = ramachandran.tables["ala"]
  assert approx_equal(t.get_score(0.0,0.0), -26.16, eps=0.01)
  assert approx_equal(t.get_score(-60,120), 10.41, eps=0.01)
  assert approx_equal(t.get_score(90,90), -7.43, eps=0.01)
  assert approx_equal(t.get_energy(0.0,0.0), 53.81, eps=0.01)
  assert approx_equal(t.get_energy(-60,120), 17.24, eps=0.01)

def exercise_lbfgs_simple (verbose=False) :
  # three peptides:
  #  1 = poly-ALA, favored
  #  2 = poly-ALA, outlier
  #  3 = poly-TRP, outlier
  #
  # Note that the ramalyze score for the first actually gets slightly worse,
  # but it's still good and we're starting from an excellent score anyway.
  #
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
ATOM     15  H   TRP A   2      -7.188   2.620   1.026  1.00  0.00           H
ATOM     16  HA  TRP A   2      -7.228   1.875   3.695  1.00  0.00           H
ATOM     17 2HB  TRP A   2      -5.686   0.023   2.715  1.00  0.00           H
ATOM     18 3HB  TRP A   2      -5.628   0.906   1.206  1.00  0.00           H
ATOM     19 1HD  TRP A   2      -7.656   0.544  -0.487  1.00  0.00           H
ATOM     20 1HE  TRP A   2      -9.932  -0.827  -0.520  1.00  0.00           H
ATOM     21 3HE  TRP A   2      -7.606  -0.756   4.438  1.00  0.00           H
ATOM     22 2HZ  TRP A   2     -11.361  -2.251   1.484  1.00  0.00           H
ATOM     23 3HZ  TRP A   2      -9.469  -2.109   5.345  1.00  0.00           H
ATOM     24 2HH  TRP A   2     -11.317  -2.853   3.892  1.00  0.00           H
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
ATOM     39  H   TRP A   3      -4.705   0.635   3.775  1.00  0.00           H
ATOM     40  HA  TRP A   3      -3.251   1.516   5.640  1.00  0.00           H
ATOM     41 2HB  TRP A   3      -1.935   1.609   2.820  1.00  0.00           H
ATOM     42 3HB  TRP A   3      -2.223   0.135   3.716  1.00  0.00           H
ATOM     43 1HD  TRP A   3       0.296   2.898   3.210  1.00  0.00           H
ATOM     44 1HE  TRP A   3       2.456   2.577   4.724  1.00  0.00           H
ATOM     45 3HE  TRP A   3      -1.667  -0.692   6.244  1.00  0.00           H
ATOM     46 2HZ  TRP A   3       3.025   0.903   6.951  1.00  0.00           H
ATOM     47 3HZ  TRP A   3      -0.264  -1.646   8.045  1.00  0.00           H
ATOM     48 2HH  TRP A   3       2.045  -0.858   8.401  1.00  0.00           H
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
ATOM     63  H   TRP A   4      -3.002   3.318   6.730  1.00  0.00           H
ATOM     64  HA  TRP A   4      -1.012   5.126   6.080  1.00  0.00           H
ATOM     65 2HB  TRP A   4      -2.422   7.293   5.781  1.00  0.00           H
ATOM     66 3HB  TRP A   4      -3.852   6.305   5.592  1.00  0.00           H
ATOM     67 1HD  TRP A   4      -4.316   5.299   3.166  1.00  0.00           H
ATOM     68 1HE  TRP A   4      -3.115   5.322   0.796  1.00  0.00           H
ATOM     69 3HE  TRP A   4       0.080   7.302   4.779  1.00  0.00           H
ATOM     70 2HZ  TRP A   4      -0.539   6.173  -0.058  1.00  0.00           H
ATOM     71 3HZ  TRP A   4       1.901   7.682   3.148  1.00  0.00           H
ATOM     72 2HH  TRP A   4       1.597   7.133   0.763  1.00  0.00           H
END
"""
  for i, peptide in enumerate([pdb1, pdb2, pdb3]) :
    pdb_in = iotbx.pdb.input(source_info="peptide",
      lines=flex.split_lines(peptide))
    mon_lib_srv = server.server()
    ener_lib = server.ener_lib()
    params = pdb_interpretation.master_params.extract()
    processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      pdb_inp=pdb_in,
      log=StringIO())
    log = StringIO()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    rama0, rama_list0 = ramalyze().analyze_pdb(hierarchy=pdb_hierarchy)
    atoms = pdb_hierarchy.atoms()
    proxies = ramachandran.extract_proxies(pdb_hierarchy, log=log)
    assert (len(proxies) == 1)
    grm = processed_pdb_file.geometry_restraints_manager()
    assert (grm is not None)
    flags = cctbx.geometry_restraints.flags.flags(default=True)
    sites_cart_1 = atoms.extract_xyz().deep_copy()
    e = grm.energies_sites(sites_cart=sites_cart_1)
    b0 = e.bond_deviations()[-1]
    a0 = e.angle_deviations()[-1]
    sites_cart_2 = sites_cart_1.deep_copy()
    lbfgs = geometry_minimization.lbfgs(
      sites_cart=sites_cart_2,
      geometry_restraints_manager=grm,
      geometry_restraints_flags=flags,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
          max_iterations=500))
    a1, b1 = lbfgs.rmsd_angles, lbfgs.rmsd_bonds
    atoms.set_xyz(sites_cart_2)
    rama1, rama_list1 = ramalyze().analyze_pdb(hierarchy=pdb_hierarchy)
    grm.set_generic_restraints(
      proxies=proxies,
      restraints_helper=ramachandran.generic_restraints_helper())
    sites_cart_3 = sites_cart_1.deep_copy()
    lbfgs = geometry_minimization.lbfgs(
      sites_cart=sites_cart_3,
      geometry_restraints_manager=grm,
      geometry_restraints_flags=flags,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
          max_iterations=500))
    a2, b2 = lbfgs.rmsd_angles, lbfgs.rmsd_bonds
    atoms.set_xyz(sites_cart_3)
    rama2, rama_list2 = ramalyze().analyze_pdb(hierarchy=pdb_hierarchy)
    phi0, psi0 = rama_list0[0][4:6]
    phi1, psi1 = rama_list1[0][4:6]
    phi2, psi2 = rama_list2[0][4:6]
    r0 = float(rama_list0[0][3])
    r1 = float(rama_list1[0][3])
    r2 = float(rama_list2[0][3])
    assert (r2 > 1.0)
    if (i > 0) : assert (r2 > r1)
    if verbose :
      print "peptide %d" % i
      print " before: rmsbonds=%-6.4f rmsangles=%-6.3f rama=%-.2f" % (b0,a0,r0)
      print "         phi=%-.1f psi=%-.1f" % (phi0, psi0)
      print " simple: rmsbonds=%-6.4f rmsangles=%-6.3f rama=%-.2f" % (b1,a1,r1)
      print "         phi=%-.1f psi=%-.1f" % (phi1, psi1)
      print " + Rama: rmsbonds=%-6.4f rmsangles=%-6.3f rama=%-.2f" % (b2,a2,r2)
      print "         phi=%-.1f psi=%-.1f" % (phi2, psi2)
      print ""

if __name__ == "__main__" :
  exercise_basic()
  exercise_lbfgs_simple(True)
  print "OK"
