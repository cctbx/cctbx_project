
from __future__ import division
from mmtbx import ramachandran
from mmtbx.validation.ramalyze import ramalyze
from mmtbx.monomer_library import server, pdb_interpretation
from mmtbx.command_line import geometry_minimization
from iotbx import file_reader
import iotbx.pdb
import cctbx.geometry_restraints
import scitbx.lbfgs
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
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
  for i, peptide in enumerate([pdb1, pdb2, pdb3]) :
    pdb_in = iotbx.pdb.input(source_info="peptide",
      lines=flex.split_lines(peptide))
    mon_lib_srv = server.server()
    ener_lib = server.ener_lib()
    params = pdb_interpretation.master_params.extract()
    #params.peptide_link_params.ramachandran_restraints = True
    processed_pdb_file = pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      params=params,
      pdb_inp=pdb_in,
      log=StringIO())
    log = StringIO()
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    atoms = pdb_hierarchy.atoms()
    proxies = ramachandran.extract_proxies(pdb_hierarchy, log=log)
    assert (len(proxies) == 1)
    sites_cart_1 = atoms.extract_xyz().deep_copy()
    gradients_fd = flex.vec3_double(sites_cart_1.size(), (0,0,0))
    gradients_an = flex.vec3_double(sites_cart_1.size(), (0,0,0))
    params = ramachandran.master_phil.fetch().extract()
    params.use_finite_differences = False
    restraints_helper = ramachandran.generic_restraints_helper(params)
    residual_fd = restraints_helper.restraints_residual_sum(
      sites_cart=sites_cart_1,
      proxies=proxies,
      gradient_array=gradients_fd)
    params.use_finite_differences = True
    #restraints_helper = ramachandran.generic_restraints_helper(params)
    residual_an = restraints_helper.restraints_residual_sum(
      sites_cart=sites_cart_1,
      proxies=proxies,
      gradient_array=gradients_an)
    assert approx_equal(residual_fd, residual_an, eps=0.01)
    if verbose :
      print "peptide %d gradients" % (i+1)
    for g1, g2 in zip(gradients_fd, gradients_an) :
      if (g1 != (0.0,0.0,0.0)) or (g2 != (0.0,0.0,0.0)) :
        if verbose :
          print ("  (%.2f, %.2f, %.2f)" % tuple(g1)), \
                ("(%.2f, %.2f, %.2f)" % tuple(g2))
        for u in range(3) :
          assert approx_equal(g1[u], g2[u], eps=0.01)
  if verbose :
    print ""
  for i, peptide in enumerate([pdb1, pdb2, pdb3]) :
    pdb_in = iotbx.pdb.input(source_info="peptide",
      lines=flex.split_lines(peptide))
    o = benchmark_structure(pdb_in, verbose)
    phi0, psi0 = o.rama_list0[0][4:6]
    phi1, psi1 = o.rama_list1[0][4:6]
    phi2, psi2 = o.rama_list2[0][4:6]
    r0 = float(o.rama_list0[0][3])
    r1 = float(o.rama_list1[0][3])
    r2 = float(o.rama_list2[0][3])
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
  pdb_in = file_reader.any_file(file_name).file_object
  o = benchmark_structure(pdb_in, verbose, 1.0)
  if verbose :
    show_results(o, "3mhk")

def show_results (o, structure_name) :
  r0 = [ row[-3] for row in o.rama_list0 ].count("OUTLIER") / len(o.rama_list0)
  r1 = [ row[-3] for row in o.rama_list1 ].count("OUTLIER") / len(o.rama_list1)
  r2 = [ row[-3] for row in o.rama_list2 ].count("OUTLIER") / len(o.rama_list2)
  rr0 = [ row[-3] for row in o.rama_list0 ].count("Favored") / len(o.rama_list0)
  rr1 = [ row[-3] for row in o.rama_list1 ].count("Favored") / len(o.rama_list1)
  rr2 = [ row[-3] for row in o.rama_list2 ].count("Favored") / len(o.rama_list2)
  print structure_name
  print " before: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b0,o.a0,r0*100,rr0*100)
  print " simple: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b1,o.a1,r1*100,rr1*100)
  print " + Rama: bonds=%-6.4f angles=%-6.3f outliers=%.1f%% favored=%.1f%%"\
    % (o.b2,o.a2,r2*100,rr2*100)
  print ""

def benchmark_structure (pdb_in, verbose=False, w=1.0) :
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
  sites_cart_1 = atoms.extract_xyz().deep_copy()
  sites_cart_2 = sites_cart_1.deep_copy()
  proxies = ramachandran.extract_proxies(pdb_hierarchy, log=log)
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
  rama1, rama_list1 = ramalyze().analyze_pdb(hierarchy=pdb_hierarchy)
  params = ramachandran.master_phil.fetch().extract()
  restraints_helper = ramachandran.generic_restraints_helper(params)
  grm.set_generic_restraints(
    proxies=proxies,
    restraints_helper=restraints_helper)
  lbfgs = geometry_minimization.lbfgs(
    sites_cart=sites_cart_2,
    geometry_restraints_manager=grm,
    geometry_restraints_flags=flags,
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=500))
  a2 = lbfgs.rmsd_angles
  b2 = lbfgs.rmsd_bonds
  atoms.set_xyz(sites_cart_2)
  rama2, rama_list2 = ramalyze().analyze_pdb(hierarchy=pdb_hierarchy)
  return group_args(
    a0=a0,
    a1=a1,
    a2=a2,
    b0=b0,
    b1=b1,
    b2=b2,
    rama_list0=rama_list0,
    rama_list1=rama_list1,
    rama_list2=rama_list2)

def exercise_other () :
  file_name = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3mku.pdb",
    test=os.path.isfile)
  if (file_name is None) :
    print "Skipping test."
    return
    mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
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

if __name__ == "__main__" :
  exercise_basic()
  exercise_lbfgs_simple(("--verbose" in sys.argv) or ("-v" in sys.argv))
  exercise_lbfgs_big(("--verbose" in sys.argv) or ("-v" in sys.argv))
  print "OK"
