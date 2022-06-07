from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from six.moves import cStringIO as StringIO
import cctbx.adp_restraints
from libtbx.test_utils import approx_equal
import libtbx.load_env
from six.moves import range

phe_pdb = """\
remark CRYST1   25.000   35.000   45.000  80.00 70.00 100.00 P 1           1
CRYST1    4.000    5.000    7.000  80.00 70.00 100.00 P 1           1
ATOM      1  CB  PHE A   1       7.767   5.853   7.671  1.00 15.00           C
ATOM      2  CG  PHE A   1       6.935   5.032   8.622  1.00 44.21           C
ANISOU    2  CG  PHE A   1     4894   5943   5958    725   -324    -10       C
ATOM      3  CD1 PHE A   1       5.918   4.176   8.140  1.00 17.66           C
ANISOU    3  CD1 PHE A   1     1321   2083   3304   1237   1214   2356       C
ATOM      4  CD2 PHE A   1       7.161   5.107  10.012  1.00 36.94           C
ANISOU    4  CD2 PHE A   1     4099   5961   3974    855     31     95       C
ATOM      5  CE1 PHE A   1       5.126   3.395   9.038  1.00 17.96           C
ATOM      6  CE2 PHE A   1       6.382   4.336  10.930  1.00 98.17           C
ATOM      7  CZ  PHE A   1       5.360   3.476  10.439  1.00 78.27           C
ATOM      8  C   PHE A   1       7.956   7.811   6.133  1.00 17.21           C
ANISOU    8  C   PHE A   1      816   2415   3307   1079   1332   2287       C
ATOM      9  O   PHE A   1       8.506   7.237   5.169  1.00 70.55           O
ANISOU    9  O   PHE A   1     8784   9152   8868     -8    122   -190       O
ATOM     10  OXT PHE A   1       8.143   9.010   6.428  1.00 27.61           O
ATOM     13  N   PHE A   1       5.875   6.461   6.183  1.00 17.13           N
ATOM     15  CA  PHE A   1       7.000   7.000   7.000  1.00 87.18           C
END
"""

def fd(xray_structure, restraints_manager, eps=1.e-2):
    uc = xray_structure.unit_cell()
    adp_ro = cctbx.adp_restraints.adp_aniso_restraints(
                                       xray_structure     = xray_structure,
                                       restraints_manager = restraints_manager,
                                       use_hd = False)
    g_aniso = adp_ro.gradients_aniso_cart
    g_iso = adp_ro.gradients_iso
    for i_seq, scatterer in enumerate(xray_structure.scatterers()):
        fl = scatterer.flags
        #print "scatterer :", i_seq, fl.use_u_iso(), fl.grad_u_iso(), \
        #                     fl.use_u_aniso(), fl.grad_u_aniso()
        if(fl.use_u_aniso()):
           for i_ind in range(6):
               xrs1 = xray_structure.deep_copy_scatterers()
               xrs2 = xray_structure.deep_copy_scatterers()
               sc1  = xrs1.scatterers()
               sc2  = xrs2.scatterers()
               us1 = sc1.extract_u_cart(uc)
               us2 = sc1.extract_u_cart(uc)
               m1 = flex.double(us1[i_seq])
               m2 = flex.double(us2[i_seq])
               m1[i_ind] = m1[i_ind] - eps
               m2[i_ind] = m2[i_ind] + eps
               us1[i_seq] = list(m1)
               us2[i_seq] = list(m2)
               sc1.set_u_cart(uc, us1)
               sc2.set_u_cart(uc, us2)
               adp_ro1 = cctbx.adp_restraints.adp_aniso_restraints(
                                       xray_structure     = xrs1,
                                       restraints_manager = restraints_manager,
                                       use_hd = False)
               adp_ro2 = cctbx.adp_restraints.adp_aniso_restraints(
                                       xray_structure     = xrs2,
                                       restraints_manager = restraints_manager,
                                       use_hd = False)
               g1 = (adp_ro2.target - adp_ro1.target)/(2*eps)
               g2 = g_aniso[i_seq][i_ind]
               #print "   fin.diff.= %10.5f anal.= %10.5f diff.= %10.5f"%(g1, g2, g1-g2)
               assert approx_equal(g1,g2,1.e-4)
        if(fl.use_u_iso()):
           sel = xray_structure.use_u_iso()
           xrs1 = xray_structure.deep_copy_scatterers()
           xrs2 = xray_structure.deep_copy_scatterers()
           sc1  = xrs1.scatterers()
           sc2  = xrs2.scatterers()
           us1 = sc1.extract_u_iso()
           us2 = sc1.extract_u_iso()
           us1[i_seq] = us1[i_seq] - eps
           us2[i_seq] = us2[i_seq] + eps
           sc1.set_u_iso(us1, sel, uc)
           sc2.set_u_iso(us2, sel, uc)
           adp_ro1 = cctbx.adp_restraints.adp_aniso_restraints(
                                      xray_structure     = xrs1,
                                      restraints_manager = restraints_manager,
                                      use_hd = False)
           adp_ro2 = cctbx.adp_restraints.adp_aniso_restraints(
                                      xray_structure     = xrs2,
                                      restraints_manager = restraints_manager,
                                      use_hd = False)
           g1 = (adp_ro2.target - adp_ro1.target)/(2*eps)
           g2 = g_iso[i_seq]
           #print "   fin.diff.= %10.5f anal.= %10.5f diff.= %10.5f"%(g1, g2, g1-g2)
           assert approx_equal(g1,g2,1.e-4)

def exercise():
  if (not libtbx.env.has_module("mmtbx")):
    print("Skipping exercise(): mmtbx module not available")
    return
  if (libtbx.env.find_in_repositories(relative_path="chem_data") is None):
    print("Skipping exercise(): chem_data directory not available")
    return
  from mmtbx.monomer_library import pdb_interpretation
  file_name = "phe_tst_adp_aniso_restraints.pdb"
  with open(file_name, "w") as f:
    f.write(phe_pdb)
  out = StringIO()
  processed_pdb_file = pdb_interpretation.run(
                                        args                     = [file_name],
                                        strict_conflict_handling = False,
                                        log                      = out)
  geo = processed_pdb_file.geometry_restraints_manager()
  xray_structure = processed_pdb_file.xray_structure()
  xray_structure.scatterers().flags_set_grads(state=False)
  xray_structure.scatterers().flags_set_grad_u_iso(
    iselection=xray_structure.use_u_iso().iselection())
  xray_structure.scatterers().flags_set_grad_u_aniso(
    iselection=xray_structure.use_u_aniso().iselection())
  adp_rm = cctbx.adp_restraints.adp_aniso_restraints(
                                           xray_structure     = xray_structure,
                                           restraints_manager = geo,
                                           use_hd = False)
  assert approx_equal(flex.mean(adp_rm.gradients_iso), 0.713756592583)
  assert approx_equal(flex.mean(adp_rm.gradients_aniso_cart.as_double()), -0.118959432097)
  assert approx_equal(adp_rm.target, 8.97112989232)
  fd(xray_structure = xray_structure, restraints_manager = geo, eps=1.e-4)

def run(args):
  assert len(args) == 0
  exercise()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
