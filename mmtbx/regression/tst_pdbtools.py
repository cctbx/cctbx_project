from __future__ import absolute_import, division, print_function
#import libtbx.load_env
import sys, os, math
from cctbx.array_family import flex
from libtbx.utils import remove_files
from mmtbx import utils
from libtbx.test_utils import approx_equal, not_approx_equal, run_command, \
  show_diff
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import adptbx
import mmtbx.model
from six.moves import cStringIO as StringIO, zip
from libtbx import easy_run

full_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
full_params.pdb_interpretation.flip_symmetric_amino_acids=False

class xray_structure_plus(object):
  def __init__(self, file_name):
    log = StringIO()
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    self.model = mmtbx.model.manager(model_input = pdb_inp, log = log)
    self.model.process(pdb_interpretation_params = full_params)
    self.xray_structure = self.model.get_xray_structure()
    uc = self.model.get_xray_structure().unit_cell()
    self.occ           = self.xray_structure.scatterers().extract_occupancies()
    self.u_iso            = self.xray_structure.scatterers().extract_u_iso()
    self.u_cart           = self.xray_structure.scatterers().extract_u_cart(uc)
    self.sites_cart       = self.xray_structure.sites_cart()
    self.use_u_iso        = self.xray_structure.use_u_iso()
    self.use_u_aniso      = self.xray_structure.use_u_aniso()
    self.u_iso_used       = self.u_iso.select(self.use_u_iso)
    self.u_cart_used      = self.u_cart.select(self.use_u_aniso)
    self.u_iso_not_used   = self.u_iso.select(~self.use_u_iso)
    self.u_cart_not_used  = self.u_cart.select(~self.use_u_aniso)

  def selection(self, selection_strings):
    return utils.get_atom_selections(iselection       = False,
                                     model = self.model,
                                     selection_strings= selection_strings)[0]

def exercise_basic():
  pdb_str = """
CRYST1   12.000   11.000   13.000  80.00 70.00 100.00 P 1           1
ATOM      1  CB  PHE A   1       7.767   5.853   7.671  1.00 17.45           C
ANISOU    1  CB  PHE A   1     1090   2307   3233   1244   1260   2218       C
ATOM      2  CG  PHE A   1       6.935   5.032   8.622  1.00 17.61           C
ANISOU    2  CG  PHE A   1     1258   2172   3262   1259   1261   2264       C
ATOM      3  CD1 PHE A   1       5.918   4.176   8.140  1.00 17.66           C
ANISOU    3  CD1 PHE A   1     1321   2083   3304   1237   1214   2356       C
ATOM      4  CD2 PHE A   1       7.161   5.107  10.012  1.00 17.82           C
ANISOU    4  CD2 PHE A   1     1374   2150   3248   1279   1308   2217       C
ATOM      5  CE1 PHE A   1       5.126   3.395   9.038  1.00 17.96           C
ANISOU    5  CE1 PHE A   1     1508   1983   3333   1224   1216   2398       C
ATOM      6  CE2 PHE A   1       6.382   4.336  10.930  1.00 18.17           C
ANISOU    6  CE2 PHE A   1     1573   2056   3276   1258   1310   2261       C
ATOM      7  CZ  PHE A   1       5.360   3.476  10.439  1.00 18.27           C
ANISOU    7  CZ  PHE A   1     1645   1979   3318   1226   1265   2351       C
ATOM      8  C   PHE A   1       7.956   7.811   6.133  1.00 17.21           C
ANISOU    8  C   PHE A   1      816   2415   3307   1079   1332   2287       C
ATOM      9  O   PHE A   1       8.506   7.237   5.169  1.00 17.39           O
ANISOU    9  O   PHE A   1      818   2555   3234   1097   1244   2232       O
ATOM     10  OXT PHE A   1       8.143   9.010   6.428  1.00 17.15           O
ANISOU   10  OXT PHE A   1      737   2422   3356   1020   1411   2296       O
ATOM     13  N   PHE A   1       5.875   6.461   6.183  1.00 17.13           N
ANISOU   13  N   PHE A   1      934   2156   3419   1076   1303   2452       N
ATOM     15  CA  PHE A   1       7.000   7.000   7.000  1.00 17.18           C
ANISOU   15  CA  PHE A   1      931   2249   3347   1119   1342   2337       C
ATOM     13  CB  PHE B   2      11.715   4.672   7.185  1.00 34.89           C
ATOM     14  CG  PHE B   2      10.876   4.117   8.301  1.00 35.22           C
ATOM     15  CD2 PHE B   2      10.127   2.966   8.118  1.00 35.30           C
ATOM     16  CD1 PHE B   2      10.836   4.746   9.534  1.00 35.63           C
ATOM     17  CE2 PHE B   2       9.355   2.454   9.143  1.00 35.88           C
ATOM     18  CE1 PHE B   2      10.066   4.239  10.563  1.00 36.30           C
ATOM     19  CZ  PHE B   2       9.324   3.091  10.367  1.00 36.47           C
ATOM     20  C   PHE B   2      11.961   6.316   5.313  1.00 34.38           C
ATOM     21  O   PHE B   2      11.902   5.976   4.132  1.00 34.72           O
ATOM     22  OXT PHE B   2      12.817   7.140   5.635  1.00 34.24           O
ATOM     23  N   PHE B   2       9.817   5.175   5.710  1.00 34.25           N
ATOM     24  CA  PHE B   2      11.002   5.735   6.347  1.00 34.35           C
ATOM      1  CB  PHE C   3      10.767   8.853   7.671  1.00 52.27           C
ANISOU    1  CB  PHE C   3     3279   6815   9766   3688   3821   6736       C
ATOM      2  CG  PHE C   3       9.935   8.032   8.622  1.00 52.60           C
ANISOU    2  CG  PHE C   3     3617   6545   9825   3719   3822   6829       C
ATOM      3  CD1 PHE C   3       8.918   7.176   8.140  1.00 52.70           C
ANISOU    3  CD1 PHE C   3     3745   6367   9912   3673   3727   7011       C
ATOM      4  CD2 PHE C   3      10.161   8.107  10.012  1.00 53.03           C
ANISOU    4  CD2 PHE C   3     3850   6502   9797   3761   3917   6735       C
ATOM      5  CE1 PHE C   3       8.126   6.395   9.038  1.00 53.33           C
ANISOU    5  CE1 PHE C   3     4119   6169   9973   3647   3734   7098       C
ATOM      6  CE2 PHE C   3       9.382   7.336  10.930  1.00 53.74           C
ANISOU    6  CE2 PHE C   3     4248   6315   9853   3723   3923   6823       C
ATOM      7  CZ  PHE C   3       8.360   6.476  10.439  1.00 53.93           C
ANISOU    7  CZ  PHE C   3     4391   6159   9943   3655   3834   7005       C
ATOM      8  C   PHE C   3      10.956  10.811   6.133  1.00 51.80           C
ANISOU    8  C   PHE C   3     2733   7032   9917   3360   3966   6875       C
ATOM      9  O   PHE C   3      11.506  10.237   5.169  1.00 52.17           O
ANISOU    9  O   PHE C   3     2737   7311   9773   3399   3791   6767       O
ATOM     10  OXT PHE C   3      11.143  12.010   6.428  1.00 51.70           O
ANISOU   10  OXT PHE C   3     2580   7045  10018   3240   4125   6894       O
ATOM     13  N   PHE C   3       8.875   9.461   6.183  1.00 51.64           N
ANISOU   13  N   PHE C   3     2968   6515  10138   3352   3906   7205       N
ATOM     15  CA  PHE C   3      10.000  10.000   7.000  1.00 51.73           C
ANISOU   15  CA  PHE C   3     2963   6697   9995   3440   3985   6975       C
END
"""
  verbose=False
  file_name = "exercise_basic_pdbtools.pdb"
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str)
  pi.write_pdb_file(file_name=file_name)
  output = "exercise_basic_pdbtools_mod.pdb"
  xrsp_init = xray_structure_plus(file_name = file_name)
  assert file_name.find('"') < 0
  base = \
      'phenix.pdbtools "%s" suffix=none output.prefix=%s '%(file_name, output.replace(".pdb",""))
  for selection_str in [None, "chain A or chain C"]:
    selection = xrsp_init.selection(selection_strings = selection_str)
    if(selection_str is None):
      assert selection.size() == selection.count(True)
    else:
      assert selection.size() == 36 and selection.count(True) == 24
    #
    cmd = base + 'adp.randomize=true modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.set_b_iso=10.0 modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_set_b_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.shift_b_iso=20.0 modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.scale_adp=2.0 modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.convert_to_iso=true modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_to_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.convert_to_aniso=true modify.selection="%s"'%str(selection_str)
    print(cmd)
    check_adp_to_aniso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    shake = 1.5
    cmd = base+'sites.shake=%s modify.selection="%s"'%(str(shake), str(selection_str))
    print(cmd)
    check_sites_shake(
      cmd, xrsp_init, output, selection, selection_str, shake, verbose)
    #
    cmd = base+'sites.rotate="1,2,3" sites.translate="4,5,6" modify.selection="%s"'%(
                                                            str(selection_str))
    print(cmd)
    check_sites_rt(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base+'occupancies.randomize=true modify.selection="%s"'%(str(selection_str))
    print(cmd)
    check_occ_randomize(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base+'occupancies.set=0.75 modify.selection="%s"'%(str(selection_str))
    print(cmd)
    check_occ_set(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    remove_selection_str = "element C"
    cmd = base+'remove="%s" modify.selection="%s"'%(
      str(remove_selection_str), str(selection_str))
    print(cmd)
    check_remove_selection(
      cmd, xrsp_init, output, selection, selection_str,
      remove_selection_str, verbose)
    #
    keep_selection_str = "element C"
    cmd = base+'keep="%s" modify.selection="%s"'%(
      str(keep_selection_str), str(selection_str))
    print(cmd)
    check_keep_selection(
      cmd, xrsp_init, output, selection, selection_str,
      keep_selection_str, verbose)
    #
  #
  cmd = base
  print(cmd)
  check_all_none(cmd, xrsp_init, output, verbose)
  #
  cmd = base
  check_keep_remove_conflict(cmd, output, verbose)

def check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,        xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.sites_cart, xrsp_init.sites_cart,tolerance)
  assert approx_equal(xrsp.use_u_iso,  xrsp_init.use_u_iso,tolerance)
  assert approx_equal(xrsp.use_u_aniso,xrsp_init.use_u_aniso,tolerance)
  assert approx_equal(xrsp.u_iso_not_used,  xrsp_init.u_iso_not_used,tolerance)
  assert approx_equal(xrsp.u_cart_not_used,xrsp_init.u_cart_not_used,tolerance)
  if(selection_str is None):
    assert not_approx_equal(xrsp.u_iso_used,  xrsp_init.u_iso_used,tolerance)
    assert not_approx_equal(xrsp.u_cart_used, xrsp_init.u_cart_used,tolerance)
  else:
    arg1 = xrsp.u_iso_used.select(selection.select(xrsp.use_u_iso))
    arg2 = xrsp_init.u_iso_used.select(selection.select(xrsp_init.use_u_iso))
    if(arg1.size() > 0): assert not_approx_equal(arg1, arg2,tolerance)
    arg1 =xrsp.u_cart_used.select(selection.select(xrsp.use_u_aniso))
    arg2 =xrsp_init.u_cart_used.select(selection.select(xrsp_init.use_u_aniso))
    if(arg1.size() > 0): assert not_approx_equal(arg1, arg2,tolerance)

def check_adp_set_b_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,        xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.sites_cart, xrsp_init.sites_cart,tolerance)
  assert approx_equal(xrsp.use_u_iso,  xrsp_init.use_u_iso,tolerance)
  assert approx_equal(xrsp.use_u_aniso,xrsp_init.use_u_aniso,tolerance)
  assert approx_equal(xrsp.u_iso_not_used,  xrsp_init.u_iso_not_used,tolerance)
  assert approx_equal(xrsp.u_cart_not_used,xrsp_init.u_cart_not_used,tolerance)
  if(selection_str is None):
    assert not_approx_equal(xrsp.u_iso_used,  xrsp_init.u_iso_used,tolerance)
    for ucart in xrsp.u_cart:
      b_iso = adptbx.u_as_b(adptbx.u_cart_as_u_iso(ucart))
      if b_iso > 0: assert approx_equal(b_iso, 10, 0.005)
      else: assert approx_equal(b_iso, -78.956, 0.005)
  else:
    arg1 = xrsp.u_iso_used.select(selection.select(xrsp.use_u_iso))
    arg2 = xrsp_init.u_iso_used.select(selection.select(xrsp_init.use_u_iso))
    if(arg1.size() > 0): assert not_approx_equal(arg1, arg2,tolerance)
    for ucart in xrsp.u_cart:
      b_iso = adptbx.u_as_b(adptbx.u_cart_as_u_iso(ucart))
      if b_iso > 0: assert approx_equal(b_iso, 10, 0.005)
      else: assert approx_equal(b_iso, -78.956, 0.005)

def check_adp_to_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,        xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.sites_cart, xrsp_init.sites_cart,tolerance)
  assert not_approx_equal(xrsp.use_u_iso,  xrsp_init.use_u_iso,tolerance)
  assert not_approx_equal(xrsp.use_u_aniso,xrsp_init.use_u_aniso,tolerance)
  assert xrsp.u_iso_not_used.size() == 0
  assert xrsp_init.u_iso_not_used.size() > 0
  assert xrsp.u_cart_used.size() == 0
  assert xrsp_init.u_cart_used.size() > 0

def check_adp_to_aniso(
      cmd, xrsp_init, output, selection, selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,        xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.sites_cart, xrsp_init.sites_cart,tolerance)
  if(selection_str is None):
    assert not_approx_equal(xrsp.use_u_iso,  xrsp_init.use_u_iso,tolerance)
    assert not_approx_equal(xrsp.use_u_aniso,xrsp_init.use_u_aniso,tolerance)
    assert xrsp.u_iso_used.size() == 0
    assert xrsp_init.u_iso_used.size() > 0
    assert xrsp.u_cart_not_used.size() == 0
    assert xrsp_init.u_cart_not_used.size() > 0
  else:
    assert approx_equal(xrsp.use_u_iso,  xrsp_init.use_u_iso,tolerance)
    assert approx_equal(xrsp.use_u_aniso,xrsp_init.use_u_aniso,tolerance)

def check_sites_shake(
      cmd, xrsp_init, output, selection, selection_str, shake, verbose,
      tolerance=1.e-2):
  remove_files(output)
  print(cmd)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,    xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.u_iso,  xrsp_init.u_iso,tolerance)
  assert approx_equal(xrsp.u_cart, xrsp_init.u_cart,tolerance)
  if(selection_str is None):
    diff = xrsp.sites_cart - xrsp_init.sites_cart
    assert approx_equal(
                      math.sqrt(flex.mean(diff.dot())), shake, 1.e-3,tolerance)
  else:
    diff = xrsp.sites_cart - xrsp_init.sites_cart
    assert approx_equal(
              math.sqrt(flex.mean(diff.select(selection).dot())), shake, 1.e-3)
    assert approx_equal(
              math.sqrt(flex.mean(diff.select(~selection).dot())),0.,tolerance)

def check_sites_rt(
      cmd, xrsp_init, output, selection, selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,    xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.u_iso,  xrsp_init.u_iso,tolerance)
  assert approx_equal(xrsp.u_cart, xrsp_init.u_cart,tolerance)
  if(selection_str is None):
    diff = xrsp.sites_cart - xrsp_init.sites_cart
    assert math.sqrt(flex.mean(diff.dot())) > 1.0
  else:
    diff = xrsp.sites_cart - xrsp_init.sites_cart
    assert math.sqrt(flex.mean(diff.select(selection).dot())) > 1.0
    assert approx_equal(
              math.sqrt(flex.mean(diff.select(~selection).dot())),0.,tolerance)

def check_occ_randomize(
      cmd, xrsp_init, output, selection,selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.sites_cart,xrsp_init.sites_cart,tolerance)
  assert approx_equal(xrsp.u_iso,     xrsp_init.u_iso,tolerance)
  assert approx_equal(xrsp.u_cart,    xrsp_init.u_cart,tolerance)
  if(selection_str is None):
    diff = flex.abs(xrsp.occ - xrsp_init.occ)
    assert flex.mean(diff) > 0.0
    assert flex.max(diff) > 0.0
  else:
    diff = flex.abs(xrsp.occ - xrsp_init.occ)
    assert flex.mean(diff) > 0.0
    assert flex.max(diff) > 0.0
    assert approx_equal(flex.mean(diff.select(~selection)),0.,tolerance)

def check_occ_set(
      cmd, xrsp_init, output, selection,selection_str, verbose,
      tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.sites_cart,xrsp_init.sites_cart,tolerance)
  assert approx_equal(xrsp.u_iso,     xrsp_init.u_iso,tolerance)
  assert approx_equal(xrsp.u_cart,    xrsp_init.u_cart,tolerance)
  if(selection_str is None):
    diff = flex.abs(xrsp.occ - xrsp_init.occ)
    assert flex.mean(diff) > 0.0
    assert flex.max(diff) > 0.0
  else:
    occ_init = xrsp_init.occ
    occ_mod  = xrsp.occ
    assert occ_init.all_eq(1.0)
    assert occ_mod.select(~selection).all_eq(1.0)
    assert occ_mod.select(selection).all_eq(0.75)

def check_remove_selection(
      cmd, xrsp_init, output, selection, selection_str,
      remove_selection_str, verbose, tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  remove_selection = ~xrsp_init.selection(
                                      selection_strings = remove_selection_str)
  assert remove_selection.size() > remove_selection.count(True)
  assert approx_equal(xrsp.sites_cart,
                      xrsp_init.sites_cart.select(remove_selection),tolerance)
  assert approx_equal(
                    xrsp.occ, xrsp_init.occ.select(remove_selection),tolerance)
  assert approx_equal(
                xrsp.u_iso, xrsp_init.u_iso.select(remove_selection),tolerance)
  assert approx_equal(
              xrsp.u_cart, xrsp_init.u_cart.select(remove_selection),tolerance)
  sct1 = xrsp_init.xray_structure.scatterers().extract_scattering_types()
  assert sct1.count("C") > 0
  sct2 = xrsp.xray_structure.scatterers().extract_scattering_types()
  assert sct2.count("C") == 0
  assert sct2.size() == remove_selection.count(True)
  assert sct1.size() > sct2.size()

def check_keep_selection(
      cmd, xrsp_init, output, selection, selection_str,
      keep_selection_str, verbose, tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  keep_selection = xrsp_init.selection(selection_strings = keep_selection_str)
  assert approx_equal(xrsp.sites_cart,
                      xrsp_init.sites_cart.select(keep_selection),tolerance)
  assert approx_equal(xrsp.occ, xrsp_init.occ.select(keep_selection),tolerance)
  assert approx_equal(
    xrsp.u_iso, xrsp_init.u_iso.select(keep_selection),tolerance)
  assert approx_equal(
    xrsp.u_cart, xrsp_init.u_cart.select(keep_selection),tolerance)
  sct1 = xrsp_init.xray_structure.scatterers().extract_scattering_types()
  assert sct1.count("C") > 0 and sct1.size() > sct1.count("C")
  sct2 = xrsp.xray_structure.scatterers().extract_scattering_types()
  assert sct2.count("C") == sct2.size()
  assert sct1.size() > keep_selection.count(True)
  assert sct1.size() > sct2.size()

def check_all_none(cmd, xrsp_init, output, verbose, tolerance=1.e-3):
  remove_files(output)
  run_command(command=cmd, verbose=verbose)
  xrsp = xray_structure_plus(file_name = output)
  assert approx_equal(xrsp.occ,         xrsp_init.occ,tolerance)
  assert approx_equal(xrsp.sites_cart,  xrsp_init.sites_cart,tolerance)
  assert approx_equal(xrsp.use_u_iso,   xrsp_init.use_u_iso,tolerance)
  assert approx_equal(xrsp.use_u_aniso, xrsp_init.use_u_aniso,tolerance)
  assert approx_equal(xrsp.u_iso,       xrsp_init.u_iso,tolerance)
  assert approx_equal(xrsp.u_cart,      xrsp_init.u_cart,tolerance)

def check_keep_remove_conflict(cmd, output, verbose):
  cmd += " keep=all remove=all "
  print(cmd)
  r = easy_run.go(cmd)
  assert r.stdout_lines[-1]==\
    "Sorry: 'keep' and 'remove' keywords cannot be used simultaneously."

def exercise_multiple():
  pdb_str = """\
ATOM      1  O   HOH A   2       5.131   5.251   5.823  1.00 10.00           O
ATOM      2  CA  LYS B  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      3  CB  LYS B  32       9.193   8.732  12.170  1.00 12.23           C
ATOM      4  CA  VAL C  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      5  CB  VAL C  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      6  O   HOH C   3       1.132   5.963   7.065  1.00 15.00           O
ATOM      7  O   HOH C   4       4.132   9.963   7.800  1.00 15.00           O
TER
END
"""
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="exercise_multiple.pdb")
  params = """\
modify{
adp {
  atom_selection = chain A
  randomize = True
}
adp {
  atom_selection = chain B
  shift_b_iso = 10
}
sites {
  atom_selection = chain B
  shake = 1.5
}
sites {
  atom_selection = chain A or chain C
  translate = 1 2 3
  rotate = 4 5 6
}
occupancies
{
  atom_selection = chain A
  randomize = True
}
occupancies
{
  atom_selection = chain C
  set = 0.1
}
}
"""
  with open("params", "w") as f:
    f.write(params)
  cmd = 'phenix.pdbtools exercise_multiple.pdb params'
  print(cmd)
  result = run_command(command=cmd, verbose=False)
  lines = [l.strip() for l in result.stdout_lines]
  expected_lines = [
    "Randomizing ADP: selected atoms: 1 of 7",
    "Adding shift = 10.00 to all ADP: selected atoms: 2 of 7",
    "Shaking sites (RMS = 1.500): selected atoms: 2 of 7",
    "Rigid body shift: selected atoms: 5 of 7",
    "Randomizing occupancies: selected atoms: 1 of 7",
    "Setting occupancies to:    0.100: selected atoms: 4 of 7"
  ]
  for line in expected_lines:
    assert line in lines

def exercise_no_cryst1():
  pdb_str = """
ATOM      1  N  AMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA AMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  AMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  AMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB AMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG AMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD AMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE AMET B  37       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  BMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA BMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  BMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  BMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB BMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG BMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD BMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE BMET B  37       8.775   5.000  10.645  1.00 10.00           C
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="exercise_no_cryst1.pdb")
  cmd = 'phenix.pdbtools exercise_no_cryst1.pdb sites.rotate="0 0 0"'
  print(cmd)
  run_command(command=cmd, verbose=False)
  lines1 = []
  with open("exercise_no_cryst1.pdb","r") as f:
    lines = f.readlines()
  for line in lines:
    line = line.strip()
    assert line.count("CRYST1") == 0
    if(line.startswith("ATOM") or line.startswith("HETATM")):
      lines1.append(line)
  lines2 = []
  with open("exercise_no_cryst1_modified.pdb","r") as f:
    lines = f.readlines()
  for line in lines:
    line = line.strip()
    assert line.count("CRYST1") == 0
    if(line.startswith("ATOM") or line.startswith("HETATM")):
      lines2.append(line)
  assert len(lines1) == len(lines2)
  for l1,l2 in zip(lines1, lines2):
    assert l1[11:70].strip() == l2[11:70].strip()

def exercise_truncate_to_polyala():
  pdb_str = """
ATOM      1  N  AMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA AMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  AMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  AMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB AMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG AMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD AMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE AMET B  37       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  BMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA BMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  BMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  BMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB BMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG BMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD BMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE BMET B  37       8.775   5.000  10.645  1.00 10.00           C
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="exercise_exercise_truncate_to_polyala.pdb")
  cmd = 'phenix.pdbtools exercise_exercise_truncate_to_polyala.pdb truncate_to_polyala=true'
  print(cmd)
  run_command(command=cmd)
  gly_atom_names = [" N  ", " CA ", " C  ", " O  ", " CB "]
  pdb_inp = iotbx.pdb.input(
    file_name="exercise_exercise_truncate_to_polyala_modified.pdb")
  for a in pdb_inp.construct_hierarchy().atoms_with_labels():
    assert a.name in gly_atom_names

def exercise_set_charge():
  input_pdb = """
ATOM      1  CL  CL  X   1       0.000   0.000   0.000  1.00 20.00          CL
END
"""
  with open("tmp_cl.pdb", "w") as f:
    f.write(input_pdb)
  cmd='phenix.pdbtools tmp_cl.pdb charge_selection="element Cl" charge=-1'
  print(cmd)
  run_command(command=cmd, verbose=False)
  pdb_in = iotbx.pdb.input("tmp_cl_modified.pdb")
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  assert (xrs.scatterers()[0].scattering_type == 'Cl1-')
  assert (hierarchy.atoms()[0].charge == '1-')

def exercise_renumber_residues():
  input_pdb = """
ATOM      1  O   GLY A   3       1.434   1.460   2.496  1.00  6.04           O
ATOM      2  O   CYS A   7       2.196   4.467   3.911  1.00  4.51           O
ATOM      3  O   CYS A   1      -1.433   4.734   5.405  1.00  7.85           O
TER
ATOM      4  O   SER B   4       0.297   0.843   7.226  1.00  7.65           O
ATOM      5  OG ASER B   4      -2.625   1.057   4.064  0.50  5.46           O
ATOM      6  OG BSER B   4      -0.885   0.189   3.843  0.50 11.74           O
ATOM      7  O   LEU B   8       3.613   4.307   6.646  1.00  5.39           O
ATOM      8  O   PRO B  -1       4.398   6.723   8.658  1.00  6.65           O
ATOM      9  O   TYR B   7       7.294   7.360   6.923  1.00  8.75           O
ATOM     10  O   CYS B   0       5.256   8.262   4.185  1.00  6.08           O
ATOM     11  O   ALA B   9       3.028  10.447   5.584  1.00  7.39           O
TER
ATOM     12  O   LEU C   0       5.613  12.448   6.864  1.00  7.32           O
TER
END
"""
  expected_output_pdb = """ATOM      1  O   GLY A   1       1.434   1.460   2.496  1.00  6.04           O
ATOM      2  O   CYS A   2       2.196   4.467   3.911  1.00  4.51           O
ATOM      3  O   CYS A   3      -1.433   4.734   5.405  1.00  7.85           O
TER
ATOM      4  O   SER B   1       0.297   0.843   7.226  1.00  7.65           O
ATOM      5  OG ASER B   1      -2.625   1.057   4.064  0.50  5.46           O
ATOM      6  OG BSER B   1      -0.885   0.189   3.843  0.50 11.74           O
ATOM      7  O   LEU B   2       3.613   4.307   6.646  1.00  5.39           O
ATOM      8  O   PRO B   3       4.398   6.723   8.658  1.00  6.65           O
ATOM      9  O   TYR B   4       7.294   7.360   6.923  1.00  8.75           O
ATOM     10  O   CYS B   5       5.256   8.262   4.185  1.00  6.08           O
ATOM     11  O   ALA B   6       3.028  10.447   5.584  1.00  7.39           O
TER
ATOM     12  O   LEU C   1       5.613  12.448   6.864  1.00  7.32           O
TER
"""
  ifn = "exercise_renumber_residues.pdb"
  with open(ifn,"w") as f:
    f.write(input_pdb)
  cmd = 'phenix.pdbtools "%s" renumber_residues=true'%ifn
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("exercise_renumber_residues_modified.pdb") as f:
    lines = f.readlines()
  for line1, line2 in zip(lines, expected_output_pdb.splitlines()):
    line1 = line1.strip()
    line2 = line2.strip()
    assert line1 == line2,[line1, line2]
  # now only a selected chain
  expected_output_pdb_2 = """\
ATOM      1  O   GLY A   3       1.434   1.460   2.496  1.00  6.04           O
ATOM      2  O   CYS A   7       2.196   4.467   3.911  1.00  4.51           O
ATOM      3  O   CYS A   1      -1.433   4.734   5.405  1.00  7.85           O
TER
ATOM      4  O   SER B   1       0.297   0.843   7.226  1.00  7.65           O
ATOM      5  OG ASER B   1      -2.625   1.057   4.064  0.50  5.46           O
ATOM      6  OG BSER B   1      -0.885   0.189   3.843  0.50 11.74           O
ATOM      7  O   LEU B   2       3.613   4.307   6.646  1.00  5.39           O
ATOM      8  O   PRO B   3       4.398   6.723   8.658  1.00  6.65           O
ATOM      9  O   TYR B   4       7.294   7.360   6.923  1.00  8.75           O
ATOM     10  O   CYS B   5       5.256   8.262   4.185  1.00  6.08           O
ATOM     11  O   ALA B   6       3.028  10.447   5.584  1.00  7.39           O
TER
ATOM     12  O   LEU C   0       5.613  12.448   6.864  1.00  7.32           O
TER
"""
  cmd = "phenix.pdbtools \"%s\" renumber_residues=true modify.selection=\"chain B\"" % ifn
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("exercise_renumber_residues_modified.pdb") as f:
    for line1, line2 in zip(f.readlines(), expected_output_pdb_2.splitlines()):
      assert (line1.strip() == line2.strip())
  cmd="phenix.pdbtools \"%s\" increment_resseq=10 modify.selection=\"chain B\"" % ifn
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("exercise_renumber_residues_modified.pdb") as f:
    pdb_new = f.read()
  expected_output_pdb_3 = """\
ATOM      1  O   GLY A   3       1.434   1.460   2.496  1.00  6.04           O
ATOM      2  O   CYS A   7       2.196   4.467   3.911  1.00  4.51           O
ATOM      3  O   CYS A   1      -1.433   4.734   5.405  1.00  7.85           O
TER
ATOM      4  O   SER B  14       0.297   0.843   7.226  1.00  7.65           O
ATOM      5  OG ASER B  14      -2.625   1.057   4.064  0.50  5.46           O
ATOM      6  OG BSER B  14      -0.885   0.189   3.843  0.50 11.74           O
ATOM      7  O   LEU B  18       3.613   4.307   6.646  1.00  5.39           O
ATOM      8  O   PRO B   9       4.398   6.723   8.658  1.00  6.65           O
ATOM      9  O   TYR B  17       7.294   7.360   6.923  1.00  8.75           O
ATOM     10  O   CYS B  10       5.256   8.262   4.185  1.00  6.08           O
ATOM     11  O   ALA B  19       3.028  10.447   5.584  1.00  7.39           O
TER
ATOM     12  O   LEU C   0       5.613  12.448   6.864  1.00  7.32           O
TER
END
"""
  assert not show_diff(pdb_new, expected_output_pdb_3)

def exercise_neutralize_scatterers():
  input_pdb = """ATOM      1  N   TYR    22      -0.813  -2.199   1.423  1.00  0.00           N
ATOM      2  CA  TYR    22       0.612  -2.082   1.127  1.00  0.00           C
ATOM      3  C   TYR    22       1.441  -2.790   2.171  1.00  0.00           C
ATOM      4  O   TYR    22       0.927  -3.375   3.128  1.00  0.00           O
ATOM      5  CB  TYR    22       1.052  -0.589   1.096  1.00  0.00           C
ATOM      6  CG  TYR    22       0.390   0.302   0.038  1.00  0.00           C
ATOM      7  CD1 TYR    22       0.807   0.254  -1.296  1.00  0.00           C
ATOM      8  CD2 TYR    22      -0.647   1.164   0.406  1.00  0.00           C
ATOM      9  CE1 TYR    22       0.188   1.057  -2.250  1.00  0.00           C
ATOM     10  CE2 TYR    22      -1.264   1.966  -0.549  1.00  0.00           C
ATOM     11  CZ  TYR    22      -0.846   1.912  -1.877  1.00  0.00           C
ATOM     12  OH  TYR    22      -1.452   2.696  -2.817  1.00  0.00           O1-
"""
  expected_output_pdb = """ATOM      1  N   TYR    22      -0.813  -2.199   1.423  1.00 10.00           N
ATOM      2  CA  TYR    22       0.612  -2.082   1.127  1.00 10.00           C
ATOM      3  C   TYR    22       1.441  -2.790   2.171  1.00 10.00           C
ATOM      4  O   TYR    22       0.927  -3.375   3.128  1.00 10.00           O
ATOM      5  CB  TYR    22       1.052  -0.589   1.096  1.00 10.00           C
ATOM      6  CG  TYR    22       0.390   0.302   0.038  1.00 10.00           C
ATOM      7  CD1 TYR    22       0.807   0.254  -1.296  1.00 10.00           C
ATOM      8  CD2 TYR    22      -0.647   1.164   0.406  1.00 10.00           C
ATOM      9  CE1 TYR    22       0.188   1.057  -2.250  1.00 10.00           C
ATOM     10  CE2 TYR    22      -1.264   1.966  -0.549  1.00 10.00           C
ATOM     11  CZ  TYR    22      -0.846   1.912  -1.877  1.00 10.00           C
ATOM     12  OH  TYR    22      -1.452   2.696  -2.817  1.00 10.00           O
"""
  ifn = "exercise_neutralize_scatterers.pdb"
  with open(ifn,"w") as f:
    f.write(input_pdb)
  cmd = 'phenix.pdbtools "%s" neutralize_scatterers=true adp.set_b_iso=10'%ifn
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open(ifn.replace(".pdb","_modified.pdb")) as f:
    lines = f.readlines()
  for line1, line2 in zip(lines, expected_output_pdb.splitlines()):
    line1 = line1.strip()
    line2 = line2.strip()
    assert line1 == line2

def exercise_remove_atoms():
  pdb_str = """
ATOM      1  N  AMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA AMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  AMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  AMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB AMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG AMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD AMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE AMET B  37       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  BMET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA BMET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  BMET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  BMET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB BMET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG BMET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD BMET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE BMET B  37       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  AMET B  38       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA AMET B  38       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  AMET B  38       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  AMET B  38       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB AMET B  38       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG AMET B  38       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD AMET B  38       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE AMET B  38       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  BMET B  38       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA BMET B  38       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  BMET B  38       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  BMET B  38       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB BMET B  38       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG BMET B  38       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD BMET B  38       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE BMET B  38       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  AMET B  39       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA AMET B  39       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  AMET B  39       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  AMET B  39       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB AMET B  39       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG AMET B  39       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD AMET B  39       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE AMET B  39       8.775   5.000  10.645  1.00 10.00           C
ATOM      1  N  BMET B  39       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA BMET B  39       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C  BMET B  39       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O  BMET B  39       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB BMET B  39       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG BMET B  39       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD BMET B  39       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE BMET B  39       8.775   5.000  10.645  1.00 10.00           C
TER
END
  """
  # Initial size - 48 atoms
  with open("exercise_remove_atoms.pdb", 'w') as f:
    f.write(pdb_str)
  cmd = " ".join([
    "phenix.pdbtools",
    "exercise_remove_atoms.pdb",
    "remove_fraction=0.1"])
  print(cmd)
  run_command(command=cmd, verbose=False)
  inp = iotbx.pdb.input(file_name="exercise_remove_atoms_modified.pdb")
  # removed 5 atoms (~10%).
  assert inp.atoms().size() == 43, inp.atoms().size()

def exercise_change_of_basis():
  with open("tmp_pdbtools_cb_op.pdb", "w") as f:
    f.write("""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O  """)
  cmd = "phenix.pdbtools tmp_pdbtools_cb_op.pdb change_of_basis='a,c,b'"
  run_command(command=cmd, verbose=False)
  with open("tmp_pdbtools_cb_op_modified.pdb") as f:
    lines = f.readlines()
  for line in lines:
    if line.startswith("CRYST1"):
      assert (line.strip() ==
        """CRYST1   21.937   23.477    4.866  90.00  90.00 107.08 P 1 1 21""")
      break

def exercise_mmcif_support():
  from libtbx.test_utils import open_tmp_file
  f = open_tmp_file(suffix="pdbtools.cif")
  f.write("""\
data_phenix
_space_group.name_H-M_alt         'C 1 2 1'
_space_group.name_Hall            ' C 2y'
_cell.length_a                    46.053
_cell.length_b                    9.561
_cell.length_c                    20.871
_cell.angle_alpha                 90.0
_cell.angle_beta                  97.43
_cell.angle_gamma                 90.0
_cell.volume                      9112.60599144
loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
  ATOM      2  CA  .  LYS  A  1  ?    7.49733  -0.62028   4.35289  1.000  10.25989  C  ?  A  ?   1  1
  ATOM     11  CA  .  LEU  A  2  ?    3.72032  -0.19320   3.89326  1.000   7.80433  C  ?  A  ?   2  1
  ATOM     19  CA  .  VAL  A  3  ?    0.78668  -0.39555   6.35234  1.000   5.03864  C  ?  A  ?   3  1
  ATOM     26  CA  .  PHE  A  4  ?   -2.75438  -0.21383   5.02429  1.000   8.93080  C  ?  A  ?   4  1
  ATOM     37  CA  .  PHE  A  5  ?   -6.05155  -0.46197   6.85390  1.000   9.57417  C  ?  A  ?   5  1
  ATOM     48  CA  .  ALA  A  6  ?   -9.57646  -0.10942   5.55847  1.000  17.73488  C  ?  A  ?   6  1
  ATOM     54  CA  .  LYS  B  1  ?   -8.86604  -5.20044   5.46515  1.000  16.15297  C  ?  B  ?   7  1
""")
  f.close()
  cmd = " ".join(["phenix.pdbtools", "\"%s\"" % f.name,
    "rename_chain_id.old_id=A",
    "rename_chain_id.new_id=C"])
  print(cmd)
  run_command(command=cmd, verbose=False)
  ofn = f.name[:].replace(".cif","_modified.cif")
  assert os.path.isfile(ofn)
  pdb_inp = iotbx.pdb.input(file_name=ofn)
  assert pdb_inp.file_type() == "mmcif"
  hierarchy = pdb_inp.construct_hierarchy()
  assert [chain.id for chain in hierarchy.chains()] == ['C', 'B']

  # check target_output_format=pdb
  cmd = " ".join(["phenix.pdbtools", "\"%s\"" % f.name,
    "rename_chain_id.old_id=A",
    "rename_chain_id.new_id=C",
    "target_output_format=pdb",])
  print(cmd)
  run_command(command=cmd, verbose=False)
  ofn = f.name[:].replace(".cif","_modified.pdb")
  assert os.path.isfile(ofn)

def exercise_mmcif_support_2(prefix="tst_pdbtools_mmcif2"):
  f = open("%s.pdb" % prefix, 'w')
  f.write("""\
HELIX    1   1 GLN A   10  GLN A   14  5                                   5
HELIX    2   2 ARG A   17  ASN A   29  1                                  13
SHEET    1   A 2 GLN A  33  GLN A  38  0
SHEET    2   A 2 GLN A  51  GLN A  56  1  O  VAL A  54   N  ILE A  35
SSBOND   1 CYS A    4    CYS A   49
CRYST1   62.654   62.654   45.906  90.00  90.00  90.00 P 43 21 2
SCALE1      0.015961  0.000000  0.000000        0.00000
SCALE2      0.000000  0.015961  0.000000        0.00000
SCALE3      0.000000  0.000000  0.021784        0.00000
ATOM      1  N   GLN A   3      23.762  12.429   1.265  1.00 45.30           N
ATOM      2  CA  GLN A   3      22.670  12.004   2.148  1.00 39.46           C
ATOM      3  C   GLN A   3      21.818  10.942   1.458  1.00 42.86           C
ATOM      4  O   GLN A   3      21.337  11.133   0.339  1.00 44.81           O
ATOM      5  CB  GLN A   3      21.794  13.195   2.571  1.00 42.53           C
ATOM      6  N   CYS A   4      21.620   9.817   2.129  1.00 36.33           N
ATOM      7  CA  CYS A   4      20.902   8.709   1.518  1.00 31.46           C
ATOM      8  C   CYS A   4      19.419   8.891   1.670  1.00 34.16           C
ATOM      9  O   CYS A   4      18.979   9.537   2.607  1.00 38.12           O
ATOM     10  CB  CYS A   4      21.307   7.409   2.205  1.00 31.53           C
ATOM     11  SG  CYS A   4      23.070   7.194   2.148  1.00 29.88           S
ATOM     12  N   SER A   5      18.657   8.290   0.760  1.00 33.83           N
ATOM     13  CA  SER A   5      17.211   8.293   0.854  1.00 38.94           C
ATOM     14  C   SER A   5      16.676   6.917   1.254  1.00 32.88           C
ATOM     15  O   SER A   5      17.295   5.879   0.981  1.00 32.57           O
ATOM     16  CB  SER A   5      16.603   8.683  -0.488  1.00 39.48           C
ATOM     17  OG  SER A   5      16.893   7.668  -1.432  1.00 47.16           O
ATOM     18  N   GLY A   6      15.520   6.924   1.901  1.00 31.69           N
ATOM     19  CA  GLY A   6      14.811   5.694   2.191  1.00 32.41           C
ATOM     20  C   GLY A   6      15.257   5.112   3.511  1.00 29.00           C
ATOM      2  CA  GLN A  10      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  11      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  12      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  13      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  14      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  33      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  34      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  35      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  36      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  37      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  38      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  51      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  52      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  53      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  54      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  55      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  56      22.670  12.004   2.148  1.00 39.46           C
""")
  f.close()
  cmd = " ".join([
      "phenix.pdb_as_cif",
      "%s.pdb" % prefix])
  run_command(command=cmd, verbose=False)
  cmd = " ".join([
      "phenix.pdbtools",
      "%s.cif" % prefix,
      "suffix=none",
      "output.prefix=%s" % prefix])
  print(cmd)
  run_command(command=cmd, verbose=False)
  cif_f = open("%s.cif" % prefix, 'r')
  cif_l = cif_f.readlines()
  cif_f.close()
  for l in ["_cell.angle_alpha                 90.000\n",
      "  _struct_conf.pdbx_PDB_helix_id\n",
      "  _struct_sheet.id\n",
      "  _struct_sheet_range.id\n",
      "  _atom_site.label_atom_id\n"]:
    assert l in cif_l, "%s not in cif file!" % l

def exercise_move_waters():
  pdb_in = """\
ATOM     16  O  AHOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM     60  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM     63  CB ALYS A  32       9.197   8.686  12.246  0.29 14.71           C
ATOM     64  CB BLYS A  32       9.193   8.732  12.170  0.71 12.23           C
ATOM     74  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM     77  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM     18  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM     19  H1  HOH A   3       1.160   5.211   6.437  1.00 15.00           H
ATOM     20  H2  HOH A   3       1.122   5.579   7.967  1.00 15.00           H
HETATM 2397  P   PO4     1      -7.520  25.376  38.369  1.00 39.37           P
HETATM 2398  O1  PO4     1      -6.610  24.262  38.967  1.00 40.00           O
HETATM 2399  O2  PO4     1      -6.901  25.919  37.049  1.00 41.07           O
HETATM 2400  O3  PO4     1      -8.894  24.741  38.097  1.00 45.09           O
HETATM 2401  O4  PO4     1      -7.722  26.556  39.350  1.00 42.48           O
ATOM     23 CL   CL  B   1       6.302   6.419   1.560  0.50 10.00          CL
HETATM 6362  O   HOH B   2      47.616  10.724 150.212  1.00 46.48       B   O
HETATM 6363  O  AHOH B   3      46.408  16.672 146.066  0.50 12.81       B   O
HETATM 6364  O   HOH B   4      29.343  12.806 185.898  1.00 35.57       B   O
HETATM 6365  O  BHOH B   5      43.786  12.615 147.734  0.50 28.43       B   O
HETATM 6366  O   HOH B   6      35.068  19.167 155.349  1.00 15.97       B   O
"""
  pdb_out = """\
ATOM      1  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      2  CB ALYS A  32       9.197   8.686  12.246  0.29 14.71           C
ATOM      3  CB BLYS A  32       9.193   8.732  12.170  0.71 12.23           C
ATOM      4  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      5  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
TER
HETATM    6  O1  PO4     1      -6.610  24.262  38.967  1.00 40.00           O
HETATM    7  O2  PO4     1      -6.901  25.919  37.049  1.00 41.07           O
HETATM    8  O3  PO4     1      -8.894  24.741  38.097  1.00 45.09           O
HETATM    9  O4  PO4     1      -7.722  26.556  39.350  1.00 42.48           O
HETATM   10  P   PO4     1      -7.520  25.376  38.369  1.00 39.37           P
ATOM     11 CL   CL  B   1       6.302   6.419   1.560  0.50 10.00          Cl
ATOM     12  O  AHOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM     13  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM     14  H1  HOH A   3       1.160   5.211   6.437  1.00 15.00           H
ATOM     15  H2  HOH A   3       1.122   5.579   7.967  1.00 15.00           H
HETATM   16  O   HOH B   2      47.616  10.724 150.212  1.00 46.48       B   O
HETATM   17  O  AHOH B   3      46.408  16.672 146.066  0.50 12.81       B   O
HETATM   18  O   HOH B   4      29.343  12.806 185.898  1.00 35.57       B   O
HETATM   19  O  BHOH B   5      43.786  12.615 147.734  0.50 28.43       B   O
HETATM   20  O   HOH B   6      35.068  19.167 155.349  1.00 15.97       B   O
END
"""
  with open("tst_pdbtools_move_waters.pdb", "w") as f:
    f.write(pdb_in)
  cmd = "phenix.pdbtools tst_pdbtools_move_waters.pdb move_waters_last=True"
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("tst_pdbtools_move_waters_modified.pdb") as f:
    pdb_new = f.read()
  assert pdb_new == pdb_out, pdb_new

def exercise_stop_for_unknowns():
  pdb_in = """\
HETATM   16  O   UNK B   2      47.616  10.724 150.212  1.00 46.48       B   O
HETATM   17  O   UNK B   3      46.408  16.672 146.066  0.50 12.81       B   O
HETATM   18  O   UNK B   4      29.343  12.806 185.898  1.00 35.57       B   O
HETATM   19  O   UNK B   5      43.786  12.615 147.734  0.50 28.43       B   O
HETATM   20  O   UNK B   6      35.068  19.167 155.349  1.00 15.97       B   O
"""
  with open("tst_pdbtools_unknown.pdb", "w") as f:
    f.write(pdb_in)
  cmd = "phenix.pdbtools tst_pdbtools_unknown.pdb set_b_iso=20"
  print(cmd)
  run_command(command=cmd, sorry_expected=True)
  cmd2 = cmd + " stop_for_unknowns=False"
  print(cmd2)
  run_command(command=cmd2)

def exercise_average_alt_confs(prefix='tst_pdbtools_average_alt_confs'):
  pdb_in = """\
ATOM     16  O  AHOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM     60  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM     63  CB ALYS A  32       9.197   8.686  12.246  0.29 14.71           C
ATOM     64  CB BLYS A  32       9.193   8.732  12.170  0.71 12.23           C
ATOM     74  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM     77  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM     18  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM     19  O  BHOH A   4       4.132   9.963   7.800  0.50 15.00           O
"""
  with open("%s.pdb" % prefix, "w") as f:
    f.write(pdb_in)
  cmd = "phenix.pdbtools %s.pdb average_alt_confs=True" % prefix
  print(cmd)
  run_command(command=cmd, verbose=False)

  with open("%s_modified.pdb" % prefix) as f:
    pdb_new = f.read()
  assert (pdb_new == '''\
ATOM      1  O  AHOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM      2  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      3  CB ALYS A  32       9.195   8.709  12.208  0.29 14.71           C
ATOM      4  CB BLYS A  32       9.195   8.709  12.208  0.71 12.23           C
ATOM      5  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      6  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      7  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM      8  O  BHOH A   4       4.132   9.963   7.800  0.50 15.00           O
TER
END
''')

def exercise_remove_alt_confs():
  pdb_in = """\
ATOM     16  O  AHOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM     60  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM     63  CB ALYS A  32       9.197   8.686  12.246  0.29 14.71           C
ATOM     64  CB BLYS A  32       9.193   8.732  12.170  0.71 12.23           C
ATOM     74  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM     77  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM     18  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM     19  O  BHOH A   4       4.132   9.963   7.800  0.50 15.00           O
"""
  with open("tst_pdbtools_alt_confs.pdb", "w") as f:
    f.write(pdb_in)
  cmd = "phenix.pdbtools tst_pdbtools_alt_confs.pdb remove_alt_confs=True"
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("tst_pdbtools_alt_confs_modified.pdb") as f:
    pdb_new = f.read()
  assert (pdb_new == """\
ATOM      1  O   HOH A   2       5.131   5.251   5.823  1.00 10.00           O
ATOM      2  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      3  CB  LYS A  32       9.197   8.686  12.246  1.00 14.71           C
ATOM      4  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      5  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      6  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
TER
END
""")
  cmd = "phenix.pdbtools tst_pdbtools_alt_confs.pdb remove_alt_confs=True " +\
    "always_keep_one_conformer=True"
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("tst_pdbtools_alt_confs_modified.pdb") as f:
    pdb_new = f.read()
  assert (pdb_new == """\
ATOM      1  O   HOH A   2       5.131   5.251   5.823  1.00 10.00           O
ATOM      2  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      3  CB  LYS A  32       9.193   8.732  12.170  1.00 12.23           C
ATOM      4  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      5  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      6  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM      7  O   HOH A   4       4.132   9.963   7.800  1.00 15.00           O
TER
END
""")
  cmd = "phenix.pdbtools tst_pdbtools_alt_confs.pdb remove_alt_confs=True altloc_to_keep='B' "
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("tst_pdbtools_alt_confs_modified.pdb") as f:
    pdb_new = f.read()

  assert (pdb_new == """\
ATOM      1  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      2  CB  LYS A  32       9.193   8.732  12.170  1.00 12.23           C
ATOM      3  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      4  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      5  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM      6  O   HOH A   4       4.132   9.963   7.800  1.00 15.00           O
TER
END
""")
  cmd = "phenix.pdbtools tst_pdbtools_alt_confs.pdb remove_alt_confs=True " +\
    "always_keep_one_conformer=True altloc_to_keep='B' "
  print(cmd)
  run_command(command=cmd, verbose=False)
  with open("tst_pdbtools_alt_confs_modified.pdb") as f:
    pdb_new = f.read()
  assert (pdb_new == """\
ATOM      1  O   HOH A   2       5.131   5.251   5.823  1.00 10.00           O
ATOM      2  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      3  CB  LYS A  32       9.193   8.732  12.170  1.00 12.23           C
ATOM      4  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      5  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
ATOM      6  O   HOH A   3       1.132   5.963   7.065  1.00 15.00           O
ATOM      7  O   HOH A   4       4.132   9.963   7.800  1.00 15.00           O
TER
END
""")


def exercise_convert_met_to_semet():
  pdb_str_met = """
ATOM      1  N   MET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA  MET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C   MET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O   MET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB  MET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG  MET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD  MET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE  MET B  37       8.775   5.000  10.645  1.00 10.00           C
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str_met)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="exercise_convert_met_to_semet.pdb")
  cmd = " ".join([
    "phenix.pdbtools",
    "exercise_convert_met_to_semet.pdb",
    "convert_met_to_semet=true"])
  print(cmd)
  run_command(command=cmd, verbose=False)
  pi_out = iotbx.pdb.input(
    file_name="exercise_convert_met_to_semet_modified.pdb"
    ).construct_hierarchy()
  for rg in pi_out.residue_groups():
      for rn in rg.unique_resnames():
        assert rn=="MSE"
  cmd = " ".join([
    "phenix.pdbtools",
    "exercise_convert_met_to_semet_modified.pdb",
    "convert_semet_to_met=true"])
  run_command(command=cmd, verbose=False)
  pi_out = iotbx.pdb.input(
    file_name="exercise_convert_met_to_semet_modified_modified.pdb"
    ).construct_hierarchy()
  for rg in pi_out.residue_groups():
      for rn in rg.unique_resnames():
        assert rn=="MET"

def exercise_switch_rotamers(prefix="exercise_switch_rotamers"):
  pdb_str_met = """
ATOM      1  N   MET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA  MET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C   MET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O   MET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB  MET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG  MET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD  MET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE  MET B  37       8.775   5.000  10.645  1.00 10.00           C
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str_met)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="%s.pdb"%prefix)
  for o in ["max_distant","min_distant","exact_match","fix_outliers"]:
    cmd = " ".join([
      "phenix.pdbtools",
      "%s.pdb"%prefix,
      "switch_rotamers=%s"%o,
      "suffix=none",
      "output.prefix=%s_%s"%(o,prefix)])
    print(cmd)
    run_command(command=cmd, verbose=False)

def exercise_segid_manipulations(prefix="tst_pdbtools_ex_segid"):
  pdb_str_met = """
ATOM      1  N   MET B  37       7.525   5.296   6.399  1.00 10.00      A    N
ATOM      2  CA  MET B  37       6.533   6.338   6.634  1.00 10.00      A    C
ATOM      3  C   MET B  37       6.175   7.044   5.330  1.00 10.00      A    C
ATOM      4  O   MET B  37       5.000   7.200   5.000  1.00 10.00      A    O
ATOM      5  CB  MET B  37       7.051   7.351   7.655  1.00 10.00      A    C
ATOM      6  CG  MET B  37       7.377   6.750   9.013  1.00 10.00      A    C
ATOM      7  SD  MET B  37       8.647   5.473   8.922  1.00 10.00      A    S
ATOM      8  CE  MET B  37       8.775   5.000  10.645  1.00 10.00      A    C
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str_met)
  ph_met_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="%s.pdb"%prefix)
  for o  in ["clear_seg_id", "set_seg_id_to_chain_id"]:
    cmd = " ".join([
        "phenix.pdbtools",
        "%s.pdb"%prefix,
        "%s=True" % o,
        "suffix=none",
        "output.prefix=%s_%s"%(o,prefix)])
    print(cmd)
    run_command(command=cmd, verbose=False)
    with open("%s_%s.pdb"%(o,prefix)) as f:
      pdb_new = f.read()
    if o == "clear_seg_id":
      assert (pdb_new == """\
ATOM      1  N   MET B  37       7.525   5.296   6.399  1.00 10.00           N
ATOM      2  CA  MET B  37       6.533   6.338   6.634  1.00 10.00           C
ATOM      3  C   MET B  37       6.175   7.044   5.330  1.00 10.00           C
ATOM      4  O   MET B  37       5.000   7.200   5.000  1.00 10.00           O
ATOM      5  CB  MET B  37       7.051   7.351   7.655  1.00 10.00           C
ATOM      6  CG  MET B  37       7.377   6.750   9.013  1.00 10.00           C
ATOM      7  SD  MET B  37       8.647   5.473   8.922  1.00 10.00           S
ATOM      8  CE  MET B  37       8.775   5.000  10.645  1.00 10.00           C
TER
END
"""), pdb_new
    else:
      assert (pdb_new == """\
ATOM      1  N   MET B  37       7.525   5.296   6.399  1.00 10.00      B    N
ATOM      2  CA  MET B  37       6.533   6.338   6.634  1.00 10.00      B    C
ATOM      3  C   MET B  37       6.175   7.044   5.330  1.00 10.00      B    C
ATOM      4  O   MET B  37       5.000   7.200   5.000  1.00 10.00      B    O
ATOM      5  CB  MET B  37       7.051   7.351   7.655  1.00 10.00      B    C
ATOM      6  CG  MET B  37       7.377   6.750   9.013  1.00 10.00      B    C
ATOM      7  SD  MET B  37       8.647   5.473   8.922  1.00 10.00      B    S
ATOM      8  CE  MET B  37       8.775   5.000  10.645  1.00 10.00      B    C
TER
END
"""), pdb_new

def exercise_flip_symmetric_amino_acids(
  prefix="exercise_flip_symmetric_amino_acids"):
  pdb_str_in = """
ATOM   3356  N   ARG H 228      20.811-103.620  24.514  1.00 56.77           N
ATOM   3357  CA  ARG H 228      21.863-104.460  23.967  1.00 64.89           C
ATOM   3358  C   ARG H 228      21.258-105.679  23.284  1.00 69.98           C
ATOM   3359  O   ARG H 228      20.179-106.141  23.657  1.00 70.23           O
ATOM   3360  CB  ARG H 228      22.810-104.925  25.075  1.00 63.55           C
ATOM   3361  CG  ARG H 228      23.354-103.810  25.936  1.00 64.01           C
ATOM   3362  CD  ARG H 228      24.570-104.271  26.712  1.00 65.11           C
ATOM   3363  NE  ARG H 228      25.126-103.194  27.523  1.00 68.09           N
ATOM   3364  CZ  ARG H 228      26.373-103.169  27.982  1.00 66.49           C
ATOM   3365  NH1 ARG H 228      26.790-102.145  28.717  1.00 64.08           N
ATOM   3366  NH2 ARG H 228      27.205-104.163  27.694  1.00 65.11           N
TER
END
  """
  pdb_str_out = """
ATOM   3355  N   ARG H 228      20.811-103.620  24.514  1.00 56.77           N
ATOM   3356  CA  ARG H 228      21.863-104.460  23.967  1.00 64.89           C
ATOM   3357  C   ARG H 228      21.258-105.679  23.284  1.00 69.98           C
ATOM   3358  O   ARG H 228      20.179-106.141  23.657  1.00 70.23           O
ATOM   3359  CB  ARG H 228      22.810-104.925  25.075  1.00 63.55           C
ATOM   3360  CG  ARG H 228      23.354-103.810  25.936  1.00 64.01           C
ATOM   3361  CD  ARG H 228      24.570-104.271  26.712  1.00 65.11           C
ATOM   3362  NE  ARG H 228      25.126-103.194  27.523  1.00 68.09           N
ATOM   3363  CZ  ARG H 228      26.373-103.169  27.982  1.00 66.49           C
ATOM   3364  NH1 ARG H 228      27.205-104.163  27.694  1.00 65.11           N
ATOM   3365  NH2 ARG H 228      26.790-102.145  28.717  1.00 64.08           N
TER
END
  """
  pi = iotbx.pdb.input(source_info=None, lines=pdb_str_in)
  ph_in = pi.construct_hierarchy()
  pi.write_pdb_file(file_name="%s.pdb"%prefix)
  cmd = " ".join([
    "phenix.pdbtools",
    "%s.pdb"%prefix,
    "flip_symmetric_amino_acids=True",
    "output.prefix=%s"%prefix])
  print(cmd)
  run_command(command=cmd, verbose=False)
  h1 = iotbx.pdb.input(source_info=None, lines=pdb_str_out).construct_hierarchy()
  h2 = iotbx.pdb.input(file_name=
    "exercise_flip_symmetric_amino_acids_modified.pdb").construct_hierarchy()
  assert h1.is_similar_hierarchy(h2)
  ph_in.flip_symmetric_amino_acids()
  assert h1.is_similar_hierarchy(ph_in)

def exercise_result_is_empty(prefix='exercise_result_is_empty'):
  pdb_str = """
HELIX    1   1 GLN A   10  GLN A   14  5                                   5
HELIX    2   2 ARG A   17  ASN A   29  1                                  13
SHEET    1   A 2 GLN A  33  GLN A  38  0
SHEET    2   A 2 GLN A  51  GLN A  56  1  O  VAL A  54   N  ILE A  35
SSBOND   1 CYS A    4    CYS A   49
CRYST1   62.654   62.654   45.906  90.00  90.00  90.00 P 43 21 2
SCALE1      0.015961  0.000000  0.000000        0.00000
SCALE2      0.000000  0.015961  0.000000        0.00000
SCALE3      0.000000  0.000000  0.021784        0.00000
ATOM      1  N   GLN A   3      23.762  12.429   1.265  1.00 45.30           N
ATOM      2  CA  GLN A   3      22.670  12.004   2.148  1.00 39.46           C
ATOM      3  C   GLN A   3      21.818  10.942   1.458  1.00 42.86           C
ATOM      4  O   GLN A   3      21.337  11.133   0.339  1.00 44.81           O
ATOM      5  CB  GLN A   3      21.794  13.195   2.571  1.00 42.53           C
ATOM      6  N   CYS A   4      21.620   9.817   2.129  1.00 36.33           N
ATOM      7  CA  CYS A   4      20.902   8.709   1.518  1.00 31.46           C
ATOM      8  C   CYS A   4      19.419   8.891   1.670  1.00 34.16           C
ATOM      9  O   CYS A   4      18.979   9.537   2.607  1.00 38.12           O
ATOM     10  CB  CYS A   4      21.307   7.409   2.205  1.00 31.53           C
ATOM     11  SG  CYS A   4      23.070   7.194   2.148  1.00 29.88           S
ATOM     12  N   SER A   5      18.657   8.290   0.760  1.00 33.83           N
ATOM     13  CA  SER A   5      17.211   8.293   0.854  1.00 38.94           C
ATOM     14  C   SER A   5      16.676   6.917   1.254  1.00 32.88           C
ATOM     15  O   SER A   5      17.295   5.879   0.981  1.00 32.57           O
ATOM     16  CB  SER A   5      16.603   8.683  -0.488  1.00 39.48           C
ATOM     17  OG  SER A   5      16.893   7.668  -1.432  1.00 47.16           O
ATOM     18  N   GLY A   6      15.520   6.924   1.901  1.00 31.69           N
ATOM     19  CA  GLY A   6      14.811   5.694   2.191  1.00 32.41           C
ATOM     20  C   GLY A   6      15.257   5.112   3.511  1.00 29.00           C
ATOM      2  CA  GLN A  10      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  11      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  12      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  13      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  14      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  33      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  34      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  35      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  36      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  37      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  38      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  51      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  52      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  53      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  54      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  55      22.670  12.004   2.148  1.00 39.46           C
ATOM      2  CA  GLN A  56      22.670  12.004   2.148  1.00 39.46           C
"""

  with open("%s.pdb"%prefix, 'w') as f:
    f.write(pdb_str)
  cmd = " ".join([
      "phenix.pdbtools",
      "%s.pdb"%prefix,
      "keep='resname ALA'",
      "output.prefix=out_%s"%(prefix)])
  print(cmd)
  run_command(command=cmd, verbose=False)
  # assert os.path.isfile("out_%s_modified.pdb"%(prefix))
  with open("out_%s_modified.pdb"%(prefix), 'r') as f:
    pdb_new = f.read()
    assert pdb_new == """\
CRYST1   62.654   62.654   45.906  90.00  90.00  90.00 P 43 21 2
SCALE1      0.015961  0.000000  0.000000        0.00000
SCALE2      0.000000  0.015961  0.000000        0.00000
SCALE3      0.000000  0.000000  0.021784        0.00000
END
"""

def exercise_renumber_and_move_waters(prefix='exercise_renumber_and_move_waters'):
  pdb_str = """
ATOM     16  O   HOH A   2       5.131   5.251   5.823  0.60 10.00           O
ATOM     18  O   HOH AAFF4       1.132   5.963   7.065  1.00 15.00           O
ATOM     19  O   HOH A   4       4.132   9.963   7.800  0.50 15.00           O
TER
ATOM     60  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM     64  CB  LYS A  32       9.193   8.732  12.170  1.00 12.23           C
ATOM     74  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM     77  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
"""

  with open("%s.pdb"%prefix, 'w') as f:
    f.write(pdb_str)
  cmd = " ".join([
      "phenix.pdbtools",
      "%s.pdb"%prefix,
      "renumber_and_move_waters=True",
      "output.prefix=out_%s"%(prefix)])
  print(cmd)
  run_command(command=cmd, verbose=False)
  # assert os.path.isfile("out_%s_modified.pdb"%(prefix))
  with open("out_%s_modified.pdb"%(prefix), 'r') as f:
    pdb_new = f.read()
    assert pdb_new == """\
ATOM      1  CA  LYS A  32      10.574   8.177  11.768  1.00 11.49           C
ATOM      2  CB  LYS A  32       9.193   8.732  12.170  1.00 12.23           C
ATOM      3  CA  VAL A  33      11.708   5.617  14.332  1.00 11.42           C
ATOM      4  CB  VAL A  33      11.101   4.227  14.591  1.00 11.47           C
TER
ATOM      5  O   HOH A 100       5.131   5.251   5.823  0.60 10.00           O
ATOM      6  O   HOH A 101       1.132   5.963   7.065  1.00 15.00           O
ATOM      7  O   HOH A 102       4.132   9.963   7.800  0.50 15.00           O
END
"""

def exercise(args):
  exercise_average_alt_confs()
  exercise_flip_symmetric_amino_acids()
  exercise_switch_rotamers()
  exercise_mmcif_support_2()
  exercise_basic()
  exercise_multiple()
  exercise_no_cryst1()
  exercise_renumber_residues()
  exercise_change_of_basis()
  exercise_move_waters()
  exercise_remove_alt_confs()
  exercise_truncate_to_polyala()
  exercise_convert_met_to_semet()
  exercise_set_charge()
  exercise_neutralize_scatterers()
  exercise_remove_atoms()
  exercise_mmcif_support()
  exercise_segid_manipulations()
  exercise_result_is_empty()
  exercise_renumber_and_move_waters()
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
