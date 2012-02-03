from libtbx import easy_run
import libtbx.load_env
import sys, os, math
from cctbx.array_family import flex
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import remove_files, search_for
from mmtbx import utils
from libtbx.test_utils import approx_equal, not_approx_equal, run_command, \
  show_diff
from libtbx.utils import show_times_at_exit
import iotbx.pdb
from scitbx.array_family import flex
from cctbx import adptbx

class xray_structure_plus(object):
  def __init__(self, file_name):
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
                                                   mon_lib_srv = mon_lib_srv,
                                                   ener_lib    = ener_lib,
                                                   file_name   = file_name)
    self.xray_structure = processed_pdb_file.xray_structure(show_summary=False)
    self.all_chain_proxies = processed_pdb_file.all_chain_proxies
    uc = self.xray_structure.unit_cell()
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
                                     all_chain_proxies= self.all_chain_proxies,
                                     selection_strings= selection_strings,
                                     xray_structure   = self.xray_structure)[0]

def exercise_basic(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  output = "modified.pdb"
  xrsp_init = xray_structure_plus(file_name = file_name)
  assert file_name.find('"') < 0
  base = \
      'phenix.pdbtools "%s" modify.output.file_name=%s '%(file_name, output)
  for selection_str in [None, "chain A or chain C"]:
    selection = xrsp_init.selection(selection_strings = selection_str)
    if(selection_str is None):
      assert selection.size() == selection.count(True)
    else:
      assert selection.size() == 36 and selection.count(True) == 24
    #
    cmd = base + 'adp.randomize=true selection="%s"'%str(selection_str)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.set_b_iso=10.0 selection="%s"'%str(selection_str)
    check_adp_set_b_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.shift_b_iso=20.0 selection="%s"'%str(selection_str)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.scale_adp=2.0 selection="%s"'%str(selection_str)
    check_adp_rand(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.convert_to_iso=true selection="%s"'%str(selection_str)
    check_adp_to_iso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base + 'adp.convert_to_aniso=true selection="%s"'%str(selection_str)
    check_adp_to_aniso(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    shake = 1.5
    cmd = base+'sites.shake=%s selection="%s"'%(str(shake), str(selection_str))
    check_sites_shake(
      cmd, xrsp_init, output, selection, selection_str, shake, verbose)
    #
    cmd = base+'sites.rotate="1 2 3" sites.translate="4 5 6" selection="%s"'%(
                                                            str(selection_str))
    check_sites_rt(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base+'occupancies.randomize=true selection="%s"'%(str(selection_str))
    check_occ_randomize(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    cmd = base+'occupancies.set=0.75 selection="%s"'%(str(selection_str))
    check_occ_set(
      cmd, xrsp_init, output, selection, selection_str, verbose)
    #
    remove_selection_str = "element C"
    cmd = base+'remove="%s" selection="%s"'%(
      str(remove_selection_str), str(selection_str))
    check_remove_selection(
      cmd, xrsp_init, output, selection, selection_str,
      remove_selection_str, verbose)
    #
    keep_selection_str = "element C"
    cmd = base+'keep="%s" selection="%s"'%(
      str(keep_selection_str), str(selection_str))
    check_keep_selection(
      cmd, xrsp_init, output, selection, selection_str,
      keep_selection_str, verbose)
    #
    test_quiet(file_name, verbose)
  #
  cmd = base
  check_all_none(cmd, xrsp_init, output, verbose)
  #
  cmd = base
  check_keep_remove_conflict(cmd, output, verbose)

def test_quiet(file_name, verbose):
  output_file_name = "shifted.pdb"
  remove_files("log")
  remove_files(output_file_name)
  assert file_name.find('"') < 0
  cmd= 'phenix.pdbtools "%s" output.file_name=%s shake=0.1 --quiet > log'%(
    file_name, output_file_name)
  print cmd
  run_command(command=cmd, verbose=verbose)
  lines = open("log","r").readlines()
  assert len(lines) == 0

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
      tolerance=1.e-3):
  remove_files(output)
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
  cmd_result = run_command(command=cmd, verbose=verbose, sorry_expected=True)
  sorry_lines = search_for(
    pattern="Sorry: Ambiguous selection:"
            " 'keep' and 'remove' keywords cannot be used simultaneously:",
    mode="==",
    lines=cmd_result.stdout_lines)
  assert len(sorry_lines) == 1

def exercise_multiple(pdb_dir, verbose):
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
}
"""
  open("params", "w").write(params)
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" modify.output.file_name=modified.pdb params' % (
    file_name)
  result = run_command(command=cmd, verbose=verbose)
  lines = result.stdout_lines
  for i,line in enumerate(lines):
    if (line.find("Performing requested model manipulations") >= 0):
      break
  else:
    raise RuntimeError("Expected output not found.")
  assert lines[i+1] == ""
  assert lines[i+6] == ""
  assert not show_diff("\n".join(lines[i+2:i+6]), """\
Randomizing ADP: selected atoms: 12 of 36
Adding shift = 10.00 to all ADP: selected atoms: 12 of 36
Shaking sites (RMS = 1.500): selected atoms: 12 of 36
Rigid body shift: selected atoms: 24 of 36""")

def exercise_no_cryst1(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "t.pdb")
  output = "modified.pdb"
  assert file_name.find('"') < 0
  base = 'phenix.pdbtools "%s" modify.output.file_name=%s '%(file_name, output)
  cmd = base+'sites.rotate="0 0 0" sites.translate="0 0 0"'
  run_command(command=cmd, verbose=verbose)
  lines1 = []
  for line in open(file_name,"r").readlines():
    line = line.strip()
    assert line.count("CRYST1") == 0
    if(line.startswith("ATOM") or line.startswith("HETATM")):
      lines1.append(line)
  lines2 = []
  for line in open(output,"r").readlines():
    line = line.strip()
    assert line.count("CRYST1") == 0
    if(line.startswith("ATOM") or line.startswith("HETATM")):
      lines2.append(line)
  assert len(lines1) == len(lines2)
  for l1,l2 in zip(lines1, lines2):
    assert l1[11:70].strip() == l2[11:70].strip()

def exercise_show_adp_statistics(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "t.pdb")
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" --show-adp-statistics'%file_name
  run_command(command=cmd, verbose=verbose)

def exercise_show_geometry_statistics(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" --show-geometry-statistics'%file_name
  run_command(command=cmd, verbose=verbose)

def exercise_show_number_of_removed(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_h.pdb")
  log = "exercise_show_number_of_removed.log"
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" remove="element H" > %s' % (file_name, log)
  remove_files(log)
  run_command(command=cmd, verbose=verbose)
  assert os.path.isfile(log)
  n_found = 0
  for line in open(log).read().splitlines():
    if (line == "Atoms to be kept: 13 of 24"):
      n_found += 1
  assert n_found == 1

def exercise_02(pdb_dir, verbose):
  from iotbx import file_reader
  file_name = os.path.join(pdb_dir, "polypro_Simon_noCRYST1.pdb")
  log = "exercise_02.log"
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" --geometry-regularization > %s' % (file_name, log)
  remove_files(log)
  run_command(command=cmd, verbose=verbose)
  pdb_1 = file_reader.any_file("polypro_Simon_noCRYST1.pdb_modified.pdb")
  assert os.path.isfile(log)
  log = "exercise_02b.log"
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" regularize_geometry=True > %s' % (file_name, log)
  remove_files(log)
  remove_files("polypro_Simon_noCRYST1.pdb_modified.pdb")
  run_command(command=cmd, verbose=verbose)
  pdb_2 = file_reader.any_file("polypro_Simon_noCRYST1.pdb_modified.pdb")
  sites_1 = pdb_1.file_object.atoms().extract_xyz()
  sites_2 = pdb_2.file_object.atoms().extract_xyz()
  assert (sites_1.size() == sites_2.size())
  assert (sites_1.rms_difference(sites_2) < 0.00001)

def exercise_truncate_to_polyala(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "enk_gbr.pdb")
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" truncate_to_polyala=true'%file_name
  run_command(command=cmd, verbose=verbose)
  ala_atom_names = [" N  ", " CA ", " C  ", " O  ", " CB "]
  pdb_inp = iotbx.pdb.hierarchy.input(file_name="enk_gbr.pdb_modified.pdb")
  counter = 0
  for a in pdb_inp.hierarchy.atoms_with_labels():
    assert a.name in ala_atom_names
    counter += 1
  assert counter == 23

def exercise_set_charge () :
  from iotbx import file_reader
  input_pdb = """
ATOM      1  CL  CL  X   1       0.000   0.000   0.000  1.00 20.00          CL
END
"""
  open("tmp_cl.pdb", "w").write(input_pdb)
  easy_run.call("phenix.pdbtools tmp_cl.pdb charge_selection='element Cl' charge=-1 --quiet")
  pdb_in = file_reader.any_file("tmp_cl.pdb_modified.pdb").file_object
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
  open(ifn,"w").write(input_pdb)
  easy_run.call('phenix.pdbtools "%s" renumber_residues=true'%ifn)
  for line1, line2 in zip(open(ifn+"_modified.pdb").readlines(), expected_output_pdb.splitlines()):
    line1 = line1.strip()
    line2 = line2.strip()
    assert line1 == line2

def exercise_set_seg_id () :
  input_pdb = """\
ATOM      1  O   GLY A   3       1.434   1.460   2.496  1.00  6.04           O
ATOM      2  O   CYS A   7       2.196   4.467   3.911  1.00  4.51           O
ATOM      3  O   CYS A   1      -1.433   4.734   5.405  1.00  7.85           O
TER
ATOM      4  O   SER B   4       0.297   0.843   7.226  1.00  7.65           O
ATOM      5  OG ASER B   4      -2.625   1.057   4.064  0.50  5.46           O
ATOM      6  OG BSER B   4      -0.885   0.189   3.843  0.50 11.74           O
TER
ATOM      7  O   LEU     0       5.613  12.448   6.864  1.00  7.32           O
TER
END
"""
  open("tmp_seg_id.pdb", "w").write(input_pdb)
  easy_run.call("phenix.pdbtools tmp_seg_id.pdb set_seg_id_to_chain_id=True --quiet")
  pdb_out_1 = open("tmp_seg_id.pdb_modified.pdb").read()
  assert (pdb_out_1 == """\
ATOM      1  O   GLY A   3       1.434   1.460   2.496  1.00  6.04      A    O
ATOM      2  O   CYS A   7       2.196   4.467   3.911  1.00  4.51      A    O
ATOM      3  O   CYS A   1      -1.433   4.734   5.405  1.00  7.85      A    O
TER
ATOM      4  O   SER B   4       0.297   0.843   7.226  1.00  7.65      B    O
ATOM      5  OG ASER B   4      -2.625   1.057   4.064  0.50  5.46      B    O
ATOM      6  OG BSER B   4      -0.885   0.189   3.843  0.50 11.74      B    O
TER
ATOM      7  O   LEU     0       5.613  12.448   6.864  1.00  7.32           O
TER
END
""")
  easy_run.call("phenix.pdbtools tmp_seg_id.pdb_modified.pdb clear_seg_id=True --quiet")
  pdb_out_2 = open("tmp_seg_id.pdb_modified.pdb_modified.pdb").read()
  assert (pdb_out_2 == input_pdb)

def exercise_remove_first_n_atoms_fraction(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "enk_gbr.pdb")
  n_atoms_start = iotbx.pdb.hierarchy.input(
    file_name=file_name).xray_structure_simple().scatterers().size()
  assert file_name.find('"') < 0
  cmd = 'phenix.pdbtools "%s" remove_first_n_atoms_fraction=0.6'%file_name
  run_command(command=cmd, verbose=verbose)
  pdb_inp = iotbx.pdb.hierarchy.input(file_name="enk_gbr.pdb_modified.pdb")
  n_atoms_final = iotbx.pdb.hierarchy.input(
    file_name="enk_gbr.pdb_modified.pdb").xray_structure_simple().scatterers().size()
  assert n_atoms_final*1./n_atoms_start == 0.4

def exercise(args):
  if ("--show-everything" in args):
    verbose = 2
  elif ("--verbose" in args):
    verbose = 1
  else:
    verbose = 0
  pdb_dir = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb", test=os.path.isdir)
  if (pdb_dir is None):
    print "Skipping exercise(): input files not available"
    return
  eargs = {"pdb_dir": pdb_dir, "verbose": verbose}
  exercise_basic(**eargs)
  exercise_multiple(**eargs)
  exercise_no_cryst1(**eargs)
  exercise_show_adp_statistics(**eargs)
  exercise_show_geometry_statistics(**eargs)
  exercise_show_number_of_removed(**eargs)
  exercise_02(**eargs)
  exercise_truncate_to_polyala(**eargs)
  exercise_renumber_residues()
  exercise_remove_first_n_atoms_fraction(**eargs)
  exercise_set_seg_id()
  exercise_set_charge()
  print "OK"

if (__name__ == "__main__"):
  show_times_at_exit()
  exercise(sys.argv[1:])
