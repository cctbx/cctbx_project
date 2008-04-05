from libtbx import easy_run
import libtbx.load_env
import sys, os, math, time
from cctbx.array_family import flex
from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import remove_files, search_for
from mmtbx import utils
from libtbx.test_utils import approx_equal, not_approx_equal, run_command, \
  show_diff
from libtbx.utils import show_times_at_exit
from iotbx import reflection_file_reader
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.f_model


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
  base = \
      "phenix.pdbtools %s output.file_name=%s "%(file_name, output)
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
  cmd= "phenix.pdbtools %s output.file_name=%s shake=0.1 --quiet > log"%(
                                        file_name, output_file_name)
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
    assert approx_equal(xrsp.u_cart, xrsp_init.u_cart,tolerance)
  else:
    arg1 = xrsp.u_iso_used.select(selection.select(xrsp.use_u_iso))
    arg2 = xrsp_init.u_iso_used.select(selection.select(xrsp_init.use_u_iso))
    if(arg1.size() > 0): assert not_approx_equal(arg1, arg2,tolerance)
    assert approx_equal(xrsp.u_cart, xrsp_init.u_cart,tolerance)

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
adp {
  selection = chain A
  randomize = True
}
adp {
  selection = chain B
  shift_b_iso = 10
}
sites {
  selection = chain B
  shake = 1.5
}
sites {
  selection = chain A or chain C
  translate = 1 2 3
  rotate = 4 5 6
}
"""
  open("params", "w").write(params)
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  cmd = "phenix.pdbtools %s output.file_name=modified.pdb params" % (
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
  base = "phenix.pdbtools %s output.file_name=%s "%(file_name, output)
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
  cmd = "phenix.pdbtools %s --show-adp-statistics"%file_name
  run_command(command=cmd, verbose=verbose)

def exercise_show_geometry_statistics(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  cmd = "phenix.pdbtools %s --show-geometry-statistics"%file_name
  run_command(command=cmd, verbose=verbose)

def exercise_f_model_option_default(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_e.pdb")
  high_resolution = 3.0
  cmd = " ".join(['phenix.pdbtools',
                  '%s'%file_name,
                  'high_resolution=%s'%str(high_resolution),
                  '--f_model'])
  run_command(command=cmd, verbose=verbose)
  reflection_file = reflection_file_reader.any_reflection_file(
    file_name          = os.path.basename(file_name)+".mtz",
    ensure_read_access = False)
  miller_arrays = reflection_file.as_miller_arrays()
  assert len(miller_arrays) == 1
  ma = miller_arrays[0]
  assert "%s"%ma.info() == \
    "%s:FMODEL,PHIFMODEL"%(os.path.basename(file_name)+".mtz")
  xray_structure = iotbx.pdb.input(file_name=file_name).xray_structure_simple()
  xray_structure.scattering_type_registry(table="n_gaussian",d_min =ma.d_min())
  assert ma.is_complex_array()
  assert approx_equal(ma.d_min(), high_resolution, 0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = ma.generate_r_free_flags(fraction = 0.1),
    target_name    = "ls_wunit_k1",
    f_obs          = abs(ma))
  r1 = abs(ma).data().min_max_mean().as_tuple()
  r2 = abs(fmodel.f_model()).data().min_max_mean().as_tuple()
  r3 = abs(fmodel.f_calc()).data().min_max_mean().as_tuple()
  r4 = abs(fmodel.f_obs).data().min_max_mean().as_tuple()
  assert approx_equal(r1, r2, 1.e-5)
  assert approx_equal(r3, r4, 1.e-5)
  assert approx_equal(r1, r4, 1.e-5)

def exercise_f_model_option_custom(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "enk_gbr.pdb")
  high_resolution = 2.0
  format = "cns"
  type = "real"
  label = "Fobs"
  low_resolution = 6.0
  algorithm = "direct"
  par = (0.35,60,3,[1,2,-3,0,0,0])
  par_str = 'k_sol=%s b_sol=%s scale=%s b_cart="%s"'%(par[0],
    par[1], par[2], " ".join([str(i) for i in par[3]]).strip())
  for type in ["real", "complex"]:
    for table in ["wk1995", "neutron"]:
      cmd = " ".join(
        ['phenix.pdbtools',
         '%s'%file_name,
         'high_resolution=%s'%str(high_resolution),
         'format=%s'%format,
         'type=%s'%type,
         'label=%s'%label,
         'low_resolution=%s'%str(low_resolution),
         'algorithm=%s'%algorithm,
         'scattering_table=%s'%table,
         '%s'%par_str,
         '--f_model'])
      run_command(command=cmd, verbose=verbose)
      xray_structure=iotbx.pdb.input(file_name=file_name).xray_structure_simple()
      reflection_file = reflection_file_reader.any_reflection_file(
        file_name          = os.path.basename(file_name)+".hkl",
        ensure_read_access = False)
      miller_arrays = reflection_file.as_miller_arrays(crystal_symmetry =
        xray_structure.crystal_symmetry())
      assert len(miller_arrays) == 1
      ma = miller_arrays[0]
      if(table == "neutron"):
        xray_structure.switch_to_neutron_scattering_dictionary()
      else:
        xray_structure.scattering_type_registry(table = table,d_min = ma.d_min())
      if(type == "real"): assert ma.is_real_array()
      if(type == "complex"): assert ma.is_complex_array()
      d_max, d_min = ma.d_max_min()
      assert approx_equal(d_min, high_resolution, 0.1)
      assert approx_equal(d_max, low_resolution, 0.1)
      assert "%s"%ma.info() == "%s:FOBS"%(os.path.basename(file_name)+".hkl")
      sf_calc_params = mmtbx.f_model.sf_and_grads_accuracy_params.extract()
      sf_calc_params.algorithm = algorithm
      fmodel = mmtbx.f_model.manager(
        xray_structure = xray_structure,
        sf_and_grads_accuracy_params = sf_calc_params,
        r_free_flags   = ma.generate_r_free_flags(fraction = 0.1),
        target_name    = "ml",
        f_obs          = abs(ma))
      fmodel.update_solvent_and_scale(verbose = -1)
      # tolerances MUST be small, otherwise ring a bell
      assert approx_equal(fmodel.r_work(),        0, 1.e-4)
      assert approx_equal(fmodel.r_free(),        0, 1.e-4)
      assert approx_equal(fmodel.k_sol(),    par[0], 1.e-2)
      assert approx_equal(fmodel.b_sol(),    par[1], 1.e-2)
      assert approx_equal(fmodel.scale_k1(), par[2], 1.e-4)
      assert approx_equal(fmodel.b_cart(),   par[3], 1.e-3)

def exercise_show_number_of_removed(pdb_dir, verbose):
  file_name = os.path.join(pdb_dir, "phe_h.pdb")
  log = "exercise_show_number_of_removed.log"
  cmd = 'phenix.pdbtools %s remove="element H" > %s' % (file_name, log)
  remove_files(log)
  run_command(command=cmd, verbose=verbose)
  assert os.path.isfile(log)
  n_found = 0
  for line in open(log).read().splitlines():
    if (line == "Atoms to be kept: 13 of 24"):
      n_found += 1
  assert n_found == 1

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
  exercise_f_model_option_default(**eargs)
  exercise_f_model_option_custom(**eargs)
  exercise_show_number_of_removed(**eargs)
  print "OK"

if (__name__ == "__main__"):
  show_times_at_exit()
  exercise(sys.argv[1:])
