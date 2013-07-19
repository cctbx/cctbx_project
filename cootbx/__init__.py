
from __future__ import division
import os.path

def write_disable_nomenclature_errors (f) :
  f.write("try :\n")
  f.write("  set_nomenclature_errors_on_read(\"ignore\")\n")
  f.write("except Exception :\n")
  f.write("  pass\n")

def create_refinement_view_script (
    mtz_file_name,
    pdb_file_name,
    coot_script_name="view_in_coot.py",
    work_dir=None,
    show_symmetry=True,
    peaks_file_name=None,
    bad_ligand_list=None,
    placed_ligand_list=None) :
  from iotbx.file_reader import any_file
  import libtbx.load_env
  have_anom_map = False
  have_anom_residual_map = False
  mtz_in = any_file(mtz_file_name).assert_file_type("hkl")
  have_anom_map = have_residual_map = False
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().labels
    if ("ANOM" in labels) :
      have_anom_map = True
    elif ("ANOMDIFF" in labels) :
      have_anom_residual_map = True
  f = open(coot_script_name, "w")
  print >> f, "import coot"
  print >> f, "import os"
  zoom_ligand_script = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/cootbx/simple_zoom_list.py",
    test=os.path.isfile)
  if (zoom_ligand_script is not None) :
    script = open(zoom_ligand_script).read()
    print >> f, script
    if (bad_ligand_list is not None) and (len(bad_ligand_list) > 0) :
      print >> f, """draw_simple_zoom_list("""
      print >> f, """  title="Residues in suspicious density","""
      print >> f, """  items=%s)""" % str(bad_ligand_list)
    if (placed_ligand_list is not None) :
      print >> f, """draw_simple_zoom_list("""
      print >> f, """  title="Placed ligands","""
      print >> f, """  items=%s)""" % str(placed_ligand_list)
  if (work_dir is not None) :
    print >> f, """os.chdir("%s")""" % work_dir
  print >> f, """read_pdb("%s")""" % pdb_file_name
  if (show_symmetry) :
    print >> f, """set_show_symmetry_master(True)"""
  if (peaks_file_name is not None) :
    print >> f, """imol = read_pdb("%s")""" % peaks_file_name
    print >> f, """set_mol_displayed(imol, False)"""
  print >> f, """auto_read_make_and_draw_maps("%s")""" % mtz_file_name
  if (have_anom_map) :
    print >> f, """imol = make_and_draw_map("%s","ANOM","PHANOM","",0,0)""" % \
      mtz_file_name
    print >> f, """set_contour_level_in_sigma(imol, 3.0)"""
    if (have_residual_map) :
      print >> f, """set_map_colour(imol, 0.0, 1.0, 1.0)""" # cyan
      print >> f, """set_map_displayed(imol, False)"""
    else :
      print >> f, """set_map_colour(imol, 1.0, 1.0, 0.0)""" # yellow
  if (have_residual_map) :
    print >> f, \
      """imol = make_and_draw_map("%s","ANOMDIFF","PHANOMDIFF","",0,1)""" % \
      mtz_file_name
    print >> f, """set_contour_level_in_sigma(imol, 3.0)"""
    print >> f, """set_map_colour(imol, 1.0, 1.0, 0.0)""" # yellow
  f.close()
