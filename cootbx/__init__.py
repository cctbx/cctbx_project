
from __future__ import absolute_import, division, print_function
import os.path

def write_disable_nomenclature_errors(f):
  f.write("try :\n")
  f.write("  set_nomenclature_errors_on_read(\"ignore\")\n")
  f.write("except Exception :\n")
  f.write("  pass\n")

def create_refinement_view_script(
    mtz_file_name,
    pdb_file_name,
    coot_script_name="view_in_coot.py",
    work_dir=None,
    show_symmetry=True,
    peaks_file_name=None,
    bad_ligand_list=None,
    placed_ligand_list=None):
  from iotbx.file_reader import any_file
  from libtbx.utils import concatenate_python_script
  import libtbx.load_env
  have_anom_map = False
  have_anom_residual_map = False
  mtz_in = any_file(mtz_file_name).assert_file_type("hkl")
  have_anom_map = have_residual_map = False
  for array in mtz_in.file_server.miller_arrays :
    labels = array.info().labels
    if ("ANOM" in labels):
      have_anom_map = True
    elif ("ANOMDIFF" in labels):
      have_anom_residual_map = True
  f = open(coot_script_name, "w")
  print("import coot", file=f)
  print("import os", file=f)
  write_disable_nomenclature_errors(f)
  load_script = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/cootbx/view_refinement.py",
    test=os.path.isfile)
  assert (load_script is not None)
  concatenate_python_script(out=f, file_name=load_script)
  zoom_ligand_script = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/cootbx/simple_zoom_list.py",
    test=os.path.isfile)
  concatenate_python_script(out=f, file_name=zoom_ligand_script)
  if (work_dir is not None):
    pdb_file_name = os.path.basename(pdb_file_name)
    mtz_file_name = os.path.basename(mtz_file_name)
  f.write("""load_refinement(\n""")
  f.write("""pdb_file="%s",\n""" % pdb_file_name)
  f.write("""map_file="%s",\n""" % mtz_file_name)
  f.write("""show_symmetry=%s,\n""" % show_symmetry)
  f.write("""have_anom_map=%s,\n""" % have_anom_map)
  f.write("""have_residual_map=%s,\n""" % have_residual_map)
  if (work_dir is not None):
    f.write("""work_dir="%s",\n""" % work_dir)
  if (peaks_file_name is not None):
    f.write("""peaks_file="%s",\n""" % peaks_file_name)
  f.write(")\n")
  if (bad_ligand_list is not None) and (len(bad_ligand_list) > 0):
    print("""draw_simple_zoom_list(""", file=f)
    print("""  title="Residues in suspicious density",""", file=f)
    print("""  items=%s)""" % str(bad_ligand_list), file=f)
  if (placed_ligand_list is not None):
    print("""draw_simple_zoom_list(""", file=f)
    print("""  title="Placed ligands",""", file=f)
    print("""  items=%s)""" % str(placed_ligand_list), file=f)
  f.close()
