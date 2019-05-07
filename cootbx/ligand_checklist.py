
from __future__ import division, print_function
from six.moves import range
import time
import os
import sys

from libtbx.utils import to_str

t_wait = 250

def start_coot_and_wait(
    pdb_file,
    map_file,
    ligand_files,
    ligand_ccs,
    cif_files=(),
    work_dir=None,
    coot_cmd="coot",
    log=None):
  from iotbx import file_reader
  from libtbx.str_utils import make_header
  from libtbx import easy_run
  import cootbx
  assert (len(ligand_files) > 0) and (len(ligand_files) == len(ligand_ccs))
  if (log is None) : log = sys.stdout
  cwd = os.getcwd()
  if (work_dir is None) : work_dir = cwd
  if (not os.path.isdir(work_dir)):
    os.makedirs(work_dir)
  os.chdir(work_dir)
  base_script = __file__.replace(".pyc", ".py")
  ligand_xyzs = []
  for pdb_file in ligand_files :
    pdb_file = to_str(pdb_file)
    pdb_in = file_reader.any_file(pdb_file, force_type="pdb")
    pdb_in.assert_file_type("pdb")
    coords = pdb_in.file_object.atoms().extract_xyz()
    ligand_xyzs.append(coords.mean())
  ligand_info = zip(ligand_files, ligand_ccs, ligand_xyzs)
  f = open("edit_in_coot.py", "w")
  f.write(open(base_script).read())
  f.write("\n")
  f.write("import coot\n")
  cootbx.write_disable_nomenclature_errors(f)
  f.write("read_pdb(\"%s\")\n" % to_str(pdb_file))
  f.write("auto_read_make_and_draw_maps(\"%s\")\n" % to_str(map_file))
  for cif_file in cif_files :
    f.write("read_cif_dictionary(\"%s\")\n" % to_str(cif_file))
  f.write("m = manager(%s)\n" % str(ligand_info))
  f.close()
  make_header("Ligand selection in Coot", log)
  rc = easy_run.call("\"%s\" --no-state-script --script edit_in_coot.py &" %
    coot_cmd)
  if (rc != 0):
    raise RuntimeError("Launching Coot failed with status %d" % rc)
  print("  Waiting for user input at %s" % str(time.asctime()), file=log)
  out_file = ".COOT_LIGANDS"
  output_files = output_ccs = None
  while (True):
    if (os.path.isfile(out_file)):
      print("  Coot editing complete at %s" % str(time.asctime()), file=log)
      ligand_indices = [ int(i) for i in open(out_file).read().split() ]
      output_files = []
      for i in ligand_indices :
        ligand_file = os.path.join(work_dir, "coot_ligand_out_%d.pdb" % (i+1))
        output_files.append(ligand_file)
      output_ccs = [ ligand_ccs[i] for i in ligand_indices ]
      break
    else :
      time.sleep(t_wait / 1000.)
  assert (output_files is not None)
  os.chdir(cwd)
  return output_files, output_ccs

class manager(object):
  def __init__(self, ligand_file_info):
    import gtk
    import coot # import dependency
    title = "Select ligand(s)"
    self.ligand_file_info = ligand_file_info
    self.ligand_imols = []
    for file_name, cc, xyz in ligand_file_info :
      self.ligand_imols.append(read_pdb(to_str(file_name)))
    self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
    self.window.set_default_size(300, 200)
    self.window.set_title(title)
    scrolled_win = gtk.ScrolledWindow()
    outside_vbox = gtk.VBox(False, 2)
    inside_vbox = gtk.VBox(False, 0)
    inside_vbox.set_border_width(2)
    self.window.add(outside_vbox)
    outside_vbox.pack_start(scrolled_win, True, True, 0) # expand fill padding
    scrolled_win.add_with_viewport(inside_vbox)
    scrolled_win.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_ALWAYS)
    frame = gtk.Frame(title)
    vbox = gtk.VBox(False, 0)
    inside_vbox.pack_start(frame, False, False, 2)
    self.model = gtk.ListStore(bool, gobject.TYPE_STRING, gobject.TYPE_STRING)
    tv = gtk.TreeView(self.model)
    cell1 = gtk.CellRendererToggle()
    cell1.connect("toggled", self.OnToggle, self.model)
    col1 = gtk.TreeViewColumn("Keep")
    col1.pack_start(cell1, False)
    col1.add_attribute(cell1, "active", 0)
    cell2 = gtk.CellRendererText()
    col2 = gtk.TreeViewColumn("File")
    col2.pack_start(cell2, False)
    col2.add_attribute(cell2, "text", 1)
    cell3 = gtk.CellRendererText()
    col3 = gtk.TreeViewColumn("CC")
    col3.pack_start(cell3, False)
    col3.add_attribute(cell3, "text", 2)
    tv.append_column(col1)
    tv.append_column(col2)
    tv.append_column(col3)
    frame.add(tv)
    for file_name, cc, xyz in ligand_file_info :
      self.model.append([False, os.path.basename(file_name), cc])
    continue_btn = gtk.Button("Close and continue")
    continue_btn.connect("clicked", self.OnContinue)
    outside_vbox.pack_end(continue_btn, False, False, 0)
    self.window.show_all()

  def OnToggle(self, cell, path, model):
    model[path][0] = not model[path][0]

  def OnContinue(self, *args):
    import gtk
    selected = []
    for i_lig in range(len(self.ligand_file_info)):
      if (self.model[i_lig][0]):
        print("  selected ligand %d" % (i_lig+1))
        pdb_out = "coot_ligand_out_%d.pdb" % (i_lig+1)
        save_coordinates(self.ligand_imols[i_lig], pdb_out)
        selected.append(str(i_lig))
    open(".COOT_LIGANDS", "w").write(" ".join(selected))
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK,
      "The selected ligands have been saved.  You may now close Coot.")
    dialog.run()
    dialog.destroy()
    gtk.main_quit()
