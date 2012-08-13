
# standard imports
import shutil
import time
import os
import sys

#-----------------------------------------------------------------------
# Phenix side
def start_coot_and_wait (
    pdb_file,
    map_file,
    data_file,
    work_dir=None,
    coot_cmd="coot",
    needs_rebuild=False,
    log=None) :
  if (log is None) : log = sys.stdout
  if (work_dir is None) : work_dir = os.getcwd()
  if (not os.path.isdir(work_dir)) :
    os.makedirs(work_dir)
  from libtbx.str_utils import make_header
  from libtbx import easy_run
  from libtbx import group_args
  base_script = __file__.replace(".pyc", ".py")
  os.chdir(work_dir)
  if (os.path.exists("coot_out_tmp.pdb")) :
    os.remove("coot_out_tmp.pdb")
  if (os.path.exists("coot_out.pdb")) :
    os.remove("coot_out.pdb")
  f = open("edit_in_coot.py", "w")
  f.write(open(base_script).read())
  f.write("\n")
  f.write("import coot\n")
  f.write("m = manager(\"%s\", \"%s\", needs_rebuild=%s)\n" %
    (pdb_file, map_file, needs_rebuild))
  f.close()
  make_header("Interactive editing in Coot", log)
  easy_run.call("\"%s\" --no-state-script --script edit_in_coot.py &" %
    coot_cmd)
  print >> log, "  Waiting for coot_out_tmp.pdb to appear at %s" % \
    str(time.asctime())
  base_dir = os.path.dirname(pdb_file)
  tmp_file = os.path.join(base_dir, "coot_out_tmp.pdb")
  while (True) :
    if (os.path.isfile(tmp_file)) :
      print >> log, "  Coot editing complete at %s" % str(time.asctime())
      break
    else :
      time.sleep(1)
  shutil.move(tmp_file, "coot_out.pdb")
  import mmtbx.maps.utils
  mmtbx.maps.utils.create_map_from_pdb_and_mtz(
    pdb_file="coot_out.pdb",
    mtz_file=data_file,
    output_file="coot_out_maps.mtz",
    fill=True,
    out=log)
  new_model = os.path.join(work_dir, "coot_out.pdb")
  new_map = os.path.join(work_dir, "coot_out_maps.mtz")
  skip_rebuild = None
  if (needs_rebuild) :
    if (os.path.isfile(os.path.join(base_dir, "NO_BUILD"))) :
      skip_rebuild = True
    else :
      skip_rebuild = False
  return group_args(
    pdb_file=new_model,
    map_file=new_map,
    skip_rebuild=skip_rebuild)

#-----------------------------------------------------------------------
# Coot side
class manager (object) :
  def __init__ (self, pdb_file, map_file, needs_rebuild=False) :
    self.file_name = pdb_file
    self.needs_rebuild = needs_rebuild
    import gtk
    import coot
    import coot_python
    auto_read_make_and_draw_maps(map_file)
    toolbar = coot_python.main_toolbar()
    assert (toolbar is not None)
    return_button = gtk.ToolButton()
    return_button.set_label("Return to PHENIX")
    return_button.set_is_important(True)
    toolbar.insert(return_button, -1)
    return_button.connect("clicked", self.OnReturn)
    return_button.show()
    self._imol = read_pdb(pdb_file)
    set_molecule_bonds_colour_map_rotation(self._imol, 30)
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK, "You may now edit your model to prepare for further "+
      "building and refinement.   When you are finished, "+
      "click the \"Save model and return\" button on the toolbar, then close "+
      "Coot.")
    dialog.run()
    dialog.destroy()

  def OnReturn (self, *args) :
    import coot
    import gtk
    dir_name = os.path.dirname(self.file_name)
    pdb_out = os.path.join(dir_name, "coot_out_tmp.pdb")
    if (self.needs_rebuild) :
      dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
        gtk.BUTTONS_YES_NO,
        "Phenix has determined that your model needs to rebuilt in "+
        "AutoBuild; however, if you have made significant changes manually "+
        "this may no longer be necessary.  Do you want to run AutoBuild now? "+
        "(Clicking \"No\" will skip to the next step.")
      response = dialog.run()
      if (dialog != gtk.RESPONSE_YES) :
        open(os.path.join(dir_name, "NO_BUILD"), "w").write("1")
      dialog.destroy()
    write_pdb_file(self._imol, pdb_out)
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK,
      "The PDB file has been saved to %s.  You may now close Coot." % pdb_out)
    dialog.run()
    dialog.destroy()
