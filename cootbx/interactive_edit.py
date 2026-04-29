
from __future__ import absolute_import, division, print_function
import shutil
import time
import os
import sys

t_wait = 250 # in milliseconds

#-----------------------------------------------------------------------
# This file is run with Coot's python, so no access to CCTBX modules
if sys.hexversion >= 0x03000000:
  unicode = str

def to_bytes(text, codec=None, errors='replace'):
  '''
  Function for handling text when it is passed to cctbx functions that expect
  bytestrings

  Changes text string type (unicode in Python 2, str in Python 3) to
  bytestring type (str or bytes in Python 2, bytes in Python 3)

  The input is returned unmodified if it is already a bytestring
  Will convert other types (e.g. int, float) to bytestring
  None is returned as None, not as 'None'

  For Linux/OS X, the default filesystem encoding is utf8.
  For Windows, the default filesystem encoding is mbcs
  This is important for handling files with basic Python functions.
  With the wrong encoding, the filesystem will not recognize the file path
  import sys; sys.getfilesystemencoding()
  '''

  if (codec is None):
    codec = 'utf8'
    if (sys.platform == 'win32'):
      codec = 'mbcs'

  if (isinstance(text, bytes)):
    return text
  elif (isinstance(text, unicode)):
    new_text = text
    try:
      new_text = text.encode(codec, errors)
    except UnicodeEncodeError: # in case errors='strict'
      raise Sorry('Unable to encode text with %s' % codec)
    return new_text
  elif (text is not None):
    return bytes(text)
  else:
    return None

#-----------------------------------------------------------------------
# Phenix side
def start_coot_and_wait(
    pdb_file,
    map_file,
    data_file,
    work_dir=None,
    coot_cmd="coot",
    needs_rebuild=False,
    log=None):
  if (log is None) : log = sys.stdout
  if (work_dir is None) : work_dir = os.getcwd()
  if (not os.path.isdir(work_dir)):
    os.makedirs(work_dir)
  import mmtbx.maps.utils
  from libtbx.str_utils import make_header
  from libtbx import easy_run
  from libtbx import group_args
  import cootbx
  base_script = __file__.replace(".pyc", ".py")
  os.chdir(work_dir)
  if (os.path.exists("coot_out_tmp.pdb")):
    os.remove("coot_out_tmp.pdb")
  if (os.path.exists("coot_out.pdb")):
    os.remove("coot_out.pdb")
  f = open("edit_in_coot.py", "w")
  f.write(open(base_script).read())
  f.write("\n")
  f.write("import coot\n")
  cootbx.write_disable_nomenclature_errors(f)
  f.write("m = manager(\"%s\", \"%s\", needs_rebuild=%s)\n" %
    (pdb_file, map_file, needs_rebuild))
  f.close()
  make_header("Interactive editing in Coot", log)
  easy_run.call("\"%s\" --no-state-script --script edit_in_coot.py &" %
    coot_cmd)
  print("  Waiting for coot_out_tmp.pdb to appear at %s" % \
    str(time.asctime()), file=log)
  base_dir = os.path.dirname(pdb_file)
  tmp_file = os.path.join(base_dir, "coot_out_tmp.pdb")
  edit_file = os.path.join(base_dir, "coot_tmp_edits.pdb")
  maps_file = os.path.join(base_dir, ".NEW_MAPS")
  while (True):
    if (os.path.isfile(tmp_file)):
      print("  Coot editing complete at %s" % str(time.asctime()), file=log)
      break
    elif (os.path.isfile(maps_file)):
      t1 = time.time()
      assert os.path.isfile(edit_file)
      mmtbx.maps.utils.create_map_from_pdb_and_mtz(
        pdb_file=edit_file,
        mtz_file=data_file,
        output_file=os.path.join(base_dir, "maps_for_coot.mtz"),
        fill=True,
        out=log)
      t2 = time.time()
      print("Calculated new map coefficients in %.1fs" % (t2-t1), file=log)
      os.remove(maps_file)
    else :
      time.sleep(t_wait/1000.)
  shutil.move(tmp_file, "coot_out.pdb")
  mmtbx.maps.utils.create_map_from_pdb_and_mtz(
    pdb_file="coot_out.pdb",
    mtz_file=data_file,
    output_file="coot_out_maps.mtz",
    fill=True,
    out=log)
  new_model = os.path.join(work_dir, "coot_out.pdb")
  new_map = os.path.join(work_dir, "coot_out_maps.mtz")
  skip_rebuild = None
  if (needs_rebuild):
    if (os.path.isfile(os.path.join(base_dir, "NO_BUILD"))):
      skip_rebuild = True
    else :
      skip_rebuild = False
  return group_args(
    pdb_file=new_model,
    map_file=new_map,
    skip_rebuild=skip_rebuild)

#-----------------------------------------------------------------------
# Coot side
class manager(object):
  def __init__(self, pdb_file, map_file, needs_rebuild=False):
    self.file_name = pdb_file
    self.needs_rebuild = needs_rebuild
    import gtk
    import coot
    import coot_python
    toolbar = coot_python.main_toolbar()
    assert (toolbar is not None)
    return_button = gtk.ToolButton()
    return_button.set_label("Return to PHENIX")
    return_button.set_is_important(True)
    toolbar.insert(return_button, -1)
    return_button.connect("clicked", self.OnReturn)
    return_button.show()
    maps_button = gtk.ToolButton()
    maps_button.set_label("Fetch new maps")
    maps_button.set_is_important(True)
    toolbar.insert(maps_button, -1)
    maps_button.connect("clicked", self.OnNewMaps)
    maps_button.show()
    self._imol = read_pdb(to_bytes(pdb_file))
    set_molecule_bonds_colour_map_rotation(self._imol, 30)
    self._map_mols = []
    self.load_maps(map_file)
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK, "You may now edit your model to prepare for further "+
      "building and refinement.   When you are finished, "+
      "click the \"Save model and return\" button on the toolbar, then close "+
      "Coot.  If you want to recalculate the map coefficients, click \"Fetch "+
      "new maps\".")
    dialog.run()
    dialog.destroy()

  def load_maps(self, map_file):
    if (len(self._map_mols) > 0):
      set_colour_map_rotation_for_map(0)
      for imol in self._map_mols :
        close_molecule(imol)
    else :
      set_colour_map_rotation_for_map(10)
    print("Loading %s" % to_bytes(map_file))
    self._map_mols = []
    map_imol = auto_read_make_and_draw_maps(to_bytes(map_file))
    if (isinstance(map_imol, int)):
      # XXX this may be dangerous, but auto_read_make_and_draw_maps only returns
      # the last imol
      self._map_mols = [ map_imol - 1, map_imol ]
    else :
      self._map_mols = map_imol
    set_imol_refinement_map(self._map_mols[0])
    set_scrollable_map(self._map_mols[0])

  def OnReturn(self, *args):
    import coot
    import gtk
    dir_name = os.path.dirname(self.file_name)
    pdb_out = os.path.join(dir_name, "coot_out_tmp.pdb")
    if (self.needs_rebuild):
      dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
        gtk.BUTTONS_YES_NO,
        "Phenix has determined that your model needs to rebuilt in "+
        "AutoBuild; however, if you have made significant changes manually "+
        "this may no longer be necessary.  Do you want to run AutoBuild now? "+
        "(Clicking \"No\" will skip to the next step.")
      response = dialog.run()
      if (dialog != gtk.RESPONSE_YES):
        open(os.path.join(dir_name, "NO_BUILD"), "w").write("1")
      dialog.destroy()
    save_coordinates(self._imol, pdb_out)
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK,
      ("The PDB file has been saved to %s.  Click the \"Okay\" button to "+
       "close Coot and continue.") % pdb_out)
    dialog.run()
    dialog.destroy()
    gtk.main_quit()

  def OnNewMaps(self, *args):
    import coot
    import gtk
    import gobject
    dir_name = os.path.dirname(self.file_name)
    pdb_out = os.path.join(dir_name, "coot_tmp_edits.pdb")
    maps_in = os.path.join(dir_name, "maps_for_coot.mtz")
    if (os.path.exists(maps_in)):
      os.remove(maps_in)
    write_pdb_file(self._imol, pdb_out)
    dialog = gtk.MessageDialog(None, gtk.DIALOG_MODAL, gtk.MESSAGE_INFO,
      gtk.BUTTONS_OK_CANCEL,
      "Phenix will generate new maps; these will be automatically loaded "+
      "when complete.")
    response = dialog.run()
    if (response == gtk.RESPONSE_OK):
      print("WRITING .NEW_MAPS")
      open(os.path.join(dir_name, ".NEW_MAPS"), "w").write("1")
      gobject.timeout_add(t_wait, self.OnWaitForMaps)
    dialog.destroy()

  def OnWaitForMaps(self, *args):
    # TODO a spinner or progress bar might be nice here...
    dir_name = os.path.dirname(self.file_name)
    maps_in = os.path.join(dir_name, "maps_for_coot.mtz")
    if (os.path.isfile(maps_in)):
      self.load_maps(maps_in)
      os.remove(maps_in)
      return False # this kills the timeout
    return True
