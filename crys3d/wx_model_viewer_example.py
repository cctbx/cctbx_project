
import sys, os
import cStringIO
from crys3d.wx_model_viewer_multi_scene import model_viewer_mixin
from crys3d.wx_selection_editor import selection_editor_mixin
from mmtbx.monomer_library import pdb_interpretation
import iotbx.pdb
import wx

########################################################################
# CLASSES AND METHODS FOR STANDALONE VIEWER
#
class App (wx.App) :
  def __init__ (self, title="crys3d.wx_model_viewer", default_size=(800,600)) :
    self.title = title
    self.default_size = default_size
    wx.App.__init__(self, 0)

  def OnInit (self) :
    self.frame = wx.Frame(None, -1, self.title, pos=wx.DefaultPosition,
      size=self.default_size)
    self.frame.CreateStatusBar()
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = selection_editor_mixin(self.frame, size=(800,600))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    return True

def run (args) :
  import cStringIO
  pdb_files = []
  cif_files = []
  for arg in args :
    if os.path.isfile(arg) :
      if iotbx.pdb.is_pdb_file(arg) :
        pdb_files.append(os.path.abspath(arg))
      elif arg.endswith(".cif") :
        cif_files.append(os.path.abspath(arg))
  if len(pdb_files) == 0 :
    print "Please specify a PDB file (and optional CIFs) on the command line."
    return
  a = App()
  out = sys.stdout
  if not "--debug" in args :
    out = cStringIO.StringIO()
  for file_name in pdb_files :
    print "Reading PDB file %s" % file_name
    processed_pdb_file = pdb_interpretation.run(args=[file_name]+cif_files,
      log=out)
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    pdb_hierarchy.atoms().reset_i_seq()
    grm = processed_pdb_file.geometry_restraints_manager()
    if grm is None or grm.shell_sym_tables is None :
      raise Sorry("Atomic bonds could not be calculated for this model. "+
        "This is probably due to a missing CRYST1 record in the PDB file.")
    atomic_bonds = grm.shell_sym_tables[0].full_simple_connectivity()
    a.view_objects.add_model(file_name, pdb_hierarchy, atomic_bonds)
  a.frame.Show()
  a.view_objects.force_update(recenter=True)
  a.MainLoop()

if __name__ == "__main__" :
  import sys
  run(sys.argv[1:])
