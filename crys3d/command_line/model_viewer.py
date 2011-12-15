# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_SIGNALS_DEFAULT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import cStringIO
from crys3d.wx_selection_editor import selection_editor_mixin
import wx
import libtbx.load_env
import sys, os, time

########################################################################
# CLASSES AND METHODS FOR STANDALONE VIEWER
#
class App (wx.App) :
  def __init__ (self, title="crys3d.wx_model_viewer", default_size=(800,600),
      viewer_class=selection_editor_mixin) :
    self.title = title
    self.default_size = default_size
    self.viewer_class = viewer_class
    wx.App.__init__(self, 0)

  def OnInit (self) :
    self.frame = wx.Frame(None, -1, self.title, pos=wx.DefaultPosition,
      size=self.default_size)
    self.frame.CreateStatusBar()
    box = wx.BoxSizer(wx.VERTICAL)
    self.view_objects = self.viewer_class(self.frame, size=(800,600))
    box.Add(self.view_objects, wx.EXPAND, wx.EXPAND)
    self.frame.SetSizer(box)
    box.SetSizeHints(self.frame)
    return True

def run (args, viewer_class=selection_editor_mixin) :
  import cStringIO
  pdb_files = []
  cif_files = []
  show_ss_restraints = False
  fast_connectivity = True
  for arg in args :
    if os.path.isfile(arg) :
      import iotbx.pdb
      if iotbx.pdb.is_pdb_file(arg) :
        pdb_files.append(os.path.abspath(arg))
      elif arg.endswith(".cif") :
        cif_files.append(os.path.abspath(arg))
    elif arg == "--ss" :
      show_ss_restraints = True
    elif arg in ["--thorough", "--slow", "--use_monomer_library"] :
      fast_connectivity = False
  if len(pdb_files) == 0 :
    print "Please specify a PDB file (and optional CIFs) on the command line."
    return
  a = App(viewer_class=viewer_class)
  a.frame.Show()
  out = sys.stdout
  if not "--debug" in args :
    out = cStringIO.StringIO()
  for file_name in pdb_files :
    print "Reading PDB file %s" % file_name
    from iotbx import file_reader
    from mmtbx.monomer_library import pdb_interpretation
    from mmtbx import secondary_structure
    t1 = time.time()
    if fast_connectivity :
      pdb_in = file_reader.any_file(file_name, force_type="pdb")
      pdb_hierarchy = pdb_in.file_object.construct_hierarchy()
      atomic_bonds = pdb_hierarchy.distance_based_simple_two_way_bond_sets()
      acp_selection = None
    else :
      processed_pdb_file = pdb_interpretation.run(args=[file_name]+cif_files,
        log=out)
      pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
      pdb_hierarchy.atoms().reset_i_seq()
      grm = processed_pdb_file.geometry_restraints_manager()
      acp_selection = processed_pdb_file.all_chain_proxies.selection
      if grm is None or grm.shell_sym_tables is None :
        raise Sorry("Atomic bonds could not be calculated for this model. "+
          "This is probably due to a missing CRYST1 record in the PDB file.")
      atomic_bonds = grm.shell_sym_tables[0].full_simple_connectivity()
    t2 = time.time()
    print "%.2fs" % (t2-t1)
    a.view_objects.add_model(file_name, pdb_hierarchy, atomic_bonds,
      mmtbx_selection_function=acp_selection)
    sec_str = secondary_structure.manager(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=None)
    sec_str.find_automatically()
    a.view_objects.set_sec_str(file_name, sec_str.selections_as_ints())
    if show_ss_restraints and acp_selection is not None :
      bonds_table = secondary_structure.process_structure(params=None,
        processed_pdb_file=processed_pdb_file,
        tmp_dir=os.getcwd(),
        log=sys.stderr)
      a.view_objects.set_noncovalent_bonds(file_name, bonds_table.bonds)
      a.view_objects.flag_show_noncovalent_bonds = True
      a.view_objects.set_model_base_color([1.0,1.0,1.0], file_name)
      a.view_objects.set_color_mode("element")
  a.view_objects.force_update(recenter=True)
  a.MainLoop()

if __name__ == "__main__" :
  if "--test" in sys.argv :
    pdb_file = libtbx.env.find_in_repositories(
      relative_path="phenix_regression/pdb/1ywf.pdb",
      test=os.path.isfile)
    run([pdb_file, "--ss"])
  else :
    run(sys.argv[1:])
