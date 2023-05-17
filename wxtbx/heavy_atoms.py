from __future__ import absolute_import, division, print_function

import wxtbx.phil_controls.floatctrl
import wx.lib.mixins.listctrl
import wx
from libtbx.utils import Sorry
import sys
from six.moves import range

# XXX TEMPORARY PYTHON 3 FIX TT
class SitesList(wx.ListCtrl): #wx.lib.mixins.listctrl.CheckListCtrlMixin):
  """
  ListCtrl for displaying and editing heavy-atom sites.  Only the occupancy
  and the selection of sites may be changed.
  """
  def __init__(self, *args, **kwds):
    wx.ListCtrl.__init__(self, *args, **kwds)
    # XXX TEMPORARY PYTHON 3 FIX TT
    #wx.lib.mixins.listctrl.CheckListCtrlMixin.__init__(self)
    self._atoms = None
    self._symm = None
    self.InsertColumn(0, "#", format=wx.LIST_FORMAT_RIGHT, width=64)
    self.InsertColumn(1, "Type", width=64)
    self.InsertColumn(2, "X", format=wx.LIST_FORMAT_RIGHT, width=80)
    self.InsertColumn(3, "Y", format=wx.LIST_FORMAT_RIGHT, width=80)
    self.InsertColumn(4, "Z", format=wx.LIST_FORMAT_RIGHT, width=80)
    self.InsertColumn(5, "Occupancy", format=wx.LIST_FORMAT_RIGHT, width=120)

  def LoadFile(self, file_name):
    """
    Read sites from a PDB file.  An error will be raised if the model contains
    other residue types (protein, NA, water, etc.).
    """
    import iotbx.pdb
    pdb_in = iotbx.pdb.input(file_name)
    self._symm = pdb_in.crystal_symmetry()
    hierarchy = pdb_in.construct_hierarchy()
    counts = hierarchy.overall_counts()
    if (counts.n_models > 1):
      raise Sorry("Multi-model PDB files not allowed.")
    elif (("common_amino_acid" in counts.resname_classes) or
          ("common_nucleic_acid" in counts.resname_classes) or
          ("common_water" in counts.resname_classes)):
      raise Sorry("This PDB file appears to contain non-heavy atoms.")
    atoms = hierarchy.atoms()
    self.LoadAtoms(atoms)

  def LoadAtoms(self, atoms):
    """
    Populate the list with contents of an atom array.  This assumes that the
    contents have already been checked for suitability.
    """
    atoms.set_chemical_element_simple_if_necessary()
    self.DeleteAllItems()
    self._atoms = atoms
    for atom in atoms :
      item = self.InsertStringItem(sys.maxunicode, atom.serial.strip())
      self.SetStringItem(item, 1, atom.element.strip())
      self.SetStringItem(item, 2, "%.3f" % atom.xyz[0])
      self.SetStringItem(item, 3, "%.3f" % atom.xyz[1])
      self.SetStringItem(item, 4, "%.3f" % atom.xyz[2])
      self.SetStringItem(item, 5, "%.2f" % atom.occ)
      if hasattr(self,'CheckItem'):  # XXX Python 3 fix
        self.CheckItem(item)

  def SetSymmetry(self, symmetry):
    self._symm = symmetry

  def SaveSites(self, file_name):
    if (self._atoms is None) or (len(self._atoms) == 0):
      raise Sorry("No atoms loaded!")
    from scitbx.array_family import flex
    selection = flex.bool(self._atoms.size(), False)
    assert (len(selection) == self.GetItemCount())
    for item in range(self.GetItemCount()):
      if (self.IsChecked(item)):
        selection[item] = True
    atoms_selected = self._atoms.select(selection)
    if (len(atoms_selected) == 0):
      raise Sorry("No atoms selected!")
    f = open(file_name, "w")
    if (self._symm is not None):
      from iotbx.pdb import format_cryst1_and_scale_records
      f.write(format_cryst1_and_scale_records(self._symm) + "\n")
    for atom in atoms_selected :
      f.write(atom.format_atom_record() + "\n")
    f.close()
    wx.MessageBox(("%d atoms saved to %s.") % (len(atoms_selected), file_name))

  def ChangeOccupancy(self):
    if (self._atoms is None) or (len(self._atoms) == 0):
      raise Sorry("No atoms loaded!")
    item = self.GetFirstSelected()
    if (item < 0):
      raise Sorry("Please select a site first.")
    assert (item < len(self._atoms))
    dlg = OccupancyDialog(self, -1, "Set occupancy")
    atom = self._atoms[item]
    dlg.SetOccupancy(atom.occ, atom.serial.strip())
    if (dlg.ShowModal() == wx.ID_OK):
      new_occ = dlg.GetOccupancy()
      assert (0. <= new_occ <= 1.)
      atom.set_occ(new_occ)
      self.SetStringItem(item, 5, "%.2f" % atom.occ)

  def OnSaveSites(self, event):
    from iotbx import file_reader
    # XXX Python 3 fix needed here (flags and wildcard not allowed)
    wildcards = file_reader.get_wildcard_strings(["pdb"])
    file_name = wx.FileSelector(
      flags=wx.FD_SAVE,
      wildcard=wildcards,
      )
    if (file_name != ""):
      self.SaveSites(file_name)

  def OnChangeOccupancy(self, event):
    self.ChangeOccupancy()

  def OnLoadSites(self, event):
    from iotbx import file_reader
    wildcards = file_reader.get_wildcard_strings(["pdb"])
    file_name = wx.FileSelector(
      wildcard=wildcards)
    if (file_name != ""):
      self.LoadFile(file_name)

class OccupancyDialog(wx.Dialog):
  def __init__(self, *args, **kwds):
    wx.Dialog.__init__(self, *args, **kwds)
    style = self.GetWindowStyle()
    style |= wx.WS_EX_VALIDATE_RECURSIVELY|wx.RAISED_BORDER|wx.CAPTION
    self.SetWindowStyle(style)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    szr2 = wx.BoxSizer(wx.VERTICAL)
    szr.Add(szr2, 1, wx.ALL|wx.EXPAND, 10)
    self.msg_text = wx.StaticText(self, -1,
      "Changing occupancy for site <unknown> (original occupancy = 0.00).")
    szr2.Add(self.msg_text, 0, wx.ALL, 5)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    szr2.Add(szr3)
    txt = wx.StaticText(self, -1, "New occupancy:")
    szr3.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.occ_ctrl = wxtbx.phil_controls.floatctrl.FloatCtrl(
      parent=self,
      size=(80,-1),
      value=0.)
    self.occ_ctrl.SetMin(0.)
    self.occ_ctrl.SetMax(1.)
    self.occ_ctrl.SetOptional(False)
    szr3.Add(self.occ_ctrl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    cancel_btn = wx.Button(self, wx.ID_CANCEL)
    ok_btn = wx.Button(self, wx.ID_OK)
    ok_btn.SetDefault()
    szr4 = wx.StdDialogButtonSizer()
    szr4.Add(cancel_btn)
    szr4.Add(ok_btn, 0, wx.LEFT, 5)
    szr2.Add(szr4, 0, wx.ALL|wx.ALIGN_RIGHT, 5)
    szr.Layout()
    self.Fit()
    self.Centre(wx.BOTH)

  def SetOccupancy(self, occ, serial):
    assert (isinstance(occ, float))
    self.occ_ctrl.SetValue(occ)
    self.msg_text.SetLabel(
      "Changing occupancy for site %s (original occupancy = %.2f)." % (serial,
      occ))
    self.Refresh()

  def GetOccupancy(self):
    return self.occ_ctrl.GetPhilValue()

  def OnOkay(self, event):
    if (not self.Validate()):
      pass
    else :
      occ = self.GetOccupancy()
      self.EndModal(wx.ID_OK)

  def OnCancel(self, event):
    self.EndModal(wx.ID_CANCEL)

class sites_panel_mixin(object):
  """
  Mixin class for panel objects which contain a SitesList and buttons to
  trigger editing functions.  Should be subclassed along with wx.Panel,
  wx.Dialog, or similar.
  """
  def create_sites_list(self,
                         sizer,
                         label="Heavy atom sites",
                         show_load_button=False):
    if (label is not None):
      txt = wx.StaticText(self, -1, label + ":")
      sizer.Add(txt, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.sites_list = SitesList(self, -1,
      size=(540,200),
      style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
    # XXX  TEMPORARY PYTHON 3 FIX TT
    sizer.Add(self.sites_list, 1, wx.ALL, 5) # |wx.EXPAND|wx.ALIGN_CENTER_VERTICAL, 5)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    sizer.Add(szr2)
    if (show_load_button):
      btn0 = wx.Button(self, -1, "Load sites...")
      szr2.Add(btn0, 0, wx.TOP|wx.LEFT|wx.BOTTOM, 5)
      self.Bind(wx.EVT_BUTTON, self.sites_list.OnLoadSites)
    btn1 = wx.Button(self, -1, "Change occupancy...")
    szr2.Add(btn1, 0, wx.TOP|wx.LEFT|wx.BOTTOM, 5)
    btn2 = wx.Button(self, -1, "Save selected...")
    szr2.Add(btn2, 0, wx.TOP|wx.LEFT|wx.BOTTOM, 5)
    self.Bind(wx.EVT_BUTTON, self.sites_list.OnChangeOccupancy, btn1)
    self.Bind(wx.EVT_BUTTON, self.sites_list.OnSaveSites, btn2)

  def load_sites(self, file_name=None, atoms=None, symmetry=None):
    assert ([file_name, atoms].count(None) == 1)
    if (file_name is not None):
      assert (symmetry is None)
      self.sites_list.LoadFile(file_name)
    else :
      self.sites_list.LoadAtoms(atoms)
      if (symmetry is not None):
        self.sites_list.SetSymmetry(symmetry)

class sites_panel(sites_panel_mixin, wx.Panel):
  def __init__(self, *args, **kwds):
    wx.Panel.__init__(self, *args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self.create_sites_list(szr, show_load_button=True)

########################################################################
# REGRESSION TESTING
def exercise():
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
CRYST1   77.324  107.578   84.396  90.00  94.23  90.00 P 1 21 1
SCALE1      0.012933  0.000000  0.000956        0.00000
SCALE2      0.000000  0.009296  0.000000        0.00000
SCALE3      0.000000  0.000000  0.011881        0.00000
ATOM      1 SE   SUB     1     -27.144 -52.994 -26.775  1.00  0.00          SE
ATOM      2 SE   SUB     2     -28.543 -78.848 -28.830  0.96  0.00          SE
ATOM      3 SE   SUB     3      -6.147  -3.557 -53.677  0.87  0.00          SE
ATOM      4 SE   SUB     4      -4.949  -0.129 -55.688  0.84  0.00          SE
ATOM      5 SE   SUB     5      -2.374  -4.774 -71.923  0.84  0.00          SE
ATOM      6 SE   SUB     6     -32.758 -55.851 -15.545  0.78  0.00          SE
ATOM      7 SE   SUB     7     -27.582 -84.707 -12.897  0.77  0.00          SE
ATOM      8 SE   SUB     8     -27.552 -80.117  -9.489  0.77  0.00          SE
ATOM      9 SE   SUB     9      -3.273 -30.674 -70.235  0.69  0.00          SE
ATOM     10 SE   SUB    10      -0.502 -86.255 -44.849  0.66  0.00          SE
ATOM     11 SE   SUB    11      -4.080 -26.592 -75.056  0.66  0.00          SE
ATOM     12 SE   SUB    12       0.931 -27.847 -57.819  0.65  0.00          SE
ATOM     13 SE   SUB    13       1.508  -2.738 -67.310  0.63  0.00          SE
ATOM     14 SE   SUB    14     -30.096 -32.213 -38.946  0.58  0.00          SE
ATOM     15 SE   SUB    15      -2.988 -25.400 -54.434  0.57  0.00          SE
ATOM     16 SE   SUB    16     -30.208 -59.010 -10.822  0.52  0.00          SE
END""")
  atoms = pdb_in.construct_hierarchy().atoms()
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Heavy-atom sites editor")
  sizer = wx.BoxSizer(wx.VERTICAL)
  frame.SetSizer(sizer)
  panel = sites_panel(frame)
  sizer.Add(panel, 1, wx.EXPAND)
  sizer.Fit(panel)
  frame.Fit()
  panel.load_sites(atoms=atoms, symmetry=pdb_in.crystal_symmetry())
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__"):
  exercise()
