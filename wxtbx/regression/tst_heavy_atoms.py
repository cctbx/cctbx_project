from __future__ import absolute_import, division, print_function

# Regression tests for wxtbx.heavy_atoms.SitesList, the heavy-atom site list
# used by the AutoSol and HySS result pages.  Under wxPython 4 it showed only
# the first site and crashed on "Save selected...".
#
# Runs headlessly: no window is shown, and the two modal wx calls that the
# widget makes (FileSelector, MessageBox) are replaced with recorders.

import os
import shutil
import sys
import tempfile

from wxtbx.regression._warnings import install_wx_deprecation_filters

install_wx_deprecation_filters()

PDB_STR = """\
CRYST1   77.324  107.578   84.396  90.00  94.23  90.00 P 1 21 1
ATOM      1 SE   SUB     1     -27.144 -52.994 -26.775  1.00  0.00          SE
ATOM      2 SE   SUB     2     -28.543 -78.848 -28.830  0.96  0.00          SE
ATOM      3 SE   SUB     3      -6.147  -3.557 -53.677  0.87  0.00          SE
ATOM      4 SE   SUB     4      -4.949  -0.129 -55.688  0.84  0.00          SE
ATOM      5 SE   SUB     5      -2.374  -4.774 -71.923  0.84  0.00          SE
END
"""
N_SITES = 5

def _sites_file(tmp_dir):
  file_name = os.path.join(tmp_dir, "ha_1.pdb")
  with open(file_name, "w") as f:
    f.write(PDB_STR)
  return file_name

def _read_serials(file_name):
  # Raw input atoms keep the serials as written; construct_hierarchy() would
  # renumber them.
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(file_name=file_name)
  return pdb_in.crystal_symmetry(), [int(a.serial) for a in pdb_in.atoms()]

def _keyword_names(func):
  """Parameter names from a wxPython builtin's signature docstring line.

  inspect.signature() rejects SIP builtins, so parse the first docstring
  line, e.g. "FileSelector(message, default_path='', ..., y=DefaultCoord) -> str".
  """
  first = func.__doc__.strip().splitlines()[0]
  inner = first[first.index("(") + 1:first.rindex(")")]
  return [p.strip().split("=")[0] for p in inner.split(",") if p.strip()]

class _recorder(object):
  """Stand-in for a modal wx call: records the call, returns a fixed value."""
  def __init__(self, result=None):
    self.result = result
    self.calls = []
  def __call__(self, *args, **kwds):
    self.calls.append((args, kwds))
    return self.result

def _make_list(frame):
  import wx
  from wxtbx.heavy_atoms import SitesList
  return SitesList(frame, -1, size=(540, 200),
    style=wx.LC_REPORT|wx.LC_SINGLE_SEL)

def exercise_load_shows_every_site(frame, tmp_dir):
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  assert sites.GetItemCount() == N_SITES, sites.GetItemCount()
  for item in range(N_SITES):
    assert sites.IsItemChecked(item), item
  assert sites.GetItemText(item=2, col=0) == "3"
  assert sites.GetItemText(item=2, col=1) == "SE"
  assert sites.GetItemText(item=2, col=5) == "0.87"
  sites.Destroy()

def exercise_save_writes_only_checked_sites(frame, tmp_dir):
  import wx
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  sites.CheckItem(1, False)
  sites.CheckItem(3, False)
  out = os.path.join(tmp_dir, "selected.pdb")
  message_box = _recorder()
  real_message_box = wx.MessageBox
  wx.MessageBox = message_box
  try:
    sites.SaveSites(out)
  finally:
    wx.MessageBox = real_message_box
  symm, serials = _read_serials(out)
  assert serials == [1, 3, 5], serials
  assert symm.unit_cell() is not None
  assert symm.unit_cell().parameters()[0] == 77.324, symm.unit_cell()
  assert len(message_box.calls) == 1
  assert message_box.calls[0][0][0] == "3 atoms saved to %s." % out
  sites.Destroy()

def exercise_save_with_nothing_checked_raises(frame, tmp_dir):
  from libtbx.utils import Sorry
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  for item in range(N_SITES):
    sites.CheckItem(item, False)
  out = os.path.join(tmp_dir, "none.pdb")
  try:
    sites.SaveSites(out)
  except Sorry as e:
    assert "No atoms selected" in str(e), str(e)
  else:
    raise AssertionError("Sorry expected")
  assert not os.path.exists(out)
  sites.Destroy()

def _assert_valid_file_selector_call(recorder):
  """The recorded call must bind against the real wx.FileSelector signature:
  one positional message string, and only keywords the builtin accepts."""
  import wx
  assert len(recorder.calls) == 1, recorder.calls
  args, kwds = recorder.calls[0]
  assert len(args) >= 1 and isinstance(args[0], str) and args[0], args
  allowed = _keyword_names(wx.FileSelector)
  assert "message" in allowed and "flags" in allowed, allowed
  unknown = sorted(set(kwds) - set(allowed[len(args):]))
  assert not unknown, unknown

def exercise_save_button_uses_file_dialog(frame, tmp_dir):
  import wx
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  out = os.path.join(tmp_dir, "from_button.pdb")
  file_selector = _recorder(result=out)
  message_box = _recorder()
  real = (wx.FileSelector, wx.MessageBox)
  wx.FileSelector, wx.MessageBox = file_selector, message_box
  try:
    sites.OnSaveSites(None)
  finally:
    wx.FileSelector, wx.MessageBox = real
  _assert_valid_file_selector_call(file_selector)
  assert file_selector.calls[0][1].get("flags") == wx.FD_SAVE
  symm, serials = _read_serials(out)
  assert serials == [1, 2, 3, 4, 5], serials
  sites.Destroy()

def exercise_save_button_cancelled_writes_nothing(frame, tmp_dir):
  import wx
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  file_selector = _recorder(result="")
  real = wx.FileSelector
  wx.FileSelector = file_selector
  try:
    sites.OnSaveSites(None)
  finally:
    wx.FileSelector = real
  assert len(file_selector.calls) == 1
  assert sorted(os.listdir(tmp_dir)) == ["ha_1.pdb"], os.listdir(tmp_dir)
  sites.Destroy()

def exercise_load_button_uses_file_dialog(frame, tmp_dir):
  import wx
  sites = _make_list(frame)
  file_selector = _recorder(result=_sites_file(tmp_dir))
  real = wx.FileSelector
  wx.FileSelector = file_selector
  try:
    sites.OnLoadSites(None)
  finally:
    wx.FileSelector = real
  _assert_valid_file_selector_call(file_selector)
  assert sites.GetItemCount() == N_SITES, sites.GetItemCount()
  sites.Destroy()

def exercise_change_occupancy_updates_row(frame, tmp_dir):
  import wx
  from wxtbx import heavy_atoms
  sites = _make_list(frame)
  sites.LoadFile(_sites_file(tmp_dir))
  sites.Select(2)
  class _dialog(object):
    seen = None
    def __init__(self, *args, **kwds):
      pass
    def SetOccupancy(self, occ, serial):
      _dialog.seen = (occ, serial)
    def ShowModal(self):
      return wx.ID_OK
    def GetOccupancy(self):
      return 0.5
  real_dialog = heavy_atoms.OccupancyDialog
  heavy_atoms.OccupancyDialog = _dialog
  try:
    sites.ChangeOccupancy()
  finally:
    heavy_atoms.OccupancyDialog = real_dialog
  assert _dialog.seen == (0.87, "3"), _dialog.seen
  assert sites.GetItemText(item=2, col=5) == "0.50"
  assert abs(sites._atoms[2].occ - 0.5) < 1e-6
  sites.Destroy()

def exercise_sites_panel_loads_every_site(frame, tmp_dir):
  from wxtbx.heavy_atoms import sites_panel
  panel = sites_panel(frame, -1)
  panel.load_sites(file_name=_sites_file(tmp_dir))
  assert panel.sites_list.GetItemCount() == N_SITES
  panel.Destroy()

def exercise():
  import wx
  try:
    app = wx.App(False)
  except (Exception, SystemExit) as e:
    # wx.App aborts on a headless box; on GTK this surfaces as a SystemExit
    # (a BaseException), which a plain "except Exception" would not catch.
    print("Skipped: cannot create wx.App ({})".format(e))
    print("OK")
    return
  frame = wx.Frame(None)
  cases = [
    exercise_load_shows_every_site,
    exercise_save_writes_only_checked_sites,
    exercise_save_with_nothing_checked_raises,
    exercise_save_button_uses_file_dialog,
    exercise_save_button_cancelled_writes_nothing,
    exercise_load_button_uses_file_dialog,
    exercise_change_occupancy_updates_row,
    exercise_sites_panel_loads_every_site,
  ]
  failed = []
  for case in cases:
    tmp_dir = tempfile.mkdtemp(prefix="tst_heavy_atoms_")
    try:
      install_wx_deprecation_filters()
      case(frame, tmp_dir)
      print("OK   {}".format(case.__name__))
    except Exception as e:
      failed.append(case.__name__)
      print("FAIL {}: {}: {}".format(case.__name__, type(e).__name__, e))
    finally:
      shutil.rmtree(tmp_dir, ignore_errors=True)
  frame.Destroy()
  if failed:
    print("failed={}".format(len(failed)))
    sys.exit(1)
  print("OK")

if __name__ == "__main__":
  exercise()
