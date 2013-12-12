
from __future__ import division
from wxtbx.phil_controls import strings
from wxtbx.phil_controls import strctrl
from wxtbx import phil_controls
import wx

WXTBX_ELEMENTS_CTRL_ALLOW_AX = 1
WXTBX_ELEMENTS_CTRL_ALLOW_CLUSTERS = 2

clusters = ['TX']

class ElementsCtrl (strings.StringsCtrl) :
  def __init__ (self, *args, **kwds) :
    elements_style = kwds.get('style', 0)
    kwds['style'] = 0
    strings.StringsCtrl.__init__(self, *args, **kwds)
    self.SetWxtbxStyle(elements_style)

  def CreateValidator (self) :
    return ElementsValidator()

class ElementsValidator (strings.StringsValidator) :
  def CheckFormat (self, value) :
    from cctbx.eltbx import chemical_elements
    allowed = chemical_elements.proper_upper_list()
    ctrl = self.GetWindow()
    style = ctrl.GetWxtbxStyle()
    elem_strings = strings.parse_strings(value)
    for elem in elem_strings :
      elem = elem.upper()
      if (elem == "AX") and (style & WXTBX_ELEMENTS_CTRL_ALLOW_AX) :
        pass
      elif (elem in clusters) and (style & WXTBX_ELEMENTS_CTRL_ALLOW_CLUSTERS):
        pass
      elif (not elem in allowed) :
        raise ValueError("'%s' is not a valid element symbol." % elem)
    return elem_strings

# TODO pop-up help
class AtomSelectionCtrl (strctrl.StrCtrl) :
  def __init__ (self, *args, **kwds) :
    size = kwds.get('size', None)
    if (size is None) :
      kwds['size'] = (400, 40)
    strctrl.StrCtrl.__init__(self, *args, **kwds)

  def CreateValidator (self) :
    return AtomSelectionValidator()

class AtomSelectionValidator (strctrl.StrValidator) :
  def CheckFormat (self, value) :
    import iotbx.pdb.hierarchy
    root = iotbx.pdb.hierarchy.root()
    sel_cache = root.atom_selection_cache()
    atom_sel = sel_cache.selection(value)
    return value

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Strings control test")
  panel = wx.Panel(frame, -1, size=(720,400))
  txt1 = wx.StaticText(panel, -1, "Elements:", pos=(20,180))
  elems_ctrl = ElementsCtrl(panel, -1, pos=(200,180), size=(120,-1),
      name="Elements", style=wx.TE_PROCESS_ENTER)
  txt2 = wx.StaticText(panel, -1, "Atom selection:", pos=(20,240))
  sel_ctrl = AtomSelectionCtrl(panel, -1, pos=(200,240), name="Selection")
  def OnUpdate (evt) :
    elems = elems_ctrl.GetPhilValue()
    print "Current elements:", elems
    sel = sel_ctrl.GetPhilValue()
    print "Atom selection:", sel
  def OnOkay (evt) :
    print "elems:", elems_ctrl.GetPhilValue()
    print "elems phil:", elems_ctrl.GetStringValue()
    print "selection:", sel_ctrl.GetPhilValue()
    print "selection phil:", sel_ctrl.GetStringValue()
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, elems_ctrl)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, sel_ctrl)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
