
from __future__ import absolute_import, division, print_function
from wxtbx.phil_controls import strings
from wxtbx.phil_controls import strctrl
from wxtbx import phil_controls
import wx

WXTBX_ELEMENTS_CTRL_ALLOW_AX = 1
WXTBX_ELEMENTS_CTRL_ALLOW_CLUSTERS = 2

clusters = set(['TX', 'XX'])

class ElementCtrl(strctrl.StrCtrl):
  """"
  Control for entering a single element symbol or scattering type, using
  phil str type
  """
  def __init__(self, *args, **kwds):
    elements_style = kwds.get('style', 0)
    kwds['style'] = 0
    strctrl.StrCtrl.__init__(self, *args, **kwds)
    self.SetWxtbxStyle(elements_style)

  def CreateValidator(self):
    return SingleElementValidator()

class ElementsCtrl(strings.StringsCtrl):
  """
  Control for entering a list of elements (comma- or space-separated) or
  scattering types, using phil strings type
  """
  def __init__(self, *args, **kwds):
    elements_style = kwds.get('style', 0)
    kwds['style'] = 0
    strings.StringsCtrl.__init__(self, *args, **kwds)
    self.SetWxtbxStyle(elements_style)

  def CreateValidator(self):
    return ElementsValidator()

class ElementsValidator(strings.StringsValidator):
  single_element = False
  def CheckFormat(self, value):
    from cctbx.eltbx import chemical_elements
    allowed = chemical_elements.proper_upper_list() + list(clusters)
    ctrl = self.GetWindow()
    style = ctrl.GetWxtbxStyle()
    elem_strings = strings.parse_strings(value)
    if (self.single_element) and (len(elem_strings) > 1):
      raise ValueError("Only a single element may be specified here!")
    for elem in elem_strings :
      elem = elem.upper() # used in Phaser EP
      if (elem == "AX" or elem == "RX") and (style & WXTBX_ELEMENTS_CTRL_ALLOW_AX):
        pass
      elif (elem in clusters) and (style & WXTBX_ELEMENTS_CTRL_ALLOW_CLUSTERS):
        pass
      elif (not elem in allowed):
        raise ValueError("'%s' is not a valid element symbol." % elem)
    if (self.single_element):
      return elem_strings[0]
    else :
      return elem_strings

class SingleElementValidator(ElementsValidator):
  single_element = True

# TODO pop-up help
class AtomSelectionCtrl(strctrl.StrCtrl):
  """Control for editing an atom selection (phil atom_selection type)"""
  def __init__(self, *args, **kwds):
    size = kwds.get('size', None)
    if (size is None):
      kwds['size'] = (400, 40)
    strctrl.StrCtrl.__init__(self, *args, **kwds)

  def CreateValidator(self):
    return AtomSelectionValidator()

class AtomSelectionValidator(strctrl.StrValidator):
  def CheckFormat(self, value):
    # FIXME this breaks on mmtbx-supported selections!
    #import iotbx.pdb.hierarchy
    #root = iotbx.pdb.hierarchy.root()
    #sel_cache = root.atom_selection_cache()
    #atom_sel = sel_cache.selection(value)
    return value

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Strings control test")
  panel = wx.Panel(frame, -1, size=(720,400))
  txt0 = wx.StaticText(panel, -1, "Scatterer:", pos=(20,80))
  sc_ctrl = ElementCtrl(panel, -1, pos=(200,80), size=(120,-1),
      name="Scatterer",
      style=WXTBX_ELEMENTS_CTRL_ALLOW_AX|WXTBX_ELEMENTS_CTRL_ALLOW_CLUSTERS)
  txt1 = wx.StaticText(panel, -1, "Elements:", pos=(20,180))
  elems_ctrl = ElementsCtrl(panel, -1, pos=(200,180), size=(120,-1),
      name="Elements", style=wx.TE_PROCESS_ENTER)
  txt2 = wx.StaticText(panel, -1, "Atom selection:", pos=(20,240))
  sel_ctrl = AtomSelectionCtrl(panel, -1, pos=(200,240), name="Selection")
  def OnUpdate(evt):
    print("OnUpdate:")
    print("  current scatterer:", sc_ctrl.GetPhilValue())
    elems = elems_ctrl.GetPhilValue()
    print("  current elements:", elems)
    sel = sel_ctrl.GetPhilValue()
    print("  atom selection:", sel)
  def OnOkay(evt):
    print("OnOkay:")
    print("  scatterer:", sc_ctrl.GetPhilValue())
    print("  scatterer phil:", sc_ctrl.GetStringValue())
    print("  elems:", elems_ctrl.GetPhilValue())
    print("  elems phil:", elems_ctrl.GetStringValue())
    print("  selection:", sel_ctrl.GetPhilValue())
    print("  selection phil:", sel_ctrl.GetStringValue())
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, sc_ctrl)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, elems_ctrl)
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, sel_ctrl)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
