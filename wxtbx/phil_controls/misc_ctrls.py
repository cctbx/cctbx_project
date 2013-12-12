
from __future__ import division
from wxtbx.phil_controls import strings
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

if (__name__ == "__main__") :
  app = wx.App(0)
  frame = wx.Frame(None, -1, "Strings control test")
  panel = wx.Panel(frame, -1, size=(600,400))
  txt1 = wx.StaticText(panel, -1, "Elements:", pos=(100,180))
  elems_ctrl = ElementsCtrl(panel, -1, pos=(300,180), size=(120,-1),
      name="Elements", style=wx.TE_PROCESS_ENTER)
  def OnUpdate (evt) :
    elems = elems_ctrl.GetPhilValue()
    print "Current elements:", elems
  def OnOkay (evt) :
    print "elems:", elems_ctrl.GetPhilValue()
    print "elems phil:", elems_ctrl.GetStringValue()
  frame.Bind(phil_controls.EVT_PHIL_CONTROL, OnUpdate, elems_ctrl)
  btn = wx.Button(panel, -1, "Process input", pos=(400, 360))
  frame.Bind(wx.EVT_BUTTON, OnOkay, btn)
  frame.Fit()
  frame.Show()
  app.MainLoop()
