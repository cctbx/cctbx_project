
from wx.lib.agw import customtreectrl
import wx
from libtbx.utils import Sorry
import re

class PhilTreeCtrl (customtreectrl.CustomTreeCtrl) :
  def __init__ (self, *args, **kwds) :
    kwds = dict(kwds)
    kwds['agwStyle'] = wx.TR_HAS_VARIABLE_ROW_HEIGHT|wx.TR_HAS_BUTTONS|wx.TR_TWIST_BUTTONS|wx.TR_HIDE_ROOT
    customtreectrl.CustomTreeCtrl.__init__(self, *args, **kwds)
    self.AddRoot("(parameters)")
    self._nodes_lookup_short = {}
    self._nodes_lookup_full = {}
    self._node_names_ordered = []
    self.ClearSearch()

  def DeleteAllItems (self) :
    customtreectrl.CustomTreeCtrl.DeleteAllItems(self)
    self._nodes_lookup_short = {}
    self._nodes_lookup_full = {}
    self._node_names_ordered = []

  def ClearSearch (self) :
    self._search_results = []
    self._current_search_text = None
    self._current_search_item = -1

  def SaveNode (self, node, phil_object) :
    name = phil_object.name
    full_path = phil_object.full_path
    self._node_names_ordered.append((name, full_path))
    node.SetData(phil_object)
    if (name in self._nodes_lookup_short) :
      self._nodes_lookup_short[name].append(node)
    else :
      self._nodes_lookup_short[name] = [node]
    if (full_path in self._nodes_lookup_full) :
      self._nodes_lookup_full[full_path].append(node)
    else :
      self._nodes_lookup_full[full_path] = [node]

  def DrawPhilObject (self, phil_root) :
    assert (phil_root.name == "")
    for phil_object in phil_root.objects :
      self._DrawPhilObject(phil_object, self.GetRootItem())

  def _DrawPhilObject (self, phil_object, current_node) :
    if (phil_object.is_definition) :
      node = self.AppendItem(current_node, phil_object.as_str().strip())
      self.SaveNode(node, phil_object)
    else :
      new_node = self.AppendItem(current_node, phil_object.name)
      self.SaveNode(new_node, phil_object)
      for object_ in phil_object.objects :
        self._DrawPhilObject(object_, new_node)
      self.Expand(new_node)

  def SearchItems (self, search_text, partial=False) :
    if (search_text != self._current_search_text) :
      results = []
      if ("." in search_text) :
        for name, path in self._node_names_ordered :
          if (partial and (search_text in path)) or (path == search_text) :
            results.extend(self._nodes_lookup_full[path])
      else :
        for name, path in self._node_names_ordered :
          if (partial and (search_text in name)) or (name == search_text) :
            results.extend(self._nodes_lookup_short[name])
      self._search_results = results
    self.OnNext(None)
    return len(self._search_results)

  def OnNext (self, event) :
    if (len(self._search_results) == 0) :
      return
    self._current_search_item += 1
    if (self._current_search_item == len(self._search_results)) :
      self._current_search_item = 0
    node = self._search_results[self._current_search_item]
    self.SelectItem(node)
    self.ScrollTo(node)

valid_text = re.compile("^[a-zA-Z]{1,}[a-zA-z0-9_]*$")
valid_text_partial = re.compile("^[a-zA-Z0-9_]*$")

class PhilTreeFrame (wx.Frame) :
  def __init__ (self, *args, **kwds) :
    wx.Frame.__init__(self, *args, **kwds)
    self.panel = wx.Panel(self)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.panel.SetSizer(szr)
    #szr.Add(self.panel, 1, wx.EXPAND)
    self.tree = PhilTreeCtrl(self.panel, -1, size=(600,400),
      style=wx.SUNKEN_BORDER)
    txt1 = wx.StaticText(self.panel, -1, "Search:")
    szr2 = wx.FlexGridSizer(cols=2)
    szr.Add(szr2)
    szr.Add(self.tree, 1, wx.EXPAND, 2)
    szr2.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr3 = wx.BoxSizer(wx.HORIZONTAL)
    szr2.Add(szr3, 0, wx.ALIGN_CENTER_VERTICAL)
    search_box = wx.SearchCtrl(self.panel, style=wx.TE_PROCESS_ENTER,
      size=(160,-1))
    search_box.ShowCancelButton(True)
    szr3.Add(search_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.Bind(wx.EVT_TEXT_ENTER, self.OnSearch, search_box)
    self.Bind(wx.EVT_SEARCHCTRL_SEARCH_BTN, self.OnSearch, search_box)
    self.Bind(wx.EVT_SEARCHCTRL_CANCEL_BTN, self.OnCancel, search_box)
    self.partial_box = wx.CheckBox(self.panel, -1, "Include partial matches")
    #szr3 = wx.BoxSizer(wx.HORIZONTAL)
    szr3.Add(self.partial_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.search_result = wx.StaticText(self.panel, -1, "", size=(300,-1))
    szr2.Add((1,1))
    szr2.Add(self.search_result, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    szr.Layout()
    szr.Fit(self.panel)
    self.Fit()

  def DrawPhilTree (self, phil_object) :
    self.tree.DrawPhilObject(phil_object)
    self.Refresh()

  def OnSearch (self, event) :
    search_text = event.GetEventObject().GetValue()
    partial = self.partial_box.GetValue()
    if ((partial and (not valid_text_partial.match(search_text))) or
        (not partial and (not valid_text.match(search_text)))) :
      self.search_result.SetLabel("Invalid search text!")
      self.search_result.SetForegroundColour((200,0,0))
      raise Sorry("Invalid search text - only alphanumeric characters ("+
        "including) underscore are allowed.  If partial matches are not "+
        "included, the search text must also begin with a letter.")
    n_items = self.tree.SearchItems(search_text, partial)
    if (n_items == 0) :
      self.search_result.SetForegroundColour((200,0,0))
    else :
      self.search_result.SetForegroundColour((0,0,0))
    self.search_result.SetLabel("%d items found" % n_items)

  def OnCancel (self, event) :
    event.GetEventObject().Clear()
    self.tree.ClearSearch()

if (__name__ == "__main__") :
  from mmtbx.command_line import fmodel
  app = wx.App(0)
  frame  = PhilTreeFrame(None, -1, "Phenix settings")
  frame.DrawPhilTree(fmodel.fmodel_from_xray_structure_master_params)
  frame.Show()
  app.MainLoop()
