
"""
Miscellaneous custom wxPython windowing objects.
"""

from __future__ import division
import wx

class ChoiceBook (wx.Panel) :
  """
  Notebook-like container with a wx.Choice control instead of tabs.
  """
  def __init__ (self, *args, **kwds) :
    super(ChoiceBook, self).__init__(*args, **kwds)
    szr = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(szr)
    self._pages = []
    self._current_page = None
    self._chooser = wx.Choice(parent=self, )#size=(200,-1))
    szr.Add(self._chooser, 0, wx.ALL|wx.ALIGN_CENTER, 5)
    self._page_sizer = wx.BoxSizer(wx.VERTICAL)
    szr.Add(self._page_sizer, 1, wx.EXPAND)
    self.Bind(wx.EVT_SHOW, self.OnShow, self)
    self.Bind(wx.EVT_CHOICE, self.OnChoose, self._chooser)

  def AddPage (self, page, label) :
    """
    Add a panel to the notebook (and the label to the wx.Choice).
    """
    page.Hide()
    self._pages.append(page)
    self.GetSizer().Add(page, 1, wx.ALL|wx.EXPAND, 0)
    items = list(self._chooser.GetItems())
    items.append(label)
    self._chooser.SetItems(items)

  def InsertPage (self, index, page, label) :
    """
    Insert a panel in the notebook (and the label to the wx.Choice) at the
    specified index.
    """
    page.Hide()
    self._pages.insert(0, page)
    self.GetSizer().Add(page, 1, wx.ALL|wx.EXPAND, 0)
    items = list(self._chooser.GetItems())
    items.insert(0, label)
    self._chooser.SetItems(items)

  def GetPage (self, i_page) :
    return self._pages[i_page]

  def SetPage (self, i_page) :
    if (self._current_page is not None) :
      self._current_page.Hide()
    self._current_page = self.GetPage(i_page)
    self._current_page.Layout()
    if hasattr(self._current_page, "SetupScrolling") : # scrolledpanel hack
      self._current_page.SetupScrolling(scrollToTop=False)
    self._current_page.Show()
    self._chooser.SetSelection(i_page)
    self.Layout()

  def SetSelection (self, i_page) :
    self.SetPage(i_page)

  def GetPageIndex (self, page) :
    for i_page, other in enumerate(self._pages) :
      if (other is page) :
        return i_page
    return -1

  def GetPageCount (self) :
    return len(self._pages)

  def OnShow (self, evt) :
    self._chooser.Layout()
    self.GetSizer().Layout()
    if (self.GetPageCount() > 0) and (self._current_page is None) :
      self.SetPage(0)

  def OnChoose (self, evt) :
    i_page = self._chooser.GetSelection()
    self.SetPage(i_page)
