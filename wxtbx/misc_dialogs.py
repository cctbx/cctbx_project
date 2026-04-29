
"""
Various convenience dialogs (and associated wrapper functions).
"""

from __future__ import absolute_import, division, print_function
from wxtbx.utils import add_ok_cancel_buttons
import wx
from libtbx.utils import Abort
import os.path

class ChoiceDialog(wx.Dialog):
  """
  Dialog containing a single wx.Choice control.  Optionally includes a button
  to flag the selected action as the default for future choices.
  """
  def __init__(self, title, message, choices,
      parent=None,
      choiceLabel=None,
      defaultChoice=-1,
      captionImage=None,
      showDefaultButton=False):
    super(ChoiceDialog, self).__init__(
      parent=parent,
      title=title)
    outer_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(outer_sizer)
    szr = wx.BoxSizer(wx.VERTICAL)
    outer_sizer.Add(szr, 1, wx.EXPAND|wx.ALL, 5)
    if (captionImage is not None):
      assert 0
    else :
      caption = wx.StaticText(parent=self, label=message)
      caption.Wrap(480)
      szr.Add(caption, 0, wx.ALL, 5)
    choice_box = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(choice_box)
    if (choiceLabel is not None):
      label = wx.StaticText(parent=self, label=choiceLabel)
      choice_box.Add(label, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    self.chooser = wx.Choice(parent=self, choices=choices)
    self.chooser.SetMinSize((200,-1))
    choice_box.Add(self.chooser, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
    if (defaultChoice >= 0):
      self.chooser.SetSelection(defaultChoice)
    self.default_btn = None
    if (showDefaultButton):
      self.default_btn = wx.CheckBox(parent=self,
        label="Always do this in the future")
      szr.Add(self.default_btn, 0, wx.ALL, 5)
    add_ok_cancel_buttons(self, szr)
    self.Centre(wx.BOTH)
    outer_sizer.Fit(self)

  def GetChoice(self):
    return self.chooser.GetSelection()

  def GetChoiceString(self):
    return self.chooser.GetStringSelection()

  def GetDefault(self):
    if (self.default_button is not None):
      return self.default_button.GetValue()
    return None

def ChoiceSelector(
    title,
    message,
    choices,
    parent=None,
    choiceLabel=None,
    defaultChoice=0,
    captionImage=None):
  """
  Convenience method for ChoiceDialog (minus default checkbox).  Returns the
  selected string value.
  """
  dlg = ChoiceDialog(
    parent=parent,
    title=title,
    message=message,
    choices=choices,
    defaultChoice=defaultChoice,
    choiceLabel=choiceLabel,
    captionImage=captionImage)
  choice = None
  if (dlg.ShowModal() == wx.ID_OK):
    choice = dlg.GetChoiceString()
  wx.CallAfter(dlg.Destroy)
  if (choice is None):
    raise Abort()
  return choice

def get_shelx_file_data_type(file_name):
  data_type = ChoiceSelector(
    title="SHELX reflection file input",
    message=("The data in this SHELX-format file (%s) may be "+
      "either amplitudes or intensities.  If you know what the data type "+
      "should be, please select from the list below and click 'Okay'.") %
      os.path.basename(file_name),
    choices=["unknown", "amplitudes", "intensities"],
    choiceLabel="Experimental data type:")
  if (data_type == "unknown"):
    data_type = None
  return data_type

if (__name__ == "__main__"):
  app = wx.App(0)
  print(get_shelx_file_data_type("data.hkl"))
