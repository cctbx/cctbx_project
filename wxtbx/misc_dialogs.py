
from __future__ import division
from wxtbx.utils import std_sizer_flags, add_ok_cancel_buttons
import wx
import os.path

class ShelxFormatDialog (wx.Dialog) :
  def __init__ (self, file_name, parent=None) :
    wx.Dialog.__init__(self, parent=parent,
      title="SHELX reflection file input")
    outer_sizer = wx.BoxSizer(wx.VERTICAL)
    self.SetSizer(outer_sizer)
    szr = wx.BoxSizer(wx.VERTICAL)
    outer_sizer.Add(szr, 1, wx.EXPAND|wx.ALL, 5)
    txt = wx.StaticText(parent=self,
      label=("The data in this SHELX-format file (%s) may be "+
      "either amplitudes or intensities.  If you know what the data type "+
      "should be, please select from the list below and click 'Okay'.") %
      os.path.basename(file_name),
      size=(480,-1))
    txt.Wrap(480)
    szr.Add(txt, 0, wx.ALL, 5)
    szr2 = wx.BoxSizer(wx.HORIZONTAL)
    szr.Add(szr2)
    label = wx.StaticText(parent=self,
      label="Experimental data type:")
    szr2.Add(label, 0, std_sizer_flags, 5)
    self.data_type_choice = wx.Choice(parent=self,
      size=(200,-1),
      choices=["unknown", "amplitudes", "intensities"])
    szr2.Add(self.data_type_choice, 0, std_sizer_flags, 5)
    add_ok_cancel_buttons(self, szr)
    self.Centre(wx.BOTH)
    outer_sizer.Fit(self)
    self.Fit()

  def GetDataType (self) :
    data_type = self.data_type_choice.GetStringSelection()
    if (data_type == "unknown") :
      data_type = None
    return data_type

def get_shelx_file_data_type (file_name) :
  dlg = ShelxFormatDialog(file_name=file_name)
  data_type = None
  if (dlg.ShowModal() == wx.ID_OK) :
    data_type = dlg.GetDataType()
  wx.CallAfter(dlg.Destroy)
  return data_type
