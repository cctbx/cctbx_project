from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME simtbx
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 12/12/2017
Last Changed: 12/12/2017
Description : SIMTBX (nanoBragg) GUI startup module.
'''

import wx
from simtbx.nanoBragg.nanoBragg_gui_init import MainWindow

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='SIMTBX (nanoBragg)')
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())

    # Find mouse position and open window there
    mx, my = wx.GetMousePosition()
    self.frame.SetPosition((mx, my))

    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
