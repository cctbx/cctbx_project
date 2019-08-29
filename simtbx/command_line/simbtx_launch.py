from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME simtbx
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 12/12/2017
Last Changed: 10/30/2018
Description : SIMTBX (nanoBragg) GUI startup module.
'''

import wx
from simtbx.nanoBragg.nanoBragg_gui_init import MainWindow

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='SIMTBX (nanoBragg)')

    # Find mouse position and open window there
    mx, my = wx.GetMousePosition()
    self.frame.SetPosition((mx, my))

    # Get display index and geometry
    n_display = wx.Display.GetFromPoint((mx, my))
    display = wx.Display(index=n_display)
    geom = display.GetGeometry()

    # Set main frame size to fraction of display geometry and center window
    self.frame.SetSize((0.5*geom[2], 0.75*geom[3]))
    self.frame.Center()

    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
