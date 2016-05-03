from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 05/02/2016
Description : IOTA GUI startup module. Version 1.0.004
'''

import wx
from iota.components.iota_gui_init import MainWindow
iota_version = '1.0.004'

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='IOTA v.{}'.format(iota_version))
    self.frame.Fit()
    self.frame.SetPosition((150, 150))
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
