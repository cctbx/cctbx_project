from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 06/01/2016
Description : IOTA GUI startup module.
'''

import wx
from iota.components.iota_gui_init import MainWindow
from iota.components.iota_misc import iota_version

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='IOTA v.{}'.format(iota_version))
    self.frame.Fit()
    self.frame.SetPosition((150, 150))
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
