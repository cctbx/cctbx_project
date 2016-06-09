from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 06/02/2016
Last Changed: 06/02/2016
Description : XFEL UI startup module.
'''

import wx
from xfel.ui.components.xfel_gui_init import MainWindow

# TODO: get actual values from initial settings window
exp = 'cxid9114'
exp_tag = 'aaron_debug'

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='CCTBX.XFEL : {} : {}'
                                            ''.format(exp, exp_tag))
    self.frame.Center()
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
