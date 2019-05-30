from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 05/19/2016
Last Changed: 06/01/2016
Description : PRIME GUI startup module.
'''

import wx
import numpy as np
assert np
from prime.postrefine.mod_gui_init import PRIMEWindow

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = PRIMEWindow(None, -1, title='PRIME')
    self.frame.place_and_size(set_size=True, center=True, set_by='mouse')
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
