from __future__ import division, print_function, absolute_import
# LIBTBX_SET_DISPATCHER_NAME iota
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 11/05/2018
Description : IOTA GUI startup module.
'''

import wx
from iota import iota_version
from iota.components.iota_ui_init import MainWindow

class MainApp(wx.App):
  ''' App for the main GUI window  '''
  def OnInit(self):
    self.frame = MainWindow(None, -1, title='IOTA v.{}'.format(iota_version))

    # Initialize Main window
    self.frame.place_and_size()
    self.frame.read_command_line_options()

    # Show main window
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
