from __future__ import division, print_function, absolute_import
# LIBTBX_SET_DISPATCHER_NAME prime
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 05/19/2016
Last Changed: 06/07/2019
Description : PRIME GUI startup module.
'''

import wx
from prime.postrefine.mod_gui_frames import PRIMEWindow as MainWindow

class MainApp(wx.App):
  ''' App for the main GUI window  '''

  def OnInit(self):

    from platform import python_version
    import matplotlib

    print('Library Versions:')
    print('  Python     : ', python_version())
    print('  wxPython   : ', wx.__version__)
    print('  MatPlotLib : ', matplotlib.__version__)

    self.frame = MainWindow(parent=None, title='PRIME')
    self.frame.place_and_size(center=True, set_by='mouse')
    self.frame.Show(True)
    self.SetTopWindow(self.frame)
    return True


if __name__ == '__main__':
  app = MainApp(useBestVisual=True)
  app.MainLoop()
