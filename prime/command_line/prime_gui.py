from __future__ import absolute_import, division, print_function

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
from prime.postrefine.mod_gui_frames import PRIMEWindow

class MainApp(wx.App):
  """ App to launch the main GUI window  """

  def OnInit(self):
    from platform import python_version
    import matplotlib

    print('Library Versions:')
    print('  Python     : ', python_version())
    print('  wxPython   : ', wx.__version__)
    print('  MatPlotLib : ', matplotlib.__version__)

    # Initialize Main window
    self.frame = PRIMEWindow(None, -1, title='PRIME')
    self.frame.place_and_size(set_size='v_default', set_by='mouse', center=True)

    # Show main window
    self.frame.Show(True)

    self.SetTopWindow(self.frame)
    return True


if __name__ == '__main__':
  app = MainApp()
  app.MainLoop()
