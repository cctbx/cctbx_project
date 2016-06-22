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
from xfel.ui.components.xfel_gui_dialogs import SettingsDialog

class MainApp(wx.App):
  ''' App for the main GUI window  '''

  def OnInit(self):

    self.frame = MainWindow(None, -1, title='CCTBX.XFEL')

    # select primary display and center on that
    display = wx.Display(0)
    x, y, w, h = display.GetGeometry()
    fw, fh = self.frame.GetSize()
    frame_position = (x + w / 2 - fw / 2, y + h / 2 - fh / 2)
    self.frame.SetMinSize(self.frame.GetEffectiveMinSize())
    self.frame.SetPosition(frame_position)
    #self.frame.Center()

    # Start with login dialog before opening main window
    self.login = SettingsDialog(self.frame, self.frame.params)
    self.login.SetTitle('CCTBX.XFEL Login')
    self.login.Center()
    if (self.login.ShowModal() == wx.ID_OK):
      if self.frame.connect_to_db(drop_tables=self.login.drop_tables):
        self.exp_tag = '| {}'.format(self.login.db_cred.ctr.GetValue())
        self.exp = '| {}'.format(self.login.experiment.ctr.GetValue())
        self.frame.SetTitle('CCTBX.XFEL {} {}'.format(self.exp, self.exp_tag))
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        self.frame.start_run_sentinel()
        return True
      else:
        return False
    else:
      return False

if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
