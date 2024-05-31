from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 06/02/2016
Last Changed: 06/02/2016
Description : XFEL UI startup module.
'''

import wx, sys
import matplotlib as mp
mp.use('PS')


from xfel.ui.components.xfel_gui_init import MainWindow
from xfel.ui.components.xfel_gui_dialogs import SettingsDialog
from xfel.ui import load_cached_settings, save_cached_settings

class MainApp(wx.App):
  ''' App for the main GUI window  '''

  def OnInit(self):
    params = load_cached_settings()
    self.frame = MainWindow(None, -1, title='CCTBX.XFEL', params=params)

    # select primary display and center on that
    self.frame.SetSize((800, -1))
    minx, miny = self.frame.GetEffectiveMinSize()
    dispx, dispy =  wx.GetDisplaySize()
    self.frame.SetMinSize((min(minx, dispx), min(miny, dispy)))
    self.frame.Center()

    # Start with login dialog before opening main window
    self.login = SettingsDialog(self.frame, params=params)
    self.login.SetTitle('CCTBX.XFEL Login')
    self.login.Center()
    if (self.login.ShowModal() == wx.ID_OK):
      save_cached_settings(params)
      if self.frame.connect_to_db(drop_tables=self.login.drop_tables):
        self.exp_tag = '| {}'.format(self.login.db_cred.ctr.GetValue())
        self.exp = '| {}'.format(self.login.experiment.ctr.GetValue())
        self.frame.SetTitle('CCTBX.XFEL {} {}'.format(self.exp, self.exp_tag))
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        self.frame.onTabChange(None)
        #self.frame.start_run_sentinel()
        #self.frame.start_job_monitor()
        #self.frame.start_prg_sentinel()
        self.frame.run_window.show_hide_tabs()
        return True
      else:
        return False
    else:
      return False

def run(args):
  import libtbx
  libtbx.mpi_import_guard.disable_mpi = True

  if '-h' in args or '--help' in args or '-c' in args:
    from xfel.ui import master_phil_str
    from libtbx.phil import parse
    print("""
Runs the cctbx.xfel GUI. For more information see:
http://cci.lbl.gov/publications/download/CCN_2019_p22_Brewster.pdf

Configuration options:

""")
    parse(master_phil_str, process_includes = True).show(attributes_level=2)
    return
  app = MainApp(0)
  app.MainLoop()

if __name__ == '__main__':
  run(sys.argv[1:])
