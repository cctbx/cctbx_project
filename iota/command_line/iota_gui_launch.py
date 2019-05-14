from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME iota
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 02/15/2019
Description : IOTA GUI startup module.
'''

import wx
import argparse

from iota import iota_version, help_message
from iota.components import iota_init as init
from iota.components.iota_ui_frames import MainWindow


def parse_command_args():
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog='iota.run',
                                   formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description=help_message,
                                   epilog=('\n{:-^70}\n'.format('')))
  parser.add_argument('path', type=str, nargs='*', default=None,
                      help='Path to data or file with IOTA parameters')
  parser.add_argument('--version', action='version',
                      version='IOTA {}'.format(iota_version),
                      help='Prints version info of IOTA')
  parser.add_argument('-d', '--default', action='store_true',
                      help='Generate default settings files and stop')
  parser.add_argument('--ha14', action='store_true',
                      help='Run IOTA with old HA14 backend')
  parser.add_argument('--random', type=int, nargs=1, default=0,
                      help='Size of randomized subset, e.g. "--random 10"')
  parser.add_argument('--range', type=str, nargs='?', default=None,
                      help='Range of images, e.g."--range 1-5,25,200-250"')
  parser.add_argument('-n', '--nproc', type=int, nargs=1, default=0,
                      help='Specify a number of cores for a multiprocessor run"')
  parser.add_argument('--analyze', type=str, nargs='?', const=None,
                      default=None,
                      help='Use for analysis only; specify run number or folder')
  parser.add_argument('--tmp', type=str, nargs=1, default=None,
                      help='Path to temp folder')
  parser.add_argument('--silent', action='store_true', default=True,
                      help='Run IOTA in silent mode')
  return parser


class MainApp(wx.App):
  """ App to launch the main GUI window  """

  def OnInit(self):
    from platform import python_version
    import matplotlib

    print('Library Versions:')
    print('  Python     : ', python_version())
    print('  wxPython   : ', wx.__version__)
    print('  MatPlotLib : ', matplotlib.__version__)

    # Initialize IOTA instance
    args, phil_args = parse_command_args().parse_known_args()
    input_dict, phil, msg = init.initialize_interface(args, phil_args, True)

    # Initialize Main window
    self.frame = MainWindow(None, -1, title='IOTA v.{}'.format(iota_version),
                            input_dict=input_dict, phil=phil, msg=msg)
    self.frame.place_and_size(set_size=True, set_by='mouse', center=True)

    # Show main window
    self.frame.Show(True)

    self.SetTopWindow(self.frame)
    return True


if __name__ == '__main__':
  app = MainApp(0)
  app.MainLoop()
