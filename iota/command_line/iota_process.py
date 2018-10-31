from __future__ import division, print_function, absolute_import
# LIBTBX_SET_DISPATCHER_NAME iota.process

'''
Author      : Lyubimov, A.Y.
Created     : 07/26/2014
Last Changed: 10/30/2018
Description : IOTA image processing submission module
'''

from iota.components.iota_base import ProcessGeneral

class SilentProcess(ProcessGeneral):
  """ Process module customized for 'silent' running (for UI and queueing) """
  def __init__(self, init, iterable, stage, abort_file):
    ProcessGeneral.__init__(self, init=init, iterable=iterable,
                            stage=stage, abort_file=abort_file)

class UIProcess(ProcessGeneral):
  """ Process module customized for 'silent' running (for UI and queueing) """
  def __init__(self, init, iterable, stage, abort_file):
    ProcessGeneral.__init__(self, init=init, iterable=iterable,
                            stage=stage, abort_file=abort_file)

  def callback(self, result):
    """ Will add object file to tmp list for inclusion in info """
    if result:
      with open(self.init.info.obj_list_file, 'a') as olf:
        olf.write('{}\n'.format(result.obj_file))


def parse_command_args():
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog='iota.process')
  parser.add_argument('init', type=str, default=None,
                      help='Path to init file')
  parser.add_argument('--files', type=str, nargs='?', const=None, default=None,
                      help='Specify input file list')
  parser.add_argument('--type', type=str, nargs='?', const=None,
                      default='image',
                      help='Specify input type')
  parser.add_argument('--stopfile', type=str, default=None,
                      help='Path to temporary hidden abort signal file')
  parser.add_argument('--mode', type=str, nargs='?', const=None, default='ui',
                      help='Specify input type')

  return parser

# ============================================================================ #
if __name__ == "__main__":
  import argparse
  args, unk_args = parse_command_args().parse_known_args()

  from libtbx import easy_pickle as ep
  init = ep.load(args.init)
  iterable = ep.load(args.files)

  if args.mode == 'ui':
    proc = UIProcess(init=init, iterable=iterable, stage='all',
                         abort_file=args.stopfile)
  elif args.mode == 'silent':
    proc = SilentProcess(init=init, iterable=iterable, stage='all',
                         abort_file=args.stopfile)
  else:
    proc = None

  if proc is not None:
    proc.start()
