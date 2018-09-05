from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.process

'''
Author      : Lyubimov, A.Y.
Created     : 07/26/2014
Last Changed: 08/29/2018
Description : IOTA image processing submission module
'''

import os
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep

import iota.components.iota_image as img
from iota.components.iota_threads import IOTATermination

class ProcessImage():
  ''' Wrapper class to do full processing of an image '''
  def __init__(self, init, input_entry, input_type = 'image'):
    self.init = init
    self.input_entry = input_entry
    self.input_type = input_type


  def run(self):
    if self.input_type == 'image':
      img_object = img.SingleImage(self.input_entry, self.init)
      img_object.import_image()
    elif self.input_type == 'object':
      img_object = self.input_entry[2]
      img_object.import_int_file(self.init)
    else:
      img_object = None

    if self.init.params.image_conversion.convert_only:
      return img_object
    else:
      img_object.process()
      return img_object

class ProcessAll(Thread):
  def __init__(self,
               init,
               iterable,
               input_type='image',
               abort_file=None):
    Thread.__init__(self)
    self.init = ep.load(init)
    self.iterable = ep.load(iterable)
    self.type = input_type
    self.abort_file = abort_file

  def run(self):
    try:
      parallel_map(iterable=self.iterable,
                   func = self.full_proc_wrapper,
                   processes=self.init.params.n_processors)
      end_filename = os.path.join(self.init.tmp_base, 'finish.cfg')
      with open(end_filename, 'w') as ef:
        ef.write('')
    except IOTATermination as e:
      aborted_file = os.path.join(self.init.int_base, '.aborted.tmp')
      with open(aborted_file, 'w') as abtf:
        abtf.write('')
      raise e

  def full_proc_wrapper(self, input_entry):
    abort = os.path.isfile(self.abort_file)
    if abort:
      print 'ABORTING ... NOW!!'
      os.remove(self.abort_file)
      raise IOTATermination('IOTA: Run aborted by user')
    else:
      print 'Processing --- {}'.format(input_entry[2])
      proc_image_instance = ProcessImage(init=self.init,
                                         input_entry=input_entry,
                                         input_type=self.type)
      proc_image_instance.run()


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
  return parser

# ============================================================================ #
if __name__ == "__main__":
  import argparse
  args, unk_args = parse_command_args().parse_known_args()

  proc = ProcessAll(init=args.init, iterable=args.files,
                    input_type=args.type, abort_file=args.stopfile)
  proc.start()
