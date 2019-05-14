from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME iota.single_image
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

'''
Author      : Lyubimov, A.Y.
Created     : 05/31/2018
Last Changed: 01/30/2019
Description : IOTA Single Image: can process single image using DIALS,
with an array of options (i.e. anything from only spotfinding, to indexing,
space group determination, refinement, integration)
'''

import os
import time

from iota.components.iota_init import initialize_single_image
from iota.components.iota_base import ProcessingBase

def parse_command_args():
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog='iota.single_image')
  parser.add_argument('path', type=str, nargs = '?', default = None,
                      help = 'Path to data file')
  # parser.add_argument('--backend', type=str, default='dials',
  #                     help='Backend for processing')
  parser.add_argument('--paramfile', type=str, default=None,
                      help='Parameter file for processing')
  parser.add_argument('--output_file', type=str, default=None,
                      help='Output filename')
  parser.add_argument('--output_dir', type=str, default=None,
                      help='Output directory (for BluIce)')
  parser.add_argument('--termfile', type=str, default='.stop',
                      help='Termination signal filename')
  parser.add_argument('--index', type=int, default=1,
                      help='Numerical index of the image')
  parser.add_argument('--min_bragg', type=int, default=10,
                      help='Minimum spots for successful spotfinding result')
  parser.add_argument('--nproc', type=int, default=1,
                      help='Number of processors')
  parser.add_argument('--action', type=str, default='spotfinding',
                      help='Code for how far to go; available codes: '
                           'spotfind, index, integrate')
  parser.add_argument('--verbose', action = 'store_true',
                      help='Print information to stdout')

  return parser

class SingleImageProcessor(ProcessingBase):
  def __init__(self, *args, **kwargs):
    ProcessingBase.__init__(self, *args, **kwargs)

  def process(self):
    file_wait_start = time.time()

    errors = []
    n_spots = 0
    n_overloads = 0
    res = (99, 99)
    n_rings = 0
    avg_I = 0
    score = 0
    sg = None
    uc = None
    lres = 999
    hres = 999
    img = self.params.input[0]

    img_object = None
    while True:
      elapsed = time.time() - file_wait_start
      if elapsed > 30:
        errors.append('{} does not exist'.format(img))
        print('DEBUG: ELAPSED = ', time.time() - file_wait_start)

        break
      if os.path.isfile(img):
        input_entry = (1, img)
        img_object = self.import_and_process(input_entry)

        n_spots = img_object.final['spots']
        score = img_object.final['indexed']
        hres = img_object.final['res']
        lres = img_object.final['lres']
        sg = img_object.final['sg']
        uc = ' '.join([
          '{:.2f}'.format(img_object.final['a']),
          '{:.2f}'.format(img_object.final['b']),
          '{:.2f}'.format(img_object.final['c']),
          '{:.2f}'.format(img_object.final['alpha']),
          '{:.2f}'.format(img_object.final['beta']),
          '{:.2f}'.format(img_object.final['gamma'])
                     ])
        errors.extend(img_object.errors)
        break

    if img_object:
      if self.verbose:
        print ('SPOTS FOUND: {}'.format(n_spots))
        print ('INDEXING: {} INDEXED SPOTS'.format(score))
        if res[0] != 999:
          print ('RESOLUTION: {:.2f} - {:.2f}'.format(lres, hres))
        if sg and uc:
          print ('BRAVAIS LATTICE: {}'.format(sg))
          print ('UNIT CELL: {}'.format(uc))
        print ('TOTAL PROCESSING TIME: {:.2f} SEC'
               ''.format(time.time() - file_wait_start))

        if errors:
          for e in errors:
            print (e)

      # info = [self.index, len(observed), self.img, sg, uc]

      if self.info.obj_list_file:
        with open(self.info.obj_list_file, 'a') as outf:
          info_line = '{} {} {} {} {}'.format(0, n_spots, img, sg, uc)
          outf.write('{}\n'.format(info_line))

    if self.verbose:
      if errors:
        err = errors[0]
        print_errors = True
      else:
        err = ''
        print_errors = False

      print ('\n__RESULTS__')
      print ('{} {} {} {:.2f} {} {} {} {} {{{}}}' .format(n_spots, n_overloads,
                                      score, res[1], n_rings, 0, avg_I, 0, err))

      if print_errors:
        print ("__ERRORS__")
        for e in errors:
          print (e)


# ============================================================================ #
if __name__ == "__main__":
  import argparse
  args, unk_args = parse_command_args().parse_known_args()

  info, iparams = initialize_single_image(img=os.path.abspath(args.path),
                                          paramfile=args.paramfile,
                                          output_file=args.output_file,
                                          output_dir=args.output_dir,
                                          min_bragg=args.min_bragg)

  interceptor = SingleImageProcessor.for_single_image(info, iparams,
                                                action_code=args.action,
                                                verbose=args.verbose)
  if args.output_dir is not None:
    if not os.path.isdir(args.output_dir):
      os.makedirs(args.output_dir)

  interceptor.start()
