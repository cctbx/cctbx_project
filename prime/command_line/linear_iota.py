from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 07/29/2015
Description : IOTA command-line module for running modules in order.
              Version 2.00
'''

iota_version = '2.00'
help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

Auto mode
Usage: prime.iota [OPTIONS] path/to/raw/images
Generates two files, parameter file for IOTA (iota.param) and
target file for cctbx.xfel (target.phil). Integrates a random
subset of images without target cell. Outputs basic analysis.
Converts raw images into pickle format and crops to ensure that
beam center is in center of image.

Script mode
Usage: prime.iota [OPTIONS] <script>.param
Run using IOTA parameter file and target PHIL file generated from
the dry run or auto mode. Make sure that IOTA parameter file has
the path to the input image folder under "input". Converts raw
images into pickle format and modifies by cropping or padding to
ensure that beam center is in center of image. Can also blank out
beam stop shadow.

"""

from libtbx.easy_mp import parallel_map
from prime.iota.iota_init import InitAll
from prime.iota.iota_analysis import Analyzer
import prime.iota.iota_image as img
import prime.iota.iota_cmd as cmd
import prime.iota.iota_misc as misc


def run_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  prog_count = input_entry[0]
  n_img = input_entry[1]
  gs_prog = cmd.ProgressBar(title='PROCESSING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()
  return run_one_image(input_entry, init)

def run_one_image(image, init):
    # Import image
    single_image = img.SingleImage(image, init)
    img_object = single_image.import_image()

    # Check / convert / triage image
    img_object = single_image.convert_image()

    # Exit if the image does not have diffraction
    if img_object.triage == 'rejected':
      misc.main_log(init.logfile, 'Image does not have diffraction!')
      return img_object

    # Grid search
    img_object = single_image.integrate('grid search')

    # Selection
    img_object = single_image.select()

    # Exit image not integrated
    if img_object.final['final'] == None:
      misc.main_log(init.logfile, 'Image not integrated!')
      return img_object

    # Final integration
    img_object = single_image.integrate('integrate')

    return img_object


# ============================================================================ #
if __name__ == "__main__":

  # Initialize IOTA parameters and log
  init = InitAll(iota_version, help_message)
  init.run()

  # Run all modules in order in multiprocessor mode
  cmd.Command.start("Processing {} images".format(len(init.input_list)))
  img_list = [[i, len(init.input_list) + 1, j] for i, j in enumerate(init.input_list, 1)]
  img_objects = parallel_map(iterable  = img_list,
                             func      = run_wrapper,
                             processes = init.params.n_processors)
  cmd.Command.end("Processing {} images -- DONE ".format(len(init.input_list)))

  final_objects = [i for i in img_objects if i.triage == 'accepted' and\
                                             i.prefilter == True and\
                                             i.final['final'] != None]

  # Analysis of integration results
  analysis = Analyzer(final_objects, init.logfile, iota_version, init.now)
  analysis.print_results()
  analysis.unit_cell_analysis(init.params.advanced.cluster_threshold,
                              init.int_base)
  analysis.print_summary(init.int_base)
  analysis.make_prime_input(init.int_base)

  misc.iota_exit(iota_version)
