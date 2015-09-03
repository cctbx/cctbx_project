from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 09/02/2015
Description : IOTA command-line module. Version 2.11
'''

iota_version = '2.11'
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

MPI mode
Usage:
prime.linear_iota iota.param --mpi init
prime.linear_iota iota.param --mpi process
prime.linear_iota iota.param --mpi analyze

Run IOTA in three separate batches (initialization, image processing,
analysis); can use MPI (mpirun) to run the image processing step.
Can run these in sequence in a shell script or any other kind of a
submission script. Useful for huge datasets.


"""
from prime.iota.iota_analysis import Analyzer
from prime.iota.iota_init import InitAll
import prime.iota.iota_image as img
import prime.iota.iota_cmd as cmd
import prime.iota.iota_misc as misc
from libtbx.easy_mp import parallel_map

def importer_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  prog_count = input_entry[0]
  n_img = input_entry[1]
  gs_prog = cmd.ProgressBar(title='IMPORTING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()
  img_object = img.SingleImage(input_entry, init)
  return img_object.import_image()

def gs_importer_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  prog_count = input_entry[0]
  n_img = input_entry[1]
  gs_prog = cmd.ProgressBar(title='IMPORTING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()
  img_object = input_entry[2]
  return img_object.import_int_file(init)

def conversion_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  prog_count = input_entry[0]
  n_img = input_entry[1]
  img_object = input_entry[2]
  gs_prog = cmd.ProgressBar(title='CONVERTING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()
  return img_object.convert_image()

def processing_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  prog_count = input_entry[0]
  n_img = input_entry[1]
  img_object = input_entry[2]
  gs_prog = cmd.ProgressBar(title='PROCESSING')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()
  return img_object.process()

# ============================================================================ #
if __name__ == "__main__":

  # Initialize IOTA parameters and log
  init = InitAll(iota_version, help_message)
  init.run()

  if init.params.selection.select_only.flag_on:
    # Generate image objects and modify with saved grid search results
    cmd.Command.start("Generating {} image objects".format(len(init.gs_img_objects)))
    img_list = [[i, len(init.gs_img_objects) + 1, j] for i, j in enumerate(init.gs_img_objects, 1)]
    img_objects = parallel_map(iterable  = img_list,
                               func      = gs_importer_wrapper,
                               processes = init.params.n_processors)
    cmd.Command.end("Generating {} image objects -- DONE ".format(len(init.gs_img_objects)))

  else:
    # Import and process raw images or image pickles
    # Make list of image objects
    cmd.Command.start("Importing {} images".format(len(init.input_list)))
    img_list = [[i, len(init.input_list) + 1, j] for i, j in enumerate(init.input_list, 1)]
    img_objects = parallel_map(iterable  = img_list,
                               func      = importer_wrapper,
                               processes = init.params.n_processors)
    cmd.Command.end("Importing {} images -- DONE ".format(len(init.input_list)))

    # Check / convert / triage images
    cmd.Command.start("Checking / converting {} images".format(len(img_objects)))
    img_list = [[i, len(img_objects) + 1, j] for i, j in enumerate(img_objects, 1)]
    img_objects = parallel_map(iterable  = img_list,
                               func      = conversion_wrapper,
                               processes = init.params.n_processors)

    # Remove rejected images from image object list
    acc_img_objects = [i.fail for i in img_objects if i.fail == None]
    cmd.Command.end("Accepted {} of {} images -- DONE "\
                    "".format(len(acc_img_objects), len(img_objects)))

    # Exit if none of the images have diffraction
    if len(acc_img_objects) == 0:
      misc.main_log(init.logfile, 'No images have diffraction!', True)
      misc.iota_exit(iota_version)
    else:
      misc.main_log(init.logfile, "{} out of {} images have diffraction"\
                                         "".format(len(acc_img_objects),
                                                   len(img_objects)))

    # Check for -c option and exit if true
    if init.params.image_conversion.convert_only:
      misc.iota_exit(iota_version)

  cmd.Command.start("Processing {} images".format(len(img_objects)))
  img_list = [[i, len(img_objects) + 1, j] for i, j in enumerate(img_objects, 1)]
  img_objects = parallel_map(iterable  = img_list,
                             func      = processing_wrapper,
                             processes = init.params.n_processors)
  cmd.Command.end("Processing {} images -- DONE ".format(len(img_objects)))

  # Analysis of integration results
  analysis = Analyzer(img_objects, init.logfile, iota_version, init.now)
  analysis.print_results()
  analysis.unit_cell_analysis(init.params.analysis.cluster_threshold,
                              init.int_base)
  analysis.print_summary(init.int_base)
  analysis.make_prime_input(init.int_base)

  # Spotfinding heatmap
  if init.params.analysis.heatmap != None:
    hm_file = "{}/{}".format(init.viz_base, 'heatmap.png')
    if init.params.analysis.heatmap == 'show':
      analysis.show_heatmap()
    elif init.params.analysis.heatmap == 'file':
      analysis.show_heatmap(show=False, hm_file=hm_file)
    elif init.params.analysis.heatmap == 'both':
      analysis.show_heatmap(hm_file=hm_file)

  misc.iota_exit(iota_version)
