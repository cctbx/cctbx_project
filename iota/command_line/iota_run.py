from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.run

from builtins import str
from builtins import object
'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 02/22/2017
Description : IOTA command-line module.
'''
from __future__ import print_function

from iota.components.iota_analysis import Analyzer
from iota.components.iota_init import InitAll
import iota.components.iota_image as img
import iota.components.iota_cmd as cmd
import iota.components.iota_misc as misc
from libtbx.easy_mp import parallel_map

iota_version = misc.iota_version

help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

Auto mode
Usage: iota.run [OPTIONS] path/to/raw/images
Generates two files, parameter file for IOTA (iota.param) and
target file for cctbx.xfel (target.phil). Integrates a random
subset of images without target cell. Outputs basic analysis.
Converts raw images into pickle format and crops to ensure that
beam center is in center of image.

Single-image mode
Usage: iota.run [OPTIONS] /path/to/single/image.[pickle/cbf/img/mccd]
Same as AUTO mode, but can accept a single image file. Will generate default
iota.param and target.phil files and integrate the single image provided. Can
also be used with iota.single_image in the same manner.

Script mode
Usage: iota.run [OPTIONS] <script>.param
Run using IOTA parameter file and target PHIL file generated from
the dry run or auto mode. Make sure that IOTA parameter file has
the path to the input image folder under "input". Converts raw
images into pickle format and modifies by cropping or padding to
ensure that beam center is in center of image. Can also blank out
beam stop shadow.

"""



class XTermIOTA(object):
  ''' Main class that will initalize and run everything'''

  def __init__(self):
    self.prog_count = 0
    self.init = InitAll(help_message)
    self.init.run()

    self.full = self.init.args.full

  def proc_wrapper(self, input_entry):
    ''' Wrapper for processing function using the image object
    @param input_entry: [image_number, total_images, image_object]
    @return: image object
    '''
    try:
      if self.stage == 'import':
        if self.init.params.cctbx.selection.select_only.flag_on:
          img_object = input_entry[2]
          img_object.import_int_file(self.init)
        else:
          img_object = img.SingleImage(input_entry, self.init)
          img_object.import_image()
      elif self.stage == 'process':
        img_object = input_entry[2]
        img_object.process()
      elif self.stage == 'all':
        if self.init.params.cctbx.selection.select_only.flag_on:
          img_object = input_entry[2]
          img_object.import_int_file(self.init)
        else:
          img_object = img.SingleImage(input_entry, self.init)
          img_object.import_image()
        img_object.process()
    except Exception as e:
      pass

    return img_object

  def callback(self, result):
    ''' To be run on completion of each step
    @param result: image_object
    @return:
    '''
    if self.prog_count < len(self.init.input_list):
      prog_step = 100 / len(self.init.input_list)
      self.gs_prog.update(self.prog_count * prog_step, self.prog_count)
      self.prog_count += 1
    else:
      self.gs_prog.finished()

  def run_all(self):
    ''' Run the full processing in multiprocessing mode '''

    # Determine whether reading in image objects or images, create image list
    if self.init.params.cctbx.selection.select_only.flag_on:
      self.img_list = [[i, len(self.init.gs_img_objects) + 1, j] for i, j in
                       enumerate(self.init.gs_img_objects, 1)]
    else:
      self.img_list = [[i, len(self.init.input_list) + 1, j] for i, j in
                       enumerate(self.init.input_list, 1)]
    cmd.Command.start("Processing {} images".format(len(self.img_list)))

    # Run processing
    self.prog_count = 0
    self.gs_prog = cmd.ProgressBar(title='PROCESSING')
    self.img_objects = parallel_map(iterable=self.img_list,
                                    func=self.proc_wrapper,
                                    callback=self.callback,
                                    processes=self.init.params.n_processors)
    cmd.Command.end("Processing {} images -- DONE "
                    "".format(len(self.img_objects)))


  def run_import(self):
    ''' Import images or image objects '''
    if self.init.params.cctbx.selection.select_only.flag_on:
      msg = "Reading {} image objects".format(len(self.init.gs_img_objects))
      title = 'READING IMAGE OBJECTS'
      self.img_list = [[i, len(self.init.gs_img_objects) + 1, j] for i, j in
                       enumerate(self.init.gs_img_objects, 1)]
    else:
      msg = "Importing {} images".format(len(self.init.input_list))
      title = 'IMPORTING IMAGES'
      self.img_list = [[i, len(self.init.input_list) + 1, j] for i, j in
                       enumerate(self.init.input_list, 1)]

    cmd.Command.start(msg)
    self.prog_count = 0
    self.gs_prog = cmd.ProgressBar(title=title)
    self.img_objects = parallel_map(iterable=self.img_list,
                                    func=self.proc_wrapper,
                                    callback=self.callback,
                                    processes=self.init.params.n_processors)

  def run_process(self):
    ''' Run indexing / integration of imported images '''
    cmd.Command.start("Processing {} images".format(len(self.img_objects)))
    self.img_list = [[i, len(self.img_objects) + 1, j] for i, j in
                      enumerate(self.img_objects, 1)]
    self.prog_count = 0
    self.gs_prog = cmd.ProgressBar(title='PROCESSING')
    self.img_objects = parallel_map(iterable=self.img_list,
                               func=self.proc_wrapper,
                               callback=self.callback,
                               processes=self.init.params.n_processors)
    cmd.Command.end("Processing {} images -- DONE "
                    "".format(len(self.img_objects)))


  def run_analysis(self):
    ''' Run analysis of integrated images '''
    cmd.Command.start("Analyzing results ")
    analysis = Analyzer(init=self.init,
                        all_objects=self.img_objects,
                        gui_mode=False)
    cmd.Command.end("Analyzing results -- DONE")
    analysis.print_results()
    analysis.unit_cell_analysis()
    analysis.print_summary()
    analysis.make_prime_input()

  def run_full_proc(self):
    ''' Run IOTA in full-processing mode (i.e. process image from import to
    integration; allows real-time tracking of output '''

    # Process images
    self.stage = 'all'
    self.run_all()

    # Analysis of integration results
    final_objects = [i for i in self.img_objects if i.fail == None]
    if len(final_objects) > 0:
      self.run_analysis()
    else:
      print('No images successfully integrated!')

    # Exit IOTA
    misc.iota_exit()


  def run(self):
    ''' Run IOTA '''

    # Import Images
    self.stage = 'import'
    self.run_import()

    # Remove rejected images from image object list
    acc_img_objects = [i.fail for i in self.img_objects if i.fail == None]
    cmd.Command.end("Accepted {} of {} images -- DONE " \
                    "".format(len(acc_img_objects), len(self.img_objects)))

    # Exit if none of the images have diffraction
    if str(self.init.params.image_triage.type).lower() != 'none':
      if len(acc_img_objects) == 0:
        misc.main_log(self.init.logfile, 'No images have diffraction!', True)
        misc.iota_exit()
      else:
        misc.main_log(self.init.logfile,
                      "{} out of {} images have diffraction "
                      "(at least {} Bragg peaks)"
                      "".format(len(acc_img_objects),
                                len(self.img_objects),
                                self.init.params.image_triage.min_Bragg_peaks))

    # Check for -c option and exit if true
    if self.init.params.image_conversion.convert_only:
      misc.iota_exit()

    # Process Images
    self.stage = 'process'
    self.run_process()

    # Analysis of integration results
    final_objects = [i for i in self.img_objects if i.fail == None]
    if len(final_objects) > 0:
      self.run_analysis()
    else:
      print('No images successfully integrated!')

    # Exit IOTA
    misc.iota_exit()


# ============================================================================ #
if __name__ == "__main__":

  iota = XTermIOTA()
  if iota.full:
    iota.run_full_proc()
  else:
    iota.run()
