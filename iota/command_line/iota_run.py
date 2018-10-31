from __future__ import division, print_function, absolute_import
# LIBTBX_SET_DISPATCHER_NAME iota.run

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 10/18/2018
Description : IOTA command-line module.
'''
import os

from iota import iota_version
from iota.components.iota_init import XInitAll
from iota.components.iota_base import ProcessGeneral
import dials.util.command_line as cmd
import iota.components.iota_utils as util

iota_version = iota_version

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

class XProcessAll(ProcessGeneral):
  """ Process module customized for command line use """
  def __init__(self, init, iterable, stage, abort_file):
    ProcessGeneral.__init__(self, init=init, iterable=iterable,
                            stage=stage, abort_file=abort_file)
    self.prog_count = 0

  def callback(self, result):
    """ To be run on completion of each step
    @param result: image_object
    """
    if self.prog_count < len(self.init.input_list):
      prog_step = 100 / len(self.init.input_list)
      self.gs_prog.update(self.prog_count * prog_step)
      self.prog_count += 1
    else:
      self.gs_prog.finished()

  def run(self):
    """ Run Process module from command line """

    if self.stage == 'import':    # Import Images
      self.create_image_iterable()

      if self.init.params.cctbx_ha14.selection.select_only.flag_on:
        msg = "Reading {} image objects".format(len(self.init.gs_img_objects))
        title = 'READING IMAGE OBJECTS'
      else:
        msg = "Importing {} images".format(len(self.init.input_list))
        title = 'IMPORTING IMAGES'
      cmd.Command.start(msg)
      self.prog_count = 0
      self.gs_prog = cmd.ProgressBar(title=title)
      self.run_process()

      # Remove rejected images from image object list
      acc_img_objects = [i.fail for i in self.img_objects if i.fail is None]
      cmd.Command.end("Accepted {} of {} images -- DONE " \
                      "".format(len(acc_img_objects), len(self.img_objects)))


      # Exit if none of the images have diffraction
      if str(self.init.params.cctbx_ha14.image_triage.type).lower() != 'none':
        if len(acc_img_objects) == 0:
          util.main_log(self.init.logfile, 'No images have diffraction!', True)
          util.iota_exit()
        else:
          util.main_log(self.init.logfile,
                        "{} out of {} images have diffraction "
                        "(at least {} Bragg peaks)"
                        "".format(len(acc_img_objects),
                                  len(self.img_objects),
                                  self.init.params.cctbx_ha14.image_triage.min_Bragg_peaks))
          self.stage = 'process'

      # Check for -c option and exit if true
      if self.init.params.cctbx_ha14.image_conversion.convert_only:
        util.iota_exit()

    # Process Images
    self.create_image_iterable()
    cmd.Command.start("Processing {} images".format(len(self.img_list)))
    self.prog_count = 0
    self.gs_prog = cmd.ProgressBar(title='PROCESSING')
    self.run_process()
    cmd.Command.end("Processing {} images -- DONE".format(len(self.img_list)))

    # Analysis of integration results
    final_objects = [i for i in self.img_objects if i.fail is None]
    if len(final_objects) > 0:
      self.run_analysis()
    else:
      print ('No images successfully integrated!')

    # Exit IOTA
    util.iota_exit()


class XTermIOTA():
  """ Main class that will initalize and run everything"""
  def __init__(self):
    self.init = XInitAll(help_message)

  def run(self):
    good_init, msg = self.init.run()  # Returns False if something goes wrong

    if not good_init:
      print (msg)
      util.iota_exit()

    if self.init.args.full:
      stage = 'all'
    else:
      stage = 'import'

    abort_file = os.path.join(self.init.int_base, '.abort')
    processor = XProcessAll(init=self.init, iterable=self.init.input_list,
                            stage=stage, abort_file=abort_file)
    processor.start()

# ============================================================================ #
if __name__ == "__main__":

  iota = XTermIOTA()
  iota.run()
