from __future__ import division, print_function, absolute_import
# LIBTBX_SET_DISPATCHER_NAME iota.run

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 11/29/2018
Description : IOTA command-line module.
'''
import os
from contextlib import contextmanager

from iota import iota_version
from iota.components.iota_init import XInitAll
from iota.components.iota_base import ProcessingThreadBase
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

@contextmanager  # Will print start / stop messages around some processes
def prog_message(msg, msg2=None):
  cmd.Command.start(msg)
  yield
  if msg2:
    cmd.Command.end('{} -- DONE'.format(msg2))
  else:
    cmd.Command.end('{} -- DONE'.format(msg))

class XProcessAll(ProcessingThreadBase):
  """ Process module customized for command line use """
  def __init__(self, init, iterable, stage, abort_file):
    ProcessingThreadBase.__init__(self, init=init, iterable=iterable,
                                  stage=stage)
    self.prog_count = 0

    # Initialize importer and processor depending on backend
    if init.params.advanced.processing_backend == 'ha14':
      from iota.components.iota_cctbx_ha14 import ImageImporter as Importer
      from iota.components.iota_cctbx_ha14 import Integrator
    else:
      from iota.components.iota_image import ImageImporter as Importer
      from iota.components.iota_processing import Integrator
    self.importer   = Importer(init=init)
    self.integrator = Integrator(init=init)

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

  def process(self):
    """ Run Process module from command line """

    if self.stage == 'import':    # Import Images
      self.create_image_iterable()

      cmd.Command.start("Importing {} images".format(len(self.init.input_list)))
      self.prog_count = 0
      self.gs_prog = cmd.ProgressBar(title='IMPORTING IMAGES')
      self.run_process()

      # Remove rejected images from image object list
      acc_img_objects = [i.fail for i in self.img_objects if i.fail is None]
      cmd.Command.end("Accepted {} of {} images -- DONE " \
                      "".format(len(acc_img_objects), len(self.img_objects)))


      # Exit if none of the images have diffraction
      if len(acc_img_objects) == 0:
        util.main_log(self.init.logfile, 'No images have diffraction!', True)
        util.iota_exit()
      else:
        util.main_log(self.init.logfile,
                      "{} out of {} images have diffraction "
                      "(at least {} Bragg peaks)"
                      "".format(len(acc_img_objects),
                                len(self.img_objects),
                                self.init.params.image_import.minimum_Bragg_peaks))
        self.stage = 'process'

    # Process Images
    self.create_image_iterable()
    with prog_message("Processing {} images".format(len(self.iterable))):
      self.prog_count = 0
      self.gs_prog = cmd.ProgressBar(title='PROCESSING')
      self.run_process()

    # Analysis of integration results
    final_objects = [i for i in self.img_objects if i.fail is None]
    if len(final_objects) > 0:
      self.run_analysis()

      # Write info object / file (TODO: Revisit this later!)
      from iota.components.iota_threads import ObjectReaderThread
      self.object_reader = ObjectReaderThread(self,
                                              info=self.init.info,
                                              source=self.img_objects,
                                              info_file=self.init.info_file)
      self.object_reader.start()

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
      if msg:
        print (msg)
      util.iota_exit()

    if self.init.args.full:
      stage = 'all'
    else:
      stage = 'import'

    # Save init and image iterable for potential UI recovery
    from libtbx import easy_pickle
    easy_pickle.dump(self.init.init_file, self.init)
    easy_pickle.dump(self.init.iter_file, self.init.input_list)

    abort_file = os.path.join(self.init.int_base, '.abort')
    processor = XProcessAll(init=self.init, iterable=self.init.input_list,
                            stage=stage, abort_file=abort_file)
    processor.start()

# ============================================================================ #
if __name__ == "__main__":

  iota = XTermIOTA()
  iota.run()
