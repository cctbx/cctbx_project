from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.single_image

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/13/2016
Description : Single-image IOTA module. Version 1.0.001G
'''

iota_version = '1.0.001G'
help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

Single image mode
Usage: iota.single_image [OPTIONS] path/to/image.file
Imports, converts, triages, indexes and integrates a single
diffraction image using CCTBX defaults. Can use either IOTA
default parameters, or individual parameters can be changed
in the command line.
"""

from iota.components.iota_init import InitAll
import iota.components.iota_image as img
import iota.components.iota_cmd as cmd
import iota.components.iota_misc as misc

def run_one_image(image, init):

  # Import image
  single_image = img.SingleImage(image, init, verbose=True)
  img_object = single_image.import_image()

  # Exit if the image does not have diffraction
  if img_object.fail == 'failed triage':
    return img_object

  # Process image
  img_object.process(single_image=True)

  return img_object

# ============================================================================ #
if __name__ == "__main__":

  # Initialize IOTA parameters and log
  init = InitAll(iota_version, help_message)
  init.run()

  # Run single image
  image = [1, 1, init.input_list[0]]

  cmd.Command.start("Processing single image")
  img_object = run_one_image(image, init)
  cmd.Command.end("Processing single image -- DONE")

  misc.iota_exit(iota_version)
