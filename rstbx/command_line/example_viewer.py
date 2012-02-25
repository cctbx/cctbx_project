# LIBTBX_SET_DISPATCHER_NAME phenix.example_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import os
import sys
from rstbx.command_line.viewer import run

def modified_image_factory (filename):
  from iotbx.detectors.generic_detector import GenericDetector
  I = GenericDetector(filename)
  I.readHeader()
  return I

def modify_the_iotbx_detector_list():
  from iotbx import detectors
  detectors.ImageFactory = modified_image_factory

if (__name__ == "__main__") :
  modify_the_iotbx_detector_list()
  file_arguments = sys.argv[1:]
  if len(file_arguments) > 0:
    run(sys.argv[1:])
  else:
    print "Use phenix.example_viewer <data file name[s]>"
