from __future__ import absolute_import, division, print_function

# Template file for user preset buttons for the HKLviewer. To be placed
# in the users HOME directory.
# Additional buttons can be defined by copying current formatted PHIL
# output into new button elements in the list below. Obtain current
# formatted PHIL output from the settings dialog in the HKLviewer
# for the way in which reflections are currently shown


buttonsdeflist = [
# Simple button displaying sigmas of an intensity dataset. Uncomment lines below to take effect
#
# ("ShowSigmaI", "Show SigI", """
#                            viewer {
#                              data_array {
#                                label = "IOBS,SIGIOBS"
#                                datatype = "Intensity"
#                              }
#                            }
#                            hkls {
#                              sigma_color_radius = True
#                            }
#  """
#  ),
]
