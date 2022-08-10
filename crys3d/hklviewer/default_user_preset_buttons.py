from __future__ import absolute_import, division, print_function

# Template file for user preset buttons for the HKLviewer. To be placed
# in the users HOME directory.
# Additional buttons can be defined by copying current formatted PHIL
# output into new button elements in the list below. Obtain current
# formatted PHIL output from the settings dialog in the HKLviewer
# for the way in which reflections are currently shown


buttonsdeflist = [
  # If user defines a twin operator where the name of the operator contains the string "twin" this
  # button will become enabled.
  # Pressing this button will expand amplitude data to P1 and slice it with twin axis perpendicular
  # to the screen. One can then step through layers of reflections with the +/- buttons in the GUI
  ("TwinAxis", "Slice perpendicular to twin axis", """
                                        clip_plane {
                                          hkldist = 0
                                          normal_vector = "twin"
                                          clip_width = 0.5
                                        }
                                        viewer {
                                          data_array {
                                            label = "F-obs,SIGF-obs"
                                            datatype = "Amplitude"
                                          }
                                          show_vector = "['twin', True]"
                                          user_vector {
                                            label = "twin"
                                          }
                                          fixorientation = *vector None
                                          angle_around_vector = "['twin', 0.0]"
                                        }
                                        hkls {
                                          expand_to_p1 = True
                                          expand_anomalous = True
                                        }

  """),

]
