from __future__ import absolute_import, division, print_function

# Template file for user preset buttons for the HKLviewer. To be placed in the users HOME directory.
# Additional buttons can be defined by copying current formatted PHIL output into new button elements in the
# list below. Obtain current formatted PHIL output from the settings dialog in the HKLviewer for your
# favouruite way of displaying certain reflections


buttonsdeflist = [
  # If user defines a twin operator where the name of the operator contains the string "twin" this
  # button will become enabled.
  # Pressing this button will expand amplitude data to P1 and slice it with twin axis perpendicular
  # to the screen. One can then step through layers of reflections with the +/- buttons in the GUI
  ("TwinAxis", "Slice perpendicular to twin axis", """

                                        clip_plane {
                                          angle_around_vector = "['twin', 0.0]"
                                          animate_rotation_around_vector = "['twin', -1.0]"
                                          hkldist = 0
                                          normal_vector = "twin"
                                          clip_width = 1
                                          auto_clip_width = False
                                        }
                                        viewer {
                                          data_array.label = "F,SIGF"
                                          data_array.datatype = "Amplitude"
                                          show_vector = "['twin', True]"
                                          user_label = "twin"
                                          fixorientation = *vector None
                                          nth_power_scale_radii = nan
                                          expand_to_p1 = True
                                          expand_anomalous = True
                                          color_scheme = *jet
                                          color_powscale = 0.5131581182
                                        }
                                        NGL {
                                          show_tooltips = *click
                                          fontsize = 10
                                        }
  """)


]
