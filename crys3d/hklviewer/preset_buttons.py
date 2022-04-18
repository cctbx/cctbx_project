from __future__ import absolute_import, division, print_function

buttonsdeflist = [
  ("H_I", "Show plane of intensities with constant H", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = "H (1,0,0)"
                                                  is_assoc_real_space_vector = True
                                                  clip_width = 0.5
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  #nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                }
                                                NGL.fontsize = 10
"""),
  ("K_I", "Show plane of intensities with constant K", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = "K (0,1,0)"
                                                  is_assoc_real_space_vector = True
                                                  clip_width = 0.650000
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                }
                                                NGL.fontsize = 10
  """),
  ("L_I", "Show plane of intensities with constant L", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = "L (0,0,1)"
                                                  is_assoc_real_space_vector = True
                                                  clip_width = 0.5
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                }
                                                NGL.fontsize = 10
"""),
  ("H_F", "Show plane of amplitudes with constant H", """
                                              clip_plane {
                                                hkldist = 0
                                                normal_vector = "H (1,0,0)"
                                                is_assoc_real_space_vector = True
                                                clip_width = 0.50000
                                              }
                                              viewer {
                                                data_array.label = "F,SIGF"
                                                data_array.datatype = "Amplitude"
                                                fixorientation = *vector None
                                                nth_power_scale_radii = nan
                                                expand_to_p1 = True
                                                expand_anomalous = True
                                              }
                                               NGL.fontsize = 10

  """),

  ("K_F", "Show plane of amplitudes with constant K", """
                                              clip_plane {
                                                hkldist = 0
                                                normal_vector = "K (0,1,0)"
                                                is_assoc_real_space_vector = True
                                                clip_width = 0.50000
                                              }
                                              viewer {
                                                data_array.label = "F,SIGF"
                                                data_array.datatype = "Amplitude"
                                                fixorientation = *vector None
                                                nth_power_scale_radii = nan
                                                expand_to_p1 = True
                                                expand_anomalous = True
                                              }
                                               NGL.fontsize = 10

  """),

  ("L_F", "Show plane of amplitudes with constant L", """
                                              clip_plane {
                                                hkldist = 0
                                                normal_vector = "L (0,0,1)"
                                                is_assoc_real_space_vector = True
                                                clip_width = 0.50000
                                              }
                                              viewer {
                                                data_array.label = "F,SIGF"
                                                data_array.datatype = "Amplitude"
                                                fixorientation = *vector None
                                                nth_power_scale_radii = nan
                                                expand_to_p1 = True
                                                expand_anomalous = True
                                              }
                                               NGL.fontsize = 10
  """),


]
