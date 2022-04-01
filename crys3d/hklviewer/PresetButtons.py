from __future__ import absolute_import, division, print_function

buttonsdeflist = [
  ("H_I", "Show Zone H for intensities", """
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
"""),
  ("K_I", "Show Zone K for intensities", """
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
  """),
  ("L_I", "Show Zone L for intensities", """
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
"""),
  ("H_F", "Show Zone H for amplitudes", """
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

  """),

  ("K_F", "Show Zone K for amplitudes", """
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

  """),

  ("L_F", "Show Zone L for amplitudes", """
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

  """),


  ("Test", "I,SIGI test", """
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                }
  """),

]
