from __future__ import absolute_import, division, print_function

buttonsdeflist = [
  ("H_I", "Show Zone H for intensities", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = 0
                                                  is_assoc_real_space_vector = True
                                                  clipwidth = 0.5
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                  color_scheme = *CMRmap
                                                  color_powscale = 0.2633312543
                                                }
  """),
  ("K_I", "Show Zone K for intensities", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = 1
                                                  is_assoc_real_space_vector = True
                                                  clipwidth = 0.650000
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                  color_scheme = *CMRmap
                                                  color_powscale = 0.2633312543
                                                }
  """),
  ("L_I", "Show Zone L for intensities", """
                                                clip_plane {
                                                  hkldist = 0
                                                  normal_vector = 2
                                                  is_assoc_real_space_vector = True
                                                  clipwidth = 0.5
                                                }
                                                viewer {
                                                  data_array.label = "I,SIGI"
                                                  data_array.datatype = "Intensity"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = nan
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                  color_scheme = *CMRmap
                                                  color_powscale = 0.2633312543
                                                }
  """),
  ("aniso", "Show Anisotropy", ""),
  ("EINFO", "Show Information in bits", ""),
  ("TNCS", "Show TNCS normal", """
                                                clip_plane {
                                                  angle_around_vector = "[3, 0]"
                                                  animate_rotation_around_vector = "[0, -1.000000]"
                                                  normal_vector = 3
                                                  normal_vector_length_scale = 0.87128
                                                  is_assoc_real_space_vector = False
                                                  clipwidth = 0.35
                                                }
                                                viewer {
                                                  data_array.label = "TEPS"
                                                  data_array.datatype = "Floating-point"
                                                  show_vector = "[3, False]"
                                                  fixorientation = *vector None
                                                  nth_power_scale_radii = 0.1
                                                  expand_to_p1 = True
                                                  expand_anomalous = True
                                                }
  """),
  ("TNCSpar", "Show TNCS paralllel", """
                                              clip_plane {
                                                angle_around_vector = "[3, 0.0]"
                                                animate_rotation_around_vector = "[3, 6.000000]"
                                                is_assoc_real_space_vector = False
                                                clipwidth = 3
                                              }
                                              viewer {
                                                data_array.label = "TEPS"
                                                data_array.datatype = "Floating-point"
                                                show_vector = "[3, True]"
                                                is_parallel = True
                                                fixorientation = *vector None
                                                nth_power_scale_radii = 0.1
                                                expand_to_p1 = True
                                                expand_anomalous = True
                                                color_scheme = *rainbow
                                              }
  """)


]
