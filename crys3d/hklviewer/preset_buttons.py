from __future__ import absolute_import, division, print_function

buttonsdeflist = [
  ("H_I", "Show plane of intensities with constant H", """
                                                          clip_plane {
                                                            normal_vector = "H (1,0,0)"
                                                            is_assoc_real_space_vector = True
                                                            clip_width = 1.184439576
                                                          }
                                                          viewer {
                                                            data_array {
                                                              label = 'I<<FSQ,SIGI<<FSQ'
                                                              datatype = 'Intensity'
                                                            }
                                                            fixorientation = *vector None
                                                          }
                                                          hkls {
                                                            expand_to_p1 = True
                                                            expand_anomalous = True
                                                          }
"""),
  ("K_I", "Show plane of intensities with constant K", """

                                                            clip_plane {
                                                              normal_vector = "K (0,1,0)"
                                                              is_assoc_real_space_vector = True
                                                              clip_width = 1.184
                                                            }
                                                            viewer {
                                                              data_array {
                                                                label = 'I<<FSQ,SIGI<<FSQ'
                                                                datatype = 'Intensity'
                                                              }
                                                              fixorientation = *vector None
                                                            }
                                                            hkls {
                                                              expand_to_p1 = True
                                                              expand_anomalous = True
                                                            }
"""),
  ("L_I", "Show plane of intensities with constant L", """
                                                            clip_plane {
                                                              normal_vector = "L (0,0,1)"
                                                              is_assoc_real_space_vector = True
                                                              clip_width = 1.184
                                                            }
                                                            viewer {
                                                              data_array {
                                                                label = 'I<<FSQ,SIGI<<FSQ'
                                                                datatype = 'Intensity'
                                                              }
                                                              fixorientation = *vector None
                                                            }
                                                            hkls {
                                                              expand_to_p1 = True
                                                              expand_anomalous = True
                                                            }
"""),
  ("H_F", "Show plane of amplitudes with constant H", """
                                                            clip_plane {
                                                              normal_vector = "H (1,0,0)"
                                                              is_assoc_real_space_vector = True
                                                              clip_width = 1.184
                                                            }
                                                            viewer {
                                                              data_array {
                                                                label = 'FP,SIGFP'
                                                                datatype = 'Amplitude'
                                                              }
                                                              fixorientation = *vector None
                                                            }
                                                            hkls {
                                                              expand_to_p1 = True
                                                              expand_anomalous = True
                                                            }
  """),

  ("K_F", "Show plane of amplitudes with constant K", """
                                                            clip_plane {
                                                              normal_vector = "K (0,1,0)"
                                                              is_assoc_real_space_vector = True
                                                              clip_width = 1.184
                                                            }
                                                            viewer {
                                                              data_array {
                                                                label = 'FP,SIGFP'
                                                                datatype = 'Amplitude'
                                                              }
                                                              fixorientation = *vector None
                                                            }
                                                            hkls {
                                                              expand_to_p1 = True
                                                              expand_anomalous = True
                                                            }

  """),

  ("L_F", "Show plane of amplitudes with constant L", """
                                                            clip_plane {
                                                              normal_vector = "L (0,0,1)"
                                                              is_assoc_real_space_vector = True
                                                              clip_width = 1.184
                                                            }
                                                            viewer {
                                                              data_array {
                                                                label = 'FP,SIGFP'
                                                                datatype = 'Amplitude'
                                                              }
                                                              fixorientation = *vector None
                                                            }
                                                            hkls {
                                                              expand_to_p1 = True
                                                              expand_anomalous = True
                                                            }
  """),




]
