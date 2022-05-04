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

("FoversigF", "F/SigF",
 """
                miller_array_operation = "('newarray._data= array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'FoverSigF2', ['FOBS,SIGFOBS', 'Amplitude'], ['', ''])"
                viewer {
                  data_array {
                    label = "FoverSigF2"
                    datatype = "Amplitude"
                  }
                }


 """),
("IoverSigI", "I/SigI",
 """
              miller_array_operation = "('newarray._data = array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'IoverSigI', ['I<<FSQ,SIGI<<FSQ', 'Intensity'], ['', ''])"
              viewer {
                data_array {
                  label = "IoverSigI"
                  datatype = "Intensity"
                }
              }
 """),
("Evalues", "E-values",
 """
                    miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = array1.normalize().sigmas()\\n', 'E-values', ['FP,SIGFP', 'Amplitude'], ['', ''])"
                    viewer {
                      data_array {
                        label = "E-values"
                        datatype = "Amplitude"
                      }
                    }
 """),




]
