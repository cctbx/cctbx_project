from __future__ import absolute_import, division, print_function

buttonsdeflist = [
  ("Intensities", "Show Intensities", """
                                            viewer {
                                              data_array {
                                                label = "I,SIGI"
                                                datatype = "Intensity"
                                              }
                                            }
"""),
  ("amplitudes", "Show Amplitudes", """
                                            viewer {
                                              data_array {
                                                label = "FOBS,SIGFOBS"
                                                datatype = "Amplitude"
                                              }
                                            }
"""),
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

("FoversigF", "F/SigF ( from miller_array.data()/miller_array.sigmas() )",
 """
                miller_array_operation = "('newarray._data= array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'FoverSigF2', ['FOBS,SIGFOBS', 'Amplitude'], ['', ''])"
                viewer {
                  data_array {
                    label = "FoverSigF2"
                    datatype = "Amplitude"
                  }
                }


 """),
("IoverSigI", "I/SigI >= 2 ( from miller_array.data()/miller_array.sigmas() )",
 """
              miller_array_operation = "('newarray._data=array1.data()/array1.sigmas()', 'IoverSigI', ['I<<FSQ,SIGI<<FSQ', 'Intensity'], ['', ''])"
              binning {
                scene_bin_thresholds = -10000 1 2 3 4 5 460 793.55 2750
                binlabel = "IoverSigI"
                bin_opacity = 0 0
                bin_opacity = 1 1
                bin_opacity = 1 2
                bin_opacity = 1 3
                bin_opacity = 1 4
                bin_opacity = 1 5
                bin_opacity = 1 6
                bin_opacity = 1 7
                bin_opacity = 1 8
                bin_opacity = 1 9
                bin_opacity = 1 10
                nbins = 9
              }
              viewer {
                data_array {
                  label = "IoverSigI"
                  datatype = "Intensity"
                }
              }
 """),
("Evalues", "E-values ( miller_array.normalize() )",
 """
                    miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = array1.normalize().sigmas()\\n', 'E-values', ['FP,SIGFP', 'Amplitude'], ['', ''])"
                    viewer {
                      data_array {
                        label = "E-values"
                        datatype = "Amplitude"
                      }
                    }
 """),
("EINFO", "eInfo < 0.35",
 """
                    binning {
                      scene_bin_thresholds = -1 0.35 0.7 1 1.25 1.85 17.59 100
                      binlabel = 'EINFO'
                      bin_opacity = 1 0
                      bin_opacity = 0 1
                      bin_opacity = 0 2
                      bin_opacity = 0 3
                      bin_opacity = 0 4
                      bin_opacity = 0 5
                      bin_opacity = 0 6
                      bin_opacity = 0 7
                      nbins = 8
                    }
                    viewer {
                      data_array {
                        label = "EINFO"
                      }
                    }
 """),
("INFO", "Info < 0.35",
 """
                    binning {
                      scene_bin_thresholds = -1 0.35 0.7 1 1.25 1.85 17.59 100
                      binlabel = 'INFO'
                      bin_opacity = 1 0
                      bin_opacity = 0 1
                      bin_opacity = 0 2
                      bin_opacity = 0 3
                      bin_opacity = 0 4
                      bin_opacity = 0 5
                      bin_opacity = 0 6
                      bin_opacity = 0 7
                      nbins = 8
                    }
                    viewer {
                      data_array {
                        label = "INFO"
                      }
                    }
 """),




]
