from __future__ import absolute_import, division, print_function

cctbx_buttonsdeflist = [
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
# ConstanAxesSliceIntens and ConstanAxesSliceAmpl rely on hard coded names "H-axis (1,0,0)", "K-axis (0,1,0)"
# and "L-axis (0,0,1)" being present in the list of vectors. The button PHIL parameter show_vector="['-axis', True]"
# will then entail comboboxes being created from where H, K or L axes can be selected. This is more compact
# than having threee separate buttons for each axes
  ("ConstanAxesSliceIntens", "Show plane of intensities with constant: ", """
          clip_plane {
            normal_vector = "-axis"
            is_assoc_real_space_vector = True
            clip_width = 1.184
          }
          viewer {
            data_array {
              label = 'I,SIGI'
              datatype = 'Intensity'
            }
            show_vector = "['-axis', True]"
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }
 """),

   ("ConstanAxesSliceAmpl", "Show plane of amplitudes with constant: ", """
          clip_plane {
            normal_vector = "-axis"
            is_assoc_real_space_vector = True
            clip_width = 1.184
          }
          viewer {
            data_array {
              label = 'FP,SIGFP'
              datatype = 'Amplitude'
            }
            show_vector = "['-axis', True]"
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
          bin_opacity = 0 1
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

]





phenix_buttonsdeflist = [
("INFO", "Reflections with information bits less than 0.35", """
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

  ("aniso4", "Rotate around one anisotropy principal axis:", """
        binning {
          binlabel = "ANISO"
          bin_opacity = 1 0
          bin_opacity = 1 1
          bin_opacity = 0 2
          bin_opacity = 0 3
          bin_opacity = 0 4
          bin_opacity = 0 5
          bin_opacity = 1 6
          bin_opacity = 1 7
          nbins = 8
        }
        viewer {
          animate_rotation_around_vector = "['Anisotropy', 5.0]"
          data_array {
            label = "ANISO"
            datatype = "Floating-point"
          }
          show_vector = "['Anisotropy', True]"
          is_parallel = True
          fixorientation = *vector None
        }
        hkls {
          expand_to_p1 = True
          expand_anomalous = True
        }
  """),

  ("aniso", "Show anisotropy principal axes", """
          real_space_unit_cell_scale_fraction = 0
          binning {
            binlabel = 'ANISO'
            bin_opacity = 1 0
            bin_opacity = 1 1
            bin_opacity = 0 2
            bin_opacity = 0 3
            bin_opacity = 0 4
            bin_opacity = 0 5
            bin_opacity = 1 6
            bin_opacity = 1 7
            nbins = 8
          }
          viewer {
            data_array {
              label = "ANISO"
              datatype = "Floating-point"
            }
            show_vector = "['Anisotropy1', True]"
            show_vector = "['Anisotropy2', True]"
            show_vector = "['Anisotropy3', True]"
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }
  """),
  ("TNCSlayer", "Slice perpendicular to TNCS vector", """
          clip_plane {
            normal_vector = "TNCS"
            clip_width = 0.380397231
          }
          viewer {
            data_array {
              label = "TEPS"
              datatype = "Floating-point"
            }
            show_vector = "['TNCS', True]"
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

"""),
  ("TNCSvecrotate", "Rotate around TNCS vector", """
            clip_plane {
              clip_width = 6
              auto_clip_width = False
            }
            viewer {
              data_array {
                label = "TEPS"
                datatype = "Floating-point"
              }
              show_vector = "['TNCS', True]"
              is_parallel = True
              fixorientation = *vector None
              animate_rotation_around_vector = "['TNCS', 5.0]"
            }
            hkls {
              expand_to_p1 = True
              expand_anomalous = True
            }
"""),

  # If user defines a twin operator where the name of the operator contains the string "twin" this
  # button will become enabled.
  # Pressing this button will expand amplitude data to P1 and slice it with twin axis perpendicular
  # to the screen. One can then step through layers of reflections with the +/- buttons in the GUI
  ("TwinAxisampl", "Slice amplitudes perpendicular to twin axis", """

          clip_plane {
            hkldist = 0.0
            normal_vector = "twin"
            clip_width = 0.5
          }
          viewer {
            data_array {
              label = "F,SIGF"
              datatype = "Amplitude"
            }
            show_vector = "['twin', True]"
            user_vector {
              label = "twin"
            }
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

  """),

  ("TwinAxisrotampl", "Rotate amplitudes around twin axis", """

          clip_plane {
            clip_width = 10
            auto_clip_width = False
          }
          viewer {
            data_array {
              label = "F,SIGF"
              datatype = "Amplitude"
            }
            show_vector = "['twin', True]"
            is_parallel = True
            user_vector {
              label = "twin"
            }
            fixorientation = *vector None
            animate_rotation_around_vector = "['twin', 2.0]"
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

  """),

  ("TwinAxisintens", "Slice intensities perpendicular to twin axis", """

          clip_plane {
            hkldist = 0.0
            normal_vector = "twin"
            clip_width = 0.5
          }
          viewer {
            data_array {
              label = "I,SIGI"
              datatype = "Intensity"
            }
            show_vector = "['twin', True]"
            user_vector {
              label = "twin"
            }
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

  """),

  ("TwinAxisintensrot", "Rotate intensities around twin axis", """

          clip_plane {
            clip_width = 10
            auto_clip_width = False
          }
          viewer {
            data_array {
              label = "I,SIGI"
              datatype = "Intensity"
            }
            show_vector = "['twin', True]"
            is_parallel = True
            user_vector {
              label = "twin"
            }
            fixorientation = *vector None
            animate_rotation_around_vector = "['twin', 2.0]"
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

  """),

]
