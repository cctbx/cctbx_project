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

  ("Amplitudes", "Show Amplitudes", """
            viewer {
              data_array {
                label = "FOBS,SIGFOBS"
                datatype = "Amplitude"
              }
            }
"""),
# ConstantAxesSliceIntens and ConstantAxesSliceAmpl rely on hard coded names "H-axis (1,0,0)", "K-axis (0,1,0)"
# and "L-axis (0,0,1)" being present in the list of vectors. The button PHIL parameter show_vector="['-axis', True]"
# will then entail comboboxes being created from where H, K or L axes can be selected. This is more compact
# than having threee separate buttons for each axes
  ("ConstantAxesSliceIntens", "Show plane of intensities with constant: ", """
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

   ("ConstantAxesSliceAmpl", "Show plane of amplitudes with constant: ", """
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

("FoversigF", "F/SigF",
 """
          miller_array_operation = "('newarray._data= array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'FoverSigF', ['FOBS,SIGFOBS', 'Amplitude'], ['', ''])"
          viewer {
            data_array {
              label = "FoverSigF"
              datatype = "Amplitude"
            }
          }

 """),

("IoverSigI", "I/SigI",
 """
        miller_array_operation = "('newarray._data=array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'IoverSigI', ['I<<FSQ,SIGI<<FSQ', 'Intensity'], ['', ''])"
        binning {
          scene_bin_thresholds = -10000 1 2 3 4 5 460 793.55 2750
          binlabel = "IoverSigI"
          bin_opacity = 1 0
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

("Evalues", "E-values",
 """
          miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = None', 'E-values', ['FP', 'Amplitude'], ['', ''])"
          viewer {
            data_array {
              label = "E-values"
              datatype = "Amplitude"
            }
          }
 """),

("Merged", "Merged Intensities",
 """
      miller_array_operation = "('from crys3d.hklviewer import display2\\nnewarray = display2.MergeData( array1, show_anomalous_pairs=False)[0]\\nfrom cctbx.xray import observation_types\\nnewarray.set_observation_type( observation_types.intensity())', 'Imerge', ['I,SIGI','Intensity'], ['', ''])"
      viewer {
        data_array {
          label = "Imerge,SigImerge"
          datatype = "Intensity"
        }
      }

 """),

("Multiplicities", "Multiplicities",
 """
      miller_array_operation = "('from crys3d.hklviewer import display2\\nmultiplicities = display2.MergeData( array1, show_anomalous_pairs=False)[1]\\n# use double to avoid being interpreted as R-free\\nnewarray._data = multiplicities.data().as_double()\\nnewarray._indices = multiplicities.indices()\\nnewarray._sigmas = None\\nfrom cctbx.xray import observation_types\\nnewarray.set_observation_type(None)', 'multiplicity', ['I,SIGI','Intensity'], ['', ''])"
      viewer {
        data_array {
          label = "multiplicity"
          datatype = "Floating-point"
        }
      }

 """),

("INFO035", "INFO < 0.35 bits", """
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

("InfoISigI05", "INFO > 0.2 bits and I/SigI < 0.5", """
        miller_array_operation = "('ISigIarray = array2.deep_copy()\\nISigIarray._data = array2.data()/array2.sigmas()\\nnewarray = array1.select(ISigIarray.data()<0.5)', 'INFO_ISigI_0.5', ['INFO', None], ['IOBS,SIGIOBS', 'Intensity'])"
        binning {
          scene_bin_thresholds = -0.1 0.2 nan nan
          binlabel = 'INFO_ISigI_0.5'
          bin_opacity = 0 0
          bin_opacity = 1 1
          bin_opacity = 1 2
          bin_opacity = 1 3
          nbins = 4
        }
        viewer {
          data_array {
            label = "INFO_ISigI_0.5"
            datatype = "Floating-point"
          }
        }
 """),
 # we omit datatype of INFO in miller_array_operation as to force validate_preset_buttons()
 # only to match an array that is labelled INFO
  ("aniso4", "Rotate around one anisotropy principal axis:", """
        real_space_unit_cell_scale_fraction = 0.9
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
          draw_real_space_unit_cell = True
          real_space_unit_cell_scale_fraction = 0.9
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
  ("TNCSlayer_xtricorder", "Slice perpendicular to tNCS_xtricorder vector", """
          clip_plane {
            normal_vector = "tNCS_xtricorder"
            clip_width = 0.380397231
          }
          viewer {
            data_array {
              label = "TEPS"
              datatype = "Floating-point"
            }
            show_vector = "['tNCS_xtricorder', True]"
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

"""),
  ("TNCSvecrotate_xtricorder", "Rotate around tNCS_xtricorder vector", """
            clip_plane {
              clip_width = 6
              auto_clip_width = False
            }
            viewer {
              data_array {
                label = "TEPS"
                datatype = "Floating-point"
              }
              show_vector = "['tNCS_xtricorder', True]"
              is_parallel = True
              fixorientation = *vector None
              animate_rotation_around_vector = "['tNCS_xtricorder', 5.0]"
            }
            hkls {
              expand_to_p1 = True
              expand_anomalous = True
            }
"""),
  ("TNCSlayer_xtriage", "Slice perpendicular to tNCS_xtriage vector", """
          clip_plane {
            normal_vector = "tNCS_xtriage"
            clip_width = 0.380397231
          }
          miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = None', 'E-values', ['FP', 'Amplitude'], ['', ''])"
          viewer {
            data_array {
              label = "E-values"
              datatype = "Amplitude"
            }
            show_vector = "['tNCS_xtriage', True]"
            user_vector {
              label = "tNCS_xtriage"
            }
            fixorientation = *vector None
          }
          hkls {
            expand_to_p1 = True
            expand_anomalous = True
          }

"""),
  ("TNCSvecrotate_xtriage_F", "Rotate around tNCS_xtriage vector", """
            clip_plane {
              clip_width = 6
              auto_clip_width = False
            }
            miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = None', 'E-values(F)', ['FP', 'Amplitude'], ['', ''])"
            viewer {
              data_array {
                label = "E-values(F)"
                datatype = "Amplitude"
              }
              show_vector = "['tNCS_xtriage', True]"
              user_vector {
                label = "tNCS_xtriage"
              }
              is_parallel = True
              fixorientation = *vector None
              animate_rotation_around_vector = "['tNCS_xtriage', 5.0]"
            }
            hkls {
              expand_to_p1 = True
              expand_anomalous = True
            }
"""),
  ("TNCSvecrotate_xtriage_I", "Rotate around tNCS_xtriage vector", """
            clip_plane {
              clip_width = 6
              auto_clip_width = False
            }
            miller_array_operation = "('newarray._data = array1.normalize().data()\\nnewarray._sigmas = None', 'E-values(I)', ['I', 'Intensity'], ['', ''])"
            viewer {
              data_array {
                label = "E-values(I)"
                datatype = "Amplitude"
              }
              show_vector = "['tNCS_xtriage', True]"
              user_vector {
                label = "tNCS_xtriage"
              }
              is_parallel = True
              fixorientation = *vector None
              animate_rotation_around_vector = "['tNCS_xtriage', 5.0]"
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

]
