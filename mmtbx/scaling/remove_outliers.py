from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import sgtbx
from cctbx import xray
from cctbx.array_family import flex
from libtbx.utils import Sorry, multi_out
import iotbx.phil
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import pdb
import mmtbx.scaling
from mmtbx.scaling import outlier_rejection
from mmtbx import f_model
from libtbx.str_utils import StringIO
import sys, os



master_params = iotbx.phil.parse("""\
outlier_utils{
  input{
    unit_cell=None
    .type=unit_cell
    space_group=None
    .type=space_group
    xray_data{
      file_name=None
      .type=path
      obs_labels=None
      .type=str
      free_flags=None
      .type=str
    }
    model{
      file_name=None
      .type=path
    }
  }
  outlier_detection{
    protocol=basic *extreme beamstop model
    .type=choice
    parameters{
      basic_wilson{
        level=1E-6
        .type=float
      }
      extreme_wilson{
        level=0.01
        .type=float
      }
      beamstop{
        level=0.001
        .type=float
        d_min=10.0
        .type=float
      }
      model_based{
        level=0.01
        .type=float
      }
    }
  }
  additional_parameters{
    free_flag_generation{
      fraction = 0.10
      .type=float
      max_number = 2000
      .type=float
      lattice_symmetry_max_delta=5.0
      .type=float
      use_lattice_symmetry=true
      .type=bool
    }
  }

  output{
    logfile=outlier_tools.log
    .type=str
    hklout=None
    .type=path
  }
}
""")

def print_help(command_name):
  print("""
usage: %(command_name)s <options>

Options are defined by the following phil scope:

outlier_utils{
  input{
    unit_cell=None                   << Unit cell. Needed for CNS and SHELX formatted Fobs/Iobs files
    space_group=None                 << space group.  Needed for CNS and SHELX formatted Fobs/Iobs files
    xray_data{
      file_name=None                 << reflection file in any format.
      obs_labels=None                << The labels of the the observed labels (CNS and MTZ)
      free_flags=None                << The labels of the Free Flags (CNS and MTZ)
    }
    model{
      file_name=None                 << A PDB file of the model corresponding to the provided xray data
    }
  }
  outlier_detection{
    protocol=basic *extreme beamstop model   << outlier protocol. See below.
    parameters{
      basic_wilson{
        level=1E-6                   << Outlier rejection level for protocol basic
      }
      extreme_wilson{
        level=0.01                   << Outlier rejection level for protocol extreme
      }
      model_based{
        level=0.01                   << outlier rejection level for protocol model
      }
    }
  }
  additional_parameters{
    free_flag_generation{            << If no free flags are provided, they will be generated
      fraction = 0.10                << Fraction of free flags
      max_number = 2000              << Maximum number of free flags
      lattice_symmetry_max_delta=5.0 << Lattice symmetry maximum delta (no need to touch this)
      use_lattice_symmetry=true      << whether or not to use the lattice symmetry in free flag generation
    }
  }

  output{
    logfile=outlier_tools.log        << The output logfile. Has a copy of the screen dump and (possibly) added plots
    hklout=None                      << an mtz file with observations that are not outliers according to the selected protocol
  }
}


A basic run looks like:

%(command_name)s data.file=my_precious_data.mtz \\
  data.obs=FOBS data.free=TEST \\
  model.file=my_baby.pdb proto=model hklout=new.mtz

A short description of the 3 outlier protocols are described below

1. protocol: basic
   Outliers are reflection for which 1 - P(E^2=>E_obs^2) < level
   For the default value, this results in (approximately) designating
   reflections with E value larger then 3.8 as an outlier.
    see Read, Acta Cryst. (1999). D55, 1759-1764 for details.

2. protocol: extreme
   Outliers are reflections for which 1 - [P(E^2=>E_obs^2)]^Nobs < level
   Nobs is the number of observations. In this manner, the size of the dataset
   is taken into account in the decision if something is an outlier.
   The basic and extrenme protocol usually result in similar outliers.

3. protocol: model
   Calculated amplitudes and estimated values of alpha and beta
   are used to compute the log-likelihood of the observed amplitude.
   The method is inspired by Read, Acta Cryst. (1999). D55, 1759-1764.
   Outliers are rejected on the basis of the assumption that the log
   likelihood differnce log[P(Fmode)]-log[P(Fobs)] is distributed
   according to a Chi-square distribution
   (see http://en.wikipedia.org/wiki/Likelihood-ratio_test ).
   The outlier threshold of relates to the p-value of the
   extreme value distribution of the chi-square distribution.

""" % vars())

def run(args, command_name="phenix.remove_outliers"):
  if (len(args)==0 or "--help" in args or "--h" in args or "-h" in args):
    print_help(command_name=command_name)
  else:
    log = multi_out()
    plot_out = None
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = master_params.command_line_argument_interpreter(
      home_scope="outlier_detection")

    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
      if arg=="--quiet":
        arg_is_processed = True
      if (os.path.isfile(arg)): ## is this a file name?
        # check if it is a phil file
        try:
          command_line_params = iotbx.phil.parse(file_name=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except Exception : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except Exception : pass

      if not arg_is_processed:
        print("##----------------------------------------------##", file=log)
        print("## Unknown file or keyword:", arg, file=log)
        print("##----------------------------------------------##", file=log)
        print(file=log)
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()
    if not os.path.exists( params.outlier_utils.input.xray_data.file_name ):
      raise Sorry("File %s can not be found"%(params.outlier_utils.input.xray_data.file_name) )
    if params.outlier_utils.input.model.file_name is not None:
      if not os.path.exists( params.outlier_utils.input.model.file_name ):
        raise Sorry("File %s can not be found"%(params.outlier_utils.input.model.file_name) )



    # now get the unit cell from the pdb file

    hkl_xs = None
    if params.outlier_utils.input.xray_data.file_name is not None:
      hkl_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.outlier_utils.input.xray_data.file_name)
    pdb_xs = None
    if params.outlier_utils.input.model.file_name is not None:
      pdb_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.outlier_utils.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.outlier_utils.input.unit_cell,
      space_group_info=params.outlier_utils.input.space_group  )

    phil_xs.show_summary()
    hkl_xs.show_summary()


    combined_xs = crystal.select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[hkl_xs])

    # inject the unit cell and symmetry in the phil scope please
    params.outlier_utils.input.unit_cell = combined_xs.unit_cell()
    params.outlier_utils.input.space_group = \
      sgtbx.space_group_info( group = combined_xs.space_group() )

    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)

    if params.outlier_utils.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.outlier_utils.input.space_group is None:
      raise Sorry("space group not specified")
    if params.outlier_utils.input.xray_data.file_name is None:
      raise Sorry("Xray data not specified")
    if params.outlier_utils.input.model.file_name is None:
      print("PDB file not specified. Basic wilson outlier rejections only.")



    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #
    phil_xs = crystal.symmetry(
      unit_cell=params.outlier_utils.input.unit_cell,
      space_group_info=params.outlier_utils.input.space_group  )

    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = phil_xs,
      force_symmetry = True,
      reflection_files=[])

    miller_array = None

    miller_array = xray_data_server.get_xray_data(
      file_name = params.outlier_utils.input.xray_data.file_name,
      labels = params.outlier_utils.input.xray_data.obs_labels,
      ignore_all_zeros = True,
      parameter_scope = 'outlier_utils.input.xray_data',
      parameter_name = 'obs_labels'
      )

    info = miller_array.info()

    miller_array = miller_array.map_to_asu()

    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))

    #we have to check if the sigma's make any sense at all
    if not miller_array.sigmas_are_sensible():
      miller_array = miller_array.customized_copy(
        data = miller_array.data(),
        sigmas=None).set_observation_type(miller_array)
    miller_array = miller_array.select(
      miller_array.data() > 0 )
    if  miller_array.sigmas() is not None:
      miller_array = miller_array.select(
        miller_array.sigmas() > 0 )

    if (miller_array.is_xray_intensity_array()):
      miller_array = miller_array.f_sq_as_f()
    elif (miller_array.is_complex_array()):
      miller_array = abs(miller_array)

    miller_array.set_info(info)
    merged_anomalous=False
    if miller_array.anomalous_flag():
      miller_array = miller_array.average_bijvoet_mates().set_observation_type(
        miller_array )
      merged_anomalous=True
    miller_array = miller_array.map_to_asu()

    # get the free reflections please
    free_flags = None
    if params.outlier_utils.input.xray_data.free_flags is None:
      free_flags = miller_array.generate_r_free_flags(
         fraction=params.outlier_utils.\
           additional_parameters.free_flag_generation.fraction,
         max_free=params.outlier_utils.\
           additional_parameters.free_flag_generation.max_number,
         lattice_symmetry_max_delta=params.outlier_utils.\
           additional_parameters.free_flag_generation.lattice_symmetry_max_delta,
         use_lattice_symmetry=params.outlier_utils.\
           additional_parameters.free_flag_generation.use_lattice_symmetry
        )
    else:
      free_flags = xray_data_server.get_xray_data(
        file_name = params.outlier_utils.input.xray_data.file_name,
        labels = params.outlier_utils.input.xray_data.free_flags,
        ignore_all_zeros = True,
        parameter_scope = 'outlier_utils.input.xray_data',
        parameter_name = 'free_flags'
        )

      if free_flags.anomalous_flag():
        free_flags = free_flags.average_bijvoet_mates()
        merged_anomalous=True
      free_flags = free_flags.customized_copy(
        data = flex.bool( free_flags.data() == 1 ))
      free_flags = free_flags.map_to_asu()
      free_flags = free_flags.common_set( miller_array )


    print(file=log)
    print("Summary info of observed data", file=log)
    print("=============================", file=log)
    miller_array.show_summary(f=log)
    if merged_anomalous:
      print("For outlier detection purposes, the Bijvoet pairs have been merged.", file=log)
    print(file=log)

    print("Constructing an outlier manager", file=log)
    print("===============================", file=log)
    print(file=log)
    outlier_manager = outlier_rejection.outlier_manager(
      miller_array,
      free_flags,
      out=log)

    basic_array = None
    extreme_array = None
    model_based_array = None

    basic_array = outlier_manager.basic_wilson_outliers(
      p_basic_wilson = params.outlier_utils.outlier_detection.\
                       parameters.basic_wilson.level,
      return_data = True)

    extreme_array = outlier_manager.extreme_wilson_outliers(
      p_extreme_wilson = params.outlier_utils.outlier_detection.parameters.\
                         extreme_wilson.level,
      return_data = True)

    beamstop_array = outlier_manager.beamstop_shadow_outliers(
      level = params.outlier_utils.outlier_detection.parameters.\
               beamstop.level,
      d_min = params.outlier_utils.outlier_detection.parameters.\
               beamstop.d_min,
      return_data=True)



    #----------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    if params.outlier_utils.input.model.file_name is not None:
      model = pdb.input(file_name=params.outlier_utils.input.model.file_name).xray_structure_simple(
        crystal_symmetry=phil_xs)
      print("Atomic model summary", file=log)
      print("====================", file=log)
      model.show_summary(f=log)
      print(file=log)


      # please make an f_model object for bulk solvent scaling etc etc

      f_model_object = f_model.manager(
        f_obs = miller_array,
        r_free_flags = free_flags,
        xray_structure = model )
      print("Bulk solvent scaling of the data", file=log)
      print("================================", file=log)
      print("Maximum likelihood bulk solvent scaling.", file=log)
      print(file=log)
      f_model_object.update_all_scales(log=log, remove_outliers=False)
      plot_out = StringIO()
      model_based_array = outlier_manager.model_based_outliers(
        f_model_object.f_model(),
        level=params.outlier_utils.outlier_detection.parameters.model_based.level,
        return_data=True,
        plot_out=plot_out)
    #check what needs to be put out please
    if params.outlier_utils.output.hklout is not None:
      if params.outlier_utils.outlier_detection.protocol == "model":
        if params.outlier_utils.input.model.file_name == None:
          print("Model based rejections requested. No model was supplied.", file=log)
          print("Switching to writing out rejections based on extreme value Wilson statistics.", file=log)
          params.outlier_utils.outlier_detection.protocol="extreme"

      output_array = None
      print(file=log)
      if params.outlier_utils.outlier_detection.protocol == "basic":
        print("Non-outliers found by the basic wilson statistics", file=log)
        print("protocol will be written out.", file=log)
        output_array = basic_array
        new_set_of_free_flags = free_flags.common_set( basic_array )

      if params.outlier_utils.outlier_detection.protocol == "extreme":
        print("Non-outliers found by the extreme value wilson statistics", file=log)
        print("protocol will be written out.", file=log)
        output_array = extreme_array
        new_set_of_free_flags = free_flags.common_set( extreme_array )

      if params.outlier_utils.outlier_detection.protocol == "model":
        print("Non-outliers found by the model based", file=log)
        print("protocol will be written out to the file:", file=log)
        print(params.outlier_utils.output.hklout, file=log)
        print(file=log)
        output_array = model_based_array
        new_set_of_free_flags = free_flags.common_set( model_based_array )

      if params.outlier_utils.outlier_detection.protocol == "beamstop":
        print("Outliers found for the beamstop shadow", file=log)
        print("problems detection protocol will be written to the file:", file=log)
        print(params.outlier_utils.output.hklout, file=log)
        print(file=log)
        output_array = model_based_array
        new_set_of_free_flags = free_flags.common_set( model_based_array )

      mtz_dataset = output_array.as_mtz_dataset(
        column_root_label="FOBS")
      mtz_dataset = mtz_dataset.add_miller_array(
        miller_array = new_set_of_free_flags,
        column_root_label = "Free_R_Flag"
        )
      mtz_dataset.mtz_object().write(
        file_name=params.outlier_utils.output.hklout)

    if (params.outlier_utils.output.logfile is not None):
      final_log = StringIO()
      print(string_buffer.getvalue(), file=final_log)
      print(file=final_log)
      if plot_out is not None:
        print(plot_out.getvalue(), file=final_log)
      outfile = open( params.outlier_utils.output.logfile, 'w' )
      outfile.write( final_log.getvalue() )
      print(file=log)
      print("A logfile named %s was created."%(
        params.outlier_utils.output.logfile), file=log)
      print("This logfile contains the screen output and", file=log)
      print("(possibly) some ccp4 style loggraph plots", file=log)
