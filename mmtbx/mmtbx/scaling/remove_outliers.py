import math
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
from cctbx import adptbx
from scitbx import lbfgs
from mmtbx import masks
import cctbx.xray.structure_factors
from cctbx.eltbx.xray_scattering import wk1995
from libtbx import adopt_init_args
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx.scaling import twin_detwin_data
from mmtbx.scaling import outlier_rejection
from mmtbx import f_model
import libtbx.phil.command_line
from mmtbx.twinning import twin_f_model
from cStringIO import StringIO
from scitbx.python_utils import easy_pickle
import sys, os

def select_crystal_symmetry(
      from_command_line,
      from_parameter_file,
      from_coordinate_files,
      from_reflection_files):
  result = crystal.symmetry(
    unit_cell=None,
    space_group_info=None)
  if (from_command_line is not None):
    result = result.join_symmetry(
      other_symmetry=from_command_line, force=False)
  if (from_parameter_file is not None):
    result = result.join_symmetry(
      other_symmetry=from_parameter_file, force=False)
  if (result.unit_cell() is None):
    for crystal_symmetry in from_reflection_files:
      unit_cell = crystal_symmetry.unit_cell()
      if (unit_cell is not None):
        result = crystal.symmetry(
          unit_cell=unit_cell,
          space_group_info=result.space_group_info(),
          assert_is_compatible_unit_cell=False)
        break
  for crystal_symmetry in from_coordinate_files:
    result = result.join_symmetry(other_symmetry=crystal_symmetry, force=False)
  if (result.space_group_info() is None):
    for crystal_symmetry in from_reflection_files:
      space_group_info = crystal_symmetry.space_group_info()
      if (space_group_info is not None):
        result = crystal.symmetry(
          unit_cell=result.unit_cell(),
          space_group_info=space_group_info,
          assert_is_compatible_unit_cell=False)
        break
  return result



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
    protocol=basic *extreme model
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

def print_help():
  print """
usage: mmtbx.remove_outliers <options>

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
    protocol=basic *extreme model    << outlier protocol. See below.
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

mmtbx.remove_outliers data.file=my_precious_data.mtz data.obs=FOBS data.free=TEST \
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

  """






def run(command_name, args):
  if len(args)==0:
    print_help()
  elif ( "--help" in args ):
    print_help()
  elif ( "--h" in args ):
    print_help()
  elif ("-h" in args ):
    print_help()
  else:
    log = multi_out()
    if (not "--quiet" in args):
      log.register(label="stdout", file_object=sys.stdout)
    string_buffer = StringIO()
    string_buffer_plots = StringIO()
    log.register(label="log_buffer", file_object=string_buffer)

    phil_objects = []
    argument_interpreter = libtbx.phil.command_line.argument_interpreter(
      master_params=master_params,
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
        except : pass
      else:
        try:
          command_line_params = argument_interpreter.process(arg=arg)
          if command_line_params is not None:
            phil_objects.append(command_line_params)
            arg_is_processed = True
        except KeyboardInterrupt: raise
        except : pass

      if not arg_is_processed:
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params = effective_params.extract()

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


    combined_xs = select_crystal_symmetry(
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
      print "PDB file not specified. Basic wilson outlier rejections only."



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


    print >> log
    print >> log, "Summary info of observed data"
    print >> log, "============================="
    miller_array.show_summary(f=log)
    if merged_anomalous:
      print >> log, "For outlier detection purposes, the Bijvoet pairs have been merged."
    print >> log

    print >> log, "Constructing an outlier manager"
    print >> log, "==============================="
    print >> log
    outlier_manager = outlier_rejection.outlier_manager(
      miller_array, out=log)

    basic_array = None
    extreme_array = None
    model_based_array = None

    basic_array = outlier_manager.basic_wilson_outliers(
      p_basic_wilson = params.outlier_utils.outlier_detection.\
                       parameters.basic_wilson.level,
      return_array = True)

    extreme_array = outlier_manager.extreme_wilson_outliers(
      p_extreme_wilson = params.outlier_utils.outlier_detection.parameters.\
                         extreme_wilson.level,
      return_array = True)

    #----------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    if params.outlier_utils.input.model.file_name is not None:
      model = xray_structure.from_pdb(
        file_name=params.outlier_utils.input.model.file_name,
        crystal_symmetry=phil_xs,
        force_symmetry=True)
      print >> log, "Atomic model summary"
      print >> log, "===================="
      model.show_summary(f=log)
      print >> log


      # please make an f_model object for bulk solvent scaling etc etc

      f_model_object = f_model.manager(
        f_obs = miller_array,
        r_free_flags = free_flags,
        xray_structure = model )
      print >> log, "Bulk solvent scaling of the data"
      print >> log, "================================"
      print >> log, "Maximum likelihood bulk solvent scaling."
      print >> log
      f_model_object.update_solvent_and_scale(out=log)
      b_cart = f_model_object.b_cart()
      k_sol = f_model_object.k_sol()
      b_sol = f_model_object.b_sol()
      ls_scale = 1.0/f_model_object.scale_k1()
      print >> log
      print >> log, "The observed data is scaled by a multiplier"
      print >> log, "equal to %5.2e"%(ls_scale)
      print >> log, "This brings the data to an approximate absolute scale."

      # update the outlier object please
      outlier_manager.apply_scale_to_original_data( ls_scale)
      free_flags = free_flags.common_set( outlier_manager.miller_obs )

      # redo the f model object please
      f_model_object = f_model.manager(
        f_obs = outlier_manager.miller_obs,
        r_free_flags = free_flags,
        xray_structure = model)
      # reset the bulk solvent parameters please
      f_model_object.update_core(b_cart=b_cart,
                                 k_sol=k_sol,
                                 b_sol=b_sol)
      f_model_data = f_model_object.f_model()

      plot_out = StringIO()
      # get alphas and betas please
      alpha,beta = f_model_object.alpha_beta()
      # get suspected outliers
      model_based_array = outlier_manager.model_based_outliers(
        f_model_data,
        alpha,
        beta,
        level=params.outlier_utils.outlier_detection.parameters.model_based.level,
        return_array=True,
        plot_out=plot_out)

    #check what needs to be put out please
    if params.outlier_utils.output.hklout is not None:
      if params.outlier_utils.outlier_detection.protocol == "model":
        if params.outlier_utils.input.model.file_name == None:
          print >> log, "Model based rejections requested. No model was supplied."
          print >> log, "Switching to writing out rejections based on extreme value Wilson statistics."
          params.outlier_utils.outlier_detection.protocol="extreme"

      output_array = None
      print >> log
      if params.outlier_utils.outlier_detection.protocol == "basic":
        print >> log, "Outliers found by the basic wilson statistics"
        print >> log, "protocol will be written out."
        output_array = basic_array
        new_set_of_free_flags = free_flags.common_set( basic_array )

      if params.outlier_utils.outlier_detection.protocol == "extreme":
        print >> log, "Outliers found by the extreme value wilson statistics"
        print >> log, "protocol will be written out."
        output_array = extreme_array
        new_set_of_free_flags = free_flags.common_set( extreme_array )

      if params.outlier_utils.outlier_detection.protocol == "model":
        print >> log, "Outliers found by the model based"
        print >> log, "protocol will be written out to the file:"
        print >> log, params.outlier_utils.output.hklout
        print >> log
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

    if params.outlier_utils.output.logfile is not None:
      final_log = StringIO()
      print >> final_log, string_buffer.getvalue()
      print >> final_log
      print >> final_log, plot_out.getvalue()
      outfile = open( params.outlier_utils.output.logfile, 'w' )
      outfile.write( final_log.getvalue() )
      print >> log
      print >> log, "A logfile named %s was created."%(
        params.outlier_utils.output.logfile)
      print >> log, "This logfile contains the screen output and"
      print >> log, "(possibly) some ccp4 style loggraph plots"


if (__name__ == "__main__" ):
  outlier_utils(sys.argv[1:])
