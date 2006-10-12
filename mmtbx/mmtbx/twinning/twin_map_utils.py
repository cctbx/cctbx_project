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
from iotbx import pdb
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews, twin_analyses
from mmtbx.scaling import basic_analyses, pair_analyses
from mmtbx.scaling import twin_detwin_data
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

  tmp = [from_command_line, from_parameter_file]+from_coordinate_files \
        +from_reflection_files
  if tmp.count(None)==len(tmp):
    raise Sorry("No unit cell or symmetry information available. Check input")

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
twin_utils{
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
      free_flag=None
      .type=str
    }
    model{
      file_name=None
      .type=path
    }
  }
  parameters{
    twinning{
      twin_law=None
      .type=None
      max_delta=3.0
      .type=float
    }
  }
  output{
    logfile=twin_tools.log
    .type=str
    map_coeffs_root=MAP_COEFFS
    .type=str
  }
}
""")

def print_help():
  print """

mmtbx.twin_map_utils
--------------------

A command line utility to compute map coefficients for twinned data.
Bulk solvent parameters and twin fractions are determined automatically.
If no twin law is specified, map coefficents for all twin laws will
be computed.

The keywords are sumarized below:

twin_utils{
  input{
    unit_cell=None
    space_group=None
    xray_data{
      file_name=None
      obs_labels=None
      free_flag=None
      }
    model{
      file_name=None
      }
  }
  parameters{
    twinning{
      twin_law=None
      max_delta=3.0
    }
  }
  output{
    logfile=twin_tools.log
    map_coeffs_root=MAP_COEFFS
  }
}

A typical run looks like this:

mmtbx.twin_map_utils data.file=mydata.mtz model.file=mymodel.pdb \\
  obs_labels=FP,SIGFP free_flag=TEST

If no unit cell is specified, the unit cell of the reflection
file or the model will be used.


"""


def run(command_name,args):
  print
  log=sys.stdout
  params=None
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
      home_scope="map_coefs")

    for arg in args:
      command_line_params = None
      arg_is_processed = False
      # is it a file?
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
    """
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    """
    # now get the unit cell from the pdb file

    hkl_xs = crystal_symmetry_from_any.extract_from(
      file_name=params.twin_utils.input.xray_data.file_name)
    pdb_xs = crystal_symmetry_from_any.extract_from(
      file_name=params.twin_utils.input.model.file_name)

    phil_xs=None
    if ([params.twin_utils.input.unit_cell,
         params.twin_utils.input.space_group]).count(None)<2:
      phil_xs = crystal.symmetry(
        unit_cell=params.twin_utils.input.unit_cell,
        space_group_info=params.twin_utils.input.space_group  )


    combined_xs = crystal.select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[hkl_xs])

    # inject the unit cell and symmetry in the phil scope please
    params.twin_utils.input.unit_cell = combined_xs.unit_cell()
    params.twin_utils.input.space_group = \
      sgtbx.space_group_info( group = combined_xs.space_group() )

    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)

    if params.twin_utils.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.twin_utils.input.space_group is None:
      raise Sorry("space group not specified")
    if params.twin_utils.input.xray_data.file_name is None:
      raise Sorry("Xray data not specified")
    if params.twin_utils.input.model.file_name is None:
      raise Sorry("pdb file with  model not specified")

    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #
    phil_xs = crystal.symmetry(
      unit_cell=params.twin_utils.input.unit_cell,
      space_group_info=params.twin_utils.input.space_group  )

    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = phil_xs,
      force_symmetry = True,
      reflection_files=[])

    miller_array = None

    miller_array = xray_data_server.get_xray_data(
      file_name = params.twin_utils.input.xray_data.file_name,
      labels = params.twin_utils.input.xray_data.obs_labels,
      ignore_all_zeros = True,
      parameter_scope = 'twin_utils.input.xray_data',
      parameter_name = 'obs_labels'
      )

    info = miller_array.info()

    miller_array = miller_array.map_to_asu()

    miller_array = miller_array.select(
      miller_array.indices() != (0,0,0))

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

    free_flags = None
    if params.twin_utils.input.xray_data.free_flag is not None:
      free_flags = xray_data_server.get_xray_data(
        file_name = params.twin_utils.input.xray_data.file_name,
        labels = params.twin_utils.input.xray_data.free_flag,
        ignore_all_zeros = False,
        parameter_scope = 'twin_utils.input.xray_data',
        parameter_name = 'free_r_flags'
        )
      if free_flags.anomalous_flag():
        free_flags = free_flags.average_bijvoet_mates()
        merged_anomalous=True
      free_flags = free_flags.customized_copy(
        data = flex.bool( free_flags.data() == 1 ))
      free_flags = free_flags.map_to_asu()
      free_flags = free_flags.common_set( miller_array )
    else:
      free_flags = miller_array.generate_r_free_flags(use_lattice_symmetry=True)




    print >> log
    print >> log, "Summary info of observed data"
    print >> log, "============================="
    miller_array.show_summary(f=log)
    print >> log


    #----------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    model = pdb.input(file_name=params.twin_utils.input.model.file_name).xray_structure_simple(
      crystal_symmetry=phil_xs)
    print >> log, "Atomic model summary"
    print >> log, "===================="
    model.show_summary(f=log)
    print >> log


    #----------------------------------------------------------------
    # step 3: get the twin laws for this xs
    twin_laws = twin_analyses.twin_laws(
      miller_array,
      lattice_symmetry_max_delta=\
         params.twin_utils.parameters.twinning.max_delta,
      out=log)

    print >> log
    print >> log, "Preliminary data analyses"
    print >> log, "=========================="
    twin_laws.show(out=log)


    #---------
    # step 3:
    # make twin model managers for all twin laws
    print >> log
    print >> log, "Overall and bulk solvent scale paranmeters and twin fraction estimation"
    print >> log, "======================================================================="
    twin_models = []
    operator_count = 0
    for twin_law in twin_laws.operators:
      operator_count += 1
      operator_hkl = sgtbx.change_of_basis_op( twin_law.operator ).as_hkl()
      twin_model = twin_f_model.twin_model_manager(
        f_obs_array=miller_array,
        free_array = free_flags,
        xray_structure=model,
        twin_law = twin_law.operator,
        out=log)
      print >> log, "--- bulk solvent scaling ---"
      twin_model.update_bulk_solvent_parameters(de_search=True,
                                                refine=False)

      twin_model.r_values()
      twin_model.target()
      twin_model.twin_fraction_scan()


      tfofc,wtfofc,grad = twin_model.map_coefficients(False)

      mtz_dataset = tfofc.as_mtz_dataset(
        column_root_label="2FOFC")
      mtz_dataset = mtz_dataset.add_miller_array(
        miller_array = wtfofc,
        column_root_label = "FWT"
        )
      mtz_dataset = mtz_dataset.add_miller_array(
        miller_array = grad,
        column_root_label = "GRAD"
        )
      name = params.twin_utils.output.map_coeffs_root+"_"+str(operator_count)+".mtz"
      print "writing %s for twin law %s"%(name,operator_hkl)
      mtz_dataset.mtz_object().write(
        file_name=name)

    print >> log
    print >> log
    print >> log, "All done "
    print >> log
    logfile = open(params.twin_utils.output.logfile,'w')
    print >> logfile,  string_buffer.getvalue()




if (__name__ == "__main__" ):
  run(sys.argv[0],sys.argv[1:])
