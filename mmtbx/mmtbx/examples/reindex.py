from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import xray
from cctbx import eltbx
import cctbx.xray.structure_factors
from cctbx.eltbx.xray_scattering import wk1995
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
import libtbx.phil.command_line
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
reindex_utils{
  input{
    unit_cell=None
    .type=unit_cell
    space_group=None
    .type=space_group
    xray_data{
      file_name=None
      .type=path
      labels=None
      .type=str
    }
    model{
      file_name=None
      .type=path
    }
  }
  parameters{
    reindex{
      standard_laws = niggli *reference_setting invert user_supplied
      .type=choice
      user_supplied_law='h,k,l'
      .type=str
    }
  }
  output{
    logfile=reindex.log
    .type=str
    hklout=reindexed.mtz
    .type=str
    xyzout=reindexmodel.pdb
    .type=str
  }
}
""")

def print_help():
  print """
mmtbx.reindex:

Allows one to quickly reindex a dataset and apply the effect on the atomic coordinates as well.

The keywords are sumarized below:

reindex_utils{
  input{
    unit_cell=None
    space_group=None
    xray_data{
      file_name=None
      labels=None
    }
    model{
      file_name=None
    }
  }
  parameters{
    reindex{
      standard_laws = niggli *reference_setting invert user_supplied
      user_supplied_law='h,k,l'
    }
  }
  output{
    logfile=reindex.log
    hklout=reindexed.mtz
    xyzout=reindexmodel.pdb
  }


  """


def reindex_utils(args):
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

    print >> log, "#phil __OFF__"
    print >> log, "================="
    print >> log, "    REINDEX      "
    print >> log, "A reindexing tool"
    print >> log, "================="
    print >> log


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
      file_name=params.reindex_utils.input.xray_data.file_name)
    pdb_xs = crystal_symmetry_from_any.extract_from(
      file_name=params.reindex_utils.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.reindex_utils.input.unit_cell,
      space_group_info=params.reindex_utils.input.space_group  )


    combined_xs = select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[hkl_xs])

    # inject the unit cell and symmetry in the phil scope please
    params.reindex_utils.input.unit_cell = combined_xs.unit_cell()
    params.reindex_utils.input.space_group = \
      sgtbx.space_group_info( group = combined_xs.space_group() )

    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    print >> log, "#phil __END__"

    if params.reindex_utils.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.reindex_utils.input.space_group is None:
      raise Sorry("space group not specified")
    if params.reindex_utils.input.xray_data.file_name is None:
      raise Sorry("Xray data not specified")
    if params.reindex_utils.input.model.file_name is None:
      raise Sorry("pdb file with  model not specified")

    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #
    phil_xs = crystal.symmetry(
      unit_cell=params.reindex_utils.input.unit_cell,
      space_group_info=params.reindex_utils.input.space_group  )

    xray_data_server =  reflection_file_utils.reflection_file_server(
      crystal_symmetry = phil_xs,
      force_symmetry = True,
      reflection_files=[])

    miller_array = None

    miller_array = xray_data_server.get_xray_data(
      file_name = params.reindex_utils.input.xray_data.file_name,
      labels = params.reindex_utils.input.xray_data.labels,
      ignore_all_zeros = True,
      parameter_scope = 'reindex_utils.input.xray_data',
      parameter_name = 'labels'
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
    print >> log
    print >> log, "Summary info of observed data"
    print >> log, "============================="
    miller_array.show_summary(f=log)
    print >> log


    #----------------------------------------------------------------
    # Step 2: get an xray structure from the PDB file
    #
    model = xray_structure.from_pdb(
      file_name=params.reindex_utils.input.model.file_name,
      crystal_symmetry=phil_xs,
      force_symmetry=True)
    print >> log, "Atomic model summary"
    print >> log, "===================="
    model.show_summary()
    print >> log


    #----------------------------------------------------------------
    # step 3: get the reindex laws

    xs           = crystal.symmetry(unit_cell = model.unit_cell(),
                                    space_group = model.space_group() )
    to_niggli    = xs.change_of_basis_op_to_niggli_cell()
    to_reference = xs.change_of_basis_op_to_reference_setting()
    to_inverse   = xs.change_of_basis_op_to_inverse_hand()
    cb_op = None
    if (params.reindex_utils.parameters.reindex.standard_laws == "niggli"):
      cb_op = to_niggli
    if (params.reindex_utils.parameters.reindex.standard_laws == "reference_setting"):
      cb_op = to_reference
    if (params.reindex_utils.parameters.reindex.standard_laws == "invert"):
      cb_op = to_inverse
    if (params.reindex_utils.parameters.reindex.standard_laws == "user_supplied"):
      cb_op = sgtbx.change_of_basis_op( params.reindex_utils.parameters.reindex.law )

    if cb_op is None:
      raise Sorry("No change of basis operation is supplied.")

    print >> log, "Supplied reindexing law:"
    print >> log, "========================"
    print >> log, cb_op.as_hkl()


    #----------------------------------------------------------------
    # step 4: do the reindexing
    #
    # step 4a: first do the miller array object
    new_miller_array = miller_array.change_basis( cb_op )
    #
    # step 4b: the xray structure
    new_model = model.change_basis( cb_op )

    #----------------------------------------------------------------
    # step 5a: write the new mtz file
    print >> log
    print >> log, "The data and model have been reindexed"
    print >> log, "--------------------------------------"
    print >> log
    print >> log, "Writing output files...."
    print >> log, "writing mtz file with name %s"%(params.reindex_utils.output.hklout)
    mtz_dataset = new_miller_array.as_mtz_dataset(
      column_root_label="FOBS")
    """
    #This is how to add (possibly) a free r flag column that has been reindexed as well
    mtz_dataset = mtz_dataset.add_miller_array(
      miller_array = new_set_of_free_flags,
      column_root_label = "Free_R_Flag"
      )
    """
    mtz_dataset.mtz_object().write(
      file_name=params.reindex_utils.output.hklout)

    #step 5b: write the new pdb file
    print >> log, "writing pdb file with name %s"%(params.reindex_utils.output.xyzout)
    pdb = new_model.as_pdb_file()
    pdb_file = open( params.reindex_utils.output.xyzout, 'w')
    print >> pdb_file, pdb
    pdb_file.close()

    #write the logfile
    logger = open( params.reindex_utils.output.logfile, 'w')
    print >> log, "writing log file with name %s"%(params.reindex_utils.output.logfile)
    print >> log
    print >> logger, string_buffer.getvalue()



if (__name__ == "__main__" ):
  reindex_utils(sys.argv[1:])
