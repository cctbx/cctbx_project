from mmtbx.xmanip import write_as_pdb_file
import iotbx.pdb
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
import iotbx.phil
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.math import matrix
from libtbx.utils import Sorry, multi_out
from cStringIO import StringIO
import sys, os

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
    action = *reindex operator manipulate_pdb
    .type=choice
    chain_id_increment = 0
    .type=int
    inverse=False
    .type=bool
    reindex{
      standard_laws = niggli *reference_setting invert user_supplied
      .type=choice
      user_supplied_law='h,k,l'
      .type=str
    }
    apply_operator{
      standard_operators = *user_supplied
      .type=choice
      user_supplied_operator="x,y,z"
      .type=str
      concatenate_model=False
      .type=bool
    }
    manipulate_pdb{
      set_b = True
      .type=bool
      b_iso = 30
      .type=float
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

Allows one to quickly reindex a dataset and apply the effect on the
atomic coordinates as well.

The keywords are sumarized below:
"""
  master_params.show()

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
    argument_interpreter = master_params.command_line_argument_interpreter(
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
        print >> log, "##----------------------------------------------##"
        print >> log, "## Unknown file or keyword:", arg
        print >> log, "##----------------------------------------------##"
        print >> log
        raise Sorry("Unknown file or keyword: %s" % arg)

    effective_params = master_params.fetch(sources=phil_objects)
    params_root = effective_params.extract()
    params = params_root.reindex_utils

    # now get the unit cell from the files
    hkl_xs = None
    pdb_xs = None
    if params.input.xray_data.file_name is not None:
      hkl_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.input.xray_data.file_name)
    if params.input.model.file_name is not None:
      pdb_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.input.unit_cell,
      space_group_info=params.input.space_group  )


    combined_xs = crystal.select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[hkl_xs])

    # inject the unit cell and symmetry in the phil scope please
    params.input.unit_cell = combined_xs.unit_cell()
    params.input.space_group = combined_xs.space_group_info()

    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params_root)
    new_params.show(out=log)
    print >> log, "#phil __END__"

    if params.input.unit_cell is None:
      raise Sorry("unit cell not specified")
    if params.input.space_group is None:
      raise Sorry("space group not specified")

    #-----------------------------------------------------------
    #
    # step 1: read in the reflection file
    #
    miller_array = None
    if  params.input.xray_data.file_name is not None:
      phil_xs = crystal.symmetry(
        unit_cell=params.input.unit_cell,
        space_group_info=params.input.space_group  )

      xray_data_server =  reflection_file_utils.reflection_file_server(
        crystal_symmetry = phil_xs,
        force_symmetry = True,
        reflection_files=[])

      miller_array = xray_data_server.get_xray_data(
        file_name = params.input.xray_data.file_name,
        labels = params.input.xray_data.labels,
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
    pdb_model = None

    if params.input.model.file_name is not None:
      pdb_model = iotbx.pdb.input(
        file_name=params.input.model.file_name)
      model = pdb_model.xray_structure_simple(crystal_symmetry=phil_xs)
      print >> log, "Atomic model summary"
      print >> log, "===================="
      model.show_summary()
      print >> log

    if params.parameters.action=="reindex":
      #----------------------------------------------------------------
      # step 3: get the reindex laws
      to_niggli    = phil_xs.change_of_basis_op_to_niggli_cell()
      to_reference = phil_xs.change_of_basis_op_to_reference_setting()
      to_inverse   = phil_xs.change_of_basis_op_to_inverse_hand()
      cb_op = None
      pr = params.parameters.reindex
      if (pr.standard_laws == "niggli"):
        cb_op = to_niggli
      elif (pr.standard_laws == "reference_setting"):
        cb_op = to_reference
      elif (pr.standard_laws == "invert"):
        cb_op = to_inverse
      else:
        cb_op = sgtbx.change_of_basis_op(pr.user_supplied_law)

      if cb_op is None:
        raise Sorry("No change of basis operation is supplied.")
      if params.parameters.inverse:
        cb_op = cb_op.inverse()

      print >> log, "Supplied reindexing law:"
      print >> log, "========================"
      print >> log, "hkl notation: ", cb_op.as_hkl()
      print >> log, "xyz notation: ", cb_op.as_xyz()
      print >> log, "abc notation: ", cb_op.as_abc()
      #----------------------------------------------------------------
      # step 4: do the reindexing
      #
      # step 4a: first do the miller array object
      new_miller_array = None
      if miller_array is not None:
        new_miller_array = miller_array.change_basis( cb_op )
      #
      # step 4b: the xray structure
      new_model = None
      if pdb_model is not None:
        new_model = model.change_basis( cb_op )

      #----------------------------------------------------------------
      # step 5a: write the new mtz file
      print >> log
      print >> log, "The data and model have been reindexed"
      print >> log, "--------------------------------------"
      print >> log
      print >> log, "Writing output files...."
      if miller_array is not None:
        print >> log, "writing mtz file with name %s" % (params.output.hklout)
        mtz_dataset = new_miller_array.as_mtz_dataset(
          column_root_label="FOBS")
        mtz_dataset.mtz_object().write(file_name=params.output.hklout)

      #step 5b: write the new pdb file
      if new_model is not None:
        pdb_file = open(params.output.xyzout, 'w')
        print >> log, "Wring pdb file to: %s" % params.output.xyzout
        write_as_pdb_file(
          input_pdb = pdb_model,
          input_xray_structure = new_model,
          out = pdb_file,
          chain_id_increment = params.parameters.chain_id_increment,
          additional_remark = "Generated by mmtbx reindex")

        print >> pdb_file, "END"
        pdb_file.close()
      if ( [miller_array,new_model]).count(None)==2:
        print >>log, "No input reflection of coordinate files have been given"

    if params.parameters.action=="operator":
      rt_mx = sgtbx.rt_mx(
        params.parameters.apply_operator.user_supplied_operator,t_den=12*8 )
      if params.parameters.inverse:
        rt_mx = rt_mx.inverse()
      print >> log
      print >> log, "Applied operator : ", rt_mx.as_xyz()
      print >> log

      sites = model.sites_frac()
      new_sites = flex.vec3_double()
      for site in sites:
        new_site = rt_mx.r()*matrix.col(site)
        new_site = flex.double(new_site)+flex.double( rt_mx.t().as_double() )
        new_sites.append( tuple(new_site) )
      new_model = model.deep_copy_scatterers()

      new_model.set_sites_frac( new_sites )
      # write the new [pdb file please
      pdb_file = open(params.output.xyzout, 'w')
      print >> log, "Wring pdb file to: %s" % params.output.xyzout
      if params.parameters.apply_operator.concatenate_model:
        write_as_pdb_file(
          input_pdb = pdb_model,
          input_xray_structure = model,
          out = pdb_file,
          chain_id_increment = 0,
          additional_remark = None,
          print_cryst_and_scale=True)

      write_as_pdb_file(
        input_pdb = pdb_model,
        input_xray_structure = new_model,
        out = pdb_file,
        chain_id_increment = params.parameters.chain_id_increment,
        additional_remark = None,
        print_cryst_and_scale=False)

      print >> pdb_file, "END"
      pdb_file.close()

    if params.parameters.action=="manipulate_pdb":
      #rest all the b values
      if params.parameters.manipulate_pdb.set_b:
        b_iso = params.reindex_utils.parameters.manipulate_pdb.b_iso
        new_model = model.set_b_iso( value = b_iso )
        print >> log
        print >> log, "All B-values have been set to %5.3f"%(b_iso)
        print >> log, "Writing PDB file %s"%(params.output.xyzout)
        print >> log

      pdb_file = open(params.output.xyzout, 'w')
      write_as_pdb_file(
        input_pdb = pdb_model,
        input_xray_structure = new_model,
        out = pdb_file,
        chain_id_increment = 0,
        additional_remark = None,
        print_cryst_and_scale=True)
      print >> pdb_file, "END"
      pdb_file.close()

    #write the logfile
    logger = open(params.output.logfile, 'w')
    print >> log, "Writing log file with name %s" % params.output.logfile
    print >> log
    print >> logger, string_buffer.getvalue()

if (__name__ == "__main__" ):
  reindex_utils(sys.argv[1:])
