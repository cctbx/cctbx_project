from cctbx import crystal
from cctbx import sgtbx
import cctbx.xray.structure_factors
from libtbx.utils import Sorry, multi_out
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from iotbx.pdb import xray_structure
from iotbx import pdb
from cStringIO import StringIO
from mmtbx import f_model
import sys, os




master_params = iotbx.phil.parse("""\
sfcalc{
  input{
    unit_cell=None
    .type=unit_cell
    space_group=None
    .type=space_group
    model{
      file_name=None
      .type=path
    }
  }
  parameters{
    d_min = 2.0
    .type=float
    overall{
      b_cart{
        b_11 = 0
        .type=float
        b_22 = 0
        .type=float
        b_33 = 0
        .type=float
        b_12 = 0
        .type=float
        b_13 = 0
        .type=float
        b_23 = 0
        .type=float
      }
      k_overall = 0.1
      .type=float
    }
    solvent{
      k_sol = 0.3
      .type=float
      b_sol = 56.0
      .type=float
    }
    output_type = complex *amplitudes intensities
    .type=choice
  }
  output{
    logfile=sfcalc.log
    .type=str
    hklout=sfcalc.mtz
    .type=str
  }
}
""")

def print_help():
  print """
No help
  """


def sfcalc(args):
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
      home_scope="sfcalc")

    print >> log, "#phil __OFF__"
    print >> log, "================="
    print >> log, "     SFCALC      "
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
    params = effective_params.extract()

    # now get the unit cell from the files
    hkl_xs = None
    pdb_xs = None

    if params.sfcalc.input.model.file_name is not None:
      pdb_xs = crystal_symmetry_from_any.extract_from(
        file_name=params.sfcalc.input.model.file_name)

    phil_xs = crystal.symmetry(
      unit_cell=params.sfcalc.input.unit_cell,
      space_group_info=params.sfcalc.input.space_group  )

    combined_xs = crystal.select_crystal_symmetry(
      None,phil_xs, [pdb_xs],[None])
    combined_xs.show_summary()
    if combined_xs.unit_cell() is None:
      raise Sorry("Unit cell not defined")
    if combined_xs.space_group() is None:
      raise Sorry("Space group not defined")


    # inject the unit cell and symmetry in the phil scope please
    params.sfcalc.input.unit_cell = combined_xs.unit_cell()
    params.sfcalc.input.space_group = \
      sgtbx.space_group_info( group = combined_xs.space_group() )

    print >> log, "#phil __ON__"
    new_params =  master_params.format(python_object=params)
    new_params.show(out=log)
    print >> log, "#phil __END__"

    pdb_model = None

    if params.sfcalc.input.model.file_name is not None:
      pdb_model = pdb.input(file_name=params.sfcalc.input.model.file_name)
      model = pdb_model.xray_structure_simple(crystal_symmetry=phil_xs)
      print >> log, "Atomic model summary"
      print >> log, "===================="
      model.show_summary()
      print >> log

      #make an f_model object please
      b_cart = params.sfcalc.parameters.overall.b_cart
      b_cart = [b_cart.b_11,
                b_cart.b_22,
                b_cart.b_33,
                b_cart.b_12,
                b_cart.b_13,
                b_cart.b_23]
      dummy = abs(model.structure_factors(
        d_min          = params.sfcalc.parameters.d_min,
        anomalous_flag = False).f_calc())

      flags = dummy.generate_r_free_flags(fraction = 0.1,
                                          max_free = 99999999)

      fmodel = f_model.manager( xray_structure   = model,
                                r_free_flags     = flags,
                                target_name      = "ls_wunit_k1",
                                f_obs            = dummy,
                                b_cart           = b_cart,
                                k_sol            = params.sfcalc.parameters.solvent.k_sol,
                                b_sol            = params.sfcalc.parameters.solvent.b_sol )

      calc_data_with_solvent_contrib = fmodel.f_model()
      calc_data_with_solvent_contrib = calc_data_with_solvent_contrib.array(
        data=calc_data_with_solvent_contrib.data()*params.sfcalc.parameters.overall.k_overall)
      result = None
      label = None
      if params.sfcalc.parameters.output_type == "complex":
        result = calc_data_with_solvent_contrib
        label="FMODEL"
      if params.sfcalc.parameters.output_type == "amplitudes":
        result = abs(calc_data_with_solvent_contrib)
        label="FMODEL"
      if params.sfcalc.parameters.output_type == "intensities":
        result = abs(calc_data_with_solvent_contrib)
        result = result.f_as_f_sq()
        label="IMODEL"


      #write an mtz file with the data
      mtz_dataset = result.as_mtz_dataset(
        column_root_label=label)

      mtz_dataset.mtz_object().write(
        file_name=params.sfcalc.output.hklout)


    #write the logfile
    logger = open( params.sfcalc.output.logfile, 'w')
    print >> log, "writing log file with name %s"%(params.sfcalc.output.logfile)
    print >> log
    print >> logger, string_buffer.getvalue()



if (__name__ == "__main__" ):
  sfcalc(sys.argv[1:])
