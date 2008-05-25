# LIBTBX_SET_DISPATCHER_NAME phenix.model_vs_data

import sys, os
from cctbx.array_family import flex
from iotbx import pdb
from cctbx import adptbx
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.f_model
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import crystal_symmetry_from_any
from libtbx import smart_open
from iotbx import reflection_file_utils
from iotbx import reflection_file_reader
from libtbx.str_utils import format_value
import iotbx
from cStringIO import StringIO
from mmtbx import utils

def get_fmodel_object(xray_structure, f_obs, r_free_flags):
  sel = f_obs.d_spacings().data() > 0.25
  f_obs = f_obs.select(sel)
  r_free_flags = r_free_flags.select(sel)
  n_outl = sel.count(False)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = r_free_flags,
    target_name    = "ls_wunit_k1",
    f_obs          = f_obs)
  fmodel.update_xray_structure(update_f_calc = True, update_f_mask = True)
  fmodel.update_solvent_and_scale(verbose = -1)
  sel = fmodel.outlier_selection()
  fmodel = fmodel.select(selection = sel)
  n_outl += sel.count(False)
  return fmodel, n_outl

def get_xray_structure(pdb_file_name):
  pdb_raw_records = smart_open.for_reading(
    file_name=pdb_file_name).read().splitlines()
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv = mon_lib_srv,
    ener_lib    = ener_lib,
    raw_records = pdb_raw_records)
  xray_structure = processed_pdb_file.xray_structure()
  if(xray_structure is None or xray_structure.scatterers().size()==0):
    raise Sorry("Cannot extract xray_structure.")
  return xray_structure

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show(file_name, fmodel, n_outl, test_flag_value):
  uc = fmodel.xray_structure.unit_cell()
  b_cart = " ".join([("%8.2f"%v).strip() for v in fmodel.b_cart()])
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*100./fmodel.r_free_flags.data().size()
  b_isos = fmodel.xray_structure.extract_u_iso_or_u_equiv()
  n_aniso = fmodel.xray_structure.use_u_aniso().count(True)
  n_not_positive_definite = \
    fmodel.xray_structure.is_positive_definite_u().count(False)
  print "Model file: "+format_value("%5s",os.path.basename(file_name))
  print "  unit cell:       ", uc
  print "  space group:     ",fmodel.xray_structure.crystal_symmetry().\
    space_group_info().symbol_and_number()
  print "  Data:"
  result = " \n    ".join([
    "high_res=       "+format_value("%-5.2f",d_min),
    "low_res=        "+format_value("%-6.2f",d_max),
    "R_work=         "+format_value("%-6.4f",fmodel.r_work()),
    "R_free=         "+format_value("%-6.4f",fmodel.r_free()),
    "hkl_size=       "+format_value("%-8d",fmodel.f_obs.data().size()),
    "test_set%=      "+format_value("%-8.1f",flags_pc),
    "test_flag=      "+format_value("%-d",test_flag_value),
    "twin=           "+format_value("%-5s",fmodel.twin_test()),
    "#Fobs_outl=     "+format_value("%-8d",n_outl),
    "anom_flag=      "+format_value("%-6s",fmodel.f_obs.anomalous_flag())])
  print "   ", result
  print "  Model:"
  result = " \n    ".join([
    "#atoms=         "+format_value("%-8d",fmodel.xray_structure.scatterers().size()),
    "k_sol=          "+format_value("%-5.2f",fmodel.k_sol()),
    "b_sol=          "+format_value("%-6.1f",fmodel.b_sol()),
    "b_cart=         "+format_value("%-s",b_cart),
    "wilsonB=        "+format_value("%-6.1f",fmodel.wilson_b()),
    "b_mean=         "+format_value("%-6.1f",adptbx.u_as_b(flex.mean(b_isos))),
    "b_min=          "+format_value("%-6.1f",adptbx.u_as_b(flex.min(b_isos))),
    "b_max=          "+format_value("%-6.1f",adptbx.u_as_b(flex.max(b_isos))),
    "#aniso=         "+format_value("%-7s",n_aniso),
    "#not pos. def.= "+format_value("%-8s",n_not_positive_definite)])
  print "   ", result

def run(args, command_name = "phenix.model_vs_data"):
  if(len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="%s reflection_file pdb_file [options]" % command_name,
    description='Example: %s data.mtz model.pdb'%command_name)
    .option(None, "--f_obs_label",
      action="store",
      default=None,
      type="string",
      help="Label for F-obs (or I-obs).")
    .option(None, "--r_free_flags_label",
      action="store",
      default=None,
      type="string",
      help="Label for free R flags.")
    ).process(args=args)
  crystal_symmetry = None
  hkl_file_name = None
  pdb_file_name = None
  reflection_file = None
  for arg in command_line.args:
    arg_is_processed = False
    if(not os.path.isfile(arg)):
      raise Sorry("The command line argument %s is not a file."%arg)
    else:
      if(pdb.is_pdb_file(file_name=arg)):
        pdb_file_name = arg
        arg_is_processed = True
        try:
          crystal_symmetry = crystal_symmetry_from_pdb.extract_from(
            file_name=arg)
        except RuntimeError, e:
          if(str(e) == "No CRYST1 record."): pass
      if(not arg_is_processed):
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          arg_is_processed = True
          crystal_symmetry = crystal_symmetry_from_any.extract_from(arg)
      if(not arg_is_processed):
        raise Sorry(
          "The command line argument %s is not a valid PDB or data file."%arg)
  if(crystal_symmetry is None): raise Sorry("Crystal symmetry is not defined.")
  if(pdb_file_name is None): raise Sorry("No PDB file given.")
  if(reflection_file is None): raise Sorry("No reflection file given.")
  xray_structure = get_xray_structure(pdb_file_name = pdb_file_name)
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = [reflection_file])
  parameters = utils.data_and_flags.extract()
  if(command_line.options.f_obs_label is not None):
    parameters.labels = command_line.options.f_obs_label
  if(command_line.options.r_free_flags_label is not None):
    parameters.r_free_flags.label = command_line.options.r_free_flags_label
  determine_data_and_flags_result = utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    parameters              = parameters,
    data_parameter_scope    = "refinement.input.xray_data",
    flags_parameter_scope   = "refinement.input.xray_data.r_free_flags",
    data_description        = "X-ray data",
    keep_going              = True,
    log                     = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=1
  if(not xray_structure.crystal_symmetry().is_similar_symmetry(
     f_obs.crystal_symmetry())):
    raise Sorry("Inconsistent crystal symmetry.")
  fmodel, n_outl = get_fmodel_object(xray_structure = xray_structure,
                                     f_obs          = f_obs,
                                     r_free_flags   = r_free_flags)
  show(file_name       = pdb_file_name,
       fmodel          = fmodel,
       n_outl          = n_outl,
       test_flag_value = test_flag_value)

if(__name__ == "__main__"):
   run(sys.argv[1:])
