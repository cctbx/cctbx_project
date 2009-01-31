import sys, os, time
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
from mmtbx import utils
from iotbx import pdb
from libtbx import easy_run
from cStringIO import StringIO
from mmtbx import model_statistics
from iotbx.pdb import extract_rfactors_resolutions_sigma
import iotbx.pdb.remark_3_interpretation
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
import iotbx.pdb
import iotbx.pdb.amino_acid_codes
import random
from cctbx import xray

def get_processed_pdb_file(pdb_file_name, cryst1, show_geometry_statistics):
  pdb_raw_records = smart_open.for_reading(
    file_name=pdb_file_name).read().splitlines()
  for rec in pdb_raw_records:
    if(rec.startswith("CRYST1 ")): pdb_raw_records.remove(rec)
  pdb_raw_records.append(cryst1)
  #############################################################################
  cif_objects = None
  if(show_geometry_statistics):
    t = time.ctime().split() # to make it safe to remove files
    time_stamp = "_"+t[4]+"_"+t[1].upper()+"_"+t[2]+"_"+t[3][:-3].replace(":","h")
    prefix = os.path.basename(pdb_file_name)
    from elbow.scripts import elbow_on_pdb_file
    from elbow.command_line import join_cif_files
    if len(sys.argv)>1: del sys.argv[1:]
    rc = elbow_on_pdb_file.run("\n".join(pdb_raw_records),
                               silent=True,
                               )
    cif_file = prefix+time_stamp+".cif"
    if rc is not None:
      hierarchy, cifs = rc
      if cifs:
        cif_lines = []
        for key in cifs:
          lines = cifs[key]
          if lines:
            cif_lines.append(lines)
        rc = join_cif_files.run(cif_lines, cif_file, no_file_access=True)
        if rc:
          f=file(cif_file, "wb")
          f.write(rc)
          f.close()

    cif_objects = None
    if(os.path.isfile(cif_file)):
      cif_objects = []
      cif_objects.append((cif_file,
        mmtbx.monomer_library.server.read_cif(file_name = cif_file)))
  #############################################################################
  pdb_ip = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  pdb_ip.clash_guard.nonbonded_distance_threshold = -1.0
  pdb_ip.clash_guard.max_number_of_distances_below_threshold = 100000000
  pdb_ip.clash_guard.max_fraction_of_distances_below_threshold = 1.0
  processed_pdb_files_srv = utils.process_pdb_file_srv(
    cif_objects               = cif_objects,
    pdb_interpretation_params = pdb_ip,
    log                       = StringIO())
  processed_pdb_file, pdb_inp = \
    processed_pdb_files_srv.process_pdb_files(raw_records = pdb_raw_records)
  return processed_pdb_file, pdb_raw_records, pdb_inp

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*1./fmodel.r_free_flags.data().size()
  print "  Data:"
  try: f_obs_labels = f_obs_labels[:f_obs_labels.index(",")]
  except ValueError: pass
  result = " \n    ".join([
    "data_label              : %s"%f_obs_labels,
    "high_resolution         : "+format_value("%-5.2f",info.d_min),
    "low_resolution          : "+format_value("%-6.2f",info.d_max),
    "completeness_in_range   : "+format_value("%-6.2f",info.completeness_in_range),
    "completeness(d_min-inf) : "+format_value("%-6.2f",info.completeness_d_min_inf),
    "completeness(6A-inf)    : "+format_value("%-6.2f",info.completeness_6_inf),
    "wilson_b                : "+format_value("%-6.1f",fmodel.wilson_b()),
    "number_of_reflections   : "+format_value("%-8d",  info.number_of_reflections),
    "test_set_size           : "+format_value("%-8.4f",flags_pc),
    "test_flag_value         : "+format_value("%-d",   test_flag_value),
    "number_of_Fobs_outliers : "+format_value("%-8d",  n_outl),
    "anomalous_flag          : "+format_value("%-6s",  fmodel.f_obs.anomalous_flag())])
  print "   ", result

def show_model_vs_data(fmodel):
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = fmodel.r_free_flags.data().count(True)*100./\
    fmodel.r_free_flags.data().size()
  if(flags_pc == 0): r_free = None
  else: r_free = fmodel.r_free()
  k_sol = format_value("%-5.2f",fmodel.k_sol())
  b_sol = format_value("%-7.2f",fmodel.b_sol())
  b_cart = " ".join([("%8.2f"%v).strip() for v in fmodel.b_cart()])
  print "  Model_vs_Data:"
  result = " \n    ".join([
    "r_work(re-computed)                : "+format_value("%-6.4f",fmodel.r_work()),
    "r_free(re-computed)                : "+format_value("%-6.4f",r_free),
    "bulk_solvent_(k_sol,b_sol)         : %s %s"%(k_sol,b_sol),
    "overall_anisotropic_scale_(b_cart) : "+format_value("%-s",b_cart)])
  print "   ", result

def detect_dummy_atom (file_name):
  """ TODO: method to look for a dummy atom placed in density, most likeley with Coot"""
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
  for model in pdb_obj.hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          for atom in ag.atoms():
              atom_names = set([" DM  ", " X  "])
              if (atom.name in atom_names):
                  print "Dummy found:"
                  print "        atom.name:  ", atom.name
                  print "        atom.xyz:  ", atom.xyz
                  print "        atom.occ:  ", atom.occ
                  print "        atom.b:    ", atom.b

  return "done"

def grow_density(fmodel, file_name, xray_structures, x_center, y_center, z_center, radius ):
  """ method to improve local density using dummy atoms """
  print "grow density method"
#  xray_structure = xray_structures[0]
#  orth = xray_structure.unit_cell().orthogonalize
#  pdb_inp = iotbx.pdb.input(file_name=file_name)
#  hierarchy = pdb_inp.construct_hierarchy()
#  #
#  f = open("tmp2.pdb","w")
#  cs = xray_structure.crystal_symmetry()
#  print >> f, iotbx.pdb.format_cryst1_record(crystal_symmetry = cs)
#  print >> f, iotbx.pdb.format_scale_records(unit_cell = cs.unit_cell())
#
#  for atom in hierarchy.atoms():
#     print '     ', atom.format_atom_record()
#     print >> f, atom.format_atom_record()

#  f.close()
#  below code left just now - sort of works but can work out how to get chain, res etc
#  for i, sc in enumerate(xray_structure.scatterers()):
#    print i, sc.site
#    a = iotbx.pdb.hierarchy.atom_with_labels()
#    label = sc.label.upper()
#    a.name = sc.scattering_type
#    a.resname = label[:3]
#    a.serial = i+1
#    a.resseq = i+1
#    a.xyz = orth(sc.site)
#    a.occ = sc.occupancy
#    a.b = adptbx.u_as_b(sc.u_iso)
#    a.element = sc.scattering_type
#    print >> f, a.format_atom_record_group()


#  rr = random.randrange
#  for k in range(1,100):
#    scatterer = xray.scatterer(
#      site = (rr(-k,k)/100., rr(-k,k)/100., rr(-k,k)/100.),
#      scattering_type = "C",
#      u = 0.2)
#    xray_structure.add_scatterer(scatterer)
#
#
#  for i, sc in enumerate(xray_structure.scatterers()):
#    a = iotbx.pdb.hierarchy.atom_with_labels()
#    a.serial = i+1
#    a.name = " C  "
#    a.resname = "DUM"
#    a.resseq = i+1
#    a.xyz = orth(sc.site)
#    a.occ = sc.occupancy
#    a.b = adptbx.u_as_b(sc.u_iso)
#    a.element = sc.scattering_type
#    print >> f, a.format_atom_record_group()



  """TODO: trying to work out how to add atom correctly to pdb, so can create models /
  need to make new atom group?  Currently adds atoms to last residue (as shown in junk.pdb) /
  which is obviously not correct.  Code below will be able to produce atom grids
  """

  x_start = float(x_center) - (float(radius)/2)
  x_end = float(x_center) + (float(radius)/2)

  y_start = float(y_center) - (float(radius)/2)
  y_end = float(y_center) + (float(radius)/2)

  z_start = float(z_center) - (float(radius)/2)
  z_end = float(z_center )+ (float(radius)/2)

  step_size = 1


  for overlap_start in [0.0, 0.2, 0.4, 0.6, 0.8]:
      """TODO: needs to be in x y and z for overlap to fill box?"""
      """TODO: Perhaps pull out into a different method so I can alternative atom fillers"""
       """TODO: e.g. random, aster, etc..."""
      f = open("tmp"+str(overlap_start)+".pdb","w")
      xray_structure = "" # trying to blank scatterers
      xray_structure = xray_structures[0]
      orth = xray_structure.unit_cell().orthogonalize
      pdb_inp = iotbx.pdb.input(file_name=file_name)
      hierarchy = pdb_inp.construct_hierarchy()

      cs = xray_structure.crystal_symmetry()
      print >> f, iotbx.pdb.format_cryst1_record(crystal_symmetry = cs)
      print >> f, iotbx.pdb.format_scale_records(unit_cell = cs.unit_cell())
      for atom in hierarchy.atoms(): print >> f, atom.format_atom_record()

      xray_structure_dummy_atoms = xray.structure(crystal_symmetry=cs)

      #
      x_coord = x_start + overlap_start
      """ TODO: need to think about which way overlap moves """
      y_coord = y_start + overlap_start
      z_coord = z_start + overlap_start

      while x_coord <= x_end:
          y_coord = y_start + overlap_start
          while y_coord <= y_end:
              z_coord = z_start + overlap_start
              while z_coord <= z_end:
                  print x_coord, y_coord, z_coord
                  scatterer = xray.scatterer(
                  site = (x_coord, y_coord, z_coord),
                  scattering_type = "N",
                  u = 0.2)
                  xray_structure_dummy_atoms.add_scatterer(scatterer)
                  #print scatterer.site



                  z_coord = z_coord + step_size
              y_coord = y_coord + step_size
          x_coord = x_coord + step_size

      for i, sc in enumerate(xray_structure_dummy_atoms.scatterers()):
         a = iotbx.pdb.hierarchy.atom_with_labels()
         a.serial = i+1
         a.name = sc.scattering_type
         a.resname = "DUM"
         a.resseq = i+1
         a.xyz = sc.site
         #print sc.site, a.xyz
         a.occ = sc.occupancy
         a.b = adptbx.u_as_b(sc.u_iso)
         a.element = sc.scattering_type
         print >> f, a.format_atom_record_group()
      f.close()
  print "   "



def run(args,
        command_name             = "mmtbx.grow_density",
        show_geometry_statistics = True):
  if(len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="%s reflection_file pdb_file [options]" % command_name,
    description='Example: %s data.mtz model.pdb --x_center=-2.76 --y_center=0.89 --z_center=9.71 --radius=10'%command_name)
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
    .option(None, "--twin_law",
      action="store",
      default=None,
      type="string",
      help="Provide twin law operator if twinned data used (for example: h,-h-k,-l).")
    .option(None, "--scattering_table",
      action="store",
      default="n_gaussian",
      type="string",
      help="Choice for scattering table: n_gaussian (default) or wk1995 or it1992 or neutron.")
     .option(None, "--x_center",
      action="store",
      default=None,
      type="string",
      help="X coordinate of poor density.")
     .option(None, "--y_center",
      action="store",
      default=None,
      type="string",
      help="Y coordinate of poor density.")
     .option(None, "--z_center",
      action="store",
      default=None,
      type="string",
      help="Z coordinate of poor density.")
     .option(None, "--radius",
      action="store",
      default=None,
      type="string",
      help="Radius of required box.")
    .option("--ignore_giant_models_and_datasets",
      action="store_true",
      help="Ignore too big models and data files to avoid potential memory problems.")
    ).process(args=args)
  if(command_line.options.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  if(command_line.options.x_center is None):
    raise Sorry("Need to specify x center.")
  if(command_line.options.y_center is None):
    raise Sorry("Need to specify y center.")
  if(command_line.options.z_center is None):
    raise Sorry("Need to specify z center.")
  if(command_line.options.radius is None):
    raise Sorry("Need to specify radius.")
  crystal_symmetry = None
  crystal_symmetry_data = None
  crystal_symmetry_model = None
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
          crystal_symmetry_model = crystal_symmetry_from_pdb.extract_from(
            file_name=arg)
        except RuntimeError, e:
          if(str(e) == "No CRYST1 record."): pass
      if(not arg_is_processed):
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          hkl_file_name = arg
          arg_is_processed = True
          crystal_symmetry_data = crystal_symmetry_from_any.extract_from(arg)
      if(not arg_is_processed):
        raise Sorry(
          "The command line argument %s is not a valid PDB or data file."%arg)
  if(hkl_file_name is None): raise Sorry("No file with Fobs is given.")
  if(pdb_file_name is None): raise Sorry("No PDB file is given.")
  print "Model and data files: %s %s"%(
    format_value("%5s",os.path.basename(pdb_file_name)),
    format_value("%5s",os.path.basename(hkl_file_name)))
  if([crystal_symmetry_model, crystal_symmetry_data].count(None)==0):
    if([crystal_symmetry_model.space_group_info(),
        crystal_symmetry_data.space_group_info()].count(None)==0):
      if(not crystal_symmetry_model.is_similar_symmetry(crystal_symmetry_data)):
        raise Sorry("Crystal symmetry mismatch between data and PDB files.")
  if(crystal_symmetry_model is not None):
    crystal_symmetry = crystal_symmetry_model
  if(crystal_symmetry_data is not None):
    crystal_symmetry = crystal_symmetry_data
  if(crystal_symmetry is None): raise Sorry("Crystal symmetry is not defined.")
  if(pdb_file_name is None): raise Sorry("No PDB file given.")
  if(reflection_file is None): raise Sorry("No reflection file given.")
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
  number_of_reflections = f_obs.indices().size()
  if(command_line.options.ignore_giant_models_and_datasets and
     number_of_reflections > data_size_max_reflections):
    raise Sorry("Too many reflections: %d"%number_of_reflections)
  #
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  #
  processed_pdb_file, pdb_raw_records, pdb_inp = get_processed_pdb_file(
    pdb_file_name = pdb_file_name,
    cryst1 = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry),
    show_geometry_statistics = False)
  #
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = command_line.options.scattering_table,
    d_min              = f_obs.d_min())
  xray_structures = xsfppf.xray_structures
  if(not xray_structures[0].crystal_symmetry().is_similar_symmetry(
     f_obs.crystal_symmetry())):
    raise Sorry("Inconsistent crystal symmetry.")
  #
  print "  Unit cell:       ", f_obs.unit_cell()
  print "  Space group:     ", f_obs.crystal_symmetry().space_group_info().\
    symbol_and_number(), "number of symmetry operators:",\
    f_obs.crystal_symmetry().space_group_info().type().group().order_z()
  print "  Unit cell volume: %-15.4f" % f_obs.unit_cell().volume()
  #
  #
  fmodel = utils.fmodel_simple(xray_structures = xray_structures,
                               f_obs           = f_obs,
                               r_free_flags    = r_free_flags,
                               twin_law        = command_line.options.twin_law)
  n_outl = f_obs.data().size() - fmodel.f_obs.data().size()
  show_data(fmodel          = fmodel,
            n_outl          = n_outl,
            test_flag_value = test_flag_value,
            f_obs_labels    = f_obs.info().label_string())
  show_model_vs_data(fmodel)
  #
  grow_density(fmodel, pdb_file_name, xray_structures,x_center=command_line.options.x_center,\
  y_center=command_line.options.y_center, z_center=command_line.options.z_center,\
  radius=command_line.options.radius )
