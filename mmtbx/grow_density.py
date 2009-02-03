from cctbx.array_family import flex
import mmtbx.f_model
from mmtbx import utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.pdb import crystal_symmetry_from_pdb
from cStringIO import StringIO
from cctbx.xray import ext
from mmtbx import utils
from scitbx.array_family import shared
from cctbx import maptbx
import iotbx.phil
from iotbx import crystal_symmetry_from_any
from cctbx import adptbx
from libtbx.utils import Sorry
import os, math, time
from cctbx import miller
from mmtbx import map_tools
from libtbx import adopt_init_args
from libtbx import Auto, group_args
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx import pdb
from libtbx.str_utils import format_value
from libtbx import smart_open
from cctbx import xray

master_params_str = """\
pdb_file_name = None
  .type = str
  .multiple = False
  .help = PDB file name.
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
high_resolution = None
  .type = float
low_resolution = None
  .type = float

x_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

y_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

z_center = 0.5
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

radius = 10
  .type = float
  .help = X coordinate center of problem density. Determine in Coot by placing a temporary atom in area.
  .expert_level = 2

step_size = 1
  .type = float
  .help = Grid step size.
  .expert_level = 2

overlap_interval = 0.2
  .type = float
  .help = Grid interval size.
  .expert_level = 2

scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
"""
master_params = iotbx.phil.parse(master_params_str, process_includes=False)

def grow_density(fmodel, file_name, xray_structures, x_center, y_center, z_center, radius, step_size=1, overlap_interval=0.2 ):
  """ method to improve local density using dummy atoms """

  """TODO: trying to work out how to add atom correctly to pdb, so can create models /
  need to make new atom group?  Currently adds atoms to last residue (as shown in junk.pdb) /
  which is obviously not correct.  Code below will be able to produce atom grids /
  """

  x_start = float(x_center) - (float(radius)/2)
  x_end = float(x_center) + (float(radius)/2)

  y_start = float(y_center) - (float(radius)/2)
  y_end = float(y_center) + (float(radius)/2)

  z_start = float(z_center) - (float(radius)/2)
  z_end = float(z_center )+ (float(radius)/2)

  step_size = step_size
  overlap_interval = overlap_interval
  overlap_list = [0.0]

  while (overlap_list[len(overlap_list)-1] < step_size - overlap_interval ):
      overlap_list.append(overlap_list[len(overlap_list)-1]+ overlap_interval)

  for overlap_start in overlap_list:
      """TODO: needs to be in x y and z for overlap to fill box?"""
      """TODO: Perhaps pull out into a different method so I can alternative atom fillers"""
      """TODO: e.g. random, asteroid, etc..."""
      """TODO: A lot of repeated code to write the pdb files - if I'm keeping this then /
      should refactor"""
      x = open("tmp_x"+str(overlap_start)+".pdb","w")
      y = open("tmp_y"+str(overlap_start)+".pdb","w")
      z = open("tmp_z"+str(overlap_start)+".pdb","w")
      xray_structure = "" # trying to blank scatterers
      xray_structure = xray_structures[0]
      orth = xray_structure.unit_cell().orthogonalize
      pdb_inp = iotbx.pdb.input(file_name=file_name)
      hierarchy = pdb_inp.construct_hierarchy()
      cs = xray_structure.crystal_symmetry()
      print >> x, iotbx.pdb.format_cryst1_record(crystal_symmetry = cs)
      print >> y, iotbx.pdb.format_cryst1_record(crystal_symmetry = cs)
      print >> z, iotbx.pdb.format_cryst1_record(crystal_symmetry = cs)
      print >> x, iotbx.pdb.format_scale_records(unit_cell = cs.unit_cell())
      print >> y, iotbx.pdb.format_scale_records(unit_cell = cs.unit_cell())
      print >> z, iotbx.pdb.format_scale_records(unit_cell = cs.unit_cell())
      for atom in hierarchy.atoms():
          print >> x, atom.format_atom_record()
          print >> y, atom.format_atom_record()
          print >> z, atom.format_atom_record()
      xray_structure_dummy_atoms_x = xray.structure(crystal_symmetry=cs)
      xray_structure_dummy_atoms_y = xray.structure(crystal_symmetry=cs)
      xray_structure_dummy_atoms_z = xray.structure(crystal_symmetry=cs)
      #
      x_coord = x_start + overlap_start
      y_coord = y_start #+ overlap_start
      z_coord = z_start #+ overlap_start
      """TODO:  Adding very simpe check to make a sphere - there is probably a lot better \
       and more efficient way of doing this, for instance look here:
       http://en.literateprograms.org/Special:Downloadcode/Generating_all_integer_lattice_points_(Python)
       This should be all be refactored - unnecessary loops.  Not time limiting step.
       """
      while x_coord <= x_end:
          y_coord = y_start + overlap_start
          while y_coord <= y_end:
              z_coord = z_start + overlap_start
              while z_coord <= z_end:
                  distance = ((float(x_center) - float(x_coord))**2 + (float(y_center) - float(y_coord))**2 + (float(z_center) - float(z_coord))**2 )**1/2
                  if(distance <= float(radius)): xray_structure_dummy_atoms_x.add_scatterer(xray.scatterer( site = (x_coord, y_coord, z_coord), scattering_type = "N",u = 0.2))
                  z_coord = z_coord + step_size
              y_coord = y_coord + step_size
          x_coord = x_coord + step_size
      #
      x_coord = x_start #+ overlap_start
      y_coord = y_start + overlap_start
      z_coord = z_start #+ overlap_start
      while y_coord <= y_end:
          z_coord = z_start + overlap_start
          while z_coord <= z_end:
              x_coord = x_start + overlap_start
              while x_coord <= x_end:
                  distance = ((float(x_center) - float(x_coord))**2 + (float(y_center) - float(y_coord))**2 + (float(z_center) - float(z_coord))**2 )**1/2
                  if(distance <= float(radius)): xray_structure_dummy_atoms_y.add_scatterer(xray.scatterer( site = (x_coord, y_coord, z_coord), scattering_type = "N",u = 0.2))
                  x_coord = x_coord + step_size
              z_coord = z_coord + step_size
          y_coord = y_coord + step_size
      #
      x_coord = x_start #+ overlap_start
      y_coord = y_start #+ overlap_start
      z_coord = z_start + overlap_start
      while z_coord <= z_end:
          x_coord = x_start + overlap_start
          while x_coord <= x_end:
              y_coord = y_start + overlap_start
              while y_coord <= y_end:
                  distance = ((float(x_center) - float(x_coord))**2 + (float(y_center) - float(y_coord))**2 + (float(z_center) - float(z_coord))**2 )**1/2
                  if(distance <= float(radius)): xray_structure_dummy_atoms_z.add_scatterer(xray.scatterer( site = (x_coord, y_coord, z_coord), scattering_type = "N",u = 0.2))
                  y_coord = y_coord + step_size
              x_coord = x_coord + step_size
          z_coord = z_coord + step_size
#
      for i, sc in enumerate(xray_structure_dummy_atoms_x.scatterers()):
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
         print "x:", a.xyz
         print >> x, a.format_atom_record_group()
#
      for i, sc in enumerate(xray_structure_dummy_atoms_y.scatterers()):
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
         print "y:", a.xyz
         print >> y, a.format_atom_record_group()
#
      for i, sc in enumerate(xray_structure_dummy_atoms_z.scatterers()):
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
         print "z:", a.xyz
         print >> z, a.format_atom_record_group()


      x.close()
      y.close()
      z.close()
  print "   "


def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

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

def cmd_run(args, command_name):
  msg = """\

Description:
Tool to compute local map correlation coefficient.

How to use:
1: Run this command: phenix.real_space_correlation;
2: Copy, save into a file and edit the parameters shown between the lines *** below;
3: Run the command with this parameters file:
   phenix.real_space_correlation parameters.txt
"""
  if(len(args) == 0):
    print msg
    print "*"*79
    master_params.show()
    print "*"*79
    return
  else :
    arg = args[-1]
    if(not os.path.isfile(arg)):
      raise Sorry("%s is not a file."%arg)
    parsed_params = None
    try: parsed_params = iotbx.phil.parse(file_name=arg)
    except KeyboardInterrupt: raise
    except RuntimeError, e:
      print e
    params, unused_definitions = master_params.fetch(
        sources = [parsed_params],
        track_unused_definitions = True)
    if(len(unused_definitions)):
      print "*"*79
      print "ERROR:",
      print "Unused parameter definitions:"
      for obj_loc in unused_definitions:
        print " ", str(obj_loc)
      print "*"*79
      raise Sorry("Fix parameters file and run again.")
      print
    run(params = params.extract())

def run(params, d_min_default=1.5, d_max_default=999.9) :
  if(params.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  if(params.x_center is None):
    raise Sorry("Need to specify x center.")
  if(params.y_center is None):
    raise Sorry("Need to specify y center.")
  if(params.z_center is None):
    raise Sorry("Need to specify z center.")
  if(params.radius is None):
    raise Sorry("Need to specify radius.")
  crystal_symmetry = None
  crystal_symmetry_data = None
  crystal_symmetry_model = None
  hkl_file_name = None
  pdb_file_name = None
  reflection_file = None
  files = [str(params.pdb_file_name), str(params.reflection_file_name)]

  print params.pdb_file_name, params.reflection_file_name
  for arg in files:
    arg_is_processed = False
    if(not os.path.isfile(arg)):
      raise Sorry(" %s is not a file."%arg)
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
    scattering_table   = params.scattering_table,
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
  twin_law = None
  fmodel = utils.fmodel_simple(xray_structures = xray_structures,
                               f_obs           = f_obs,
                               r_free_flags    = r_free_flags,
                               twin_law        = twin_law)
  n_outl = f_obs.data().size() - fmodel.f_obs.data().size()
  #
  grow_density(fmodel, pdb_file_name, xray_structures,x_center=params.x_center,\
  y_center=params.y_center, z_center=params.z_center, radius=params.radius, step_size=params.radius, overlap_interval=params.overlap_interval )
