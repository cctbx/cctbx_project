from cctbx.array_family import flex
import mmtbx.f_model
from cctbx.array_family import flex
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
from cctbx.array_family import flex
import time
from cctbx import adptbx
import mmtbx.tls.tools
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.refinement
import scitbx.lbfgs
import random


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

atom_gap = 1.2
  .type = float
  .help = Thje gap between the atoms in the grid.
  .expert_level = 2

overlap_interval = 0.4
  .type = float
  .help = Grid interval size - the gap between two grids.
  .expert_level = 2

atom_type = N
    .type = str
    .help = Atom to use to make grid.  Nitrogen is a good choice


scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations

cycles = 5
    .type = int
    .help = Number of refinement cycles

iterations = 20
    .type = int
    .help = Number of refinement iterations

bfac_dummy_atoms = 25
  .type = float
  .help = The initial B factor for the dummy atoms. Should be around protein average.
  .expert_level = 2

bfac_cutoff = 40.0
  .type = float
  .help = The value above which to remove atoms for the final model.
  .expert_level = 2

occ_cutoff = 1.0
  .type = float
  .help = The value below which to remove atoms for the final model.
  .expert_level = 2



"""
master_params = iotbx.phil.parse(master_params_str, process_includes=False)

def grow_density(f_obs, r_free_flags, scattering_table, file_name, xray_structures, x_center, y_center, z_center, radius,

                overlap_interval=0.2,
                step_size=1,
                atom_type = "N",
                number_of_iterations = 20,
                number_of_cycles     = 5,

                bfac_dummy_atoms = 25,
                bfac_cutoff = 60,
                occ_cutoff = 1.0,

                sf_algorithm         = "fft", # XXX use "fft" in prectice
                free_r_fraction      = 0.05,
                target_name          = "ls_wunit_kunit",

                 ):
    """ method to improve local density using dummy atoms """
    """TODO: trying to work out how to add atom correctly to pdb, so can create models /
    need to make new atom group?  Currently adds atoms to last residue (as shown in junk.pdb) /
    which is obviously not correct.  Code below will be able to produce atom grids /
    """
    step = overlap_interval
    atom_gap = step_size

    kept_atoms = []
    x_start = float(x_center) - (float(radius))
    x_end = float(x_center) + (float(radius))
    y_start = float(y_center) - (float(radius))
    y_end = float(y_center) + (float(radius))
    z_start = float(z_center) - (float(radius))
    z_end = float(z_center )+ (float(radius))
    step_size = step_size
    overlap_interval = overlap_interval
    overlap_list = [0.0]

    kept_atoms = []

    center = [x_center,y_center,z_center]
    xray_structure = xray_structures[0]
    orth = xray_structure.unit_cell().orthogonalize
    pdb_inp = iotbx.pdb.input(file_name=file_name)
    hierarchy = pdb_inp.construct_hierarchy()
    cs = xray_structure.crystal_symmetry()
    atom_count = count_atoms(file_name)
    """still not sure this is right for making grid of grids"""
    x_range = frange(center[0],center[0]+atom_gap,step)
    y_range = frange(center[1],center[1]+atom_gap,step)
    z_range = frange(center[2],center[2]+atom_gap,step)
    all_atoms = []
    number_of_grids = len(x_range) * len(y_range) * len(z_range)
    print "Creating %s grids with atom spacing %s, each grid is %s apart" %(str(number_of_grids),str(atom_gap), str(step))
    print "Criteria for removing final atoms is b factor above %s, or occupancy below %s" %(str(bfac_cutoff), str(occ_cutoff))

    for x_start in x_range:
        for y_start in y_range:
            for z_start in z_range:
                new_center = [x_start, y_start, z_start ]
                print "Number of grids left: ", number_of_grids
                number_of_grids = number_of_grids - 1
                atom_grid = make_grid(new_center, radius, atom_gap,1.0,bfac_dummy_atoms)
                add_dummy_atoms = add_atoms_with_occ_b(file_name, atom_type, atom_grid, "tmp.pdb")

                processed_pdb_file, pdb_raw_records, pdb_inp = get_processed_pdb_file(
                pdb_file_name = "tmp.pdb",
                cryst1 = pdb.format_cryst1_record(crystal_symmetry = cs),
                show_geometry_statistics = False)

                twin_law = None
                xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
                processed_pdb_file = processed_pdb_file,
                scattering_table   = scattering_table,
                d_min              = f_obs.d_min())
                xray_structures = xsfppf.xray_structures
                fmodel = utils.fmodel_simple(xray_structures = xray_structures,
                               f_obs           = f_obs,
                               r_free_flags    = r_free_flags)
                new_model, new_r_factor, new_rfree = refine_atoms(fmodel, number_of_iterations, number_of_cycles )
                orth = new_model.unit_cell().orthogonalize
                #
                current_atom = 0
                for i_scatterer, sc in enumerate(new_model.scatterers()):
                    if (sc.occupancy > occ_cutoff and adptbx.u_as_b(sc.u_iso)< bfac_cutoff) and current_atom > atom_count :
                        a = iotbx.pdb.hierarchy.atom_with_labels()
                        a.xyz = orth(sc.site)
                        #dummy_atom = "HETATM %s  N   HET X%s    %8.3f%8.3f%8.3f%6.2f%6.2f           N"%(i_scatterer+1000,i_scatterer+1000,a.xyz[0] ,a.xyz[1] ,a.xyz[2] , sc.occupancy, adptbx.u_as_b(sc.u_iso))
                        dummy_atom_tuple = [a.xyz[0] ,a.xyz[1] ,a.xyz[2] , sc.occupancy, adptbx.u_as_b(sc.u_iso)]
                        kept_atoms.append(dummy_atom_tuple)
                    current_atom = current_atom + 1


    add_dummy_atoms = add_atoms_with_occ_b(file_name, atom_type, kept_atoms, "final.pdb")
    print "Finished"


def refine_atoms(fmodel, number_of_iterations, number_of_cycles):
  fmodels = mmtbx.fmodels(fmodel_xray = fmodel)
  size = fmodel.xray_structure.scatterers().size()
  selection = flex.size_t(xrange(size))
  lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
    max_iterations = number_of_iterations)
  # XXX Can we refine occupancy and B-factor sumultaneously ?
  for i in xrange(number_of_cycles):
    if 1: # XXX refine q and B separately
      fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = selection)
      minimized = mmtbx.refinement.minimization.lbfgs(
        fmodels                  = fmodels,
        lbfgs_termination_params = lbfgs_termination_params,
        collect_monitor          = False)
      print "Refined B-factors     Rwork = %8.6f Rfree = %8.6f"%(fmodel.r_work(),
        fmodel.r_free())
      fmodel.xray_structure.scatterers().flags_set_grad_occupancy(
        iselection = selection)
      minimized = mmtbx.refinement.minimization.lbfgs(
        fmodels                  = fmodels,
        lbfgs_termination_params = lbfgs_termination_params,
        collect_monitor          = False)
      print "Refined occupancies   Rwork = %8.6f Rfree = %8.6f"%(fmodel.r_work(),
        fmodel.r_free())
    if 0: # XXX refine q and B sumultaneously... This might be unstable....
      for sc in fmodel.xray_structure.scatterers():
        fl = sc.flags
        fl.set_grad_u_iso(True)
        fl.set_grad_occupancy(True)
      minimized = mmtbx.refinement.minimization.lbfgs(
        fmodels                  = fmodels,
        lbfgs_termination_params = lbfgs_termination_params,
        collect_monitor          = False)
      print "Refined q and B       Rwork = %8.6f Rfree = %8.6f"%(fmodel.r_work(),
        fmodel.r_free())
  assert minimized.xray_structure is fmodels.fmodel_xray().xray_structure
  assert minimized.xray_structure is fmodel.xray_structure
  structure = minimized.xray_structure
  orth = structure.unit_cell().orthogonalize
  return structure, fmodel.r_work(), fmodel.r_free()

#########################################
#########################################

def count_atoms(input_pdb_file_name):
    """A very simple class to count atoms - probably a better was in pdb toolbox or similiar"""
    count = 0
    try:
        input_file = open(input_pdb_file_name)
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        print "could not open file - are you sure directory and filename is correct:", input_pdb_file_name
        exit()
    #
    for line in input_file.readlines():
        #if (line.find("TER") or line.find("END"))== -1 :
        if ("ATOM" in line or "HETATOM" in line or "HEATOM" in line  ) :
            count = count +1
    return count

def add_atoms_with_occ_b(input_pdb_file_name, final_atom_type, dummy_atoms, output_pdb_file_name):
    """A class menthod to add the atoms - occ and b in array"""
    try:
        input_file = open(input_pdb_file_name)
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        print "could not open file - are you sure directory and filename is correct:", input_pdb_file_name
        exit()
    try:
        output_file = open(output_pdb_file_name,"w")
    except IOError, (errno, strerror):
        print "I/O error(%s): %s" % (errno, strerror)
        print "could not output file - are the permissions for this folder correct?"
        exit()
    #
    final_line = ""
    for line in input_file.readlines():
        #if (line.find("TER") or line.find("END"))== -1 :
        if not ("TER" in line or "END" in line ) :
            output_file.write(line)
            final_line = line
    #
    final_atom_id  = int(final_line[7:11].replace(' ', '')) + 1
    final_residue_id  = int(final_line[22:26].replace(' ', '')) + 1
    #
    for dummy_coords in dummy_atoms:
        #dummy_atom = "HETATM %4i  %s   HET X%4i    %8.3f%8.3f%8.3f%6.2f%6.2f           %s"\ #HOH doesn't make ugly lines...
        dummy_atom = "HETATM%5i  %s   HOH X%5i   %8.3f%8.3f%8.3f%6.2f%6.2f           %s"\
        %(final_atom_id,final_atom_type, final_residue_id, float(dummy_coords[0]), \
        float(dummy_coords[1]),float(dummy_coords[2]), float(dummy_coords[3]), float(dummy_coords[4]), final_atom_type )
        output_file.write(dummy_atom+"\n")
        final_atom_id  = final_atom_id + 1
        final_residue_id  = final_residue_id + 1
    #
    output_file.close()
    input_file.close()
    print "New file created"
    return " "


def make_grid(center, radius, gap,occ,bfac):
    """Make the grid of atoms"""
    print "Making grid"
    atom_list = []
    x_start = float(center[0]) - (float(radius))
    x_end = float(center[0]) + (float(radius))
    y_start = float(center[1]) - (float(radius))
    y_end = float(center[1]) + (float(radius))
    z_start = float(center[2]) - (float(radius))
    z_end = float(center[2])+ (float(radius))
    #
    x_range = frange(x_start,x_end,gap)
    y_range = frange(y_start,y_end,gap)
    z_range = frange(z_start,z_end,gap)
    #
    for x in x_range:
        for y in y_range:
            for z in z_range:
                atom_list.append([x,y,z,occ,bfac])
    print "Finished making grid"
    return atom_list


def frange(start, end=None, inc=None):
    """A range function, that does accept float increments..."""
    if end == None:
        end = start + 0.0
        start = 0.0
    if inc == None:
        inc = 1.0
    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
    return L



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
   phenix.grow_density parameters.txt
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
  if(params.atom_gap is None):
    raise Sorry("Need to specify atom_gap.")
  if(params.overlap_interval is None):
    raise Sorry("Need to specify overlap_interval.")

  if(params.atom_type is None):
    raise Sorry("Need to specify atom_type.")
  if(params.cycles is None):
    raise Sorry("Need to specify cycles.")
  if(params.iterations is None):
    raise Sorry("Need to specify iterations.")
  if(params.bfac_dummy_atoms is None):
    raise Sorry("Need to specify bfac_dummy_atoms.")
  if(params.bfac_cutoff is None):
    raise Sorry("Need to specify bfac_cutoff.")
  if(params.occ_cutoff is None):
    raise Sorry("Need to specify iterations.")

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
  scattering_table   = params.scattering_table
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = scattering_table,
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
                               r_free_flags    = r_free_flags)
  n_outl = f_obs.data().size() - fmodel.f_obs.data().size()
  #
  grow_density(f_obs, r_free_flags, scattering_table, pdb_file_name, xray_structures,x_center=params.x_center,\
  y_center=params.y_center, z_center=params.z_center, radius=params.radius, step_size=params.atom_gap, overlap_interval=params.overlap_interval, \
  atom_type= params.atom_type, number_of_cycles = params.cycles, number_of_iterations = params.iterations, \
  bfac_dummy_atoms = params.bfac_dummy_atoms, bfac_cutoff = params.bfac_cutoff , occ_cutoff = params.occ_cutoff  )


