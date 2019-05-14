from __future__ import absolute_import, division, print_function
import os,sys
from iotbx import pdb
from iotbx import reflection_file_reader
from iotbx import file_reader
from mmtbx.refinement.real_space import individual_sites
import mmtbx
import libtbx.phil.command_line

master_phil = libtbx.phil.parse("""
flip_base {
  pdb_file = None
    .type = path
    .help = '''input PDB file'''
  reflection_file = None
    .type = path
    .help = '''Reflection file'''
  out_pdb_file = None
    .type = str
    .help = '''input PDB file'''
  chain = None
    .type = str
    .help = '''Chain of the residue that is to be flipped'''
  alt_loc = None
    .type = str
    .help = '''Alternate location of the residue that is to be flipped'''
  res_num = None
    .type = int
    .help = '''Residue number of the residue that is to be flipped'''
  n_refine_cycles = 3
    .type = int
    .help = '''Number of real-space refinement cycles'''
  help = False
    .type = bool
    .help = '''Show help message'''
}
""", process_includes=True)

def usage(msg='', log=sys.stderr):
  s = '''
******************************************************************************
Usage :
  python.phenix flipbase.py xxxx.mtz yyyy.pdb chain=A res_num=1

Will flip base of chain A residue 1 of yyyy.pdb and do a real-space
refinement using xxxx.mtz.

Required :
  pdb_file           input PDB file
  reflection_file    Reflection file
  chain              Chain of the residue that is to be flipped
  res_num            Residue number of the residue that is to be flipped

Options :
  out_pdb_file       input PDB file
  alt_loc            Alternate location of the residue that is to be flipped
  n_refine_cycles    Number of real-space refinement cycles
  help               Show help message
******************************************************************************

'''
  if msg != '' :
    s = '*'*79 + '\n\n!!!!!  %s  !!!!!\n' % msg + s
  print(s);sys.exit()

base_rotation_axes = {
  "A" : ["C1'", "N9"],
  "G" : ["C1'", "N9"],
  "C" : ["C1'", "N1"],
  "T" : ["C1'", "N1"],
  "U" : ["C1'", "N1"],
}
base_rotatable_atoms = {
  "A" : ["N1", "C2", "H2", "N3", "C4", "C5", "C6", "N6", "H61", "H62", "N7",
         "C8", "H8"],
  "G" : ["N1", "H1", "C2", "N2", "H21", "H22", "N3", "C4", "C5", "C6", "O6",
         "N7", "C8", "H8"],
  "C" : ["C2", "O2", "N3", "C4", "N4", "H41", "H42", "C5", "H5", "C6", "H6"],
  "T" : ["C2", "O2", "N3", "H3", "C4", "O4", "C5", "C7", "H71", "H72", "H73",
         "C6", "H6"],
  "U" : ["C2", "O2", "N3", "H3", "C4", "O4", "C5", "H5", "C6", "H6"],
}

def flip_base(atom_group, angle=180):
  import scitbx.matrix
  axis_point_1 = axis_point_2 = None
  rotateable_atoms = []
  base_name = atom_group.resname.strip()
  if ("r" in base_name):
    base_name = base_name.replace("r")
  elif (base_name.startswith("D") and len(base_name) == 2):
    base_name = base_name[1]
  assert base_name in base_rotation_axes.keys(), base_name
  for atom in atom_group.atoms():
    atom_name = atom.name.strip()
    if (atom_name == base_rotation_axes[base_name][0]):
      axis_point_1 = atom.xyz
    elif (atom_name == base_rotation_axes[base_name][1]):
      axis_point_2 = atom.xyz
    elif (atom_name in base_rotatable_atoms[base_name]):
      rotateable_atoms.append(atom)
  if (None in [axis_point_1, axis_point_2]):
    raise RuntimeError("Missing atom(s) for rotateable axis.")
  elif (len(rotateable_atoms) == 0):
    raise RuntimeError("Missing nucleotide base.")
  for atom in rotateable_atoms :
    atom.xyz = scitbx.matrix.rotate_point_around_axis(
      axis_point_1=axis_point_1,
      axis_point_2=axis_point_2,
      point=atom.xyz,
      angle=angle,
      deg=True)

def get_target_map(reflection_file_name,  log=sys.stderr):
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    reflection_file_name).as_miller_arrays()
  ma = miller_arrays[0]
  fft_map = ma.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  print("\nUsing sigma scaled map.\n", file=log)
  target_map = fft_map.real_map_unpadded()
  return target_map

def flip_and_refine(pdb_hierarchy,
                    xray_structure,
                    target_map,
                    geometry_restraints_manager,
                    chain,
                    res_num,
                    alt_loc = None,
                    n_refine_cycles = 3,
                    log = sys.stdout):
  sites_cart = xray_structure.sites_cart()
  ero = False
  for ch in pdb_hierarchy.chains():
    if ch.id.strip() != chain : continue
    for rg in ch.residue_groups():
      if rg.resseq_as_int() != res_num : continue
      if rg.have_conformers() and not alt_loc :
        s = 'Specified residue has alternate conformations. Please specify '
        raise RuntimeError(s + 'alt_loc on the command line')
      for residue in rg.atom_groups():
        if alt_loc and alt_loc != residue.altloc.strip():
          continue
        flip_base(residue, angle=180)

        sites_cart.set_selected(residue.atoms().extract_i_seq(),
          residue.atoms().extract_xyz())
        xray_structure = xray_structure.replace_sites_cart(sites_cart)
        sele = residue.atoms().extract_i_seq()
        print('real-space refinement BEGIN'.center(79,'*'), file=log)
        for i in range(n_refine_cycles):
          print('real-space refinement cycle %i...' % (i + 1), file=log)
          ero = individual_sites.easy(
            map_data                    = target_map,
            xray_structure              = xray_structure,
            pdb_hierarchy               = pdb_hierarchy,
            geometry_restraints_manager = geometry_restraints_manager,
            selection                   = sele)
        print('real-space refinement FINISHED'.center(79,'*'), file=log)
  if not ero : raise RuntimeError('Specified residue not found')
  return ero.pdb_hierarchy

def run(args):

  # phil parsing----------------------------------------------------------
  interpreter = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  sources = []
  for arg in args:
    if os.path.isfile(arg): #Handles loose filenames
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb"):
        sources.append(interpreter.process(arg="pdb_file=\"%s\"" % arg))
      if (input_file.file_type == "hkl"):
        sources.append(interpreter.process(arg="reflection_file=\"%s\"" % arg))
      elif (input_file.file_type == "phil"):
        sources.append(input_file.file_object)
    else: #Handles arguments with xxx=yyy formatting
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  params = work_params.flip_base
  if work_params.flip_base.pdb_file == None :
    usage('PDB file not provided!')
  if work_params.flip_base.reflection_file == None :
    usage('Reflection file not provided!')
  if work_params.flip_base.chain == None :
    usage('chain not provided!')
  if work_params.flip_base.res_num == None :
    usage('res_num file not provided!')
  if work_params.flip_base.out_pdb_file == None :
    fn = work_params.flip_base.pdb_file.replace('.pdb','_baseflip.pdb')
    work_params.flip_base.out_pdb_file = fn
    #usage('out_pdb_file file not provided!')
  params = work_params.flip_base

  if params.help:
    usage()
    sys.exit()
  # end phil parsing ------------------------------------------------------

  pdb_file_name = params.pdb_file
  reflection_file_name = params.reflection_file
  log = sys.stdout
  print('\ngettinsg target_map...\n', file=log)
  target_map = get_target_map(reflection_file_name, log)
  ppf = mmtbx.utils.process_pdb_file_srv(log=False).process_pdb_files(
    [pdb_file_name])[0]
  grm = mmtbx.restraints.manager(
      geometry      = ppf.geometry_restraints_manager(show_energies = False),
      normalization = True)
  pdb_hierarchy  = ppf.all_chain_proxies.pdb_hierarchy
  pdb_hierarchy.atoms().reset_i_seq()
  xray_structure = ppf.xray_structure(show_summary = False)
  flip_hierarchy = flip_and_refine(pdb_hierarchy,
                  xray_structure,
                  target_map = target_map,
                  geometry_restraints_manager = grm,
                  chain = params.chain,
                  res_num = params.res_num,
                  alt_loc = params.alt_loc,
                  n_refine_cycles = params.n_refine_cycles,
                  log= log)

  flip_hierarchy.write_pdb_file(params.out_pdb_file)
  print('\nOut written to %s' % params.out_pdb_file, file=log)

if __name__ == "__main__":
  run(sys.argv[1:])
