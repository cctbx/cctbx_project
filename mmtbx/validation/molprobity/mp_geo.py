
from __future__ import absolute_import, division, print_function
import os
from mmtbx.validation import restraints
from mmtbx.monomer_library import server, pdb_interpretation
import iotbx.phil
from six.moves import cStringIO as StringIO

def get_master_phil():
  return iotbx.phil.parse(
    input_string="""
  mp_geo {
    pdb = None
      .type = path
      .multiple = True
    out_file = None
      .type = path
    bonds_and_angles = False
      .type = bool
    kinemage = False
      .type = bool
    rna_backbone = False
      .type = bool
    outliers_only = False
      .type = bool
    cdl = False
      .type = bool
  }
  """,process_includes=True)

def get_bond_and_angle_outliers(
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      use_segids,
      outliers_only=False,
      type=None):
  rc = restraints.combined(
         pdb_hierarchy=pdb_hierarchy,
         xray_structure=xray_structure,
         geometry_restraints_manager=geometry_restraints_manager,
         ignore_hd=True,
         outliers_only=outliers_only,
         use_segids_in_place_of_chainids=use_segids)
  return rc

def get_atoms_str(atoms_info):
  return_str = ""
  connector = ""
  # bond
  if len(atoms_info) == 2:
    connector = '--'
  # angle
  elif len(atoms_info) == 3:
    connector = '-'
  for atom_info in atoms_info:
    return_str = return_str+atom_info.name.strip()+connector
  return return_str.strip('-')

def get_altloc(atoms_info):
  altloc = ' '
  for atom_info in atoms_info:
    if altloc == ' ':
      if atom_info.altloc != '':
        altloc = atom_info.altloc
  return altloc

def run(args):
  """
  I suggest adding here:
  cctbx_project/mmtbx/validation/regression/tst_mp_geo.py
  test cases with just .pdb, without arguments, etc.
  """
  master_phil = get_master_phil()
  import iotbx.phil
  input_objects = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    pdb_file_def="mp_geo.pdb")
  work_params = input_objects.work.extract()
  assert len(work_params.mp_geo.pdb) == 1, "Need a model file to run"
  file_name = work_params.mp_geo.pdb[0]
  out_file = None
  if work_params.mp_geo.out_file != None:
    out_file = work_params.mp_geo.out_file
  do_bonds_and_angles = work_params.mp_geo.bonds_and_angles
  do_kinemage = work_params.mp_geo.kinemage
  do_rna_backbone = work_params.mp_geo.rna_backbone
  outliers_only = work_params.mp_geo.outliers_only
  use_cdl = work_params.mp_geo.cdl
  log = StringIO()
  basename = os.path.basename(file_name)
  if out_file == None:
    import sys
    out = sys.stdout
  else:
    if do_bonds_and_angles:
      out = open(out_file, 'w')
    elif do_kinemage:
      out = open(out_file, 'a')
    elif do_rna_backbone:
      out = open(out_file, 'w')
  restraints_loading_flags = {}
  restraints_loading_flags["use_neutron_distances"]=False
  from mmtbx.validation import utils
  params = pdb_interpretation.master_params.extract()
  params.restraints_library.cdl = use_cdl
  params.clash_guard.nonbonded_distance_threshold = None
  params.allow_polymer_cross_special_position=True
  params.flip_symmetric_amino_acids=False
  processed_pdb_file = pdb_interpretation.process(
    params                   = params,
    mon_lib_srv              = server.server(),
    ener_lib                 = server.ener_lib(),
    file_name                = file_name,
    strict_conflict_handling = True,
    restraints_loading_flags = restraints_loading_flags,
    force_symmetry           = True,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log                      = log)
  grm = processed_pdb_file.geometry_restraints_manager()
  use_segids = utils.use_segids_in_place_of_chainids(
                 hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  if do_bonds_and_angles or do_kinemage:
    rc = get_bond_and_angle_outliers(
           pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
           xray_structure=processed_pdb_file.xray_structure(),
           geometry_restraints_manager=grm,
           use_segids=use_segids,
           outliers_only=outliers_only)
    #get chain types
    chain_types = {}
    for chain in processed_pdb_file.all_chain_proxies.\
                   pdb_hierarchy.models()[0].chains():
      if use_segids:
        chain_id = utils.get_segid_as_chainid(chain=chain)
      else:
        chain_id = chain.id
      main_conf = chain.conformers()[0]
      if chain_types.get(chain_id) not in ["NA", "PROTEIN"]:
        if (main_conf.is_na()):
          chain_types[chain_id] = "NA"
        elif (main_conf.is_protein()):
          chain_types[chain_id] = "PROTEIN"
        else:
          chain_types[chain_id] = "UNK"
    outliers = []
    #bonds
    #for result in rc.bonds.results:
    for result in sorted(rc.bonds.results, key=lambda x: (x.atoms_info[1].resseq, get_altloc(atoms_info=x.atoms_info), get_atoms_str(atoms_info=x.atoms_info))):
      atom_info = result.atoms_info[1]
      # label:chain:number:ins:alt:type:measure:value:sigmas:class
      atoms_str = get_atoms_str(atoms_info=result.atoms_info)
      altloc = get_altloc(atoms_info=result.atoms_info)
      chain_id = atom_info.chain_id
      outliers.append( [chain_id,
                        atom_info.resseq,
                        atom_info.icode,
                        altloc,
                        atom_info.resname,
                        atoms_str,
                        result.model,
                        result.score,
                        chain_types[atom_info.chain_id]] )
    #angles
    #for result in rc.angles.results:
    for result in sorted(rc.angles.results, key=lambda x: (x.atoms_info[1].resseq, get_altloc(atoms_info=x.atoms_info), get_atoms_str(atoms_info=x.atoms_info))):
      atom_info = result.atoms_info[1]
      # label:chain:number:ins:alt:type:measure:value:sigmas:class
      atoms_str = get_atoms_str(atoms_info=result.atoms_info)
      altloc = get_altloc(atoms_info=result.atoms_info)
      chain_id = atom_info.chain_id
      outliers.append( [chain_id,
                        atom_info.resseq,
                        atom_info.icode,
                        altloc,
                        atom_info.resname,
                        atoms_str,
                        result.model,
                        result.score,
                        chain_types[atom_info.chain_id]] )

    for result in sorted(rc.chiralities.results, key=lambda x: (x.atoms_info[0].resseq, get_altloc(atoms_info=x.atoms_info), get_atoms_str(atoms_info=x.atoms_info))):
      atom_info = result.atoms_info[0]
      # label:chain:number:ins:alt:type:measure:value:sigmas:class
      #atoms_str = get_atoms_str(atoms_info=result.atoms_info)
      atoms_str = result.atoms_info[0].name.strip() #report chiral center instead of atom list
      altloc = get_altloc(atoms_info=result.atoms_info)
      chain_id = atom_info.chain_id
      outliers.append( [chain_id,
                        atom_info.resseq,
                        atom_info.icode,
                        altloc,
                        atom_info.resname,
                        atoms_str,
                        result.model,
                        result.score,
                        chain_types[atom_info.chain_id]] )

    if do_bonds_and_angles:
      for outlier in outliers:
        print("%s:%2s:%s:%s:%s:%s:%s:%.3f:%.3f:%s" % (
          basename, outlier[0], outlier[1], outlier[2], outlier[3],
          outlier[4], outlier[5], outlier[6], outlier[7], outlier[8]), file=out)
    elif do_kinemage:
      print(rc.bonds.kinemage_header, file=out)
      for result in rc.bonds.results:
        print(result.as_kinemage(), file=out)
      print(rc.angles.kinemage_header, file=out)
      for result in rc.angles.results:
        print(result.as_kinemage(), file=out)
      print(rc.chiralities.kinemage_header, file=out)
      for result in rc.chiralities.results:
        print(result.as_kinemage(), file=out)
    out.close()
  elif do_rna_backbone:
    from mmtbx.validation import utils
    rna_bb = utils.get_rna_backbone_dihedrals(processed_pdb_file)
    print(rna_bb, file=out)
    if out_file is not None:
      out.close()
