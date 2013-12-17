
from __future__ import division
import os
from mmtbx.validation import restraints
from mmtbx.monomer_library import server, pdb_interpretation
from libtbx.utils import Usage
from cStringIO import StringIO

def get_bond_and_angle_outliers(
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      type=None):
  rc = restraints.combined(
         pdb_hierarchy=pdb_hierarchy,
         xray_structure=xray_structure,
         geometry_restraints_manager=geometry_restraints_manager,
         ignore_hd=True,
         outliers_only=False) #get them all
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
  if len(args) != 2:
    raise Usage(
      "mmtbx.mp_geo input.pdb out_file")
  file_name = args[0]
  out_file = args[1]
  log = StringIO()
  basename = os.path.basename(file_name)
  out = file(out_file, 'w')
  use_neutron_distances = False
  processed_pdb_file = pdb_interpretation.process(
    mon_lib_srv              = server.server(),
    ener_lib                 = server.ener_lib(),
    file_name                = file_name,
    strict_conflict_handling = True,
    use_neutron_distances    = use_neutron_distances,
    force_symmetry           = True,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log                      = log)
  grm = processed_pdb_file.geometry_restraints_manager()
  rc = get_bond_and_angle_outliers(
         pdb_hierarchy=processed_pdb_file.all_chain_proxies.pdb_hierarchy,
         xray_structure=processed_pdb_file.xray_structure(),
         geometry_restraints_manager=grm)

  #get chain types
  chain_types = {}
  for chain in processed_pdb_file.all_chain_proxies.\
                 pdb_hierarchy.models()[0].chains() :
    chain_id = chain.id
    main_conf = chain.conformers()[0]
    if chain_types.get(chain_id) not in ["NA", "PROTEIN"]:
      if (main_conf.is_na()) :
        chain_types[chain_id] = "NA"
      elif (main_conf.is_protein()):
        chain_types[chain_id] = "PROTEIN"
      else:
        chain_types[chain_id] = "UNK"
  outliers = []
  #bonds
  for result in rc.bonds.results:
    atom_info = result.atoms_info[0]
    # label:chain:number:ins:alt:type:measure:value:sigmas:class
    atoms_str = get_atoms_str(atoms_info=result.atoms_info)
    altloc = get_altloc(atoms_info=result.atoms_info)
    outliers.append( [atom_info.chain_id,
                      atom_info.resseq,
                      atom_info.icode,
                      altloc,
                      atom_info.resname,
                      atoms_str,
                      result.model,
                      result.score,
                      chain_types[atom_info.chain_id]] )
  #angles
  for result in rc.angles.results:
    #print dir(result)
    atom_info = result.atoms_info[0]
    # label:chain:number:ins:alt:type:measure:value:sigmas:class
    atoms_str = get_atoms_str(atoms_info=result.atoms_info)
    altloc = get_altloc(atoms_info=result.atoms_info)
    outliers.append( [atom_info.chain_id,
                      atom_info.resseq,
                      atom_info.icode,
                      altloc,
                      atom_info.resname,
                      atoms_str,
                      result.model,
                      result.score,
                      chain_types[atom_info.chain_id]] )
  #print output
  #print >> out, "#bonds:%d:%d" % (rc.bonds.n_outliers, rc.bonds.n_total)
  #print >> out, "#angles:%d:%d" % (rc.angles.n_outliers, rc.angles.n_total)
  for outlier in outliers:
    print >> out, "%s:%2s:%s:%s:%s:%s:%s:%.3f:%.3f:%s" % (
      basename, outlier[0], outlier[1], outlier[2], outlier[3],
      outlier[4], outlier[5], outlier[6], outlier[7], outlier[8])
  out.close()
