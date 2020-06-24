from __future__ import absolute_import, division, print_function
import sys
import time

from libtbx.utils import Sorry
from mmtbx.conformation_dependent_library.rdl_database import \
  get_rdl_database
from mmtbx.validation import rotalyze

substitute_residue_lookup = {
  "MSE" : "MET",
  }

rdl_database = get_rdl_database()

def get_rotamer_data(atom_group,
                     sa,
                     rotamer_evaluator,
                     rotamer_id,
                     all_dict,
                     sites_cart,
                     ):
  resname=substitute_residue_lookup.get(atom_group.resname,
                                        atom_group.resname,
                                        )
  if resname not in rdl_database.keys(): return None
  if resname in ["PRO", "GLY"]: return None
  model_rot, m_chis, value = rotalyze.evaluate_rotamer(
    atom_group=atom_group,
    sidechain_angles=sa,
    rotamer_evaluator=rotamer_evaluator,
    rotamer_id=rotamer_id,
    all_dict=all_dict,
    sites_cart=sites_cart,
    )
  return model_rot, m_chis, value

def setup_restraints():
  return None

def generate_rotamer_data(pdb_hierarchy,
                          exclude_outliers=True,
                          data_version="8000",
                          ):
  assert data_version == "8000","data_version not recognized."
  r = rotalyze.rotalyze(pdb_hierarchy=pdb_hierarchy,
                        data_version=data_version)
  for rot in r.results:
    #if rot.outlier:
    #  print (rot.resname,
    #       rot.id_str(),
    #       rot.chain_id,
    #       rot.altloc,
    #       rot.rotamer_name)
    #  assert 0
    if exclude_outliers and rot.outlier: continue
    yield (rot.resname,
           rot.id_str(),
           rot.chain_id,
           rot.altloc,
           rot.rotamer_name)

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      rdl_proxies=None,
                      esd_factor=1.,
                      exclude_backbone=False,
                      assert_rotamer_found=False,
                      rdl_selection=None,
                      log=None,
                      verbose=False,
                      data_version="8000",
                      ):
  assert data_version == "8000","data_version not recognized."
  assert not exclude_backbone
  loud=True and False
  if loud:
    verbose=1
  from mmtbx.rotamer.sidechain_angles import SidechainAngles
  from mmtbx.rotamer import rotamer_eval
  from mmtbx.validation import rotalyze
  #
  def _set_or_reset_angle_restraints(geometry,
                                     lookup,
                                     use_rdl_weights=False,
                                     verbose=False,
                                     ):
    count = 0
    for angle_proxy in geometry.angle_proxies:
      if angle_proxy.i_seqs in lookup:
        if verbose: print(" i_seqs %-15s initial %12.3f %12.3f" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          ))
        assert angle_proxy.angle_ideal<181
        angle_proxy.angle_ideal = lookup[angle_proxy.i_seqs][0]
        if use_rdl_weights: angle_proxy.weight = esd_factor/lookup[angle_proxy.i_seqs][1]**2
        if verbose: print(" i_seqs %-15s final   %12.3f %12.3f" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          ))
        count += 1
    return count
  #
  def _set_or_reset_dihedral_restraints(geometry,
                                        lookup,
                                        use_rdl_weights=False,
                                        verbose=False,
                                        ):
    count = 0
    for angle_proxy in geometry.dihedral_proxies:
      if angle_proxy.i_seqs in lookup:
        if verbose: print(" i_seqs %-15s initial %12.3f %12.3f %s %d" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          angle_proxy.alt_angle_ideals,
          angle_proxy.periodicity,
          ))
        angle_proxy.angle_ideal = lookup[angle_proxy.i_seqs][0]
        if use_rdl_weights: angle_proxy.weight = esd_factor/lookup[angle_proxy.i_seqs][1]**2
        angle_proxy.alt_angle_ideals=None
        angle_proxy.periodicity = lookup[angle_proxy.i_seqs][2]
        if verbose: print(" i_seqs %-15s final   %12.3f %12.3f %s %d" % (
          angle_proxy.i_seqs,
          angle_proxy.angle_ideal,
          angle_proxy.weight,
          angle_proxy.alt_angle_ideals,
          angle_proxy.periodicity,
          ))
        count += 1
    return count
  #
  t0=time.time()
  sa = SidechainAngles(False)
  rotamer_id = rotamer_eval.RotamerID()
  rotamer_evaluator = rotamer_eval.RotamerEval(data_version=data_version)
  sites_cart = None
  if current_geometry:
    sites_cart = current_geometry.sites_cart()
  #
  if rdl_proxies is None:
    rdl_proxies = []
  i_seqs_restraints = {}
  i_seqs_restraints_reverse = {}
  #
  def _alt_loc_atom_generator(residue_group, atom_group):
    atoms = []
    for ag in residue_group.atom_groups():
      if ag.altloc.strip()=="" or ag.altloc.strip()==atom_group.altloc.strip():
        for atom in ag.atoms(): yield atom
  #
  for model in hierarchy.models():
    #if verbose: print 'model: "%s"' % model.id
    for chain in model.chains():
      #if verbose: print 'chain: "%s"' % chain.id
      for residue_group in chain.residue_groups():
        all_dict = rotalyze.construct_complete_sidechain(residue_group)
        for atom_group in residue_group.atom_groups():
          if rdl_selection is None: pass
          elif 'all' in rdl_selection: pass
          elif atom_group.resname in rdl_selection: pass
          else: continue
          rc = get_rotamer_data(atom_group,
                                sa,
                                rotamer_evaluator,
                                rotamer_id,
                                all_dict=all_dict,
                                sites_cart=sites_cart,
            )
          if rc is None: continue
          rotamer_name, chis, value = rc
          if verbose:
            chis_str = "["
            if chis:
              for chi in chis:
                chis_str += " %6.1f," % chi
              chis_str = chis_str[:-1]+']'
            try:
              print("  %s %s %s %-5s %-60s %0.1f" % (
                chain.id,
                atom_group.resname,
                residue_group.resseq,
                rotamer_name,
                chis_str,
                value,
              ), file=log)
            except TypeError as e:
              print("  %s %s %s %-5s %-60s %s" % (
                chain.id,
                atom_group.resname,
                residue_group.resseq,
                rotamer_name,
                chis_str,
                value,
              ), file=log)
          if loud: print('exclude_backbone',exclude_backbone)
          if rotamer_name in ["OUTLIER"]: continue
          resname_lookup = substitute_residue_lookup.get(atom_group.resname,
                                                         atom_group.resname,
                                                         )
          if rotamer_name not in rdl_database[resname_lookup]:
            if assert_rotamer_found:
              raise Sorry("rotamer %s not found in db" % rotamer_name)
            else:
              continue
          if loud: print('not outlier')
          if rotamer_name not in rdl_database[resname_lookup]:
            if loud: print("rotamer_name %s not in RDL db" % rotamer_name)
            continue
          if loud: print('rotamer_name %s found' % rotamer_name)
          restraints = rdl_database[resname_lookup][rotamer_name]
          defaults = rdl_database[resname_lookup]["default"]
          rdl_proxies.append(atom_group.id_str())
          for names, values in restraints.items():
            i_seqs = []
            if exclude_backbone and names[1] in ["N", "CA", "C"]: continue
            for name in names:
              for atom in _alt_loc_atom_generator(residue_group, atom_group):
                if name.strip()==atom.name.strip():
                  i_seqs.append(atom.i_seq)
                  break
            if len(i_seqs)!=len(names): continue
            i_seqs_restraints[tuple(i_seqs)] = values
            i_seqs.reverse()
            i_seqs_restraints[tuple(i_seqs)] = values
            if names in defaults:
              i_seqs_restraints_reverse[tuple(i_seqs)] = defaults[names]
              i_seqs.reverse()
              i_seqs_restraints_reverse[tuple(i_seqs)] = defaults[names]
  #
  if loud:
    for i, atom in enumerate(hierarchy.atoms()):
      print(i, atom.quote())
  count = _set_or_reset_dihedral_restraints(geometry,
                                            i_seqs_restraints_reverse,
                                            #verbose=loud,
                                            )
  count_d = _set_or_reset_dihedral_restraints(geometry,
                                            i_seqs_restraints,
                                            verbose=loud,
                                            )
  count = _set_or_reset_angle_restraints(geometry,
                                         i_seqs_restraints_reverse,
                                         #verbose=loud,
                                         )
  count_a = _set_or_reset_angle_restraints(geometry,
                                         i_seqs_restraints,
                                         verbose=loud,
                                         )
  #
  print("    Number of angles, dihedrals RDL adjusted : %d, %d" % (
    count_a,
    count_d,
    ), file=log)
  print("    Time to adjust                           : %0.1fs" % (
    time.time()-t0), file=log)
  return rdl_proxies

def run(pdb_filename):
  from iotbx import pdb
  print(pdb_filename)
  pdb_inp = pdb.input(pdb_filename)
  pdb_hieratchy = pdb_inp.construct_hierarchy()
  for residue_group in pdb_hieratchy.residue_groups():
    rotamer_name, chis, value = get_rotamer_data(residue_group,
                                                 sa,
                                                 rotamer_evaluator,
                                                 rotamer_id,
      )
    print(rotamer_name)

if __name__=="__main__":
  run(sys.argv[1])
