from __future__ import division
from string import ascii_uppercase

import iotbx
from mmtbx.chemical_components import get_bond_pairs

def add_side_chain_hydrogen(atom_group):
  def _get_side_chain_hydrogen_xyz(atom_names):
    from elbow.chemistry.xyzClass import xyzClass
    vec_oo = xyzClass(atom_names[0].xyz) - xyzClass(atom_names[2].xyz)
    mid = xyzClass(atom_names[0].xyz) - vec_oo/2
    vec_cm = mid - xyzClass(atom_names[1].xyz)
    h_xyz = xyzClass(atom_names[2].xyz) + vec_cm*1/abs(vec_cm)
    return h_xyz
  #
  heavy_name = None
  if atom_group.resname == "ASP":
    atom_names = ["OD1", "CG", "OD2", "HD2"]
  elif atom_group.resname == "GLU":
    atom_names = ["OE1", "CD", "OE2", "HE2"]
  else:
    raise Sorry("residue %s has no side chain acid group" % atom_group.resname)
  count=0
  for i, name in enumerate(atom_names):
    for atom in atom_group.atoms():
      if atom.name.strip()==name:
        atom_names[i]=atom
        count+=1
        break
  if count!=3: return None
  h_xyz = _get_side_chain_hydrogen_xyz(atom_names)
  new_atom = iotbx.pdb.hierarchy.atom()
  new_atom.element = "H"
  new_atom.name = atom_names[3]
  new_atom.xyz = h_xyz
  atom_group.insert_atom(-1, new_atom)
  return new_atom

def neutron_exchange_hydrogens(hierarchy,
                               cifs=None,
                               exchange_sites_only=True,
                               perdeuterate=False,
                               only_chain_id=None,
                               only_resseq=None,
                               side_chain_acids=False,
                               verbose=False,
                               ):
  if verbose:
    print("""neutron_exchange_hydrogens
    hierarchy : %s
    cifs : %s
    exchange_sites_only : %s
    perdeuterate : %s
    """ % (hierarchy,
           cifs,
           exchange_sites_only,
           perdeuterate,
           ))
  def _get_exchange_sites(atom_group, side_chain_acids=False, verbose=False):
    exchanges = []
    for atom in atom_group.atoms():
      # if atom.hetero: continue
      if atom.element.strip() == "":
        raise Sorry("Need element columns in input PDB")
      if atom.element.strip() in ["H"]:
        for a1, a2 in bonds1:
          other = None
          if atom.name.strip()==a1:
            other = a2
            break
          elif atom.name.strip()==a2:
            other = a1
            break
        else:
          for a1, a2 in bonds2:
            if atom.name.strip()==a1:
              other = a2
              break
            elif atom.name.strip()==a2:
              other = a1
              break
        if other:
          for heavy in atom_group.atoms():
            if heavy.element.strip() in ["H"]: continue
            if heavy.name.strip() == other.strip():
              break
          else:
            continue
          if verbose: print('heavy :%s: :%s:' % (heavy.quote(),
                                                 heavy.element.strip()))
          if heavy.element.strip() in ["C"]: continue
          exchanges.append(atom)
    if side_chain_acids:
      rc = add_side_chain_hydrogen(atom_group)
      if rc:
        exchanges.append(rc)
    return exchanges
  ##############################################
  def _get_hydrogens(atom_group, verbose=False):
    exchanges = []
    for atom in atom_group.atoms():
      # if atom.hetero: continue
      if atom.element.strip() == "":
        raise Sorry("Need element columns in input PDB")
      if not atom.element.strip() in ["H", "D"]: continue
      exchanges.append(atom)
    return exchanges
  ############################
  if verbose: hierarchy.show()
  exchange_count = 0
  for model in hierarchy.models():
    if verbose: print('model: "%s"' % model.id)
    for chain in model.chains():
      if only_chain_id is not None and chain.id!=only_chain_id: continue
      if verbose: print('chain: "%s"' % chain.id)
      for residue_group in chain.residue_groups():
        if ( only_resseq is not None and
             residue_group.resseq.strip()!=only_resseq
             ):
          continue
        if verbose: print('  residue_group: resseq="%s" icode="%s"' % (
          residue_group.resseq, residue_group.icode))
        altlocs = []
        for atom_group in residue_group.atom_groups():
          altlocs.append(atom_group.altloc)
        for atom_group in residue_group.atom_groups():
          if verbose: print('  atom_group: resname="%s" altloc="%s"' % (
            atom_group.resname, atom_group.altloc))
          bonds1 = get_bond_pairs(atom_group.resname)
          bonds2 = get_bond_pairs(atom_group.resname,
                                  alternate=True)
          if 0: # PVA: This triggers a bug. Added 13-MAR-2025, modules/phenix_regression/refinement/neutron/tst_ready_set_all_d.py
          #if cifs:
            bonds3 = get_bond_pairs_from_cif(cifs.get(atom_group.resname, None))
            if bonds3:
              bonds1 = bonds3
              bonds2 = bonds3
            else:
              if bonds1 is None: continue
              if bonds2 is None: continue
          else:
            if bonds1 is None: continue
            if bonds2 is None: continue
          if verbose: print('    atom_group: altloc="%s" resname="%s"' % (
            atom_group.altloc, atom_group.resname))

          # for deuteriums
          deuteriums=True
          for atom in atom_group.atoms():
            if atom.hetero: continue
            if atom.element.strip() in ["D"]:
              break
          else:
            deuteriums=False
          if deuteriums: break

          if exchange_sites_only:
            exchanges = _get_exchange_sites(atom_group,
                                            side_chain_acids=side_chain_acids,
                                            verbose=verbose,
            )
          else:
            exchanges = _get_hydrogens(atom_group, verbose=verbose)

          for u, l in enumerate(ascii_uppercase):
            if l not in altlocs:
              break
          else:
            assert 0

          if perdeuterate:
            for atom in atom_group.atoms():
              if atom.element.strip() not in ["H", "D"]: continue
              if atom in exchanges: continue
              atom.element = " D"
              atom.name = atom.name.replace("H","D", 1)

          if exchanges and not atom_group.altloc.strip():
            for i_exchange, exchange in enumerate(exchanges):
              if verbose:
                print('Exchanging',exchange.quote())
              exchange.occ=0.5
              atom_group.remove_atom(exchange)
              if i_exchange==0:
                new_atom_group = iotbx.pdb.hierarchy.atom_group()
                new_atom_group.altloc = l
                new_atom_group.resname = atom_group.resname
              new_atom_group.append_atom(exchange)
              exchange_count += 1
            residue_group.append_atom_group(new_atom_group)
            other = new_atom_group.detached_copy()
            other.altloc = ascii_uppercase[u+1]
            for exchange in other.atoms():
              exchange.element = " D"
              exchange.name = exchange.name.replace("H","D", 1)
            residue_group.append_atom_group(other)
  hierarchy.atoms().reset_serial()
  # hierarchy.overall_counts().show()
  if verbose: print("    Deuterium exchange count : %d" % exchange_count)
  return hierarchy
