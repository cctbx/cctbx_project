import os, sys
import time

from mmtbx.conformation_dependent_library.rdl_database import rdl_database

from mmtbx.validation import rotalyze

def generate_rotamer_data(pdb_hierarchy,
                          exclude_outliers=True,
                          ):
  r = rotalyze.rotalyze(pdb_hierarchy=pdb_hierarchy)
  for rot in r.results:
    if exclude_outliers and rot.outlier: continue    
    yield (rot.resname,
           rot.id_str(),
           rot.chain_id,
           rot.altloc,
           rot.rotamer_name)

def adjust_rotomer_restraints(pdb_hierarchy,
                              geometry_restraints_manager,
                              i_seqs_restraints=None,
                              log=None,
                              verbose=False,
                              ):
  if log is None: log = sys.stdout
  t0=time.time()
  print >> log, "  Rotamer Conformation Library"
  def _get_i_seqs_restraints(pdb_hierarchy):
    name_restraints = {}
    chains=[]
    residues=[]
    for item in generate_rotamer_data(pdb_hierarchy=pdb_hierarchy):
      resname, id_str, chain, altloc, rotomer = item
      if resname not in rdl_database: continue
      if rotomer not in rdl_database[resname]: continue
      restraints = rdl_database[resname][rotomer]
      defaults = rdl_database[resname]["default"]
      name_restraints[id_str] = (chain,
                                resname,
                                altloc,
                                restraints,
                                defaults,
        )
      chains.append(chain)
      residues.append(resname)
    #
    i_seqs_restraints = {}
    i_seqs_restraints_reversal = {}
    #
    for model in pdb_hierarchy.models():
      if verbose: print 'model: "%s"' % model.id
      for chain in model.chains():
        if chain.id.strip() not in chains: continue
        if verbose: print 'chain: "%s"' % chain.id
        for conformer in chain.conformers():
          if verbose: print '  conformer: altloc="%s" only_residue="%s"' % (
            conformer.altloc, conformer.only_residue)
          for residue in conformer.residues():
            if residue.resname not in residues: continue
            if verbose:
              if residue.resname not in ["HOH"]:
                print '    residue: resname="%s" resid="%s"' % (
                  residue.resname, residue.resid())
                for atom in residue.atoms():
                  if verbose: print '         atom: name="%s"' % (atom.name)
            for id_str, values in name_restraints.items():
              chain_id, resname, altloc, restraints, defaults = values
              if chain_id.strip()!=chain.id: continue
              if resname.strip() !=residue.resname.strip(): continue
              if altloc.strip()  !=conformer.altloc.strip(): continue
              for names, values in restraints.items():
                i_seqs = []
                for name in names:
                  for atom in residue.atoms():
                    if name.strip()==atom.name.strip():
                      i_seqs.append(atom.i_seq)
                      break
                i_seqs_restraints[tuple(i_seqs)] = values
                i_seqs_restraints_reversal[tuple(i_seqs)] = defaults[names]
                i_seqs.reverse()
                i_seqs_restraints[tuple(i_seqs)] = values
                i_seqs_restraints_reversal[tuple(i_seqs)] = defaults[names]
    return i_seqs_restraints, i_seqs_restraints_reversal
  #
  i_seqs_restraints_reversal = None
  if i_seqs_restraints is None:
    i_seqs_restraints, i_seqs_restraints_reversal = _get_i_seqs_restraints(
      pdb_hierarchy,
      )
  count = 0
  for angle_proxy in geometry_restraints_manager.angle_proxies:
    if angle_proxy.i_seqs in i_seqs_restraints:
      if verbose: print " i_seqs %-15s initial %12.3f %12.3f" % (
        angle_proxy.i_seqs,
        angle_proxy.angle_ideal,
        angle_proxy.weight,
        )
      angle_proxy.angle_ideal = i_seqs_restraints[angle_proxy.i_seqs][0]
      angle_proxy.weight = 1/i_seqs_restraints[angle_proxy.i_seqs][1]**2
      if verbose: print " i_seqs %-15s final   %12.3f %12.3f" % (
        angle_proxy.i_seqs,
        angle_proxy.angle_ideal,
        angle_proxy.weight,
        )
      count += 1
  print >> log, "    Number of angles RDL adjusted : %d" % count
  print >> log, "    Time to adjust                : %0.3f" % (time.time()-t0)
  return i_seqs_restraints, i_seqs_restraints_reversal
    
