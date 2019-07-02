from __future__ import absolute_import, division, print_function
# (jEdit options) :folding=explicit:collapseFolds=1:
#
#cablam_fingerprints
#Author: Christopher Williams, contact christopher.j.williams@duke.edu
#
#cablam_fingerprints holds objects and methods related to the identification of
#  secondary structure and other motifs in protein structures.
#It is called primarily by cablam_training, but may become called elsewhere as
#  motifs move into validation and correction routines
#
#2013-09-17: Initial Upload
#2014-02-07:
#Found motifs now stored as objects rather than as labels in the linked_residue
#  objects.  This allows more robust and varied handling downstream (eg output),
#  but required substantial rewrite of the functions in this file.

import os, sys
import libtbx.load_env
from libtbx import easy_pickle
from libtbx.utils import Sorry

#Checks return either "False" or a debug string describing the reason for failure.
#  This will let me get debug with maybe-cleaner code

#{{{ class objects
#-------------------------------------------------------------------------------
class motif(object):
  def __init__(self,
    motif_name="",residue_names={},superpose_order = {}):
    self.motif_name = motif_name
    self.residue_names = residue_names #labels for printing keyed by indices for residues in the motif
    self.superpose_order = superpose_order # format as {'a':['CA','O'],'b':['CA']}
    self.residues = []

  def add_residue(self, allowed_resname=[], banned_resname=[],cis_or_trans=None,
      sequence_move=None, bond_move='',
      end_of_motif=False, index=''):
    new_residue = residue(allowed_resname, banned_resname, cis_or_trans,
      sequence_move, bond_move, end_of_motif, index)
    self.residues.append(new_residue)
    return new_residue

#-------------------------------------------------------------------------------
class residue(object):
  def check_resname(self,residue):
    if residue.resname in reslist:
      return True

  def add_bond(self, required=True,banned=False,allow_bifurcated=False,
      src_atom='', trg_index=''):
    new_bond = bond(required, banned, allow_bifurcated, src_atom, trg_index)
    self.bonds.append(new_bond)
    return new_bond

  def __init__(self,
    allowed_resname=[], banned_resname=[], cis_or_trans=None,
    sequence_move=None, bond_move='',
    end_of_motif=False, index=''):

    self.allowed_resname=allowed_resname
    self.banned_resname=banned_resname
    self.cis_or_trans=cis_or_trans
    self.sequence_move=sequence_move
    self.bond_move=bond_move
    self.end_of_motif=end_of_motif
    self.index=index
    self.bonds=[]

#-------------------------------------------------------------------------------
class bond(object):
  def add_target_atom(self, atomname=None, anyatom=False, seqdist=None,
    anyseqdist=False):
    self.trg_atoms.append(target_atom(atomname, anyatom, seqdist, anyseqdist))

  def __init__(self,
    required=True, banned=False, allow_bifurcated=False,
    src_atom='', trg_index=''):

    self.required=required
    self.banned=banned
    self.allow_bifurcated=allow_bifurcated
    self.src_atom=src_atom
    self.trg_atoms=[] #what actually goes here? list of target_atom objects
    self.trg_index=trg_index

#-------------------------------------------------------------------------------
class target_atom(object):
  def __init__(self,
    atomname = None, anyatom=False, seqdist=None, anyseqdist=False):
    self.atomname = atomname
    self.anyatom = anyatom
    self.seqdist = seqdist
    self.anyseqdist = anyseqdist
#}}}

#{{{ found motif object
class motif_instance(object):
  def is_complete(self):
    if len(self.residues) == self.needed_length: return True
    else: return False

  #a check before printing that all residues have all values needed to print
  def has_all_measures(self, kinorder):
    for residue in self.residues.values():
      for needed_measure in kinorder:
        if needed_measure not in residue.measures:
          return False
    return True

  def add_to_instance(self, index, residue):
    self.residues[index] = residue

  def build_superposition(self, motif):
    superpose_list = []
    for index in motif.superpose_order:
      atoms = motif.superpose_order[index]
      residue = self.residues[index]
      chain = residue.chain
      resseq = str(residue.resnum)
      for atom in atoms:
        #altloc " " breaks down on alternates, probably do care about this.
        if len(residue.alts) > 1 and residue.firstalt(atom):
          atom_order = "(chain "+chain+" and resseq "+resseq+" and name "+atom+ " and altloc "+residue.firstalt(atom)+")"
        else:
          atom_order = "(chain "+chain+" and resseq "+resseq+" and name "+atom+")"
        superpose_list.append(atom_order)
    self.superpose_thus = "\"" + " or ".join(superpose_list) + "\""

  def __init__(self, fingerprint):
    self.needed_length = len(fingerprint.residues)
    self.residues = {}
    self.names = fingerprint.residue_names
    self.superpose_thus = ""
#}}}

#{{{ check protein
#-------------------------------------------------------------------------------
def check_protein(protein, motif_name_list):
  #debug = True
  debug = False
  found_motifs = {}
  motif_list = fetch_fingerprints(motif_name_list)
  # TODO: is protein a dict, if so put a hint comment
  reslist = list(protein.keys())
  reslist.sort()
  for motif in motif_list:
    found_motifs[motif.motif_name] = []
  #for resid in protein:
  for resid in reslist:
    residue = protein[resid]
    for motif in motif_list:
      candidate = motif_instance(motif)
      failed = check_for_motif(motif,residue,candidate)#returns either False or debug info
      if failed:
        candidate = None
        if debug==True: print(failed)
        continue
      else:
        #print candidate.needed_length, candidate.residues, len(candidate.residues)
        if candidate.is_complete():
          #print candidate.needed_length, candidate.residues, len(candidate.residues)
#1b16FH_A.pdb
#finished motif wide_helix_turn
#7 {'': <class 'mmtbx.cablam.cablam_fingerprints.residue'>} 1
          candidate.build_superposition(motif)
          found_motifs[motif.motif_name].append(candidate)
  return found_motifs
#-------------------------------------------------------------------------------
#}}}

#{{{ check for motif function
#-------------------------------------------------------------------------------
def check_for_motif(motif, active_residue, candidate):
  motif_index = 0
  motif_complete = False
  while not motif_complete:
    motif_residue = motif.residues[motif_index]
    if not active_residue:
      return "Missing residue for index " + motif_residue.index
    fail_residue_checks = fail_residue_check(motif_residue, active_residue)
    if fail_residue_checks:
      return fail_residue_checks
    for bond in motif_residue.bonds:
      fail_forbidden_bond_checks = fail_forbidden_bond_check(bond, active_residue)
      if fail_forbidden_bond_checks:
        return fail_forbidden_bond_checks
      fail_required_bond_checks = fail_required_bond_check(bond, active_residue, candidate)
      if fail_required_bond_checks:
        return fail_required_bond_checks
    fail_index_checks = fail_index_match(motif_residue.index, active_residue, candidate)
    if fail_index_checks:
      return fail_index_checks
    candidate.add_to_instance(motif_residue.index, active_residue)
    if motif_residue.end_of_motif:
      motif_complete = True
    else:
      motif_index += 1
      active_residue = do_move(motif_residue, active_residue, candidate)
  sys.stderr.write("finished motif "+ motif.motif_name+"\n")
  return 0
#-------------------------------------------------------------------------------
#}}}

#{{{ do move (keep)
def do_move(check_residue, active_residue, motif_in_progress):
  new_res = active_residue
  if check_residue.sequence_move:
    move_dist = check_residue.sequence_move
    if move_dist > 0:
      while move_dist > 0:
        if not new_res.nextres:
          return None
        else:
          new_res = new_res.nextres
          move_dist -= 1
    elif move_dist < 0:
      while move_dist < 0:
        if not new_res.prevres:
          return None
        else:
          new_res = new_res.prevres
          move_dist += 1
    return new_res
  elif check_residue.bond_move:
    move_index = check_residue.bond_move
    try:
      return motif_in_progress.residues[move_index]
    except KeyError:
      sys.stderr.write('\nTried to move to a residue not yet indexed at index \"'+check_residue.bond_move+'\"\n')
      sys.stderr.write('Please check fingerprint definition. Exiting . . .\n')
      sys.exit()
  else:
    sys.stderr.write('\nNo move specified for residue indexed as\"'+check_residue.index+'\"\n')
    sys.stderr.write('Please check fingerprint definition. Exiting . . .\n')
    sys.exit()
#}}}

#{{{ fail residue check
#-------------------------------------------------------------------------------
def fail_residue_check(motif_residue, active_residue):
  resname = active_residue.id_with_resname()[0:3]
  if motif_residue.allowed_resname:
    if resname.upper() not in motif_residue.allowed_resname:
      return active_residue.id_with_resname() + " resname not in allowed list"
  if motif_residue.banned_resname:
    if resname.upper() in motif_residue.banned_resname:
      return active_residue.id_with_resname() + " resname in banned list"
  if motif_residue.cis_or_trans:
    fail_cis_or_trans_check = fail_cis_or_trans_check(motif_residue.cis_or_trans, active_residue)
    if fail_cis_or_trans_check:
      return fail_cis_or_trans_check
  return 0

###minor function for fail residue check
def fail_cis_or_trans_check(cis_or_trans, active_residue):
  if 'omega' not in active_residue.measures:
    return active_residue.id_with_resname() + " has no value for omega dihedral"
  omega = active_residue.measures['omega']
  if cis_or_trans=='trans' and not (omega >= 120 or omega <= -120):
    return active_residue.id_with_resname() + " omega dihedral is not trans"
  if cis_or_trans=='cis' and not (omega >= -60 and omega <= 60):
    return active_residue.id_with_resname() + " omega dihedral is not cis"
  return 0
#-------------------------------------------------------------------------------
#}}}

#{{{ fail required bond check function
#-------------------------------------------------------------------------------
def fail_required_bond_check(check_bond, active_residue, candidate):
  if check_bond.banned or not check_bond.required:
    return 0
  if not active_residue.probe:
    return active_residue.id_with_resname() + " Missing all probe information"
  if check_bond.src_atom not in active_residue.probe:
    return active_residue.id_with_resname() + " Missing probe information for src_atom "+check_bond.src_atom
  if len(active_residue.probe[check_bond.src_atom]) > 1 and not check_bond.allow_bifurcated:
    return active_residue.id_with_resname() + " Disallowed bifucated bond at src_atom "+check_bond.src_atom
  #-----------------------------------------------------------------------------
  for src_atom in active_residue.probe:
    if src_atom == check_bond.src_atom:
      for bond_name in active_residue.probe[src_atom]:#ie for each bond partner
        observed_bond = active_residue.probe[src_atom][bond_name]
        for trg_atom in check_bond.trg_atoms:
          if (trg_atom.anyatom or trg_atom.atomname == observed_bond.atom) and (trg_atom.anyseqdist or trg_atom.seqdist == observed_bond.seqdist):
            index = check_bond.trg_index
            trg_residue = observed_bond.residue
            fail_index_matching = fail_index_match(index, trg_residue, candidate)
            if fail_index_matching:
              return fail_index_matching
            else:
              if index: candidate.add_to_instance(index, trg_residue)
              return 0
  return active_residue.id_with_resname() + " No matching bond for " + check_bond.src_atom
#-------------------------------------------------------------------------------
#}}}

#{{{ fail forbidden bond check
#-------------------------------------------------------------------------------
#At the moment, does not check for generic bans against indexed residues w/o
#  specific sequence relationship
def fail_forbidden_bond_check(check_bond, active_residue):
  if not check_bond.banned:
    return 0
  # TODO: is active residie.prob a dict ? if so put a hint comment
  if check_bond.src_atom not in active_residue.probe.keys():
    return 0
  src_atom = check_bond.src_atom
  for bond_name in active_residue.probe[src_atom]:
    observed_bond = active_residue.probe[src_atom][bond_name]
    for trg_atom in check_bond.trg_atoms:
      if ((trg_atom.anyatom)or(trg_atom.atomname == observed_bond.atom)) and ((trg_atom.anyseqdist)or(trg_atom.seqdist == observed_bond.seqdist)):
        return active_residue.id_with_resname() + " Residue has forbidden bond"
      else: pass
    else: pass
  else: return 0
#-------------------------------------------------------------------------------
#}}}

#{{{ bonding info storage reference
#-------------------------------------------------------------------------------
#Reference for how bond info gets stored into residue.probe
#    recordkey = trgResidue.id_with_resname() + trgAtomname
#    record = group_args(residue  = trgResidue,
#                        atom     = trgAtomname,
#                        dotcount = dotcount,
#                        mingap   = mingap,
#                        seqdist  = srcResidue.seq_dist(trgResidue))
#    srcResidue.probe[srcAtomname][recordkey] = record
#for srcatom in residue.probe:
 # for record in residue.probe[srcatom]:
  #  check_bond = residue.probe[srcatom][record]
#-------------------------------------------------------------------------------
#}}}

#{{{ fail index match
#-------------------------------------------------------------------------------
def fail_index_match(index, residue, motif_in_progress):
  if not index: return 0 #i.e. the dafault index of ''
  for observed_index in motif_in_progress.residues.keys():
    if index == observed_index and motif_in_progress.residues[index] != residue:
      return "Index "+index+" already identified, does not match new residue "+residue.id_with_resname()
    if motif_in_progress.residues[observed_index] == residue and observed_index != index:
      return "Residue "+residue.id_with_resname()+" already indexed, does not match new index "+index
  else: return 0
#-------------------------------------------------------------------------------
#}}}

#{{{ make pickle function
#This .pickles and stores a motif fingerprint object
#-------------------------------------------------------------------------------
def make_pickle(motif):
  pwd = os.getcwd()
  fingerprints_dir = libtbx.env.find_in_repositories(
    "cctbx_project/mmtbx/cablam/fingerprints")
  if fingerprints_dir is None:
    raise Sorry("""\
Problem locating cablam fingerprints dir""")
  os.chdir(fingerprints_dir)
  filename = motif.motif_name + ".pickle"
  print("Converting", motif.motif_name, "to pickle file . . .")
  easy_pickle.dump(file_name=filename,obj=motif)
  print(". . . Done")
  os.chdir(pwd)
#-------------------------------------------------------------------------------
#}}}

#{{{ fetch fingerprints function
#Given a list of string which match .pickled motif names, returns a list of
#  unpickled motif fingerprint objects
#-------------------------------------------------------------------------------
def fetch_fingerprints(motif_list):
  fingerprint_list = []
  for motif in motif_list:
    path = "cctbx_project/mmtbx/cablam/fingerprints/"+motif+".pickle"
    picklefile = libtbx.env.find_in_repositories(
      relative_path=path, test=os.path.isfile)
    if picklefile is None:
      raise Sorry("\nCould not find a needed pickle file for motif "+motif+" in chem_data.\nExiting.\n")
    else:
      fingerprint_list.append(easy_pickle.load(file_name=picklefile))
  return fingerprint_list
#-------------------------------------------------------------------------------
#}}}

#{{{ get all labels function
#Returns a list of all the residues labels for all the motifs in a list
#May be needed for some output functions
#-------------------------------------------------------------------------------
def get_all_labels(motif_name_list):
  motifs = fetch_fingerprints(motif_name_list)
  label_list = []
  for motif in motifs:
    for residue_name in motif.residue_names.values():
      label_list.append(residue_name)
  return label_list
#-------------------------------------------------------------------------------
#}}}
