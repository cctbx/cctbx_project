from __future__ import division
# (jEdit options) :folding=explicit:collapseFolds=1:
import os, sys
import libtbx.load_env
from libtbx import easy_pickle
from libtbx.utils import Sorry

#{{{ class objects
class motif(object):
  def __init__(self,
    motif_name="",residue_names={}):
    self.motif_name = motif_name
    self.residue_names = residue_names #labels for printing keyed by indices for residues in the motif
    self.residues = []

  def add_residue(self, allowed_resname=[], banned_resname=[],
      sequence_move=None, bond_move='',
      end_of_motif=False, index=''):
    new_residue = residue(allowed_resname, banned_resname, sequence_move,
      bond_move, end_of_motif, index)
    self.residues.append(new_residue)
    return new_residue

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
    allowed_resname=[], banned_resname=[],
    sequence_move=None, bond_move='',
    end_of_motif=False, index=''):

    self.allowed_resname=allowed_resname
    self.banned_resname=banned_resname
    self.sequence_move=sequence_move
    self.bond_move=bond_move
    self.end_of_motif=end_of_motif
    self.index=index
    self.bonds=[]

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

class target_atom(object):
  def __init__(self,
    atomname = None, anyatom=False, seqdist=None, anyseqdist=False):
    self.atomname = atomname
    self.anyatom = anyatom
    self.seqdist = seqdist
    self.anyseqdist = anyseqdist
#}}}

def check_protein(protein, motif_name_list):
  motif_list = fetch_fingerprints(motif_name_list)
  for resid in protein:
    residue = protein[resid]
    for motif in motif_list:
      check_for_motif(motif,residue)

def check_for_motif(motif, active_residue):
  #This checks a single residue for whether it qualifies as the start point for a motif
  #Its 'return' as such is an addition to the residue.motifs {}
  #Format {'motifname':True} for reasons of output formatting, I guess
  motif_in_progress = {}
  for motif_residue in motif.residues:
    if not active_residue: #probably due to do_move returning None
      return False
    if not do_residue_check(motif_residue,active_residue):
      return False
    for bond in motif_residue.bonds:
      if bond.banned:
        if not check_for_forbidden_bonds(bond,active_residue):
          return False
      elif bond.required:
        bonding_partner = do_bond_check(bond,active_residue, motif_in_progress) #returns either False or a residue
        if bonding_partner:
          #if not check_index_match(bond.trg_index, bonding_partner, motif_in_progress):
          #  return False #move this inside the check loop so bifur can be handled
          if bond.trg_index: #i.e. not the default '' index
            motif_in_progress[bond.trg_index] = bonding_partner
        else:
          return False #return to previous function
    #add current residue to the motif in progress . . .
    if not check_index_match(motif_residue.index, active_residue, motif_in_progress):
      return False
    motif_in_progress[motif_residue.index] = active_residue
    # . . . move to the next residue . . .
    if not motif_residue.end_of_motif:
      #print active_residue
      active_residue = do_move(motif_residue, active_residue, motif_in_progress)
    else:
      sys.stderr.write("finished motif "+ motif.motif_name+"\n")
    # . . . and go back to top of 'for' loop
  #---------
  #This loop only occurs if the checking loop above has complete without return
  for index in motif_in_progress:
    if index == '' or index not in motif.residue_names.keys():
      continue
    residue = motif_in_progress[index]
    residue.motifs[motif.residue_names[index]] = True

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
      return motif_in_progress[move_index]
    except KeyError:
      sys.stderr.write('\nTried to move to a residue not yet indexed at index \"'+check_residue.bond_move+'\"\n')
      sys.stderr.write('Please check fingerprint definition. Exiting . . .\n')
      sys.exit()
  else:
    sys.stderr.write('\nNo move specified for residue indexed as\"'+check_residue.index+'\"\n')
    sys.stderr.write('Please check fingerprint definition. Exiting . . .\n')
    sys.exit()

def do_residue_check(motif_residue, active_residue):
  resname = active_residue.id_with_resname()[0:3] #should slice out the resname from this id
  if motif_residue.allowed_resname: #default is all-inclusive
    #for allowed in motif_residue.allowed_resname:
    #  if resname.upper() == allowed:
    #    return True
    #else:
    #  return False
    if resname.upper() not in motif_residue.allowed_resname:
      #sys.stderr.write(active_residue.id_with_resname() + '\n')
      return False
  if resname.upper() in motif_residue.banned_resname:
    return False
  return True

def do_bond_check(check_bond, active_residue, motif_in_progress):
  if check_bond.required:
    if not active_residue.probe or check_bond.src_atom not in active_residue.probe:
      return False
    if len(active_residue.probe[check_bond.src_atom]) > 1 and not check_bond.allow_bifurcated:
      return False
    for src_atom in active_residue.probe:
      if src_atom == check_bond.src_atom:
        for bond_name in active_residue.probe[src_atom]:
          real_bond = active_residue.probe[src_atom][bond_name]
          for trg_atom in check_bond.trg_atoms:
            if (trg_atom.anyatom or trg_atom.atomname == real_bond.atom) and (trg_atom.anyseqdist or trg_atom.seqdist == real_bond.seqdist):
              #print active_residue.id_with_resname(), src_atom,'|',real_bond.seqdist,'|',real_bond.residue.id_with_resname(), bond_name
              if check_index_match(check_bond.trg_index, real_bond.residue, motif_in_progress):
                return real_bond.residue
    #bonding_partner = check_for_required_bonds(check_bond, active_residue)
  return False

#At the moment, does not check for generic bans against indexed residues w/o
#  specific sequence relationship
def check_for_forbidden_bonds(check_bond, active_residue):
  if not active_residue.probe:
    return True #no problem with an empty set
  if check_bond.src_atom in active_residue.probe.keys():
    src_atom = check_bond.src_atom
    for bond_name in active_residue.probe[src_atom]:
      real_bond = active_residue.probe[src_atom][bond_name]
      for trg_atom in check_bond.trg_atoms:
        if trg_atom.anyatom and trg_atom.anyseqdist:
          return False
        elif trg_atom.anyatom and (trg_atom.seqdist == real_bond.seqdist):
          return False
        elif  trg_atom.anyseqdist and (trg_atom.atomname == real_bond.atom):
          return False
        elif (trg_atom.atomname == real_bond.atom) and (trg_atom.seqdist == real_bond.seqdist):
          return False
        else:
          continue
        #I think I can do that in one line.  Should I?:
        #if (trg_atom.anyatom or trg_atom.atomname == active_residue.probe[real_bond].atom) and (trg_atom.anyseqdist or trg_atom.seqdist == active_residue.probe[real_bond].seqdist):
        #  return False
  return True #only if the for loop completes without issue

###def check_for_required_bonds(check_bond, active_residue):
###  sys.stderr.write('checking')
####  if not active_residue.probe:
####    return False
####  if len(active_residue.probe[check_bond.src_atom]) > 1 and not check_bond.allow_bifurcated:
####    return False
###  for src_atom in active_residue.probe:
###    for bond_name in active_residue.probe[src_atom]:
###      real_bond = active_residue.probe[src_atom][bond_name]
###      for trg_atom in check_bond.trg_atoms:
###        #Four possible cases:
###        #1: a bond is needed, and any will do
###        #2: a bond is needed in certain sequence relationship
###        #3: a bond is needed to certain atom type
###        #4: a bond is needed to certain atom type in certain sequence relationship
###        if (trg_atom.anyatom or trg_atom.atomname == real_bond.atom) and (trg_atom.anyseqdist or trg_atom.seqdist == real_bond.seqdist):
###          print 'found'
###          #Then check if the matching bond is bifurcated when it should not be
###          if len(active_residue.probe[src_atom]) > 1 and not check_bond.allow_bifurcated:
###            return False
###          else:
###            #print real_bond
###            #print real_bond.residue
###            print active_residue.id_with_resname(),'|',realbond.seqdist,'|',real_bond.residue.id_with_resname()
###            return real_bond.residue #need data for indexing etc.
###  return False

def check_index_match(index, residue, motif_in_progress):
  if not index: #i.e. the default index of ''
    return True
  for test_index in motif_in_progress:
    if index == test_index:
      if motif_in_progress[index] != residue: #index present w/o residue
        return False
    if motif_in_progress[test_index] == residue:
      if test_index != index: #residue present x/o index
        return False
  return True #either residue/index match or both totally absent

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

def make_pickle(motif):
  pwd = os.getcwd()
  fingerprints_dir = libtbx.env.find_in_repositories(
    "cctbx_project/mmtbx/cablam/fingerprints")
  if fingerprints_dir is None:
    raise Sorry("""\
Problem locating cablam fingerprints dir""")
  os.chdir(fingerprints_dir)
  filename = motif.motif_name + ".pickle"
  print "Converting", motif.motif_name, "to pickle file . . ."
  easy_pickle.dump(file_name=filename,obj=motif)
  print ". . . Done"
  os.chdir(pwd)

def fetch_fingerprints(motif_list):
  #give this a list of strings; it will return a list of motif objects
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

def get_all_labels(motif_name_list):
  motifs = fetch_fingerprints(motif_name_list)
  label_list = []
  for motif in motifs:
    for residue_name in motif.residue_names.values():
      label_list.append(residue_name)
  return label_list
