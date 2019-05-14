from __future__ import absolute_import, division, print_function
# (jEdit options) :folding=explicit:collapseFolds=1:
#This module contains the linked_residue class and the functions needed to build
#  and access instances of it.
#2012-09-05:
#  prunerestype() moved to this module from cablam_training
#  linked_residue.id_with_resname() changed to return pdb-column-formatted ids
#2012-10-09:
#  self.resid is now stored in each linked_residue object
#2013_02_01: Added a step in linked_residue.__init__() that flags HETATOMs and a
#  step in construct_linked_residues() that skips adding them to the
#  resdata/protein object. Will this be a problem for synthetic or other
#  non-standard residues?
#2013-09-17: changed formatting for id_with_resname to rg.atom.pdb_label_columns
#  Tried to add handling for HETATOMS in __init__. May or may not fully work.
#Next: added mp_id() to produce MolPribity-friendly resids


import sys
from mmtbx.cablam.cablam_math import veclen, vectorize

#{{{ linked_residue class
#This class holds information on protein residues (as defined by residue_group
#  in pdb.hierarchy). See the __init__ for details. It also maintains sequence
#  relationships through a connectivity check and references "prevres" and
#  "nextres" to the class instances of sequence-adjacent residues. The intent is
#  to allow easy stepping forward and backward through the residue sequence.
#For example: self.resnum returns the current residue number,
#  self.nextres.resnum returns the residue number of the next residue in
#  sequence, if such a residue exists. This is not protected by a .get or
#  anything, so you have to do your own error catching, usually of the form:
#  if self.nextres: return self.nextres.resnum
#The class was designed for use with cablam, which requires forward and backward
#  sequence awareness, but it could be used elsewhere. It contains a generic
#  results={} dictionary  usable by anyone willing to make unique keys for their
#  data.
#It is recommended that all instances of this class for a protein be members of
#  some iterable (I use a dictionary keyed with the cablam_key() function in
#  this module). Do not rely on sequence connectivity alone (which breaks due to
#  multiple chains or missing atoms) to let you find all residues.
#The class should deal with insertion codes well, since different icodes seem to
#  be in different rg's in the hierarchy.
#The class does not deal with alts well, however, as an rg may contain several
#  ag's. The firstalt() function provides some relief, but is really a dodge
#  rather than a fix. I do not expect this to be a problem at low resolution
#  (where alts are rare), but it does present a problem for high-res training.
#-------------------------------------------------------------------------------
class linked_residue(object):

  #Prints kinemage point-style output for a list of measures given in kinorder
  def printtokin(self, kinorder, writeto=sys.stdout):
    outline = ['{'+self.pdbid+' '+self.id_with_resname()+'}']
    for order in kinorder:
      try:
        outline.append(str(self.measures[order]))
      except KeyError:
        return
    writeto.write(' '.join(outline)+'\n')

  #Prints comma-separated output for a list of measures given in kinorder
  def printtocsv(self, kinorder, doconnections=False, writeto=sys.stdout):
    outline = [self.id_with_resname()]
    for order in kinorder:
      try:
        outline.append(str(self.measures[order]))
      except KeyError:
        outline.append('NULL')
    if doconnections:
      if self.prevres: outline.append(self.prevres.id_with_resname())
      else: outline.append('NULL')
      if self.nextres: outline.append(self.nextres.id_with_resname())
      else: outline.append('NULL')
    writeto.write(','.join(outline)+'\n')  #note the ',' in the join here

  #id_to_string and id_with_resname return string concatenations of residue
  #  identifiers.  The identifier order should be standard with other RLab
  def id_to_str(self, sep=' '):
    resid_string = sep.join(
      [self.pdbid, self.model, self.chain, str(self.resnum), self.icode])
    return resid_string

  def id_with_resname(self):
    # Formatted as: 'ALA A####I'
    resid_string = self.rg.atoms()[0].pdb_label_columns()[5:]
    return resid_string

  def mp_id(self):
    #An id consistent with MolProbity 'cnit' ids
    #Formatted as: ccnnnnilttt
    #  c: 2-char Chain ID, space for none
    #  n: sequence number, right justified, space padded
    #  i: insertion code, space for none
    #  l: alternate ID, space for none
    #  t: residue type (ALA, LYS, etc.), all caps left justified, space padded
    #(Not sure about the 2-char Chain IDs just yet)
    #(alternates are not going to be handled properly yet)
    resid_string = self.id_with_resname()
    resname = resid_string[0:3]
    chain   = resid_string[3:5]
    resnum  = resid_string[5:9]
    ins     = resid_string[9:10]
    alt     = self.firstalt("CA")
    if alt is None or alt == '':
      alt = " "
    mpid_string = chain + resnum + ins + alt + resname
    return mpid_string

  #Removes the references that sequence-adjacent linked_residue class instances
  #  have to this instance. Helps maintain correct sequence connectivity and may
  #  allow this instance to be removed from memory.
  #Used in cablam_training.stripB() and cablam_training.prunerestype()
  def removelinks(self):
    if self.prevres:
      self.prevres.nextres = None
    if self.nextres:
      self.nextres.prevres = None

  #Returns the first alt index in alts that has an atom of the requested name
  #  Removes some guesswork from atom lookup, but really just acrobatics around
  #  the problem of how to store and access alts usefully
  def firstalt(self, atomname):
    for alt in self.alts:
      try:
        if self.atomxyz[alt][atomname]:
          return alt
        else:
          continue
      except KeyError:
        continue
    else:
      return None

  #simplified retrieval around firstalt for the common case of atom coords
  def getatomxyz(self, atomname):
    firstalt = self.firstalt(atomname)
    try:
      return self.atomxyz[firstalt][atomname]
    except KeyError:
      return None

  #There needs to be a CA-only consecutive check. Adding one is a high priority.
  def consecutive(self, res1, res2):
    if res1 and res2: #check against empties
      try:
        C = res1.atomxyz[res1.firstalt('C')]['C']
        N = res2.atomxyz[res2.firstalt('N')]['N']
        #firstalt returns None if it can't find an atom,
        #  and a key of None gives a KeyError here
      except KeyError:
        #if there aren't a C and an N, assume the two are not properly bonded
        return False
      bondlen = veclen(vectorize(C,N))
      if bondlen <= 2.0:
        #2.0A is the peptide bond cutoff used by O for model building and
        #  potentially-fragmented chains. O's generous cutoff seemed appropriate
        #  since I expect to process in-progress models with this program
        #RLab (probably) uses 1.4 +- 0.3, official textbook is about 1.33
        #O's Ca-Ca distance is 4.5A
        return True
      else:
        return False
    else:
      return False

  def seq_dist(self, otherres):
    ############################################
    # Returns distance in sequence with insertion codes accounted for.
    # The return value is negative if res is N-terminal to this.
    # The return value is positive if res is C-terminal to this.
    # The return value is None if res couldn't be found due to chain break etc.
    ############################################
    if self is otherres:
      return 0
    if self.chain != otherres.chain:# or self.model != otherres.model:
      return None
    #guess which direction to look
    #  the "<=" and ">=" should let this look back or forward from within an
    #  insertion
    if self.resnum <= otherres.resnum:
      delta = 0
      cur = self
      while cur != None:
        delta += 1
        if cur.nextres is otherres: return delta
        cur = cur.nextres
    if self.resnum >= otherres.resnum:
      delta = 0
      cur = self
      while cur != None:
        delta -= 1
        if cur.prevres is otherres: return delta
        cur = cur.prevres
    return None

  def __init__(self,
    rg, prevres=None, pdbid='pdbid', modelid='', chainid='',
    targetatoms=["CA","O","C","N"]
    ):

    self.rg = rg #the source residue group is preserved for additional data and
    #ease in transfering back to hierarchy mode

    self.pdbid = pdbid
    self.model = modelid
    self.chain = chainid
    self.resnum = int(rg.resseq.strip())
    #self.resseq = rg.resid[:-1]
    self.icode = rg.icode
    self.resid = cablam_key(self.model, self.chain, self.resnum, self.icode)
    self.hetero = False #marks whether this is a HETATOM

    #alts: 'alt' and 'resname' keyed by ag.altloc in the form of '','A','B' etc.
    #atomxyz: xyz coords, indexed by ag.altloc, and atom.name within each alt
    #  e.g. atomxyz['']['CA'] returns the coordinates of a non-alt Calpha
    #atomb: atomic b, indexed by ag.altloc, and atom.name within each alt
    self.alts = {}
    self.atomxyz = {}
    self.atomb   = {}
    ### What about anisou? Not handled yet.

    #hierachy looping and data extraction
    for ag in rg.atom_groups():
      #if not ag.is_protein(): Need sopmething like this that works
      #  self.is_protein=True
      self.alts[ag.altloc] = {'alt':ag.altloc, 'resname':ag.resname}
      self.atomxyz[ag.altloc] = {}
      self.atomb[ag.altloc]   = {}
      for atom in ag.atoms():
        if atom.hetero and ag.resname.upper() != 'MSE':
           self.hetero=True
        for targetatom in targetatoms:
          if atom.name.strip() == targetatom:
            self.atomxyz[ag.altloc][targetatom] = atom.xyz
            self.atomb[ag.altloc][targetatom] = atom.b

    #Note that a reference to the related residue is stored, not a dictionary
    #  key for the wrapper dictionary
    #Someone clever may want to teach me how to use weakref() if the mutual
    #  references that result from this cause memory problems
    if prevres and self.consecutive(prevres, self):
      self.prevres = prevres     #Connect this residue to previous
      prevres.nextres = self     #And the previous residue to this one
    else:
      self.prevres = None #Adjacency is handled in an outside function
    self.nextres = None

    self.probe = {'O':{},'H':{}}   #holder for hydrogen bonding, indexed by 'target' residue+atom, see cablam_training.add_probe_data()
    self.probeH = []
    self.probeO = []

    self.measures = {} #Holder for cablam-space geometric measures
    self.motifs = {} #Holder for identified motifs from fingerprints/probetrain
    self.results = {} #Generic holder for calcuated values of interest

#-------------------------------------------------------------------------------
#}}}

#{{{ prunerestype function
#Deletes all members of a given residue type from a dictionary of residues
#"Residue type" is determined from pdb.hierarchy's ag.resname, the upshot being
#  that non-residues like "HOH" waters and het groups can also be pruned to
#  improve performance if their three-letter codes are known.
#-------------------------------------------------------------------------------
def prunerestype(resdata, restype):
  # TODO is resdata a dict ? put a hint
  reslist = list(resdata.keys())
  for residue in reslist:
    for alt in resdata[residue].alts:
      if resdata[residue].alts[alt]['resname'].strip() == restype:
        resdata[residue].removelinks()
        trash = resdata.pop(residue)
        break
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_key function
#The "protein" or "resdata" dictionary returned by construct_linked_residues()
#  below uses a particular key construction to access (and order) its contents
#This function provides that construction so that residues in resdata may be
#  accessed from anywhere. The string returned is .sort()able
#-------------------------------------------------------------------------------
def cablam_key(modelid=None, chainid=None, resnum=None, icode=None):
  if None not in [modelid, chainid, resnum, icode]:
    resid_string = ' '.join([modelid, chainid, '%04i' % resnum, icode])
    #The bit of string formatting here ('%04i' % resnum) helps .sort() later by
    #  adding 0's to the left side of resnum until resnum is 4 characters long.
    #  May or may not be compatible with Hybrid36 or other numbering schemes.
    return resid_string
  else:
    sys.stderr.write("""
Missing value for cablam_res.cablam_key(pdbid, modelid, chainid, resnum, icode)
Please pass complete information
      """)
    sys.exit()
#-------------------------------------------------------------------------------
#}}}

#{{{ construct_linked_residues function
#-------------------------------------------------------------------------------
#This function returns a dictionary of linked_residue objects with keys as
#  defined by cablam_key() above. It is responsible for iterating over part of
#  the hierarchy, but most of the work is done by the __init__() in the
#  linked_residue class.
#targetatoms handling is likely to see some change as I add compatibility for
#  CA-only mainchain traces.  Currently "C" and "N" are necessary for the
#  sequency adjacency check linked_residue.consecutive(), which checks for
#  appropriate peptide bond length
#targetatoms is a list of strings that identify pdb atoms, e.g. "CA" or "CD1"
def construct_linked_residues(
  hierarchy, targetatoms=["CA","O","C","N"], pdbid='pdbid'
  ):
  protein = {}
  for model in hierarchy.models():
    for chain in model.chains():
      prevres = None
      for rg in chain.residue_groups():
        residue = linked_residue(
          rg, prevres=prevres, pdbid=pdbid, modelid=model.id, chainid=chain.id,
          targetatoms=targetatoms)
        resid_string = cablam_key(model.id, chain.id, residue.resnum, rg.icode)
        if not residue.hetero: #automatically skip het atoms
          protein[resid_string] = residue
          prevres = residue   #important update for determining connectivity
  return protein
#-------------------------------------------------------------------------------
#}}}
