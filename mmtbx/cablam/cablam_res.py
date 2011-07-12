# (jEdit options) :folding=explicit:collapseFolds=1:
#This module contains the linked_residue class and the functions needed to build
#  and access instances of it.

import sys
from cjw_vectormath import veclen, vectorize

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
#  sequence awareness. But it could be used elsewhere. It contains a generic
#  results={} dictionary used by cablam, but usable by anyone willing to make
#  unique keys for their data.
#It is recommended that all instances of this class for a protein be members of
#  some iterable (I use a dictionary keyed with the cablam_key() function in
#  this module). Do not rely on sequence connectivity (which breaks due to
#  multiple chains or missing atoms) to let you find all residues.
#The class should deal with insertion codes well, since different icodes seem to
#  be in different rg's in the hierarchy.
#The class does not deal with alts well, however, as an rg may contain several
#  ag's. The firstalt() function provides some relief, but is really a dodge
#  rather than a fix. I do not expect this to be a problem at low resolution
#  (where alts are rare), but it does present a problem for high-res training.
#-------------------------------------------------------------------------------
class linked_residue():

  #printrec is a simple debugging function to print the contents of a pdb to
  #  screen to give an idea of where and how values are stored
  def printrec(self):
    print self.pdbid, self.model, self.chain, self.resnum, self.icode
    for alt in self.alts:
      print self.alts[alt]['resname'], self.alts[alt]['alt']
      print ' ', self.atomxyz[alt], self.atomb[alt]

  def printtokin(self, kinorder, writeto=sys.stdout):
    outline = ['{'+self.id_with_resname(sep=' ')+' '+self.dssp+'}']
    for order in kinorder:
      try:
        outline.append(str(self.results[order]))
      except KeyError:
        return
    writeto.write(' '.join(outline)+'\n')

  def printtocsv(self, kinorder, writeto=sys.stdout):
    outline = [self.id_with_resname(sep=':')]
    for order in kinorder:
      try:
        outline.append(str(self.results[order]))
      except KeyError:
        return
    outline.append(self.dssp)
    writeto.write(','.join(outline)+'\n')  #note the ',' in the join here

  def id_to_str(self, sep=' '):
    resid_string = sep.join(
      [self.pdbid, self.model, self.chain, str(self.resnum), self.icode])
    return resid_string

  def id_with_resname(self, alt='', sep=' '):
    try:
      resname = self.alts[alt]['resname']
    except KeyError:
      resname = 'XXX'
    resid_string = sep.join(
      [self.pdbid, self.model, self.chain, resname+alt,
      str(self.resnum), self.icode])
    return resid_string

  #Removes the references that sequence-adjacent linked_residue class instances
  #  have to this instance. Helps maintain correct sequence connectivity and may
  #  allow this instance to be removed from memory.
  #Used in cablam_training.stripB() and cablam_training.prunerestype()
  def removelinks(self):
    if self.prevres:
      self.prevres.nextres = None
    if self.nextres:
      self.nextres.prevres = None

  #Returns the first alt in alts that has an atom of the requested name
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

  def dssprob_sort(self):
    self.dssprob_decending = self.dssprob.keys()  #list of dsspkeys
    self.dssprob_decending.sort(reverse=True, key=lambda x: self.dssprob[x])

  def __init__(self,
    rg, prevres=None, pdbid='pdbid', modelid='', chainid='',
    targetatoms=["CA","O","C","N"]
    ):

    self.pdbid = pdbid
    self.model = modelid
    self.chain = chainid
    self.resnum = int(rg.resseq.strip())
    self.icode = rg.icode

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
      self.alts[ag.altloc] = {'alt':ag.altloc, 'resname':ag.resname}
      self.atomxyz[ag.altloc] = {}
      self.atomb[ag.altloc]   = {}
      for atom in ag.atoms():
        for targetatom in targetatoms:
          if atom.name.strip() == targetatom:
            self.atomxyz[ag.altloc][targetatom] = atom.xyz
            self.atomb[ag.altloc][targetatom] = atom.b

    #Note that a reference to the related residue is stored, not a dictionary
    #  key for the wrapper dictionary
    #Someone clever may want to teach me how to use weakref() if the mutual
    #  references that reault from this cause memory leaking or something
    if prevres and self.consecutive(prevres, self):
      self.prevres = prevres     #Connect this residue to previous
      prevres.nextres = self     #And the previous residue to this one
    else:
      self.prevres = None #Adjacency is handled in an outside function
    self.nextres = None

    #self.resname = self.alts[self.firstalt('CA')]['resname']

    self.dssp = ''    #space fo DSSP code if needed

    self.results = {} #Holder for calcuated values of interest
    self.dssprob = {} #Holder for probabilities of various secondary structures
    #Yes, dssprob and dsspath are the sorts of over-clever names for which I
    #  should probably be scolded.
    self.dssprob_decending = []

#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_key function
#The "protein" or "resdata" dictionary returned by rollongRecords() below uses
#a   particular key construction to access (and order) its contents
#This function provides that construction so that residues in resdata may be
#  accessed from anywhere. The string returned is .sort()able
#-------------------------------------------------------------------------------
def cablam_key(modelid=None, chainid=None, resnum=None, icode=None):
  if None not in [modelid, chainid, resnum, icode]:
    resid_string = ' '.join([modelid, chainid, '%04i' % resnum, icode])
    #The bit of string formatting here ('%04i' % resnum) helps .sort() later by
    #  adding 0's to the left side of resnum until resnum is 4 characters long.
    #  May or may not be campatible with Hybrid36 or other numbering schemes.
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
        protein[resid_string] = residue
        prevres = residue   #important update for determining connectivity
  return protein
#-------------------------------------------------------------------------------
#}}}
