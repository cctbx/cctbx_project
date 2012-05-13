# (jEdit options) :folding=explicit:collapseFolds=1:
import sys

#{{{ fingerprint superclass, hold top-level methods and data
class ssfingerprint():

  #Returns bool True if two residues are in the same model and chain
  def in_same_chain(self, srcres, trgres):
    if (trgres.model == srcres.model) and (trgres.chain == srcres.chain):
      return True
    else:
      return False

  #Returns sequence separation of two residues
  #  Does not check chain id or handle insertion codes
  def find_bond_jump(self, srcres, trgres):
    return trgres.resnum - srcres.resnum

  #Wraps residue-level checking for this fingerprint
  #Good luck if you're trying to suss out the logic here.
  def checkall(self, startres, debug=False):
    foundres = {}
    if debug: sys.stderr.write('checking '+str(startres.resnum)+'\n')
    #adds to residue.motifs True or False, keyed by the string name of the fingerprint
    currentres = startres
    for memberres in self.members:
      if debug: print ' member '
      if currentres: #should catch a stray None
        if currentres == 'Done': break
        if memberres.checkthis(currentres, foundres):
          currentres = memberres.do_move(currentres,foundres)
        else:
          if debug: print 'checkthis returned false'
          break
      else:
        if debug: print 'currentres is None'
        break

    if currentres == 'Done':
      for index in foundres:
        if index in self.labelhash:
          foundres[index].motifs[self.labelhash[index]] = True
        else:
          continue

  def __init__(self, name):
    self.members = []
    self.name = name
    self.labelhash = {} #holds output labels for assignment to each residue
    #  of interest
    self.labellist = [] #holds member residue labels in output-relevant order
#}}}

#{{{ member_residue subclass, holds residue-level methods and data
class member_residue():

  def in_same_chain(self, srcres, trgres):
    if (trgres.model == srcres.model) and (trgres.chain == srcres.chain):
      return True
    else:
      return False

  def find_bond_jump(self, srcres, trgres):
    return trgres.resnum - srcres.resnum

  #returns True if the residue type for this residue is acceptible
  def resnamecheck(self, srcres):
    if self.resaccept:
      firstalt = srcres.firstalt('CA')
      resname = srcres.alts[firstalt]['resname']
      if resname.upper() not in resaccept:
        return False
    elif self.resbanned:
      firstalt = srcres.firstalt('CA')
      resname = srcres.alts[firstalt]['resname']
      if resname.upper() in resbanned:
        return False
    else:
      #This happens if the residue passes the check or if no check was specified
      return True

  #{{{ Check for correct index
  def indexcheck(self, residue, index, foundres):
    #Four cases:
    #1 no index, no res = pass
    #2 correctly paired index and res = pass
    #3 unique index w/wrong res = fail
    if index and index in foundres.keys():
      if foundres[index] != residue: return False
      else: pass #okay
    else: pass
    #4 nonunique res w/wrong index
    if index and residue in foundres.values():
      for key in foundres:
        if foundres[key] == residue and key != index: return False
        else: pass
    else:
      pass
    return True #if everything checks out
  #}}}

  #{{{ checkthis - all checks for current residue
  def checkthis(self, srcres, foundres, debug=False):
    if not self.indexcheck(srcres, self.index, foundres): return False
    if debug: print '  index pass '

    if self.resaccept or self.resbanned:
      if not resnamecheck(srcres): return False
      else: pass
    else: pass

    if self.Obonding:
      if debug: print '  try Obond '
      #verify Obonding
      if len(self.Obonding) == 1:
        bondjump = self.Obonding[0]
        if bondjump == 'any': #0 or 1 bonds, no preference which
          if len(srcres.probeO) > 1: return False
          else: pass
        elif bondjump == "!":
          if len(srcres.probeO) > 0: return False
          else: pass
        elif len(srcres.probeO) != 1:
          if debug: print '  FAIL: Wrong bond count', len(srcres.probeO)
          return False
        else:
          #single trgres confirmed, procede with testing
          trgres = srcres.probeO[0]
          if str(bondjump).isalpha():
            #check that bondjump is an okay index for trgres
            if not self.indexcheck(trgres, bondjump, foundres):
              if debug: print '  FAIL: bond has wrong index'
              return False
            #if it gets through that, then update foundres
            foundres[self.index] = srcres
            foundres[bondjump] = trgres #bondjump *is* the trgres index, as specified in the fingerprint definition
            pass #I'm using 'pass' to help me keep tabs on where bits of logic end
          else: #(if bondjump is an int)
            if self.find_bond_jump(srcres,trgres) != bondjump: return False
            else:
              if not self.indexcheck(srcres, self.index, foundres): return False #Redundant?
              else:
                foundres[self.index] = srcres
                #Don't know the index for trgres
                pass
      elif len(self.Obonding) == 2:
        if len(srcres.probeO) > 2: return False
        else: pass #The above checks against too many bonds
        requestsfilled = [] #Will hold T/F values to be checked later
        requestidx = -1
        for bondjump in self.Obonding:
          requestidx += 1 #(First idx used will be 0)
          requestsfilled.append(False)
          if bondjump == 'any': #an 'any' is an automatic pass
            requestsfilled[requestidx] = True
            continue
          elif str(bondjump).isalpha():
            for trgres in srcres.probeO:
              if not self.indexcheck(trgres,bondjump,foundres): continue
              else:
                foundres[bondjump] = trgres
                requestsfilled[requestidx] = True
                break
            pass
          else: #(if bondjump is an int)
            for trgres in srcres.probeO:
              if self.find_bond_jump(srcres,trgres) != bondjump:
                continue
              else: # more checks
                if not self.indexcheck(srcres, self.index, foundres):
                  continue
                else:
                  requestsfilled[requestidx] = True
                  break #exit the for loop
        for request in requestsfilled:
          if not request: return False
        foundres[self.index] = srcres
        pass #end of logic block

      else:
        pass

    if self.Hbonding:
      if debug: print '  try Hbond '
      #verify Hbonding
      if len(self.Hbonding) == 1:
        bondjump = self.Hbonding[0]
        if bondjump == 'any': #0 or 1 bonds, no preference which
          if len(srcres.probeH) > 1: return False
          else: pass
        elif bondjump == "!":
          if len(srcres.probeO) > 0: return False
          else: pass
        elif len(srcres.probeH) != 1:
          if debug: print '  FAIL: Wrong bond count', len(srcres.probeH)
          return False
        else:
          #single trgres confirmed, procede with testing
          trgres = srcres.probeH[0]
          if str(bondjump).isalpha():
            #check that bondjump is an okay index for trgres
            if not self.indexcheck(trgres, bondjump, foundres):
              if debug: print '  FAIL: bond has wrong index'
              return False
            #if it gets through that, then update foundres
            foundres[self.index] = srcres
            foundres[bondjump] = trgres #bondjump *is* the trgres index, as specified in the fingerprint definition
            pass #I'm using 'pass' to help me keep tabs on where bits of logic end
          else: #(if bondjump is an int)
            if debug: print '  **doing check'
            if debug: print self.find_bond_jump(srcres,trgres)
            if self.find_bond_jump(srcres,trgres) != bondjump: return False
            else:
              if not self.indexcheck(srcres, self.index, foundres): return False #Redundant?
              else:
                foundres[self.index] = srcres
                #Don't know the index for trgres
                pass
      elif len(self.Hbonding) == 2:
        if len(srcres.probeH) > 2: return False
        else: pass #The above checks against too many bonds
        requestsfilled = [] #Will hold T/F values to be checked later
        requestidx = -1
        for bondjump in self.Hbonding:
          requestidx += 1 #(First idx used will be 0)
          requestsfilled.append(False)
          if bondjump == 'any': #an 'any' is an automatic pass
            requestsfilled[requestidx] = True
            continue
          elif str(bondjump).isalpha():
            for trgres in srcres.probeH:
              if not self.indexcheck(trgres,bondjump,foundres): continue
              else:
                foundres[bondjump] = trgres
                requestsfilled[requestidx] = True
                break
            pass
          else: #(if bondjump is an int)
            for trgres in srcres.probeH:
              if self.find_bond_jump(srcres,trgres) != bondjump:
                continue
              else: # more checks
                if not self.indexcheck(srcres, self.index, foundres):
                  continue
                else:
                  requestsfilled[requestidx] = True
                  break #exit the for loop
        for request in requestsfilled:
          if not request: return False
        foundres[self.index] = srcres
        pass #end of logic block

      else:
        pass

    #Obanned and Hbanned list hydrogen bonding to be specifically disallowed
    if self.Obanned:
      for banned_bond in self.Obanned:
        if banned_bond.isalpha():
          if banned_bond == 'any' and len(srcres.probeO) > 0: return False
          else: pass
          for trgres in srcres.probeO:
            if banned_bond in foundres and foundres[banned_bond] == trgres:
              return False
            else:
              pass
        else:
          for trgres in srcres.probeO:
            if self.find_bond_jump(srcres,trgres) == banned_bond: return False
            else: pass

    if self.Hbanned:
      for banned_bond in self.Hbanned:
        if banned_bond.isalpha():
          if banned_bond == 'any' and len(srcres.probeH) > 0: return False
          else: pass
          for trgres in srcres.probeH:
            if banned_bond in foundres and foundres[banned_bond] == trgres:
              return False
            else:
              pass
        else:
          for trgres in srcres.probeH:
            if self.find_bond_jump(srcres,trgres) == banned_bond: return False
            else: pass

    if debug: sys.stderr.write(' check pass \n')
    foundres[self.index] = srcres
    return True #if it gets through all of the above without hitting a return False, it must be okay
    #Egad, the whitespace in this function is a nightmare!

  #}}}

  #{{{ do_move method
  def do_move(self, srcres, foundres, debug=False):
    if not self.move: #An empty or None move indicated motif's end
      if debug: print "   ***Done***"
      return "Done"
    elif str(self.move).isalpha():
      return foundres[self.move]
    else:
      seqdist = self.move
      trgres = srcres
      if seqdist > 0:
        while seqdist > 0:
          if not trgres.nextres:
            if debug: sys.stderr.write('  move fail: no nextres\n')
            return None
          else:
            trgres = trgres.nextres
            seqdist -= 1
      if seqdist < 0:
        while seqdist < 0:
          if not trgres.prevres:
            if debug: sys.stderr.write('  move fail: no prevres\n')
            return None
          else:
            trgres = trgres.prevres
            seqdist += 1
      return trgres
  #}}}

  #{{{ __init__ method
  def __init__(self):
    self.resaccept = []
    self.resbanned = []
    #^formatted as GLY,PRO,etc. (3-letter codes for flexibility) case should be set to .upper() on use
    #negtives are held in resbanned
    #Empty is all-permissive
    self.Obonding = []
    self.Hbonding = []
    #^required bonding patterns, may gain explicit strength requirements later
    #formatted as AtomnameBondingindex, so O4 is oxygen bonding to i+4, N-3 is amide H bonding to i-3
    #this bonding nomenclature matches what's currently used in my probe interpretter
    self.Obanned = []
    self.Hbanned = []
    #^bonds that are NOT permitted in the pattern
    self.geometry = []
    #this is a placeholder for geometic information about this motif
    #formatting unknown, it may end up as a link to a contour file
    self.distance = []
    #^this is a placeholder for spacial, non-sequence relationships among residues
    #formatting undetermined
    self.move = []
    #^this describes the move necessary to get to the next residue to be considered
    self.index = None
    #Important note! All indices must pass string.isalpha() as True!
    ##do_move() and checkthis() (maybe others) are dependent on this formatting!
  #}}}
#}}}

### I'm bailing on handling bifurcation in the program,
### for now, write fingerprints carefully, and the indices will handle the rest

#Labels should be unique and tied to their parent fingerprint.

#{{{ alpha helix definitions
#-------------------------------------------------------------------------------
alpha_helix_single = ssfingerprint('alpha_helix_single')
alpha_helix_single.labelhash = {'a':'alpha_helix_single'}
alpha_helix_single.labellist = ['alpha_helix_single']
alpha_helix_single.members.append(member_residue())
alpha_helix_single.members[0].Obonding = [4]
alpha_helix_single.members[0].Hbonding = [-4]
alpha_helix_single.members[0].index = 'a'
alpha_helix_single.members[0].move = None

alpha_helix_3 = ssfingerprint('alpha_helix_3')
alpha_helix_3.labelhash = {
'a':'alpha_helix_3_nside',
'b':'alpha_helix_3_center',
'c':'alpha_helix_3_cside'}
alpha_helix_3.labellist = [
'alpha_helix_3_nside',
'alpha_helix_3_center',
'alpha_helix_3_cside']
alpha_helix_3.members.append(member_residue())
alpha_helix_3.members[0].Obonding = [4]
alpha_helix_3.members[0].move = 1 #a number means sequence, a string means index
alpha_helix_3.members[0].index = 'a'
alpha_helix_3.members.append(member_residue())
alpha_helix_3.members[1].Obonding = [4]
alpha_helix_3.members[1].move = 1 #a number means sequence, a string means index
alpha_helix_3.members[1].index = 'b'
alpha_helix_3.members.append(member_residue())
alpha_helix_3.members[2].Obonding = [4]
alpha_helix_3.members[2].move = None #a number means sequence, a string means index
alpha_helix_3.members[2].index = 'c'

alpha_turn = ssfingerprint('alpha_turn')
alpha_turn.labelhash = {
'a':'alpha_turn',
'e':'alpha_turn'}
alpha_turn.labellist = [
'alpha_turn']
alpha_turn.members.append(member_residue())
alpha_turn.members[0].Obonding = [4]
alpha_turn.members[0].Hbonding = []
alpha_turn.members[0].move = 1
alpha_turn.members[0].index = 'a'
alpha_turn.members.append(member_residue())
alpha_turn.members[1].Obonding = [4]
alpha_turn.members[1].Hbonding = []
alpha_turn.members[1].move = 1
alpha_turn.members[1].index = 'b'
alpha_turn.members.append(member_residue())
alpha_turn.members[2].Obonding = []
alpha_turn.members[2].Hbonding = []
alpha_turn.members[2].move = 1
alpha_turn.members[2].index = 'c'
alpha_turn.members.append(member_residue())
alpha_turn.members[3].Obonding = []
alpha_turn.members[3].Hbonding = [-4]
alpha_turn.members[3].move = None
alpha_turn.members[3].index = 'd'
alpha_turn.members.append(member_residue())
alpha_turn.members[4].Obonding = []
alpha_turn.members[4].Hbonding = [-4]
alpha_turn.members[4].move = None
alpha_turn.members[4].index = 'e'

reg_alpha = ssfingerprint('reg_alpha')
reg_alpha.labelhash = {
'a':'reg_alpha',
'b':'reg_alpha',
'c':'reg_alpha',
'd':'reg_alpha'}
reg_alpha.labellist = ['reg_alpha']
#There seems at present no reason to distinguish among these residues
reg_alpha.members.append(member_residue())
reg_alpha.members[0].Obonding = [4]
reg_alpha.members[0].Hbonding = [-4]
reg_alpha.members[0].move = 1
reg_alpha.members[0].index = 'a'
reg_alpha.members.append(member_residue())
reg_alpha.members[1].Obonding = [4]
reg_alpha.members[1].Hbonding = [-4]
reg_alpha.members[1].move = 1
reg_alpha.members[1].index = 'b'
reg_alpha.members.append(member_residue())
reg_alpha.members[2].Obonding = [4]
reg_alpha.members[2].Hbonding = [-4]
reg_alpha.members[2].move = 1
reg_alpha.members[2].index = 'c'
reg_alpha.members.append(member_residue())
reg_alpha.members[3].Obonding = [4]
reg_alpha.members[3].Hbonding = [-4]
reg_alpha.members[3].move = None
reg_alpha.members[3].index = 'd'
#-------------------------------------------------------------------------------
#}}}

#{{{ alpha helix interrupts
#-------------------------------------------------------------------------------
wide_helix_turn = ssfingerprint('wide_helix_turn')
wide_helix_turn.labelhash = {
'a':'wide_helix_turn_helix_in',
'b':'wide_helix_turn_bifur',
'c':'wide_helix_turn_after_bifur',
'd':'wide_helix_turn_miss_bond_1',
'e':'wide_helix_turn_miss_bond_2',
'f':'wide_helix_turn_alpha_partner',
'g':'wide_helix_turn_pi_partner'}
wide_helix_turn.labellist = [
'wide_helix_turn_helix_in',
'wide_helix_turn_bifur',
'wide_helix_turn_after_bifur',
'wide_helix_turn_miss_bond_1',
'wide_helix_turn_miss_bond_2',
'wide_helix_turn_alpha_partner',
'wide_helix_turn_pi_partner']
wide_helix_turn.members.append(member_residue())
wide_helix_turn.members[0].Obonding = [4]
wide_helix_turn.members[0].move = 1
wide_helix_turn.members[0].index = 'a'
wide_helix_turn.members.append(member_residue())
wide_helix_turn.members[1].Obonding = [4,5]
wide_helix_turn.members[1].move = 1
wide_helix_turn.members[1].index = 'b'
wide_helix_turn.members.append(member_residue())
wide_helix_turn.members[2].Obonding = [5]
wide_helix_turn.members[2].move = 1
wide_helix_turn.members[2].index = 'c'
wide_helix_turn.members.append(member_residue())
#should there be an Obanned here?
wide_helix_turn.members[3].move = 1
wide_helix_turn.members[3].index = 'd'
wide_helix_turn.members.append(member_residue())
#should there be an Obanned here?
wide_helix_turn.members[4].move = 1
wide_helix_turn.members[4].index = 'e'
wide_helix_turn.members.append(member_residue())
wide_helix_turn.members[5].Obonding = [4]
wide_helix_turn.members[5].Hbonding = [-4]
wide_helix_turn.members[5].move = 1
wide_helix_turn.members[5].index = 'f'
wide_helix_turn.members.append(member_residue())
wide_helix_turn.members[6].Obonding = [4]
wide_helix_turn.members[6].Hbonding = [-5]
wide_helix_turn.members[6].move = None
wide_helix_turn.members[6].index = 'g'
#-------------------------------------------------------------------------------
#}}}

#{{{ helix caps
#-------------------------------------------------------------------------------
#Obonding pattern: 4,5,3
dans_c_cap = ssfingerprint('dans_c_cap')
#need to get the actual names for these from Dan
dans_c_cap.labelhash = {
'a':'dans_c_cap_reg_helix',
'b':'dans_c_cap_5_bonding',
'c':'dans_c_cap_3_bonding',
'd':'dans_c_cap_intervening',
'e':'dans_c_cap_last_alpha_bond',
'f':'dans_c_cap_3_bonded_cap',
'g':'dans_c_cap_5_bonded_cap'}
dans_c_cap.labellist = [
'dans_c_cap_reg_helix',
'dans_c_cap_5_bonding',
'dans_c_cap_3_bonding',
'dans_c_cap_intervening',
'dans_c_cap_last_alpha_bond',
'dans_c_cap_3_bonded_cap',
'dans_c_cap_5_bonded_cap']
dans_c_cap.members.append(member_residue())
dans_c_cap.members[0].Obonding = [4]
dans_c_cap.members[0].Hbonding = [-4]
dans_c_cap.members[0].index = 'a'
dans_c_cap.members[0].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[1].Obonding = [5,'any']
dans_c_cap.members[1].Hbonding = [-4]
dans_c_cap.members[1].index = 'b'
dans_c_cap.members[1].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[2].Obonding = [3]
dans_c_cap.members[2].Hbonding = [-4]
dans_c_cap.members[2].index = 'c'
dans_c_cap.members[2].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[3].Obanned = ['any']
dans_c_cap.members[3].index = 'd'
dans_c_cap.members[3].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[4].Obanned = ['any']
dans_c_cap.members[4].Hbonding = [-4]
dans_c_cap.members[4].index = 'e'
dans_c_cap.members[4].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[5].Hbonding = [-3,'any']
dans_c_cap.members[5].index = 'f'
dans_c_cap.members[5].move = 1
dans_c_cap.members.append(member_residue())
dans_c_cap.members[6].Hbonding = [-5]
dans_c_cap.members[6].index = 'g'
dans_c_cap.members[6].move = None
#-------------------------------------------------------------------------------
#}}}

#{{{ other helix
#-------------------------------------------------------------------------------
threeten_helix_single = ssfingerprint('threeten_helix_single')
threeten_helix_single.labelhash = {'a':'threeten_helix_single'}
threeten_helix_single.labellist = ['threeten_helix_single']
threeten_helix_single.members.append(member_residue())
threeten_helix_single.members[0].Obonding = [3]
threeten_helix_single.members[0].Hbonding = [-3]
threeten_helix_single.members[0].index = 'a'
threeten_helix_single.members[0].move = None

pi_helix_single = ssfingerprint('pi_helix_single')
pi_helix_single.labelhash = {'a':'pi_helix_single'}
pi_helix_single.labellist = ['pi_helix_single']
pi_helix_single.members.append(member_residue())
pi_helix_single.members[0].Obonding = [5]
pi_helix_single.members[0].Hbonding = [-5]
pi_helix_single.members[0].index = 'a'
pi_helix_single.members[0].move = None
#-------------------------------------------------------------------------------
#}}}

#{{{ beta definitions
#-------------------------------------------------------------------------------

#Two strands:
# g (h) i (j) k
# r (q) p (o) n
#Might want to include h,j,q,o in .labels at some point. Depends.
antiparallel_beta_close = ssfingerprint('antiparallel_beta_close')
antiparallel_beta_close.labelhash = {
'i':'antiparallel_beta_close',
'p':'antiparallel_beta_close'}
antiparallel_beta_close.labellist = ['antiparallel_beta_close']
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[0].Obonding = ['p']
antiparallel_beta_close.members[0].Hbonding = ['p']
antiparallel_beta_close.members[0].move = 'p'
antiparallel_beta_close.members[0].index = 'i'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[1].Obonding = ['i']
antiparallel_beta_close.members[1].Hbonding = ['i']
antiparallel_beta_close.members[1].move = 2
antiparallel_beta_close.members[1].index = 'p'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[2].Hbonding = ['g']
antiparallel_beta_close.members[2].move = 'g'
antiparallel_beta_close.members[2].index = 'r'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[3].Obonding = ['r']
antiparallel_beta_close.members[3].move = 4
antiparallel_beta_close.members[3].index = 'g'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[4].Hbonding = ['n']
antiparallel_beta_close.members[4].move = 'n'
antiparallel_beta_close.members[4].index = 'k'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[5].Obonding = ['k']
antiparallel_beta_close.members[5].move = 2
antiparallel_beta_close.members[5].index = 'n'
antiparallel_beta_close.members.append(member_residue())
antiparallel_beta_close.members[6].Obonding = ['i']
antiparallel_beta_close.members[6].Hbonding = ['i']
antiparallel_beta_close.members[6].move = None
antiparallel_beta_close.members[6].index = 'p'

#Two strands:
# (g) h i j (k)
# (r) q p o (n)
antiparallel_beta_wide = ssfingerprint('antiparallel_beta_wide')
antiparallel_beta_wide.labelhash = {
'i':'antiparallel_beta_wide',
'q':'antiparallel_beta_wide'}
antiparallel_beta_wide.labellist = ['antiparallel_beta_wide']
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[0].move = 1
antiparallel_beta_wide.members[0].index = 'i'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[1].Obonding = ['o']
antiparallel_beta_wide.members[1].Hbonding = ['o']
antiparallel_beta_wide.members[1].move = 'o'
antiparallel_beta_wide.members[1].index = 'j'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[2].Obonding = ['j']
antiparallel_beta_wide.members[2].Hbonding = ['j']
antiparallel_beta_wide.members[2].move = 1
antiparallel_beta_wide.members[2].index = 'o'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[3].move = 1
antiparallel_beta_wide.members[3].index = 'p'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[4].Obonding = ['h']
antiparallel_beta_wide.members[4].Hbonding = ['h']
antiparallel_beta_wide.members[4].move = 'h'
antiparallel_beta_wide.members[4].index = 'q'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[5].Obonding = ['q']
antiparallel_beta_wide.members[5].Hbonding = ['q']
antiparallel_beta_wide.members[5].move = 1
antiparallel_beta_wide.members[5].index = 'h'
antiparallel_beta_wide.members.append(member_residue())
antiparallel_beta_wide.members[6].move = None
antiparallel_beta_wide.members[6].index = 'i'

#Two strands:
# (g) h i j (k)
# (n) o p q (r)
parallel_beta = ssfingerprint('parallel_beta')
parallel_beta.labelhash = {
'i':'parallel_beta',
'q':'parallel_beta'}
parallel_beta.labellist = ['parallel_beta']
parallel_beta.members.append(member_residue())
parallel_beta.members[0].Obonding = ['q']
parallel_beta.members[0].Hbonding = ['o']
parallel_beta.members[0].move = 2
parallel_beta.members[0].index = 'i'
parallel_beta.members.append(member_residue())
parallel_beta.members[1].Hbonding = ['q']
parallel_beta.members[1].move = 'q'
parallel_beta.members[1].index = 'k'
parallel_beta.members.append(member_residue())
parallel_beta.members[2].Obonding = ['k']
parallel_beta.members[2].Hbonding = ['i']
parallel_beta.members[2].move = -1
parallel_beta.members[2].index = 'q'
parallel_beta.members.append(member_residue())
parallel_beta.members[3].move = -1
parallel_beta.members[3].index = 'p'
parallel_beta.members.append(member_residue())
parallel_beta.members[4].Obonding = ['i']
parallel_beta.members[4].Hbonding = ['g']
parallel_beta.members[4].move = 'g'
parallel_beta.members[4].index = 'o'
parallel_beta.members.append(member_residue())
parallel_beta.members[5].Obonding = ['o']
parallel_beta.members[5].move = 2
parallel_beta.members[5].index = 'g'
parallel_beta.members.append(member_residue())
#parallel_beta.members[6].Obonding = ['r']
#parallel_beta.members[6].Hbonding = ['p']
#Don't need to re-check the bonds if the index is right
parallel_beta.members[6].move = None
parallel_beta.members[6].index = 'i'
#-------------------------------------------------------------------------------
#}}}

#{{{ beta interrupts
#-------------------------------------------------------------------------------
bulge1 = ssfingerprint('bulge1')
#Needs a better name - is this classic (or will it be once properly defined?)
bulge1.labelhash = {
'p':'bulge1_bifur',
'i':'bulge1_distort',
'j':'bulge1_other'}
bulge1.labellist = ['bulge1_bifur','bulge1_distort','bulge1_other']
bulge1.members.append(member_residue())
#bulge1.members[0].Obonding = []
bulge1.members[0].Hbonding = ['p']
bulge1.members[0].move = 1
bulge1.members[0].index = 'i'
bulge1.members.append(member_residue())
bulge1.members[1].Obonding = ['p']
bulge1.members[1].Hbonding = ['p']
bulge1.members[1].move = 'p'
bulge1.members[1].index = 'j'
bulge1.members.append(member_residue())
bulge1.members[2].Obonding = ['i','j']
bulge1.members[2].Hbonding = ['j']
bulge1.members[2].move = None
bulge1.members[2].index = 'p'
#-------------------------------------------------------------------------------
#}}}

####-------------------------------------------------------------------------------
####{{{Single-sided beta structure (commented out at the moment)
####antiparallel_beta_close_oneside = ssfingerprint('antiparallel_beta_close_oneside')
####antiparallel_beta_close_oneside.labels = {'antiparallel_beta_close_oneside':['i','p']}
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[0].Obonding = ['p']
####antiparallel_beta_close_oneside.members[0].Hbonding = ['p']
####antiparallel_beta_close_oneside.members[0].move = 'p'
####antiparallel_beta_close_oneside.members[0].index = 'i'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[1].Obonding = ['i']
####antiparallel_beta_close_oneside.members[1].Hbonding = ['i']
####antiparallel_beta_close_oneside.members[1].move = 1
####antiparallel_beta_close_oneside.members[1].index = 'p'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[2].Obanned = ['any']
####antiparallel_beta_close_oneside.members[2].Hbanned = ['any']
####antiparallel_beta_close_oneside.members[2].move = 1
####antiparallel_beta_close_oneside.members[2].index = 'q'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[3].Hbonding = ['g']
####antiparallel_beta_close_oneside.members[3].move = 'g'
####antiparallel_beta_close_oneside.members[3].index = 'r'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[4].Obonding = ['r']
####antiparallel_beta_close_oneside.members[4].move = 1
####antiparallel_beta_close_oneside.members[4].index = 'g'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[5].Obanned = ['any']
####antiparallel_beta_close_oneside.members[5].Hbanned = ['any']
####antiparallel_beta_close_oneside.members[5].move = 2
####antiparallel_beta_close_oneside.members[5].index = 'h'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[6].Obanned = ['any']
####antiparallel_beta_close_oneside.members[6].Hbanned = ['any']
####antiparallel_beta_close_oneside.members[6].move = 1
####antiparallel_beta_close_oneside.members[6].index = 'j'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[7].Hbonding = ['n']
####antiparallel_beta_close_oneside.members[7].move = 'n'
####antiparallel_beta_close_oneside.members[7].index = 'k'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[8].Obonding = ['k']
####antiparallel_beta_close_oneside.members[8].move = 1
####antiparallel_beta_close_oneside.members[8].index = 'n'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[9].Obanned = ['any']
####antiparallel_beta_close_oneside.members[9].Hbanned = ['any']
####antiparallel_beta_close_oneside.members[9].move = 1
####antiparallel_beta_close_oneside.members[9].index = 'o'
####antiparallel_beta_close_oneside.members.append(member_residue())
####antiparallel_beta_close_oneside.members[10].Obonding = ['i']
####antiparallel_beta_close_oneside.members[10].Hbonding = ['i']
####antiparallel_beta_close_oneside.members[10].move = None
####antiparallel_beta_close_oneside.members[10].index = 'p'
#####Beats me if I've got this right yet . . .
#####  ghijk
#####  nopqr
####
####antiparallel_beta_wide_oneside = ssfingerprint('antiparallel_beta_wide_oneside')
####antiparallel_beta_wide_oneside.labels = {'antiparallel_beta_wide_oneside':['i','q']}
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[0].move = 1
####antiparallel_beta_wide_oneside.members[0].index = 'i'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[1].Obonding = ['p']
####antiparallel_beta_wide_oneside.members[1].Hbonding = ['p']
####antiparallel_beta_wide_oneside.members[1].move = 'p'
####antiparallel_beta_wide_oneside.members[1].index = 'j'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[2].Obonding = ['j']
####antiparallel_beta_wide_oneside.members[2].Hbonding = ['j']
####antiparallel_beta_wide_oneside.members[2].move = 1
####antiparallel_beta_wide_oneside.members[2].index = 'p'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[3].move = 1
####antiparallel_beta_wide_oneside.members[3].index = 'q'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[4].Obonding = ['h']
####antiparallel_beta_wide_oneside.members[4].Hbonding = ['h']
####antiparallel_beta_wide_oneside.members[4].move = 'h'
####antiparallel_beta_wide_oneside.members[4].index = 'r'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[5].Obonding = ['r']
####antiparallel_beta_wide_oneside.members[5].Hbonding = ['r']
####antiparallel_beta_wide_oneside.members[5].move = 1
####antiparallel_beta_wide_oneside.members[5].index = 'h'
####antiparallel_beta_wide_oneside.members.append(member_residue())
####antiparallel_beta_wide_oneside.members[6].move = None
####antiparallel_beta_wide_oneside.members[6].index = 'i'
####
####parallel_beta_oneside = ssfingerprint('parallel_beta_oneside')
####parallel_beta_oneside.labels = {'parallel_beta_oneside':['i','q']}
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[0].Obonding = ['r']
####parallel_beta_oneside.members[0].Hbonding = ['p']
####parallel_beta_oneside.members[0].move = 2
####parallel_beta_oneside.members[0].index = 'i'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[1].Hbonding = ['r']
####parallel_beta_oneside.members[1].move = 'r'
####parallel_beta_oneside.members[1].index = 'k'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[2].Obonding = ['k']
####parallel_beta_oneside.members[2].Hbonding = ['i']
####parallel_beta_oneside.members[2].move = -1
####parallel_beta_oneside.members[2].index = 'r'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[3].move = -1
####parallel_beta_oneside.members[3].index = 'q'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[4].Obonding = ['i']
####parallel_beta_oneside.members[4].Hbonding = ['g']
####parallel_beta_oneside.members[4].move = 'g'
####parallel_beta_oneside.members[4].index = 'p'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[5].Obonding = ['p']
####parallel_beta_oneside.members[5].move = 2
####parallel_beta_oneside.members[5].index = 'g'
####parallel_beta_oneside.members.append(member_residue())
####parallel_beta_oneside.members[6].move = None
####parallel_beta_oneside.members[6].index = 'i'

#}}}
####-------------------------------------------------------------------------------

#####{{{Single-residue helix transitions (commented out)
####h310_to_alpha = ssfingerprint('h310_to_alpha')
####h310_to_alpha.labels = {'h310_to_alpha':['a']}
####h310_to_alpha.members.append(member_residue())
####h310_to_alpha.members[0].Obonding = [4]
####h310_to_alpha.members[0].Hbonding = [-3]
####h310_to_alpha.members[0].move = None
####h310_to_alpha.members[0].index = 'a'
####
####alpha_to_310 = ssfingerprint('alpha_to_310')
####alpha_to_310.labels = {'alpha_to_310':['a']}
####alpha_to_310.members.append(member_residue())
####alpha_to_310.members[0].Obonding = [3]
####alpha_to_310.members[0].Hbonding = [-4]
####alpha_to_310.members[0].move = None
####alpha_to_310.members[0].index = 'a'
####
####pi_to_alpha = ssfingerprint('pi_to_alpha')
####pi_to_alpha.labels = {'pi_to_alpha':['a']}
####pi_to_alpha.members.append(member_residue())
####pi_to_alpha.members[0].Obonding = [4]
####pi_to_alpha.members[0].Hbonding = [-5]
####pi_to_alpha.members[0].move = None
####pi_to_alpha.members[0].index = 'a'
####
####alpha_to_pi = ssfingerprint('alpha_to_pi')
####alpha_to_pi.labels = {'alpha_to_pi':['a']}
####alpha_to_pi.members.append(member_residue())
####alpha_to_pi.members[0].Obonding = [5]
####alpha_to_pi.members[0].Hbonding = [-4]
####alpha_to_pi.members[0].move = None
####alpha_to_pi.members[0].index = 'a'
#####}}}

#{{{ Other fingerprints
#-------------------------------------------------------------------------------
narrow_3 = ssfingerprint('narrow_3')
#A close-beta-type bonding pair 3 residues apart
narrow_3.labelhash = {
'a':'narrow_3_nsidebonding',
'b':'narrow_3_loop1',
'c':'narrow_3_loop2',
'd':'narrow_3_csidebonding'}
narrow_3.labellist = [
'narrow_3_nsidebonding',
'narrow_3_loop1',
'narrow_3_loop2',
'narrow_3_csidebonding']
narrow_3.members.append(member_residue())
narrow_3.members[0].Obonding = [3]
narrow_3.members[0].Hbonding = [3]
narrow_3.members[0].index = 'a'
narrow_3.members[0].move = 1
narrow_3.members.append(member_residue())
narrow_3.members[1].index = 'b'
narrow_3.members[1].move = 1
narrow_3.members.append(member_residue())
narrow_3.members[2].index = 'c'
narrow_3.members[2].move = 1
narrow_3.members.append(member_residue())
narrow_3.members[3].Obonding = [-3]
narrow_3.members[3].Hbonding = [-3]
narrow_3.members[3].index = 'd'
narrow_3.members[3].move = None

narrow_4 = ssfingerprint('narrow_4')
#A close-beta-type bonding pair 4 residues apart
narrow_4.labelhash = {
'a':'narrow_4_nsidebonding',
'b':'narrow_4_loop1',
'c':'narrow_4_loop2',
'd':'narrow_4_loop3',
'e':'narrow_4_csidebonding'}
narrow_4.labellist = [
'narrow_4_nsidebonding',
'narrow_4_loop1',
'narrow_4_loop2',
'narrow_4_loop3',
'narrow_4_csidebonding']
narrow_4.members.append(member_residue())
narrow_4.members[0].Obonding = [4]
narrow_4.members[0].Hbonding = [4]
narrow_4.members[0].index = 'a'
narrow_4.members[0].move = 1
narrow_4.members.append(member_residue())
narrow_4.members[1].index = 'b'
narrow_4.members[1].move = 1
narrow_4.members.append(member_residue())
narrow_4.members[2].index = 'c'
narrow_4.members[2].move = 1
narrow_4.members.append(member_residue())
narrow_4.members[3].index = 'd'
narrow_4.members[3].move = 1
narrow_4.members.append(member_residue())
narrow_4.members[4].Obonding = [-4]
narrow_4.members[4].Hbonding = [-4]
narrow_4.members[4].index = 'e'
narrow_4.members[4].move = None

narrow_5 = ssfingerprint('narrow_5')
#A close-beta-type bonding pair 5 residues apart
narrow_5.labelhash = {
'a':'narrow_5_nsidebonding',
'b':'narrow_5_loop1',
'c':'narrow_5_loop2',
'd':'narrow_5_loop3',
'e':'narrow_5_loop4',
'f':'narrow_5_csidebonding'}
narrow_5.labellist = [
'narrow_5_nsidebonding',
'narrow_5_loop1',
'narrow_5_loop2',
'narrow_5_loop3',
'narrow_5_loop4',
'narrow_5_csidebonding']
narrow_5.members.append(member_residue())
narrow_5.members[0].Obonding = [5]
narrow_5.members[0].Hbonding = [5]
narrow_5.members[0].index = 'a'
narrow_5.members[0].move = 1
narrow_5.members.append(member_residue())
narrow_5.members[1].index = 'b'
narrow_5.members[1].move = 1
narrow_5.members.append(member_residue())
narrow_5.members[2].index = 'c'
narrow_5.members[2].move = 1
narrow_5.members.append(member_residue())
narrow_5.members[3].index = 'd'
narrow_5.members[3].move = 1
narrow_5.members.append(member_residue())
narrow_5.members[4].index = 'e'
narrow_5.members[4].move = 1
narrow_5.members.append(member_residue())
narrow_5.members[5].Obonding = [-5]
narrow_5.members[5].Hbonding = [-5]
narrow_5.members[5].index = 'f'
narrow_5.members[5].move = None
#-------------------------------------------------------------------------------
#}}}


#Isolated ribbons (that is, just a pair of strands) often have very pronounced
#  curl that might help us get a handle on how cablam space shows curl
#With Hbanned and Obanned, I can ensure one-sided beta motifs
#Also: the overlap of residues that have wide *and* close beta motifs would be interesting

fingerprints = {'alpha_helix_3':alpha_helix_3,
'alpha_helix_single':alpha_helix_single, 'threeten_helix_single':threeten_helix_single, 'pi_helix_single':pi_helix_single,
'antiparallel_beta_close':antiparallel_beta_close,'antiparallel_beta_wide':antiparallel_beta_wide,
'parallel_beta':parallel_beta,
'wide_helix_turn':wide_helix_turn,
'reg_alpha':reg_alpha,
'alpha_turn':alpha_turn,
#'h310_to_alpha':h310_to_alpha,'alpha_to_310':alpha_to_310,
#'pi_to_alpha':pi_to_alpha,'alpha_to_pi':alpha_to_pi,
'bulge1':bulge1,
'narrow_3':narrow_3,'narrow_4':narrow_4,'narrow_5':narrow_5,
'dans_c_cap':dans_c_cap}

def get_all_labels(fingerprint_list):
  #Given a list of fingerprints of interest (args.probe.split(',') from
  #  cablam_training for example, this function returns a list of all the
  #  residue labels used by all those fingerprints.  That label_list is used for
  #  printing output to kin or csv.)
  label_list = []
  for fingerprint_name in fingerprint_list:
    if fingerprint_name in fingerprints:
      for label in fingerprints[fingerprint_name].labellist:
        label_list.append(label)
    else:
      continue
  return label_list

def annote_motif_residue(residue, fingerprint_name):
  #The idea here is to check a residue against a fingerprint
  #probably checking against the i residue, at least for a start
  #residue objects contain links to sequence-adjacent residues, but space-adjacent ones will be a problem
  #
  #Found motif information will be stored in residue.motifs for the time being
  fingerprint = fingerprints[fingerprint_name]
  fingerprint.checkall(residue)

def annote_motif_protein(resdata, fingerprint_name):
  #Wraps the above 'annote_motif_residue' function for a whole protein
  reslist = resdata.keys()
  reslist.sort()
  for resid in reslist:
    residue = resdata[resid]
    annote_motif_residue(residue, fingerprint_name)

###cablam_access.py
###def make_cablam_object(hierarchy):
###  prune HOH?
###  return resdata
###
###def add_geometry_to_cablam(cablam_object, desired_measures):
###  return resdata
###
###def refresh_probe_data(cablam_object, probe_data?):
###  return resdata
###
###def force_connected_residues(res1,res2):
###  add/ensures nextres/prevres relationship
