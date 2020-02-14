from __future__ import absolute_import, division, print_function
# (jEdit options) :folding=explicit:collapseFolds=1:
from mmtbx.validation import residue, validation, atom
from cctbx import geometry_restraints
from mmtbx.rotamer.n_dim_table import NDimTable #handles contours
from libtbx import easy_pickle #NDimTables are stored as pickle files
import libtbx.load_env
import os, sys
from iotbx.pdb.hybrid_36 import hy36decode

#{{{ global constants
#-------------------------------------------------------------------------------
MAINCHAIN_ATOMS = [" N  "," C"  ," O  "," CA "]
CA_PSEUDOBOND_DISTANCE = 4.5
PEPTIDE_BOND_DISTANCE = 2.0
#2.0A is the peptide bond cutoff used by O for model building and potentially-
#  fragmented chains. O's generous cutoff seemed appropriate since I expect to
#  process in-progress models with this program
#RLab (probably) uses 1.4 +- 0.3, official textbook is about 1.33
#O's Ca-Ca distance is 4.5A
#Tom T suggests 2.5A as low-end cutoff for Ca-Ca distance
#--He also suggests using the pdb interpreter's internal chain break def

CABLAM_OUTLIER_CUTOFF = 0.01
CABLAM_DISFAVORED_CUTOFF = 0.05
CA_GEOM_CUTOFF = 0.005
#These cutoffs are used to identify outliers in CaBLAM parameter spaces
#These values were set heuristically through inspection of known outliers

ALPHA_CUTOFF = 0.001
BETA_CUTOFF = 0.0001
THREETEN_CUTOFF = 0.001
#These cutoffs are used to identify probable secondary structure residues
#Individual secondary structure residues are assembled into complete secondary
#  structure elements
#These values were set heuristically through inspection of known secondary
#  structure in low-resolution models
#-------------------------------------------------------------------------------
#}}}

#{{{ interpretation
#-------------------------------------------------------------------------------
def interpretation():
  #prints a brief guide to interpreting CaBLAM results
  sys.stderr.write("""
---------------- *** CaBLAM validation data interpretation *** -----------------

CaBLAM uses CA-trace geometry to validate protein backbone. Since the CA trace
  is more reliably modeled from low-resolution density than the full backbone
  trace, CaBLAM can identify modeling errors and intended secondary structure
  elements in low-resolution models.

Text:
  Text output is provided in colon-separated format. The columns are as follows:
    residue : A residue identifier
    outlier_type : If the residue is an outlier by CaBLAM metrics, the type will
      appear here. There are 3 types of outliers:
        CaBLAM Outlier is similar to Ramachandran 'Outlier' in severity
        CaBLAM Disfavored is similar to Ramachandran 'Allowed'
        CA Geom Outlier indicates a severe error in the CA trace
    contour_level : The CaBLAM contour percentile score for the residue.
      0.05 (5%) or lower is Disfavored.  0.01 (1%) or lower is Outlier.
    ca_contour_level : The CA geometry contour percentile score for the residue.
      Serves as a sanity check for CaBLAM's other metrics.
      0.005 (0.5%) or lower is Outlier.
    sec struc recommendation : Secondary structure identification for this
      residue. There are 3 possible identifications - try alpha helix,
      try beta sheet, and try three-ten.
    alpha score : The alpha helix contour percentile score for this residue.
      0.001 (0.1%) or better is probable helix.
    beta score : The beta strand contour percentile score for this residue.
      0.0001 (0.01%) or better is probable strand.
    three-ten score : The 3-10 helix contour percentile score for this residue.
      0.001 (0.1%) or better is probable helix.

Kinemage:
  Kinemage output is available for visual markup of structures in KiNG.
""")
#-------------------------------------------------------------------------------
#}}}

#{{{ math functions
#-------------------------------------------------------------------------------
def perptersect(a1, a2, b1):
  #Finds the line from a1 to a2, drops a perpendicular to it from b1, and returns
  #  the point of intersection.
  A = [a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]]
  #Find the slope of line A in each direction, A is in vector notation
  t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / ((A[0]**2)+(A[1]**2)+(A[2]**2))
  #Solve the parametric equations (dot of perpendiculars=0). . .
  b2 = [a1[0]+A[0]*t, a1[1]+A[1]*t, a1[2]+A[2]*t]
  # . . . and use the result to find the new point b2 on the line
  return b2

def calculate_mu(CA1,CA2,CA3,CA4):
  #dihedral calculation for CA trace
  if None in [CA1,CA2,CA3,CA4]:
    return None
  return geometry_restraints.dihedral(sites=[CA1.xyz,CA2.xyz,CA3.xyz,CA4.xyz],
    angle_ideal=180, weight=1).angle_model

def calculate_ca_virtual_angle(CA1,CA2,CA3):
  #angle calculation for CA trace
  if None in [CA1,CA2,CA3]:
    return None
  return geometry_restraints.angle(sites=[CA1.xyz,CA2.xyz,CA3.xyz],
    angle_ideal=120, weight=1).angle_model

def calculate_nu(CA1,CA2,CA3,O1,O2):
  #dihedral calculation for peptide plane orientations
  if None in [CA1,CA2,CA3,O1,O2]:
    return None
  X1 = perptersect(CA1.xyz,CA2.xyz,O1.xyz)
  X2 = perptersect(CA2.xyz,CA3.xyz,O2.xyz)
  return geometry_restraints.dihedral(sites=[O1.xyz,X1,X2,O2.xyz],
    angle_ideal=180, weight=1).angle_model

def calculate_omega(CA1,C1,N2,CA2):
  #dihedral calculation for peptide bond
  if None in [CA1,C1,N2,CA2]:
    return None
  return geometry_restraints.dihedral(sites=[CA1.xyz,C1.xyz,N2.xyz,CA2.xyz],
    angle_ideal=180, weight=1).angle_model
#-------------------------------------------------------------------------------
#}}}

#{{{ contour fetching
#-------------------------------------------------------------------------------
def fetch_peptide_expectations():
  #This function finds, unpickles, and returns N-Dim Tables of expected residue
  #  behavior for use in determining cablam outliers
  #The return object is a dict keyed by residue type:
  #  'general','gly','transpro','cispro'
  #This set of contours defines peptide behavior (mu_in,mu_out,nu)
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled

def fetch_ca_expectations():
  #This function finds, unpickles, and returns N-Dim Tables of expected residue
  #  behavior for use in determining ca geometry outliers
  #The return object is a dict keyed by residue type:
  #  'general','gly','transpro','cispro'
  #This set of contours defines CA trace quality (mu_in,mu_d_out,ca_virtual)
  categories = ['general','gly','transpro','cispro']
  unpickled = {}
  for category in categories:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.expected."+category+"_CA.pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for category "+
        category+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[category] = ndt
  return unpickled

def fetch_motif_contours():
  #This function finds, unpickles, and returns N-Dim Tables of secondary structure
  #  structure behavior for use in determining what an outlier residue might
  #  really be
  #The return object is a dict keyed by secondary structure type:
  #  'regular_beta','loose_alpha','loose_threeten'
  motifs = ['regular_beta','loose_alpha','loose_threeten']
  unpickled = {}
  for motif in motifs:
    picklefile = libtbx.env.find_in_repositories(
      relative_path=(
        "chem_data/cablam_data/cablam.8000.motif."+motif+".pickle"),
      test=os.path.isfile)
    if (picklefile is None):
      sys.stderr.write("\nCould not find a needed pickle file for motif "+
        motif+" in chem_data.\nExiting.\n")
      sys.exit()
    ndt = easy_pickle.load(file_name=picklefile)
    unpickled[motif] = ndt
  return unpickled
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam data storage classes
#-------------------------------------------------------------------------------
class cablam_geometry():
  #holds cablam's geometry parameters for one residue
  def __init__(self, mu_in=None, mu_out=None, nu=None, ca_virtual=None, omega=None):
    self.mu_in = mu_in
    self.mu_out = mu_out
    self.nu = nu
    self.ca_virtual = ca_virtual
    self.omega=omega

class cablam_score():
  #holds cablam contour scores for one residue
  def __init__(self,cablam=None,c_alpha_geom=None,alpha=None,beta=None,threeten=None):
    self.cablam = cablam
    self.c_alpha_geom = c_alpha_geom
    self.alpha = alpha
    self.beta = beta
    self.threeten = threeten

class cablam_feedback():
  #holds outliers status and secondary structure identifications for one residue
  def __init__(self):
    self.cablam_outlier = None
    self.cablam_disfavored = None
    self.c_alpha_geom_outlier = None
    self.alpha=None
    self.beta=None
    self.threeten=None

#cablam results are stored by chains and by conformers within chains, in
#  parallel with the conformer organization of the pdb hierarchy.  These classes
#  handle that organization of the results
class cablam_chain():
  #chain-level organization for cablam results
  def __init__(self):
    self.conf_names = []
    self.confs = {}

class cablam_conf():
  #conformer-level organization for cablam results
  def __init__(self):
    self.conf_name = None
    self.results = {}
    self.sec_struc_records = []

class secondary_structure_segment():
  #holds a secondary structure element identified by cablam
  def __init__(self, start, end, segment_type, segment_length):
    self.start = start
    self.end = end
    self.segment_type = segment_type
    self.segment_length = segment_length

class wheel_wedge():
  #a wedge of space generated from rotating a peptide plane, and the corresponding cablam score
  #for drawing kinemage markup
  def __init__(self, start, end, offset, cablam_score):
    self.start = [start[0]-offset[0], start[1]-offset[1], start[2]-offset[2]]
    self.end = [end[0]-offset[0], end[1]-offset[1], end[2]-offset[2]]
    self.cablam_score = cablam_score
#-------------------------------------------------------------------------------
#}}}

#{{{ cablam_result class
#-------------------------------------------------------------------------------
class cablam_result(residue):
  """
  Result class for cablam analysis
  """
  __slots__ = residue.__slots__ + [
    "residue",
    "prevres",
    "nextres",
    "_outlier_i_seqs",
    "has_ca",
    "has_mc",
    "has_any_alts",
    "has_mc_alts",
    "mc_alts",
    "alts",
    "measures",
    "scores",
    "feedback"
    ]

  @staticmethod
  def header():
    return "%-12s %-19s %-7s %-7s %-17s %-7s %-7s %-7s" % (
      'Residue','Outlier type','Cablam','CA geom',' Sec structure','Alpha','Strand','3-10')

  #{{{ as_string
  #-----------------------------------------------------------------------------
  def as_string(self):
    #print text output
    outlist = self.as_list_for_text(outliers_only=False)
    if outlist is None:
      return ""
    else:
      return "%s %s %s %s %s %s %s %s" % (
        outlist[0], outlist[1], outlist[2], outlist[3], outlist[4], outlist[5], outlist[6], outlist[7])
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ mp_id
  #-----------------------------------------------------------------------------
  def mp_id(self):
    #Returns an id consistent with MolProbity 'cnit' ids
    #Formatted as: ccnnnnilttt
    #  c: 2-char Chain ID, space for none
    #  n: sequence number, right justified, space padded
    #  i: insertion code, space for none
    #  l: alternate ID, space for none
    #  t: residue type (ALA, LYS, etc.), all caps left justified, space padded
    ##id_str = self.residue.id_str()
    ##| A  75 |
    ##chain = id_str[0:2]
    chain = self.chain_id
    resnum = self.residue.resseq
    ins = self.residue.icode
    resname = self.resname
    alt = self.altloc
    if alt == '':
      alt = ' '
    return chain + resnum + ins + alt + resname
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ sorting_id
  #-----------------------------------------------------------------------------
  def sorting_id(self):
    #Returns an id used for sorting residues
    #Formatted as: ccnnnni
    #  c: 2-char Chain ID, space for none
    #  n: sequence number, right justified, space padded
    #  i: insertion code, space for none
    chain = self.chain_id
    resnum = self.residue.resseq
    ins = self.residue.icode
    return chain + resnum + ins
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get_atom
  #-----------------------------------------------------------------------------
  #returns the desired atom from a hierarchy residue object
  def get_atom(self, atom_name):
    for atom in self.residue.atoms():
      if atom.name == atom_name: return atom
    else: return None
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ get_cablam_atoms
  #-----------------------------------------------------------------------------
  def get_cablam_atoms(self):
    #finds and returns all the atoms necessary for the CaBLAM calculations for
    #  this residue
    #returned atoms are:
    #res0_CA
    #res1_CA, res1_O, res1_C
    #res2_CA, res2_O, res2_N (res2 is the current residue)
    #res3_CA
    #res4_CA
    #returns None for missing/unfound atoms
    atom_set = {}
    atom_set["res2_CA"]=self.get_atom(" CA ")
    atom_set["res2_O"]= self.get_atom(" O  ")#for nu
    atom_set["res2_N"]= self.get_atom(" N  ")#for omega
    if self.prevres:
      atom_set["res1_CA"]=self.prevres.get_atom(" CA ")
      atom_set["res1_O"]= self.prevres.get_atom(" O  ")#for nu
      atom_set["res1_C"]= self.prevres.get_atom(" C  ")#for omega
      if self.prevres.prevres:
        atom_set["res0_CA"]=self.prevres.prevres.get_atom(" CA ")
      else:
        atom_set["res0_CA"]=None
    else:
      atom_set["res1_CA"]=None
      atom_set["res1_O"]= None
      atom_set["res1_C"]= None
      atom_set["res0_CA"]=None
    if self.nextres:
      atom_set["res3_CA"]=self.nextres.get_atom(" CA ")
      if self.nextres.nextres:
        atom_set["res4_CA"]=self.nextres.nextres.get_atom(" CA ")
      else:
        atom_set["res4_CA"]=None
    else:
      atom_set["res3_CA"]=None
      atom_set["res4_CA"]=None
    return atom_set
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ check_atoms
  #-----------------------------------------------------------------------------
  #checks whether atoms necessary for cablam calculations are present in this
  #  residue.
  #Sets has_ca = True if the residue has a CA atom (minimum cablam requirement)
  #Sets has_mc = True if the residue has all 4 mainchain heavy atoms
  def check_atoms(self):
    if self.get_atom(' CA ') is None: pass
    else:
      self.has_ca = True
      for atom_name in [' N  ',' C  ',' O  ']:
        if self.get_atom(atom_name) is None:
          break
      else:
        self.has_mc = True
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ link_residues
  #-----------------------------------------------------------------------------
  def link_residues(self, previous_result):
    #CaBLAM calculations depend on traversing sequential residues both backwards
    #  and forwards in sequence.  This function creates "links" between
    #  sequential residues.
    #Links are established if the residues are within probable bonding distance.
    #prevres and nextres values default to None during cablam_data.__init__()
    if previous_result is None:
      #no previous residue to link to
      return #Default link is None
    elif not previous_result.has_ca or not self.has_ca:
      #previous residue's proximity cannot be checked
      return
    elif not previous_result.has_mc or not self.has_mc:
      #CA-trace-only: use CA to check proximity
      ca1 = previous_result.get_atom(' CA ').xyz
      ca2 = self.get_atom(' CA ').xyz
      cadist = ((ca1[0]-ca2[0])**2 + (ca1[1]-ca2[1])**2 + (ca1[2]-ca2[2])**2)**0.5
      if cadist > CA_PSEUDOBOND_DISTANCE:
        #CA atoms are too far apart
        return
      else:
        #
        previous_result.nextres = self
        self.prevres = previous_result
    else: #has full set of mc heavy atoms available
      # use previous ' C  ' and current ' N  ' to check proximity
      c = previous_result.get_atom(' C  ').xyz
      n = self.get_atom(' N  ').xyz
      peptidedist = ((c[0]-n[0])**2 + (c[1]-n[1])**2 + (c[2]-n[2])**2)**0.5
      if peptidedist > PEPTIDE_BOND_DISTANCE:
        #atoms that would peptide bond are too far apart
        return
      else:
        previous_result.nextres = self
        self.prevres = previous_result
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ calculate_cablam_geomtery
  #-----------------------------------------------------------------------------
  def calculate_cablam_geometry(self):
    #populates self.measures with geometric measures relevant to CaBLAM
    #these measures are: mu_in, mu_out, nu, ca_virtual, omega
    atom_set = self.get_cablam_atoms()
    self.measures = cablam_geometry(
      mu_in  = calculate_mu(atom_set['res0_CA'],atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA']),
      mu_out = calculate_mu(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA'],atom_set['res4_CA']),
      nu = calculate_nu(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA'],atom_set['res1_O'],atom_set['res2_O']),
      ca_virtual = calculate_ca_virtual_angle(atom_set['res1_CA'],atom_set['res2_CA'],atom_set['res3_CA']),
      omega = calculate_omega(atom_set['res1_CA'],atom_set['res1_C'],atom_set['res2_N'],atom_set['res2_CA'])
      )
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ contour_category
  #-----------------------------------------------------------------------------
  def contour_category(self):
    #determines the category of the current residue so that it can be paired
    #  with the correct contours
    #these categories are: 'general', 'gly', 'transpro', 'cispro'
    resname = self.resname.upper()
    if resname == "GLY": return "gly"
    elif resname == "PRO":
      if self.measures.omega is not None and self.measures.omega < 90 and self.measures.omega > -90:
        return "cispro"
      else: return "transpro"
    else: return "general"
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ calculate_contour_values
  #-----------------------------------------------------------------------------
  def calculate_contour_values(self, cablam_contours, ca_contours, motif_contours):
    #populates self.scores[alt] with contour values for the current residue
    #these contour values are: cablam, c_alpha_geom, alpha, beta, threeten
    category = self.contour_category()
    cablam_point = [self.measures. mu_in,self.measures.mu_out, self.measures.nu]
    ca_point = [self.measures.mu_in, self.measures.mu_out, self.measures.ca_virtual]
    motif_point = [self.measures.mu_in, self.measures.mu_out]
    if None in cablam_point: cablam=None
    else: cablam = cablam_contours[category].valueAt(cablam_point)
    if None in ca_point: c_alpha_geom=None
    else: c_alpha_geom = ca_contours[category].valueAt(ca_point)
    if None in motif_point: alpha, beta, threeten = 0, 0, 0
    else:
      alpha =    motif_contours['loose_alpha'].valueAt(motif_point)
      beta =     motif_contours['regular_beta'].valueAt(motif_point)
      threeten = motif_contours['loose_threeten'].valueAt(motif_point)
    self.scores = cablam_score(
      cablam=cablam,
      c_alpha_geom=c_alpha_geom,
      alpha=alpha,
      beta=beta,
      threeten=threeten)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ set_cablam_feedback
  #-----------------------------------------------------------------------------
  def set_cablam_feedback(self):
    #populates self.feedback with True/False outlier status for easy access
    #these statuses are: cablam_outlier, cablam_disfavored, c_alpha_geom_outlier
    self.feedback = cablam_feedback()
    #outlier flags
    if self.scores.cablam is not None and self.scores.cablam < CABLAM_OUTLIER_CUTOFF:
      self.feedback.cablam_outlier = True
    else:
      self.feedback.cablam_outlier = False
    if self.scores.cablam is not None and self.scores.cablam < CABLAM_DISFAVORED_CUTOFF:
      self.feedback.cablam_disfavored = True
    else:
      self.feedback.cablam_disfavored = False
    if self.scores.c_alpha_geom is not None and self.scores.c_alpha_geom < CA_GEOM_CUTOFF:
      self.feedback.c_alpha_geom_outlier = True
    else:
      self.feedback.c_alpha_geom_outlier = False
    if self.feedback.cablam_outlier or self.feedback.cablam_disfavored or self.feedback.c_alpha_geom_outlier:
      self.outlier = True
    #secondary structure
    #This is semi-duplicated from assemble_secondary_structure
    if not self.prevres or not self.nextres:
      #alpha, beta, and threeten defaults are already None
      return
    if ((self.scores.alpha >= ALPHA_CUTOFF or self.scores.threeten >=THREETEN_CUTOFF)
      and (self.prevres.scores.alpha >= ALPHA_CUTOFF or self.prevres.scores.threeten >= THREETEN_CUTOFF)
      and (self.nextres.scores.alpha >= ALPHA_CUTOFF or self.nextres.scores.threeten >= THREETEN_CUTOFF)):
      if (self.scores.threeten > self.scores.alpha and
        (self.nextres.scores.threeten > self.nextres.scores.alpha or
          self.prevres.scores.threeten > self.prevres.scores.alpha)):
        self.feedback.threeten=True
      else:
        self.feedback.alpha=True
    if self.scores.beta >= BETA_CUTOFF and self.prevres.scores.beta >= BETA_CUTOFF and self.nextres.scores.beta >= BETA_CUTOFF:
      self.feedback.beta=True
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ find_single_outlier_type, find_single_structure_suggestion
  #-----------------------------------------------------------------------------
  def find_single_outlier_type(self):
    if self.feedback.c_alpha_geom_outlier:
      outlier_type = ' CA Geom Outlier    '
    elif self.feedback.cablam_outlier:
      outlier_type = ' CaBLAM Outlier     '
    elif self.feedback.cablam_disfavored:
      outlier_type = ' CaBLAM Disfavored  '
    else:
      outlier_type = '                    '
    return outlier_type

  def find_single_structure_suggestion(self):
    if self.feedback.c_alpha_geom_outlier:
      suggestion = '                 '
    elif self.feedback.threeten:
      suggestion = ' try three-ten   '
    elif self.feedback.alpha:
      suggestion = ' try alpha helix '
    elif self.feedback.beta:
      suggestion = ' try beta sheet  '
    else:
      suggestion = '                 '
    return suggestion
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_kinemage
  #-----------------------------------------------------------------------------
  def as_kinemage(self, mode=None, out=sys.stdout):
    #prints kinemage markup for this residue
    #has separate output modes for cablam outliers and for ca geom outliers
    if mode == 'ca_geom':
      if self.feedback.c_alpha_geom_outlier is not None:
        stats = self.mp_id() + " ca_geom=%.2f alpha=%.2f beta=%.2f three-ten=%.2f" %(self.scores.c_alpha_geom*100, self.scores.alpha*100, self.scores.beta*100, self.scores.threeten*100)
        CA_1 = self.prevres.get_atom(' CA ').xyz
        CA_2 = self.get_atom(' CA ').xyz
        CA_3 = self.nextres.get_atom(' CA ').xyz
        out.write('\n{'+stats+'} P '+str(CA_2[0]-(CA_2[0]-CA_1[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_1[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_1[2])*0.9))
        out.write('\n{'+stats+'} '+str(CA_2[0])+' '+str(CA_2[1])+' '+str(CA_2[2]))
        out.write('\n{'+stats+'} '+str(CA_2[0]-(CA_2[0]-CA_3[0])*0.9)+' '+str(CA_2[1]-(CA_2[1]-CA_3[1])*0.9)+' '+str(CA_2[2]-(CA_2[2]-CA_3[2])*0.9))
    elif mode == 'cablam':
      if self.feedback.cablam_outlier is not None:
        stats = self.mp_id() + " cablam=%.2f alpha=%.2f beta=%.2f three-ten=%.2f" %(self.scores.cablam*100, self.scores.alpha*100, self.scores.beta*100, self.scores.threeten*100)
        CA_1, O_1 = self.prevres.get_atom(' CA ').xyz,self.prevres.get_atom(' O  ').xyz
        CA_2, O_2 = self.get_atom(' CA ').xyz,self.get_atom(' O  ').xyz
        CA_3      = self.nextres.get_atom(' CA ').xyz
        X_1 = perptersect(CA_1,CA_2,O_1)
        X_2 = perptersect(CA_2,CA_3,O_2)
        midpoint = [ (X_1[0]+X_2[0])/2.0 , (X_1[1]+X_2[1])/2.0 , (X_1[2]+X_2[2])/2.0 ]
        out.write('\n{'+stats+'} P '+ str(O_1[0]) +' '+ str(O_1[1]) +' '+ str(O_1[2]))
        out.write('\n{'+stats+'} '+ str(X_1[0]) +' '+ str(X_1[1]) +' '+ str(X_1[2]))
        out.write('\n{'+stats+'} '+ str(midpoint[0]) +' '+ str(midpoint[1]) +' '+ str(midpoint[2]))
        out.write('\n{'+stats+'} '+ str(X_2[0]) +' '+ str(X_2[1]) +' '+ str(X_2[2]))
        out.write('\n{'+stats+'} '+ str(O_2[0]) +' '+ str(O_2[1]) +' '+ str(O_2[2]))
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ calculate_kinemage_wheels
  #-----------------------------------------------------------------------------
  def calculate_kinemage_wheels(self, cablam_contours):
    from scitbx.matrix import rotate_point_around_axis
    category = self.contour_category()
    CA_1 = self.prevres.get_atom(' CA ').xyz
    O_1  = self.prevres.get_atom(' O  ').xyz
    CA_2 = self.get_atom(' CA ').xyz
    O_2  = self.get_atom(' O  ').xyz
    CA_3 = self.nextres.get_atom(' CA ').xyz
    if None in [CA_1, CA_2, CA_3, O_1, O_2]: return
    #markup wheels should extend to less than the full CO length
    #moving the used CO position is an easy way to propagate this across calculations
    #along the X-O line used in the dihedral so as not to change the cablam results
    scaling = 0.75
    X1 = perptersect(CA_1,CA_2,O_1)
    X2 = perptersect(CA_2,CA_3,O_2)
    O_1 = ( (O_1[0]-X1[0])*scaling+X1[0], (O_1[1]-X1[1])*scaling+X1[1], (O_1[2]-X1[2])*scaling+X1[2])
    O_2 = ( (O_2[0]-X2[0])*scaling+X2[0], (O_2[1]-X2[1])*scaling+X2[1], (O_2[2]-X2[2])*scaling+X2[2])
    #-----------------------------
    #markup wheels are offset from the carbonyl oxygen position for visual clarity
    #offset is a move along the CA-CA line
    #calculate unit vector alone CA-CA line, offset is some fraction of that
    offset = 0.15
    CA_2_1 = (CA_1[0]-CA_2[0], CA_1[1]-CA_2[1], CA_1[2]-CA_2[2])
    CA_2_1_len = (CA_2_1[0]**2 + CA_2_1[1]**2 + CA_2_1[2]**2)**0.5
    CA_2_1_offset = (CA_2_1[0]/CA_2_1_len*offset, CA_2_1[1]/CA_2_1_len*offset, CA_2_1[2]/CA_2_1_len*offset)
    CA_2_3 = (CA_3[0]-CA_2[0], CA_3[1]-CA_2[1], CA_3[2]-CA_2[2])
    CA_2_3_len = (CA_2_3[0]**2 + CA_2_3[1]**2 + CA_2_3[2]**2)**0.5
    CA_2_3_offset = (CA_2_3[0]/CA_2_3_len*offset, CA_2_3[1]/CA_2_3_len*offset, CA_2_3[2]/CA_2_3_len*offset)

    #Each CaBLAM outlier is based on the relative positions of *two* peptide planes, represented by CO positions
    #Test rotation of each peptide plane independently, and draw a wheel for each
    #Starting with O_2 is arbitrary, but O_2 is the CO in the residue named as an outlier by CaBLAM convention
    angle = -10
    #prev_O_2_xyz = O_2
    wheel1_center = (X2[0]-CA_2_3_offset[0], X2[1]-CA_2_3_offset[1], X2[2]-CA_2_3_offset[2])
    wheel1 = []
    while angle < 360:
      angle += 10
      new_xyz = rotate_point_around_axis(
        axis_point_1 = CA_2,
        axis_point_2 = CA_3,
        point        = O_2,
        angle        = angle,
        deg          = True)
      new_nu = geometry_restraints.dihedral(sites=[O_1, X1, X2, new_xyz],
                                            angle_ideal=180, weight=1).angle_model
      cablam_point = [self.measures.mu_in, self.measures.mu_out, new_nu]
      cablam_score = cablam_contours[category].valueAt(cablam_point)
      if cablam_score >= 0.05:
        wheel1.append(None)
        continue
      wedge_start = rotate_point_around_axis(
        axis_point_1 = CA_2,
        axis_point_2 = CA_3,
        point        = O_2,
        angle        = angle-5,
        deg          = True)
      wedge_end = rotate_point_around_axis(
        axis_point_1=CA_2,
        axis_point_2=CA_3,
        point=O_2,
        angle=angle+5,
        deg=True)
      wedge = wheel_wedge(wedge_start, wedge_end, CA_2_3_offset, cablam_score)
      wheel1.append(wedge)

    angle = -10
    #prev_O_1_xyz = O_1
    wheel2_center = (X1[0]-CA_2_1_offset[0], X1[1]-CA_2_1_offset[1], X1[2]-CA_2_1_offset[2])
    wheel2 = []
    while angle <= 360:
      angle += 10
      new_xyz = rotate_point_around_axis(
        axis_point_1 = CA_1,
        axis_point_2 = CA_2,
        point        = O_1,
        angle        = angle,
        deg          = True)
      new_nu = geometry_restraints.dihedral(sites=[new_xyz,X1,X2,O_2],
        angle_ideal=180, weight=1).angle_model
      cablam_point = [self.measures.mu_in,self.measures.mu_out, new_nu]
      cablam_score = cablam_contours[category].valueAt(cablam_point)
      if cablam_score >= 0.05:
        wheel2.append(None)
        continue
      wedge_start = rotate_point_around_axis(
        axis_point_1=CA_1,
        axis_point_2=CA_2,
        point=O_1,
        angle=angle-5,
        deg=True)
      wedge_end = rotate_point_around_axis(
        axis_point_1=CA_1,
        axis_point_2=CA_2,
        point=O_1,
        angle=angle+5,
        deg=True)
      wedge = wheel_wedge(wedge_start, wedge_end, CA_2_1_offset, cablam_score)
      wheel2.append(wedge)
    return [(wheel1, wheel1_center), (wheel2, wheel2_center)]
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_kinemage_point
  #-----------------------------------------------------------------------------
  def as_kinemage_point(self, out=sys.stdout):
    #printing for pointcloud kinemage output
    if not (self.measures.mu_in and self.measures.mu_out and self.measures.nu):
      return
    point_name = "{"+self.mp_id()+"}"
    print(point_name, "%.2f %.2f %.2f" % (self.measures.mu_in, self.measures.mu_out, self.measures.nu), file=out)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ is_same_as_other_result
  #-----------------------------------------------------------------------------
  def is_same_as_other_result(self,other_result):
    #Compare this result object to another to see if they are effectively the
    #  same.  Identical cablam geometry (mu_in, mu_out, and nu) is assumed to
    #  mean identical residues:
    #This method will probably change
    if (self.measures.mu_in != other_result.measures.mu_in
      or self.measures.mu_out != other_result.measures.mu_out
      or self.measures.nu != other_result.measures.nu):
      return False
    return True
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_list_for_text
  #-----------------------------------------------------------------------------
  def as_list_for_text(self, outliers_only=False):
    if not self.has_ca:
      return None
    if outliers_only:
      if not (self.feedback.cablam_disfavored or self.feedback.c_alpha_geom_outlier):
        return None

    if self.feedback.c_alpha_geom_outlier:
      outlier_type = ' CA Geom Outlier    '
    elif self.feedback.cablam_outlier:
      outlier_type = ' CaBLAM Outlier     '
    elif self.feedback.cablam_disfavored:
      outlier_type = ' CaBLAM Disfavored  '
    else:
      outlier_type = '                    '

    if self.scores.cablam is not None:
      cablam_level = '%.5f' %self.scores.cablam
    else:
      cablam_level = '       ' #default printing for CA-only models
    if self.scores.c_alpha_geom is not None:
      ca_geom_level = '%.5f' %self.scores.c_alpha_geom
    else:
      return None #if this is missing, there's nothing

    if self.feedback.c_alpha_geom_outlier:
      suggestion = '                 '
    elif self.feedback.threeten:
      suggestion = ' try three-ten   '
    elif self.feedback.alpha:
      suggestion = ' try alpha helix '
    elif self.feedback.beta:
      suggestion = ' try beta sheet  '
    else:
      suggestion = '                 '

    outlist = [self.mp_id() ,outlier_type, cablam_level, ca_geom_level, suggestion, '%.5f' %self.scores.alpha, '%.5f' %self.scores.beta, '%.5f' %self.scores.threeten]
    return outlist
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_table_row_phenix
  #-----------------------------------------------------------------------------
  def as_table_row_phenix(self):
    return [ self.chain_id,
      "%s %s" % (self.resname, self.resid),
      self.find_single_outlier_type().strip(),
      self.scores.cablam,
      self.scores.c_alpha_geom,
      self.find_single_structure_suggestion().strip(),
      self.scores.alpha,
      self.scores.beta,
      self.scores.threeten ]
  #-----------------------------------------------------------------------------
  #}}}
#-------------------------------------------------------------------------------
#}}}

#{{{ cablamalyze class
#-------------------------------------------------------------------------------
class cablamalyze(validation):
  """
  Frontend for calculating cablam statistics for a model
  """
  __slots__ = validation.__slots__ + [
    "residue_count",
    "outlier_count",
    "out",
    "pdb_hierarchy",
    "all_results",
    "summary_stats"
    ]

  program_description = "Analyze protein CA geometry for secondary structure identification - recommended for low-resolution structures"

  gui_list_headers = ["Chain","Residue","Evaluation","CaBLAM Score","CA Geometry Score","Secondary Structure","Helix Score","Beta Score","3-10 Score"]
  gui_formats = ["%s", "%s", "%s", "%.5f", "%.5f", "%s", "%.5f", "%.5f", "%.5f"]
  wx_column_widths = [125]*9

  def get_result_class(self): return cablam_result

  #{{{ __init__
  #-----------------------------------------------------------------------------
  def __init__(self,
    pdb_hierarchy,
      outliers_only,
      out,
      quiet,
      cablam_contours=None,
      ca_contours=None,
      motif_contours=None):
    validation.__init__(self)
    from mmtbx.validation import utils
    from scitbx.array_family import flex
    #self._outlier_i_seqs = flex.size_t()
    self.out = out
    if cablam_contours is None:
      cablam_contours = fetch_peptide_expectations()
    if ca_contours is None:
      ca_contours = fetch_ca_expectations()
    if motif_contours is None:
      motif_contours = fetch_motif_contours()
    self.pdb_hierarchy = pdb_hierarchy.deep_copy()
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if all_i_seqs.all_eq(0):
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)

    self.all_results = {}
    ordered_keys = {}
    all_keys = []
    confs = []
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if not chain.is_protein():
          continue
        for conf in chain.conformers():
          if conf.is_protein(): break #at least one conformer must be protein
        else: continue
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain).rjust(2)
        else:
          chain_id = chain.id.rjust(2)
        #The above .rjust(2)'s are to force 2-char chain ids
        current_chain = cablam_chain()
        self.all_results[chain_id] = current_chain
        previous_result = None
        for conf in chain.conformers():
          if not conf.is_protein():
            continue
          current_conf = cablam_conf()
          current_conf.conf_name = conf.altloc
          current_chain.confs[conf.altloc] = current_conf
          current_chain.conf_names.append(conf.altloc)
          for residue in conf.residues():
            result = cablam_result(
              residue=residue,
              resseq=residue.resseq,
              icode=residue.icode,
              altloc=conf.altloc,
              chain_id=chain_id,
              resname = residue.resname,
              #formatting note: residue.id_str() = 'pdbres="THR B 182 "'
              #chain=residue.id_str()[11:13],
              #residue.id_str() turned out to break on some segid formatting
              prevres=None,
              nextres=None,
              has_ca=False,
              has_mc=False,
              outlier = False,
              measures=None,
              scores=None
              )
            result.check_atoms()
            result.link_residues(previous_result)
            #Occasionally a conformer may have more than one "residue" that has
            # the same sorting_id (sorting_id is just chain+resseq+icode)
            # see phenix_regression/pdb/lysozyme_nohoh_plus6H.pdb, where WAT 14
            # and ARG 14 have the same sorting_id
            #This check my not be the correct behavior to catch such a
            # formatting error.
            if result.sorting_id() not in current_conf.results:
              current_conf.results[result.sorting_id()] = result
            previous_result = result
    for chain in self.all_results:
      for conf in self.all_results[chain].confs:
        for result in self.all_results[chain].confs[conf].results:
          if self.all_results[chain].confs[conf].results[result].has_ca:
            self.all_results[chain].confs[conf].results[result].calculate_cablam_geometry()
            self.all_results[chain].confs[conf].results[result].calculate_contour_values(cablam_contours, ca_contours, motif_contours)
            self.all_results[chain].confs[conf].results[result].xyz = self.all_results[chain].confs[conf].results[result].get_atom(' CA ').xyz
    self.outlier_count = 0
    for chain in self.all_results:
      for conf in self.all_results[chain].confs:
        for result in self.all_results[chain].confs[conf].results:
          if self.all_results[chain].confs[conf].results[result].has_ca:
            self.all_results[chain].confs[conf].results[result].set_cablam_feedback()
            if self.all_results[chain].confs[conf].results[result].outlier:
              self.outlier_count += 1
    self.assemble_secondary_structure()
    self.make_single_results_object(confs, all_keys)
    self.residue_count = len(self.results)
    self.summary_stats = self.make_summary_stats()
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ assemble_secondary_structure
  #-----------------------------------------------------------------------------
  def assemble_secondary_structure(self):
    #assembles complete secondary structure elements (alpha helices, three-ten
    #  helices, and beta strands) from individual residue
    #have to assemble from scores to handle helix transitions
    for chain in self.all_results:
      for conf_id in self.all_results[chain].confs:
        conf = self.all_results[chain].confs[conf_id]
        records = []
        record_start = None
        helix_in_progress = False
        result_ids = list(conf.results.keys())
        #result_ids.sort(key=lambda k: (k[0:2], int(hy36decode(4,k[2:6])), k[6:7])) #this broke for non 2-char segids
        result_ids.sort(key=lambda k: (conf.results[k].chain_id, int(hy36decode(len(conf.results[k].resseq),conf.results[k].resseq)), conf.results[k].icode))
        for result_id in result_ids:
          result = conf.results[result_id]
          #is it evaluable?
          if not result.prevres:
            continue
          if not result.has_ca or not result.nextres:
            if helix_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
              helix_in_progress = False
            continue

          #helix building
          #This requires that the residues be the center of three in any combination of helix types
          #threeten segments of only 2 are lost in this method relative to the previous
          if ((result.scores.alpha >= ALPHA_CUTOFF or result.scores.threeten >=THREETEN_CUTOFF)
            and (result.prevres.scores.alpha >= ALPHA_CUTOFF or result.prevres.scores.threeten >= THREETEN_CUTOFF)
            and (result.nextres.scores.alpha >= ALPHA_CUTOFF or result.nextres.scores.threeten >= THREETEN_CUTOFF)):
            #now determine which helix type the current residue should be identified as
            #if at least two residues in a row have higher threeten scores than alpha scores, they can be considered threeten
            if (result.scores.threeten > result.scores.alpha and
              (result.nextres.scores.threeten > result.nextres.scores.alpha or
                result.prevres.scores.threeten > result.prevres.scores.alpha)):
              thisres = 'threeten'
            else:
              thisres = 'alpha'
            if helix_in_progress:
              if thisres == helix_in_progress: #is it same as previous residue
                record_length += 1
                continue
              else: #or has it changed helix types
                records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
                helix_in_progress = thisres
                record_start = result
                record_length = 1
            else:
              helix_in_progress = thisres
              record_start = result
              record_length = 1
          else: #(current residue is not helix)
            if helix_in_progress:
              #might fail on chain breaks?
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type=helix_in_progress, segment_length=record_length))
              helix_in_progress = False
              record_start = None
              record_length = 0
            else:
              continue
        #helix building end

        #beta strands require another, separate pass
        strand_in_progress = False
        record_start = None
        for result_id in result_ids:
          result = conf.results[result_id]
          if not result.prevres:
            continue
          if not result.has_ca or not result.nextres:
            if strand_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type='beta', segment_length=record_length))
              strand_in_progress = False
            continue
          if result.scores.beta >= BETA_CUTOFF and result.prevres.scores.beta >= BETA_CUTOFF and result.nextres.scores.beta >= BETA_CUTOFF:
            if strand_in_progress:
              record_length += 1
              continue
            else:
              strand_in_progress = True
              record_start = result
              record_length = 1
          else:
            if strand_in_progress:
              records.append(secondary_structure_segment(start=record_start, end=result.prevres, segment_type='beta', segment_length=record_length))
              strand_in_progress = False
              record_start = None
              record_length = 0
            else:
              continue
        #beta strand building end
        #NOTE: Each strand is currently treated as its own sheet
        #Developing or implementing strand-to-sheet assembly that does not rely on
        #  H-bonds is a major future goal
        conf.sec_struc_records = records
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ make_single_results_object
  #-----------------------------------------------------------------------------
  def make_single_results_object(self, confs, all_keys):
    #should work without any arguments
    #populates self.results
    self.results = []
    chains = list(self.all_results.keys())
    chains.sort()
    for chain_id in chains:
      chain = self.all_results[chain_id]
      #take the first conformer as the basis for comparison
      conf  = chain.confs[chain.conf_names[0]]
      if len(chain.conf_names) == 0:
        for result_id in conf.results:
          result = conf.results[result_id]
          result.altloc = ''
          #set self.results id
        continue #go to next chain
      #else, combine results into single list
      result_ids = list(conf.results.keys())
      ###result_ids.sort()
      #for result_id in conf.results:
      for result_id in result_ids:
        result = conf.results[result_id]
        if not result.has_ca: continue
        #results without CAs have measures=None and break the
        #  is_same_as_other_result check. Also, they aren't evaluable residues.
        self.results.append(result)
        found_meaningful_alt = False
        for other_conf in chain.conf_names[1:]:
          if result.sorting_id() in chain.confs[other_conf].results:
            other_result = chain.confs[other_conf].results[result.sorting_id()]
            if not other_result.has_ca: continue
            if not result.is_same_as_other_result(other_result):
              self.results.append(other_result)
              found_meaningful_alt = True
        if not found_meaningful_alt:
          result.altloc = ''
          #set self.results id
          pass
    self.results.sort(key=lambda r: (r.chain_id, int(hy36decode(len(r.resseq),r.resseq)), r.icode, r.altloc))
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_text
  #-----------------------------------------------------------------------------
  def as_text(self, outliers_only=False):
    #prints colon-separated text for CaBLAM validation
    #one line per residue (or alternate)
    #Output is formatted to be human-readable, and is also used by MolProbity
    #This is the default output for running this script from commandline
    self.out.write('residue : outlier_type : contour_level : ca_contour_level : sec struc recommendation : alpha score : beta score : three-ten score')
    for result in self.results:
      if not result.has_ca:
        continue
      #if not result.feedback:
      #  continue
      if outliers_only:
        if not (result.feedback.cablam_disfavored or result.feedback.c_alpha_geom_outlier):
          continue

      if result.feedback.c_alpha_geom_outlier:
        outlier_type = ' CA Geom Outlier    '
      elif result.feedback.cablam_outlier:
        outlier_type = ' CaBLAM Outlier     '
      elif result.feedback.cablam_disfavored:
        outlier_type = ' CaBLAM Disfavored  '
      else:
        outlier_type = '                    '

      if result.scores.cablam is not None:
        cablam_level = '%.5f' %result.scores.cablam
      else:
        cablam_level = '       ' #default printing for CA-only models
      if result.scores.c_alpha_geom is not None:
        ca_geom_level = '%.5f' %result.scores.c_alpha_geom
      else:
        continue #if this is missing, there's nothing

      if result.feedback.c_alpha_geom_outlier:
        suggestion = '                 '
      elif result.feedback.threeten:
        suggestion = ' try three-ten   '
      elif result.feedback.alpha:
        suggestion = ' try alpha helix '
      elif result.feedback.beta:
        suggestion = ' try beta sheet  '
      else:
        suggestion = '                 '

      outlist = [result.mp_id() ,outlier_type, cablam_level, ca_geom_level, suggestion, '%.5f' %result.scores.alpha, '%.5f' %result.scores.beta, '%.5f' %result.scores.threeten]
      self.out.write('\n'+':'.join(outlist))
    self.out.write('\n')
    self.show_summary(out=self.out,prefix="")
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ make_summary_stats
  #-----------------------------------------------------------------------------
  def make_summary_stats(self):
    #calculates whole-model stats used by show_summary, gui_summary, and
    #  as_oneline
    residue_count = 0
    ca_residue_count = 0
    cablam_outliers = 0
    cablam_disfavored = 0
    ca_geom_outliers = 0
    prev_result_id = None
    is_residue = 0
    is_cablam_outlier = 0
    is_cablam_disfavored = 0
    is_ca_geom_outlier = 0
    alpha_count, beta_count, threeten_count = 0,0,0
    is_alpha, is_beta, is_threeten = 0,0,0
    correctable_helix_count, correctable_beta_count = 0,0
    is_correctable_helix, is_correctable_beta = 0,0
    for result in self.results:
      if not result.has_ca:
        continue
      is_residue = 1
      if result.scores.cablam is None:
        is_residue = 0
      if result.sorting_id() != prev_result_id:
        #new residue; update counts
        residue_count    += is_residue
        if result.scores.c_alpha_geom is not None:
          ca_residue_count += 1
        cablam_outliers  += is_cablam_outlier
        cablam_disfavored+= is_cablam_disfavored
        ca_geom_outliers += is_ca_geom_outlier
        is_cablam_outlier    = 0
        is_cablam_disfavored = 0
        is_ca_geom_outlier   = 0
        alpha_count    += is_alpha
        beta_count     += is_beta
        threeten_count += is_threeten
        is_alpha, is_beta, is_threeten = 0,0,0
        correctable_helix_count += is_correctable_helix
        correctable_beta_count += is_correctable_beta
        is_correctable_helix, is_correctable_beta = 0,0
      if result.scores.cablam is not None and result.scores.cablam < CABLAM_OUTLIER_CUTOFF and is_residue:
        is_cablam_outlier    = 1
        if result.feedback.alpha or result.feedback.threeten:
          is_correctable_helix = 1
        elif result.feedback.beta:
          is_correctable_beta = 1
      if result.scores.cablam is not None and result.scores.cablam < CABLAM_DISFAVORED_CUTOFF and is_residue:
        is_cablam_disfavored = 1
      if result.scores.c_alpha_geom is not None and result.scores.c_alpha_geom < CA_GEOM_CUTOFF:
        is_ca_geom_outlier = 1
      #---parse secondary structure---
      if result.feedback.alpha and not is_beta and not is_threeten:
        is_alpha = 1
      elif result.feedback.beta and not is_alpha and not is_threeten:
        is_beta = 1
      elif result.feedback.threeten and not is_alpha and not is_beta:
        is_threeten = 1
      #---endparse secondary structure---
      prev_result_id = result.sorting_id()
    residue_count    += is_residue
    cablam_outliers  += is_cablam_outlier
    cablam_disfavored+= is_cablam_disfavored
    ca_geom_outliers += is_ca_geom_outlier
    alpha_count    += is_alpha
    beta_count     += is_beta
    threeten_count += is_threeten
    return {'residue_count':residue_count,'ca_residue_count':ca_residue_count,
      'cablam_outliers':cablam_outliers,'cablam_disfavored':cablam_disfavored,'ca_geom_outliers':ca_geom_outliers,
      'alpha_count':alpha_count, 'beta_count':beta_count,'threeten_count':threeten_count,
      'correctable_helix_count':correctable_helix_count,'correctable_beta_count':correctable_beta_count}
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_oneline
  #-----------------------------------------------------------------------------
  def as_oneline(self,pdbid='pdbid'):
    #prints a one-line summary of cablam statistics for a structure
    #for oneline purposes, alternates are collapsed: each residue contributes up
    #  to 1 to each outlier count, regarless of how many outlier alternates it
    #  may contain
    if self.count_residues() == 0:
      self.out.write(pdbid+':0:0:0:0\n')
    else:
      self.out.write('%s:%i:%.1f:%.1f:%.2f\n' %(pdbid,self.count_residues(), self.percent_outliers(), self.percent_disfavored(), self.percent_ca_outliers()) )
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ Summary retrieval functions
  #-----------------------------------------------------------------------------
  def count_residues(self):
    return self.summary_stats['residue_count']
  def count_ca_residues(self):
    return self.summary_stats['ca_residue_count']
  def count_outliers(self):
    return self.summary_stats['cablam_outliers']
  def percent_outliers(self):
    if self.count_residues() == 0:
      return 0
    return self.count_outliers()/self.count_residues()*100
  def count_disfavored(self):
    return self.summary_stats['cablam_disfavored']
  def percent_disfavored(self):
    if self.count_residues() == 0:
      return 0
    return self.count_disfavored()/self.count_residues()*100
  def count_ca_outliers(self):
    return self.summary_stats['ca_geom_outliers']
  def percent_ca_outliers(self):
    if self.count_ca_residues() == 0:
      return 0
    return self.count_ca_outliers()/self.count_ca_residues()*100
  def count_helix(self):
    return self.summary_stats['alpha_count']+self.summary_stats['threeten_count']
  def count_beta(self):
    return self.summary_stats['beta_count']
  def percent_helix(self):
    if self.count_ca_residues() == 0:
      return 0
    return self.count_helix()/self.count_ca_residues()*100
  def percent_beta(self):
    if self.count_ca_residues() == 0:
      return 0
    return self.count_beta()/self.count_ca_residues()*100
  def count_correctable_helix(self):
    return self.summary_stats['correctable_helix_count']
  def count_correctable_beta(self):
    return self.summary_stats['correctable_beta_count']
  def percent_correctable_helix(self):
    if self.count_residues() == 0:
      return 0
    return self.count_correctable_helix()/self.count_residues()*100
  def percent_correctable_beta(self):
    if self.count_residues() == 0:
      return 0
    return self.count_correctable_beta()/self.count_residues()*100
  #-----------------------------------------------------------------------------
  #}}}

  def cablam_wheel_triangle(self, wheel_center, wedge, color):
    self.out.write('\n{} P X %s %.3f %.3f %.3f' % (color, wheel_center[0], wheel_center[1], wheel_center[2]))
    self.out.write('\n{} %s %.3f %.3f %.3f' % (color, wedge.start[0], wedge.start[1], wedge.start[2]))
    self.out.write('\n{} %s %.3f %.3f %.3f' % (color, wedge.end[0], wedge.end[1], wedge.end[2]))

  def cablam_wheel_edge(self, wheel_center, wedge, prevwedge):
    pass

  #{{{ as_kinemage
  #-----------------------------------------------------------------------------
  def as_kinemage(self):
    #output cablam validation as standalone kinemage markup for viewing in KiNG
    self.out.write('\n@subgroup {cablam disfavored} dominant\n')
    self.out.write('@vectorlist {cablam disfavored} color= purple width= 4 master={cablam disfavored} off') #default off
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.cablam_disfavored:
        result.as_kinemage(mode="cablam", out=self.out)
    self.out.write('\n@subgroup {cablam outlier} dominant\n')
    self.out.write('@vectorlist {cablam outlier} color= magenta width= 4 master={cablam outlier}') #default on
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.cablam_outlier:
        result.as_kinemage(mode="cablam", out=self.out)
    self.out.write('\n@subgroup {ca geom outlier} dominant\n')
    self.out.write('@vectorlist {ca geom outlier} color= red width= 4 master={ca geom outlier}') #default on
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.c_alpha_geom_outlier:
        result.as_kinemage(mode="ca_geom", out=self.out)
    #----------------------------
    #"wheels" show favorable and unfavorabe regions for each peptide plane involved in a cablam outlier
    #Some additional calculations are required to generate these wheels
    cablam_contours = fetch_peptide_expectations()
    wheels_list = []
    for result in self.results:
      if not result.has_ca:
        continue
      if result.feedback.cablam_disfavored:
        wheels_list.extend(result.calculate_kinemage_wheels(cablam_contours=cablam_contours))
    #wheels are made of 10-degree wedges, each wedge drawn as a triangle
    self.out.write('\n@subgroup {cablam_wheels} dominant master={cablam wheels}\n')
    self.out.write('@trianglelist {cablam_wheels} alpha=0.75')
    for cablam_wheel in wheels_list:
      wheel = cablam_wheel[0]
      wheel_center = cablam_wheel[1]
      for wedge in wheel:
        if wedge is None:
          continue
        elif wedge.cablam_score < 0.01:
          color = 'magenta'
        else:
          color = 'purple'
        self.cablam_wheel_triangle(wheel_center, wedge, color)
        #self.out.write('\n{} P X %s %.3f %.3f %.3f' % (color, wheel_center[0], wheel_center[1], wheel_center[2]))
        #self.out.write('\n{} %s %.3f %.3f %.3f' % (color, wedge.start[0], wedge.start[1], wedge.start[2]))
        #self.out.write('\n{} %s %.3f %.3f %.3f' % (color, wedge.end[0], wedge.end[1], wedge.end[2]))
    #a thin black line outlining the wheel greatly aids visual interpretation
    self.out.write('\n@vectorlist {cablam_wheels_lines} color=deadblack width= 1 alpha=0.75')
    for cablam_wheel in wheels_list:
      wheel = cablam_wheel[0]
      wheel_center = cablam_wheel[1]
      prevwedge = wheel[-1]
      new_poly = ' P' #starts a new polyline in kinemage format, print this for each new wheel
      for wedge in wheel:
        if wedge and prevwedge:
          self.out.write('\n{}%s %.3f %.3f %.3f' % (new_poly, wedge.start[0], wedge.start[1], wedge.start[2]))
          self.out.write('\n{} %.3f %.3f %.3f' % (wedge.end[0], wedge.end[1], wedge.end[2]))
        elif wedge and not prevwedge:
          self.out.write('\n{}%s %.3f %.3f %.3f' % (new_poly, wheel_center[0], wheel_center[1], wheel_center[2]))
          self.out.write('\n{} %.3f %.3f %.3f' % (wedge.start[0], wedge.start[1], wedge.start[2]))
          self.out.write('\n{} %.3f %.3f %.3f' % (wedge.end[0], wedge.end[1], wedge.end[2]))
        elif prevwedge and not wedge:
          self.out.write('\n{}%s %.3f %.3f %.3f' % (new_poly, prevwedge.end[0], prevwedge.end[1], prevwedge.end[2]))
          self.out.write('\n{} %.3f %.3f %.3f' % (wheel_center[0], wheel_center[1], wheel_center[2]))
        else:
          prevwedge = wedge
          continue
        prevwedge = wedge
        new_poly=''
    #----------------------------
    self.out.write('\n')
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_full_kinemage
  #-----------------------------------------------------------------------------
  def as_full_kinemage(self,pdbid=''):
    #output cablam validation as kinemage markup on pdb model. Future version
    #  of this will also print ribbons, but standalone ribbon code must be
    #  developed first.
    #Pdb-to-kinemage printing has been hijacked from mmtbx.kinemage.validation
    #  That code was not meant to be used outside its original context, so this
    #  may be fragile.
    from mmtbx.kinemage.validation import get_kin_lots, build_name_hash
    from mmtbx import monomer_library
    from mmtbx.monomer_library import pdb_interpretation
    i_seq_name_hash = build_name_hash(pdb_hierarchy=self.pdb_hierarchy)
    sites_cart=self.pdb_hierarchy.atoms().extract_xyz()
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
    pdb_io=self.pdb_hierarchy.as_pdb_input()
    processed_pdb_file = pdb_interpretation.process(
        mon_lib_srv=mon_lib_srv,
        ener_lib=ener_lib,
        pdb_inp=pdb_io,
        #params=work_params.kinemage.pdb_interpretation,
        substitute_non_crystallographic_unit_cell_if_necessary=True)
    geometry = processed_pdb_file.geometry_restraints_manager()
    flags = geometry_restraints.flags.flags(default=True)
    #angle_proxies = geometry.angle_proxies
    pair_proxies = geometry.pair_proxies(flags=flags, sites_cart=sites_cart)
    bond_proxies = pair_proxies.bond_proxies
    quick_bond_hash = {}
    for bp in bond_proxies.simple:
      if (i_seq_name_hash[bp.i_seqs[0]][9:14] == i_seq_name_hash[bp.i_seqs[1]][9:14]):
        if quick_bond_hash.get(bp.i_seqs[0]) is None:
          quick_bond_hash[bp.i_seqs[0]] = []
        quick_bond_hash[bp.i_seqs[0]].append(bp.i_seqs[1])

    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        self.out.write(get_kin_lots(chain, bond_hash=quick_bond_hash, i_seq_name_hash=i_seq_name_hash, pdbID=pdbid, index=0, show_hydrogen=True))
    self.as_kinemage()
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_pointcloud_kinemage
  #-----------------------------------------------------------------------------
  def as_pointcloud_kinemage(self):#, out=self.out):
    print("@group {cablam-space points} dominant", file=self.out)
    print("@dotlist (cablam-space points)", file=self.out)
    for result in self.results:
      result.as_kinemage_point()
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_secondary_structure
  #-----------------------------------------------------------------------------
  def as_secondary_structure(self, conf_request=None):
    #returns CaBLAM secondary structure identification as phenix's preferred
    #  iotbx.secondary_structure objects
    #each individual beta strand is currently represented as a whole sheet
    #Proper sheet reporting will depend on developing or implementing a
    #  strand-to-sheet assembly that does not require H-bonds
    from iotbx.pdb import secondary_structure
    helix_i = 0
    sheet_i = 0
    helix_records = []
    strand_records = []

    chain_list = list(self.all_results.keys())
    chain_list.sort()
    for chain_id in chain_list:
      chain = self.all_results[chain_id]
      if conf_request in chain.conf_names:
        conf = conf_request
      else:
        conf = chain.conf_names[0]
      for record in chain.confs[conf].sec_struc_records:
        if record.segment_type == 'alpha' or record.segment_type == 'threeten':
          helix_i += 1
          if record.segment_type == 'alpha':
            helix_class = 1
          elif record.segment_type == 'threeten':
            helix_class = 5
          return_record = secondary_structure.pdb_helix(
            serial = helix_i,
            helix_id = helix_i,
            start_resname  = record.start.resname,
            start_chain_id = record.start.chain_id,
            #start_chain_id = " A",
            start_resseq   = record.start.resseq,
            start_icode    = record.start.icode,
            end_resname    = record.end.resname,
            end_chain_id   = record.end.chain_id,
            end_resseq     = record.end.resseq,
            end_icode      = record.end.icode,
            helix_class    = helix_class,
            comment = "",
            length = record.segment_length)
          helix_records.append(return_record)
        if record.segment_type == 'beta':
          sheet_i += 1
          strand_record = secondary_structure.pdb_strand(
            sheet_id       = sheet_i,
            strand_id      = 1,
            start_resname  = record.start.resname,
            start_chain_id = record.start.chain_id,
            start_resseq   = record.start.resseq,
            start_icode    = record.start.icode,
            end_resname    = record.end.resname,
            end_chain_id   = record.end.chain_id,
            end_resseq     = record.end.resseq,
            end_icode      = record.end.icode,
            sense          = 1
            )
          return_record=secondary_structure.pdb_sheet(
            sheet_id = sheet_i,
            n_strands = 1,
            strands = [strand_record],
            registrations = [None],
            hbond_list = []
            )
          strand_records.append(return_record)
    return secondary_structure.annotation(helices=helix_records,sheets=strand_records)
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_records
  #-----------------------------------------------------------------------------
  def as_records(self, conf_request=None):
    #outputs pdb-style HELIX and SHEET secondary structure records
    #uses the iotbx.secondary_structure object
    #By default, this returns the first conformation (alt) in each chain
    #Other conformations can be accessed with conf_request
    records = self.as_secondary_structure(conf_request=conf_request)
    self.out.write(records.as_pdb_str()+'\n')
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_records_and_pdb
  #-----------------------------------------------------------------------------
  def as_records_and_pdb(self, conf_request=None):
    #outputs pdb-style HELIX and SHEET secondary structure records, followed by
    #  the pdb file
    self.as_records(conf_request=conf_request)
    self.out.write(self.pdb_hierarchy.as_pdb_string())
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ as_coot_data
  #-----------------------------------------------------------------------------
  def as_coot_data(self):
    data = []
    for result in self.results:
      if result.is_outlier:
        data.append((result.chain_id, result.resid, result.resname, result.scores.cablam, result.xyz))
    return data
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ show_summary
  #-----------------------------------------------------------------------------
  def show_summary(self, out=sys.stdout, prefix="  "):
    #print whole-model cablam summary
    # double percent sign %% is how to get an escape-char'd % in string formatting
    if self.count_residues() == 0 and self.count_ca_residues == 0:
      out.write("SUMMARY: CaBLAM found no evaluable protein residues.\n")
    else:
      out.write(prefix+"SUMMARY: Note: Regardless of number of alternates, each residue is counted as having at most one outlier.\n")
      out.write(prefix+"SUMMARY: CaBLAM found %s full protein residues and %s CA-only residues\n" % (self.count_residues(),self.count_ca_residues()-self.count_residues()))
      out.write(prefix+"SUMMARY: %s residues (%.1f%%) have disfavored conformations. (<=5%% expected).\n" % (self.count_disfavored(),self.percent_disfavored()))
      out.write(prefix+"SUMMARY: %s residues (%.1f%%) have outlier conformations. (<=1%% expected)\n" % (self.count_outliers(),self.percent_outliers()))
      out.write(prefix+"SUMMARY: %s residues (%.2f%%) have severe CA geometry outliers. (<=0.5%% expected)\n" % (self.count_ca_outliers(),self.percent_ca_outliers()))
      out.write(prefix+"SUMMARY: %s residues (%.2f%%) are helix-like, %s residues (%.2f%%) are beta-like\n" % (self.count_helix(),self.percent_helix(),self.count_beta(),self.percent_beta()))
      out.write(prefix+"SUMMARY: %s residues (%.2f%%) are correctable to helix, %s residues (%.2f%%) are correctable to beta\n" % (self.count_correctable_helix(),self.percent_correctable_helix(),self.count_correctable_beta(),self.percent_correctable_beta()))

    #if self.count_residues() == 0:
    #  print >> out, "SUMMARY: CaBLAM found no fully evaluable protein residues."
    #  if self.count_ca_residues() > 0:
    #    print >> out,prefix+"SUMMARY: CaBLAM found %s CA geometry evaluable residues." % (self.count_ca_residues())
    #    print >> out,prefix+"SUMMARY: %.2f" % (self.percent_ca_outliers())+"% of these residues have severe CA geometry outliers. (<=0.5% expected)"
    #    print >> out,prefix+"SUMMARY: %.2f%% helix, %.2f%% beta" % (self.percent_helix(),self.percent_beta())
    #else:
    #  print >> out,prefix+"SUMMARY: Note: Regardless of number of alternates, each residue is counted as having at most one outlier."
    #  print >> out,prefix+"SUMMARY: CaBLAM found %s fully evaluable residues and %s CA geometry evaluable residues." % (self.count_residues(),self.count_ca_residues()-self.count_residues())
    #  print >> out,prefix+"SUMMARY: %.1f" % (self.percent_disfavored())+"% of these residues have disfavored conformations. (<=5% expected)"
    #  print >> out,prefix+"SUMMARY: %.1f" % (self.percent_outliers())+"% of these residues have outlier conformations. (<=1% expected)"
    #  print >> out,prefix+"SUMMARY: %.2f" % (self.percent_ca_outliers())+"% of these residues have severe CA geometry outliers. (<=0.5% expected)"
    #  print >> out,prefix+"SUMMARY: %.2fpct helix, %.2fpct beta" % (self.percent_helix(),self.percent_beta())
  #-----------------------------------------------------------------------------
  #}}}

  #{{{ gui_summary
  #-----------------------------------------------------------------------------
  def gui_summary(self):
    output = []
    return "GUI summary goes here!"
  #-----------------------------------------------------------------------------
  #}}}
#-------------------------------------------------------------------------------
#}}}
