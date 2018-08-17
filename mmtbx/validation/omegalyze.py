from __future__ import division
from mmtbx.validation import residue, validation, atom
import iotbx.phil
import os.path
from libtbx import slots_getstate_setstate
import numpy
import os, sys

################################################################################
# omegalyze.py
# This is a class to assess the omega (peptide bond) dihedral in protein
# backbone.  It originated following concerns that cis-peptides, especially
# non-proline cis-peptides were not being flagged by MolProbity, and so
# structures with an improbable over-abundance of cis-peptides were passing
# validation.
#
# This code reuses existing ramalyze code, structure, and naming conventions
# where possible.  Interfacing with this code should be either the same as or
# parallel to interfacing with ramalyze code.
#
################################################################################

def get_master_phil():
  return iotbx.phil.parse(input_string="""
    model = None
      .type = path
      .help = "Model file (PDB or mmCIF)"
    nontrans_only = True
      .type = bool
      .help = "Controls whether trans peptides are stored and printed"
    text = True
      .type = bool
      .help = "Prints verbose, colon-delimited text output and summary"
    kinemage = False
      .type = bool
      .help = "Prints kinemage markup for cis-peptides"
    oneline = False
      .type = bool
      .help = "Prints oneline-style summary statistics"
    help = False
      .type = bool
      .help = "Prints this help message if true"
""", process_includes=True)
prog = os.getenv('LIBTBX_DISPATCHER_NAME')
usage_string = """
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB or mmCIF file
  nontrans_only=True    only print nontrans residues (does not affect kinemage)
  text=True          verbose colon-delimited text output (default)
  kinemage=False        Create kinemage markup (overrides text output)
  help = False          Prints this help message if true

  text output is colon-delimited and follows the format:
    residue:type:omega:conformation
      'residue' is a unique residue identifier
      'type' is either proline or general case
      'omega' is the calculated omega dihedral for the peptide between this
        residue and the preceeding residue
      'conformation' is: cis for omega within 30 degrees of planar cis
                         trans for omega within 30 degrees of planar trans
                         twisted for omega not within 30 degrees of planar

  SUMMARY statistics provide:
    counts of cis prolines and twisted prolines relative to total prolines with
      measurable omega dihedrals across all chains
    counts of non-proline cis and twisted peptides relative to total non-proline
      peptides with measurable omega dihedrals across all chains

  Cis Prolines occur in ~5%% of prolines (1 in 20) at high resolution
  Non-Proline Cis residues occur in ~0.05%% of residues (1 in 2000) and require
    clear support from experimental data or homology.
  Twisted peptides are even less frequent and are highly suspect without
    high-resolution data.

Example:

  %(prog)s model=1ubq.pdb kinemage=True
""" % locals()

#{{{ XXX Use these constants internally, not the strings
OMEGA_GENERAL = 0
OMEGA_PRO = 1

OMEGALYZE_TRANS   =0
OMEGALYZE_CIS     =1
OMEGALYZE_TWISTED =2

res_types = ["non-proline", "proline"] #used in GUI table
res_type_labels = ["non-Pro", "Pro    "] #used in text output for MolProbity
res_type_kin = ["nonPro", "Pro"]
omega_types = ["Trans", "Cis", "Twisted"]
#}}}

class kin_atom(slots_getstate_setstate):
  """Container class used in generation of kinemages."""
  __slots__ = ['id_str','xyz']
  def __init__(self, id_str, xyz):
    self.id_str = id_str
    self.xyz = xyz

def dist(p1,p2):
  return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)**0.5

def single_offset(p1,p2,offset):
  d = dist(p1,p2)
  return (p1[0]-(p1[0]-p2[0])/d*offset, p1[1]-(p1[1]-p2[1])/d*offset, p1[2]-(p1[2]-p2[2])/d*offset)

class omega_result(residue):
  """
  Result class for protein backbone omega angle analysis (molprobity.omegalyze).
  """
  __omega_attr__ = [
    "res_type",
    "omega_type",
    "omega",
    "is_nontrans",
    "markup_atoms",
    "highest_mc_b",
    "prev_resseq",
    "prev_icode",
    "prev_resname",
    "prev_altloc"
  ]
  __slots__ = residue.__slots__ + __omega_attr__

  @staticmethod
  def header():
    return "%-31s %-12s %6s %-13s %6s" % ("Residues", "Type", "Omega","Conformation","MC_high_b")

  def residue_type(self):
    return res_type_labels[self.res_type]

  def omegalyze_type(self):
    return omega_types[self.omega_type]

  def residue_type_kin(self):
    return res_type_kin[self.res_type]

  def prev_id_str(self):
    return "%2s%4s%1s%1s%3s" % (
      self.chain_id, self.prev_resseq, self.prev_icode, self.prev_altloc,
      self.prev_resname)

  def as_string(self):
    return "%-12s to %-15s %-12s %6.2f %-13s %6.2f" % (
      self.prev_id_str(),self.id_str(),self.residue_type(),self.omega, self.omegalyze_type(),self.highest_mc_b)

  #For backwards compatibility
  def id_str_old(self):
    return "%s%4s%1s %1s%s" % (self.chain_id, self.resseq, self.icode,
      self.altloc, self.resname)

  def format_old(self):
    return "%s to %s: %s :%s:%s:%s" % (self.prev_id_str(),self.id_str(), self.residue_type(),
      ('%.2f'%self.omega).rjust(7), self.omegalyze_type().ljust(8),self.highest_mc_b)

  def as_kinemage(self, triangles=False, vectors=False):
    ca1,c,n,ca2  = self.markup_atoms[0].xyz, self.markup_atoms[1].xyz, self.markup_atoms[2].xyz, self.markup_atoms[3].xyz
    o = 0.1

    d1 = dist(ca1,ca2)
    d2 = dist(n,c)
    d3 = dist(ca1,c)
    d4 = dist(n,ca2)
    d_diag = dist(ca1,n)
    c_n_vec = (n[0]-c[0], n[1]-c[1], n[2]-c[2])
    c_ca1_vec = (ca1[0]-c[0], ca1[1]-c[1], ca1[2]-c[2])
    ca2_ca1_vec =  (ca1[0]-ca2[0], ca1[1]-ca2[1], ca1[2]-ca2[2])
    ca2_n_vec = (n[0]-ca2[0], n[1]-ca2[1], n[2]-ca2[2])
    diag_vec = (n[0]-ca1[0], n[1]-ca1[1], n[2]-ca1[2])

    theta = numpy.arccos(numpy.dot(c_ca1_vec,diag_vec)/d_diag/d3)
    ca1_offset_len = o/numpy.sin(theta)
    ca1_offset = single_offset(ca1,n,ca1_offset_len)

    theta = numpy.arccos(numpy.dot(ca2_n_vec,diag_vec)/d_diag/d4)
    n_offset_len = o/numpy.sin(theta)
    n_offset = single_offset(n,ca1,n_offset_len)

    theta_ca2 = numpy.arccos(numpy.dot(ca2_n_vec,ca2_ca1_vec)/d4/d1)
    ca2_to_ca1_offset_len = o/numpy.sin(theta_ca2)
    ca2_to_ca1_offset = single_offset(ca2,ca1,ca2_to_ca1_offset_len)
    vec_from_offset = (ca2_to_ca1_offset[0]+n[0]-ca2[0], ca2_to_ca1_offset[1]+n[1]-ca2[1], ca2_to_ca1_offset[2]+n[2]-ca2[2])
    theta = -numpy.arccos(numpy.dot(diag_vec,ca2_ca1_vec)/d_diag/d1)
    v1 = -dist(ca1,ca1_offset)*numpy.sin(theta)
    x = v1/numpy.sin(theta_ca2)
    ca2_offset = single_offset(ca2_to_ca1_offset,vec_from_offset, x)

    theta_n = numpy.arccos(numpy.dot(diag_vec,c_n_vec)/d_diag/d2)
    n_to_c_offset_len = o/numpy.sin(theta_n)
    n_to_c_offset = single_offset(c,n,n_to_c_offset_len)

    v2 = dist(n,n_offset)*numpy.sin(theta_n)
    theta_c = numpy.arccos(numpy.dot(c_ca1_vec,c_n_vec)/d2/d3)
    c_to_n_offset_len = o/numpy.sin(theta_c)
    c_to_n_offset = single_offset(c,n,c_to_n_offset_len)
    x2 = v2/numpy.sin(theta_c)
    vec_from_offset = (c_to_n_offset[0]+ca1[0]-c[0], c_to_n_offset[1]+ca1[1]-c[1], c_to_n_offset[2]+ca1[2]-c[2])

    c_offset = single_offset(c_to_n_offset,vec_from_offset,x2)

    #This commented block was the first pass at offseting the triangles
    #  it was abandoned due to moving the drawn ca1-n diagonal off of the actual diagonal
    #o = 0.1 #this is the offset from the backbone in angstrom
    #d1 = ((ca1[0]-ca2[0])**2 + (ca1[1]-ca2[1])**2 + (ca1[2]-ca2[2])**2)**0.5
    #d2 = ((n[0]-c[0])**2 + (n[1]-c[1])**2 + (n[2]-c[2])**2)**0.5
    #d3 = ((ca1[0]-c[0])**2 + (ca1[1]-c[1])**2 + (ca1[2]-c[2])**2)**0.5
    #d4 = ((n[0]-ca2[0])**2 + (n[1]-ca2[1])**2 + (n[2]-ca2[2])**2)**0.5
    #ca1_offset = (ca1[0]+(ca2[0]-ca1[0])/d1*o - (ca1[0]-c[0])/d3*o, ca1[1]+(ca2[1]-ca1[1])/d1*o - (ca1[1]-c[1])/d3*o, ca1[2]+(ca2[2]-ca1[2])/d1*o - (ca1[2]-c[2])/d3*o)
    #ca2_offset = (ca2[0]-(ca2[0]-ca1[0])/d1*o - (ca2[0]-n[0])/d4*o, ca2[1]-(ca2[1]-ca1[1])/d1*o - (ca2[1]-n[1])/d4*o, ca2[2]-(ca2[2]-ca1[2])/d1*o - (ca2[2]-n[2])/d4*o)
    #c_offset = (c[0]+(n[0]-c[0])/d2*o+(ca1[0]-c[0])/d3*o, c[1]+(n[1]-c[1])/d2*o+(ca1[1]-c[1])/d3*o, c[2]+(n[2]-c[2])/d2*o+(ca1[2]-c[2])/d3*o)
    #n_offset = (n[0]-(n[0]-c[0])/d2*o+(ca2[0]-n[0])/d4*o, n[1]-(n[1]-c[1])/d2*o+(ca2[1]-n[1])/d4*o, n[2]-(n[2]-c[2])/d2*o+(ca2[2]-n[2])/d4*o)
    if triangles:
      triangle1_line1 = "{%s CA  (%s %s, omega= %.2f)} P X %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % ca1_offset)
      triangle1_line2 =      "{%s C  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[1].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % c_offset)
      triangle1_line3 =      "{%s N  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[2].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % n_offset)
      triangle2_line1 = "{%s CA  (%s %s, omega= %.2f)} P X %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % ca1_offset)
      triangle2_line2 =      "{%s N  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[2].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % n_offset)
      triangle2_line3 =     "{%s CA  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[3].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % ca2_offset)
      out_this = triangle1_line1 + triangle1_line2 + triangle1_line3 + triangle2_line1 + triangle2_line2 + triangle2_line3
    elif vectors:
      vector_line1 = "{%s CA  (%s %s, omega= %.2f)} P %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % ca1_offset)
      if self.omega_type == OMEGALYZE_CIS:
        vector_line2 = "{%s CA  (%s %s, omega= %.2f)} %s\n" % (
          self.markup_atoms[3].id_str, self.omegalyze_type(),
          self.residue_type_kin(), self.omega,
          "%.3f %.3f %.3f" % ca2_offset)
      elif self.omega_type == OMEGALYZE_TWISTED:
        vector_line2 = "{%s N  (%s %s, omega= %.2f)} %s\n" % (
          self.markup_atoms[2].id_str, self.omegalyze_type(),
          self.residue_type_kin(), self.omega,
          "%.3f %.3f %.3f" % n_offset)
      else:
        return ""
      out_this = vector_line1 + vector_line2
    return out_this

  #def as_table_row_phenix(self):
  #  return [ self.chain_id, "%s %s" % (self.resname, self.resid),
  #           res_types[self.res_type], self.omega, omega_types[self.omega_type] ]

  def as_table_row_phenix(self):
    #'%4s%1s' string formatting for previous residue matched string formatting within self.resid
    return [ self.chain_id, "%1s%s %4s%1s to %1s%s %s" % (self.prev_altloc, self.prev_resname, self.prev_resseq, self.prev_icode, self.altloc, self.resname, self.resid),
             res_types[self.res_type], self.omega, omega_types[self.omega_type] ]

#the ramachandran_ensemble class is only called in mmtbx/validation/ensembles
# and does not seem to provide functionality useful to omega analysis
#So it is omitted for the moment

class omegalyze(validation):
  """
  Frontend for calculating omega angle statistics for a model.
  """
  __slots__ = validation.__slots__ + [
    "residue_count",
    "omega_count",
    "_outlier_i_seqs"
    ]

  program_description = "Analyze protein backbone peptide dihedrals (omega)"
  output_header = "residues:type:omega:conformation:mc_bmax"

  gui_list_headers = ["Chain","Residues","Residue type","omega","conformation"]
  gui_formats = ["%s", "%s", "%s", "%.2f", "%s"]
  wx_column_widths = [75,200,125,125,125]

  def get_result_class(self): return omega_result

  def __init__(self,
      pdb_hierarchy,
      nontrans_only=False,
      out=sys.stdout,
      quiet=True):
    validation.__init__(self)
    self.residue_count = [0, 0]
    #[OMEGA_GENERAL, OMEGA_PRO]
    self.omega_count = [[0,0,0], [0,0,0]]
    #[OMEGA_GENERAL, OMEGA_PRO], then
    #[OMEGALYZE_TRANS, OMEGALYZE_CIS, OMEGALYZE_TWISTED]

    from mmtbx.validation import utils
    from scitbx.array_family import flex
    self._outlier_i_seqs = flex.size_t()
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if all_i_seqs.all_eq(0):
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)

    prev_rezes, next_rezes = None, None
    prev_resid = None
    cur_resseq = None
    next_resseq = None
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        prev_rezes, next_rezes = None, None
        prev_resid = None
        cur_resseq = None
        next_resseq = None
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain)
        else:
          chain_id = chain.id
        residues = list(chain.residue_groups())
        for i, residue_group in enumerate(residues):
          # The reason I pass lists of atom_groups to get_phi and get_psi is to
          # deal with the particular issue where some residues have an A alt
          # conf that needs some atoms from a "" alt conf to get calculated
          # correctly.  See 1jxt.pdb for examples.  This way I can search both
          # the alt conf atoms and the "" atoms if necessary.
          prev_atom_list, next_atom_list, atom_list = None, None, None
          if cur_resseq is not None:
            prev_rezes = rezes
            prev_resseq = cur_resseq
          rezes = construct_residues(residues[i])
          cur_resseq = residue_group.resseq_as_int()
          cur_icode = residue_group.icode.strip()
          if (i > 0):
            #check for insertion codes
            if (cur_resseq == residues[i-1].resseq_as_int()) :
              if (cur_icode == '') and (residues[i-1].icode.strip() == '') :
                continue
            elif (cur_resseq != (residues[i-1].resseq_as_int())+1):
              continue
          for atom_group in residue_group.atom_groups():
            alt_conf = atom_group.altloc
            if rezes is not None:
              atom_list = rezes.get(alt_conf)
            if prev_rezes is not None:
              prev_atom_list = prev_rezes.get(alt_conf)
              if (prev_atom_list is None):
                prev_keys = sorted(prev_rezes.keys())
                prev_atom_list = prev_rezes.get(prev_keys[0])
            omega=get_omega(prev_atom_list, atom_list)
            highest_mc_b = get_highest_mc_b(prev_atom_list, atom_list)
            if omega is not None:
              resname = atom_group.resname[0:3]
              coords = get_center(atom_group)
              if resname == "PRO":
                res_type = OMEGA_PRO
              else:
                res_type = OMEGA_GENERAL
              self.residue_count[res_type] += 1
              omega_type = find_omega_type(omega)
              is_nontrans = False
              if omega_type == OMEGALYZE_CIS or omega_type == OMEGALYZE_TWISTED:
                self.n_outliers += 1
                is_nontrans = True
              self.omega_count[res_type][omega_type] += 1
              markup_atoms = [None, None, None, None] #for kinemage markup
              if is_nontrans:
                for a in prev_atom_list:
                  if a is None: continue
                  a_ = atom(pdb_atom=a)
                  if a.name.strip() == "CA":
                    markup_atoms[0] = kin_atom(
                      id_str=a_.atom_group_id_str(),xyz=a_.xyz)
                  elif a.name.strip() == "C":
                    markup_atoms[1] = kin_atom(
                      id_str=a_.atom_group_id_str(),xyz=a_.xyz)
                for a in atom_list:
                  if a is None: continue
                  a_ = atom(pdb_atom=a)
                  if a.name.strip() == "N":
                    markup_atoms[2] = kin_atom(
                      id_str=a_.atom_group_id_str(),xyz=a_.xyz)
                  elif a.name.strip() == "CA":
                    markup_atoms[3] = kin_atom(
                      id_str=a_.atom_group_id_str(),xyz=a_.xyz)
                #------------
              #prevres=residues[i-1]
              #find prev res identities for printing
              prev_alts = []
              prev_resnames = {}
              for ag in residues[i-1].atom_groups():
                prev_alts.append(ag.altloc)
                prev_resnames[ag.altloc] = ag.resname
              if alt_conf in prev_alts:
                prev_altloc = alt_conf
              else:
                if len(prev_alts) > 1:
                  prev_altloc = prev_alts[1]
                else:
                  prev_altloc = prev_alts[0]
              prev_resname = prev_resnames[prev_altloc]
              #done finding prev res identities
              result = omega_result(
                chain_id=chain_id,
                resseq=residue_group.resseq,
                icode=residue_group.icode,
                resname=atom_group.resname,
                altloc=atom_group.altloc,
                prev_resseq=residues[i-1].resseq,
                prev_icode=residues[i-1].icode,
                prev_resname=prev_resname,
                prev_altloc=prev_altloc,
                segid=None,
                omega=omega,
                omega_type=omega_type,
                res_type=res_type,
                is_nontrans=is_nontrans,
                outlier=is_nontrans,
                highest_mc_b=highest_mc_b,
                xyz=coords,
                markup_atoms=markup_atoms)
              if is_nontrans or not nontrans_only: #(not nontrans_only or is_nontrans)
                self.results.append(result)
              if is_nontrans:
                i_seqs = atom_group.atoms().extract_i_seq()
                assert (not i_seqs.all_eq(0)) #This assert copied from ramalyze
                self._outlier_i_seqs.extend(i_seqs)

  def _get_count_and_fraction(self, res_type, omega_type):
    total = self.residue_count[res_type]
    if total == 0:
      return 0, 0.0
    else:
      count = self.omega_count[res_type][omega_type]
      fraction = float(count) / total
      return count, fraction

  def as_kinemage(self):
    outlist = []
    cisprolist = []
    cisnonprolist = []
    cisprovectorlist = []
    cisnonprovectorlist = []
    twistlist = []
    twistvectorlist = []
    cisprohead = ["@subgroup {Cis proline} dominant master= {Cis proline} off\n",
      "@trianglelist {cis pro omega triangles} color= sea\n"]
    cisnonprohead = [
      "@subgroup {Cis peptides} dominant master= {Cis non-proline}\n",
      "@trianglelist {cis nonpro omega triangles} color= lime\n"]
    twisthead = [
      "@subgroup {Twisted peptides} dominant master= {Twisted peptides}\n",
      "@trianglelist {twisted omega triangles} color= yellow\n"]
    cisprovectorhead = ["@vectorlist {cis pro omega vectors} color= sea width=3\n"]
    cisnonprovectorhead = ["@vectorlist {cis nonpro omega vectors} color= lime width=3\n"]
    twistvectorhead=[
    "@vectorlist {twisted omega vectors} color= yellow width=3\n"]
    for result in self.results:
      if result.omega_type == OMEGALYZE_CIS:
        if result.res_type == OMEGA_PRO:
          cisprolist.append(result.as_kinemage(triangles=True))
          cisprovectorlist.append(result.as_kinemage(vectors=True))
        else:
          cisnonprolist.append(result.as_kinemage(triangles=True))
          cisnonprovectorlist.append(result.as_kinemage(vectors=True))
      elif result.omega_type == OMEGALYZE_TWISTED:
        twistlist.append(result.as_kinemage(triangles=True))
        twistvectorlist.append(result.as_kinemage(vectors=True))
    if cisprolist:
      outlist = outlist + cisprohead + cisprolist + cisprovectorhead + cisprovectorlist
    if cisnonprolist:
      outlist = outlist + cisnonprohead + cisnonprolist + cisnonprovectorhead + cisnonprovectorlist
    if twistlist:
      outlist = outlist + twisthead + twistlist + twistvectorhead + twistvectorlist
    return "".join(outlist)
    #it's my understanding that .join(list) is more efficient than string concat

  def as_coot_data(self):
    data = []
    for result in self.results:
      if result.is_nontrans:
        data.append((result.chain_id, result.resid, result.resname, result.score, result.xyz))
    return data

  def show_summary(self, out=sys.stdout, prefix=""):
    print >> out, prefix + 'SUMMARY: %i cis prolines out of %i PRO' % (
      self.n_cis_proline(),
      self.n_proline())
    print >> out, prefix + 'SUMMARY: %i twisted prolines out of %i PRO' % (
      self.n_twisted_proline(),
      self.n_proline())
    print >> out, prefix + 'SUMMARY: %i other cis residues out of %i nonPRO' % (
      self.n_cis_general(),
      self.n_general())
    print >> out, prefix + 'SUMMARY: %i other twisted residues out of %i nonPRO' % (
      self.n_twisted_general(),
      self.n_general())

  def summary_only(self, out=sys.stdout, pdbid="pdbid"):
    out.write(os.path.basename(pdbid) + ":")
    if self.n_cis_proline() == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.n_cis_proline()/self.n_proline()*100)+":")
    if self.n_twisted_proline() == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.n_twisted_proline()/self.n_proline()*100)+":")
    out.write("%i" % self.n_proline() + ":")
    if self.n_cis_general() == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.n_cis_general()/self.n_general()*100)+":")
    if self.n_twisted_general() == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.n_twisted_general()/self.n_general()*100)+":")
    out.write("%i" % self.n_general() + "\n")

  def gui_summary(self):
    output = []
    if self.n_cis_proline() or self.n_proline():
      output.append('%i cis prolines out of %i PRO' % (self.n_cis_proline(),self.n_proline()))
    if self.n_twisted_proline():
      output.append('%i twisted prolines out of %i PRO' % (self.n_twisted_proline(),self.n_proline()))
    if self.n_cis_general():
      output.append('%i cis residues out of %i nonPRO' % (self.n_cis_general(),self.n_general()))
    if self.n_twisted_general():
      output.append('%i twisted residues out of %i nonPRO' % (self.n_twisted_general(),self.n_general()))
    return "\n".join(output)

  def n_proline(self):
    return self.residue_count[OMEGA_PRO]
  def n_trans_proline(self):
    return self.omega_count[OMEGA_PRO][OMEGALYZE_TRANS]
  def n_cis_proline(self):
    return self.omega_count[OMEGA_PRO][OMEGALYZE_CIS]
  def n_twisted_proline(self):
    return self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED]

  def n_general(self):
    return self.residue_count[OMEGA_GENERAL]
  def n_trans_general(self):
    return self.omega_count[OMEGA_GENERAL][OMEGALYZE_TRANS]
  def n_cis_general(self):
    return self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS]
  def n_twisted_general(self):
    return self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED]

def write_header(writeto=sys.stdout):
  writeto.write("residue:omega:evaluation\n")

def find_omega_type(omega):
  if (omega > -30) and (omega < 30): omega_type = OMEGALYZE_CIS
  elif (omega < -150) or (omega > 150): omega_type = OMEGALYZE_TRANS
  else: omega_type = OMEGALYZE_TWISTED
  return omega_type

#construct_residues was adpated from ramalyze.construct_complete_residues
#main change invloved accepting incomplete residues, since complete ones were
#  not necessary for calculation of omega
#This required checks against None-valued atoms elsewhere in the code, but
#  allows assessment of all models omega dihedrals
def construct_residues(res_group):
  if (res_group is not None):
    complete_dict = {}
    nit, ca, co, oxy = None, None, None, None
    atom_groups = res_group.atom_groups()
    reordered = []
    # XXX always process blank-altloc atom group first
    for ag in atom_groups :
      if (ag.altloc == '') :
        reordered.insert(0, ag)
      else :
        reordered.append(ag)
    for ag in reordered :
      changed = False
      for atom in ag.atoms():
        if (atom.name == " N  "): nit = atom
        if (atom.name == " CA "): ca = atom
        if (atom.name == " C  "): co = atom
        if (atom.name == " O  "): oxy = atom
        if (atom.name in [" N  ", " CA ", " C  ", " O  "]) :
          changed = True
      if (changed):
        complete_dict[ag.altloc] = [nit, ca, co, oxy]
    if len(complete_dict) > 0:
      return complete_dict
  return None

def get_omega(prev_atoms, atoms):
  import mmtbx.rotamer
  prevCA, prevC, thisN, thisCA = None, None, None, None
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if atom is not None:
        if (atom.name == " CA "): prevCA = atom
        if (atom.name == " C  "): prevC = atom
  if (atoms is not None):
    for atom in atoms:
      if atom is not None:
        if (atom.name == " N  "): thisN = atom
        if (atom.name == " CA "): thisCA = atom
  if (prevCA is not None and prevC is not None and thisN is not None and thisCA is not None):
    return mmtbx.rotamer.omega_from_atoms(prevCA, prevC, thisN, thisCA)
  else:
    return None

def get_highest_mc_b(prev_atoms, atoms):
  highest_mc_b = 0
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if atom is not None and atom.name in [" CA "," C  "," N  "," O  ","CB"]:
        if atom.b > highest_mc_b:
          highest_mc_b = atom.b
  if (atoms is not None):
    for atom in atoms:
      if atom is not None and atom.name in [" CA "," C  "," N  "," O  ","CB"]:
        if atom.b > highest_mc_b:
          highest_mc_b = atom.b
  return highest_mc_b

def get_center(ag):
  for atom in ag.atoms():
    if (atom.name == " N  "):
      return atom.xyz
  return ag.atoms().extract_xyz().mean()

def run (args, out=sys.stdout, quiet=False) :
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  #if (params.model is None or params.help) :
    #help printing is handled in iotbx.phil.process_command_line_with_files()
  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.hierarchy
  hierarchy.atoms().reset_i_seq()
  result = omegalyze(
    pdb_hierarchy=hierarchy,
    nontrans_only=params.nontrans_only,
    out=out,
    quiet=quiet)
  if params.kinemage:
    print >> out, result.as_kinemage()
  elif params.oneline:
    result.summary_only(out=out, pdbid=params.model)
  elif params.text:
    result.show_old_output(out=out, verbose=True)

if (__name__ == "__main__") :
  run(sys.argv[1:])
