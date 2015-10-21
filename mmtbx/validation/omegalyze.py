from __future__ import division
from mmtbx.validation import residue, validation, atom
import iotbx.phil
import os.path
from libtbx import slots_getstate_setstate
import sys

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

usage_string = """
phenix.omegalyze file.pdb [params.eff] [options ...]

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

  Cis Prolines occur in ~5% of prolines (1 in 20) at high resolution
  Non-Proline Cis residues occur in ~0.05% of residues (1 in 2000) and require
    clear support from experimental data or homology.
  Twisted peptides are even less frequent and are highly suspect without
    high-resolution data.

Example:

  phenix.omegalyze model=1ubq.pdb kinemage=True
"""

#{{{ XXX Use these constants internally, not the strings
OMEGA_GENERAL = 0
OMEGA_PRO = 1

OMEGALYZE_TRANS   =0
OMEGALYZE_CIS     =1
OMEGALYZE_TWISTED =2

res_types = ["general", "proline"]
res_type_labels = ["General", "Pro"]
res_type_kin = ["nonPro", "Pro"]
omega_types = ["Trans", "Cis", "Twisted"]
#}}}

class kin_atom(slots_getstate_setstate):
  """Container class used in generation of kinemages."""
  __slots__ = ['id_str','xyz']
  def __init__(self, id_str, xyz):
    self.id_str = id_str
    self.xyz = xyz

class omega_result(residue):
  """
  Result class for protein backbone omega angle analysis (phenix.omegalyze).
  """
  __omega_attr__ = [
    "res_type",
    "omega_type",
    "omega",
    "is_nontrans",
    "markup_atoms"
  ]
  __slots__ = residue.__slots__ + __omega_attr__

  @staticmethod
  def header():
    return "%-20s %-12s %6s %-20s" % ("Residue", "Type", "Omega","Conformation")

  def residue_type(self):
    return res_type_labels[self.res_type]

  def omegalyze_type(self):
    return omega_types[self.omega_type]

  def residue_type_kin(self):
    return res_type_kin[self.res_type]

  def as_string(self):
    return "%-20s %-12s %6.2f %-20s" % (self.id_str(), self.residue_type(),
      self.omega, self.omegalyze_type())

  #For backwards compatibility
  def id_str_old(self):
    return "%s%4s%1s %1s%s" % (self.chain_id, self.resseq, self.icode,
      self.altloc, self.resname)

  def format_old(self):
    return "%s:%s:%.2f:%s" % (self.id_str(), self.residue_type(),
      self.omega, self.omegalyze_type())

  def as_kinemage(self, triangles=False, vectors=False):
    if triangles:
      triangle1_line1 = "{%s CA  (%s %s, omega= %.2f)} P X %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[0].xyz)
      triangle1_line2 =      "{%s C  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[1].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[1].xyz)
      triangle1_line3 =      "{%s N  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[2].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[2].xyz)
      triangle2_line1 = "{%s CA  (%s %s, omega= %.2f)} P X %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[0].xyz)
      triangle2_line2 =      "{%s N  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[2].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[2].xyz)
      triangle2_line3 =     "{%s CA  (%s %s, omega= %.2f)} %s\n" % (
        self.markup_atoms[3].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[3].xyz)
      out_this = triangle1_line1 + triangle1_line2 + triangle1_line3 + triangle2_line1 + triangle2_line2 + triangle2_line3
    elif vectors:
      vector_line1 = "{%s CA  (%s %s, omega= %.2f)} P %s\n" % (
        self.markup_atoms[0].id_str, self.omegalyze_type(),
        self.residue_type_kin(), self.omega,
        "%.3f %.3f %.3f" % self.markup_atoms[0].xyz)
      if self.omega_type == OMEGALYZE_CIS:
        vector_line2 = "{%s CA  (%s %s, omega= %.2f)} %s\n" % (
          self.markup_atoms[3].id_str, self.omegalyze_type(),
          self.residue_type_kin(), self.omega,
          "%.3f %.3f %.3f" % self.markup_atoms[3].xyz)
      elif self.omega_type == OMEGALYZE_TWISTED:
        vector_line2 = "{%s N  (%s %s, omega= %.2f)} %s\n" % (
          self.markup_atoms[2].id_str, self.omegalyze_type(),
          self.residue_type_kin(), self.omega,
          "%.3f %.3f %.3f" % self.markup_atoms[2].xyz)
      else:
        return ""
      out_this = vector_line1 + vector_line2
    return out_this

  def as_table_row_phenix(self):
    return [ self.chain_id, "%s %s" % (self.resname, self.resid),
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
  output_header = "residue:type:omega:conformation"
  gui_list_headers = ["Chain","Residue","Residue type","omega","conformation"]
  gui_formats = ["%s", "%s", "%s", "%.2f", "%s"]
  wx_column_widths = [125]*5

  def get_result_class(self): return omega_result

  def __init__(self,
      pdb_hierarchy,
      nontrans_only,
      out,
      quiet):
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
              result = omega_result(
                chain_id=chain_id,
                resseq=residue_group.resseq,
                icode=residue_group.icode,
                resname=atom_group.resname,
                altloc=atom_group.altloc,
                segid=None,
                omega=omega,
                omega_type=omega_type,
                res_type=res_type,
                is_nontrans=is_nontrans,
                outlier=is_nontrans,
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
    cislist = []
    cisvectorlist = []
    twistlist = []
    twistvectorlist = []
    cishead = ["@subgroup {Cis peptides} dominant master= {Cis peptides}\n",
      "@trianglelist {cis omega triangles} color= sea\n"]
    twisthead = [
      "@subgroup {Twisted peptides} dominant master= {Twisted peptides}\n",
      "@trianglelist {twisted omega triangles} color= lime\n"]
    cisvectorhead = ["@vectorlist {cis omega vectors} color= sea width=3\n"]
    twistvectorhead=[
    "@vectorlist {twisted omega vectors} color= lime width=3\n"]
    for result in self.results:
      if result.omega_type == OMEGALYZE_CIS:
        cislist.append(result.as_kinemage(triangles=True))
        cisvectorlist.append(result.as_kinemage(vectors=True))
      elif result.omega_type == OMEGALYZE_TWISTED:
        twistlist.append(result.as_kinemage(triangles=True))
        twistvectorlist.append(result.as_kinemage(vectors=True))
    if cislist:
      outlist = outlist + cishead + cislist + cisvectorhead + cisvectorlist
    if twistlist:
      outlist = outlist + twisthead + twistlist + twistvectorhead + twistvectorlist
    return "".join(outlist)
    #it's my understanding that .join(list) is more efficient than string concat

##################################################
#gui_list_headers = ["Chain","Residue","Residue type","omega","conformation"]
# GUI output
###  def as_table_row_phenix(self):
###    return [ self.chain_id, "%s %s" % (self.resname, self.resid),
###             res_types[self.res_type], self.omega, omega_types[self.omega_type] ]

#from validation
##  def as_gui_table_data (self, outliers_only=True, include_zoom=False) :
##    """
##    Format results for display in the Phenix GUI.
##    """
##    table = []
##    for result in self.iter_results(outliers_only) :
##      extra = []
##      if (include_zoom) :
##        extra = result.zoom_info()
##      row = result.as_table_row_phenix()
##      assert (len(row) == len(self.gui_list_headers) == len(self.gui_formats))
##      table.append(row + extra)
##    return table

#from entity
##def zoom_info (self) :
##    """
##    Returns data needed to zoom/recenter the graphics programs from the Phenix
##    GUI.
##    """
##    return [ self.as_selection_string(), self.xyz ]

  def as_coot_data(self):
    data = []
    for result in self.results:
      if result.is_nontrans:
        data.append((result.chain_id, result.resid, result.resname, result.score, result.xyz))
    return data

  def show_summary(self, out=sys.stdout, prefix=""):
    print >> out, prefix + 'SUMMARY: %i cis prolines out of %i PRO' % (
      self.omega_count[OMEGA_PRO][OMEGALYZE_CIS],
      self.residue_count[OMEGA_PRO])
    print >> out, prefix + 'SUMMARY: %i twisted prolines out of %i PRO' % (
      self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED],
      self.residue_count[OMEGA_PRO])
    print >> out, prefix + 'SUMMARY: %i other cis residues out of %i nonPRO' % (
      self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS],
      self.residue_count[OMEGA_GENERAL])
    print >> out, prefix + 'SUMMARY: %i other twisted residues out of %i nonPRO' % (
      self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED],
      self.residue_count[OMEGA_GENERAL])

  def summary_only(self, out=sys.stdout, pdbid="pdbid"):
    out.write(os.path.basename(pdbid) + ":")
    if self.omega_count[OMEGA_PRO][OMEGALYZE_CIS] == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.omega_count[OMEGA_PRO][OMEGALYZE_CIS]/self.residue_count[OMEGA_PRO]*100)+":")
    if self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED] == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED]/self.residue_count[OMEGA_GENERAL]*100)+":")
    out.write("%i" % self.residue_count[OMEGA_PRO] + ":")
    if self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS] == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS]/self.residue_count[OMEGA_GENERAL]*100)+":")
    if self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED] == 0:
      out.write("0:")
    else:
      out.write('%.3f' % (self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED]/self.residue_count[OMEGA_GENERAL]*100)+":")
    out.write("%i" % self.residue_count[OMEGA_GENERAL] + "\n")

  def gui_summary(self):
    output = []
    
    if self.omega_count[OMEGA_PRO][OMEGALYZE_CIS]:
      output.append('%i cis prolines out of %i PRO' % (self.omega_count[OMEGA_PRO][OMEGALYZE_CIS],self.residue_count[OMEGA_PRO]))
    if self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED]:
      output.append('%i twisted prolines out of %i PRO' % (self.omega_count[OMEGA_PRO][OMEGALYZE_TWISTED],self.residue_count[OMEGA_PRO]))
    if self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS]:
      output.append('%i cis residues out of %i nonPRO' % (self.omega_count[OMEGA_GENERAL][OMEGALYZE_CIS],self.residue_count[OMEGA_GENERAL]))
    if self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED]:
      output.append('%i twisted residues out of %i nonPRO' % (self.omega_count[OMEGA_GENERAL][OMEGALYZE_TWISTED],self.residue_count[OMEGA_GENERAL]))
    return "\n".join(output)

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

def get_center(ag):
  coords = None
  for atom in ag.atoms():
    if (atom.name == " CA "):
      coords = atom.xyz
  return coords

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
    result.summary_only(pdbid=params.model)
  elif params.text:
    result.show_old_output(out=out, verbose=True)


if (__name__ == "__main__") :
  run(sys.argv[1:])
