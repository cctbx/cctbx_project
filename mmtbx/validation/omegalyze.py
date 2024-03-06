from __future__ import absolute_import, division, print_function
from mmtbx.validation import residue, validation, atom
import os.path
from libtbx import slots_getstate_setstate
import numpy
import os, sys
import collections
import json
from iotbx.pdb.hybrid_36 import hy36decode

from mmtbx.conformation_dependent_library import generate_protein_fragments

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
    "prev_altloc",
    "model_id"
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

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if s != 'markup_atoms' and hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    slots_as_dict["omega_type"] = omega_types[slots_as_dict["omega_type"]]
    return json.dumps(slots_as_dict, indent=2)

#{'': {'A': {'  41 ': {'': {'altloc': '',
#                           'atom_selection': None,
#                           'chain_id': 'A',
#                           'highest_mc_b': 28.76,
#                           'icode': ' ',
#                           'is_nontrans': True,
#                           'model_id': '',
#                           'occupancy': None,
#                           'omega': -14.27418253081719,
#                           'omega_type': 'Cis',
#                           'outlier': True,
#                           'prev_altloc': '',
#                           'prev_icode': ' ',
#                           'prev_resname': 'PHE',
#                           'prev_resseq': '  40',
#                           'res_type': 1,
#                           'resname': 'PRO',
#                           'resseq': '  41',
#                           'score': None,
#                           'segid': None,
#                           'xyz': [-9.124,
#                                   3.291,
#                                   40.391]}},

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

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
    "residue_count_by_model",
    "omega_count_by_model",
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
    self.residue_count_by_model = {}
    self.omega_count_by_model = {}

    from mmtbx.validation import utils
    from scitbx.array_family import flex
    self._outlier_i_seqs = flex.size_t()
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if all_i_seqs.all_eq(0):
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)

    first_conf_altloc = None
    prev_chain_id = None
    for twores in generate_protein_fragments(
        pdb_hierarchy,
        length=2,
        geometry=None,
        include_non_standard_peptides=True,
        include_d_amino_acids=True):
      main_residue = twores[1] #this is the relevant residue for id-ing cis-Pro
      conf_altloc = get_conformer_altloc(twores)
      omega_atoms = get_omega_atoms(twores)
      # omega_atoms is the list [CA1 C1 N2 CA2], with None for missing atoms
      if None in omega_atoms:
        continue
      twores_altloc = local_altloc_from_atoms(omega_atoms)

      model_id = twores[0].parent().parent().parent().id
      if model_id not in self.residue_count_by_model:
        self.residue_count_by_model[model_id] = [0,0]
        self.omega_count_by_model[model_id] = [[0,0,0],[0,0,0]]
      chain = main_residue.parent().parent()
      if use_segids:
        chain_id = utils.get_segid_as_chainid(chain=chain)
      else:
        chain_id = chain.id

      if chain_id != prev_chain_id: #if we've moved to a new chain...
        first_conf_altloc = conf_altloc #...reset reference altloc
        prev_chain_id = chain_id
      if (conf_altloc != first_conf_altloc) and twores_altloc == '':
        #skip non-alternate residues unless this is the first time thru a chain
        continue

      omega = get_omega(omega_atoms)
      if omega is None: continue
      omega_type = find_omega_type(omega)
      if omega_type == OMEGALYZE_TRANS:
        is_nontrans = False
      else:
        is_nontrans = True
        self.n_outliers += 1
      if main_residue.resname == "PRO": res_type = OMEGA_PRO
      else:                             res_type = OMEGA_GENERAL
      self.residue_count[res_type] += 1
      self.residue_count_by_model[model_id][res_type] += 1
      self.omega_count[res_type][omega_type] += 1
      self.omega_count_by_model[model_id][res_type][omega_type] += 1
      highest_mc_b = get_highest_mc_b(twores[0].atoms(),twores[1].atoms())
      coords = get_center(main_residue)
      markup_atoms = []
      for omega_atom in omega_atoms:
        markup_atoms.append(kin_atom(omega_atom.parent().id_str(), omega_atom.xyz))

      result = omega_result(
                model_id=model_id,
                chain_id=chain_id,
                resseq=main_residue.resseq,
                icode=main_residue.icode,
                resname=main_residue.resname,
                altloc=local_altloc_from_atoms(omega_atoms[2:]),
                prev_resseq=twores[0].resseq,
                prev_icode=twores[0].icode,
                prev_resname=twores[0].resname,
                prev_altloc=local_altloc_from_atoms(omega_atoms[0:2]),
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
        i_seqs = main_residue.atoms().extract_i_seq()
        assert (not i_seqs.all_eq(0)) #This assert copied from ramalyze
        self._outlier_i_seqs.extend(i_seqs)
    self.results.sort(key=lambda r: (r.model_id, r.chain_id, int(hy36decode(len(r.resseq), r.resseq)), r.icode, r.altloc))

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

  def as_JSON(self, addon_json={}):
    # self.chain_id, "%1s%s %4s%1s to %1s%s %s" % (self.prev_altloc, self.prev_resname, self.prev_resseq, self.prev_icode, self.altloc, self.resname, self.resid),
    #         res_types[self.res_type], self.omega, omega_types[self.omega_type] ]
    # keep names roughly the same
    #check name in program template
    # {model: {1: {chain: {A: {residue_group: {1A: results}}}}}}
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "omegalyze"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    #hierarchical_results_dict = collections.defaultdict(dict)
    #hierarchical_results_dict['model_id']['chain_id']['resseq+icode']['resname'] = 'result'
    #print(hierarchical_results_dict)
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    for model_id in self.residue_count_by_model.keys():
      summary_results[model_id] = {
        "num_cis_proline" : self.n_cis_proline_by_model(model_id),
        "num_twisted_proline" : self.n_twisted_proline_by_model(model_id),
        "num_proline" : self.n_proline_by_model(model_id),
        "num_cis_general" : self.n_cis_general_by_model(model_id),
        "num_twisted_general" : self.n_twisted_general_by_model(model_id),
        "num_general" : self.n_general_by_model(model_id),
      }
    data['summary_results'] = summary_results

    return json.dumps(data, indent=2)

  def show_summary(self, out=sys.stdout, prefix=""):
    print(prefix + 'SUMMARY: %i cis prolines out of %i PRO' % (
      self.n_cis_proline(),
      self.n_proline()), file=out)
    print(prefix + 'SUMMARY: %i twisted prolines out of %i PRO' % (
      self.n_twisted_proline(),
      self.n_proline()), file=out)
    print(prefix + 'SUMMARY: %i other cis residues out of %i nonPRO' % (
      self.n_cis_general(),
      self.n_general()), file=out)
    print(prefix + 'SUMMARY: %i other twisted residues out of %i nonPRO' % (
      self.n_twisted_general(),
      self.n_general()), file=out)

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

  def n_proline_by_model(self, model_id):
    return self.residue_count_by_model[model_id][OMEGA_PRO]
  def n_trans_proline_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_PRO][OMEGALYZE_TRANS]
  def n_cis_proline_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_PRO][OMEGALYZE_CIS]
  def n_twisted_proline_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_PRO][OMEGALYZE_TWISTED]

  def n_general_by_model(self, model_id):
    return self.residue_count_by_model[model_id][OMEGA_GENERAL]
  def n_trans_general_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_GENERAL][OMEGALYZE_TRANS]
  def n_cis_general_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_GENERAL][OMEGALYZE_CIS]
  def n_twisted_general_by_model(self, model_id):
    return self.omega_count_by_model[model_id][OMEGA_GENERAL][OMEGALYZE_TWISTED]

def write_header(writeto=sys.stdout):
  writeto.write("residue:omega:evaluation\n")

def find_omega_type(omega):
  if (omega > -30) and (omega < 30): omega_type = OMEGALYZE_CIS
  elif (omega < -150) or (omega > 150): omega_type = OMEGALYZE_TRANS
  else: omega_type = OMEGALYZE_TWISTED
  return omega_type

def get_omega(omega_atoms):
  #omega_atoms = [CA1 C1 N2 CA2]
  if None in omega_atoms:
    return None
  import mmtbx.rotamer
  return mmtbx.rotamer.omega_from_atoms(omega_atoms[0], omega_atoms[1], omega_atoms[2], omega_atoms[3])

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

def get_conformer_altloc(twores):
  return twores[0].parent().altloc #go to conformer level

def local_altloc_from_atoms(atom_list):
  for atom in atom_list:
    if atom is not None:
      altloc = atom.parent().altloc #go to atom_group level
      if altloc != '':
        return altloc
  return ''

def get_omega_atoms(twores):
  atomlist = [twores[0].find_atom_by(name=" CA "),
              twores[0].find_atom_by(name=" C  "),
              twores[1].find_atom_by(name=" N  "),
              twores[1].find_atom_by(name=" CA ")]
  return atomlist
