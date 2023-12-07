from __future__ import division
from __future__ import nested_scopes, generators, absolute_import
from __future__ import with_statement, print_function

from mmtbx.validation import residue
from mmtbx.suitename.suitename import compute
from mmtbx.suitename.suitenamedefs import findBase
from mmtbx.validation import rna_geometry
import json
import sys

#reasons is a simple list version of an Enum from suitenamedefs
reasons = ["","delta-1", "epsilon-1", "zeta-1", "alpha", "beta", "gamma", "delta"]

def usage(out=sys.stdout):
  print("""This is documentation for the suitealyze suitename interface
You should see this message if you ran the program with no input file
""",file=out)

class suite(residue):
  """
    RNA backbone "suite", analyzed using 'suitename' validation.
  """
  __slots__ = residue.__slots__ + [
    "suite_id",
    #"suite",
    "suiteness",
    #"triaged_angle",
    "model_id",
    "valid",
    "angle","angles",
    "deltaMinus","epsilon","zeta","alpha","beta","gamma","delta",
    "bin", "cluster", "distance", "suiteness", "situation", "issue", "comment", "pointMaster", "pointColor",
    "outlier",
    "atoms",
    "base"
  ]
  def __init__(self,resid,angles,atoms):
    residue.__init__(self)
    self.model_id = resid["model"]
    self.chain_id = resid["chain"]
    self.resseq = resid["resseq"]
    self.icode = resid["icode"]
    self.altloc = resid["alt"]
    self.resname = resid["resname"]
    self.suite_id = self.altloc + self.resname + self.chain_id + self.resseq + self.icode
    self.angles = angles
    assert len(self.angles) == 7, "wrong # of dihedrals passed to suite"
    self.valid = self.validate()
    self.angle = [0]+angles+[0] #chi dihedrals, needed for compatibility with existing suitename functions
    self.deltaMinus = self.angles[0]
    self.epsilon = self.angles[1]
    self.zeta = self.angles[2]
    self.alpha = self.angles[3]
    self.beta = self.angles[4]
    self.gamma = self.angles[5]
    self.delta = self.angles[6]
    #print(self.suiteclustername)
    self.outlier = False #False by default, set True in suitealyze after computations if needed
    self.atoms = atoms
    assert len(self.atoms) == 3, "wrong # of atom positions passed to suite"
    self.xyz = atoms[1]
    self.base = findBase("%3s" % self.resname)
    if self.base is None:
      # covers the case where there is a DNA residue in the chain
      self.base = "?"

  def validate(self):
    # make sure that angles deltaMinus through delta exist and are in range
    if None in self.angles:
      return False
    for angle in self.angles:
      if angle < 0 or angle > 360:
        return False
    return True

  @property
  def suite(self):
    return self.cluster.name

  @property
  def triaged_angle(self):
    if self.issue:
      return reasons[self.issue.value]
    else:
      return None

  @staticmethod
  def header():
    return "%-20s  %8s  %9s  %12s" % ("Suite ID", "suite", "suiteness",
                                      "triaged angle")

  def format_values(self):
    return "%-20s  %8s  %9s  %8s" % (self.suite_id, self.suite, self.suiteness,
                                     self.triaged_angle)

  def as_string(self, prefix=""):
    return prefix + self.format_values()

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if hasattr(self, s) and s != "bin" and s != "cluster" and s != "issue"]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    slots_as_dict["bin"] = self.bin.name
    slots_as_dict["cluster"] = self.cluster.name
    reason = ""
    note = ""
    if self.issue:
      reason = reasons[self.issue.value]
    if self.cluster.status == "wannabe":
      note = "wannabe"
    slots_as_dict["reason"] = reason
    slots_as_dict["note"] = note
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  def as_table_row_phenix(self):
    return [self.suite_id, self.suite, self.suiteness, self.triaged_angle]

  def as_text_output(self,filename=" ",show_causes=None):
    #"{} {} {} {:5.3f}{}{}\n".format(outIDs, suite.bin.name, suite.cluster.name, float(suite.suiteness), reason, note)
    #:1: C:   1: : :  A inc  __ 0.000
    #:1: C:   4: : :  C 33 p 1a 0.717
    #:1: C:   5: : :  U trig !! 0.000 delta
    resid = ":".join([filename, self.model_id.strip() or '1', "%2s" % self.chain_id, self.resseq, self.icode, "%1s" % self.altloc, self.resname])
    #self.model_id.strip() or '1' defaults to printing the model# as 1 if unnumbered
    reason = ""
    note = ""
    if self.issue:
      reason = " " + reasons[self.issue.value]
    elif self.comment:
      reason = " " + self.comment
    elif self.situation and show_causes:
      reason += " " + self.situation
    if self.cluster.status == "wannabe":
      note = " wannabe"
    outline = "{} {} {} {:5.3f}{}{}".format(resid, self.bin.name, self.cluster.name, float(self.suiteness), reason, note)
    return outline

  def as_kinemage_markup(self):
    #return a vectorlist bar, similar to the Ramachandran markup, for a suite outlier
    #centered at P atom, extending towards ends of suite but not quite touching
    import numpy
    if None in self.atoms:
      return ""
    offset = 0.95
    p = numpy.array(self.atoms[0])
    c = numpy.array(self.atoms[1])
    o = numpy.array(self.atoms[2])
    vec1 = p-c
    vec2 = p-o
    c_offset = p - vec1*offset
    o_offset = p - vec2*offset
    return """{%s}P %.3f %.3f %.3f
{%s} %.3f %.3f %.3f
{%s} %.3f %.3f %.3f""" % (self.suite_id, c_offset[0],c_offset[1],c_offset[2],
                          self.suite_id, p[0], p[1], p[2],
                          self.suite_id, o_offset[0],o_offset[1],o_offset[2])

  def as_kinemage_label(self):
    #return a kinemage label with the suite name offset from the P atom of that suite
    #self.atoms[0] is the P atom location
    if not self.atoms[0]:
      return ""
    return "{%4s} %.3f %.3f %.3f" % (self.suite, self.atoms[0][0], self.atoms[0][1], self.atoms[0][2])

def setOptions(optionsIn):
  """optionsIn may be the result of parseOptions above
     or the result of an argparse parse_args operation"""
  from mmtbx.suitename.suitename import loadOptions
  from suitenamedefs import globals
  from mmtbx.suitename import suites
  if not optionsIn:
    optionsIn = suites.parseOptions("report=true")
    globals.options = optionsIn
    loadOptions(optionsIn)
  else:
    globals.options = optionsIn.suitename
    loadOptions(optionsIn.suitename)

class suitealyze(rna_geometry):
  output_header = "#suiteID:suite:suiteness:triaged_angle"
  label = "Backbone torsion suites"
  gui_list_headers = ["Suite ID", "Suite", "Suiteness", "Triaged angles",]
  gui_formats = ["%s"] * 4
  wx_column_widths = [200] * 4
  def __init__(self,
               pdb_hierarchy = None,
               model = None,
               options = None,
               outliers_only = False):
    assert (pdb_hierarchy is not None) or (model is not None),"no model or hierarchy passed to suitealyze"
    rna_geometry.__init__(self)
    if pdb_hierarchy is None:
      pdb_hierarchy = model.get_hierarchy()
    setOptions(options)
    self.results = self.make_suite_dihedrals(pdb_hierarchy)
    self.results = compute(self.results)
    for result in self.results:
      if result.suite == "!!":
        result.outlier=True
    if outliers_only:
      outlier_results = []
      for outlier in self.select_results(include_suites=["!!"]):
        outlier_results.append(outlier)
      self.results = outlier_results
    self.model_list = list(set([result.model_id for result in self.results]))
    self.chain_list = list(set([result.chain_id for result in self.results]))
    self.n_outliers = self.count_outliers()
    self.n_total    = self.count_suites()

  def get_result_class(self) : return suite

  def select_results(self, model=None, chain=None, include_suites=[], exclude_suites=[]):
    for result in self.results:
      if model is not None and result.model_id != model:
        continue
      elif chain is not None and result.chain_id != chain:
        continue
      elif include_suites and result.suite not in include_suites:
        continue
      elif result.suite in exclude_suites:
        continue
      yield result

  def average_suiteness(self, model=None, chain=None, include_suites=[], exclude_suites=[]):
    #Return the average suiteness of the structure or a specific chain
    suitenesses = []
    for result in self.select_results(model=model, chain=chain,
                                      include_suites=include_suites,
                                      exclude_suites=exclude_suites):
      suitenesses.append(result.suiteness)
    if len(suitenesses) == 0:
      return 0
    else:
      return sum(suitenesses)/len(suitenesses)

  def average_suiteness_all(self):
    return self.average_suiteness(exclude_suites=["__"])
  def average_suiteness_1a(self):
    return self.average_suiteness(include_suites=['1a'])
  def average_suiteness_non1a(self):
    return self.average_suiteness(exclude_suites=['1a',"__"])

  def count_suites(self, model=None, chain=None, include_suites=[], exclude_suites=[]):
    count = 0
    for result in self.select_results(model=model, chain=chain,
                                      include_suites=include_suites,
                                      exclude_suites=exclude_suites):
      count += 1
    return count

  def count_outliers(self, model=None, chain=None):
    return self.count_suites(model=model, chain=chain, include_suites=["!!"])

  def count_suites_in_suiteness_range(self, s_min, s_max, model=None, chain=None, include_suites=[], exclude_suites=[]):
    count = 0
    for result in self.select_results(model=model, chain=chain,
                                      include_suites=include_suites,
                                      exclude_suites=exclude_suites):
      if result.is_outlier():
        continue
      elif s_min <= result.suiteness < s_max:
        count += 1
    return count

  def standard_suiteness_bin_counts(self, model=None, chain=None, include_suites=[], exclude_suites=[]):
    #text output uses the following 12 categories or bins for suiteness.
    #  Outliers, exactly 0, then advancing by 0.1 until the max of 1.0
    #  index 0 is outliers
    #  index 1 is suiteness 0
    #  other indices correspond to the range 0.index-2 to 0.index-1
    #See suiteness_summary_block function below for reference text output
    return [self.count_outliers(model, chain),
            self.count_suites_in_suiteness_range(0, 0.000001, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.000001, 0.1, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.1, 0.2, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.2, 0.3, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.3, 0.4, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.4, 0.5, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.5, 0.6, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.6, 0.7, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.7, 0.8, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.8, 0.9, model, chain, include_suites, exclude_suites),
            self.count_suites_in_suiteness_range(0.9, 1.1, model, chain, include_suites, exclude_suites)
            ]

  def suiteness_summary_block(self, model=None, chain=None, include_suites=[], exclude_suites=[]):
    # For all 35 suites: average suiteness== 0.550 (power==3.00) #this includes !!
    #      4 suites are  outliers
    #      0 suites have suiteness == 0
    #      3 suites have suiteness >  0 <.1
    #      1 suites have suiteness >=.1 <.2
    #      2 suites have suiteness >=.2 <.3
    #      2 suites have suiteness >=.3 <.4
    #      2 suites have suiteness >=.4 <.5
    #      2 suites have suiteness >=.5 <.6
    #      1 suites have suiteness >=.6 <.7
    #      5 suites have suiteness >=.7 <.8
    #      9 suites have suiteness >=.8 <.9
    #      4 suites have suiteness >=.9
    outlist = []
    bins = self.standard_suiteness_bin_counts(model=model, chain=chain, include_suites=include_suites,
                                              exclude_suites=exclude_suites)
    if include_suites == [] and exclude_suites == ["__"]:
      outlist.append("For all %i suites: average suiteness== %.3f (power==3.00)" % (self.count_suites(model=model, chain=chain, exclude_suites=exclude_suites), self.average_suiteness(model=model, chain=chain, exclude_suites=["__"])))
      outlist.append("%6i suites are  outliers" % bins[0])
    elif include_suites == ['1a']:
      outlist.append(" A form (1a) %i suites: average suiteness== %.3f (power==3.00)" % (self.count_suites(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites), self.average_suiteness(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites)))
    elif exclude_suites == ['1a',"__"]:
      outlist.append(" non-1a  has %i suites: average suiteness== %.3f (power==3.00)" % (self.count_suites(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites), self.average_suiteness(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites)))
    else:
      outlist.append(" selection has %i suites: average suiteness== %.3f (power==3.00)" % (self.count_suites(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites), self.average_suiteness(model=model, chain=chain, include_suites=include_suites, exclude_suites=exclude_suites)))
    outlist.append("%6i suites have suiteness == 0" % bins[1])
    outlist.append("%6i suites have suiteness >  0 <.1" % bins[2])
    for i in [3,4,5,6,7,8,9,10]:
      outlist.append("%6i suites have suiteness >=.%i <.%i" % (bins[i], i-2, i-1))
    outlist.append("%6i suites have suiteness >=.9" % bins[11])
    return "\n".join(outlist)

  def count_triage(self,model=None,chain=None):
    count = 0
    for result in self.select_results(model=model,chain=chain):
      if result.triaged_angle:
        count += 1
    return count

  def show_summary(self, out=sys.stdout, prefix=""):
    #Found 38 complete suites derived from 39 entries
    #3 suites were triaged, leaving 35 assigned to bins

    print("Found %i complete suites derived from %i entries" % (sum(1 for r in self.select_results(exclude_suites=["__"])), len(self.results)), file=out)
    print("%i suites were triaged, leaving %i assigned to bins" % (self.count_triage(), len(self.results)), file=out)
    print(self.suiteness_summary_block(exclude_suites=["__"]), file=out)
    print(self.suiteness_summary_block(include_suites=["1a"]), file=out)
    print(self.suiteness_summary_block(exclude_suites=["1a","__"]), file=out)

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "rna_suites"
    data = addon_json
    flat_results = []
    hierarchical_results = {}
    summary_results = {}
    for result in self.results:
      flat_results.append(json.loads(result.as_JSON()))
      hier_result = json.loads(result.as_hierarchical_JSON())
      hierarchical_results = self.merge_dict(hierarchical_results, hier_result)

    data['flat_results'] = flat_results
    data['hierarchical_results'] = hierarchical_results
    for mod_id in self.model_list:
      summary_results[mod_id] = {"num_outliers": self.count_outliers(mod_id),
                                 "num_suites": self.count_suites(mod_id)}
    data['summary_results'] = summary_results
    data['suitestrings'] = self.as_suitestrings()
    return json.dumps(data, indent=2)

  def show_old_output(self, out=sys.stdout, verbose=False):
    for result in self.results:
      print(result.as_text_output(), file=out)
    if verbose:
      self.show_summary(out=out)

  def as_kinemage_markup(self):
    outlist = []
    outlist.append("@subgroup {suites} nobutton")
    outlist.append("@labellist {suitelabels} color=gold off nobutton master={suite labels}")
    for result in self.results:
      outlist.append(result.as_kinemage_label())
    outlist.append("@vectorlist {suites} color=gold width=4 nobutton master={suite outliers}")
    for result in self.results:
      if result.is_outlier():
        outlist.append(result.as_kinemage_markup())
    return "\n".join(outlist)

  def suitestring_for_chain(self, model=None, chain=None, alt='A'):
    allowed_alts = ['',alt]
    outlist = []
    for result in self.select_results(model=model, chain=chain):
      if result.altloc not in allowed_alts:
        continue
      outlist.append(str(result.suite) + str(result.base))
    return ("".join(outlist))

  def suitestring_to_block(self, suitestring):
    from textwrap import wrap
    return("\n".join(wrap(suitestring,30)))

  def as_suitestrings(self,alt='A'):
    #blockform is a more human-readable format that adds a return character every 10 suites
    #the default is a single long string (for each chain)
    suitestrings = {}
    for model_id in self.model_list:
      suitestrings[model_id] ={}
      for chain_id in self.chain_list:
        suitestrings[model_id][chain_id] = self.suitestring_for_chain(model=model_id, chain=chain_id, alt=alt)
    return suitestrings

  def display_suitestrings(self, alt='A', blockform=True, out=sys.stdout):
    suitestrings = self.as_suitestrings(alt=alt)
    outlist = []
    for model_id in sorted(suitestrings):
      if len(self.model_list)>1:
        outlist.append("Model "+model_id)
      for chain_id in sorted(suitestrings[model_id]):
        outlist.append("Chain "+chain_id)
        suitestring = suitestrings[model_id][chain_id]
        if blockform:
          suitestring = self.suitestring_to_block(suitestring)
        outlist.append(suitestring)
        outlist.append("")
    print("\n".join(outlist), file=out)

  def make_suite_dihedrals(self, pdb_hierarchy):
    from mmtbx.conformation_dependent_library import generate_dna_rna_fragments
    first_altloc_by_chain = {}
    results = []
    prev_chain = None
    for residue_pair in generate_dna_rna_fragments(pdb_hierarchy,
                                           geometry=None,
                                           length=2,
                                           #include_non_linked=False,
                                           include_non_linked=True,
                                           backbone_only=False,
                                           include_non_standard_bases=True,
                                           retain_selection='all',
                                           verbose=False,
                                           ):
      chainid = residue_pair[1].parent().parent().id
      modelid = residue_pair[1].parent().parent().parent().id
      first_altloc = first_altloc_by_chain.get(chainid)
      if first_altloc is None:
        first_altloc = [conf.altloc for conf in residue_pair[1].parent().parent().conformers()][0]
        first_altloc_by_chain[chainid] = first_altloc
      conformer_altloc = residue_pair[1].parent().altloc
      if chainid != prev_chain or modelid != prev_model:
        #The first residue of a chain is handled separately
        #It will be incomplete, but it's important for counts and suitestring representation
        res0atoms = {" C5'": None, " C4'": None, " C3'": None, " O3'": None}
        if conformer_altloc == '':
          local_altloc = ''
        else:
          local_altloc = self.local_altloc_from_atoms(list(self.get_res1atoms(residue_pair[0]).values()))
        if not (local_altloc == '' and local_altloc != conformer_altloc and conformer_altloc != first_altloc):
          res1atoms = self.get_res1atoms(residue_pair[0])
          if res1atoms[" O2'"] is None:
            # this residue is not RNA and should be skipped
            prev_chain = chainid
            prev_model = modelid
            continue
          deltaMinus, epsilon, zeta, alpha, beta, gamma, delta = self.get_7_dihedrals(res0atoms, res1atoms)
          resid = {"model": modelid, "chain": chainid, "resseq": residue_pair[0].resseq, "icode": residue_pair[0].icode,
                   "alt": local_altloc, "resname": residue_pair[0].resname}
          atoms = [res1atoms[" P  "] and res1atoms[" P  "].xyz,  # if the atom is None, you get None
                   res0atoms[" C5'"] and res0atoms[" C5'"].xyz,  # if the atom exists, you get its xyz
                   res1atoms[" O3'"] and res1atoms[" O3'"].xyz]
          results.append(suite(resid, [deltaMinus, epsilon, zeta, alpha, beta, gamma, delta], atoms))
          prev_chain = chainid
          prev_model = modelid
      #first residue done, resume normal path
      if residue_pair.are_linked():
        res0atoms = self.get_res0atoms(residue_pair[0])
        if res0atoms[" O2'"] is None:
          # preceding residue is not RNA, pair should be treated as non-linked
          res0atoms = {" C5'": None, " C4'": None, " C3'": None, " O3'": None, " O2'": None}
      else:
        res0atoms = {" C5'": None, " C4'": None, " C3'": None, " O3'": None, " O2'":None}
      res1atoms = self.get_res1atoms(residue_pair[1])
      if res1atoms[" O2'"] is None:
        #this residue is not RNA and should be skipped
        continue

      if conformer_altloc == '':
        local_altloc = ''
      else:
        local_altloc = self.local_altloc_from_atoms(list(res0atoms.values()) + list(res1atoms.values()))

      if local_altloc == '' and local_altloc != conformer_altloc and conformer_altloc != first_altloc:
        #if all true, this residue pair is a duplicate of one already seen
        continue
      #resid_str = residue_pair[1].id_str()
      resseq = residue_pair[1].resseq
      icode = residue_pair[1].icode
      resname = residue_pair[1].resname
      #resid = ":".join([chainid, resname, local_altloc, resseq, icode])
      resid = {"model":modelid, "chain":chainid, "resseq":resseq, "icode":icode, "alt":local_altloc, "resname":resname}
      deltaMinus, epsilon, zeta, alpha, beta, gamma, delta = self.get_7_dihedrals(res0atoms, res1atoms)
      atoms = [res1atoms[" P  "] and res1atoms[" P  "].xyz, #if the atom is None, you get None
               res0atoms[" C5'"] and res0atoms[" C5'"].xyz, #if the atom exists, you get its xyz
               res1atoms[" O3'"] and res1atoms[" O3'"].xyz]
      results.append(suite(resid, [deltaMinus, epsilon, zeta, alpha, beta, gamma, delta], atoms))
    return results

  def local_altloc_from_atoms(self, atom_list):
    for atom in atom_list:
      if atom is not None:
        altloc = atom.parent().altloc
        if altloc != '':
          return altloc
    return ''

  def get_res0atoms(self,residue):
    return {" C5'": residue.find_atom_by(name=" C5'"),
            " C4'": residue.find_atom_by(name=" C4'"),
            " C3'": residue.find_atom_by(name=" C3'"),
            " O3'": residue.find_atom_by(name=" O3'"),
            " O2'": residue.find_atom_by(name=" O2'"), }

  def get_res1atoms(self,residue):
    return {" P  ": residue.find_atom_by(name=" P  "),
            " O5'": residue.find_atom_by(name=" O5'"),
            " C5'": residue.find_atom_by(name=" C5'"),
            " C4'": residue.find_atom_by(name=" C4'"),
            " C3'": residue.find_atom_by(name=" C3'"),
            " O3'": residue.find_atom_by(name=" O3'"),
            " O2'": residue.find_atom_by(name=" O2'"), }

  def get_7_dihedrals(self, res0atoms, res1atoms):
    deltaMinus = self.get_dihedral([res0atoms[" C5'"], res0atoms[" C4'"], res0atoms[" C3'"], res0atoms[" O3'"]])
    epsilon = self.get_dihedral([res0atoms[" C4'"], res0atoms[" C3'"], res0atoms[" O3'"], res1atoms[" P  "]])
    zeta = self.get_dihedral([res0atoms[" C3'"], res0atoms[" O3'"], res1atoms[" P  "], res1atoms[" O5'"]])
    alpha = self.get_dihedral([res0atoms[" O3'"], res1atoms[" P  "], res1atoms[" O5'"], res1atoms[" C5'"]])
    beta = self.get_dihedral([res1atoms[" P  "], res1atoms[" O5'"], res1atoms[" C5'"], res1atoms[" C4'"]])
    gamma = self.get_dihedral([res1atoms[" O5'"], res1atoms[" C5'"], res1atoms[" C4'"], res1atoms[" C3'"]])
    delta = self.get_dihedral([res1atoms[" C5'"], res1atoms[" C4'"], res1atoms[" C3'"], res1atoms[" O3'"]])
    return deltaMinus, epsilon, zeta, alpha, beta, gamma, delta

  def get_dihedral(self, four_atom_list):
    from cctbx import geometry_restraints
    if None in four_atom_list:
      return None
    dh = geometry_restraints.dihedral(
      sites=[atom.xyz for atom in four_atom_list],
      angle_ideal=-40,
      weight=1).angle_model
    if dh < 0:
      #suitename expects all dihedrals to be between 0 and 360
      dh += 360
    return dh
