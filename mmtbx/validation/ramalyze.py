
from __future__ import absolute_import, division, print_function
from mmtbx.validation import residue, validation, atom
from mmtbx.validation import graphics
from iotbx import data_plots
from libtbx import slots_getstate_setstate
from mmtbx.conformation_dependent_library import generate_protein_threes
from mmtbx.conformation_dependent_library.cdl_utils import get_c_ca_n
import sys
from scitbx.array_family import flex
from mmtbx.validation import utils
from mmtbx.validation.rotalyze import get_center
import mmtbx.rotamer
from mmtbx.rotamer import ramachandran_eval
from mmtbx.validation.fav_lists import fav_tables
from six.moves import range
import json

# XXX Use these constants internally, never strings!
RAMA_GENERAL = 0
RAMA_GLYCINE = 1
RAMA_CISPRO = 2
RAMA_TRANSPRO = 3
RAMA_PREPRO = 4
RAMA_ILE_VAL = 5

RAMALYZE_OUTLIER = 0
RAMALYZE_ALLOWED = 1
RAMALYZE_FAVORED = 2
RAMALYZE_ANY = 3
RAMALYZE_NOT_FAVORED = 4

# Values are rama_keys, from 0 to 5
res_types_dict = {"general": 0, "glycine":1, "cis-proline":2, "trans-proline":3,
             "pre-proline":4, "isoleucine or valine":5}

res_types = ["general", "glycine", "cis-proline", "trans-proline",
             "pre-proline", "isoleucine or valine"]
res_type_labels = ["General", "Gly", "cis-Pro", "trans-Pro", "pre-Pro",
                   "Ile/Val"]
res_type_plot_labels = ["all non-Pro/Gly residues", "Glycine", "cis-Proline",
  "trans-Proline", "pre-Proline residues", "Ile or Val"]
rama_types = ["OUTLIER", "Allowed", "Favored", "Any", "Allowed/Outlier"]
rama_type_labels = ["Outlier", "Allowed", "Favored", "Any", "Allowed/Outlier"]

class c_alpha(slots_getstate_setstate):
  """Container class used in the generation of kinemages."""
  __slots__ = ['id_str', 'xyz']
  def __init__(self, id_str, xyz):
    self.id_str = id_str
    self.xyz = xyz

class ramachandran(residue):
  """
  Result class for protein backbone Ramachandran analysis (phenix.ramalyze).
  """
  __rama_attr__ = [
    "res_type",
    "rama_type",
    "score",
    "phi",
    "psi",
    #"c_alphas",
    "markup",
    "model_id",
  ]
  __slots__ = residue.__slots__ + __rama_attr__

  @staticmethod
  def header():
    return "%-20s %-12s %10s %6s %-20s" % ("Residue", "Type", "Region", "Score",
      "Phi/Psi")

  def residue_type(self):
    return res_type_labels[self.res_type]

  def ramalyze_type(self):
    return rama_types[self.rama_type]

  def as_string(self):
    return "%-20s %-12s %10s %6.2f %10s" % (self.id_str(), self.residue_type(),
      self.ramalyze_type(), self.score,
      ",".join([ "%.1f" % x for x in [self.phi, self.psi] ]))

  def as_csv(self):
    rc = [self.chain_id.strip(), self.resseq.strip(),
      # self.icode, self.altloc,
      self.resname.strip()]
    rc += [self.residue_type(), self.ramalyze_type(), '%0.1f' % (self.score)]
    rc += [ "%.1f" % x for x in [self.phi, self.psi]]
    return ','.join(rc)

  # Backwards compatibility
  def id_str_old(self):
    return "%s%4s%1s %1s%s" % (self.chain_id, self.resseq, self.icode,
      self.altloc, self.resname)

  def format_old(self):
    return "%s:%.2f:%.2f:%.2f:%s:%s" % (self.id_str(), self.score,
      self.phi, self.psi, self.ramalyze_type(),
      res_types[self.res_type].capitalize())

  def as_kinemage(self):
    assert self.is_outlier()
#    ram_out = "{%s CA}P %s\n" % (self.c_alphas[0].id_str, "%.3f %.3f %.3f" %
#      self.c_alphas[0].xyz)
#    ram_out += "{%s CA} %s\n" % (self.c_alphas[1].id_str, "%.3f %.3f %.3f" %
#      self.c_alphas[1].xyz)
#    ram_out += "{%s CA} %s\n" % (self.c_alphas[2].id_str, "%.3f %.3f %.3f" %
#      self.c_alphas[2].xyz)
#    return ram_out
    return self.markup

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if s != 'markup' and hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    slots_as_dict["rama_type"] = rama_types[slots_as_dict["rama_type"]]
    res_type_index = slots_as_dict['res_type']
    slots_as_dict['res_type'] = res_types[res_type_index]
    slots_as_dict['res_type_label'] = res_type_labels[res_type_index]
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  # GUI output
  def as_table_row_phenix(self):
    return [ self.chain_id, "%1s%s %s" % (self.altloc,self.resname,self.resid),
             self.residue_type(), self.score, self.phi, self.psi ]

class ramachandran_ensemble(residue):
  """Container for results for an ensemble of residues"""
  __slots__ = ramachandran.__slots__
  def __init__(self, all_results):
    self._copy_constructor(all_results[0])
    self.res_type = all_results[0].res_type
    self.rama_type = [ r.rama_type for r in all_results ]
    self.phi = flex.double([ r.phi for r in all_results ])
    self.psi = flex.double([ r.psi for r in all_results ])
    self.score = flex.double([ r.score for r in all_results ])

  def phi_min_max_mean(self):
    return self.phi.min_max_mean()

  def psi_min_max_mean(self):
    return self.psi.min_max_mean()

  def score_statistics(self):
    return self.score.min_max_mean()

  def phi_range(self):
    pass

class ramalyze(validation):
  """
  Frontend for calculating Ramachandran statistics for a model.  Can directly
  generate the corresponding plots.
  """
  __slots__ = validation.__slots__ + ["out_percent", "fav_percent",
    "n_allowed", "n_favored", "n_allowed_by_model", "n_favored_by_model", "n_type", "_outlier_i_seqs" ]
  program_description = "Analyze protein backbone ramachandran"
  output_header = "residue:score%:phi:psi:evaluation:type"
  gui_list_headers = ["Chain","Residue","Residue type","Score","Phi","Psi"]
  gui_formats = ["%s", "%s", "%s", "%.2f", "%.1f", "%.1f"]
  wx_column_widths = [75, 125, 125, 100, 125, 125]

  def get_result_class(self) : return ramachandran

  def __init__(self,
      pdb_hierarchy,
      outliers_only=False,
      show_errors=False,
      out=sys.stdout,
      quiet=False):
    # Optimization hint: make it possible to pass
    # ramachandran_eval.RamachandranEval() from outside.
    # Better - convert this to using mmtbx.model.manager where
    # RamachandranEval is already available.
    validation.__init__(self)
    self.n_allowed = 0
    self.n_favored = 0
    self.n_allowed_by_model = {}
    self.n_favored_by_model = {}
    self.n_type = [ 0 ] * 6
    self._outlier_i_seqs = flex.size_t()
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if (all_i_seqs.all_eq(0)):
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)
    analysis = ""
    output_list = []
    count_keys = []
    uniqueness_keys = []
    r = ramachandran_eval.RamachandranEval()
    ##if use_segids:
    ##      chain_id = utils.get_segid_as_chainid(chain=chain)
    ##    else:
    ##      chain_id = chain.id
    for three in generate_protein_threes(hierarchy=pdb_hierarchy, geometry=None):
      main_residue = three[1]
      phi_psi_atoms = three.get_phi_psi_atoms()
      if phi_psi_atoms is None:
        continue
      phi_atoms, psi_atoms = phi_psi_atoms
      phi = get_dihedral(phi_atoms)
      psi = get_dihedral(psi_atoms)
      coords = get_center(main_residue) #should find the CA of the center residue

      if (phi is not None and psi is not None):
        model_id = main_residue.parent().parent().parent().id
        if model_id not in self.n_total_by_model:
          self.n_total_by_model[model_id] = 0
          self.n_outliers_by_model[model_id] = 0
          self.n_allowed_by_model[model_id] = 0
          self.n_favored_by_model[model_id] = 0
        res_type = RAMA_GENERAL
        #self.n_total += 1
        if (main_residue.resname[0:3] == "GLY"):
          res_type = RAMA_GLYCINE
        elif (main_residue.resname[0:3] == "PRO"):
          is_cis = is_cislike_peptide(three)
          if is_cis:
            res_type = RAMA_CISPRO
          else:
            res_type = RAMA_TRANSPRO
        elif (three[2].resname == "PRO"):
          res_type = RAMA_PREPRO
        elif (main_residue.resname[0:3] == "ILE" or \
              main_residue.resname[0:3] == "VAL"):
          res_type = RAMA_ILE_VAL
        #self.n_type[res_type] += 1
        value = r.evaluate(res_types[res_type], [phi, psi])
        ramaType = self.evaluateScore(res_type, value)
        is_outlier = ramaType == RAMALYZE_OUTLIER

        c_alphas = None
        # XXX only save kinemage data for outliers
        if is_outlier :
          c_alphas = get_cas_from_three(three)
          assert (len(c_alphas) == 3)
          markup = self.as_markup_for_kinemage(c_alphas)
        else:
          markup = None
        result = ramachandran(
          model_id=main_residue.parent().parent().parent().id,
          chain_id=main_residue.parent().parent().id,
          resseq=main_residue.resseq,
          icode=main_residue.icode,
          resname=main_residue.resname,
          #altloc=main_residue.parent().altloc,
          altloc=get_altloc_from_three(three),
          segid=None, # XXX ???
          phi=phi,
          psi=psi,
          rama_type=ramaType,
          res_type=res_type,
          score=value*100,
          outlier=is_outlier,
          xyz=coords,
          markup=markup)
        #if result.chain_id+result.resseq+result.icode not in count_keys:
        result_key = result.model_id+result.chain_id+result.resseq+result.icode
        if result.altloc in ['','A'] and result_key not in count_keys:
          self.n_total += 1
          self.n_type[res_type] += 1
          self.n_total_by_model[model_id] += 1
          self.add_to_validation_counts(ramaType, model_id)
          count_keys.append(result_key)
        if (not outliers_only or is_outlier):
          if (result.altloc != '' or
            result_key not in uniqueness_keys):
            #the threes/conformers method results in some redundant result
            #  calculations in structures with alternates. Using the
            #  uniqueness_keys list prevents redundant results being added to
            #  the final list
            self.results.append(result)
            uniqueness_keys.append(result_key)
        if is_outlier :
          i_seqs = main_residue.atoms().extract_i_seq()
          assert (not i_seqs.all_eq(0))
          self._outlier_i_seqs.extend(i_seqs)
    #self.results.sort(key=lambda r: (r.model_id,r. chain_id, r.resseq_as_int(), r.icode, r.altloc))
    self.results.sort(key=lambda r: (r.model_id,r. chain_id, r.resseq_as_int(), r.icode, r.altloc))
    out_count, out_percent = self.get_outliers_count_and_fraction()
    fav_count, fav_percent = self.get_favored_count_and_fraction()
    self.out_percent = out_percent * 100.0
    self.fav_percent = fav_percent * 100.0

  def __add__(self, other):
    self.results += other.results
    return self

  def get_plots(self, show_labels=True, point_style='bo',
      markersize=10,markeredgecolor="black", dpi=100,markerfacecolor="white",
      show_filling=True, show_contours=True):
    """
    Create a dictionary of six PNG images representing the plots for each residue type.
    :param out: log filehandle
    """
    result = {}
    for pos in range(6):
      stats = utils.get_rotarama_data(
        pos_type=res_types[pos],
        convert_to_numpy_array=True)
      points, coords = self.get_plot_data(position_type=pos)
      result[pos] = draw_ramachandran_plot(
        points=points,
        rotarama_data=stats,
        position_type=pos,
        title=format_ramachandran_plot_title(pos, '*'),
        show_labels=show_labels,
        markeredgecolor=markeredgecolor,
        markerfacecolor=markerfacecolor,
        show_filling=show_filling,
        show_contours=show_contours,
        point_style=point_style,
        markersize=markersize)
    return result

  def write_plots(self, plot_file_base, out, show_labels=True, point_style='bo',
    markersize=10,markeredgecolor="black", show_filling=True, show_contours=True,
    dpi=100,markerfacecolor="white"):
    """
    Write a set of six PNG images representing the plots for each residue type.

    :param plot_file_base: file name prefix
    :param out: log filehandle
    """
    print("", file=out)
    print("Creating images of plots...", file=out)
    plots = self.get_plots(
        show_labels=show_labels,
        point_style=point_style,
        markersize=markersize,
        markeredgecolor=markeredgecolor,
        show_filling=show_filling,
        show_contours=show_contours,
        dpi=dpi,
        markerfacecolor=markerfacecolor)
    for pos in range(6):
      file_label = res_type_labels[pos].replace("/", "_")
      plot_file_name = plot_file_base + "_rama_%s.png" % file_label
      plots[pos].save_image(plot_file_name, dpi=dpi)
      print("  wrote %s" % plot_file_name, file=out)

  def display_wx_plots(self, parent=None,
      title="MolProbity - Ramachandran plots"):
    import wxtbx.plots.molprobity     # causes GUI error when moved to top?
    frame = wxtbx.plots.molprobity.ramalyze_frame(
      parent=parent, title=title, validation=self)
    frame.Show()
    return frame

  def show_summary(self, out=sys.stdout, prefix=""):
    print(prefix + 'SUMMARY: %i Favored, %i Allowed, %i Outlier out of %i residues (altloc A where applicable)' % (self.n_favored, self.n_allowed, self.n_outliers, self.n_total), file=out)
    print(prefix + 'SUMMARY: %.2f%% outliers (Goal: %s)' % \
      (self.out_percent, self.get_outliers_goal()), file=out)
    print(prefix + 'SUMMARY: %.2f%% favored (Goal: %s)' % \
      (self.fav_percent, self.get_favored_goal()), file=out)

  def get_plot_data(self, position_type=RAMA_GENERAL, residue_name="*",
      point_type=RAMALYZE_ANY):
    assert isinstance(position_type, int) and (0 <= position_type <= 5), \
      position_type
    points, coords = [], []
    for i, residue in enumerate(self.results):
      if ((residue.res_type == position_type) and
          ((residue_name == '*') or (residue_name == residue.resname))):
        if ((point_type == RAMALYZE_ANY) or
            (point_type == residue.rama_type) or
            ((residue.rama_type in [RAMALYZE_ALLOWED,RAMALYZE_OUTLIER]) and
             (point_type == RAMALYZE_NOT_FAVORED))):
          points.append((residue.phi, residue.psi, residue.simple_id(),
            residue.is_outlier()))
          coords.append(residue.xyz)
    return (points, coords)

  @staticmethod
  def evalScore(resType, value):
    if (value >= 0.02):
      return RAMALYZE_FAVORED
    if (resType == RAMA_GENERAL):
      if (value >= 0.0005):
        return RAMALYZE_ALLOWED
      else:
        return RAMALYZE_OUTLIER
    elif (resType == RAMA_CISPRO):
      if (value >=0.0020):
        return RAMALYZE_ALLOWED
      else:
        return RAMALYZE_OUTLIER
    else:
      if (value >= 0.0010):
        return RAMALYZE_ALLOWED
      else:
        return RAMALYZE_OUTLIER

  def evaluateScore(self, resType, value):
    ev = ramalyze.evalScore(resType, value)
    assert ev in [RAMALYZE_FAVORED, RAMALYZE_ALLOWED, RAMALYZE_OUTLIER]
    #if ev == RAMALYZE_FAVORED:
    #  self.n_favored += 1
    #elif ev == RAMALYZE_ALLOWED:
    #  self.n_allowed += 1
    #elif ev == RAMALYZE_OUTLIER:
    #  self.n_outliers += 1
    return ev

  def add_to_validation_counts(self, ev, model_id=""):
    if ev == RAMALYZE_FAVORED:
      self.n_favored += 1
      if model_id in self.n_favored_by_model:
        self.n_favored_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in ramalyze favored count dictionary, make sure you are calling this function with the correct model ID")
    elif ev == RAMALYZE_ALLOWED:
      self.n_allowed += 1
      if model_id in self.n_allowed_by_model:
        self.n_allowed_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in ramalyze allowed count dictionary, make sure you are calling this function with the correct model ID")
    elif ev == RAMALYZE_OUTLIER:
      self.n_outliers += 1
      if model_id in self.n_outliers_by_model:
        self.n_outliers_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in ramalyze outliers count dictionary, make sure you are calling this function with the correct model ID")

  def get_outliers_goal(self):
    return "< 0.2%"

  def _get_count_and_fraction(self, res_type):
    if (self.n_total != 0):
      count = self.n_type[res_type]
      fraction = float(count) / self.n_total
      return count, fraction
    return 0, 0.

  @property
  def percent_favored(self):
    n_favored, frac_favored = self.get_favored_count_and_fraction()
    return frac_favored * 100.

  @property
  def percent_allowed(self):
    n_allowed, frac_allowed = self.get_allowed_count_and_fraction()
    return frac_allowed * 100.

  def get_allowed_count_and_fraction(self):
    if (self.n_total != 0):
      fraction = self.n_allowed / self.n_total
      return self.n_allowed, fraction
    return 0, 0.

  def get_allowed_goal(self):
    return "> 99.8%"

  def get_favored_count_and_fraction(self):
    if (self.n_total != 0):
      fraction = self.n_favored / self.n_total
      return self.n_favored, fraction
    return 0, 0.

  def get_favored_goal(self):
    return "> 98%"

  def get_general_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_GENERAL)

  def get_gly_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_GLYCINE)

  def get_cis_pro_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_CISPRO)

  def get_trans_pro_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_TRANSPRO)

  def get_prepro_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_PREPRO)

  def get_ileval_count_and_fraction(self):
    return self._get_count_and_fraction(RAMA_ILE_VAL)

  def get_phi_psi_residues_count(self):
    return self.n_total

  def as_markup_for_kinemage(self,c_alphas):
    #atom.id_str() returns 'pdb=" CA  LYS    16 "'
    #The [9:-1] slice gives ' LYS    16 '
    if None in c_alphas: return ''
    ram_out = "{%s CA}P %s\n" % (c_alphas[0].id_str()[9:-1], "%.3f %.3f %.3f" %
      c_alphas[0].xyz)
    ram_out += "{%s CA} %s\n" % (c_alphas[1].id_str()[9:-1], "%.3f %.3f %.3f" %
      c_alphas[1].xyz)
    ram_out += "{%s CA} %s\n" % (c_alphas[2].id_str()[9:-1], "%.3f %.3f %.3f" %
      c_alphas[2].xyz)
    return ram_out

  def as_kinemage(self):
    ram_out = "@subgroup {Rama outliers} master= {Rama outliers}\n"
    ram_out += "@vectorlist {bad Rama Ca} width= 4 color= green\n"
    for rama in self.results :
      if rama.is_outlier():
        ram_out += rama.as_kinemage()
    return ram_out

  def as_coot_data(self):
    data = []
    for result in self.results :
      if result.is_outlier():
        data.append((result.chain_id, result.resid, result.resname,
          result.score, result.xyz))
    return data

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "ramalyze"
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
    data['summary_results'] = summary_results
    for model_id in self.n_total_by_model.keys():
      summary_results[model_id] = {"num_favored" : self.n_favored_by_model[model_id],
        "num_allowed" : self.n_allowed_by_model[model_id],
        "num_outliers" : self.n_outliers_by_model[model_id],
        "num_residues" : self.n_total_by_model[model_id],
        "outlier_percentage" : self.n_outliers_by_model[model_id]/self.n_total_by_model[model_id]*100,
        "outlier_goal" : self.get_outliers_goal(),
        "favored_percentage" : self.n_favored_by_model[model_id]/self.n_total_by_model[model_id]*100,
        "favored_goal" : self.get_favored_goal()
      }
    data['summary_results'] = summary_results
    #{summary: {"summary_text": "74 Favored, 0 Allowed, 0 Outlier out of 74 residues (altloc A where applicable)"}
    #           "info": {"res_count": 74,
    #                    "favored": 74,}
                #}
    #filebase = os.path.basename(self.data_manager.get_model_names()[0])
    return json.dumps(data, indent=2)

def get_matching_atom_group(residue_group, altloc):
  match = None
  if (residue_group != None):
    for ag in residue_group.atom_groups():
      if (ag.altloc == "" and match == None): match = ag
      if (ag.altloc == altloc): match = ag
  return match

def get_dihedral(four_atom_list):
  from cctbx import geometry_restraints
  if None in four_atom_list:
    return None
  return geometry_restraints.dihedral(
    sites=[atom.xyz for atom in four_atom_list],
    angle_ideal=-40,
    weight=1).angle_model

def get_phi(prev_atoms, atoms):
  import mmtbx.rotamer
  prevC, resN, resCA, resC = None, None, None, None;
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if (atom.name == " C  "): prevC = atom
  if (atoms is not None):
    for atom in atoms:
      if (atom.name == " N  "): resN = atom
      if (atom.name == " CA "): resCA = atom
      if (atom.name == " C  "): resC = atom
  if (prevC is not None and resN is not None and resCA is not None and resC is not None):
    return mmtbx.rotamer.phi_from_atoms(prevC, resN, resCA, resC)

def get_psi(atoms, next_atoms):
  import mmtbx.rotamer
  resN, resCA, resC, nextN = None, None, None, None
  if (next_atoms is not None):
    for atom in next_atoms:
      if (atom.name == " N  "): nextN = atom
  if (atoms is not None):
    for atom in atoms:
      if (atom.name == " N  "): resN = atom
      if (atom.name == " CA "): resCA = atom
      if (atom.name == " C  "): resC = atom
  if (nextN is not None and resN is not None and resCA is not None and resC is not None):
    return mmtbx.rotamer.psi_from_atoms(resN, resCA, resC, nextN)

def get_omega_atoms(three):
  ccn1, outl1 = get_c_ca_n(three[0])
  ccn2, outl2 = get_c_ca_n(three[1])
  if ccn1: ca1, c = ccn1[1], ccn1[0]
  else: ca1, c = None, None
  if ccn2: n, ca2 = ccn2[2], ccn2[1]
  else: n, ca2 = None, None
  #ca1, c, n, ca2 = ccn1[1], ccn1[0], ccn2[2], ccn2[1]
  omega_atoms = [ca1, c, n, ca2]
  return omega_atoms

def is_cis_peptide(three): #Depricated in favor of is_cislike_peptide below
  omega_atoms = get_omega_atoms(three)
  omega = get_dihedral(omega_atoms)
  if omega is None:
    return False
  if(omega > -30 and omega < 30):
    return True
  else:
    return False

def is_cislike_peptide(three):
  omega_atoms = get_omega_atoms(three)
  omega = get_dihedral(omega_atoms)
  if omega is None:
    return False
  if(omega > -90 and omega < 90):
    return True
    #Splitting omega at 90 accounts for twisted peptides and placed them in the cis/trans
    #  category closest to their modeled angle
    #This split is currently only relevant for twisted Proline, since cis nonPro is extremely rare
    #  and we do not define a separate Rama category for it
    #Users are advised to inspect every twisted Proline manually, since twisted peptides are higly
    #  unusual and should not be assumed to be close to correct. But they have to go somewhere on
    #  the Ramachandran plot
  else:
    return False

def get_cas_from_three(three):
  cas = []
  for residue in three:
    for atom in residue.atoms():
      if atom.name == " CA ":
        cas.append(atom)
        break
    else:
      cas.append(None)
  return cas
  ##  c_ca_n = get_c_ca_n(residue)
  ##  if c_ca_n[0] is None:
  ##    cas.append(None)
  ##  else:
  ##    cas.append(c_ca_n[0][1])
  ##return cas

def get_altloc_from_three(three):
  #look at all atoms in the main residues, plus the atoms used in calculations
  #  from adjacent residues to find if any have altlocs
  ##mc_atoms = (" N  ", " CA ", " C  ", " O  ")
  for atom in three[1].atoms():
    altchar = atom.parent().altloc
    if altchar != '':
      return altchar
  for atom in three[0].atoms():
    if atom.name != ' C  ':
      continue
    altchar = atom.parent().altloc
    if altchar != '':
      return altchar
  for atom in three[2].atoms():
    if atom.name != ' N  ':
      continue
    altchar = atom.parent().altloc
    if altchar != '':
      return altchar
  return ''

def construct_complete_residues(res_group):
  if (res_group is not None):
    complete_dict = {}
    nit, ca, co, oxy = None, None, None, None
    atom_groups = res_group.atom_groups()
    reordered = []
    # XXX always process blank-altloc atom group first
    for ag in atom_groups :
      if (ag.altloc == ''):
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
        if (atom.name in [" N  ", " CA ", " C  ", " O  "]):
          changed = True
      if (not None in [nit, ca, co, oxy]) and (changed):
        # complete residue backbone found
        complete_dict[ag.altloc] = [nit, ca, co, oxy]
    if len(complete_dict) > 0:
      return complete_dict
  return None

def isPrePro(residues, i):
  if (i < 0 or i >= len(residues) - 1): return False
  else:
    next = residues[i+1]
    for ag in next.atom_groups():
      if (ag.resname[0:3] == "PRO"): return True
  return False

def get_favored_regions(rama_key):
  """
  Returns list of tuples (phi, psi) inside separate favorable regions on
  particula Ramachandran plot.
  It is not the best idea to use strings, but it is not clear how
  conviniently use constants defined in the beginning of the file.
  """
  assert rama_key in range(6)

  if rama_key == RAMA_GENERAL:
    return [(-99, 119), (-63, -43), (53, 43), (60,-120)]
  if rama_key == RAMA_GLYCINE:
    return [(63, 41), (-63, -41), (79, -173), (-79, 173)]
  if rama_key == RAMA_CISPRO:
    return [(-75, 155), (-89, 5)]
  if rama_key == RAMA_TRANSPRO:
    # return [(-56, -55), (-55, 135)]
    return [(-57, -37), (-59, 143), (-81, 65)]
  if rama_key == RAMA_PREPRO:
    return [(-57, -45), (-100, 120), (49, 57)]
  if rama_key == RAMA_ILE_VAL:
    return [(-63, -45), (-119, 127)]
  return None

def get_favored_peaks(rama_key):
  """
  returns exact favored peaks with their score value
  """
  assert rama_key in range(6)

  if rama_key == RAMA_GENERAL:
    return [((-115.0, 131.0), 0.57068),
            ((-63.0, -43.0), 1.0),
            ((53.0, 43.0), 0.323004),
            ((53.0, -127.0), 0.0246619)]
  if rama_key == RAMA_GLYCINE:
    return [((63.0, 41.0), 1.0),
            ((-63.0, -41.0), 1.0),
            ((79.0, -173.0), 0.553852),
            # ((-79.0, 173.0), 0.553852),
            ]
  if rama_key == RAMA_CISPRO:
    return [((-75.0, 155.0), 1.0),
            ((-89.0, 5.0), 0.701149)]
  if rama_key == RAMA_TRANSPRO:
    return [((-57.0, -37.0), 0.99566),
            ((-59.0, 143.0), 1.0),
            ((-81.0, 65.0), 0.0896269)]
  if rama_key == RAMA_PREPRO:
    return [((-57.0, -45.0), 1.0),
            ((-67.0, 147.0), 0.992025),
            ((49.0, 57.0), 0.185259)]
  if rama_key == RAMA_ILE_VAL:
    return [((-63.0, -45.0), 1.0),
            ((-121.0, 129.0), 0.76163)]
  return None

def find_region_max_value(rama_key, phi, psi, allow_outside=False):
  def normalize(angle):
    a = int(angle)
    while a >= 180:
      a -= 360
    while a <= -180:
      a += 360
    return a
  from mmtbx.rotamer import ramachandran_eval
  from collections import Counter
  r = ramachandran_eval.RamachandranEval()
  value = r.evaluate(rama_key, [phi, psi])
  ev = ramalyze.evalScore(rama_key, value)
  if ev != RAMALYZE_FAVORED and not allow_outside:
    return None
  ph = int(phi)
  ps = int(psi)
  peaks = get_favored_peaks(rama_key)
  v = fav_tables[rama_key][normalize(ph)+180][normalize(ps)+180]
  values = []
  if v == 0:
    # look around, rounding problems
    for i in [-1,0,1]:
      for j in [-1,0,1]:
        values.append(fav_tables[rama_key][normalize(ph+i)+180][normalize(ps+j)+180])
    for e in Counter(values).elements():
      if e != 0:
        return peaks[e-1]
  if allow_outside:
    # do more comprehensive search, basically looking for the nearest
    # favorite region
    c = 1
    flag = True
    while flag:
      for i in range(-c,c):
        for j in range(-c,c):
          reg_number = fav_tables[rama_key][normalize(ph+i)+180][normalize(ps+j)+180]
          if reg_number != 0:
            flag = False
            return peaks[reg_number-1]
      c += 2
    return peaks[reg_number-1]
  if v == 0:
    return None
  else:
    return peaks[v-1]

#-----------------------------------------------------------------------
# GRAPHICS OUTPUT
def format_ramachandran_plot_title(position_type, residue_type):
  if (residue_type == '*'):
    title = "Ramachandran plot for " + res_type_plot_labels[position_type]
  else :
    title = "Ramachandran plot for " + residue_type
  return title

class ramachandran_plot_mixin(graphics.rotarama_plot_mixin):
  extent = [-180,180,-180,180]

  def set_labels(self, y_marks=()):
    self.plot.set_xlabel("Phi")
    self.plot.set_xticks([-120,-60,0,60,120])
    # self.plot.set_xticks([-160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160])
    self.plot.set_ylabel("Psi")
    self.plot.set_yticks([-120,-60,0,60,120])
    # self.plot.set_yticks([-160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160])
    self.plot.set_ylim((-182,182))
    self.plot.set_xlim((-182,182))
    # self.plot.grid(which='both', color='lime', linestyle='-', linewidth=2)

class ramachandran_plot(data_plots.simple_matplotlib_plot,
                         ramachandran_plot_mixin):
  def __init__(self, *args, **kwds):
    data_plots.simple_matplotlib_plot.__init__(self, *args, **kwds)
    ramachandran_plot_mixin.__init__(self, *args, **kwds)

def get_contours(position_type):
  '''
  Function for determining the contours in a Ramachandran plot

  Parameters
  ----------
  position_type: int, defined in beginning of file (e.g. RAMA_GENERAL)

  Returns
  -------
  list containing contours (2 numbers)

  data for plotting is being "scaled" in
  mmtbx/validation/utils.py: export_ramachandran_distribution():
    return npz ** scale_factor, # scale_factor = 0.25
  Therefore to calculate contours we need to look at
  mmtbx/validation/ramalyze.py: evalScore() for the logic and
  put the cutoff numbers to the power of 0.25
  '''
  general_contours = [0.1495, 0.376] # [0.0005**0.25, 0.02**0.25]
  cispro_contours = [0.21147, 0.376] # [0.002**0.25, 0.02**0.25]
  default_contours = [0.1778, 0.376] # [0.001**0.25, 0.02**0.25]

  contours = default_contours
  if position_type == RAMA_GENERAL :
    contours = general_contours
  elif position_type == RAMA_CISPRO :
    contours = cispro_contours
  return contours

def draw_ramachandran_plot(points,
                            rotarama_data,
                            position_type,
                            title,
                            show_labels=True,
                            markerfacecolor="white",
                            markeredgecolor="black",
                            show_filling=True,
                            show_contours=True,
                            markersize=10,
                            point_style='bo'):
  p = ramachandran_plot()
  contours = None
  if show_contours:
    contours = get_contours(position_type)
  p.draw_plot(
    stats=rotarama_data,
    title=title,
    points=points,
    show_labels=show_labels,
    colormap="Blues",
    contours=contours,
    markerfacecolor=markerfacecolor,
    markeredgecolor=markeredgecolor,
    show_filling=show_filling,
    markersize=markersize,
    point_style=point_style)
  return p
