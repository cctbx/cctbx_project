from __future__ import absolute_import, division, print_function
from mmtbx.validation import residue, validation
from mmtbx.validation import graphics
from iotbx import data_plots
from iotbx.pdb import aa_utils
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
import operator
import json
import os, sys
from scitbx.array_family import flex

OUTLIER_THRESHOLD = 0.003
ALLOWED_THRESHOLD = 0.02

class rotamer(residue):
  """
  Result class for protein sidechain rotamer analysis (molprobity.rotalyze).
  """
  __rotamer_attr__ = [
    "score",
    "evaluation",
    "rotamer_name",
    "chi_angles",
    "incomplete",
    "model_id"
  ]
  __slots__ = residue.__slots__ + __rotamer_attr__

  def __hash__(self):
    return self.as_string().__hash__()

  @staticmethod
  def header():
    return "%-20s %8s %6s   %-20s" % ("Residue", "Rotamer", "Score",
      "Chi angles")

  def get_chi1_chi2(self):
    if (len(self.chi_angles) < 2):
      raise ValueError("Less than 2 chi angles for this residue (%s)" %
        self.id_str())
    return self.chi_angles[0], self.chi_angles[1]

  def as_string(self):
    return "%-20s %8s %6.2f   %-20s" % (self.id_str(), self.rotamer_name,
      self.score, self.format_chi_angles())

  # Backwards compatibility for scripts that expect old rotalyze output
  def id_str_old(self):
    return "%s%4s%1s %s" % (self.chain_id, self.resseq, self.icode,
      self.altloc + self.resname.strip())

  def format_chi_angles(self, pad=False, sep=","):
    formatted = []
    for chi in self.chi_angles :
      if pad or (chi is not None):
        formatted.append(format_value("%.1f", chi,
          replace_none_with="").strip())
    return sep.join(formatted)

  # Old output
  def format_old(self):
    s_occ = format_value("%.2f", self.occupancy)
    s_score = format_value("%.1f", self.score)
    chis = list(self.chi_angles)
    return "%s:%s:%s:%s:%s:%s" % (self.id_str(), s_occ, s_score,
      self.format_chi_angles(pad=True, sep=":"),
      self.evaluation,self.rotamer_name)

  def as_JSON(self):
    serializable_slots = [s for s in self.__slots__ if hasattr(self, s)]
    slots_as_dict = ({s: getattr(self, s) for s in serializable_slots})
    #slots_as_dict["rama_type"] = rama_types[slots_as_dict["rama_type"]]
    return json.dumps(slots_as_dict, indent=2)

  def as_hierarchical_JSON(self):
    hierarchical_dict = {}
    hierarchy_nest_list = ['model_id', 'chain_id', 'resid', 'altloc']
    return json.dumps(self.nest_dict(hierarchy_nest_list, hierarchical_dict), indent=2)

  # GUI output
  def as_table_row_phenix(self):
    return [ self.chain_id, "%1s%s %s" % (self.altloc,self.resname,self.resid),
             self.score ] + self.format_chi_angles(pad=True).split(",")#list(self.chi_angles)

class rotamer_ensemble(residue):
  """Container for validation results for an ensemble of residues."""
  __slots__ = rotamer.__slots__
  def __init__(self, all_results):
    self._copy_constructor(all_results[0])
    self.score = [ r.score for r in all_results ]
    self.rotamer_name = [ r.rotamer_name for r in all_results ]
    self.chi_angles = [ r.chi_angles for r in all_results ]

  def rotamer_frequencies(self):
    rotamers = []
    for rot_id in set(self.rotamer_name):
      n_rotamer = self.rotamer_name.count(rot_id)
      rotamers.append((rot_id, n_rotamer))
    return sorted(rotamers, key=operator.itemgetter(1), reverse=True)

  def as_string(self):
    rotamers = self.rotamer_frequencies()
    rot_out = [ "%s (%d)" % (rid, n_rot) for (rid, n_rot) in rotamers ]
    return "%-20s %s" % (self.id_str(), ", ".join(rot_out))

class rotalyze(validation):
  __slots__ = validation.__slots__ + ["n_allowed", "n_allowed_by_model",
    "n_favored", "n_favored_by_model", "out_percent",
    "outlier_threshold", "data_version", "percent_favored","percent_allowed",
    "_outlier_i_seqs"]
  program_description = "Analyze protein sidechain rotamers"
  output_header = "residue:occupancy:score%:chi1:chi2:chi3:chi4:"
  output_header+= "evaluation:rotamer"
  gui_list_headers = ["Chain","Residue","Score","Chi1","Chi2","Chi3","Chi4"]
  gui_formats = ["%s", "%s", "%.2f", "%.1f", "%.1f", "%.1f", "%.1f"]
  wx_column_widths = [75, 120, 100, 100, 100, 100, 100]

  def get_result_class(self) : return rotamer

  def __init__(self, pdb_hierarchy,
      data_version="8000",
      outliers_only=False,
      show_errors=False,
      use_parent=False,
      out=sys.stdout,
      quiet=False):
    validation.__init__(self)
    self.n_allowed = 0
    self.n_favored = 0
    self.n_allowed_by_model = {}
    self.n_favored_by_model = {}
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    from mmtbx.rotamer.rotamer_eval import RotamerID
    from mmtbx.validation import utils
    self.data_version = data_version
    self._outlier_i_seqs = flex.size_t()
#   if self.data_version == "500":    self.outlier_threshold = 0.01
    if self.data_version == "8000": self.outlier_threshold = OUTLIER_THRESHOLD
    else: raise ValueError(
      "data_version given to RotamerEval not recognized (%s)." % data_version)
    sidechain_angles = SidechainAngles(show_errors)
    rotamer_evaluator = rotamer_eval.RotamerEval(
                             data_version=data_version)
    rotamer_id = rotamer_eval.RotamerID() # loads in the rotamer names
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=pdb_hierarchy)
    current_rotamers = {}
    for model in pdb_hierarchy.models():
      self.n_allowed_by_model[model.id] = 0
      self.n_favored_by_model[model.id] = 0
      self.n_outliers_by_model[model.id] = 0
      self.n_total_by_model[model.id] = 0
      for chain in model.chains():
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain)
        else:
          chain_id = chain.id
        for rg in chain.residue_groups():
          all_dict = construct_complete_sidechain(rg)
          for atom_group in rg.atom_groups():
            coords = get_center(atom_group)
            resname = atom_group.resname
            occupancy = get_occupancy(atom_group)
            kwargs = {
              "model_id" : model.id,
              "chain_id" : chain_id,
              "resseq" : rg.resseq,
              "icode" : rg.icode,
              "altloc" : atom_group.altloc,
              "resname" : resname,
              "xyz" : coords,
              "occupancy" : occupancy,
            }
            atom_dict = all_dict.get(atom_group.altloc)
            if use_parent:
              parent_name = aa_utils.get_aa_parent(atom_group.resname)
              if parent_name!=atom_group.resname: atom_group.resname=parent_name
            res_key = get_residue_key(atom_group=atom_group)
            try:
              chis = sidechain_angles.measureChiAngles(
                       atom_group,
                       atom_dict)#.get(conformer.altloc))
            except AttributeError:
              if show_errors:
                kwargs['incomplete'] = True
                result = rotamer(**kwargs)
                print('%s is missing some sidechain atoms' % \
                  result.id_str(), file=out)
                self.results.append(result)
              continue
            if (chis is not None):
              if None in chis:
                continue
              cur_res = resname.lower().strip()
              if cur_res == 'mse':
                cur_res = 'met'
              elif use_parent: cur_res=parent_name
              value = rotamer_evaluator.evaluate(cur_res, chis)
              if value is not None:
                self.n_total += 1
                self.n_total_by_model[model.id] += 1
                kwargs['score'] = value * 100
                wrap_chis = rotamer_id.wrap_chis(resname.strip(), chis,
                  symmetry=False)
                sym_chis = wrap_chis[:]
                sym_chis = rotamer_id.wrap_sym(resname.strip(), sym_chis)
                evaluation = self.evaluateScore(value, model.id)
                kwargs['evaluation'] = evaluation
                if evaluation == "OUTLIER":
                  kwargs['outlier'] = True
                  kwargs['rotamer_name'] = evaluation
                  self._outlier_i_seqs.extend(atom_group.atoms().extract_i_seq())
                else:
                  kwargs['outlier'] = False
                  if use_parent: resname=parent_name
                  kwargs['rotamer_name'] = rotamer_id.identify(resname,
                    wrap_chis)
                  #deal with unclassified rotamers
                  if kwargs['rotamer_name'] == '':
                    kwargs['rotamer_name'] = "UNCLASSIFIED"
                while (len(wrap_chis) < 4):
                  wrap_chis.append(None)
                kwargs['chi_angles'] = wrap_chis
                result = rotamer(**kwargs)
                if (result.is_outlier()) or (not outliers_only):
                  self.results.append(result)
    out_count, out_percent = self.get_outliers_count_and_fraction()
    self.out_percent = out_percent * 100.0
    assert abs(self.percent_outliers - self.out_percent) < 1.e-6
    self.percent_favored, self.percent_allowed = 0, 0
    if(self.n_total>0):
      self.percent_favored  = self.n_favored * 100. / self.n_total
      self.percent_allowed  = self.n_allowed * 100. / self.n_total
      # Checksum assert
      assert abs(self.percent_favored+self.percent_allowed+
                 self.percent_outliers-100.) < 1.e-6

  def get_rotamer_count(self):
    rotamer_count = {'total': {}}
    for residue in self.results:
      if residue.is_outlier(): continue
      current = rotamer_count['total'].setdefault(residue.resname, {})
      current.setdefault(residue.rotamer_name, 0)
      current[residue.rotamer_name]+=1
    return rotamer_count

  def evaluateScore(self, value, model_id=""):
    if value >= ALLOWED_THRESHOLD :
      self.n_favored += 1
      if model_id in self.n_favored_by_model:
        self.n_favored_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in favored count dictionary, make sure you are calling this function with the correct model ID")
      return "Favored"
    elif value >= self.outlier_threshold:
      self.n_allowed += 1
      if model_id in self.n_allowed_by_model:
        self.n_allowed_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in allowed count dictionary, make sure you are calling this function with the correct model ID")
      return "Allowed"
    else:
      self.n_outliers += 1
      if model_id in self.n_outliers_by_model:
        self.n_outliers_by_model[model_id] += 1
      else:
        raise Sorry("Model ID not found in outliers count dictionary, make sure you are calling this function with the correct model ID")
      return "OUTLIER"

  def show_summary(self, out=sys.stdout, prefix=""):
    print(prefix + 'SUMMARY: %.2f%% outliers (Goal: %s)' % \
      (self.out_percent, self.get_outliers_goal()), file=out)

  def get_outliers_goal(self):
#   if self.data_version == '500' : return "< 1%"
    return "< 0.3%"

  def get_favored_fraction_for_model(self, model_id):
    if (self.n_total_by_model[model_id] != 0):
      fraction = float(self.n_favored_by_model[model_id]) / self.n_total_by_model[model_id]
      assert fraction <= 1.0
      return fraction
    return 0.

  def get_favored_goal(self):
    return "> 98%"

  def coot_todo(self):
    return ""

  def get_plot_data(self, residue_name, point_type):
    assert (point_type in ["All", "None", "Outlier"])
    points = []
    coords = []
    for i, residue in enumerate(self.results):
      if (residue.resname == residue_name):
        if ((point_type == "All") or (residue.is_outlier())):
          chi1, chi2 = residue.get_chi1_chi2()
          points.append((chi1, chi2, residue.simple_id(), residue.is_outlier()))
          coords.append(residue.xyz)
    return (points, coords)

  def display_wx_plots(self, parent=None,
      title="MolProbity - Sidechain Chi1/Chi2 plots"):
    import wxtbx.plots.molprobity
    frame = wxtbx.plots.molprobity.rotalyze_frame(
      parent=parent, title=title, validation=self)
    frame.Show()
    return frame

  def as_JSON(self, addon_json={}):
    if not addon_json:
      addon_json = {}
    addon_json["validation_type"] = "rotalyze"
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
    for model_id in self.n_total_by_model:
      summary_results[model_id] = {
        "num_favored" : self.n_favored_by_model[model_id],
        "num_allowed" : self.n_allowed_by_model[model_id],
        "num_outliers" : self.n_outliers_by_model[model_id],
        "num_residues" : self.n_total_by_model[model_id],
        "outlier_percentage" : self.get_outliers_fraction_for_model(model_id) * 100,
        "outlier_goal" : self.get_outliers_goal(),
        "favored_percentage" : self.get_favored_fraction_for_model(model_id) * 100,
        "favored_goal" : self.get_favored_goal()
      }
    data['summary_results'] = summary_results
    return json.dumps(data, indent=2)

  def as_coot_data(self):
    data = []
    for result in self.results :
      if result.is_outlier():
        data.append((result.chain_id, result.resid, result.resname,
          result.score, result.xyz))
    return data

def evaluate_rotamer(
    atom_group,
    sidechain_angles,
    rotamer_evaluator,
    rotamer_id,
    all_dict,
    outlier_threshold=0.003,
    sites_cart=None):
  atom_dict = all_dict.get(atom_group.altloc)
  resname = atom_group.resname
  try:
    chis = sidechain_angles.measureChiAngles(atom_group, atom_dict, sites_cart)
    value = rotamer_evaluator.evaluate(
              atom_group.resname.lower().strip(),
              chis)
  except Exception:
    return None, None, None
  wrap_chis = rotamer_id.wrap_chis(resname.strip(), chis, symmetry=False)
  rotamer_name = rotamer_id.identify(resname.strip(), wrap_chis)
  if (value is None):
    return None, None, None
  elif (value < outlier_threshold):
    return 'OUTLIER', chis, value
  else:
    return rotamer_name, chis, value

def split_rotamer_names(rotamer):
  split_rotamer = []
  multi = ""
  if rotamer in ['UNCLASSIFIED', 'OUTLIER', 'Cg_exo', 'Cg_endo']:
    split_rotamer.append(rotamer)
    return split_rotamer
  for i, c in enumerate(rotamer):
    if c in ['t', 'p', 'm']:
      split_rotamer.append(c)
    elif (c in ['-', '?']) or (c >= '0' and c<='9'):
      multi += c
  if len(multi) > 0:
    split_rotamer.append(multi)
  return split_rotamer

def get_residue_key(atom_group):
  altloc = atom_group.altloc
  if altloc == "":
    altloc = " "
  key = None
  for atom in atom_group.atoms():
    cur_label = atom.pdb_label_columns()+atom.segid
    cur_altloc = cur_label[4:5]
    if key is None:
      if altloc == cur_altloc:
        key = cur_label[4:]
    else:
      if altloc == cur_altloc:
        if (key != cur_label[4:]):
          if os.getenv('CCP4'):
            outl = 'a model editor'
          else:
            outl = 'phenix.pdbtools or phenix.pdb_editor'
          raise Sorry("""\
Incompatible identifiers for one or more atoms in a residue:
%s
This is usually caused by atoms with a different segid from the rest of the
residue.  You can use %s to reset the
segid.""" % (atom.format_atom_record(), outl))
        assert key == cur_label[4:]
  return key

def evaluate_residue(
      residue_group,
      sa,
      r,
      all_dict,
      sites_cart=None):
  is_outlier = False
  data_version = "8000"
  outlier_threshold = OUTLIER_THRESHOLD
  for ag in residue_group.atom_groups():
    atom_dict = all_dict.get(ag.altloc)
    try:
      chis = sa.measureChiAngles(ag, atom_dict)
      value = r.evaluate(ag.resname.lower().strip(), chis, sites_cart)
    except Exception:
      #print ag.resname.lower()+residue_group.resseq+" is missing some sidechain atoms"
      value = None;
      is_outlier = None;
      return is_outlier, value
    if (value is None):
      is_outlier = False
      return is_outlier, value
    elif (value < outlier_threshold):
      is_outlier = True
      return is_outlier, value
    else:
      return is_outlier, value

class residue_evaluator(object):
  def __init__(self):
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    self.sa = SidechainAngles(False)
    self.data_version = "8000"
    self.r = rotamer_eval.RotamerEval(data_version=self.data_version)

  def evaluate_residue(self, residue_group):
    all_dict = construct_complete_sidechain(residue_group)
    return evaluate_residue(
      residue_group=residue_group,
      sa=self.sa,
      r=self.r,
      all_dict=all_dict,
      data_version=self.data_version)

  def __call__(self, *args, **kwds):
    return self.evaluate_residue(*args, **kwds)

def get_center(residue):
  for atom in residue.atoms():
    if atom.name == " CA ":
      return atom.xyz
  return residue.atoms().extract_xyz().mean()

def has_heavy_atoms(atoms):
  for atom in atoms:
    if not atom.element_is_hydrogen():
      return True
  return False

def construct_complete_sidechain(residue_group):
  if (residue_group is not None):
    complete_dict = {}
    atom_dict = {}
    for ag in residue_group.atom_groups():
      if not has_heavy_atoms(ag.atoms()):
        continue
      for atom in ag.atoms():
        #handle hydrogen/deuterium swaps
        if atom_dict.get(atom.name) == None:
          if atom_dict.get(atom.name.replace("H","D",1)) != None:
            del(atom_dict[atom.name.replace("H","D",1)])
          elif atom_dict.get(atom.name.replace("D","H",1)) != None:
            del(atom_dict[atom.name.replace("D","H",1)])
        atom_dict[atom.name] = atom
      clone_dict = {}
      clone_dict.update(atom_dict)
      complete_dict[ag.altloc] = clone_dict
    if len(complete_dict) > 0:
      return complete_dict
  return {}

# XXX does this need to be smarter?
def get_occupancy(atom_group):
  max_partial_occ = 0.
  for atom in atom_group.atoms():
    if (atom.occ > max_partial_occ) and (atom.occ < 1):
      max_partial_occ = atom.occ
  if (max_partial_occ == 0.):
    return max([ atom.occ for atom in atom_group.atoms() ])
  else :
    return max_partial_occ

#-----------------------------------------------------------------------
# GRAPHICS
class rotamer_plot_mixin(graphics.rotarama_plot_mixin):
  def set_labels(self, y_marks=(60,180,300)):
    self.plot.set_xlabel("Chi1")
    self.plot.set_xticks([60,180,300])
    self.plot.set_ylabel("Chi2")
    self.plot.set_yticks(list(y_marks))
    self.plot.grid(True, color="0.75")

class rotamer_plot(data_plots.simple_matplotlib_plot, rotamer_plot_mixin):
  def __init__(self, *args, **kwds):
    data_plots.simple_matplotlib_plot.__init__(self, *args, **kwds)
    rotamer_plot_mixin.__init__(self, *args, **kwds)

def draw_rotamer_plot(rotalyze_data,
                       rotarama_data,
                       residue_name,
                       file_name,
                       show_labels=True):
  points, coords = get_residue_rotamer_data(
    rotalyze_data=rotalyze_data,
    residue_name=residue_name,
    point_type="All")
  p = rotamer_plot()
  title = "Chi1-Chi2 plot for %s" % residue_name
  p.draw_plot(
    stats=rotarama_data,
    contours=[0.05477, 0.14142],
    title=title,
    points=points,
    xyz=coords,
    colormap="Blues")
  p.save_image(file_name)
