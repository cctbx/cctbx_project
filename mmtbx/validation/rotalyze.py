
from __future__ import division
from mmtbx.validation import residue, validation
from mmtbx.validation import graphics
from iotbx import data_plots
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
import sys

class rotamer (residue) :
  """
  Result class for protein sidechain rotamer analysis (phenix.rotalyze).
  """
  __rotamer_attr__ = [
    "score",
    "rotamer_name",
    "chi_angles",
    "incomplete",
  ]
  __slots__ = residue.__slots__ + __rotamer_attr__

  @staticmethod
  def header () :
    return "%-20s %8s %6s   %-20s" % ("Residue", "Rotamer", "Score",
      "Chi angles")

  def get_chi1_chi2 (self) :
    if (len(self.chi_angles) < 2) :
      raise ValueError("Less than 2 chi angles for this residue (%s)" %
        self.id_str())
    return self.chi_angles[0], self.chi_angles[1]

  def as_string (self) :
    return "%-20s %8s %6.2f   %-20s" % (self.id_str(), self.rotamer_name,
      self.score, self.format_chi_angles())

  # Backwards compatibility for scripts that expect old rotalyze output
  def id_str_old (self) :
    return "%s%4s%1s %s" % (self.chain_id, self.resseq, self.icode,
      self.altloc + self.resname.strip())

  def format_chi_angles (self, pad=False, sep=",") :
    formatted = []
    for chi in self.chi_angles :
      if pad or (chi is not None) :
        formatted.append(format_value("%.1f", chi,
          replace_none_with="").strip())
    return sep.join(formatted)

  # Old output
  def format_old (self) :
    s_occ = format_value("%.2f", self.occupancy)
    s_score = format_value("%.1f", self.score)
    chis = list(self.chi_angles)
    return "%s:%s:%s:%s:%s" % (self.id_str_old(), s_occ, s_score,
      self.format_chi_angles(pad=True, sep=":"), self.rotamer_name)

class rotalyze (validation) :
  __slots__ = validation.__slots__ + ["out_percent"]
  program_description = "Analyze protein sidechain rotamers"
  output_header = "residue:occupancy:score%:chi1:chi2:chi3:chi4:rotamer"

  def get_result_class (self) : return rotamer

  def __init__ (self, pdb_hierarchy,
      outliers_only=False,
      show_errors=False,
      out=sys.stdout,
      quiet=False) :
    validation.__init__(self)
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    from mmtbx.rotamer.rotamer_eval import RotamerID
    from mmtbx.validation import utils
    sidechain_angles = SidechainAngles(show_errors)
    rotamer_evaluator = rotamer_eval.RotamerEval()
    rotamer_id = rotamer_eval.RotamerID() # loads in the rotamer names
    use_segids = utils.use_segids_in_place_of_chainids(
                   hierarchy=pdb_hierarchy)
    current_rotamers = {}
    for model in pdb_hierarchy.models():
      for chain in model.chains():
        if use_segids:
          chain_id = utils.get_segid_as_chainid(chain=chain)
        else:
          chain_id = chain.id
        for rg in chain.residue_groups():
          all_dict = construct_complete_sidechain(rg)
          #print all_dict
          for atom_group in rg.atom_groups() :
            coords = get_center(atom_group)
            resname = atom_group.resname
            occupancy = get_occupancy(atom_group)
            kwargs = {
              "chain_id" : chain_id,
              "resseq" : rg.resseq,
              "icode" : rg.icode,
              "altloc" : atom_group.altloc,
              "resname" : resname,
              "xyz" : coords,
              "occupancy" : occupancy,
            }
            atom_dict = all_dict.get(atom_group.altloc)
            res_key = get_residue_key(atom_group=atom_group)
            try:
              chis = sidechain_angles.measureChiAngles(
                       atom_group,
                       atom_dict)#.get(conformer.altloc))
            except AttributeError:
              if show_errors:
                kwargs['incomplete'] = True
                result = rotamer(**kwargs)
                print >> out, '%s is missing some sidechain atoms' % \
                  result.id_str()
                self.results.append(result)
              continue
            if (chis is not None):
              if None in chis:
                continue
              value = rotamer_evaluator.evaluate(resname.lower().strip(), chis)
              if value is not None:
                self.n_total += 1
                kwargs['score'] = value * 100
                wrap_chis = rotamer_id.wrap_chis(resname.strip(), chis,
                  symmetry=False)
                sym_chis = wrap_chis[:]
                sym_chis = rotamer_id.wrap_sym(resname.strip(), sym_chis)
                if value < 0.01:
                  self.n_outliers += 1
                  kwargs['outlier'] = True
                  kwargs['rotamer_name'] = "OUTLIER"
                else:
                  kwargs['outlier'] = False
                  kwargs['rotamer_name'] = rotamer_id.identify(resname,
                    wrap_chis)
                while (len(wrap_chis) < 4) :
                  wrap_chis.append(None)
                kwargs['chi_angles'] = wrap_chis
                result = rotamer(**kwargs)
                if (result.is_outlier()) or (not outliers_only) :
                  self.results.append(result)
    out_count, out_percent = self.get_outliers_count_and_fraction()
    self.out_percent = out_percent * 100.0

  def show_summary (self, out=sys.stdout, prefix="") :
    print >> out, prefix + 'SUMMARY: %.2f%% outliers (Goal: %s)' % \
      (self.out_percent, self.get_outliers_goal())

  def get_outliers_goal(self):
    return "< 1%"

  def coot_todo (self):
    return ""

  def get_plot_data (self, residue_name, point_type) :
    assert (point_type in ["All", "Outlier"])
    points = []
    coords = []
    for i, residue in enumerate(self.results) :
      if (residue.resname == residue_name) :
        if ((point_type == "All") or (residue.is_outlier())) :
          chi1, chi2 = residue.get_chi1_chi2()
          points.append((chi1, chi2, residue.simple_id(), residue.is_outlier()))
          coords.append(residue.xyz)
    return (points, coords)

  def display_wx_plots (self, parent=None,
      title="MolProbity - Sidechain Chi1/Chi2 plots") :
    import wxtbx.plots.molprobity
    frame = wxtbx.plots.molprobity.rotalyze_frame(
      parent=parent, title=title, validation=self)
    frame.Show()

  def as_coot_data (self) :
    data = []
    for result in self.results :
      if result.is_outlier() :
        data.append((result.chain_id, result.resid, result.resname,
          result.score, result.xyz))
    return data

def evaluate_rotamer(
    atom_group,
    sidechain_angles,
    rotamer_evaluator,
    rotamer_id,
    all_dict,
    sites_cart=None) :
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
  elif (value < 0.01):
    return 'OUTLIER', chis, value
  else:
    return rotamer_name, chis, value

def split_rotamer_names(rotamer):
  split_rotamer = []
  multi = ""
  if rotamer in ['OUTLIER', 'Cg_exo', 'Cg_endo']:
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
        if (key != cur_label[4:]) :
          raise Sorry("""\
Incompatible identifiers for one or more atoms in a residue:
%s
This is usually caused by atoms with a different segid from the rest of the
residue.  You can use phenix.pdbtools or phenix.pdb_editor to reset the
segid.""" % atom.format_atom_record())
        assert key == cur_label[4:]
  return key

def evaluate_residue(
      residue_group,
      sa,
      r,
      all_dict,
      sites_cart=None):
  is_outlier = False
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
    elif (value < 0.01):
      is_outlier = True
      return is_outlier, value
    else:
      return is_outlier, value

class residue_evaluator (object) :
  def __init__ (self) :
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    self.sa = SidechainAngles(False)
    self.r = rotamer_eval.RotamerEval()

  def evaluate_residue (self, residue_group) :
    all_dict = construct_complete_sidechain(residue_group)
    return evaluate_residue(
      residue_group=residue_group,
      sa=self.sa,
      r=self.r,
      all_dict=all_dict)

  def __call__ (self, *args, **kwds) :
    return self.evaluate_residue(*args, **kwds)

def get_center (residue):
  for atom in residue.atoms():
    if atom.name == " CA ":
      return atom.xyz
  return None

def construct_complete_sidechain(residue_group):
  if (residue_group is not None):
    complete_dict = {}
    atom_dict = {}
    for ag in residue_group.atom_groups():
      for atom in ag.atoms():
        #if atom.name not in atom_dict:
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
def get_occupancy (atom_group) :
  max_partial_occ = 0.
  for atom in atom_group.atoms() :
    if (atom.occ > max_partial_occ) and (atom.occ < 1) :
      max_partial_occ = atom.occ
  if (max_partial_occ == 0.) :
    return max([ atom.occ for atom in atom_group.atoms() ])
  else :
    return max_partial_occ

#-----------------------------------------------------------------------
# GRAPHICS
class rotamer_plot_mixin (graphics.rotarama_plot_mixin) :
  def set_labels (self, y_marks=(60,180,300)) :
    axes = self.plot.get_axes()
    axes.set_xlabel("Chi1")
    axes.set_xticks([60,180,300])
    axes.set_ylabel("Chi2")
    axes.set_yticks(list(y_marks))
    axes.grid(True, color="0.75")

class rotamer_plot (data_plots.simple_matplotlib_plot, rotamer_plot_mixin) :
  def __init__ (self, *args, **kwds) :
    data_plots.simple_matplotlib_plot.__init__(self, *args, **kwds)
    rotamer_plot_mixin.__init__(self, *args, **kwds)

def draw_rotamer_plot (rotalyze_data,
                       rotarama_data,
                       residue_name,
                       file_name,
                       show_labels=True) :
  points, coords = get_residue_rotamer_data(
    rotalyze_data=rotalyze_data,
    residue_name=residue_name,
    point_type="All")
  p = rotamer_plot()
  title = "Chi1-Chi2 plot for %s" % residue_name
  p.draw_plot(
    stats=rotarama_data,
    title=title,
    points=points,
    xyz=coords,
    colormap="Blues",
    contours=None)
  p.save_image(file_name)
