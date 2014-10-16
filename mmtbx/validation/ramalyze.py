
from __future__ import division
from mmtbx.validation import residue, validation, atom
from mmtbx.validation import graphics
from iotbx import data_plots
from libtbx import slots_getstate_setstate
import sys

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

res_types = ["general", "glycine", "cis-proline", "trans-proline",
             "pre-proline", "isoleucine or valine"]
res_type_labels = ["General", "Gly", "cis-Pro", "trans-Pro", "pre-Pro",
                   "Ile/Val"]
res_type_plot_labels = ["all non-Pro/Gly residues", "Glycine", "cis-Proline",
  "trans-Proline", "pre-Proline residues", "Ile or Val"]
rama_types = ["OUTLIER", "Allowed", "Favored", "Any", "Allowed/Outlier"]
rama_type_labels = ["Outlier", "Allowed", "Favored", "Any", "Allowed/Outlier"]

class c_alpha (slots_getstate_setstate) :
  """Container class used in the generation of kinemages."""
  __slots__ = ['id_str', 'xyz']
  def __init__ (self, id_str, xyz) :
    self.id_str = id_str
    self.xyz = xyz

class ramachandran (residue) :
  """
  Result class for protein backbone Ramachandran analysis (phenix.ramalyze).
  """
  __rama_attr__ = [
    "res_type",
    "rama_type",
    "score",
    "phi",
    "psi",
    "c_alphas",
  ]
  __slots__ = residue.__slots__ + __rama_attr__

  @staticmethod
  def header () :
    return "%-20s %-12s %10s %6s %-20s" % ("Residue", "Type", "Region", "Score",
      "Phi/Psi")

  def residue_type (self) :
    return res_type_labels[self.res_type]

  def ramalyze_type (self) :
    return rama_types[self.rama_type]

  def as_string (self) :
    return "%-20s %-12s %10s %6.2f %10s" % (self.id_str(), self.residue_type(),
      self.ramalyze_type(), self.score,
      ",".join([ "%.1f" % x for x in [self.phi, self.psi] ]))

  # Backwards compatibility
  def id_str_old (self) :
    return "%s%4s%1s %1s%s" % (self.chain_id, self.resseq, self.icode,
      self.altloc, self.resname)

  def format_old (self) :
    return "%s:%.2f:%.2f:%.2f:%s:%s" % (self.id_str(), self.score,
      self.phi, self.psi, self.ramalyze_type(),
      res_types[self.res_type].capitalize())

  def as_kinemage (self) :
    assert self.is_outlier()
    ram_out = "{%s CA}P %s\n" % (self.c_alphas[0].id_str, "%.3f %.3f %.3f" %
      self.c_alphas[0].xyz)
    ram_out += "{%s CA} %s\n" % (self.c_alphas[1].id_str, "%.3f %.3f %.3f" %
      self.c_alphas[1].xyz)
    ram_out += "{%s CA} %s\n" % (self.c_alphas[2].id_str, "%.3f %.3f %.3f" %
      self.c_alphas[2].xyz)
    return ram_out

  # GUI output
  def as_table_row_phenix (self) :
    return [ self.chain_id, "%s %s" % (self.resname, self.resid),
             self.residue_type(), self.score, self.phi, self.psi ]

class ramachandran_ensemble (residue) :
  """Container for results for an ensemble of residues"""
  __slots__ = ramachandran.__slots__
  def __init__ (self, all_results) :
    self._copy_constructor(all_results[0])
    self.res_type = all_results[0].res_type
    self.rama_type = [ r.rama_type for r in all_results ]
    from scitbx.array_family import flex
    self.phi = flex.double([ r.phi for r in all_results ])
    self.psi = flex.double([ r.psi for r in all_results ])
    self.score = flex.double([ r.score for r in all_results ])

  def phi_min_max_mean (self) :
    return self.phi.min_max_mean()

  def psi_min_max_mean (self) :
    return self.psi.min_max_mean()

  def score_statistics (self) :
    return self.score.min_max_mean()

  def phi_range (self) :
    pass

class ramalyze (validation) :
  """
  Frontend for calculating Ramachandran statistics for a model.  Can directly
  generate the corresponding plots.
  """
  __slots__ = validation.__slots__ + ["out_percent", "fav_percent",
    "n_allowed", "n_favored", "n_type", "_outlier_i_seqs" ]
  program_description = "Analyze protein backbone ramachandran"
  output_header = "residue:score%:phi:psi:evaluation:type"
  gui_list_headers = ["Chain","Residue","Residue type","Score","Phi","Psi"]
  gui_formats = ["%s", "%s", "%s", "%.2f", "%.1f", "%.1f"]
  wx_column_widths = [125]*6

  def get_result_class (self) : return ramachandran

  def __init__ (self,
      pdb_hierarchy,
      outliers_only=False,
      show_errors=False,
      out=sys.stdout,
      quiet=False) :
    validation.__init__(self)
    self.n_allowed = 0
    self.n_favored = 0
    self.n_type = [ 0 ] * 6
    from mmtbx.validation import utils
    import mmtbx.rotamer
    from mmtbx.rotamer import ramachandran_eval
    from scitbx.array_family import flex
    self._outlier_i_seqs = flex.size_t()
    pdb_atoms = pdb_hierarchy.atoms()
    all_i_seqs = pdb_atoms.extract_i_seq()
    if (all_i_seqs.all_eq(0)) :
      pdb_atoms.reset_i_seq()
    use_segids = utils.use_segids_in_place_of_chainids(
      hierarchy=pdb_hierarchy)
    analysis = ""
    output_list = []
    r = ramachandran_eval.RamachandranEval()
    prev_rezes, next_rezes = None, None
    prev_resid = None
    cur_resseq = None
    next_resseq = None
    for model in pdb_hierarchy.models():
      for chain in model.chains():
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
          rezes = construct_complete_residues(residues[i])
          cur_resseq = residue_group.resseq_as_int()
          cur_icode = residue_group.icode.strip()
          if (i > 0):
            #check for insertion codes
            if (cur_resseq == residues[i-1].resseq_as_int()) :
              if (cur_icode == '') and (residues[i-1].icode.strip() == '') :
                continue
            elif (cur_resseq != (residues[i-1].resseq_as_int())+1):
              continue
          if (i < len(residues)-1):
            #find next residue
            if residue_group.resseq_as_int() == \
               residues[i+1].resseq_as_int():
              if (cur_icode == '') and (residues[i+1].icode.strip() == '') :
                continue
            elif residue_group.resseq_as_int() != \
               (residues[i+1].resseq_as_int())-1:
              continue
            next_rezes = construct_complete_residues(residues[i+1])
            next_resid = residues[i+1].resseq_as_int()
          else:
            next_rezes = None
            next_resid = None
          for atom_group in residue_group.atom_groups():
            alt_conf = atom_group.altloc
            if rezes is not None:
              atom_list = rezes.get(alt_conf)
            if prev_rezes is not None:
              prev_atom_list = prev_rezes.get(alt_conf)
              if (prev_atom_list is None):
                prev_keys = sorted(prev_rezes.keys())
                prev_atom_list = prev_rezes.get(prev_keys[0])
            if next_rezes is not None:
              next_atom_list = next_rezes.get(alt_conf)
              if (next_atom_list is None):
                next_keys = sorted(next_rezes.keys())
                next_atom_list = next_rezes.get(next_keys[0])
            phi = get_phi(prev_atom_list, atom_list)
            psi = get_psi(atom_list, next_atom_list)
            coords = get_center(atom_group)
            if (phi is not None and psi is not None):
              res_type = RAMA_GENERAL
              self.n_total += 1
              if (atom_group.resname[0:3] == "GLY"):
                res_type = RAMA_GLYCINE
              elif (atom_group.resname[0:3] == "PRO"):
                is_cis = is_cis_peptide(prev_atom_list, atom_list)
                if is_cis:
                  res_type = RAMA_CISPRO
                else:
                  res_type = RAMA_TRANSPRO
              elif (isPrePro(residues, i)):
                res_type = RAMA_PREPRO
              elif (atom_group.resname[0:3] == "ILE" or \
                    atom_group.resname[0:3] == "VAL"):
                res_type = RAMA_ILE_VAL
              self.n_type[res_type] += 1
              value = r.evaluate(res_types[res_type], [phi, psi])
              ramaType = self.evaluateScore(res_type, value)
              is_outlier = isOutlier(res_type, value)
              c_alphas = None
              # XXX only save kinemage data for outliers
              if is_outlier :
                c_alphas = []
                for atoms in [prev_atom_list, atom_list, next_atom_list] :
                  for a in atoms :
                    if (a.name.strip() == "CA") :
                      a_ = atom(pdb_atom=a)
                      c_alphas.append(c_alpha(
                        id_str=a_.atom_group_id_str(),
                        xyz=a_.xyz))
                assert (len(c_alphas) == 3)
              result = ramachandran(
                chain_id=chain_id,
                resseq=residue_group.resseq,
                icode=residue_group.icode,
                resname=atom_group.resname,
                altloc=atom_group.altloc,
                segid=None, # XXX ???
                phi=phi,
                psi=psi,
                rama_type=ramaType,
                res_type=res_type,
                score=value*100,
                outlier=is_outlier,
                xyz=coords,
                c_alphas=c_alphas)
              if (not outliers_only or is_outlier) :
                self.results.append(result)
              if is_outlier :
                i_seqs = atom_group.atoms().extract_i_seq()
                assert (not i_seqs.all_eq(0))
                self._outlier_i_seqs.extend(i_seqs)
    out_count, out_percent = self.get_outliers_count_and_fraction()
    fav_count, fav_percent = self.get_favored_count_and_fraction()
    self.out_percent = out_percent * 100.0
    self.fav_percent = fav_percent * 100.0

  def write_plots (self, plot_file_base, out) :
    """
    Write a set of six PNG images representing the plots for each residue type.

    :param plot_file_base: file name prefix
    :param out: log filehandle
    """
    from mmtbx.validation import utils
    print >> out, ""
    print >> out, "Creating images of plots..."
    for pos in range(6) :
      stats = utils.get_rotarama_data(
        pos_type=res_types[pos],
        convert_to_numpy_array=True)
      file_label = res_type_labels[pos].replace("/", "_")
      plot_file_name = plot_file_base + "_rama_%s.png" % file_label
      points, coords = self.get_plot_data(position_type=pos)
      draw_ramachandran_plot(
        points=points,
        rotarama_data=stats,
        position_type=pos,
        title=format_ramachandran_plot_title(pos, '*'),
        file_name=plot_file_name)
      print >> out, "  wrote %s" % plot_file_name

  def display_wx_plots (self, parent=None,
      title="MolProbity - Ramachandran plots") :
    import wxtbx.plots.molprobity
    frame = wxtbx.plots.molprobity.ramalyze_frame(
      parent=parent, title=title, validation=self)
    frame.Show()
    return frame

  def show_summary (self, out=sys.stdout, prefix="") :
    print >> out, prefix + 'SUMMARY: %.2f%% outliers (Goal: %s)' % \
      (self.out_percent, self.get_outliers_goal())
    print >> out, prefix + 'SUMMARY: %.2f%% favored (Goal: %s)' % \
      (self.fav_percent, self.get_favored_goal())

  def get_plot_data (self, position_type=RAMA_GENERAL, residue_name="*",
      point_type=RAMALYZE_ANY) :
    assert isinstance(position_type, int) and (0 <= position_type <= 5), \
      position_type
    points, coords = [], []
    for i, residue in enumerate(self.results) :
      if ((residue.res_type == position_type) and
          ((residue_name == '*') or (residue_name == residue.resname))) :
        if ((point_type == RAMALYZE_ANY) or
            (point_type == residue.rama_type) or
            ((residue.rama_type in [RAMALYZE_ALLOWED,RAMALYZE_OUTLIER]) and
             (point_type == RAMALYZE_NOT_FAVORED))) :
          points.append((residue.phi, residue.psi, residue.simple_id(),
            residue.is_outlier()))
          coords.append(residue.xyz)
    return (points, coords)

  def evaluateScore(self, resType, value):
    if (value >= 0.02):
      self.n_favored += 1
      return RAMALYZE_FAVORED
    if (resType == RAMA_GENERAL):
      if (value >= 0.0005):
        self.n_allowed += 1
        return RAMALYZE_ALLOWED
      else:
        self.n_outliers += 1
        return RAMALYZE_OUTLIER
    elif (resType == RAMA_CISPRO):
      if (value >=0.0020):
        self.n_allowed += 1
        return RAMALYZE_ALLOWED
      else:
        self.n_outliers += 1
        return RAMALYZE_OUTLIER
    else:
      if (value >= 0.0010):
        self.n_allowed += 1
        return RAMALYZE_ALLOWED
      else:
        self.n_outliers += 1
        return RAMALYZE_OUTLIER

  def get_outliers_goal(self):
    return "< 0.2%"

  def _get_count_and_fraction (self, res_type) :
    if (self.n_total != 0) :
      count = self.n_type[res_type]
      fraction = float(count) / self.n_total
      return count, fraction
    return 0, 0.

  @property
  def percent_favored (self) :
    n_favored, frac_favored = self.get_favored_count_and_fraction()
    return frac_favored * 100.

  @property
  def percent_allowed (self) :
    n_allowed, frac_allowed = self.get_allowed_count_and_fraction()
    return frac_allowed * 100.

  def get_allowed_count_and_fraction(self):
    if (self.n_total != 0) :
      fraction = self.n_allowed / self.n_total
      return self.n_allowed, fraction
    return 0, 0.

  def get_allowed_goal(self):
    return "> 99.8%"

  def get_favored_count_and_fraction(self):
    if (self.n_total != 0) :
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

  def as_kinemage (self) :
    ram_out = "@subgroup {Rama outliers} master= {Rama outliers}\n"
    ram_out += "@vectorlist {bad Rama Ca} width= 4 color= green\n"
    for rama in self.results :
      if rama.is_outlier() :
        ram_out += rama.as_kinemage()
    return ram_out

  def as_coot_data (self) :
    data = []
    for result in self.results :
      if result.is_outlier() :
        data.append((result.chain_id, result.resid, result.resname,
          result.score, result.xyz))
    return data

def get_matching_atom_group(residue_group, altloc):
  match = None
  if (residue_group != None):
    for ag in residue_group.atom_groups():
      if (ag.altloc == "" and match == None): match = ag
      if (ag.altloc == altloc): match = ag
  return match

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

def get_omega(prev_atoms, atoms):
  import mmtbx.rotamer
  prevCA, prevC, thisN, thisCA = None, None, None, None
  if (prev_atoms is not None):
    for atom in prev_atoms:
      if (atom.name == " CA "): prevCA = atom
      if (atom.name == " C  "): prevC = atom
  if (atoms is not None):
    for atom in atoms:
      if (atom.name == " N  "): thisN = atom
      if (atom.name == " CA "): thisCA = atom
  if (prevCA is not None and prevC is not None and thisN is not None and thisCA is not None):
    return mmtbx.rotamer.omega_from_atoms(prevCA, prevC, thisN, thisCA)

def is_cis_peptide(prev_atoms, atoms):
  omega = get_omega(prev_atoms, atoms)
  if(omega > -30 and omega < 30):
    return True
  else:
    return False

def construct_complete_residues(res_group):
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
      if (not None in [nit, ca, co, oxy]) and (changed) :
        # complete residue backbone found
        complete_dict[ag.altloc] = [nit, ca, co, oxy]
    if len(complete_dict) > 0:
      return complete_dict
  return None

def get_center(ag):
  coords = None
  for atom in ag.atoms():
    if (atom.name == " CA "):
      coords = atom.xyz
  return coords

def isPrePro(residues, i):
  if (i < 0 or i >= len(residues) - 1): return False
  else:
    next = residues[i+1]
    for ag in next.atom_groups():
      if (ag.resname[0:3] == "PRO"): return True
  return False

def isOutlier (resType, value) :
  if (resType == RAMA_GENERAL):
    if (value < 0.0005): return True
    else: return False
  elif (resType == RAMA_CISPRO):
    if (value < 0.0020): return True
    else: return False
  else:
    if (value < 0.0010): return True
    else: return False

#-----------------------------------------------------------------------
# GRAPHICS OUTPUT
def format_ramachandran_plot_title (position_type, residue_type) :
  if (residue_type == '*') :
    title = "Ramachandran plot for " + res_type_plot_labels[position_type]
  else :
    title = "Ramachandran plot for " + residue_type
  return title

class ramachandran_plot_mixin (graphics.rotarama_plot_mixin) :
  extent = [-179,179,-179,179]
  def set_labels (self, y_marks=()) :
    axes = self.plot.get_axes()
    axes.set_xlabel("Phi")
    axes.set_xticks([-120,-60,0,60,120])
    axes.set_ylabel("Psi")
    axes.set_yticks([-120,-60,0,60,120])

class ramachandran_plot (data_plots.simple_matplotlib_plot,
                         ramachandran_plot_mixin) :
  def __init__ (self, *args, **kwds) :
    data_plots.simple_matplotlib_plot.__init__(self, *args, **kwds)
    ramachandran_plot_mixin.__init__(self, *args, **kwds)

def draw_ramachandran_plot (points,
                            rotarama_data,
                            position_type,
                            title,
                            file_name,
                            show_labels=True) :
  p = ramachandran_plot()
  # XXX where do these numbers come from?
  if position_type == RAMA_GENERAL :
    contours = [0.1495, 0.376]
  else :
    contours = [0.2115, 0.376]
  p.draw_plot(
    stats=rotarama_data,
    title=title,
    points=points,
    colormap="Blues",
    contours=contours)
  p.save_image(file_name)
