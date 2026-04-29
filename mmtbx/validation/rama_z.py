from __future__ import absolute_import, division, print_function

from libtbx.utils import null_out
from libtbx import easy_pickle
import libtbx.load_env
import iotbx.phil
import iotbx.pdb
from mmtbx.secondary_structure import manager as ss_manager
from mmtbx.secondary_structure import sec_str_master_phil_str
from mmtbx.conformation_dependent_library import generate_protein_threes
from libtbx.str_utils import format_value
from scitbx.array_family import flex
from libtbx import adopt_init_args
from libtbx import group_args
from scitbx.math import linear_interpolation_2d
import numpy as np
import math
import json
import os
import sys

master_phil_str = """
rama_z {

}
"""

class result(object):
  def __init__(self, whole, helix, sheet, loop):
    adopt_init_args(self, locals())

  def as_string(self, prefix=''):
    f = format_value
    p = prefix
    w, h, s, l = self.whole, self.helix, self.sheet, self.loop
    d = "%5.2f"
    i = "%d"
    strs = [
      "\n%sRama-Z values with (uncertainties):"%p,
      "%sInterpretation: poor |Rama-Z| > 3; suspicious 2 < |Rama-Z| < 3; good |Rama-Z| < 2." % p,
      "%sScores below are scaled independently, so they are not related in a simple way." % p,
      "%s  whole: %s (%s), residues: %s"%(p, f(d,w.value),f(d,w.std).strip(),f(i,w.n)),
      "%s  helix: %s (%s), residues: %s"%(p, f(d,h.value),f(d,h.std).strip(),f(i,h.n)),
      "%s  sheet: %s (%s), residues: %s"%(p, f(d,s.value),f(d,s.std).strip(),f(i,s.n)),
      "%s  loop : %s (%s), residues: %s"%(p, f(d,l.value),f(d,l.std).strip(),f(i,l.n))
    ]
    return "\n".join(strs)

  def as_json(self):
    data = {}
    for name, obj in [('whole', self.whole), ('helix', self.helix),
                      ('sheet', self.sheet), ('loop', self.loop)]:
      data[name] = {'value':obj.value,
                    'std':obj.std,
                    'n_residues':obj.n}
    return json.dumps(data, indent=2)

class z_score_mixins(object):
  def _get_mean_one_table(self, table, verbose=False):
    if verbose:
      for i, a in enumerate(table):
        outl=''
        for j, b in enumerate(a):
          outl+='%3d '%b
        print(outl)
    # Origianl paper calc:
    reg_sum = 0
    sq_sum = 0
    for i in table:
      for j in i:
        reg_sum += j
        sq_sum += j*j
    if reg_sum > 0:
      mean = sq_sum / reg_sum
    else:
      mean = 0
    if verbose: print('mean',mean)
    return mean

  def _get_min_max_one_table(self, table):
    min_info=[1e9,0,0]
    max_info=[-1e9,0,0]
    for x, i in enumerate(table):
      for y, j in enumerate(i):
        if j<min_info[0]:
          min_info[0]=j
          min_info[1]=x
          min_info[2]=y
        if j>max_info[0]:
          max_info[0]=j
          max_info[1]=x
          max_info[2]=y
    return min_info, max_info

  def _get_std_one_table(self, table, mean):
    # Origianl paper calc:
    ch = 0
    zn = 0
    for i in table:
      for j in i:
        ch += j * (j-mean)**2
        zn += j
    zn -= 1
    if zn == 0:
      return 0
    std = math.sqrt(ch/zn)
    return std

  def _get_z_score_point_one_table(self, table, phi, psi, step=4):
    half_step=step//2
    vmin = -180+half_step
    if phi < -180+half_step:
      i = -1
      x1 = -180-half_step
      x2 = -180+half_step
    elif phi >= 180-half_step:
      i = -1
      x1 = 180-half_step
      x2 = 180+half_step
    else:
      i = int(abs(-180 + half_step - phi) // step)
      nsteps = abs(vmin - phi) // step
      x1 = vmin + nsteps * step
      x2 = x1 + step

    if psi < -180+half_step:
      j = -1
      y1 = -180-half_step
      y2 = -180+half_step
    elif psi >= 180-half_step:
      j = -1
      y1 = 180-half_step
      y2 = 180+half_step
    else:
      j = int(abs(-180 + half_step - psi) // step)
      nsteps = abs(vmin - psi) // step
      y1 = vmin + nsteps * step
      y2 = y1 + step

    xx = phi
    yy = psi
    v1 = table[i][j]
    v2 = table[i+1][j+1]
    v3 = table[i][j+1]
    v4 = table[i+1][j]

    int_sc = linear_interpolation_2d(x1,y1,x2,y2,v1,v2,v3,v4,xx,yy)
    # print('>>>',phi,psi,x1,y1,x2,y2,v1,v2,v3,v4,int_sc)
    assert phi>=x1 and phi<=x2
    assert psi>=y1 and psi<=y2, f'{y1} {psi} {y2}'
    return int_sc

class rama_z(z_score_mixins):
  def __init__(self, models, log):
    db_path = libtbx.env.find_in_repositories(
        relative_path="chem_data/rama_z/top8000_rama_z_dict.pkl",
        test=os.path.isfile)
    self.log = log
    # this takes ~0.15 seconds, so I don't see a need to cache it somehow.
    self.db = easy_pickle.load(db_path)

    # Python 3 pickle fix
    # =========================================================================
    if sys.version_info.major == 3:
      self.db = easy_pickle.fix_py2_pickle(self.db)
    # =========================================================================

    self.calibration_values = {
        'H': (-0.045355950779513175, 0.1951165524439217),
        'S': (-0.0425581278436754, 0.20068584887814633),
        'L': (-0.018457764754231075, 0.15788374669456848),
        'W': (-0.016806654295023003, 0.12044960331869274)}
    self.residue_counts = {"H": 0, "S": 0, "L":0}
    self.z_score = {"H": None, "S": None, "L":None, 'W': None}
    self.means = {"H": {}, "S": {}, "L": {}}
    self.stds = {"H": {}, "S": {}, "L": {}}

    self.phi_step = 4
    self.psi_step = 4
    self.n_phi_half = 45
    self.n_psi_half = 45

    # this is needed to disable e.g. selection functionality when
    # multiple models are present
    self.n_models = len(models)
    self.res_info = []
    for model in models:
      if model.get_hierarchy().models_size() > 1:
        hierarchy = iotbx.pdb.hierarchy.root()
        m = model.get_hierarchy().models()[0].detached_copy()
        hierarchy.append_model(m)
        asc = hierarchy.atom_selection_cache()
      else:
        hierarchy = model.get_hierarchy()
        asc = model.get_atom_selection_cache()
      sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)
      ss_params = sec_str_master_phil.fetch().extract()
      ss_params.secondary_structure.protein.search_method = "from_ca"
      ss_params.secondary_structure.from_ca_conservative = True

      ssm = ss_manager(hierarchy,
          atom_selection_cache=asc,
          geometry_restraints_manager=None,
          sec_str_from_pdb_file=None,
          # params=None,
          params = ss_params.secondary_structure,
          was_initialized=False,
          mon_lib_srv=None,
          verbose=-1,
          log=null_out(),
          # log=sys.stdout,
          )

      filtered_ann = ssm.actual_sec_str.deep_copy()
      filtered_ann.remove_short_annotations(
          helix_min_len=4, sheet_min_len=4, keep_one_stranded_sheets=True)
      self.helix_sel = asc.selection(filtered_ann.overall_helices_selection())
      self.sheet_sel = asc.selection(filtered_ann.overall_sheets_selection())

      used_atoms = set()
      for three in generate_protein_threes(hierarchy=hierarchy, geometry=None):
        main_residue = three[1]
        phi_psi_atoms = three.get_phi_psi_atoms()
        if phi_psi_atoms is None:
          continue
        phi_atoms, psi_atoms = phi_psi_atoms
        key = [x.i_seq for x in phi_atoms]+[psi_atoms[-1].i_seq]
        key = "%s" % key
        if key not in used_atoms:
          phi, psi = three.get_phi_psi_angles()
          if None in (phi, psi):
            continue
          rkey = three.get_ramalyze_key()
          resname = main_residue.resname
          ss_type = self._figure_out_ss(three)
          self.res_info.append( ["", rkey, resname, ss_type, phi, psi] )
          self.residue_counts[ss_type] += 1
          used_atoms.add(key)
    self.residue_counts["W"] = self.residue_counts["H"] + self.residue_counts["S"] + self.residue_counts["L"]

  def get_residue_counts(self):
    return self.residue_counts

  def get_result(self):
    r  = self.z_score
    if(r["W"] is None): self.get_z_scores() # XXX Odd. This should not be necessary!
    rc = self.get_residue_counts()
    def nov(x,i):
      if(x is None): return None
      else:          return x[i]
    return result(
      whole = group_args(value=nov(r["W"],0), std=nov(r["W"],1), n=rc["W"]),
      helix = group_args(value=nov(r["H"],0), std=nov(r["H"],1), n=rc["H"]),
      sheet = group_args(value=nov(r["S"],0), std=nov(r["S"],1), n=rc["S"]),
      loop  = group_args(value=nov(r["L"],0), std=nov(r["L"],1), n=rc["L"]))

  def get_z_scores(self):
    for k in ['H', 'S', 'L', 'W']:
      if k != 'W':
        element_points = [p for p in self.res_info if p[3] == k]
      else:
        element_points = self.res_info
      c = None
      try:
        c = self._get_z_score_points(element_points)
      except ZeroDivisionError:
        c = None
      if c is not None:
        zs = (c - self.calibration_values[k][0]) / self.calibration_values[k][1]
        zs_std = None
        if len(element_points) > 1:
          zs_std = self._get_z_score_accuracy(element_points, k)
        self.z_score[k] = (zs, zs_std)
    return self.z_score

  def get_detailed_values(self):
    return self.res_info

  def _get_z_score_accuracy(self, points, part):
    scores = []
    values = [x[-1] for x in points]
    sum_values = np.sum(values)
    for v in values:
      s = (sum_values - v)/(len(values)-1)
      scores.append( ( s-self.calibration_values[part][0]) / self.calibration_values[part][1] )
    return np.std(scores) * ((len(points)-1) ** 0.5)

  def get_ss_selections(self):
    if self.n_models > 1:
      raise NotImplementedError
    self.loop_sel = flex.bool([True]*self.helix_sel.size())
    self.loop_sel &= ~self.helix_sel
    self.loop_sel &= ~self.sheet_sel
    return self.helix_sel, self.sheet_sel, self.loop_sel

  def _figure_out_ss(self, three):
    iseq = three.get_phi_psi_atoms()[0][-1].i_seq
    if self.helix_sel[iseq]: return "H"
    elif self.sheet_sel[iseq]: return "S"
    else: return "L"

  def _get_z_score_points(self, points):
    score = 0
    for entry in points:
      if len(entry) == 6:
        sc = self._get_z_score_point(entry)
        entry.append(sc)
      score += entry[-1]
    return score/len(points)

  def _get_z_score_point(self, entry):
    fname, rama_type, resname, ss_type, phi, psi = entry
    phi = round(phi, 10)
    psi = round(psi, 10)
    resname = self._get_resname(rama_type, resname)
    if resname == 'cisPRO':
      ss_type = 'L'
    table = self.db[ss_type][resname]
    int_sc = self._get_z_score_point_one_table(table, phi, psi)
    if self.means[ss_type].get(resname, None) is None:
      self.means[ss_type][resname] = self._get_mean(ss_type, resname)
    if self.stds[ss_type].get(resname, None) is None:
      self.stds[ss_type][resname] = self._get_std(ss_type, resname, self.means[ss_type][resname])
    return (int_sc - self.means[ss_type][resname]) / self.stds[ss_type][resname]

  def _get_mean(self, ss_type, resname):
    return self._get_mean_one_table(self.db[ss_type][resname])

  def _get_std(self, ss_type, resname, mean):
    return self._get_std_one_table(self.db[ss_type][resname], mean)

  def _get_resname(self, rama_type, resname):
    rn = resname
    if resname == "MSE":
      rn = "MET"
    if rama_type == 2:
      rn = 'cisPRO'
    if rama_type == 3:
      rn = 'transPRO'
    if rama_type == 4:
      rn = 'prePRO'
    return rn
