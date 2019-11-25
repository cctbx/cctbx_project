from __future__ import absolute_import, division, print_function

from libtbx.utils import null_out
from libtbx import easy_pickle
import libtbx.load_env
import iotbx.phil
from mmtbx.secondary_structure import manager as ss_manager
from mmtbx.secondary_structure import sec_str_master_phil_str
from mmtbx.conformation_dependent_library import generate_protein_threes
from scitbx.array_family import flex
from scipy import interpolate
import numpy as np
import copy
import math
import os

master_phil_str = """
rama_z {

}
"""

class rama_z(object):
  def __init__(self, model, log):
    db_path = libtbx.env.find_in_repositories(
        relative_path="chem_data/rama_z/top8000_rama_z_dict.pkl",
        test=os.path.isfile)
    rmsd_path = libtbx.env.find_in_repositories(
        relative_path="chem_data/rama_z/rmsd.pkl",
        test=os.path.isfile)
    self.log = log
    # this takes ~0.15 seconds, so I don't see a need to cache it somehow.
    self.db = easy_pickle.load(db_path)
    self.rmsd_estimator = easy_pickle.load(rmsd_path)
    self.calibration_values = {
        'H': (-0.045355950779513175, 0.1951165524439217),
        'S': (-0.0425581278436754, 0.20068584887814633),
        'L': (-0.018457764754231075, 0.15788374669456848),
        'W': (-0.016806654295023003, 0.12044960331869274)}
    self.residue_counts = {"H": 0, "S": 0, "L":0}
    self.z_score = {"H": None, "S": None, "L":None, 'W': None}
    self.interpolation_fs = {"H": {}, "S": {}, "L": {}}
    self.means = {"H": {}, "S": {}, "L": {}}
    self.stds = {"H": {}, "S": {}, "L": {}}

    self.phi_step = 4
    self.psi_step = 4
    self.n_phi_half = 45
    self.n_psi_half = 45

    self.res_info = []
    asc = model.get_atom_selection_cache()
    sec_str_master_phil = iotbx.phil.parse(sec_str_master_phil_str)
    ss_params = sec_str_master_phil.fetch().extract()
    ss_params.secondary_structure.protein.search_method = "from_ca"
    ss_params.secondary_structure.from_ca_conservative = True

    self.ssm = ss_manager(model.get_hierarchy(),
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

    filtered_ann = self.ssm.actual_sec_str.deep_copy()
    filtered_ann.remove_short_annotations(
        helix_min_len=4, sheet_min_len=4, keep_one_stranded_sheets=True)
    self.helix_sel = asc.selection(filtered_ann.overall_helices_selection())
    self.sheet_sel = asc.selection(filtered_ann.overall_sheets_selection())

    used_atoms = set()
    for three in generate_protein_threes(hierarchy=model.get_hierarchy(), geometry=None):
      main_residue = three[1]
      phi_psi_atoms = three.get_phi_psi_atoms()
      if phi_psi_atoms is None:
        continue
      phi_atoms, psi_atoms = phi_psi_atoms
      key = [x.i_seq for x in phi_atoms]+[psi_atoms[-1].i_seq]
      key = "%s" % key
      if key not in used_atoms:
        phi, psi = three.get_phi_psi_angles()
        rkey = three.get_ramalyze_key()
        resname = main_residue.resname
        ss_type = self._figure_out_ss(three)
        self.res_info.append( ["", rkey, resname, ss_type, phi, psi] )
        self.residue_counts[ss_type] += 1
        used_atoms.add(key)
    self.residue_counts["W"] = self.residue_counts["H"] + self.residue_counts["S"] + self.residue_counts["L"]
    for i in self.res_info:
      print(i, file=self.log)

  def get_residue_counts(self):
    return self.residue_counts

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
        zs_std = self._get_z_score_accuracy(element_points, k)
        self.z_score[k] = (zs, zs_std)
    return self.z_score

  def _get_z_score_accuracy(self, points, part, n_shuffles=50, percent_to_keep=50):
    return np.interp(len(points), self.rmsd_estimator[0], self.rmsd_estimator[1])
    # tmp = copy.deepcopy(points)
    # scores = []
    # n_res = int(len(tmp) * percent_to_keep / 100)
    # if n_res == len(tmp):
    #   n_res -= 1
    # for i in range(n_shuffles):
    #   np.random.shuffle(tmp)
    #   c = self._get_z_score_points(tmp[:n_res])
    #   if c is not None:
    #     c = (c - self.calibration_values[part][0]) / self.calibration_values[part][1]
    #     scores.append(c)
    #   c = self._get_z_score_points(tmp[n_res:])
    #   if c is not None:
    #     c = (c - self.calibration_values[part][0]) / self.calibration_values[part][1]
    #     scores.append(c)
    # return np.std(scores)

  def get_ss_selections(self):
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
    # if len(points) < 10:
    #   return None
    score = 0
    for entry in points:
      if len(entry) == 6:
        sc = self._get_z_score_point(entry)
        entry.append(sc)
      score += entry[-1]
    return score/len(points)

  def _get_z_score_point(self, entry):
    fname, rama_type, resname, ss_type, phi, psi = entry
    resname = self._get_resname(rama_type, resname)
    if resname == 'cisPRO':
      ss_type = 'L'

    if self.interpolation_fs[ss_type].get(resname, None) is None:
      self.interpolation_fs[ss_type][resname] = self._set_interpolation_f(self.db[ss_type][resname])
    int_sc = self.interpolation_fs[ss_type][resname]([phi], [psi])[0]

    if self.means[ss_type].get(resname, None) is None:
      self.means[ss_type][resname] = self._get_mean(ss_type, resname)
    if self.stds[ss_type].get(resname, None) is None:
      self.stds[ss_type][resname] = self._get_std(ss_type, resname, self.means[ss_type][resname])
    return (int_sc - self.means[ss_type][resname]) / self.stds[ss_type][resname]

  def _get_mean(self, ss_type, resname):
    # return np.average(self.g)

    # Nonzero regular:
    # nz = []
    # for i in self.g:
    #   for j in i:
    #     if j != 0:
    #       nz.append(j)
    # self.mean = np.average(nz)
    # Paper calc:
    reg_sum = 0
    sq_sum = 0
    # calculate as in paper
    for i in self.db[ss_type][resname]:
      for j in i:
        reg_sum += j
        sq_sum += j*j
    if reg_sum > 0:
      mean = sq_sum / reg_sum
    else:
      mean = 0
    return mean

  def _get_std(self, ss_type, resname, mean):
    # return np.std(self.g)

    # Nonzero regular:
    # nz = []
    # for i in self.g:
    #   for j in i:
    #     if j != 0:
    #       nz.append(j)
    # self.std = np.std(nz)
    # Paper calc:
    ch = 0
    zn = 0
    for i in self.db[ss_type][resname]:
      for j in i:
        ch += j * (j-mean)**2
        zn += j
    zn -= 1
    if zn == 0:
      return 0
    std = math.sqrt(ch/zn)
    return std

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

  def _calc_ij(self, p):
    i = int(p[0] // self.phi_step) + self.n_phi_half
    j = int(p[1] // self.psi_step) + self.n_psi_half
    # assert  0 <= i < len(self.g), i
    # assert  0 <= j < len(self.g[1]), j
    return i,j

  def _set_interpolation_f(self, grid):
    x = range(-180-self.phi_step // 2, 180 +self.phi_step + self.phi_step //2, self.phi_step)
    y = range(-180-self.psi_step // 2, 180 +self.psi_step + self.psi_step //2, self.psi_step)
    z = []
    # print "x,y", x, y
    for i in range(len(grid)+2):
      z.append([0] * (len(grid)+2) )
    for i in range(len(z)):
      for j in range(len(z)):
        # figure out where to get value
        ii = i-1
        jj = j-1
        if i == 0:
          ii = len(grid)-1
        if i == len(z) - 1:
          ii = 0
        if j == 0:
          jj = len(grid)-1
        if j == len(z) - 1:
          jj = 0
        z[i][j] = grid[jj][ii]
    return interpolate.interp2d(x,y,z, kind='linear')

