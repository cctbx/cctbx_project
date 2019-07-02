from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math, sys
from mmtbx.utils import rotatable_bonds
import mmtbx.geometry_restraints.torsion_restraints.utils as torsion_utils
from scitbx.matrix import rotate_point_around_axis
import mmtbx.monomer_library.server
from mmtbx.validation import rotalyze
import iotbx.phil
from libtbx import adopt_init_args
from six.moves import zip

torsion_search_params_str = """\
torsion_search
  .style = box auto_align
{
  min_angle_between_solutions = 5
    .type = float
    .short_caption = Min. angle between solutions
  range_start = -40
    .type = float
  range_stop = 40
    .type = float
  step = 2
    .type = float
}
"""

def generate_range(start, stop, step):
  assert abs(start) <= abs(stop)
  inc = start
  result = []
  while abs(inc) <= abs(stop):
    result.append(inc)
    inc += step
  return result

def target(sites_cart_residue, unit_cell, m):
  sites_frac_residue = unit_cell.fractionalize(sites_cart_residue)
  result = 0
  for rsf in sites_frac_residue:
    result += m.eight_point_interpolation(rsf)
  return result

def torsion_search_params():
  return iotbx.phil.parse(input_string = torsion_search_params_str)

class rotamer_evaluator(object):
  def __init__(self, sites_cart_start,
                     unit_cell,
                     two_mfo_dfc_map,
                     mfo_dfc_map):
    adopt_init_args(self, locals())
    t1 = target(self.sites_cart_start, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(self.sites_cart_start, self.unit_cell, self.mfo_dfc_map)
    self.t1 = t1
    self.t2 = t2
    self.t_start = t1+t2
    self.t_best = self.t_start
    self.t1_start = t1
    self.t1_best = self.t1_start
    self.t2_start = t2
    self.t2_best = self.t2_start

  def is_better(
                self,
                sites_cart,
                percent_cutoff=0.0,
                verbose=False):
    t1 = target(sites_cart, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(sites_cart, self.unit_cell, self.mfo_dfc_map)
    t = t1+t2#*3 # XXX very promising thing to do, but reaaly depends on resolution
    result = False
    size = sites_cart.size()
    if 1:
      if(t > self.t_best):
        if percent_cutoff > 0.0 and self.t_best > 0.0:
          percent = (t - self.t_best) / self.t_best
          if percent < percent_cutoff:
            return False
        if((t2 > 0 and self.t2_best > 0 and t2 > self.t2_best) or
           (t2 < 0 and self.t2_best < 0 and abs(t2)<abs(self.t2_best)) or
           (t2 > 0 and self.t2_best < 0)):
          result = True
          self.t2_best = t2
          self.t1_best = t1
          self.t_best = t
    return result

def torsion_search(residue_evaluator,
                   cluster_evaluators,
                   axes_and_atoms_to_rotate,
                   rotamer_sites_cart,
                   rotamer_id_best,
                   residue_sites_best,
                   params = None,
                   rotamer_id = None,
                   include_ca_hinge = False):
  if(params is None):
    params = torsion_search_params().extract().torsion_search
    params.range_start = 0
    params.range_stop = 360
    params.step = 1.0
  rotamer_sites_cart_ = rotamer_sites_cart.deep_copy()
  n_clusters = len(axes_and_atoms_to_rotate)
  c_counter = 0
  for cluster_evaluator, aa in zip(cluster_evaluators,axes_and_atoms_to_rotate):
    #account for CA hinge at beginning of search
    if include_ca_hinge and c_counter == 0:
      cur_range_start = -6.0
      cur_range_stop = 6.0
    else:
      cur_range_start = params.range_start
      cur_range_stop = params.range_stop
    c_counter += 1
    axis = aa[0]
    atoms = aa[1]
    angle_deg_best = None
    angle_deg_good = None
    for angle_deg in generate_range(start = cur_range_start, stop =
                                    cur_range_stop, step = params.step):
      if(c_counter != n_clusters):
        if include_ca_hinge and c_counter == 1:
          new_xyz = flex.vec3_double()
          for atom in atoms:
            new_xyz.append(rotate_point_around_axis(
              axis_point_1 = rotamer_sites_cart[axis[0]],
              axis_point_2 = rotamer_sites_cart[axis[1]],
              point  = rotamer_sites_cart[atom],
              angle = angle_deg, deg=True))
        else:
          point_local = rotamer_sites_cart[atoms[0]]
          new_xyz = flex.vec3_double([rotate_point_around_axis(
            axis_point_1 = rotamer_sites_cart[axis[0]],
            axis_point_2 = rotamer_sites_cart[axis[1]],
            point  = point_local,
            angle = angle_deg, deg=True)])
      else:
        new_xyz = flex.vec3_double()
        for atom in atoms:
          new_xyz.append(rotate_point_around_axis(
            axis_point_1 = rotamer_sites_cart[axis[0]],
            axis_point_2 = rotamer_sites_cart[axis[1]],
            point  = rotamer_sites_cart[atom],
            angle = angle_deg, deg=True))
      if(cluster_evaluator.is_better(sites_cart = new_xyz)):
        if(angle_deg_best is not None and
           abs(abs(angle_deg_best)-abs(angle_deg))>
           params.min_angle_between_solutions):
          angle_deg_good = angle_deg_best
        angle_deg_best = angle_deg
    if(angle_deg_best is not None):
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = rotamer_sites_cart[axis[0]],
          axis_point_2 = rotamer_sites_cart[axis[1]],
          point  = rotamer_sites_cart[atom],
          angle = angle_deg_best, deg=True)
        rotamer_sites_cart[atom] = new_xyz
    if(angle_deg_good is not None):
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = rotamer_sites_cart_[axis[0]],
          axis_point_2 = rotamer_sites_cart_[axis[1]],
          point  = rotamer_sites_cart_[atom],
          angle = angle_deg_best, deg=True)
        rotamer_sites_cart_[atom] = new_xyz
  for rsc in [rotamer_sites_cart, rotamer_sites_cart_]:
    if(residue_evaluator.is_better(sites_cart = rsc)):
      rotamer_id_best = rotamer_id
      residue_sites_best = rsc.deep_copy()
  return residue_sites_best, rotamer_id_best

class manager(object):

  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               range_start=-10.0,
               range_stop=10.0,
               step=1.0,
               min_angle_between_solutions=0.5,
               name_hash=None,
               selection=None,
               log=None):
    if(log is None): log = sys.stdout
    self.log = log
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.unit_cell = xray_structure.unit_cell()
    self.exclude_free_r_reflections = False
    self.torsion_params = torsion_search_params().extract().torsion_search
    self.torsion_params.range_start = range_start
    self.torsion_params.range_stop = range_stop
    self.torsion_params.step = step
    self.torsion_params.min_angle_between_solutions = \
      min_angle_between_solutions
    self.selection=selection
    if self.selection is None:
      sites_cart = pdb_hierarchy.atoms().extract_xyz()
      self.selection = flex.bool(len(sites_cart), True)
    self.c_alpha_hinges = torsion_utils.get_c_alpha_hinges(
                            pdb_hierarchy=pdb_hierarchy,
                            xray_structure=xray_structure,
                            selection=self.selection)
    self.name_hash = name_hash
    if self.name_hash is None:
      self.name_hash = build_name_hash(pdb_hierarchy)
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    self.sa = SidechainAngles(False)
    self.rotamer_id = rotamer_eval.RotamerID()
    self.rotamer_evaluator = rotamer_eval.RotamerEval(mon_lib_srv=self.mon_lib_srv)
    self.target_map_data = None
    self.residual_map_data = None

  def prepare_map(
        self,
        fmodel):
    target_map_data, residual_map_data = \
      torsion_utils.prepare_map(fmodel=fmodel)
    self.target_map_data = target_map_data
    self.residual_map_data = residual_map_data

  def search(
        self,
        atom_group,
        all_dict,
        m_chis,
        r_chis,
        rotamer,
        sites_cart_moving,
        xray_structure,
        key):
    include_ca_hinge = False
    axis_and_atoms_to_rotate, tardy_labels= \
      rotatable_bonds.axes_and_atoms_aa_specific(
          residue=atom_group,
          mon_lib_srv=self.mon_lib_srv,
          remove_clusters_with_all_h=True,
          include_labels=True,
          log=None)
    if (axis_and_atoms_to_rotate is None):
      print("Skipped %s rotamer (TARDY error)" % key, file=self.log)
      return False
    assert len(m_chis) == len(r_chis)
    #exclude H-only clusters if necessary
    while len(axis_and_atoms_to_rotate) > len(m_chis):
      axis_and_atoms_to_rotate = \
        axis_and_atoms_to_rotate[:-1]
    assert len(m_chis) == len(axis_and_atoms_to_rotate)
    counter = 0
    residue_iselection = atom_group.atoms().extract_i_seq()
    cur_ca = None
    ca_add = None
    ca_axes = []
    for atom in atom_group.atoms():
      if atom.name == " CA ":
        cur_ca = atom.i_seq
    if cur_ca is not None:
      cur_c_alpha_hinges = self.c_alpha_hinges.get(cur_ca)
      if cur_c_alpha_hinges is not None:
        residue_length = len(tardy_labels)
        for ca_pt in cur_c_alpha_hinges[0]:
          residue_iselection.append(ca_pt)
          tardy_labels.append(self.name_hash[ca_pt][0:4])
        for bb_pt in cur_c_alpha_hinges[1]:
          residue_iselection.append(bb_pt)
          tardy_labels.append(self.name_hash[bb_pt][0:4])
        end_pts = (residue_length, residue_length+1)
        group = []
        for i, value in enumerate(tardy_labels):
          if i not in end_pts:
            group.append(i)
        ca_add = [end_pts, group]
        ca_axes.append(ca_add)
        for ax in axis_and_atoms_to_rotate:
          ca_axes.append(ax)
    sites_cart_residue = \
      sites_cart_moving.select(residue_iselection)
    sites_cart_residue_start = sites_cart_residue.deep_copy()
    selection = flex.bool(
                  len(sites_cart_moving),
                  residue_iselection)
    rev_first_atoms = []

    rev_start = rotamer_evaluator(
      sites_cart_start = sites_cart_residue_start,
      unit_cell        = self.unit_cell,
      two_mfo_dfc_map  = self.target_map_data,
      mfo_dfc_map      = self.residual_map_data)

    sidechain_only_iselection = flex.size_t()
    for i_seq in residue_iselection:
      atom_name = self.name_hash[i_seq][0:4]
      if atom_name not in [' N  ', ' CA ', ' C  ', ' O  ']:
        sidechain_only_iselection.append(i_seq)
    sites_cart_sidechain = \
      sites_cart_moving.select(sidechain_only_iselection)
    sites_frac_residue = self.unit_cell.fractionalize(sites_cart_sidechain)
    sigma_cutoff = 1.0
    sigma_residue = []
    for rsf in sites_frac_residue:
      if self.target_map_data.eight_point_interpolation(rsf) < sigma_cutoff:
        sigma_residue.append(False)
      else:
        sigma_residue.append(True)
    sigma_count_start = 0
    for sigma_state in sigma_residue:
      if sigma_state:
        sigma_count_start += 1
      else:
        break

    for aa in axis_and_atoms_to_rotate:
      axis = aa[0]
      atoms = aa[1]
      new_xyz = flex.vec3_double()
      angle_deg = r_chis[counter] - m_chis[counter]
      #skip angle rotations that are close to zero
      if math.fabs(angle_deg) < 0.01:
        counter += 1
        continue
      if angle_deg < 0:
        angle_deg += 360.0
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
                    axis_point_1=sites_cart_residue[axis[0]],
                    axis_point_2=sites_cart_residue[axis[1]],
                    point=sites_cart_residue[atom],
                    angle=angle_deg, deg=True)
        sites_cart_residue[atom] = new_xyz
      counter += 1

    #***** TEST *****
    sites_cart_moving.set_selected(
      residue_iselection, sites_cart_residue)
    cur_rotamer, cur_chis, cur_value = rotalyze.evaluate_rotamer(
      atom_group=atom_group,
      sidechain_angles=self.sa,
      rotamer_evaluator=self.rotamer_evaluator,
      rotamer_id=self.rotamer_id,
      all_dict=all_dict,
      sites_cart=sites_cart_moving)
    assert rotamer == cur_rotamer
    #****************

    if len(ca_axes) == 0:
      eval_axes = axis_and_atoms_to_rotate
    else:
      eval_axes = ca_axes
      include_ca_hinge = True
    for i_aa, aa in enumerate(eval_axes):
      if(i_aa == len(eval_axes)-1):
        sites_aa = flex.vec3_double()
        for aa_ in aa[1]:
          sites_aa.append(sites_cart_residue[aa_])
      elif i_aa == 0 and include_ca_hinge:
        sites_aa = flex.vec3_double()
        for aa_ in aa[1]:
          sites_aa.append(sites_cart_residue[aa_])
      else:
        sites_aa = flex.vec3_double([sites_cart_residue[aa[1][0]]])
      rev_i = rotamer_evaluator(
        sites_cart_start = sites_aa,
        unit_cell        = self.unit_cell,
        two_mfo_dfc_map  = self.target_map_data,
        mfo_dfc_map      = self.residual_map_data)
      rev_first_atoms.append(rev_i)

    rev = rotamer_evaluator(
      sites_cart_start = sites_cart_residue,
      unit_cell        = self.unit_cell,
      two_mfo_dfc_map  = self.target_map_data,
      mfo_dfc_map      = self.residual_map_data)

    residue_sites_best = sites_cart_residue.deep_copy()
    residue_sites_best, rotamer_id_best = \
      torsion_search(
        residue_evaluator=rev,
        cluster_evaluators=rev_first_atoms,
        axes_and_atoms_to_rotate=eval_axes,
        rotamer_sites_cart=sites_cart_residue,
        rotamer_id_best=rotamer,
        residue_sites_best=residue_sites_best,
        params = self.torsion_params,
        rotamer_id = rotamer,
        include_ca_hinge = include_ca_hinge)
    sites_cart_moving.set_selected(
        residue_iselection, residue_sites_best)
    xray_structure.set_sites_cart(sites_cart_moving)
    cur_rotamer, cur_chis, cur_value = rotalyze.evaluate_rotamer(
      atom_group=atom_group,
      sidechain_angles=self.sa,
      rotamer_evaluator=self.rotamer_evaluator,
      rotamer_id=self.rotamer_id,
      all_dict=all_dict,
      sites_cart=sites_cart_moving)
    rotamer_match = (cur_rotamer == rotamer)
    if rev_start.is_better(sites_cart=residue_sites_best,
                           percent_cutoff=0.15, verbose=True) and \
       rotamer_match:
      sidechain_only_iselection = flex.size_t()
      for i_seq in residue_iselection:
        atom_name = self.name_hash[i_seq][0:4]
        if atom_name not in [' N  ', ' CA ', ' C  ', ' O  ',
                             ' OXT', ' H  ', ' HA ']:
          sidechain_only_iselection.append(i_seq)
      selection = flex.bool(
                    len(sites_cart_moving),
                    sidechain_only_iselection)
      selection_within = xray_structure.selection_within(
      radius    = 1.0,
      selection = selection)
      #check for bad steric clashes
      created_clash = False
      for i, state in enumerate(selection_within):
        if state:
          if i not in sidechain_only_iselection:
            #print >> self.log, "atom clash: ", self.name_hash[i]
            created_clash = True
      if created_clash:
        sites_cart_moving.set_selected(
          residue_iselection, sites_cart_residue_start)
        xray_structure.set_sites_cart(sites_cart_moving)
        return False

      sidechain_only_iselection = flex.size_t()
      for i_seq in residue_iselection:
        atom_name = self.name_hash[i_seq][0:4]
        if atom_name not in [' N  ', ' CA ', ' C  ', ' O  ']:
          sidechain_only_iselection.append(i_seq)
      sites_cart_sidechain = \
        sites_cart_moving.select(sidechain_only_iselection)
      sites_frac_residue = self.unit_cell.fractionalize(sites_cart_sidechain)
      sigma_cutoff = 1.0
      sigma_residue = []
      for rsf in sites_frac_residue:
        if self.target_map_data.eight_point_interpolation(rsf) < sigma_cutoff:
          sigma_residue.append(False)
        else:
          sigma_residue.append(True)
      sigma_count = 0
      for sigma_state in sigma_residue:
        if sigma_state:
          sigma_count += 1
        else:
          break
      if sigma_count < sigma_count_start:
        sites_cart_moving.set_selected(
          residue_iselection, sites_cart_residue_start)
        xray_structure.set_sites_cart(sites_cart_moving)
        return False
      return True
    else:
      sites_cart_moving.set_selected(
        residue_iselection, sites_cart_residue_start)
      xray_structure.set_sites_cart(sites_cart_moving)
      return False
