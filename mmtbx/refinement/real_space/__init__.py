from __future__ import division
from libtbx import adopt_init_args
from scitbx.array_family import flex
from scitbx.matrix import rotate_point_around_axis
import time
from cctbx import maptbx
import mmtbx.utils
from mmtbx.rotamer.rotamer_eval import RotamerEval
from libtbx.utils import user_plus_sys_time
import iotbx.pdb
from cctbx import miller
from libtbx.str_utils import format_value

class residue_monitor(object):
  def __init__(self,
               residue,
               selection_sidechain,
               selection_backbone,
               selection_all,
               map_cc_sidechain=None,
               map_cc_backbone=None,
               map_cc_all=None,
               rotamer_status=None):
    adopt_init_args(self, locals())

  def format_info_string(self):
    return "%s %s %6s %6s %6s %9s"%(
      self.residue.resname,
      self.residue.resseq,
      format_value("%6.3f",self.map_cc_all),
      format_value("%6.3f",self.map_cc_backbone),
      format_value("%6.3f",self.map_cc_sidechain),
      self.rotamer_status)

class structure_monitor(object):
  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               target_map_object,
               geometry_restraints_manager):
    adopt_init_args(self, locals())
    self.unit_cell = self.xray_structure.unit_cell()
    self.xray_structure = xray_structure.deep_copy_scatterers()
    #XXX assert xrs are identical !!!
    self.states_collector = mmtbx.utils.states(
      pdb_hierarchy  = self.pdb_hierarchy,
      xray_structure = self.xray_structure)
    self.rotamer_manager = RotamerEval()
    #
    self.map_cc_whole_unit_cell = None
    self.map_cc_around_atoms = None
    self.map_cc_per_atom = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = 0
    self.number_of_rotamer_outliers = 0
    self.residue_monitors = None
    #
    self.initialize()

  def initialize(self):
    #global time_initialize_structure_monitor
    #timer = user_plus_sys_time()
    # residue monitors
    self.residue_monitors = []
    backbone_atoms = ["N","CA","C","O","CB"]
    get_class = iotbx.pdb.common_residue_names_get_class
    sites_cart = self.xray_structure.sites_cart()
    current_map = self.compute_map(xray_structure = self.xray_structure)
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue in chain.only_conformer().residues():
          if(get_class(residue.resname) == "common_amino_acid"):
            residue_i_seqs_backbone  = flex.size_t()
            residue_i_seqs_sidechain = flex.size_t()
            residue_i_seqs_all       = flex.size_t()
            for atom in residue.atoms():
              an = atom.name.strip()
              bb = an in backbone_atoms
              residue_i_seqs_all.append(atom.i_seq)
              if(bb): residue_i_seqs_backbone.append(atom.i_seq)
              else:   residue_i_seqs_sidechain.append(atom.i_seq)
            sca = sites_cart.select(residue_i_seqs_all)
            scs = sites_cart.select(residue_i_seqs_sidechain)
            scb = sites_cart.select(residue_i_seqs_backbone)
            if(scs.size()==0): ccs = None
            else: ccs = self.map_cc(sites_cart=scs, other_map = current_map)
            self.residue_monitors.append(residue_monitor(
              residue             = residue,
              selection_sidechain = residue_i_seqs_sidechain,
              selection_backbone  = residue_i_seqs_backbone,
              selection_all       = residue_i_seqs_all,
              map_cc_sidechain = ccs,
              map_cc_backbone  = self.map_cc(sites_cart=scb, other_map = current_map),
              map_cc_all       = self.map_cc(sites_cart=sca, other_map = current_map),
              rotamer_status   = self.rotamer_manager.evaluate_residue(residue)))
    # globals
    self.map_cc_whole_unit_cell = self.map_cc(other_map = current_map)
    self.map_cc_around_atoms = self.map_cc(other_map = current_map,
      sites_cart = sites_cart)
    self.map_cc_per_atom = self.map_cc(other_map = current_map,
      sites_cart = sites_cart, per_atom = True)
    es = self.geometry_restraints_manager.energies_sites(sites_cart=sites_cart)
    self.rmsd_a = es.angle_deviations()[2]
    self.rmsd_b = es.bond_deviations()[2]
    #self.dist_from_start = flex.mean(self.xray_structure_start.distances(
    #  other = self.xray_structure))
    #assert self.dist_from_start < 1.e-6
    self.number_of_rotamer_outliers = 0
    for r in self.residue_monitors:
      if(r.rotamer_status == "OUTLIER"):
        self.number_of_rotamer_outliers += 1
    #
    #time_initialize_structure_monitor += timer.elapsed()

  def compute_map(self, xray_structure):
    #global time_compute_map
    #timer = user_plus_sys_time()
    mc = self.target_map_object.miller_array.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.target_map_object.crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    #time_compute_map += timer.elapsed()
    return fft_map.real_map_unpadded()

  def map_cc(self, other_map, sites_cart=None, atom_radius=2, per_atom=False):
    #global time_map_cc
    #timer = user_plus_sys_time()
    if(sites_cart is not None):
      if(per_atom):
        result = flex.double()
        for site_cart in sites_cart:
          sel = maptbx.grid_indices_around_sites(
            unit_cell  = self.unit_cell,
            fft_n_real = other_map.focus(),
            fft_m_real = other_map.all(),
            sites_cart = flex.vec3_double([site_cart]),
            site_radii = flex.double(1, atom_radius))
          result.append(flex.linear_correlation(
            x=other_map.select(sel).as_1d(),
            y=self.target_map_object.data.select(sel).as_1d()).coefficient())
      else:
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = self.unit_cell,
          fft_n_real = other_map.focus(),
          fft_m_real = other_map.all(),
          sites_cart = sites_cart,
          site_radii = flex.double(sites_cart.size(), atom_radius))
        result = flex.linear_correlation(
          x=other_map.select(sel).as_1d(),
          y=self.target_map_object.data.select(sel).as_1d()).coefficient()
    else:
      result = flex.linear_correlation(
        x=other_map.as_1d(),
        y=self.target_map_object.data.as_1d()).coefficient()
    #time_map_cc += timer.elapsed()
    return result

  def show(self, prefix=""):
    fmt="%s: cc(cell,atoms): %6.3f %6.3f rmsd(b,a): %6.4f %5.2f moved: %6.3f rota: %3d"
    print fmt%(
      prefix,
      self.map_cc_whole_unit_cell,
      self.map_cc_around_atoms,
      self.rmsd_b,
      self.rmsd_a,
      self.dist_from_start,
      self.number_of_rotamer_outliers)

  def show_residues(self, map_cc_all=None):
    print "resid    CC(sc) CC(bb) CC(sb)   Rotamer"
    for r in self.residue_monitors:
      if(map_cc_all is None or (map_cc_all is not None and map_cc_all>r.map_cc_all)):
        print r.format_info_string()


  def update(self, xray_structure, accept_as_is=False):
    self.xray_structure = xray_structure
    self.pdb_hierarchy.adopt_xray_structure(xray_structure)
    self.initialize()
    #global time_update
    #timer = user_plus_sys_time()
    #unit_cell = xray_structure.unit_cell()
    #current_map = self.compute_map(xray_structure = xray_structure)
    #sites_cart  = xray_structure.sites_cart()
    #sites_cart_ = self.xray_structure.sites_cart()
    #for r in self.residues:
    #  sca = sites_cart.select(r.selection_all)
    #  scs = sites_cart.select(r.selection_sidechain)
    #  scb = sites_cart.select(r.selection_backbone)
    #  map_cc_all       = self.map_cc(sites_cart = sca, map = current_map)
    #  map_cc_sidechain = self.map_cc(sites_cart = scs, map = current_map)
    #  map_cc_backbone  = self.map_cc(sites_cart = scb, map = current_map)
    #  map_value_sidechain = target_simple(target_map=current_map,
    #    sites_cart=scs, unit_cell=self.unit_cell)
    #  map_value_backbone = target_simple(target_map=current_map,
    #    sites_cart=scb, unit_cell=self.unit_cell)
    #  flag = map_cc_all      > r.map_cc_all and \
    #         map_cc_backbone > r.map_cc_backbone and \
    #         map_cc_backbone > map_cc_sidechain
    #  if(r.map_value_backbone > r.map_value_sidechain):
    #    if(map_value_backbone < map_value_sidechain):
    #      flag = False
    #  if(accept_any): flag=True
    #  if(flag):
    #    residue_sites_cart_new = sites_cart.select(r.selection_all)
    #    sites_cart_.set_selected(r.selection_all, residue_sites_cart_new)
    #    r.pdb_hierarchy_residue.atoms().set_xyz(residue_sites_cart_new)
    #    rotamer_status = self.rotamer_manager.evaluate_residue(
    #      residue=r.pdb_hierarchy_residue)
    #    r.rotamer_status = rotamer_status
    #self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_)
    #self.set_globals()
    #self.set_rsr_residue_attributes()
    #self.states_collector.add(sites_cart = sites_cart_)
    #time_update += timer.elapsed()


def negate_map_around_selected_atoms_except_selected_atoms(
      xray_structure,
      map_data,
      iselection,
      selection_within_radius,
      atom_radius,
      select_good=None):
  # XXX time and memory inefficient
  #print "1:",list(iselection)
  sel_around = xray_structure.selection_within(
    radius    = selection_within_radius,
    selection = flex.bool(xray_structure.scatterers().size(), iselection)
    ).iselection()
  #print "2:",list(sel_around)
  negate_selection = flex.size_t(tuple(
    set(sel_around).difference(set(iselection))))

  #print "3:",list(negate_selection)
  if(select_good is not None):
    nsb = flex.bool(xray_structure.scatterers().size(), negate_selection)
    negate_selection = nsb & select_good
  #print "4:",list(negate_selection.iselection())

  #
  sites_cart_p1 = xray_structure.select(negate_selection).expand_to_p1(
      sites_mod_positive=True).sites_cart()
  around_atoms_selections = maptbx.grid_indices_around_sites(
    unit_cell  = xray_structure.unit_cell(),
    fft_n_real = map_data.focus(),
    fft_m_real = map_data.all(),
    sites_cart = sites_cart_p1,
    site_radii = flex.double(sites_cart_p1.size(), atom_radius))
  sel_ = flex.bool(size=map_data.size(), iselection=around_atoms_selections)
  sel_.reshape(map_data.accessor())
  md = map_data.deep_copy().set_selected(sel_, -2)
  md = md.set_selected(~sel_, 1)
  return map_data*md

class score(object):
  def __init__(self,
               unit_cell,
               target_map,
               rotamer_eval,
               residue,
               vector = None,
               slope_decrease_factor = 3,
               tmp = None):
    adopt_init_args(self, locals())
    self.target = None
    self.sites_cart = None

  def compute_target(self, sites_cart, selection=None):
    sites_frac = self.unit_cell.fractionalize(sites_cart)
    result = 0
    if(selection is None):
      for site_frac in sites_frac:
        result += self.target_map.eight_point_interpolation(site_frac)
    else:
      for sel in selection:
        result += self.target_map.eight_point_interpolation(sites_frac[sel])
    return result

  def update(self, sites_cart, selection=None, tmp=None):
    target = self.compute_target(sites_cart = sites_cart, selection=selection)
    slope = None
    if(self.vector is not None):
      slope = self.get_slope(sites_cart = sites_cart)
      #slope, y, p = self.get_slope(sites_cart = sites_cart)
      #print "   ",target, list(y), slope, p, flex.linear_regression(
      #  flex.double(tuple(self.vector)),y).slope()
      #slope=True
    if(self.target is not None):
      if(target > self.target and (slope is None or (slope is not None and slope))):
        self.residue.atoms().set_xyz(sites_cart)
        rotameric = self.rotamer_eval.evaluate_residue(residue = self.residue)
        if(rotameric):
          self.target = target
          self.sites_cart = sites_cart
          self.tmp = tmp,slope,target
    else:
      self.residue.atoms().set_xyz(sites_cart)
      rotameric = self.rotamer_eval.evaluate_residue(residue = self.residue)
      if(rotameric):
        self.target = target
        self.sites_cart = sites_cart
        self.tmp = tmp, slope, target

  def get_slope(self, sites_cart):
    sites_frac = self.unit_cell.fractionalize(sites_cart)
    y = flex.double()
    for v in self.vector:
      y.append(self.target_map.eight_point_interpolation(sites_frac[v]))
    result = True
    #for i, y_ in enumerate(y):
    #  if(i<y.size()-1):
    #    d = y[i+1]-y[i]
    #    p = d/((abs(y[i+1])+abs(y[i]))/2)
    #    if(p>0 and abs(d)>(abs(y[i+1])+abs(y[i]))/2/self.slope_decrease_factor):
    #      result = False
    #      break
    #return result, y, (d,(abs(y[i+1])+abs(y[i]))/2/self.slope_decrease_factor)
    ynew = flex.double()
    for y_ in y:
      if(y_<1): ynew.append(0)
      else:     ynew.append(1)
    found = False
    for y_ in ynew:
      if(not found and y_==0): found=True
      if(found and y_==1): return False#,list(y),list(ynew)
    return True#,list(y),list(ynew)

  def reset_with(self, sites_cart, selection=None):
    self.target = self.compute_target(sites_cart = sites_cart,
      selection = selection)
    self.sites_cart = sites_cart

def torsion_search(
      clusters,
      scorer,
      sites_cart,
      start = -50,
      stop  = 50,
      step  = 1):
  def generate_range(start, stop, step):
    assert abs(start) <= abs(stop)
    inc = start
    result = []
    while abs(inc) <= abs(stop):
      result.append(inc)
      inc += step
    return result
  for i_cl, cl in enumerate(clusters):
    if(i_cl == 0):
      scorer.reset_with(sites_cart=sites_cart.deep_copy(),
        selection=cl.selection)
    else:
      scorer.reset_with(sites_cart=scorer.sites_cart.deep_copy(),
        selection=cl.selection)
    sites_cart_ = scorer.sites_cart.deep_copy()
    for angle_deg in generate_range(start = start, stop = stop, step = step):
      xyz_moved = sites_cart_.deep_copy()
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = sites_cart_[cl.axis[0]],
          axis_point_2 = sites_cart_[cl.axis[1]],
          point        = sites_cart_[atom],
          angle        = angle_deg, deg=True)
        xyz_moved[atom] = new_xyz
      scorer.update(sites_cart = xyz_moved, selection = cl.selection)
  return scorer

def torsion_search_nested(
      clusters,
      scorer,
      sites_cart,
      start = -6,
      stop  = 6):
  n_angles = len(clusters)
  #nested_loop = flex.nested_loop(begin=[start]*n_angles, end=[stop]*n_angles,
  #  open_range=False)
  #if(n_angles==3):
  #  r1 = [-5,-5,-5]
  #  r2 = [5,5,5]
  #else: return scorer
  x = 5
  #print n_angles
  r1 = range(-x, -n_angles-x, -1)
  r2 = range(x, n_angles+x)
  #print list(r1), list(r2)
  nested_loop = flex.nested_loop(begin=r1, end=r2, open_range=False)
  selection = clusters[0].atoms_to_rotate
  scorer.reset_with(sites_cart = sites_cart, selection = selection)
  for angles in nested_loop:
    xyz_moved = sites_cart.deep_copy()
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle, deg=True)
        xyz_moved[atom] = new_xyz
    scorer.update(sites_cart = xyz_moved, selection = selection)
  return scorer
