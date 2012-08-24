from __future__ import division
import iotbx.phil
from cctbx.maptbx import real_space_target_and_gradients
from libtbx import adopt_init_args
from libtbx.str_utils import format_value
import scitbx.lbfgs
from cctbx import maptbx
from cctbx.array_family import flex
from mmtbx import utils
from mmtbx.rotamer.rotamer_eval import RotamerEval
from cctbx import miller
from libtbx.test_utils import approx_equal
import mmtbx.monomer_library
from libtbx.utils import user_plus_sys_time

master_params_str = """\
real_space_refinement
  .short_caption = Real-space refinement
  .style = menu_item auto_align
{
    mode = simple diff_map *lockit
      .type = choice(multi=False)
      .help = Real space refinement method (diff_map is much slower but might \
              have larger convergence radius especially at low resolution)
    target_map_name = 2mFo-DFc
      .type = str
      .help = Target map type to refine against
    model_map_name = Fc
      .type = str
      .help = Current model map type (used only if diff_map is selected)
    grid_resolution_factor = 1./4
      .type = float
      .help = Defines coarseness of the map
    target_weights = value grid_search *gradients_ratio
      .type = choice(multi=False)
      .help = Method for relative target weight determination
    real_space_target_weight = None
      .type = float
      .help = Weight for the real-space target (used only if diff_map is selected)
    restraints_target_weight = None
      .type = float
      .help = Weight for the restraints target (used only if diff_map is selected)
    number_of_cycles = 1
      .type = int
      .help = Number of minimization cycles
    rmsd_max_bonds = 0.07
      .type = float
      .short_caption = Max. RMSD(bonds)
      .help = Refinement result is ignored if this max allowable deviation exceeded
    rmsd_max_angles = 5.0
      .type = float
      .short_caption = Max. RMSD(angles)
      .help = Refinement result is ignored if this max allowable deviation exceeded
    grid_search_scales = 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.5 3.0
      .type = floats
      .help = Defines the range of values for the weight grid search
    verbose = 1
      .type = int
      .help = All output is supressed if it is negative
    lockit_parameters
      .expert_level=3
    {
      include scope \
        mmtbx.command_line.lockit.coordinate_refinement_export_phil_str
    }
  }
"""

time_target_simple = 0
time_compute_map   = 0
time_map_cc        = 0
time_update        = 0
time_states        = 0

def show_time(external):
  total = 0
  print "Detailed:"
  print "  time_target_simple: %6.3f" % time_target_simple
  print "  time_compute_map  : %6.3f" % time_compute_map
  print "  time_map_cc       : %6.3f" % time_map_cc
  print "  time_update       : %6.3f" % time_update
  print "  time_states       : %6.3f" % time_states
  total = time_target_simple+\
          time_compute_map+\
          time_map_cc+\
          time_update+\
          time_states
  for e in external:
    print e[0]%e[1]
    total += e[1]
  print "    sub-total       : %6.3f" % total

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

class simple(object):
  def __init__(self,
               target_map,
               selection,
               real_space_gradients_delta,
               geometry_restraints_manager=None,
               max_iterations=150):
    adopt_init_args(self, locals())
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = max_iterations)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.refined = None

  def refine(self, weight, xray_structure):
    assert self.selection.size() == xray_structure.scatterers().size()
    self.refined = maptbx.real_space_refinement_simple.lbfgs(
      selection_variable              = self.selection,
      sites_cart                      = xray_structure.sites_cart(),
      density_map                     = self.target_map,
      geometry_restraints_manager     = self.geometry_restraints_manager,
      real_space_target_weight        = weight,
      real_space_gradients_delta      = self.real_space_gradients_delta,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)

  def sites_cart(self):
    assert self.refined is not None
    sites_cart = self.refined.sites_cart
    sites_cart.set_selected(self.selection, self.refined.sites_cart_variable)
    return sites_cart

class diff_map(object):
  def __init__(self,
               miller_array,
               crystal_gridding,
               map_target,
               geometry_restraints_manager,
               restraints_target_weight = 1,
               max_iterations = 500,
               min_iterations = 500):
    adopt_init_args(self, locals())
    self.step = miller_array.d_min()/4.
    self.refined = None

  def refine(self, weight, sites_cart=None, xray_structure=None):
    assert xray_structure is not None and [sites_cart,xray_structure].count(None)==1
    self.refined = real_space_target_and_gradients.minimization(
      xray_structure              = xray_structure,
      miller_array                = self.miller_array,
      crystal_gridding            = self.crystal_gridding,
      map_target                  = self.map_target,
      max_iterations              = self.max_iterations,
      min_iterations              = self.min_iterations,
      step                        = self.step,
      real_space_target_weight    = weight,
      restraints_target_weight    = self.restraints_target_weight,
      geometry_restraints_manager = self.geometry_restraints_manager,
      target_type                 = "diff_map")

  def sites_cart(self):
    assert self.refined is not None
    return self.refined.xray_structure.sites_cart()

class refinery(object):
  def __init__(self,
               refiner,
               xray_structure,
               start_trial_weight_value = 50.,
               weight_sample_rate = 10,
               rms_bonds_limit = 0.03,
               rms_angles_limit = 3.0,
               optimize_weight = True):
    self.rms_angles_start = None
    self.rms_bonds_start = None
    self.refiner = refiner
    self.weight_start=start_trial_weight_value
    self.rms_bonds_start, self.rms_angles_start  = \
      self.rmsds(sites_cart=xray_structure.sites_cart())
    self.weight_sample_rate = weight_sample_rate
    # results
    self.weight_final = None
    self.sites_cart_result = None
    self.rms_bonds_final,self.rms_angles_final = None,None
    #
    pool = {}
    bonds = flex.double()
    angles = flex.double()
    weights = flex.double()
    #
    weight = start_trial_weight_value
    weight_last = weight
    self.adjust_weight_sample_rate(weight=weight)
    if(optimize_weight):
      while True:
        self.adjust_weight_sample_rate(weight=weight_last)
        refiner.refine(
          xray_structure = xray_structure.deep_copy_scatterers(), # XXX
          weight     = weight)
        sites_cart_result = refiner.sites_cart()
        bd, ad = self.rmsds(sites_cart=sites_cart_result)
        bonds.append(bd)
        angles.append(ad)
        weights.append(weight)
        pool.setdefault(weight,[]).append([sites_cart_result.deep_copy(),bd,ad])
        if(refiner.geometry_restraints_manager is None): break
        weight_last = weight
        if(ad>rms_angles_limit or bd > rms_bonds_limit):
          weight -= self.weight_sample_rate
        else:
          weight += self.weight_sample_rate
        #print ">>> ", "%6.2f %6.2f"%(weight, weight_last), "%6.4f %5.2f"%(bd, ad)
        if((weight<0 or weight>1000) or weight in weights): break
    else:
      refiner.refine(
        xray_structure = xray_structure.deep_copy_scatterers(), # XXX
        weight     = weight)
      sites_cart_result = refiner.sites_cart()
    # select results
    sel  = bonds <= rms_bonds_limit
    sel &= angles <= rms_angles_limit
    bonds   = bonds  .select(sel)
    angles  = angles .select(sel)
    weights = weights.select(sel)
    if(sel.count(True)>0):
      bond_max = flex.max(bonds)
      ind = None
      for i, b in enumerate(bonds):
        if(b==bond_max):
          ind = i
          break
      assert ind is not None
      self.weight_final = weights[ind]
      self.sites_cart_result = pool[self.weight_final][0][0]
      self.rms_bonds_final,self.rms_angles_final = \
        self.rmsds(sites_cart=self.sites_cart_result)
      assert approx_equal(pool[self.weight_final][0][2], angles[ind])
      assert approx_equal(pool[self.weight_final][0][1], bonds[ind])
      assert approx_equal(self.rms_angles_final, angles[ind])
      assert approx_equal(self.rms_bonds_final, bonds[ind])

  def rmsds(self, sites_cart):
    b,a = None,None
    if(self.refiner.geometry_restraints_manager is not None):
      es = self.refiner.geometry_restraints_manager.energies_sites(
        sites_cart = sites_cart)
      a = es.angle_deviations()[2]
      b = es.bond_deviations()[2]
    return b,a

  def adjust_weight_sample_rate(self, weight):
    if(  weight <= 1000.): self.weight_sample_rate=100.
    if(  weight <= 100.):  self.weight_sample_rate=10.
    if(  weight <= 10.):   self.weight_sample_rate=1.
    elif(weight <= 1.0):   self.weight_sample_rate=0.1
    elif(weight <= 0.1):   self.weight_sample_rate=0.01

class states(object):
  def __init__(self, xray_structure, pdb_hierarchy):
    adopt_init_args(self, locals())
    self.counter = 0
    self.root = iotbx.pdb.hierarchy.root()
    self.sites_carts = []

  def add(self, sites_cart):
    global time_states
    timer = user_plus_sys_time()
    self.sites_carts.append(sites_cart)
    ph = self.pdb_hierarchy.deep_copy()
    xrs = self.xray_structure.replace_sites_cart(new_sites = sites_cart)
    ph.adopt_xray_structure(xrs)
    models = ph.models()
    md = models[0].detached_copy()
    md.id = str(self.counter)
    self.root.append_model(md)
    self.counter += 1
    time_states += timer.elapsed()

  def write(self, file_name):
    self.root.write_pdb_file(file_name = file_name)

def target_simple(target_map, sites_cart=None, sites_frac=None, unit_cell=None):
  global time_target_simple
  timer = user_plus_sys_time()
  assert [sites_cart, sites_frac].count(None) == 1
  if(sites_frac is None):
    assert unit_cell is not None
    sites_frac = unit_cell.fractionalize(sites_cart)
  result = 0
  for site_frac in sites_frac:
    result += target_map.eight_point_interpolation(site_frac)
  scale = 1
  if(sites_cart.size() > 0): scale = 1./sites_cart.size()
  time_target_simple += timer.elapsed()
  return result*scale

class evaluator(object):
  def __init__(self, sites_cart, target_map, unit_cell, target):
    adopt_init_args(self, locals())
    if(target is not None):
      self.target = -9999
    else:
      self.target = target_simple(
        target_map = self.target_map,
        sites_cart = self.sites_cart,
        unit_cell  = self.unit_cell)

  def evaluate(self, sites_cart):
    t = target_simple(
      target_map = self.target_map,
      sites_cart = sites_cart,
      unit_cell  = self.unit_cell)
    #print "  ", t
    if(t >= self.target):
      self.sites_cart = sites_cart
      self.target = t

class rsr_residue(object):
  def __init__(self,
               pdb_hierarchy_residue,
               selection_sidechain,
               selection_backbone,
               selection_all,
               map_cc_sidechain=None,
               map_cc_backbone=None,
               map_cc_all=None,
               map_value_sidechain=None,
               map_value_backbone=None,
               distance_to_closest_rotamer=None,
               rotamer_status=None):
    adopt_init_args(self, locals())

def rotamer_fit(residue, target_map, mon_lib_srv, unit_cell, rotamer_manager):
  sites_cart_result = residue.atoms().extract_xyz()
  sc_start = sites_cart_result.deep_copy() #XXX DEBUG
  rotamer_iterator = get_rotamer_iterator(
    mon_lib_srv = mon_lib_srv, residue = residue)
  if(rotamer_iterator is not None):
    rotamer_status = rotamer_manager.evaluate_residue(residue=residue)
    target = None
    if(rotamer_status == "OUTLIER"): target = -9999
    e = evaluator(
      sites_cart = residue.atoms().extract_xyz(),
      target_map = target_map,
      unit_cell  = unit_cell,
      target = target)
    #print residue.resname,residue.resseq, e.target, rotamer_status
    for rotamer, rotamer_sites_cart in rotamer_iterator:
      e.evaluate(sites_cart=rotamer_sites_cart.deep_copy())
      residue.atoms().set_xyz(rotamer_sites_cart) # XXX DEBUG
      rotamer_status = rotamer_manager.evaluate_residue(residue=residue) # XXX DEBUG
      #print dir(residue)
      #STOP()
      #print "      ",rotamer.id
      assert rotamer_status != "OUTLIER"
    sites_cart_result = e.sites_cart
    #print "  final", e.target
  residue.atoms().set_xyz(sc_start) # XXX DEBUG
  return sites_cart_result

def get_rotamer_iterator(mon_lib_srv, residue):
  get_class = iotbx.pdb.common_residue_names_get_class
  rotamer_iterator = None
  if(get_class(residue.resname) == "common_amino_acid"):
    rotamer_iterator = mon_lib_srv.rotamer_iterator(
      fine_sampling = True,
      comp_id       = residue.resname,
      atom_names    = residue.atoms().extract_name(),
      sites_cart    = residue.atoms().extract_xyz())
    if(rotamer_iterator is None or
       rotamer_iterator.problem_message is not None or
       rotamer_iterator.rotamer_info is None):
      rotamer_iterator = None
  return rotamer_iterator

class monitor(object):
  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               target_map,
               geometry_restraints_manager,
               xray_structure_reference = None):
    adopt_init_args(self, locals())
    self.unit_cell = self.xray_structure.unit_cell()
    self.xray_structure_start = xray_structure.deep_copy_scatterers()
    self.states_collector = states(
      pdb_hierarchy  = self.pdb_hierarchy,
      xray_structure = self.xray_structure)
    self.cc = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = None
    self.dist_from_ref = None
    self.number_of_rotamer_outliers = 0
    self.rotamer_manager = RotamerEval()
    self.residues = self.get_residues()
    self.set_globals()
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.set_rsr_residue_attributes()

  def compute_map(self, xray_structure):
    global time_compute_map
    timer = user_plus_sys_time()
    mc = self.target_map.miller_array.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.target_map.crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    time_compute_map += timer.elapsed()
    return fft_map.real_map_unpadded()

  def map_cc(self, map, sites_cart = None):
    global time_map_cc
    timer = user_plus_sys_time()
    result = None
    if(sites_cart is not None):
      sel = maptbx.grid_indices_around_sites(
        unit_cell  = self.unit_cell,
        fft_n_real = map.focus(),
        fft_m_real = map.all(),
        sites_cart = sites_cart,
        site_radii = flex.double(sites_cart.size(), 2))
      result = flex.linear_correlation(
        x=map.select(sel).as_1d(),
        y=self.target_map.data.select(sel).as_1d()).coefficient()
    else:
      result = flex.linear_correlation(
        x=map.as_1d(),
        y=self.target_map.data.as_1d()).coefficient()
    time_map_cc += timer.elapsed()
    return result

  def show_residues(self):
    print
    print "resid    CC(sc) CC(bb) CC(sb)   Rotamer  Dist. to nearest rotamer"
    fmt="%s %s %6.3f %6.3f %6.3f %9s %s"
    for r in self.residues:
      ms=r.map_cc_sidechain
      mb=r.map_cc_backbone
      ma=r.map_cc_all
      if(ms<-1): ms = -1
      if(mb<-1): mb = -1
      if(ma<-1): ma = -1
      print fmt % (
        r.pdb_hierarchy_residue.resname,
        r.pdb_hierarchy_residue.resseq,
        ms,
        mb,
        ma,
        str(r.rotamer_status),
        format_value("%6.3f",r.distance_to_closest_rotamer))

  def set_rsr_residue_attributes(self):
    small = -1.e9
    get_class = iotbx.pdb.common_residue_names_get_class
    unit_cell = self.xray_structure.unit_cell()
    current_map = self.compute_map(xray_structure = self.xray_structure)
    sites_cart = self.xray_structure.sites_cart()
    for r in self.residues:
      sca = sites_cart.select(r.selection_all)
      scs = sites_cart.select(r.selection_sidechain)
      scb = sites_cart.select(r.selection_backbone)
      r.rotamer_status = self.rotamer_manager.evaluate_residue(
        residue=r.pdb_hierarchy_residue)
      if(r.rotamer_status != "OUTLIER"):
        r.map_cc_all       = self.map_cc(sites_cart = sca, map = current_map)
        r.map_cc_sidechain = self.map_cc(sites_cart = scs, map = current_map)
        r.map_cc_backbone  = self.map_cc(sites_cart = scb, map = current_map)
        r.map_value_sidechain = target_simple(target_map=current_map,
          sites_cart=scs, unit_cell=self.unit_cell)
        r.map_value_backbone = target_simple(target_map=current_map,
          sites_cart=scb, unit_cell=self.unit_cell)
      else:
        r.map_cc_all          = small
        r.map_cc_sidechain    = small
        r.map_cc_backbone     = small
        r.map_value_sidechain = small
        r.map_value_backbone  = small
      rotamer_iterator = get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = r.pdb_hierarchy_residue)
      if(rotamer_iterator is not None):
        dist_min = 1.e+9
        for ro, rotamer_sites_cart in rotamer_iterator:
          d=flex.mean(flex.sqrt((sites_cart.select(r.selection_all) - rotamer_sites_cart).dot()))
          if(d < dist_min):
            dist_min = d
        r.distance_to_closest_rotamer = dist_min

  def get_residues(self):
    from mmtbx.command_line import lockit
    residue_groups = self.pdb_hierarchy.residue_groups()
    backbone_atoms = ["N","CA","C","O","CB"]
    selections_all = []
    result = []
    for residue_group in residue_groups:
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          residue_i_seqs_backbone  = flex.size_t()
          residue_i_seqs_sidechain = flex.size_t()
          residue_i_seqs_all = flex.size_t()
          for atom in residue.atoms():
            an = atom.name.strip()
            bb = an in backbone_atoms
            residue_i_seqs_all.append(atom.i_seq)
            if(bb): residue_i_seqs_backbone.append(atom.i_seq)
            else:   residue_i_seqs_sidechain.append(atom.i_seq)
          selections_all.append(residue_i_seqs_all)
          result.append(rsr_residue(
            pdb_hierarchy_residue = residue,
            selection_sidechain   = residue_i_seqs_sidechain,
            selection_backbone    = residue_i_seqs_backbone,
            selection_all         = residue_i_seqs_all))
    sa = flex.size_t()
    for s in selections_all: sa.extend(s)
    assert self.xray_structure.select(sa).scatterers().size() == \
           self.xray_structure.scatterers().size()
    return result

  def set_globals(self):
    self.cc = self.map_cc(map=self.compute_map(
      xray_structure = self.xray_structure))
    es = self.geometry_restraints_manager.energies_sites(
      sites_cart = self.xray_structure.sites_cart())
    self.rmsd_a = es.angle_deviations()[2]
    self.rmsd_b = es.bond_deviations()[2]
    self.dist_from_start = flex.mean(self.xray_structure_start.distances(
      other = self.xray_structure))
    if(self.xray_structure_reference is not None):
      self.dist_from_ref = flex.mean(self.xray_structure_reference.distances(
        other = self.xray_structure))
    self.number_of_rotamer_outliers = 0
    for r in self.residues:
      rotamer_status = self.rotamer_manager.evaluate_residue(
        residue=r.pdb_hierarchy_residue)
      if(rotamer_status == "OUTLIER"):
        self.number_of_rotamer_outliers += 1

  def update(self, xray_structure, accept_any=False):
    global time_update
    timer = user_plus_sys_time()
    unit_cell = xray_structure.unit_cell()
    current_map = self.compute_map(xray_structure = xray_structure)
    sites_cart  = xray_structure.sites_cart()
    sites_cart_ = self.xray_structure.sites_cart()
    for r in self.residues:
      sca = sites_cart.select(r.selection_all)
      scs = sites_cart.select(r.selection_sidechain)
      scb = sites_cart.select(r.selection_backbone)
      map_cc_all       = self.map_cc(sites_cart = sca, map = current_map)
      map_cc_sidechain = self.map_cc(sites_cart = scs, map = current_map)
      map_cc_backbone  = self.map_cc(sites_cart = scb, map = current_map)
      map_value_sidechain = target_simple(target_map=current_map,
        sites_cart=scs, unit_cell=self.unit_cell)
      map_value_backbone = target_simple(target_map=current_map,
        sites_cart=scb, unit_cell=self.unit_cell)
      flag = map_cc_all      > r.map_cc_all and \
             map_cc_backbone > r.map_cc_backbone and \
             map_cc_backbone > map_cc_sidechain
      if(r.map_value_backbone > r.map_value_sidechain):
        if(map_value_backbone < map_value_sidechain):
          flag = False
      if(accept_any): flag=True
      if(flag):
        residue_sites_cart_new = sites_cart.select(r.selection_all)
        sites_cart_.set_selected(r.selection_all, residue_sites_cart_new)
        r.pdb_hierarchy_residue.atoms().set_xyz(residue_sites_cart_new)
        rotamer_status = self.rotamer_manager.evaluate_residue(
          residue=r.pdb_hierarchy_residue)
        r.rotamer_status = rotamer_status
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_)
    self.set_globals()
    self.set_rsr_residue_attributes()
    self.states_collector.add(sites_cart = sites_cart_)
    time_update += timer.elapsed()

  def show(self, suffix=""):
    if(self.dist_from_ref is None):
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f rota: %3d"
      print f%(self.cc,self.rmsd_b,self.rmsd_a,self.dist_from_start,
               self.number_of_rotamer_outliers),suffix
    else:
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f d(ref): %6.3f"
      print f%(self.map_cc,self.rmsd_b,self.rmsd_a,self.dist_from_start,
        self.dist_from_ref),suffix

def run(target_map,
        pdb_hierarchy,
        xray_structure,
        geometry_restraints_manager,
        xray_structure_reference = None,
        rms_bonds_limit  = 0.03,
        rms_angles_limit = 3.0,
        max_iterations   = 500,
        macro_cycles     = 20,
        minimization     = True,
        expload          = True,
        rotamer_search   = True,
        verbose          = True):
  sel = flex.bool(xray_structure.scatterers().size(), True)
  d_min = target_map.miller_array.d_min()
  #
  #geometry_restraints_manager.geometry.remove_dihedrals_in_place(sel)
  #
  rsr_simple_refiner = simple(
    target_map                  = target_map.data,
    selection                   = sel,
    real_space_gradients_delta  = d_min/4,
    max_iterations              = max_iterations,
    geometry_restraints_manager = geometry_restraints_manager.geometry)

  rsr_diff_map_refiner = diff_map(
    miller_array                = target_map.miller_array,
    crystal_gridding            = target_map.crystal_gridding,
    map_target                  = target_map.data,
    geometry_restraints_manager = geometry_restraints_manager.geometry,
    restraints_target_weight    = 1,
    max_iterations              = max_iterations,
    min_iterations              = max_iterations)
  #
  if(xray_structure_reference is not None):
    xray_structure_reference = xray_structure_reference.deep_copy_scatterers()
  monitor_object = monitor(
    pdb_hierarchy = pdb_hierarchy,
    xray_structure = xray_structure.deep_copy_scatterers(),
    target_map = target_map,
    geometry_restraints_manager = geometry_restraints_manager.geometry,
    xray_structure_reference = xray_structure_reference)
  if(verbose):
    monitor_object.show(suffix="start")
    monitor_object.show_residues()
  #
  tmp = monitor_object.xray_structure.deep_copy_scatterers()
  weight_d = 50
  weight_s = 50
  #
  for i in range(macro_cycles):
    if(expload and i>1 and i%2==0):
      tmp_dc = tmp.deep_copy_scatterers()
      tmp.shake_sites_in_place(mean_distance=3) # reverse back if refinement failed
    if(minimization):
      target_type = "simple"
      refined = refinery(
        refiner          = rsr_simple_refiner,
        xray_structure   = tmp,
        start_trial_weight_value = weight_s,
        rms_bonds_limit  = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      if(refined.sites_cart_result is not None):
        tmp = tmp.replace_sites_cart(refined.sites_cart_result)
        weight_s = refined.weight_final
        monitor_object.update(xray_structure=tmp)#, accept_any=True)
        if(verbose):
          monitor_object.show(suffix=" weight: %s"%str(weight_s))
        tmp = monitor_object.xray_structure.deep_copy_scatterers()
      else:
        tmp = tmp_dc
        print "Refinement failed."
        #
    if(rotamer_search or ((not expload or (expload and minimization)) and i>macro_cycles/2)):
      sites_cart = tmp.sites_cart()
      for r in monitor_object.residues:
        sites_cart_ = rotamer_fit(
          residue     = r.pdb_hierarchy_residue,
          target_map  = target_map.data,
          mon_lib_srv = monitor_object.mon_lib_srv,
          unit_cell   = tmp.unit_cell(),
          rotamer_manager = monitor_object.rotamer_manager)
        sites_cart.set_selected(r.selection_all, sites_cart_)
      tmp = tmp.replace_sites_cart(sites_cart)
      monitor_object.update(xray_structure=tmp, accept_any=True)
      if(verbose):
        monitor_object.show(suffix=" weight: %s"%str(None))
      tmp = monitor_object.xray_structure.deep_copy_scatterers()
        #

  #
  if(minimization):
    refined = refinery(
      refiner          = rsr_simple_refiner,
      xray_structure   = tmp,
      start_trial_weight_value = weight_s,
      rms_bonds_limit  = 0.02,
      rms_angles_limit = 2.5)
    if(verbose):
      print "FINAL:", refined.rms_bonds_final,refined.rms_angles_final
    if(refined.sites_cart_result is not None):
      tmp = tmp.replace_sites_cart(refined.sites_cart_result)
      weight_s = refined.weight_final
      monitor_object.update(xray_structure=tmp, accept_any=True) # XXX ???
      if(verbose): monitor_object.show(suffix=" weight: %s"%str(weight_s))
  #
  if(verbose): monitor_object.show_residues()
  monitor_object.states_collector.write(file_name = "all.pdb")
  return monitor_object.xray_structure


class box_refinement_manager(object):
  def __init__(self,
               xray_structure,
               pdb_hierarchy,
               target_map,
               geometry_restraints_manager,
               real_space_gradients_delta=1./4,
               max_iterations = 50):
    self.xray_structure = xray_structure
    self.sites_cart = xray_structure.sites_cart()
    self.pdb_hierarchy = pdb_hierarchy
    self.target_map = target_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.max_iterations=max_iterations
    self.real_space_gradients_delta = real_space_gradients_delta

  def update_xray_structure(self, new_xray_structure):
    self.xray_structure = new_xray_structure

  def update_target_map(self, new_target_map):
    self.target_map = new_target_map

  def refine(self,
             selection,
             selection_buffer_radius=5,
             box_cushion=2,
             monitor_clashscore=False):
    sites_cart_moving = self.sites_cart
    selection_within = self.xray_structure.selection_within(
      radius    = selection_buffer_radius,
      selection = selection)
    sel = selection.select(selection_within)
    iselection = flex.size_t()
    for i, state in enumerate(selection):
      if state:
        iselection.append(i)
    if monitor_clashscore:
      pdb_string = utils.write_pdb_file(
                     xray_structure=self.xray_structure,
                     pdb_hierarchy=self.pdb_hierarchy,
                     write_cryst1_record = False,
                     selection = selection_within,
                     return_pdb_string = True)
      csm = clashscore.probe_clashscore_manager(
              pdb_string=pdb_string)
      self.clashscore = csm.clashscore
    box = utils.extract_box_around_model_and_map(
            xray_structure   = self.xray_structure,
            pdb_hierarchy    = self.pdb_hierarchy,
            map_data         = self.target_map,
            selection        = selection_within,
            box_cushion      = box_cushion)
    new_unit_cell = box.xray_structure_box.unit_cell()
    geo_box = \
      self.geometry_restraints_manager.select(box.selection_within)
    geo_box.discard_symmetry(new_unit_cell=new_unit_cell)
    map_box = box.map_box
    sites_cart_box = box.xray_structure_box.sites_cart()
    selection = flex.bool(sites_cart_box.size(), True)
    rsr_simple_refiner = simple(
      target_map                  = map_box,
      selection                   = sel,
      real_space_gradients_delta  = self.real_space_gradients_delta,
      max_iterations              = self.max_iterations,
      geometry_restraints_manager = geo_box)
    real_space_result = refinery(
      refiner          = rsr_simple_refiner,
      xray_structure   = box.xray_structure_box)
    sites_cart_box_refined = real_space_result.sites_cart_result
    sites_cart_box_refined_shifted_back = \
      sites_cart_box_refined + box.shift_to_map_boxed_sites_back
    sites_cart_refined = sites_cart_box_refined_shifted_back.select(
                           sel)
    sites_cart_moving = sites_cart_moving.set_selected(
      iselection, sites_cart_refined)
    self.xray_structure.set_sites_cart(sites_cart_moving)
    self.sites_cart = self.xray_structure.sites_cart()
    if monitor_clashscore:
      pdb_string = utils.write_pdb_file(
                     xray_structure=self.xray_structure,
                     pdb_hierarchy=self.pdb_hierarchy,
                     write_cryst1_record = False,
                     selection = selection_within,
                     return_pdb_string = True)
      csm = clashscore.probe_clashscore_manager(
              pdb_string=pdb_string)
      self.clashscore_refined = csm.clashscore
