from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.utils
from scitbx.array_family import flex
import mmtbx.refinement.real_space.explode_and_refine
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import mmtbx.refinement.real_space.individual_sites
from libtbx import group_args
import iotbx.map_manager
import sys
import iotbx.map_model_manager
import inspect
from libtbx.utils import user_plus_sys_time

class run(object):
  """
  Build and real-space-refine ligand into map. Ligand must be a linear molecule
  such as ATP or PEG.

  map_model_manager:
    - contains masked ligand map in the box;
    - origin is zero;
    - unit cell and symmetry correspond to box cell and symmetry;
    - map values outside ligand region are set to zero;
    - map values inside ligand region are actual map values;
    - contains model with the ideal liagnd from the library positioned
      arbitrarily in space;
  """

  def __init__(self, map_model_manager, d_min, log = sys.stdout):
    # Timing
    self.total_time = 0
    self.time_strings = []
    # Aliases
    self.log = log
    self.d_min = d_min
    self.mmm = map_model_manager
    self.cs = self.mmm.crystal_symmetry()
    self.uc = self.cs.unit_cell()
    # Initiate states collector
    self.states = self.caller(func = self._init_states_accumulator)
    # Modify maps
    self.map_data_ref, self.mc_ref, self.map_data = self.caller(
      func=self._get_modified_maps)
    # Spread atoms evenly inside the blob
    self.dummy_atoms_in_map = self.caller(func = self._initial_trace)
    # pdb hierarchies corresponding to two initial placements
    self.hierarchies_placed = []
    # Place ligand atoms along the trace (forward)
    self.hierarchies_placed.append( self.caller(func = self._place_forward) )
    # Place ligand atoms along the trace (backward)
    self.hierarchies_placed.append( self.caller(func = self._place_backward) )
    # Choose one best placement, set result to mmm's model
    self.caller(func = self._choose_one_placement)
    # Model fragment analysis
    self.fragments = self.caller(func = self._fragments)
    # Explode and refine
    self.caller(func = self._ear)
    # Write final model and all states
    self.caller(func = self._write_final_model)
    self._show_timing()

  def caller(self, func):
    timer = user_plus_sys_time()
    doc = inspect.getdoc(func)
    #if(doc is not None): broadcast(m = doc, log = self.log)
    result = func()
    t = timer.elapsed()
    self.total_time += t
    fmt = "  %s: %s"%(doc, str("%8.3f"%t).strip())
    self.time_strings.append(fmt)
    #self.log.flush()
    return result

  def _print(self, m):
    if(self.log is None): return
    print(m, file=self.log)

  def _show_cc(self, sites_cart):
    xrs = self.mmm.model().get_xray_structure()
    xrs.set_sites_cart(sites_cart)
    fc = self.mc_ref.structure_factors_from_scatterers(
      xray_structure = xrs).f_calc()
    return fc.map_correlation(other=self.mc_ref)

  def _write_final_model(self):
    """
    Write final model and all states
    """
    self.states.write(file_name = "all.pdb")
    self.mmm.write_model(file_name = "final.pdb")

  def _ear(self):
    """
    Explode-and-refine
    """
    model = self.mmm.model()
    xrs = model.get_xray_structure()
    ph = model.get_hierarchy()
    states = mmtbx.utils.states(pdb_hierarchy=ph)
    states.add(sites_cart = xrs.sites_cart())
    ear = mmtbx.refinement.real_space.explode_and_refine.run(
      xray_structure     = model.get_xray_structure(),
      pdb_hierarchy      = model.get_hierarchy(),
      map_data           = self.map_data,
      restraints_manager = model.get_restraints_manager(),
      states             = states,
      resolution         = self.d_min,
      map_data_ref       = self.map_data_ref,
      mode               = "thorough",
      score_method       = ["cc","geometry"],
      fragments          = self.fragments,
      number_of_trials   = 25,
      nproc              = 1,
      log=null_out())
    #
    cc = self._show_cc(sites_cart = ear.pdb_hierarchy.atoms().extract_xyz())
    self._print("Final, CC: %6.4f"%cc)
    self.states.add(hierarchy = ear.pdb_hierarchy)

  def _fragments(self):
    """
    Split into fragments
    """
    return get_fragments(model = self.mmm.model())

  def _t_search(self, model):
    """
    Translational grid search
    """
    CCBEST=-1
    sites_cart = model.get_sites_cart()
    SCBEST=None
    for x in [-2,0,2]:
      for y in [-2,0,2]:
        for z in [-2,0,2]:
          shift = flex.vec3_double([[x,y,z],]*model.size())
          model.set_sites_cart(sites_cart = sites_cart + shift)
          model = refine(model=model, map_data=self.map_data)
          CC = get_cc(
            xrs      = model.get_xray_structure(),
            d_min    = self.d_min,
            map_data = self.map_data_ref)
          if(CC>CCBEST):
            CCBEST=CC
            SCBEST = model.get_sites_cart().deep_copy()
    model.set_sites_cart(sites_cart=SCBEST)
    self.states.add(hierarchy = model.get_hierarchy())

  def _choose_one_placement(self):
    """
    Choose one placement: forward vs backward
    """
    cc_best = -1
    sites_cart_best = None
    for ih, ph in enumerate(self.hierarchies_placed):
      self._print("Orientation %d:"%ih)
      model = get_model_adhoc(crystal_symmetry = self.cs, ph = ph)
      # Pre-refine original placement
      for it in [1,2]:
        # XXX This invalidates grm with no way to get it back!
        #if(it==1): model.set_nonbonded_weight(value=1)
        #else:      model.set_nonbonded_weight(value=1000)
        model = sa_simple(model=model, map_data=self.map_data, log=None)
        model = refine(model=model, map_data=self.map_data)
      self.states.add(hierarchy = model.get_hierarchy())
      cc = self._show_cc(sites_cart = model.get_sites_cart())
      self._print(" SA and minimization, CC: %6.4f"%cc)
      #
      self._t_search(model = model)
      # Evaluate
      cc = self._show_cc(sites_cart = model.get_sites_cart())
      self._print(" Translaton search, CC: %6.4f"%cc)
      if(cc>cc_best):
        cc_best = cc
        sites_cart_best = model.get_sites_cart()
    model.set_sites_cart(sites_cart = sites_cart_best)
    self.states.add(hierarchy = model.get_hierarchy())
    # Atom order is different, so relly need to do this
    self.mmm.set_model(model)

  def _show_timing(self):
    print("Detailed timing:")
    max_len = int(flex.max(flex.double(
      [len(" ".join(it.split()[:-1])) for it in self.time_strings])))
    fmt = "%-"+str(max_len)+"s"
    cntr = 0
    for it in self.time_strings:
      p = it.split()
      p1, p2 = " ".join(p[:-1]), p[-1]
      print("  "+fmt%p1+p2)
      cntr += float(p2)
    print("Total: %8.3f"%self.total_time)
    assert approx_equal(cntr, self.total_time) # total time and the sum match!

  def _init_states_accumulator(self):
    """
    Initialize states accumulator
    """
    return mmtbx.utils.states(pdb_hierarchy = self.mmm.model().get_hierarchy())

  def _get_modified_maps(self):
    """
    Make masked and negated maps
    """
    map_data_ref = self.mmm.map_data().deep_copy()
    map_data_ref = map_data_ref.set_selected(map_data_ref < 1, 0)
    map_data = self.mmm.map_data().deep_copy()
    map_data = map_data.set_selected(map_data<3, -10)
    fc = self.mmm.model().get_xray_structure().structure_factors(
      d_min=self.d_min).f_calc()
    mc_ref = fc.structure_factors_from_map(map=map_data_ref, use_sg=False)
    return map_data_ref, mc_ref, map_data

  def _place_forward(self):
    """
    Replace dummy atoms with ligand atoms (forward)
    """
    hierarchy = trace_to_hierarchy(
      ph         = self.mmm.model().get_hierarchy(),
      trace_cart = self.dummy_atoms_in_map,
      uc         = self.uc)
    self.states.add(hierarchy = hierarchy)
    cc = self._show_cc(sites_cart = hierarchy.atoms().extract_xyz())
    self._print("Trace forward, CC: %6.4f"%cc)
    return hierarchy

  def _place_backward(self):
    """
    Replace dummy atoms with ligand atoms (backward)
    """
    hierarchy = trace_to_hierarchy(
      ph         = self.mmm.model().get_hierarchy(),
      trace_cart = self.dummy_atoms_in_map,
      uc         = self.uc,
      reverse    = True)
    self.states.add(hierarchy = hierarchy)
    cc = self._show_cc(sites_cart = hierarchy.atoms().extract_xyz())
    self._print("Trace backward, CC: %6.4f"%cc)
    return hierarchy

  def _initial_trace(self):
    """
    Place dummy atoms into map
    """
    sites_cart = self.mmm.map_manager().trace_atoms_in_map(
      dist_min = 2, n_atoms  = self.mmm.model().size())
    #
    fmt = "HETATM %4d   O  HOH %5d    %8.3f%8.3f%8.3f  1.00 30.00           O"
    lines = "\n".join(
      [fmt%(i,i,sc[0],sc[1],sc[2]) for i, sc in enumerate(sites_cart)])
    pdb_inp = iotbx.pdb.input(source_info=None, lines = lines)
    model = get_model_adhoc(crystal_symmetry = self.cs, pdb_inp = pdb_inp)

    self.states.add(hierarchy = model.get_hierarchy())
    model = sa_simple(model=model, map_data=self.map_data, log=null_out())
    #
    sites_cart = model.get_sites_cart()
    self.states.add(hierarchy = model.get_hierarchy())
    #
    start, end = get_se(coords = sites_cart, uc = self.uc)
    distances = flex.double([dist(start, s, self.uc) for s in sites_cart])
    sel = flex.sort_permutation(distances)
    sites_cart = sites_cart.select(sel)
    #
    return list(sites_cart)














































def dist(p1,p2, uc):
  return uc.distance(uc.fractionalize(p1), uc.fractionalize(p2))

def sa_simple(
        model,
        map_data,
        log):
  tmp_xrs = model.get_xray_structure().deep_copy_scatterers()
#  ro = mmtbx.refinement.real_space.individual_sites.easy(
#    map_data                    = map_data,
#    xray_structure              = tmp_xrs,
#    pdb_hierarchy               = model.get_hierarchy().deep_copy(),
#    geometry_restraints_manager = model.get_restraints_manager(),
#    rms_bonds_limit             = 0.01,
#    rms_angles_limit            = 1.0,
#    selection                   = None, #TODO
#    log                         = log)
#  weight = ro.w
  weight=50
  #
  from mmtbx.dynamics import simulated_annealing as sa
  tmp = model.get_xray_structure().deep_copy_scatterers()
  params = sa.master_params().extract()
  params.start_temperature=5000
  params.cool_rate=500
  sa.run(
    params             = params,
    xray_structure     = tmp,
    real_space         = True,
    target_map         = map_data,
    restraints_manager = model.get_restraints_manager(),
    wx                 = weight,
    wc                 = 1.,
    verbose            = False,
    log                = log)
  model.set_sites_cart(sites_cart=tmp.sites_cart())
  return model

def refine(model, map_data, start=None, end=None):
  #rm = model.get_restraints_manager().geometry
  #print dir(rm)
  #from mmtbx.geometry_restraints import reference
  #sites_cart_reference = flex.vec3_double([end,start])
  #selection = flex.size_t([0,30])
  #rm.adopt_reference_coordinate_restraints_in_place(
  #  reference.add_coordinate_restraints(
  #      sites_cart = sites_cart_reference,
  #      selection = selection,
  #      sigma = 0.02))

  ro = mmtbx.refinement.real_space.individual_sites.easy(
    map_data                    = map_data,
    xray_structure              = model.get_xray_structure(),
    pdb_hierarchy               = model.get_hierarchy(),
    geometry_restraints_manager = model.get_restraints_manager(),
    rms_bonds_limit             = 0.01,
    rms_angles_limit            = 1,
    selection                   = None, #TODO
    selection_real_space        = flex.bool(model.size(),True),
    w                           = 10,#None,
    log                         = None)
  model.set_sites_cart(sites_cart = ro.xray_structure.sites_cart())
  return model

def merge_common(lists):
    from collections import defaultdict
    neigh = defaultdict(set)
    visited = set()
    for each in lists:
        for item in each:
            neigh[item].update(each)
    def comp(node, neigh = neigh, visited = visited, vis = visited.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            vis(node)
            nodes |= neigh[node] - visited
            yield node
    for node in neigh:
        if node not in visited:
            yield sorted(comp(node))

def get_se(coords, uc):
  d = -1
  start, end = None,None
  for i, c1 in enumerate(coords):
    for j, c2 in enumerate(coords):
      d_ = uc.distance(c1,c2)
      if(d_>d):
        start = c1
        end   = c2
        d     = d_
  return start, end

def get_fragments(model):
  rm = model.get_restraints_manager()
  atoms = model.get_hierarchy().atoms()
  all_selection = list(range(atoms.size()))
  # Planes
  planes = []
  planes_all = []
  for p in rm.geometry.planarity_proxies:
    planes.append(list(p.i_seqs))
    planes_all.extend(list(p.i_seqs))
  planes_all = list(set(planes_all))
#  print "planes    :", planes
#  print "planes_all:", planes_all
#  print
  # Chiral
  chirals = []
  chirals_unique = []
  for p in rm.geometry.chirality_proxies:
#    print "chiral:", p.i_seqs, [atoms[i].name for i in p.i_seqs]
    chirals.append(list(p.i_seqs))
    tmp=[]
    for i in p.i_seqs:
      if(not i in planes_all):
        tmp.append(i)
    chirals_unique.append(tmp)
  chirals_unique = list(merge_common(chirals_unique))[0]
  chirals_all    = list(merge_common(chirals))[0]
#  print "chirals_unique:", chirals_unique
#  print "chirals_all   :", chirals_all
  chirals_mapping = [chirals_all.index(u) for u in chirals_unique]
#  print "mapping:", chirals_mapping
  # Dihedral
  dihedrals = []
  dihedrals_unique = []
  for p in rm.geometry.dihedral_proxies:
#    print "dihedral:", p.i_seqs, [atoms[i].name for i in p.i_seqs]
    dihedrals.append(list(p.i_seqs))
    tmp=[]
    for i in p.i_seqs:
      if(not i in planes_all+chirals_all):
        tmp.append(i)
    dihedrals_unique.append(tmp)
  dihedrals_unique = list(merge_common(dihedrals_unique))[0]
  dihedrals_all    = list(merge_common(dihedrals))[0]
#  print "dihedrals_unique:", dihedrals_unique
#  print "dihedrals_all   :", dihedrals_all
  # Finalize
  pcd = dihedrals_all + chirals_all + planes_all
  #
  tmp = []
  for s in all_selection:
    if s in pcd: continue
    tmp.append(s)
  left = tmp[:]
#  print "pcd :", pcd
#  print "left:", left
  #
  angles = []
  for p in rm.geometry.angle_proxies:
    present = False
    for i in p.i_seqs:
      if i in left:
        present=True
        break
    if not present: continue
#    print list(p.i_seqs)
    angles.append(list(p.i_seqs))
  angles = list(merge_common(angles))
  #
  dihedrals_unique = [dihedrals_unique] + angles
  dihedrals_all    = [dihedrals_all   ] + angles
  dihedrals_unique = list(merge_common(dihedrals_unique))[0]
  dihedrals_all    = list(merge_common(dihedrals_all))[0]
#  print "dihedrals_unique:", dihedrals_unique
#  print "dihedrals_all   :", dihedrals_all
  dihedrals_mapping = [dihedrals_all.index(u) for u in dihedrals_unique]
#  print "mapping:", dihedrals_mapping

  # check
  pcd = dihedrals_unique + chirals_unique + planes_all
  pcd.sort()
  assert approx_equal(pcd, all_selection)
  #
  return group_args(
    planes_all       = planes_all,

    chirals_unique    = chirals_unique,
    chirals_all       = chirals_all,
    chirals_mapping   = flex.size_t(chirals_mapping),

    dihedrals_unique  = dihedrals_unique,
    dihedrals_all     = dihedrals_all,
    dihedrals_mapping = flex.size_t(dihedrals_mapping))

def get_cc(xrs, d_min, map_data):
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  fo = fc.structure_factors_from_map(map=map_data, use_sg=False)
  return fc.map_correlation(other=fo)

def trace_to_hierarchy(ph, trace_cart, uc, reverse=False):
  if(reverse): trace_cart.reverse()
  ph_dc = ph.deep_copy()
  xyz = flex.vec3_double(trace_cart)
  ss =  flex.std_string(ph_dc.as_pdb_string().splitlines())
  sites_cart = ph_dc.atoms().extract_xyz()
  start, end = get_se(coords=sites_cart, uc = uc)
  distances = flex.double([dist(start, s, uc) for s in sites_cart])
  sel = flex.sort_permutation(distances)
  ss = ss.select(sel)
  pi = iotbx.pdb.input(source_info=None, lines=ss)
  ph_dc = pi.construct_hierarchy(sort_atoms=False)
  ph_dc.atoms().set_xyz(xyz)
  return ph_dc

def get_model_adhoc(crystal_symmetry, ph=None, pdb_inp=None):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  params.pdb_interpretation.nonbonded_weight=1000
  if(pdb_inp is not None):
    model = mmtbx.model.manager(
      model_input      = pdb_inp,
      crystal_symmetry = crystal_symmetry,
      log              = null_out())
  else:
    model = mmtbx.model.manager(
       model_input      = None,
       crystal_symmetry = crystal_symmetry,
       log              = null_out(),
       pdb_hierarchy    = ph)
  model.process(pdb_interpretation_params=params, make_restraints=True)
  return model
