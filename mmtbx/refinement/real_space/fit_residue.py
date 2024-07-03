from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.refinement.real_space
from mmtbx.refinement.real_space import individual_sites
import math, sys
from cctbx import maptbx
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager
import collections
from libtbx import group_args
from collections import OrderedDict
from scitbx.matrix import rotate_point_around_axis

import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("mmtbx_rotamer_fit_ext")

def flatten(l):
  if l is None: return None
  return sum(([x] if not (isinstance(x, list) or isinstance(x, flex.size_t))
    else flatten(x) for x in l), [])

class monitor(object):
  def __init__(self, id_str, selection, map_data, unit_cell, weights, pairs,
               cmv, rotamer_evaluator, log):
    adopt_init_args(self, locals())
    self.states = collections.OrderedDict()

  def add(self, residue, state):
    vals = collections.OrderedDict()
    target     = 0
    target_neg = 0
    exceed_map_max_value = False
    for i in self.selection:
      atom = residue.atoms()[i]
      key = "%s_%s_%s"%(
        atom.parent().parent().parent().id, atom.parent().resname,
        atom.name.strip())
      name = atom.name.strip().upper()
      element = atom.element.strip().upper()
      if(element in ["H","D"]): continue
      mv = self.map_data.eight_point_interpolation(
        self.unit_cell.fractionalize(atom.xyz))
      vals[name] = mv
      if(mv > self.cmv[key]*3 and not element in ["S","SE"]):
        exceed_map_max_value = True
      target += mv
      if(mv < 0): target_neg += mv
    #
    rot = self.rotamer_evaluator.evaluate_residue(residue)
    self.states[state] = group_args(
      vals       = vals,
      sites_cart = residue.atoms().extract_xyz(),
      target     = target,
      target_neg = target_neg,
      exceed_map_max_value = exceed_map_max_value,
      rot        = rot)

  def finalize(self, residue):
    if(len(self.states.keys())==1 or self.selection.size()==0): return
    state = "start"
    S     = self.states["start"]
    F     = self.states["fitting"]
    try:
      T     = self.states["tuneup"]
    except KeyError:
      T = None
    #
    if((S.rot=="OUTLIER" and F.rot!="OUTLIER" and not F.exceed_map_max_value) or
       (F.target>S.target and F.target_neg>=S.target_neg and not F.exceed_map_max_value) or
       (F.target_neg>S.target_neg) or
       (S.exceed_map_max_value and not F.exceed_map_max_value)):
      state = "fitting"
    N = self.states[state]
    if(T is not None):
      if((N.rot=="OUTLIER" and T.rot!="OUTLIER" and not T.exceed_map_max_value) or
         (T.target>N.target and T.target_neg>=S.target_neg and not T.exceed_map_max_value) or
         (T.target_neg>N.target_neg) or
         (N.exceed_map_max_value and not T.exceed_map_max_value)):
        state = "tuneup"
    #
    residue.atoms().set_xyz(self.states[state].sites_cart)
    #
    if(T is not None and state != "tuneup"):
      self.add(residue = residue, state = "revert")
    #

  def show(self):
    if(len(self.states.keys())==1 or self.selection.size()==0): return
    print(self.id_str, file=self.log)
    for k,v in zip(self.states.keys(), self.states.values()):
      vals = " ".join(["%s: %5.2f"%(k_, v_)
        for k_,v_ in zip(v.vals.keys(), v.vals.values())])
      print("  %7s: score: %7.3f %s %s"%(k, v.target, vals, v.rot), file=self.log)

def find_peak(angles, values, threshold):
  peaks = []
  s = len(values)
  #print("looking at peaks", s)
  for i, a in enumerate(angles):
    v = values[i]
    if v < threshold: continue
    ip1 = i+1 if i+1 < s else abs(s-(i+1))
    ip2 = i+2 if i+2 < s else abs(s-(i+2))
    ip3 = i+3 if i+3 < s else abs(s-(i+3))
    vip1 = values[ip1]
    vip2 = values[ip2]
    vip3 = values[ip3]
    im1 = i-1
    im2 = i-2
    im3 = i-3
    vim1 = values[im1]
    vim2 = values[im2]
    vim3 = values[im3]
    #print("   ", a, v, "|", vim1,vim2,vim3,  vip1,vip2,vip3)
    stop = False
    for vv in [vip1,vip2,vip3, vim1,vim2,vim3]:
      if vv < threshold: stop=True
    if stop: continue
    fl=vim1<v and vim2<vim1 and vim3<vim2 and  vip1<v and vip2<vip1 and vip3<vip2
    #print(fl, "<>", vim1<v , vim2<vim1 , vim3<vim2 ,  vip1<v , vip2<vip1 , vip3<vip2)
    if(fl):
      peaks.append(a)
  return peaks

class find_all_conformers(object):

  def __init__(self, residue, map_data, mon_lib_srv, unit_cell, threshold,
                     rotamer_evaluator):
    adopt_init_args(self, locals())
    sites_cart = residue.atoms().extract_xyz()
    self.side_chain_selection = flex.size_t()
    for a in residue.atoms():
      if a.name.strip().upper() in ["CA","O","C","N","CB"]: continue
      if a.element_is_hydrogen(): continue
      self.side_chain_selection.append(a.i_seq)
    #
    self.target_start, _ = self._get_target(sites_cart)
    #
    co = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = residue,
      mon_lib_srv     = mon_lib_srv,
      backbone_sample = False)
    clusters = co.clusters
    self.rsr_eval_selection = co.rsr_eval_selection

    #print()
    #print(residue.resname, residue.resseq)
    #print("hd_selection", list(self.hd_selection))
    #print("rsr_eval_selection", list(co.rsr_eval_selection))
    #for c in clusters:
    #  print(list(c.axis))
    #  print(list(c.atoms_to_rotate))
    #
    self.sites_cart_conformers_all = []
    last=False
    for i, cl in enumerate(clusters):
      if i == len(clusters)-1: last=True
      if i==0:
        self.func(cl=cl, sites_cart=sites_cart,  last=last)
      else:
        for sites_cart_ in self.sites_cart_conformers_all[:]:
          self.func(cl=cl, sites_cart=sites_cart_, last=last)
    #
    self.rotamers = {}
    tbest = None
    tmp = []
    for sites_cart in self.sites_cart_conformers_all:
      residue.atoms().set_xyz(sites_cart)
      rotamer_status = rotamer_evaluator.evaluate_residue(residue)
      if rotamer_status == "OUTLIER": continue
      t, vals = self._get_target(sites_cart=sites_cart)
      #print(rotamer_status, " ".join(["%6.2f"%v for v in vals]))
      #
      sel = vals < threshold
      if sel.count(True)>0: continue
      #if t<=self.target_start: continue
      #print(rotamer_status, " ".join(["%6.2f"%v for v in vals]), "<<<<<<<<", t)
      self.rotamers.setdefault(rotamer_status, []).append([t,sites_cart])
      tmp.append(sites_cart)
    self.sites_cart_conformers_all = tmp
    #
    #print("rotamers")
    #for k,v in self.rotamers.items():
    #  print("      ", k,[vv[0] for vv in v])

    self.unique = {}
    for k,v in self.rotamers.items():
      vbest  = -1.e9
      ssbest = None
      for it in v:
        vv, ss = it
        if vv > vbest:
          vbest = vv
          ssbest = ss
      self.unique[k] = [vbest, ssbest]
    #print("unique")
    #for k, v in self.unique.items():
    #  print("   ", k, v)
    #
    self.unique_sorted = None

  def _get_target(self, sites_cart):
    vals = flex.double()
    for sc in sites_cart.select(self.side_chain_selection):
      vals.append(self.map_data.tricubic_interpolation(
        self.unit_cell.fractionalize(sc)))
    return flex.mean(vals), vals

  def sorted_by_map_value(self, n=None):
    tbest=-9999
    scbest = None
    if self.unique=={}: return None
    for k, v in self.unique.items():
      t = v[0]
      if t > tbest:
        tbest = t
        scbest = v[1].deep_copy()
    return [scbest]

    #if self.unique_sorted is not None: return self.unique_sorted
    #self.unique_sorted = OrderedDict()
    #for k,v in self.unique.items():
    #  self.unique_sorted[v[0]] = v[1]
    #self.unique_sorted = OrderedDict(reversed(list(self.unique_sorted.items())))
    ##for k,v in self.unique_sorted.items():
    ##  print("LOOK", k,v)
    #if n is None: return self.unique_sorted
    #else:
    #  rv = list(self.unique_sorted.values())[:n]
    #  return rv

  def func(self, cl, sites_cart, last):
    if not last: points = [sites_cart[cl.atoms_to_rotate[0]]]
    else:
      points = []
      for i in cl.atoms_to_rotate:
        if not i in self.rsr_eval_selection: continue
        points.append(sites_cart[i])
    angles = flex.double()
    values = flex.double()

    #print("using axis", list(cl.axis))

    if last and len(points)>1: end=180
    else:                      end=360
    for angle in range(0, end, 1):
      vals = flex.double()
      for point in points:
        new_point = rotate_point_around_axis(
          axis_point_1 = sites_cart[cl.axis[0]],
          axis_point_2 = sites_cart[cl.axis[1]],
          point        = point,
          angle        = angle, deg=True)
        value = self.map_data.eight_point_interpolation(
          self.unit_cell.fractionalize(new_point))
        vals.append(value)
      #print(angle, flex.mean(vals))
      angles.append(angle)
      values.append(flex.mean(vals))
    peaks = find_peak(angles=angles, values=values, threshold=self.threshold)
    #print(peaks)
    for peak in peaks:
      sites_cart_ = sites_cart.deep_copy()
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = sites_cart_[cl.axis[0]],
          axis_point_2 = sites_cart_[cl.axis[1]],
          point        = sites_cart_[atom],
          angle        = peak, deg=True)
        sites_cart_[atom] = new_xyz
      self.sites_cart_conformers_all.append(sites_cart_)

class run(object):
  def __init__(self,
               residue,
               mon_lib_srv,
               rotamer_manager,
               sin_cos_table,
               cmv,
               unit_cell,
               rotatable_hd=None,
               vdw_radii=None,
               xyzrad_bumpers=None,
               target_map=None,
               target_map_for_cb=None,
               backbone_sample=False,
               accept_only_if_max_shift_is_smaller_than=None,
               trust_map_values_real=True,
               log=None):
    adopt_init_args(self, locals())
    if(self.log is None): self.log = sys.stdout
    self.co = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = True,
      log             = self.log)
    self.m = None
    if(self.target_map is not None and len(self.co.clusters)>0):
      # Set weights
      AN = {"S":16, "O":8, "N":7, "C":6, "SE":34, "H":1, "D":5}
      #AN = {"S":1, "O":1, "N":1, "C":1, "SE":1, "H":1}
      self.weights = flex.double()
      for atom in self.residue.atoms():
        self.weights.append(AN[atom.element.strip().upper()])
      # Bonded pairs
      exclude = ["C","N","O","CA"]
      reference = exclude + ["CB"]
      atoms = self.residue.atoms()
      self.pairs = []
      for i, ai in enumerate(atoms):
        if(ai.name.strip() in reference):
          mv = self.target_map.eight_point_interpolation(
            self.unit_cell.fractionalize(ai.xyz))
        if(ai.name.strip() in exclude): continue
        if(ai.element.strip().upper() in ["H","S","SE"]): continue
        for j, aj in enumerate(atoms):
          if i==j: continue
          if(aj.name.strip() in exclude): continue
          if(aj.element.strip().upper() in ["H","S","SE"]): continue
          d = ai.distance(aj)
          if d < 1.6:
            pair = [i,j]
            pair.sort()
            if(not pair in self.pairs):
              self.pairs.append(pair)
      # Set monitor
      id_str=""
      if(self.residue.parent() is not None and
         self.residue.parent().parent() is not None):
         id_str+="chain: %s"%(self.residue.parent().parent().id)
      id_str+=" residue: %s %s"%(self.residue.resname, self.residue.resseq.strip())
      if(len(self.co.clusters)>1):
        msel = flex.size_t(flatten(self.co.clusters[1:][0].vector))
      else:
        msel = flex.size_t()
      self.m = monitor(
        id_str    = id_str,
        selection = msel,
        map_data  = self.target_map,
        unit_cell = self.unit_cell,
        weights   = self.weights,
        pairs     = self.pairs,
        cmv       = self.cmv,
        rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
        log       = self.log)
      self.m.add(residue = self.residue, state = "start")
    if(self.target_map is None):
      assert not backbone_sample
    # Actual calculations
    if(self.residue.resname == "PRO"): # Special case
      self.fit_proline()
    else:
      self.chi_angles = self.rotamer_manager.get_chi_angles(
        resname = self.residue.resname)
      if(len(self.co.clusters)>0):
        if(backbone_sample):
          self.fit_c_beta(c_beta_rotation_cluster = self.co.clusters[0])
        self.fit_side_chain(clusters = self.co.clusters[1:])
    if(self.m is not None):
      self.m.finalize(residue = self.residue)
      # Too bulky, but very useful. Use for debugging only.
      #self.m.show()

  def get_target_value(self, sites_cart, selection=None, target_map=None):
    if(target_map is None): target_map = self.target_map
    if(selection is None):
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = target_map,
        sites_cart  = sites_cart)
    else:
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = target_map,
        sites_cart  = sites_cart,
        selection   = selection)

  def get_rotamer_iterator(self):
    return mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
      mon_lib_srv = self.mon_lib_srv,
      residue     = self.residue)

  def fit_proline(self):
    """
    PRO is a special case. Just sample two possible rotamers.
    Skip if map isn't available: can't do much in this case!
    """
    if(self.target_map is None): return
    rotamer_iterator = self.get_rotamer_iterator()
    if rotamer_iterator is None:
      id_str="chain: %s"%(self.residue.parent().parent().id)
      id_str+=" residue: %s %s"%(
        self.residue.resname, self.residue.resseq.strip())
      print("Corrupt residue: %s >>> skipping"%id_str, file=self.log)
      return None
    scorer = mmtbx.refinement.real_space.score3(
      unit_cell    = self.unit_cell,
      target_map   = self.target_map,
      residue      = self.residue,
      rotamer_eval = self.rotamer_manager.rotamer_evaluator)
    scorer.reset(sites_cart = self.residue.atoms().extract_xyz())
    for rotamer, sites_cart in rotamer_iterator:
      scorer.update(sites_cart = sites_cart)
    self.residue.atoms().set_xyz(new_xyz=scorer.sites_cart)
    self.m.add(residue = self.residue, state = "fitting")

  def fit_side_chain(self, clusters):
    rotamer_iterator = self.get_rotamer_iterator()
    if(rotamer_iterator is None): return
    selection_clash = self.co.clash_eval_selection
    selection_rsr   = self.co.rsr_eval_selection
    if(self.target_map is not None):
      start_target_value = self.get_target_value(
        sites_cart = self.residue.atoms().extract_xyz(),
        selection  = selection_rsr)
    sites_cart_start = self.residue.atoms().extract_xyz()
    sites_cart_first_rotamer = list(rotamer_iterator)[0][1]
    # From this point on the coordinates in residue are to initial rotamer!
    self.residue.atoms().set_xyz(sites_cart_first_rotamer)
    axes = []
    atr = []
    for i, angle in enumerate(self.chi_angles[0]):
      cl = clusters[i]
      axes.append(flex.size_t(cl.axis))
      atr.append(flex.size_t(cl.atoms_to_rotate))
    #
    if(self.target_map is not None and self.xyzrad_bumpers is not None):
      # Get reference map values
      ref_map_vals = flex.double()
      for a in self.residue.atoms():
        key = "%s_%s_%s"%(
          a.parent().parent().parent().id, a.parent().resname,
          a.name.strip())
        ref_map_vals.append(self.cmv[key])
      # Get radii
      radii = mmtbx.refinement.real_space.get_radii(
        residue = self.residue, vdw_radii = self.vdw_radii)
      # Exclude rotatable H from clash calculation
      tmp = flex.size_t()
      for i in selection_clash:
        if(self.rotatable_hd[self.residue.atoms()[i].i_seq]): continue
        tmp.append(i)
      selection_clash = tmp[:]
      # Ad hoc: S or SE have larger peaks!
      if(self.residue.resname in ["MET","MSE"]): scale=100
      else:                                      scale=3
      if self.trust_map_values_real:
        scaleup=1
        scaledown=1
      else:
        scaleup=10000
        scaledown=0
      moving = ext.moving(
        sites_cart       = self.residue.atoms().extract_xyz(),
        sites_cart_start = sites_cart_start,
        radii            = radii,
        weights          = self.weights,
        bonded_pairs     = self.pairs,
        ref_map_max      = ref_map_vals * scale*scaleup,
        ref_map_min      = ref_map_vals / 10   *scaledown)
      #
      ro = ext.fit(
        fixed                    = self.xyzrad_bumpers,
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        density_map              = self.target_map,
        moving                   = moving,
        unit_cell                = self.unit_cell,
        selection_clash          = selection_clash,
        selection_rsr            = selection_rsr, # select atoms to compute map target
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    elif(self.target_map is not None and self.xyzrad_bumpers is None):
      ro = ext.fit(
        target_value             = start_target_value,
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        density_map              = self.target_map,
        all_points               = self.residue.atoms().extract_xyz(),
        unit_cell                = self.unit_cell,
        selection                = selection_rsr,
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    else:
      ro = ext.fit(
        sites_cart_start         = sites_cart_start.deep_copy(),
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        all_points               = self.residue.atoms().extract_xyz(),
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    sites_cart_result = ro.result()
    if(sites_cart_result.size()>0):
      dist = None
      if(self.accept_only_if_max_shift_is_smaller_than is not None):
        dist = flex.max(flex.sqrt((sites_cart_start - sites_cart_result).dot()))
      if(dist is None):
        self.residue.atoms().set_xyz(sites_cart_result)
      else:
        if(dist is not None and
           dist < self.accept_only_if_max_shift_is_smaller_than):
          self.residue.atoms().set_xyz(sites_cart_result)
        else:
          self.residue.atoms().set_xyz(sites_cart_start)
    else:
      self.residue.atoms().set_xyz(sites_cart_start)
    if(self.m): self.m.add(residue = self.residue, state = "fitting")
#    # tune up
    if(self.target_map is not None):
      tune_up(
        target_map           = self.target_map,
        residue              = self.residue,
        mon_lib_srv          = self.mon_lib_srv,
        rotamer_manager      = self.rotamer_manager.rotamer_evaluator,
        unit_cell            = self.unit_cell,
        monitor = self.m,
        torsion_search_start = -30,
        torsion_search_stop  = 30,
        torsion_search_step  = 1)

  def fit_c_beta(self, c_beta_rotation_cluster):
    selection = flex.size_t(c_beta_rotation_cluster.selection)
    sites_cart = self.residue.atoms().extract_xyz()
    sites_cart_start = sites_cart.deep_copy() # XXX
    start_target_value = self.get_target_value(
      sites_cart = sites_cart,
      selection  = selection,
      target_map = self.target_map_for_cb)
    ro = ext.fit(
      target_value             = start_target_value+1.e-6,
      axes                     = [c_beta_rotation_cluster.axis],
      rotatable_points_indices = [c_beta_rotation_cluster.atoms_to_rotate],
      angles_array             = [[i*math.pi/180] for i in range(-20,21,1)],
      density_map              = self.target_map_for_cb,
      all_points               = sites_cart,
      unit_cell                = self.unit_cell,
      selection                = selection,
      sin_table                = self.sin_cos_table.sin_table,
      cos_table                = self.sin_cos_table.cos_table,
      step                     = self.sin_cos_table.step,
      n                        = self.sin_cos_table.n)
    sites_cart_result = ro.result()
    if(sites_cart_result.size()>0):
      self.residue.atoms().set_xyz(sites_cart_result)
    else:
      self.residue.atoms().set_xyz(sites_cart_start)

class run_with_minimization(object):
  def __init__(self,
               target_map,
               residue,
               vdw_radii,
               xray_structure,
               mon_lib_srv,
               rotamer_manager,
               # This is cctbx.geometry_restraints.manager.manager
               geometry_restraints_manager,
               real_space_gradients_delta,
               selection_radius = 5,
               rms_bonds_limit = 0.03, # XXX probably needs to be much lower
               rms_angles_limit = 3.0, # XXX
               backbone_sample_angle=None,
               cmv = None,
               allow_modified_residues=False):
    adopt_init_args(self, locals())
    # load rotamer manager
    self.rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load(
      rotamers="favored")
    # pre-compute sin and cos tables
    self.sin_cos_table = scitbx.math.sin_cos_table(n=10000)
    self.backbone_atom_names = ["N", "CA", "O", "CB", "C"]
    self.residue_iselection = self.residue.atoms().extract_i_seq()
    assert (not self.residue_iselection.all_eq(0))
    self.residue_selection = flex.bool(
      xray_structure.scatterers().size(), self.residue_iselection)
    self.residue_backbone_selection = flex.size_t()
    for atom in self.residue.atoms():
      if(atom.name.strip() in self.backbone_atom_names):
        self.residue_backbone_selection.append(atom.i_seq)
    self.residue_backbone_selection = flex.bool(
      xray_structure.scatterers().size(), self.residue_backbone_selection)
    self.target_map_work = target_map
    self.target_map_orig = target_map.deep_copy()
    self.fit_backbone()
    negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
      xray_structure          = self.xray_structure,
      selection_within_radius = self.selection_radius,
      iselection              = self.residue.atoms().extract_i_seq())
    self.target_map_work = mmtbx.refinement.real_space.\
      negate_map_around_selected_atoms_except_selected_atoms(
        xray_structure   = self.xray_structure,
        map_data         = target_map,
        negate_selection = negate_selection,
        atom_radius      = 1.5)
    self.fit_rotamers()

  def fit_backbone(self):
    # move in place (pure geometry regularizaition of residue in question)
    self.real_space_refine(optimize_weight=False, start_trial_weight_value=0)
    # fit n-c-o-ca-cb only (ignore side chain!). XXX BAD: amino-acid specific!
    self.grid_sample_around_c_n_axis()
    # fine-tune
    self.real_space_refine(optimize_weight=True, start_trial_weight_value=50)

  def fit_rotamers(self):
    sps = self.xray_structure.special_position_settings()
    mmtbx.refinement.real_space.fit_residue.run(
      vdw_radii         = self.vdw_radii,
      target_map        = self.target_map_work,
      target_map_for_cb = self.target_map_orig,
      mon_lib_srv       = self.mon_lib_srv,
      unit_cell         = self.xray_structure.unit_cell(),
      residue           = self.residue,
      sin_cos_table     = self.sin_cos_table,
      cmv               = self.cmv,
      rotamer_manager   = self.rotamer_manager)
    sites_cart_poor = self.xray_structure.sites_cart()
    sites_cart_poor.set_selected(self.residue_iselection,
      self.residue.atoms().extract_xyz())
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_poor)

  def grid_sample_around_c_n_axis(self):
    sps = self.xray_structure.special_position_settings()
    scorer = mmtbx.refinement.real_space.score(
      target_map = self.target_map_work,
      residue    = self.residue,
      unit_cell  = self.xray_structure.unit_cell())
    def get_cluster(self):
      axis=[]
      atoms_to_rotate=[]
      use_in_target_selection = flex.size_t()
      counter = 0
      for atom in self.residue.atoms():
        if(atom.name.strip() in ["N", "C"]):
          axis.append(counter)
        else:
          atoms_to_rotate.append(counter)
        if(atom.name.strip() in self.backbone_atom_names):
          use_in_target_selection.append(counter)
        counter += 1
      return mmtbx.refinement.real_space.cluster(
        axis            = axis,
        atoms_to_rotate = atoms_to_rotate,
        selection       = use_in_target_selection)
    cl = get_cluster(self)
    residue_sites_cart = self.residue.atoms().extract_xyz()
    scorer.reset(
      sites_cart = residue_sites_cart,
      selection  = cl.selection)
    angle_start = 0
    angle_end = 360
    if (self.backbone_sample_angle is not None):
      assert (self.backbone_sample_angle > 0)
      angle_start = - self.backbone_sample_angle
      angle_end = self.backbone_sample_angle
    mmtbx.refinement.real_space.torsion_search(
      clusters   = [cl],
      sites_cart = residue_sites_cart,
      scorer     = scorer,
      start      = 0,
      stop       = 360,
      step       = 1)
    self.residue.atoms().set_xyz(new_xyz=scorer.sites_cart)
    selection = self.residue.atoms().extract_i_seq()
    sites_cart_poor = self.xray_structure.sites_cart()
    sites_cart_poor.set_selected(selection, scorer.sites_cart)
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_poor)

  def real_space_refine(self, optimize_weight, start_trial_weight_value):
    brm = individual_sites.box_refinement_manager(
      xray_structure              = self.xray_structure,
      target_map                  = self.target_map_work,
      geometry_restraints_manager = self.geometry_restraints_manager,
      real_space_gradients_delta  = 1./4,
      max_iterations              = 500)
    brm.refine(
      selection                = self.residue_selection,
      optimize_weight          = optimize_weight,
      start_trial_weight_value = start_trial_weight_value,
      selection_buffer_radius  = self.selection_radius,
      box_cushion              = 2,
      rms_bonds_limit          = self.rms_bonds_limit,
      rms_angles_limit         = self.rms_angles_limit)
    self.xray_structure = brm.xray_structure
    self.residue.atoms().set_xyz(brm.sites_cart.select(self.residue_iselection))

def get_rotamer_iterator(mon_lib_srv, residue):
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    fine_sampling = True,
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  if (rotamer_iterator is None):
    return None
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

class tune_up(object):
  def __init__(self,
               target_map,
               residue,
               mon_lib_srv,
               rotamer_manager,
               unit_cell,
               exclude_hd           = True,
               monitor              = None,
               torsion_search_start = -20,
               torsion_search_stop  = 20,
               torsion_search_step  = 2):
    adopt_init_args(self, locals())
    self.clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = False).clusters
    self.score_residue = mmtbx.refinement.real_space.score3(
      unit_cell    = self.unit_cell,
      target_map   = self.target_map,
      residue      = self.residue,
      rotamer_eval = self.rotamer_manager,
      exclude_hd   = self.exclude_hd)
    mmtbx.refinement.real_space.torsion_search(
      clusters   = self.clusters,
      sites_cart = self.residue.atoms().extract_xyz(),
      scorer     = self.score_residue,
      start      = self.torsion_search_start,
      stop       = self.torsion_search_stop,
      step       = self.torsion_search_step)
    self.residue.atoms().set_xyz(new_xyz=self.score_residue.sites_cart)
    if(monitor is not None):
      monitor.add(residue = self.residue, state = "tuneup")

#
# These functions are not used anywhere. And not tested anymore.
# They are here as an example of correct backrub move, according to
# original paper https://doi.org/10.1016/j.str.2005.10.007
# Unfortunately, for proper backrub move we need previous and next residues,
# but current code is build under assumption that one residue is enough for
# rotamer fitting. One will have to reconsider this idea and do some changes
# to make it possible to do proper backrub move.
#
def _find_theta(ap1, ap2, cur_xyz, needed_xyz):
  from mmtbx.building.loop_closure.ccd import ccd_python
  f, s_home, r_norm, r_home = ccd_python._get_f_r_s(
      axis_point_1=ap1,
      axis_point_2=ap2,
      moving_coor=cur_xyz,
      fixed_coor=needed_xyz)
  b = list(2*r_norm*(f.dot(r_home)))[0]
  c = list(2*r_norm*(f.dot(s_home)))[0]
  znam = math.sqrt(b*b+c*c)
  sin_alpha = c/znam
  cos_alpha = b/znam
  alpha = math.atan2(sin_alpha, cos_alpha)
  return math.degrees(alpha)

def backrub_move(
    prev_res,
    cur_res,
    next_res,
    angle,
    move_oxygens=False,
    accept_worse_rama=False,
    rotamer_manager=None,
    rama_manager=None):
  import boost_adaptbx.boost.python as bp
  ext = bp.import_ext("mmtbx_validation_ramachandran_ext")
  from mmtbx_validation_ramachandran_ext import rama_eval
  from scitbx.matrix import rotate_point_around_axis
  from mmtbx.conformation_dependent_library.multi_residue_class import ThreeProteinResidues, \
      RestraintsRegistry

  if abs(angle) < 1e-4:
    return
  if prev_res is None or next_res is None:
    return
  saved_res = [{},{},{}]
  for i, r in enumerate([prev_res, cur_res, next_res]):
    for a in r.atoms():
      saved_res[i][a.name.strip()] = a.xyz
  if rotamer_manager is None:
    rotamer_manager = RotamerEval()
  prev_ca = prev_res.find_atom_by(name=" CA ")
  cur_ca = cur_res.find_atom_by(name=" CA ")
  next_ca = next_res.find_atom_by(name=" CA ")
  if prev_ca is None or next_ca is None or cur_ca is None:
    return
  atoms_to_move = []
  atoms_to_move.append(prev_res.find_atom_by(name=" C  "))
  atoms_to_move.append(prev_res.find_atom_by(name=" O  "))
  for atom in cur_res.atoms():
    atoms_to_move.append(atom)
  atoms_to_move.append(next_res.find_atom_by(name=" N  "))
  for atom in atoms_to_move:
    assert atom is not None
    new_xyz = rotate_point_around_axis(
        axis_point_1 = prev_ca.xyz,
        axis_point_2 = next_ca.xyz,
        point        = atom.xyz,
        angle        = angle,
        deg          = True)
    atom.xyz = new_xyz
  if move_oxygens:
    registry = RestraintsRegistry()
    if rama_manager is None:
      rama_manager = rama_eval()
    tpr = ThreeProteinResidues(geometry=None, registry=registry)
    tpr.append(prev_res)
    tpr.append(cur_res)
    tpr.append(next_res)
    phi_psi_angles = tpr.get_phi_psi_angles()
    rama_key = tpr.get_ramalyze_key()
    ev_before = rama_manager.evaluate_angles(rama_key, phi_psi_angles[0], phi_psi_angles[1])
    theta1 = _find_theta(
        ap1 = prev_ca.xyz,
        ap2 = cur_ca.xyz,
        cur_xyz = prev_res.find_atom_by(name=" O  ").xyz,
        needed_xyz = saved_res[0]["O"])
    theta2 = _find_theta(
        ap1 = cur_ca.xyz,
        ap2 = next_ca.xyz,
        cur_xyz = cur_res.find_atom_by(name=" O  ").xyz,
        needed_xyz = saved_res[1]["O"])
    for a in [prev_res.find_atom_by(name=" C  "),
        prev_res.find_atom_by(name=" O  "),
        cur_res.find_atom_by(name=" C  ")]:
      new_xyz = rotate_point_around_axis(
              axis_point_1 = prev_ca.xyz,
              axis_point_2 = cur_ca.xyz,
              point        = a.xyz,
              angle        = theta1,
              deg          = True)
      a.xyz = new_xyz
    for a in [cur_res.find_atom_by(name=" C  "),
        cur_res.find_atom_by(name=" O  "),
        next_res.find_atom_by(name=" N  ")]:
      new_xyz = rotate_point_around_axis(
              axis_point_1 = cur_ca.xyz,
              axis_point_2 = next_ca.xyz,
              point        = a.xyz,
              angle        = theta2,
              deg          = True)
      a.xyz = new_xyz
    phi_psi_angles = tpr.get_phi_psi_angles()
    rama_key = tpr.get_ramalyze_key()
    ev_after = rama_manager.evaluate_angles(rama_key, phi_psi_angles[0], phi_psi_angles[1])
    if ev_before > ev_after and not accept_worse_rama:
      for a in [prev_res.find_atom_by(name=" C  "),
          prev_res.find_atom_by(name=" O  "),
          cur_res.find_atom_by(name=" C  ")]:
        new_xyz = rotate_point_around_axis(
                axis_point_1 = prev_ca.xyz,
                axis_point_2 = cur_ca.xyz,
                point        = a.xyz,
                angle        = -theta1,
                deg          = True)
        a.xyz = new_xyz
      for a in [cur_res.find_atom_by(name=" C  "),
          cur_res.find_atom_by(name=" O  "),
          next_res.find_atom_by(name=" N  ")]:
        new_xyz = rotate_point_around_axis(
                axis_point_1 = cur_ca.xyz,
                axis_point_2 = next_ca.xyz,
                point        = a.xyz,
                angle        = -theta2,
                deg          = True)
        a.xyz = new_xyz
