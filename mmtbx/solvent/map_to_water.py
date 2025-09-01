from __future__ import division, print_function
from cctbx.array_family import flex
import iotbx.pdb
import mmtbx.model
import iotbx.mrcfile
import boost_adaptbx.boost.python as bp
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
import time
from cctbx import maptbx
from libtbx import adopt_init_args
from libtbx.utils import user_plus_sys_time
from cctbx import crystal, maptbx, xray, adptbx
from iotbx.pdb.utils import all_chain_ids
import scitbx.math
from libtbx.utils import Sorry

import mmtbx.maps.correlation
from cctbx import miller
import mmtbx.refinement.real_space.adp

from libtbx.utils import null_out

elements=["N","O","S","Se","Fe","Mn","Mg","Zn","Ca","Cd","Cu","Co"]
selection_string_interaction = " or ".join(["element %s"%e for e in elements])

hoh_str = "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f  1.00  0.00           O"

def _debug_show_all_plots(
      map_data,
      unit_cell,
      center_cart,
      radius,
      plot_number,
      s_angle_sampling_step = 60,
      t_angle_sampling_step = 60):
  rho = flex.double()
  vecs = []
  one_d = flex.double()
  for s in range(0,360,s_angle_sampling_step):
    for t in range(0,360,t_angle_sampling_step):
      xc,yc,zc = scitbx.math.point_on_sphere(r=radius, s_deg=s, t_deg=t,
        center=center_cart)
      xf,yf,zf = unit_cell.fractionalize([xc,yc,zc])
      rho.append(map_data.tricubic_interpolation([xf,yf,zf]))
      ###
      v = flex.double()
      r = flex.double()
      for p in [p/20 for p in range(0,21)]:
        xp = center_cart[0] + p * (xc-center_cart[0])
        yp = center_cart[1] + p * (yc-center_cart[1])
        zp = center_cart[2] + p * (zc-center_cart[2])
        xpf,ypf,zpf = unit_cell.fractionalize([xp,yp,zp])
        v.append(map_data.tricubic_interpolation([xpf,ypf,zpf]))
        r.append(p)
      vecs.append([r,v])
      one_d.extend(v)
      ###
  import matplotlib.pyplot as plt
  import matplotlib as mpl
  mpl.use('Agg')
  fig, axes = plt.subplots(ncols=6, nrows=6)
  for i, ax in enumerate(axes.flatten()):
    ax.plot(vecs[i][0], vecs[i][1])
    ax.set_ylim(flex.min(one_d),flex.max(one_d))
    ax.set_xticks([])
    ax.set_yticks([])
  fig.savefig("fig_%s.png"%plot_number)
  ###
  del plt
  del mpl
  return None

class msg_accumulator(object):
  def __init__(self, log=None):
    self.log = log
    self.messages = []
    self.prefix=""

  def set_prefix(self, prefix):
    self.prefix=prefix

  def add(self, msg):
    msg = "%s%s"%(self.prefix, msg)
    self.messages.append(msg)
    if(self.log is not None): print(msg, file=self.log)

  def show(self, log=None):
    assert [self.log, log].count(None) == 1
    if(log is not None):
      assert self.log is None
    else:
      assert self.log is not None
      log = self.log
    for msg in self.messages:
      print(msg, file=log)

def models_as_chains(model):
  ph = model.get_hierarchy()
  r = iotbx.pdb.hierarchy.root()
  m = iotbx.pdb.hierarchy.model()
  for m_ in ph.models():
    for c_ in m_.chains():
      c_ = c_.detached_copy()
      m.append_chain(c_)
  r.append_model(m)
  #
  wat_sel = r.atom_selection_cache().selection("water")
  mm = r.select(~wat_sel)
  for m in r.select(wat_sel).models():
    for c in m.chains():
      c_ = c.detached_copy()
      mm.models()[0].append_chain(c_)
  r = mm
  #
  new_model = mmtbx.model.manager(
    model_input = None,
    pdb_hierarchy = r,
    crystal_symmetry = model.crystal_symmetry())
  new_model.set_shift_cart(shift_cart = model.shift_cart())
  return new_model

def remove_clashing(model, dist_min):
  h = model.get_hierarchy()
  # Remove overlapping water-water
  wsel = h.atom_selection_cache().selection(string="water")
  wh = h.select(wsel)
  sel = []
  for i, ai in enumerate(wh.atoms()):
    for j, aj in enumerate(wh.atoms()):
      if(j>=i or ai.i_seq==aj.i_seq): continue
      d = ai.distance(aj)
      if(d<dist_min):
        sel_ = [ai.i_seq, aj.i_seq]
        sel_.sort()
        sel.append(sel_)
  sel_clash = flex.size_t([s[0] for s in sel])
  # Remove overlapping water-nonwater
  nonwa = h.select(~wsel).atoms()
  sel = flex.size_t()
  for aw in wh.atoms():
    for ap in nonwa:
      d = aw.distance(ap)
      if(d<dist_min):
        sel.append(aw.i_seq)
        break
  sel_clash.extend(sel)
  #
  return model.select(~flex.bool(h.atoms_size(), sel_clash))

def attribute_water_to_chains(model):
  h = model.get_hierarchy()
  wsel = h.atom_selection_cache().selection(string="water")
  wh = h.select(wsel)
  nonw = h.select(~wsel)
  nonwa = nonw.atoms()
  #
  # chain ID <> max residue number mapping
  can = {}
  for c in nonw.chains():
    last = list(c.residue_groups())[-1]
    can.setdefault(c.id, []).append(last.resseq_as_int())
  for k in can.keys():
    can[k] = max(can[k])
  #
  dic = {}
  for aw in wh.atoms():
    d_min    = 1.e9
    #b_iso    = None # This un-does ADP refinement!!!
    chain_id = None
    for ap in nonwa:
      d = aw.distance(ap)
      if(d<d_min):
        d_min = d
        #b_iso    = ap.b
        chain_id = ap.parent().parent().parent().id
    #aw.set_b(b_iso)
    dic.setdefault(chain_id, []).append(aw)
  #
  pdb_model = nonw.only_model()
  for chain_id, atoms in zip(dic.keys(), dic.values()):
    new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
    for i_seq, new_atom in enumerate(atoms):
      new_atom_group = iotbx.pdb.hierarchy.atom_group(altloc="",
        resname="HOH")
      new_atom_group.append_atom(atom=new_atom.detached_copy())
      i_seq_ = i_seq + 1 + can[chain_id]
      new_residue_group = iotbx.pdb.hierarchy.residue_group(
        resseq=iotbx.pdb.resseq_encode(value=i_seq_), icode=" ")
      new_residue_group.append_atom_group(atom_group=new_atom_group)
      new_chain.append_residue_group(residue_group=new_residue_group)
    pdb_model.append_chain(chain=new_chain)
  #
  new_model = mmtbx.model.manager(
    model_input = None,
    pdb_hierarchy = nonw,
    crystal_symmetry = model.crystal_symmetry())
  new_model.set_shift_cart(model.shift_cart())
  new_model.get_hierarchy().atoms().reset_i_seq()
  return new_model

def write_map_file(map_data, cs, file_name): # FOR DEBUGGING XXX
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name   = file_name,
    unit_cell   = cs.unit_cell(),
    space_group = cs.space_group(),
    map_data    = map_data,
    labels      = flex.std_string([""]))

class run(object):
  def __init__(self, map_model_manager, params, log):
    self.mmm        = map_model_manager
    self.log        = log
    self.params     = params
    self.ma         = msg_accumulator(log = self.log)
    self.total_time = 0
    if(not self.params.keep_input_water):
      self._call(self._remove_input_water, "Removing water in input model")
    if("per_chain" in self.params.mode):
      self._call(self._run_per_chain, "Process per chain (mode=per_chain)")
    elif("whole" in params.mode):
      self._call(self._run_whole, "Process whole (mode=whole)")
    else: raise Sorry("Invalid mode:", params.mode)
    self._call(self._attribute_water_to_chains, "Attribute water to chains")
    self._call(self._summary, "Summary")

  def _call(self, func, msg):
    timer = user_plus_sys_time()
    self.ma.add(msg)
    func()
    t = timer.elapsed()
    self.total_time += t
    self.ma.add("  time (s): %s (total time: %s)"%(("%8.3f"%t).strip(),
      ("%8.3f"%self.total_time).strip()))

  def _summary(self):
    wsel = self.mmm.model().selection(string="water")
    self.ma.add("  Number of water added: %d"%wsel.count(True))

  def _remove_input_water(self):
    wsel = self.mmm.model().selection(string="water")
    nwater = wsel.count(True)
    self.ma.add("  number of water in input model to remove: %d"%nwater)
    self.mmm.set_model(model = self.mmm.model().select(~wsel))

  def _attribute_water_to_chains(self):
    model = attribute_water_to_chains(model = self.mmm.model())
    self.mmm.set_model(model = model)

  def _run_whole(self):
    self.ma.set_prefix(prefix="  ")
    self.mmm = run_one(
      ma                       = self.ma,
      map_model_manager        = self.mmm,
      resolution               = self.params.resolution,
      dist_min                 = self.params.dist_min,
      dist_max                 = self.params.dist_max,
      step                     = self.params.step,
      water_residues_ratio     = self.params.water_residues_ratio,
      map_threshold_scale      = self.params.map_threshold_scale,
      scc                      = self.params.scc,
      sphericity_filter        = self.params.sphericity_filter,
      cc_mask_filter           = self.params.cc_mask_filter,
      cc_mask_filter_threshold = self.params.cc_mask_filter_threshold,
      cc_mask_threshold_interacting_atoms = self.params.cc_mask_threshold_interacting_atoms,
      atom_radius              = self.params.atom_radius,
      scattering_table         = self.params.scattering_table,
      debug                    = self.params.debug,
      total_time               = self.total_time,
      log                      = self.log).map_model_manager
    self.ma.set_prefix(prefix="  ")

  def _run_per_chain(self):
    box_info = self.mmm.split_up_map_and_model_by_chain(
      mask_around_unselected_atoms=False, box_cushion=4)
    for i_box, box_mmm in enumerate(box_info.mmm_list):
      self.ma.set_prefix(prefix="  ")
      chain_ids = []
      for c in box_mmm.model().get_hierarchy().chains():
        chain_ids.append(c.id)
      chain_ids = list(set(chain_ids))
      assert len(chain_ids) == 1
      self.ma.add("Working on chain %s:"%chain_ids[0])
      self.ma.set_prefix(prefix="    ")
      box_mmm = run_one(
        ma                       = self.ma,
        map_model_manager        = box_mmm,
        resolution               = self.params.resolution,
        dist_min                 = self.params.dist_min,
        dist_max                 = self.params.dist_max,
        step                     = self.params.step,
        water_residues_ratio     = self.params.water_residues_ratio,
        map_threshold_scale      = self.params.map_threshold_scale,
        sphericity_filter        = self.params.sphericity_filter,
        cc_mask_filter           = self.params.cc_mask_filter,
        cc_mask_filter_threshold = self.params.cc_mask_filter_threshold,
        cc_mask_threshold_interacting_atoms = self.params.cc_mask_threshold_interacting_atoms,
        scc                      = self.params.scc,
        atom_radius              = self.params.atom_radius,
        scattering_table         = self.params.scattering_table,
        debug                    = self.params.debug,
        total_time               = self.total_time,
        log                      = self.log).map_model_manager
      self.ma.set_prefix(prefix="")
    self.mmm.merge_split_maps_and_models(
      box_info                   = box_info,
      replace_coordinates        = True,
      allow_changes_in_hierarchy = True,
      replace_u_aniso            = False)
    model = models_as_chains(model = self.mmm.model())
    model = remove_clashing(
      model    = model,
      dist_min = self.params.dist_min)
    self.mmm.set_model(model = model)

class run_one(object):
  def __init__(self,
               ma, # message_accumulator
               map_model_manager,
               resolution=None,
               dist_min=2.0,
               dist_max=3.2,
               step=0.3,
               water_residues_ratio=1.5,
               map_threshold_scale=1.5,
               sphericity_filter=True,
               cc_mask_filter = False,
               cc_mask_filter_threshold = 0.5,
               scc=0.5,
               atom_radius=2,
               cc_mask_threshold_interacting_atoms=0.5,
               scattering_table="electron",
               debug=False,
               log=None,
               total_time=0,
               prefix=""):
    # Common parameters / variables
    self.model = map_model_manager.model()
    self.n_max_water = int(self.model.get_hierarchy().overall_counts().
                           n_residues*water_residues_ratio)
    map_data = map_model_manager.map_manager().map_data()
    adopt_init_args(self, locals())
    self.total_time = total_time
    if(map_data.accessor().origin()!=(0,0,0)):
      raise Sorry("Map origin must be zero.")
    self.ma = ma
    self.map_data = self.map_data.deep_copy() # will be changed in-place!
    self.map_data_resampled = None
    self.unit_cell = self.model.crystal_symmetry().unit_cell()
    self.interaction_selection = self.model.selection(
      selection_string_interaction)
    self.sites_frac_interaction = self.model.get_sites_frac().select(
      self.interaction_selection)
    self.xrs_water = None
    self.cutoff = None
    self._wrote_file=0
    # Calculations
    self._call(self._massage_map         , "Normalize and re-sample input map")
    self._call(self._get_cutoff          , "Get cutoff")
    self._call(self._mask_out            , "Mask out molecule/solvent regions")
    self._call(self._find_peaks          , "Find peaks")
    self._call(self._filter_by_distance  , "Filter peaks by distance")
    if(self.sphericity_filter):
      self._call(self._filter_by_sphericity, "Filter peaks by sphericity")
      self._call(self._filter_by_distance  , "Filter peaks by distance")
    if(self.cc_mask_filter and
       (self.model.solvent_selection().iselection().size()>0 or
        self.xrs_water is not None and self.xrs_water.scatterers().size()>0)):
      self._call(self._refine_water_adp     , "Refine ADP")
      self._call(self._filter_by_map_model_cc, "Filter peaks by CC_mask")
      self._call(self._filter_by_distance  , "Filter peaks by distance")
      self._call(self._filter_by_map_model_cc_and_maxwater, "Filter peaks by CC_mask and number")
      self._call(self._filter_by_distance  , "Filter peaks by distance")
    self._call(self._append_to_model     , "Add to model and reset internals")
    if(self.debug):
      self._call(self._show_peak_profiles, "Show peak profiles")

  def _assert_same_model(self):
    assert self.model.is_same_model(other = self.map_model_manager.model())

  def _call(self, func, msg):
    timer = user_plus_sys_time()
    self.ma.add(msg)
    self._assert_same_model()
    func()
    t = timer.elapsed()
    self.total_time += t
    self.ma.add("  time (s): %s (total time: %s)"%(("%8.3f"%t).strip(),
      ("%8.3f"%self.total_time).strip()))

  def _map_mmm_str(self):
    return " ".join([("%8.3f"%m).strip() for m in [
      flex.min(self.map_data),
      flex.max(self.map_data),
      flex.mean(self.map_data)]])

  def _write_map(self, file_name):
    iotbx.mrcfile.write_ccp4_map(
      file_name   = file_name,
      unit_cell   = self.unit_cell,
      space_group = self.model.crystal_symmetry().space_group(),
      map_data    = self.map_data.as_double(),
      labels      = flex.std_string([""]))

  def _massage_map(self):
    self.ma.add("  input map (min,max,mean): %s"%self._map_mmm_str())
    #
    self.map_data = self.map_data - flex.mean(self.map_data)
    self.map_data = self.map_data.set_selected(self.map_data < 0, 0)
    sd = self.map_data.sample_standard_deviation()
    assert sd != 0
    self.map_data = self.map_data / sd
    self.ma.add("  re-scaled map (min,max,mean): %s"%self._map_mmm_str())
    if(self.debug):
      self._write_map(file_name = "rescaled.mrc")
    #
    a,b,c = self.unit_cell.parameters()[:3]
    n_real = self.map_data.accessor().all()
    sx, sy, sz = a/n_real[0], b/n_real[1], c/n_real[2]
    self.ma.add("  input map dimenstions: %d %d %d"%n_real)
    self.ma.add("  input map grid steps (A): %s %s %s"%(
      ("%6.3f"%sx).strip(), ("%6.3f"%sy).strip(), ("%6.3f"%sz).strip()))
    self.ma.add("  input map (min,max,mean): %s"%self._map_mmm_str())
    if(self.debug):
      self._write_map(file_name = "resampled.mrc")
    self.map_data_resampled = self.map_data.deep_copy()

  def _write_pdb_file(self, file_name_prefix):
    return
    file_name = file_name_prefix + "_%d.pdb"%self._wrote_file
    with open(file_name, "w") as fo:
      for i, site_frac in enumerate(sites_frac):
        site_cart = self.unit_cell.orthogonalize(site_frac)
        print(hoh_str%(i,i,site_cart[0],site_cart[1],site_cart[2]), file=fo)
    self._wrote_file += 1

  def _get_cutoff(self):
    map_at_atoms = flex.double()
    for site_frac in self.model.get_sites_frac():
      map_at_atoms.append( self.map_data_resampled.tricubic_interpolation(site_frac) )
    mean = flex.mean(map_at_atoms)
    self.cutoff = mean * self.map_threshold_scale
    self.ma.add("  mean density at atom centers: %8.3f"%mean)
    self.ma.add("  density cutoff for peak search: %8.3f"%self.cutoff)

  def _mask_out(self, rad_inside = 1.0, rad_outside = 6.0):
    # Zero inside
    radii = flex.double(self.model.size(), rad_inside)
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = self.model.get_sites_frac(),
      unit_cell                   = self.unit_cell,
      n_real                      = self.map_data.all(),
      mask_value_inside_molecule  = 0,
      mask_value_outside_molecule = 1,
      radii                       = radii)
    self.map_data = self.map_data * mask
    # Zero outside
    radii = flex.double(self.sites_frac_interaction.size(), rad_outside)
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = self.sites_frac_interaction,
      unit_cell                   = self.unit_cell,
      n_real                      = self.map_data.all(),
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = radii)
    self.map_data = self.map_data * mask
    #self._write_map(file_name = "masked.mrc")

  def _find_peaks(self, min_peak_peak_distance = 1.5):
    cg = maptbx.crystal_gridding(
      space_group_info      = self.model.crystal_symmetry().space_group_info(),
      symmetry_flags        = maptbx.use_space_group_symmetry,
      unit_cell             = self.unit_cell,
      pre_determined_n_real = self.map_data.accessor().all())
    cgt = maptbx.crystal_gridding_tags(gridding = cg)
    peak_search_parameters = maptbx.peak_search_parameters(
      peak_search_level      = 3,
      max_peaks              = 0,
      peak_cutoff            = self.cutoff,
      interpolate            = True,
      min_distance_sym_equiv = 0,
      general_positions_only = False, # ???????XXXXXXXXXXXXXX
      min_cross_distance     = min_peak_peak_distance,
      min_cubicle_edge       = 5)
    psr = cgt.peak_search(parameters = peak_search_parameters,
      map = self.map_data).all(max_clusters = 99999999)
    # Convert peaks into water xray.structure
    mean_b = flex.mean(self.model.get_hierarchy().atoms().extract_b())
    sp = crystal.special_position_settings(self.model.crystal_symmetry())
    scatterers = flex.xray_scatterer()
    for site_frac in psr.sites():
      sc = xray.scatterer(
        label="O", site=site_frac, u=adptbx.b_as_u(mean_b), occupancy=1.0)
      scatterers.append(sc)
    self.xrs_water = xray.structure(sp, scatterers)
    #
    self.ma.add("  total peaks found: %d"%self.xrs_water.scatterers().size())
    self.ma.add("  B factors set to : %8.3f"%mean_b)
    if(self.debug): self._write_pdb_file(file_name_prefix="hoh_all_peaks")

  def _filter_by_distance(self):
    self.ma.add(
      "  distance (A), min: %4.2f max: %4.2f"%(self.dist_min, self.dist_max))
    self.ma.add("  start: %d"%self.xrs_water.scatterers().size())
    sel = mmtbx.utils.filter_water(
      interaction_selection  = self.interaction_selection,
      sites_frac_other       = self.model.get_sites_frac(),
      sites_frac_water       = self.xrs_water.sites_frac(),
      dist_min               = self.dist_min,
      dist_max               = self.dist_max,
      unit_cell              = self.unit_cell)
    self.xrs_water = self.xrs_water.select(sel)
    self.ma.add("  final: %d"%self.xrs_water.scatterers().size())
    if(self.debug): self._write_pdb_file(file_name_prefix="hoh_dist")

  def _filter_by_map_model_cc_and_maxwater(self):
    self._filter_by_map_model_cc(limit_number_of_water=True)

  def _filter_by_map_model_cc(self, limit_number_of_water=False):
    self.ma.add("  start: %d"%self.xrs_water.scatterers().size())
    # Compute model Fourier map
    m_all = self.model.deep_copy()
    m_all.add_solvent(
      solvent_xray_structure = self.xrs_water,
      atom_name              = "O",
      residue_name           = "HOH",
      chain_id               = "S",
      refine_adp             = "isotropic")
    m_all.setup_scattering_dictionaries(scattering_table=self.scattering_table)
    xrs_all = m_all.get_xray_structure()
    f_calc = xrs_all.structure_factors(d_min=self.resolution).f_calc()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xrs_all.unit_cell(),
      space_group_info      = xrs_all.space_group_info(),
      pre_determined_n_real = self.map_data.accessor().all())
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = f_calc)
    map_calc = fft_map.real_map_unpadded()
    sites_cart = xrs_all.sites_cart()

    ccs = flex.double()
    # Remove water by CC
    wsel = m_all.selection(string="water")
    for i, s in enumerate(wsel):
      if not s: continue
      cc = mmtbx.maps.correlation.from_map_map_atom(
        map_1     = self.map_data_resampled,
        map_2     = map_calc,
        site_cart = sites_cart[i],
        unit_cell = xrs_all.unit_cell(),
        radius    = self.atom_radius)
      ccs.append(cc)

    if(limit_number_of_water and
       self.xrs_water.scatterers().size() > self.n_max_water):
      sel = flex.sort_permutation(ccs, reverse=True)
      ccs_s = ccs.select(sel)
      cc_min = max(ccs_s[self.n_max_water], self.cc_mask_filter_threshold)
      fmt = "  updated cc_mask cutoff (to limit %d water): %6.4f (was: %6.4f)"
      self.ma.add(fmt%(self.n_max_water, cc_min, self.cc_mask_filter_threshold))
    else:
      cc_min = self.cc_mask_filter_threshold

    cntr=0
    for i, s in enumerate(wsel):
      if not s: continue
      cc = ccs[cntr]
      if cc < cc_min: wsel[i] = False
      cntr+=1
    self.xrs_water = m_all.select(wsel).get_xray_structure()

    # Exclude poor macro-molecule atoms from interaction analysis
    sites_cart = self.model.get_sites_cart()
    for i, s in enumerate(self.interaction_selection):
      if not s: continue # XXX
      cc = mmtbx.maps.correlation.from_map_map_atom(
        map_1     = self.map_data_resampled,
        map_2     = map_calc,
        site_cart = sites_cart[i],
        unit_cell = xrs_all.unit_cell(),
        radius    = self.atom_radius)
      if cc < self.cc_mask_threshold_interacting_atoms:
        self.interaction_selection[i] = False
    self.sites_frac_interaction = self.model.get_sites_frac().select(
      self.interaction_selection)
    #
    self.ma.add("  final: %d"%self.xrs_water.scatterers().size())

  def _filter_by_sphericity(self):
    sel = flex.bool(self.xrs_water.scatterers().size(), False)
    for i, site_frac in enumerate(self.xrs_water.sites_frac()):
      site_cart = self.unit_cell.orthogonalize(site_frac)
      o = maptbx.sphericity_by_heuristics(
        map_data    = self.map_data,
        unit_cell   = self.unit_cell,
        center_cart = site_cart,
        radius      = 1.0)
      mi,ma,me = o.rho
      #if(self.debug):
      #  _debug_show_all_plots(
      #    map_data    = self.map_data,
      #    unit_cell   = self.unit_cell,
      #    center_cart = site_cart,
      #    radius      = 1.0,
      #    plot_number = str("%d_%4.2f_%4.2f"%(i, o.ccs[0],o2.ccs[0])))
      if(o.ccs[0]>self.scc):
        sel[i] = True
        #print i, "%6.3f %6.3f %6.3f"%mv.min_max_mean().as_tuple(),\
        #  "%6.3f %6.3f %6.3f"%ccs.min_max_mean().as_tuple(),\
        #  "%6.3f %6.3f %6.3f"%ccs2.min_max_mean().as_tuple()
      #  print(hoh_str%(i,i,sc[0],sc[1],sc[2]), file=of)
      #else:
      #  print(hoh_str%(i,i,sc[0],sc[1],sc[2]), file=of2)
      #  #print i, "%6.3f %6.3f %6.3f"%mv.min_max_mean().as_tuple(),\
      #  #  "%6.3f %6.3f %6.3f"%ccs.min_max_mean().as_tuple(),\
      #  #  "%6.3f %6.3f %6.3f"%ccs2.min_max_mean().as_tuple(), "<<< REJECTED"
    self.xrs_water = self.xrs_water.select(sel)
    self.ma.add("  peaks left: %d"%self.xrs_water.scatterers().size())

  def _show_peak_profiles(self):
    for i, site_frac in enumerate(self.xrs_water.sites_frac()):
      site_cart = self.unit_cell.orthogonalize(site_frac)
      o = maptbx.sphericity_by_heuristics(
        map_data    = self.map_data,
        unit_cell   = self.unit_cell,
        center_cart = site_cart,
        radius      = 0.5)
      o2 = maptbx.sphericity_by_heuristics(
        map_data    = self.map_data,
        unit_cell   = self.unit_cell,
        center_cart = site_cart,
        radius      = 1.0)
      _debug_show_all_plots(
        map_data    = self.map_data,
        unit_cell   = self.unit_cell,
        center_cart = site_cart,
        radius      = 1.0,
        plot_number = str("%d_%4.2f_%4.2f"%(i+1, o.ccs[0],o2.ccs[0])))

  def _refine_water_adp(self):
    # Make water-only mmm from existing objects: an ugly set of manipulations
    m = self.model.deep_copy()
    m.add_solvent(
      solvent_xray_structure = self.xrs_water,
      atom_name              = "O",
      residue_name           = "HOH",
      chain_id               = "S",
      refine_adp             = "isotropic")
    m_w = m.select( m.selection(string="water") )
    self.ma.add("  B (min/max/mean) start: %8.3f %8.3f %8.3f"%
      m_w.get_b_iso().min_max_mean().as_tuple())
    m_w.setup_scattering_dictionaries(scattering_table = self.scattering_table)
    self.map_model_manager.set_model(m_w)
    # Isotropic ADP refinement: water only
    o = mmtbx.refinement.real_space.adp.\
      ncs_aware_refinement(
        map_model_manager = self.map_model_manager,
        d_min             = self.resolution,
        nproc             = 1,
        atom_radius       = self.atom_radius,
        individual        = True,
        restraints_weight = None,
        log               = null_out())
    self.xrs_water = self.map_model_manager.model().get_xray_structure()
    # Reset mmm back to original state
    self.map_model_manager.set_model(self.model)
    self.ma.add("  B (min/max/mean) final: %8.3f %8.3f %8.3f"%
      m_w.get_b_iso().min_max_mean().as_tuple())
    self.xrs_water.set_b_iso(values=m_w.get_b_iso())

  def _append_to_model(self):
    #
    chain_ids_taken = []
    for chain in self.model.get_hierarchy().chains():
      chain_ids_taken.append(chain.id)
    unique_taken = list(set(chain_ids_taken))
    if(len(unique_taken)==1):
      solvent_chain = unique_taken[0]
    else:
      for solvent_chain in all_chain_ids():
        if(not solvent_chain in chain_ids_taken):
          break
    self.ma.add("  new water chain ID: '%s'"%solvent_chain)
    bisos = self.xrs_water.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1)
    if bisos.size()>0:
      self.ma.add("  B (min/max/mean): %8.3f %8.3f %8.3f"%
        bisos.min_max_mean().as_tuple())
    #
    self.model.add_solvent(
      solvent_xray_structure = self.xrs_water,
      conformer_indices = None,
      atom_name    = "O",
      residue_name = "HOH",
      chain_id     = solvent_chain,
      refine_adp   = "isotropic")
