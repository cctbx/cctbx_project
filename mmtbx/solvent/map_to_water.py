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

elements=["N","O","S","Se","Fe","Mn","Mg","Mn","Zn","Ca","Cd","Cu","Co"]
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
  dic = {}
  for aw in wh.atoms():
    d_min    = 1.e9
    b_iso    = None
    chain_id = None
    for ap in nonwa:
      d = aw.distance(ap)
      if(d<d_min):
        d_min = d
        b_iso    = ap.b
        chain_id = ap.parent().parent().parent().id
    aw.set_b(b_iso)
    dic.setdefault(chain_id, []).append(aw)
  #
  pdb_model = nonw.only_model()
  for chain_id, atoms in zip(dic.keys(), dic.values()):
    new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
    for i_seq, new_atom in enumerate(atoms):
      new_atom_group = iotbx.pdb.hierarchy.atom_group(altloc="",
        resname="HOH")
      new_atom_group.append_atom(atom=new_atom.detached_copy())
      new_residue_group = iotbx.pdb.hierarchy.residue_group(
        resseq=iotbx.pdb.resseq_encode(value=i_seq), icode=" ")
      new_residue_group.append_atom_group(atom_group=new_atom_group)
      new_chain.append_residue_group(residue_group=new_residue_group)
    pdb_model.append_chain(chain=new_chain)
  #
  new_model = mmtbx.model.manager(
    model_input = None,
    pdb_hierarchy = nonw,
    crystal_symmetry = model.crystal_symmetry())
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
    #
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
      ma                = self.ma,
      map_model_manager = self.mmm,
      dist_min          = self.params.dist_min,
      dist_max          = self.params.dist_max,
      step              = self.params.step,
      map_threshold     = self.params.map_threshold,
      scc05             = self.params.scc05,
      scc1              = self.params.scc1,
      sphericity_filter = self.params.sphericity_filter,
      debug             = self.params.debug,
      total_time        = self.total_time,
      log               = self.log).map_model_manager
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
        ma                = self.ma,
        map_model_manager = box_mmm,
        dist_min          = self.params.dist_min,
        dist_max          = self.params.dist_max,
        step              = self.params.step,
        map_threshold     = self.params.map_threshold,
        sphericity_filter = self.params.sphericity_filter,
        scc05             = self.params.scc05,
        scc1              = self.params.scc1,
        debug             = self.params.debug,
        total_time        = self.total_time,
        log               = self.log).map_model_manager
      self.ma.set_prefix(prefix="")
    self.mmm.merge_split_maps_and_models(
      box_info                   = box_info,
      replace_coordinates        = True,
      allow_changes_in_hierarchy = True,
      replace_u_aniso            = False)
    model = models_as_chains(model = self.mmm.model())
    model = remove_clashing(model = model, dist_min = self.params.dist_min)
    self.mmm.set_model(model = model)

class run_one(object):
  def __init__(self,
               ma, # message_accumulator
               map_model_manager,
               dist_min=2.0,
               dist_max=3.2,
               step=0.3,
               map_threshold=1.5,
               sphericity_filter=True,
               scc05=0.97,
               scc1=0.9,
               debug=False,
               log=None,
               total_time=0,
               prefix=""):
    # Common parameters / variables
    self.model = map_model_manager.model()
    map_data = map_model_manager.map_manager().map_data()
    adopt_init_args(self, locals())
    self.total_time = total_time
    if(map_data.accessor().origin()!=(0,0,0)):
      raise Sorry("Map origin must be zero.")
    self.ma = ma
    self.map_data = self.map_data.deep_copy() # will be changed in-place!
    self.unit_cell = self.model.crystal_symmetry().unit_cell()
    self.sites_frac = self.model.get_xray_structure().sites_frac()
    self.interaction_selection = self.model.selection(
      selection_string_interaction)
    self.sites_frac_interaction = self.sites_frac.select(
      self.interaction_selection)
    self.sites_frac_water = None
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
    self._call(self._append_to_model     , "Add to model and reset internals")
    if(self.debug):
      self._call(self._show_peak_profiles, "Show peak profiles")

  def _call(self, func, msg):
    timer = user_plus_sys_time()
    self.ma.add(msg)
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
    if(max(sx, sy, sz) > self.step):
      n_real_fine = (int(a/self.step), int(b/self.step), int(c/self.step))
      self.ma.add("  re-sampled map dimenstions: %d %d %d"%n_real_fine)
      map_fine = flex.double(flex.grid(n_real_fine), 0)
      maptbx.resample(
        map_data     = self.map_data,
        map_data_new = map_fine,
        unit_cell    = self.unit_cell)
      self.map_data = map_fine
      self.ma.add("  input map (min,max,mean): %s"%self._map_mmm_str())
    if(self.debug):
      self._write_map(file_name = "resampled.mrc")

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
    for site_frac in self.sites_frac:
      map_at_atoms.append( self.map_data.tricubic_interpolation(site_frac) )
    mean = flex.mean(map_at_atoms)
    self.cutoff = self.map_data.as_1d().select(self.map_data.as_1d()>0
      ).standard_deviation_of_the_sample() * self.map_threshold
    self.ma.add("  mean density at atom centers: %8.3f"%mean)
    self.ma.add("  density cutoff for peak search: %8.3f"%self.cutoff)

  def _mask_out(self, rad_inside = 1.0, rad_outside = 6.0):
    # Zero inside
    radii = flex.double(self.sites_frac.size(), rad_inside)
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = self.sites_frac,
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
    self.sites_frac_water = psr.sites()
    self.ma.add("  total peaks found: %d"%self.sites_frac_water.size())
    if(self.debug): self._write_pdb_file(file_name_prefix="hoh_all_peaks")

  def _filter_by_distance(self):
    sel = mmtbx.utils.filter_water(
      sites_frac_interaction = self.sites_frac_interaction,
      sites_frac_other       = self.sites_frac.select(~self.interaction_selection),
      sites_frac_water       = self.sites_frac_water,
      dist_min               = self.dist_min,
      dist_max               = self.dist_max,
      unit_cell              = self.unit_cell)
    self.sites_frac_water = self.sites_frac_water.select(sel)
    self.ma.add(
      "  distance (A), min: %4.2f max: %4.2f"%(self.dist_min, self.dist_max))
    self.ma.add("  peaks left: %d"%self.sites_frac_water.size())
    if(self.debug): self._write_pdb_file(file_name_prefix="hoh_dist")

  def _filter_by_sphericity(self):
    tmp = flex.vec3_double()
    for i, site_frac in enumerate(self.sites_frac_water):
      site_cart = self.unit_cell.orthogonalize(site_frac)
      # XXX make it so it is called once! XXX
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
      mi,ma,me = o2.rho
      #if(self.debug):
      #  _debug_show_all_plots(
      #    map_data    = self.map_data,
      #    unit_cell   = self.unit_cell,
      #    center_cart = site_cart,
      #    radius      = 1.0,
      #    plot_number = str("%d_%4.2f_%4.2f"%(i, o.ccs[0],o2.ccs[0])))
      fl = (mi>0 and ma>0 and me>0 and o.ccs[0]>self.scc05 and o2.ccs[0]>self.scc1)
      if(fl):
        tmp.append(site_frac)
        #print i, "%6.3f %6.3f %6.3f"%mv.min_max_mean().as_tuple(),\
        #  "%6.3f %6.3f %6.3f"%ccs.min_max_mean().as_tuple(),\
        #  "%6.3f %6.3f %6.3f"%ccs2.min_max_mean().as_tuple()
      #  print(hoh_str%(i,i,sc[0],sc[1],sc[2]), file=of)
      #else:
      #  print(hoh_str%(i,i,sc[0],sc[1],sc[2]), file=of2)
      #  #print i, "%6.3f %6.3f %6.3f"%mv.min_max_mean().as_tuple(),\
      #  #  "%6.3f %6.3f %6.3f"%ccs.min_max_mean().as_tuple(),\
      #  #  "%6.3f %6.3f %6.3f"%ccs2.min_max_mean().as_tuple(), "<<< REJECTED"
    self.sites_frac_water = tmp
    self.ma.add("  peaks left: %d"%self.sites_frac_water.size())

  def _show_peak_profiles(self):
    for i, site_frac in enumerate(self.sites_frac_water):
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

  def _append_to_model(self):
    mean_b = flex.mean(self.model.get_hierarchy().atoms().extract_b())
    self.ma.add("  new water B factors will be set to mean B: %8.3f"%mean_b)
    sp = crystal.special_position_settings(self.model.crystal_symmetry())
    scatterers = flex.xray_scatterer()
    for site_frac in self.sites_frac_water:
      sc = xray.scatterer(
        label="O", site=site_frac, u=adptbx.b_as_u(mean_b), occupancy=1.0)
      scatterers.append(sc)
    xrs_water = xray.structure(sp, scatterers)
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
    self.ma.add("  new water will have chain ID: '%s'"%solvent_chain)
    #
    self.model.add_solvent(
      solvent_xray_structure = xrs_water,
      atom_name    = "O",
      residue_name = "HOH",
      chain_id     = solvent_chain,
      refine_adp   = "isotropic")
