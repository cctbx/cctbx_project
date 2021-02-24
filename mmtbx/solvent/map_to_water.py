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

  def add(self, msg):
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


def run_split(map_model_manager, params, log):
  #
  mmm = map_model_manager
#  box_info = mmm.split_up_map_and_model_by_chain(
#    mask_around_unselected_atoms=False)
#  for box_mmm in box_info.mmm_list:
#    box_mmm = run(
#      map_model_manager = box_mmm,
#      dist_min          = params.dist_min,
#      dist_max          = params.dist_max,
#      step              = params.step,
#      mean_scale        = params.mean_scale,
#      debug             = params.debug,
#      log               = log)
#  mmm.merge_split_maps_and_models(
#    box_info = box_info,
#    replace_coordinates = True,
#    replace_u_aniso = False)
#  return mmm
  #
  o = run(
    map_model_manager = mmm,
    dist_min          = params.dist_min,
    dist_max          = params.dist_max,
    step              = params.step,
    mean_scale        = params.mean_scale,
    debug             = params.debug,
    log               = log)
  return o.map_model_manager

class run(object):
  def __init__(self,
               #model,
               #map_data,
               map_model_manager,
               dist_min=2.0,
               dist_max=3.2,
               step=0.3,
               mean_scale=0.7,
               debug=False,
               log=None):
    # Common parameters / variables
    self.model = map_model_manager.model()
    map_data = map_model_manager.map_manager().map_data()
    adopt_init_args(self, locals())
    self.total_time = 0
    if(map_data.accessor().origin()!=(0,0,0)):
      raise Sorry("Map origin must be zero.")
    self.ma = msg_accumulator(log = log)
    self.map_data = self.map_data.deep_copy() # will be changed in-place!
    self.unit_cell = self.model.crystal_symmetry().unit_cell()
    self.sites_frac = self.model.get_xray_structure().sites_frac()
    self.sites_frac_interaction = self.sites_frac.select(
      self.model.selection(selection_string_interaction))
    self.sites_frac_water = None
    self.cutoff = None
    self._wrote_file=0
    # Calculations
    self._call(self._massage_map         , "Normalize and re-sample input map")
    self._call(self._get_cutoff          , "Get cutoff")
    self._call(self._mask_out            , "Mask out molecule/solvent regions")
    self._call(self._find_peaks          , "Find peaks")
    self._call(self._filter_by_distance  , "Filter peaks by distance")
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
    file_name = file_name_prefix + "_%d.pdb"%self._wrote_file
    assert self.sites_frac_water is not None
    with open(file_name, "w") as fo:
      for i, site_frac in enumerate(self.sites_frac_water):
        site_cart = self.unit_cell.orthogonalize(site_frac)
        print(hoh_str%(i,i,site_cart[0],site_cart[1],site_cart[2]), file=fo)
    self._wrote_file += 1

  def _get_cutoff(self):
    map_at_atoms = flex.double()
    for site_frac in self.sites_frac:
      map_at_atoms.append( self.map_data.tricubic_interpolation(site_frac) )
    mean = flex.mean(map_at_atoms)
    sel = map_at_atoms > mean/3 # so it is not spoiled by too low/high values
    self.cutoff = flex.mean( map_at_atoms.select(sel) )*self.mean_scale
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
      sites_frac       = self.sites_frac_interaction,
      sites_frac_water = self.sites_frac_water,
      dist_min         = self.dist_min,
      dist_max         = self.dist_max,
      unit_cell        = self.unit_cell)
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
      mi,ma,me = o.rho
      if(self.debug):
        _debug_show_all_plots(
          map_data    = self.map_data,
          unit_cell   = self.unit_cell,
          center_cart = site_cart,
          radius      = 1.0,
          plot_number = str("%d_%4.2f_%4.2f"%(i, o.ccs[0],o2.ccs[0])))
      fl = (mi>0 and ma>0 and me>0 and o.ccs[0]>0.95 and o2.ccs[0]>0.90)
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
        radius      = 0.5,
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
