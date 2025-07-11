"""Compute averaged radial density distribution"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_model_histogram

import sys, math
import iotbx.pdb
import iotbx.phil
import mmtbx.f_model
from iotbx import reflection_file_reader
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import mmtbx.masks
from six.moves import zip
from six.moves import range
from iotbx import extract_xtal_data

master_params_str = """\
bulk_solvent_mode = *fast slow
  .type=choice(multi=False)
remove_outliers = True
  .type = bool
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
grid_step = 0.5
  .type = float
map_type = Fo-Fc
  .type = str
apply_sigma_scaling = False
  .type = bool
apply_volume_scaling = True
  .type = bool
neutron = False
  .type = bool
n_radial_shells = 1
  .type = int
include_f000 = True
  .type = bool
use_exact_phases = False
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)


def show_histogram(data, n_slots):
  print(flex.min(data), flex.max(data), flex.mean(data))
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print("%10.3f - %-10.3f : %10.2f" % (lc_1, hc_1, float(n_1)/(data.size())*100.))
    lc_1 = hc_1

def get_f_obs_and_flags(reflection_file_name,
                        crystal_symmetry,
                        f_obs_label = None,
                        r_free_flags_label = None,
                        log = None):
  reflection_files = []
  reflection_files.append(reflection_file_reader.any_reflection_file(
    file_name = reflection_file_name, ensure_read_access = False))
  reflection_file_server = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = reflection_files,
    err              = log)
  parameters = extract_xtal_data.data_and_flags_master_params().extract()
  if(f_obs_label is not None):
    parameters.labels = command_line.options.f_obs_label
  if(r_free_flags_label is not None):
    parameters.r_free_flags.label = command_line.options.r_free_flags_label
  determine_data_and_flags_result = extract_xtal_data.run(
    reflection_file_server = reflection_file_server,
    parameters             = parameters,
    keep_going             = True)
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  return f_obs, r_free_flags

def expand_to_p1(xrs):
  xrsp1 = xrs.expand_to_p1()
  for sc in xrsp1.scatterers():
    site = sc.site
    x,y,z = site
    for t in range(3):
      if x < 0: x = x + 1.
      if x > 1: x = x - 1.
      #
      if y < 0: y = y + 1.
      if y > 1: y = y - 1.
      #
      if z < 0: z = z + 1.
      if z > 1: z = z - 1.
    #
    sc.site = [x,y,z]
    #
    site = sc.site
    assert site[0] >= 0 and site[0] <= 1
    assert site[1] >= 0 and site[1] <= 1
    assert site[2] >= 0 and site[2] <= 1
  return xrsp1

def map_stat(distances, map_values):
  result = []
  #
  n_points_max = -1
  nn=20
  x = [[i/100,i/100+nn/100.] for i in range(0,800, nn)]
  for x_ in x:
    l,r = x_
    sel  = distances >= l
    sel &= distances < r
    mv = map_values.select(sel)
    if(mv.size()>n_points_max): n_points_max = mv.size()
  #
  for x_ in x:
    l,r = x_
    sel  = distances >= l
    sel &= distances < r
    mv = map_values.select(sel)
    if(mv.size()>0):
      sz = mv.size()
      rms = math.sqrt( flex.sum(mv*mv)/sz )
      #fr = sz*100./map_values.size()
      fr = sz*1./n_points_max
      result.append([l, r, flex.mean(mv), rms, sz, fr])
  return result

def get_map_values_and_grid_sites_frac(
      fmodel,
      map_type,
      grid_step,
      d_min,
      apply_sigma_scaling,
      apply_volume_scaling,
      include_f000,
      sel_bb,
      use_exact_phases):
  #
  resolution_factor = grid_step/d_min
  mp = mmtbx.masks.mask_master_params.extract()
  mp.grid_step_factor = 1./resolution_factor
  mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
    xray_structure = fmodel.xray_structure,
    d_min          = d_min,
    mask_params    = mp)
  bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
  sel = bulk_solvent_mask > 0
  bulk_solvent_mask = bulk_solvent_mask.set_selected(sel, 1)
  cr_gr = maptbx.crystal_gridding(
    unit_cell             = fmodel.xray_structure.unit_cell(),
    space_group_info      = fmodel.f_obs().space_group_info(),
    pre_determined_n_real = bulk_solvent_mask.focus())
  from mmtbx import map_tools
  from cctbx import miller
  #
  #mc = map_tools.electron_density_map(fmodel = fmodel).map_coefficients(
  #  map_type = map_type,
  #  acentrics_scale = 1.0,
  #  centrics_pre_scale = 1.0)
  if not use_exact_phases:
    k = fmodel.k_isotropic()*fmodel.k_anisotropic()
    print("flex.mean(k):", flex.mean(k))
    f_model = fmodel.f_model()
    mc_data = abs(fmodel.f_obs()).data()/k - abs(f_model).data()/k

    tmp = miller.array(miller_set = f_model,
      data = flex.double(f_model.indices().size(), 1)
      ).phase_transfer(phase_source = f_model)
    mc = miller.array(miller_set = tmp,
      data = mc_data * tmp.data())
  else:
    fmodel.update_all_scales(fast=True, remove_outliers=False)
    k = fmodel.k_isotropic()*fmodel.k_anisotropic()
    fo = fmodel.f_obs().customized_copy(data = fmodel.f_obs().data()/k)
    fo = fo.phase_transfer(phase_source = fmodel.f_model())
    fc = fmodel.f_calc().customized_copy(data = fmodel.f_calc().data())
    mc = miller.array(miller_set = fo,
      data = fo.data()-fc.data())




  ######## XXX
  fft_map = miller.fft_map(
    crystal_gridding     = cr_gr,
    fourier_coefficients = mc)
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()

  xrs = fmodel.xray_structure
  sites_cart = xrs.sites_cart().select(sel_bb)
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = xrs.unit_cell(),
    fft_n_real = map_data.focus(),
    fft_m_real = map_data.all(),
    sites_cart = sites_cart,
    site_radii = flex.double(sites_cart.size(), 0.5))
  map_in  = map_data.select(sel)
  mm = flex.mean(map_in)
  print("mean in (1):", mm)
  #
  #sites_frac = xrs.sites_frac().select(sel_bb)
  #mm = 0
  #for sf in sites_frac:
  #  mm += map_data.eight_point_interpolation(sf)
  #mm = mm/sites_frac.size()
  #print "mean in (2):", mm
  ########

  #
  # Add F000
  #reg = fmodel.xray_structure.scattering_type_registry(table = "wk1995")
  #f_000 = reg.sum_of_scattering_factors_at_diffraction_angle_0() +\
  #  0.4*fmodel.xray_structure.unit_cell().volume()
  if(include_f000):
    #f_000 = include_f000*fmodel.xray_structure.unit_cell().volume()*0.3
    #f_000 = None # XXX
    f_000 = abs(mm * xrs.unit_cell().volume())
    #f_000 = 0.626*fmodel.xray_structure.unit_cell().volume()*0.35
  else:
    f_000 = None
  print("f_000:", f_000)
  #print "XXX", include_f000*fmodel.xray_structure.unit_cell().volume()*0.3
  #
  fft_map = miller.fft_map(
    crystal_gridding     = cr_gr,
    fourier_coefficients = mc,
    f_000 = f_000)
  #
  assert [apply_sigma_scaling, apply_volume_scaling].count(True) == 1
  if(apply_sigma_scaling):    fft_map.apply_sigma_scaling()
  elif(apply_volume_scaling): fft_map.apply_volume_scaling()
  else: assert RuntimeError
  nx,ny,nz = fft_map.n_real()
  map_data = fft_map.real_map_unpadded()

  #map_data = map_data * bulk_solvent_mask
  print("n_real:", nx,ny,nz, map_data.size())
  grid_sites_frac = flex.vec3_double()
  map_values = flex.double()
  for ix in range(nx):
    for iy in range(ny):
      for iz in range(nz):
        mv = map_data[(ix,iy,iz)]
        if 1: #if(mv != 0):
          xf,yf,zf = ix/float(nx), iy/float(ny), iz/float(nz)
          grid_sites_frac.append([xf,yf,zf])
          map_at_ixiyiz = map_data[(ix,iy,iz)]
          map_values.append(map_at_ixiyiz)
  return map_values, grid_sites_frac

def show_fmodel(fmodel, prefix=""):
  print("%s r_work=%6.4f r_free=%6.4f d_min=%6.4f nref=%d"%(prefix,
    fmodel.r_work(),fmodel.r_free(),fmodel.f_obs().d_min(),
    fmodel.f_obs().data().size()))

def get_stats(pdb_file_name,
              f_obs,
              r_free_flags,
              map_type,
              grid_step,
              apply_volume_scaling,
              apply_sigma_scaling,
              mode,
              neutron,
              n_radial_shells,
              include_f000,
              remove_outliers,
              use_exact_phases):
  results = []
  for bulk_solvent in [True, False]:
    print("bulk_solvent:", bulk_solvent, "-"*50)
    pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
    pdb_hierarchy = pdb_inp.construct_hierarchy()
    sstring = """ pepnames and (name ca or name n or name c) and altloc " " """
    sel_bb = pdb_hierarchy.atom_selection_cache().selection(string = sstring)
    xrs = pdb_inp.xray_structure_simple()
    if(neutron):
      xrs.switch_to_neutron_scattering_dictionary()
    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.n_radial_shells = n_radial_shells
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      xray_structure = xrs,
      r_free_flags   = r_free_flags,
      mask_params    = mask_params)
    #
    show_fmodel(fmodel=fmodel, prefix="start:")
    if(mode=="fast"): fast=True
    elif(mode=="slow"): fast=False
    else: assert 0
    fmodel.update_all_scales(fast=fast, remove_outliers=remove_outliers)
    if(not bulk_solvent):
      k_mask = flex.double(fmodel.f_obs().data().size(), 0)
      fmodel.update(k_mask = k_mask)
    show_fmodel(fmodel=fmodel, prefix="final:")
    #
    map_values, grid_sites_frac = get_map_values_and_grid_sites_frac(
      fmodel               = fmodel,
      map_type             = map_type,
      grid_step            = grid_step,
      d_min                = f_obs.d_min(),
      apply_sigma_scaling  = apply_sigma_scaling,
      apply_volume_scaling = apply_volume_scaling,
      include_f000         = include_f000,
      sel_bb               = sel_bb,
      use_exact_phases     = use_exact_phases)
    #
    res = mmtbx.utils.density_distribution_per_atom(
      sites_frac_atoms = expand_to_p1(xrs = xrs).sites_frac(),
      sites_frac_peaks = grid_sites_frac,
      density_values   = map_values,
      unit_cell        = xrs.unit_cell())
    distances = res.distances()
    map_values = res.map_values()
    #
    result = map_stat(distances=distances, map_values=map_values)
    results.append(result)
  for result in zip(results[0],results[1]):
    l,r,p = result[0][0], result[0][1], result[0][5]
    m1,r1 = result[0][2], result[0][3]
    m2,r2 = result[1][2], result[1][3]
    print("%4.2f-%4.2f %10.4f | %8.4f %8.4f | %8.4f %8.4f" % (l, r, p, m1,r1, m2,r2))

def run(args, log = None):
  if(log is None): log = sys.stdout
  if(len(args)==0):
    print("Usage:\n")
    print("phenix.map_to_model_histogram model.pdb data.mtz [parameters]")
  processed_args = mmtbx.utils.process_command_line_args(args = args, log = log,
    master_params = master_params())
  print("-"*79, file=log)
  print("\nParameters:\n", file=log)
  processed_args.params.show(out = log, prefix=" ")
  if(len(args)==0): return
  params = processed_args.params.extract()
  #
  print("-"*79, file=log)
  print("\nData:\n", file=log)
  f_obs, r_free_flags = get_f_obs_and_flags(
    reflection_file_name = processed_args.reflection_file_names[0],
    crystal_symmetry     = processed_args.crystal_symmetry,
    log                  = log)
  print("-"*79, file=log)
  print("\nCalculations:\n", file=log)
  get_stats(pdb_file_name        = processed_args.pdb_file_names[0],
            f_obs                = f_obs,
            r_free_flags         = r_free_flags,
            grid_step            = params.grid_step,
            map_type             = params.map_type,
            apply_volume_scaling = params.apply_volume_scaling,
            apply_sigma_scaling  = params.apply_sigma_scaling,
            mode                 = params.bulk_solvent_mode,
            neutron              = params.neutron,
            n_radial_shells      = params.n_radial_shells,
            include_f000         = params.include_f000,
            remove_outliers      = params.remove_outliers,
            use_exact_phases     = params.use_exact_phases)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

