# LIBTBX_SET_DISPATCHER_NAME phenix.map_to_model_histogram

import sys, math
import iotbx.pdb
import iotbx.phil
import mmtbx.f_model
from iotbx import reflection_file_reader
from scitbx.array_family import flex
from cctbx import maptbx
import mmtbx.utils
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils

master_params_str = """\
d_min = None
  .type = float
bulk_solvent = False
  .type = bool
map_filter = 6.0
  .type = float
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
grid_step = 0.1
  .type = float
map_type = Fo-Fc
  .type = str
apply_sigma_scaling = False
  .type = bool
apply_volume_scaling = True
  .type = bool
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)


def show_histogram(data, n_slots):
  print flex.min(data), flex.max(data), flex.mean(data)
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print "%10.3f - %-10.3f : %10.2f" % (lc_1, hc_1, float(n_1)/(data.size())*100.)
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
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  if(f_obs_label is not None):
    parameters.labels = command_line.options.f_obs_label
  if(r_free_flags_label is not None):
    parameters.r_free_flags.label = command_line.options.r_free_flags_label
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = reflection_file_server,
    parameters              = parameters,
    keep_going              = True,
    log                     = log)
  f_obs = determine_data_and_flags_result.f_obs
  r_free_flags = determine_data_and_flags_result.r_free_flags
  return f_obs, r_free_flags


def get_stats(pdb_file_name,
              f_obs,
              r_free_flags,
              bulk_solvent,
              map_filter,
              map_type,
              grid_step,
              apply_volume_scaling,
              apply_sigma_scaling,
              d_min):
  xrs = iotbx.pdb.input(file_name = pdb_file_name).xray_structure_simple()
  #
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    xray_structure = xrs,
    r_free_flags   = r_free_flags)
  fmodel = fmodel.resolution_filter(d_min=d_min)
  fmodel.remove_outliers()
  #
  fmodel.info().show_rfactors_targets_scales_overall()
  #params = bss.master_params.extract()
  #params.bulk_solvent=bulk_solvent
  fmodel.update_solvent_and_scale()#(params = params)
  fmodel.info().show_rfactors_targets_scales_overall()
  fmodel.remove_outliers()
  if(not bulk_solvent):
    fmodel.update(k_sol=0, b_sol=0)
  fmodel.info().show_rfactors_targets_scales_overall()
  #
  print
  fft_map = fmodel.electron_density_map().fft_map(
    map_type          = map_type,
    resolution_factor = grid_step/d_min,
    acentrics_scale   = 1.,
    symmetry_flags    = maptbx.use_space_group_symmetry)
  #
  assert [apply_sigma_scaling, apply_volume_scaling].count(True) == 1
  if(apply_sigma_scaling):    fft_map.apply_sigma_scaling()
  elif(apply_volume_scaling): fft_map.apply_volume_scaling()
  else: assert RuntimeError
  nx,ny,nz = fft_map.n_real()
  print "n_real:", nx,ny,nz
  map_data = fft_map.real_map_unpadded()
  unit_cell = xrs.unit_cell()
  a,b,c = unit_cell.parameters()[:3]
  print "unit cell:", a,b,c
  #
  grid_sites_frac = flex.vec3_double()
  map_values = flex.double()
  for ix in xrange(nx):
    for iy in xrange(ny):
      for iz in xrange(nz):
        xf,yf,zf = ix/float(nx), iy/float(ny), iz/float(nz)
        grid_sites_frac.append([xf,yf,zf])
        map_at_ixiyiz = map_data[(ix,iy,iz)]
        map_values.append(map_at_ixiyiz)
  #
  xrsp1 = xrs.expand_to_p1()
  for sc in xrsp1.scatterers():
    site = sc.site
    x,y,z = site
    for t in xrange(3):
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
  #
  print
  show_histogram(data = map_values, n_slots=10)
  if(map_filter is not None):
    #for i in xrange(1000):
    #  map_values.set_selected(map_values == flex.max(map_values), 0)
    #  map_values.set_selected(map_values == flex.min(map_values), 0)
    mva = flex.abs(map_values)
    mmax =  flex.mean(mva)*map_filter
    mmin = -flex.mean(mva)*map_filter
    selp = map_values >= mmax
    seln = map_values <= mmin
    map_values.set_selected(selp, mmax)
    map_values.set_selected(seln, mmin)
    show_histogram(data = map_values, n_slots=10)
    print
  #
  print "ok up to here..."
  res = mmtbx.utils.density_distribution_per_atom(
    sites_frac_atoms = xrsp1.sites_frac(),
    sites_frac_peaks = grid_sites_frac,
    density_values   = map_values,
    unit_cell        = unit_cell)
  distances = res.distances()
  map_values = res.map_values()
  #
  print
  print "number of points", distances.size()
  #
  inc = 0.1
  n_points_max = -1
  for d in xrange(101):
    d = d/10.
    l = d
    r = d+inc
    sel  = distances >= l
    sel &= distances < r
    mv = map_values.select(sel)
    if(mv.size()>n_points_max): n_points_max = mv.size()
  #
  inc = 0.1
  for d in xrange(101):
    d = d/10.
    l = d
    r = d+inc
    sel  = distances >= l
    sel &= distances < r
    mv = map_values.select(sel)
    if(mv.size()>0):
      rms = math.sqrt( flex.sum(mv*mv)/mv.size() )
      print "%4.2f-%4.2f %8.4f %8.4f %10d %10.4f" % (l, r, flex.mean(mv), rms,
        mv.size(), mv.size()*1./n_points_max)

def run(args, log = None):
  if(log is None): log = sys.stdout
  if(len(args)==0):
    print "Usage:\n"
    print "phenix.map_to_model_histogram model.pdb data.mtz [parameters]"
  processed_args = mmtbx.utils.process_command_line_args(args = args, log = log,
    master_params = master_params())
  print >> log, "-"*79
  print >> log, "\nParameters:\n"
  processed_args.params.show(out = log, prefix=" ")
  if(len(args)==0): return
  params = processed_args.params.extract()
  #
  print >> log, "-"*79
  print >> log, "\nData:\n"
  f_obs, r_free_flags = get_f_obs_and_flags(
    reflection_file_name = processed_args.reflection_file_names[0],
    crystal_symmetry     = processed_args.crystal_symmetry,
    log                  = log)
  print >> log, "-"*79
  print >> log, "\nCalculations:\n"
  get_stats(pdb_file_name        = processed_args.pdb_file_names[0],
            f_obs                = f_obs,
            r_free_flags         = r_free_flags,
            bulk_solvent         = params.bulk_solvent,
            map_filter           = params.map_filter,
            grid_step            = params.grid_step,
            map_type             = params.map_type,
            apply_volume_scaling = params.apply_volume_scaling,
            apply_sigma_scaling  = params.apply_sigma_scaling,
            d_min                = params.d_min)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
