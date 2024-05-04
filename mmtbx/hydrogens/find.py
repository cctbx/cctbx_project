from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx import find_peaks
from mmtbx import utils
import iotbx.phil
from scitbx import matrix
from libtbx import adopt_init_args
from mmtbx.refinement import print_statistics
import mmtbx.utils
from cctbx import sgtbx
import cctbx
from six.moves import zip,range

master_params_part1 = iotbx.phil.parse("""\
map_type = mFobs-DFmodel
  .type = str
  .help = Map type to be used to find hydrogens
map_cutoff = 2.0
  .type = float
  .help = Map cutoff
secondary_map_type = 2mFobs-DFmodel
  .type = str
  .help = Map type to be used to validate peaks in primary map
secondary_map_cutoff = 1.4
  .type = float
  .help = Secondary map cutoff
angular_step = 3.0
  .type = float
  .help = Step in degrees for 6D rigid body search for best fit
dod_and_od = False
  .type = bool
  .help = Build DOD/OD/O types of waters for neutron models
filter_dod = False
  .type = bool
  .help = Filter DOD/OD/O by correlation
min_od_dist = 0.60
  .type = float
  .help = Minimum O-D distance when building water from peaks
max_od_dist = 1.35
  .type = float
  .help = Maximum O-D distance when building water from peaks
min_dod_angle = 85.0
  .type = float
  .help = Minimum D-O-D angle when building water from peaks
max_dod_angle = 170.0
  .type = float
  .help = Maximum D-O-D angle when building water from peaks
h_bond_min_mac = 1.8
  .type = float
  .short_caption = H-bond minimum for DOD solvent-model
  .expert_level = 1
h_bond_max = 3.9
  .type = float
  .short_caption = Maximum H-bond length in DOD solvent model
  .expert_level = 1
""")

master_params_part2 = find_peaks.master_params.fetch(iotbx.phil.parse("""\
use_sigma_scaled_maps = True
resolution_factor = 1./4.
map_next_to_model
{
  min_model_peak_dist = 0.7
  max_model_peak_dist = 1.05
  min_peak_peak_dist = 0.7
  use_hydrogens = False
}
peak_search
{
  peak_search_level = 1
  min_cross_distance = 1.0
  general_positions_only=True
}
"""))

def all_master_params():
  return iotbx.phil.parse("""\
    include scope mmtbx.hydrogens.find.master_params_part1
    include scope mmtbx.hydrogens.find.master_params_part2
""", process_includes=True)

class h_peak(object):
  def __init__(self, site_frac,
                     height,
                     dist,
                     scatterer_o,
                     pdb_atom_o,
                     i_seq_o):
    self.site_frac = site_frac
    self.height = height
    self.dist = dist
    self.scatterer_o = scatterer_o
    self.pdb_atom_o = pdb_atom_o
    self.i_seq_o = i_seq_o

class water_and_peaks(object):
  def __init__(self, i_seq_o,
                     i_seq_h1,
                     i_seq_h2,
                     peaks_sites_frac):
    assert [i_seq_o,i_seq_h1,i_seq_h2,peaks_sites_frac].count(None) == 0
    adopt_init_args(self, locals())

def water_bond_angle(o,h1,h2):
  result = None
  a = h1[0]-o[0], h1[1]-o[1], h1[2]-o[2]
  b = h2[0]-o[0], h2[1]-o[1], h2[2]-o[2]
  a = matrix.col(a)
  b = matrix.col(b)
  return a.angle(b, deg=True)

def find_hydrogen_peaks(fmodel,
                        pdb_atoms,
                        params,
                        log):
  #
  from cctbx import maptbx
  e_map = fmodel.electron_density_map()
  crystal_symmetry = fmodel.xray_structure.crystal_symmetry()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = crystal_symmetry.unit_cell(),
    space_group_info = crystal_symmetry.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = 0.6)
  coeffs = e_map.map_coefficients(
    map_type     = params.map_type,
    fill_missing = False,
    isotropize   = True)
  fft_map = coeffs.fft_map(crystal_gridding = crystal_gridding)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  #
  fp_manager = find_peaks.manager(map_data       = map_data,
                                  xray_structure = fmodel.xray_structure,
                                  map_cutoff = params.map_cutoff,
                                  params     = params,
                                  log        = log)
  result = fp_manager.peaks_mapped()
  fp_manager.show_mapped(pdb_atoms = pdb_atoms)
  return result

def make_peak_dict(peaks, selection, obs_map, cutoff):
  result = {}
  for i in flex.sort_permutation(data=peaks.iseqs_of_closest_atoms):
    s = peaks.sites[i]
    h = peaks.heights[i]
    obsh = obs_map.eight_point_interpolation(s)
    if obsh<cutoff:
      continue
    i_seq = peaks.iseqs_of_closest_atoms[i]
    if(selection[i_seq]):
      if i_seq in result:
        result[i_seq].extend( [(h, s)] )
      else:
        result[i_seq] = [(h, s)]
  return result

def extract_hoh_peaks(peaks, pdb_hierarchy, pdb_atoms, xray_structure):
  scatterers = xray_structure.scatterers()
  assert scatterers.size() == pdb_atoms.size()
  assert peaks.sites.size() == peaks.heights.size()
  assert peaks.heights.size() == peaks.iseqs_of_closest_atoms.size()
  sentinel = pdb_atoms.reset_tmp(first_value=0, increment=0)
  for i_seq in peaks.iseqs_of_closest_atoms:
    pdb_atoms[i_seq].tmp = 1
  get_class = iotbx.pdb.common_residue_names_get_class
  o_i_seq_ag = {}
  result = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          is_water = (get_class(name=ag.resname) == "common_water")
          for atom in ag.atoms():
            if (atom.tmp == 0): continue
            assert atom.element.strip() not in ['H','D']
            if (is_water):
              assert atom.element.strip() == 'O'
              o_i_seq_ag[atom.i_seq] = ag
              break
  del sentinel
  o_i_seq_i_peak = {}
  for i_seq in o_i_seq_ag:
    o_i_seq_i_peak[i_seq] = []
  for i_peak,i_seq in enumerate(peaks.iseqs_of_closest_atoms):
    if(o_i_seq_i_peak.get(i_seq) is not None):
      o_i_seq_i_peak[i_seq].append(i_peak)
  for i_seq,ag in o_i_seq_ag.items():
    ag_atoms = ag.atoms()
    assert ag_atoms.size() == 3
    o_atom = None
    h_atoms = []
    for atom in ag_atoms:
      el = atom.element.strip()
      if (el == 'O'):
        o_atom = atom
      else:
        assert el in ['H','D']
        h_atoms.append(atom)
    assert o_atom is not None
    assert len(h_atoms) == 2
    result.append(water_and_peaks(
      i_seq_o = o_atom.i_seq,
      i_seq_h1 = h_atoms[0].i_seq,
      i_seq_h2 = h_atoms[1].i_seq,
      peaks_sites_frac = [peaks.sites[i_peak]
        for i_peak in o_i_seq_i_peak[i_seq]]))
  return result

def fit_water(water_and_peaks, xray_structure, params, log):
  scatterers = xray_structure.scatterers()
  uc = xray_structure.unit_cell()
  site_frac_o  = scatterers[water_and_peaks.i_seq_o ].site
  site_frac_h1 = scatterers[water_and_peaks.i_seq_h1].site
  site_frac_h2 = scatterers[water_and_peaks.i_seq_h2].site
  peak_sites_frac = water_and_peaks.peaks_sites_frac
  if(len(peak_sites_frac) == 1):
    sc1 = scatterers[water_and_peaks.i_seq_h1]
    sc2 = scatterers[water_and_peaks.i_seq_h2]
    if(sc1.occupancy < sc2.occupancy and sc2.occupancy > 0.3):
      site_frac_h2 = sc2.site
    elif(sc1.occupancy > sc2.occupancy and sc1.occupancy > 0.3):
      site_frac_h2 = scatterers[water_and_peaks.i_seq_h1].site
    else:
      site_frac_h2 = peak_sites_frac[0]
    result = mmtbx.utils.fit_hoh(
      site_frac_o     = site_frac_o,
      site_frac_h1    = site_frac_h1,
      site_frac_h2    = site_frac_h2,
      site_frac_peak1 = peak_sites_frac[0],
      site_frac_peak2 = site_frac_h2,
      angular_shift   = params.angular_step,
      unit_cell       = uc)
    d_best = result.dist_best()
    o = uc.fractionalize(result.site_cart_o_fitted)
    h1 = uc.fractionalize(result.site_cart_h1_fitted)
    h2 = uc.fractionalize(result.site_cart_h2_fitted)
  else:
    peak_pairs = []
    for i, s1 in enumerate(peak_sites_frac):
      for j, s2 in enumerate(peak_sites_frac):
        if i < j:
          peak_pairs.append([s1,s2])
    d_best = 999.
    for pair in peak_pairs:
      result = mmtbx.utils.fit_hoh(
        site_frac_o     = site_frac_o,
        site_frac_h1    = site_frac_h1,
        site_frac_h2    = site_frac_h2,
        site_frac_peak1 = pair[0],
        site_frac_peak2 = pair[1],
        angular_shift   = params.angular_step,
        unit_cell       = uc)
      if(result.dist_best() < d_best):
        d_best = result.dist_best()
        o = uc.fractionalize(result.site_cart_o_fitted)
        h1 = uc.fractionalize(result.site_cart_h1_fitted)
        h2 = uc.fractionalize(result.site_cart_h2_fitted)
  if(d_best < 1.0):
    # do not move HOH located on special position
    skip = False
    sites = [scatterers[water_and_peaks.i_seq_o ].site,
             scatterers[water_and_peaks.i_seq_h1].site,
             scatterers[water_and_peaks.i_seq_h2].site]
    for site in sites:
      site_symmetry = sgtbx.site_symmetry(
        xray_structure.unit_cell(),
        xray_structure.space_group(),
        site,
        0.5,
        True)
      if(site_symmetry.n_matrices() != 1):
        skip = True
        break
    #
    if(not skip):
      scatterers[water_and_peaks.i_seq_o ].site = o
      scatterers[water_and_peaks.i_seq_h1].site = h1
      scatterers[water_and_peaks.i_seq_h2].site = h2
  print("%6.3f"%d_best, file=log)

def run(fmodel, model, log, params = None):
  print_statistics.make_header("Fit water hydrogens into residual map",
    out = log)
  if(params is None):
    params = all_master_params().extract()
  print_statistics.make_sub_header("find peak-candidates", out = log)
  peaks = find_hydrogen_peaks(
    fmodel = fmodel,
    pdb_atoms = model.get_atoms(),
    params = params,
    log = log)
  waters_and_peaks = extract_hoh_peaks(
    peaks = peaks,
    pdb_hierarchy = model.get_hierarchy(),
    pdb_atoms = model.get_atoms(),
    xray_structure = model.get_xray_structure())
  print_statistics.make_sub_header("6D rigid body fit of HOH", out = log)
  print("Fit quality:", file=log)
  for water_and_peaks in waters_and_peaks:
    fit_water(water_and_peaks = water_and_peaks,
              xray_structure  = model.get_xray_structure(),
              params          = params,
              log             = log)
  # adjust ADP for H
  # TODO mrt: probably H bfactors should be equal to those
  # of the bonded atom
  u_isos = model.get_xray_structure().extract_u_iso_or_u_equiv()
  u_iso_mean = flex.mean(u_isos)
  sel_big = u_isos > u_iso_mean*2
  hd_sel = model.get_hd_selection()
  sel_big.set_selected(~hd_sel, False)
  model.get_xray_structure().set_u_iso(value = u_iso_mean, selection = sel_big)

def build_dod_and_od(model, fmodels, log=None, params=None):
  fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
    update_f_calc  = True,
    update_f_mask  = True)
  fmodels.show_short()
  hbparams = params.hydrogens.build
  if hbparams is None:
    hbparams = all_master_params().extract()
  if log is None:
    log = fmodels.log
  if fmodels.fmodel_n is not None:
    fmodel = fmodels.fmodel_neutron()
  else:
    fmodel = fmodels.fmodel_xray()
  title = "xray"
  if fmodel.xray_structure.guess_scattering_type_neutron():
    title="neutron"
  print_statistics.make_header("Build water hydrogens into "+title+" difference"
    " map", out=log)
  model = build_water_hydrogens_from_map2(model, fmodel, params=hbparams, log=log)
  #
  if hbparams.filter_dod:
    bmax=params.ordered_solvent.b_iso_max
    water_map_correlations(model, fmodels, log)
    model = remove_zero_occupancy(model,0.01,bmax)
    model.reprocess_pdb_hierarchy_inefficient()
  fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
    update_f_calc  = True,
    update_f_mask  = True)
  fmodels.show_short()
  mmtbx.utils.assert_model_is_consistent(model)
  sol_sel = model.solvent_selection()
  hd_sel = sol_sel & model.get_hd_selection()
  print("Final number of water hydrogens: ", hd_sel.count(True), file=log)
  return model, fmodels

def remove_zero_occupancy(model, min_occupancy=0.05, max_b_iso=80.):
  atoms = model.get_xray_structure().scatterers()
  occ = atoms.extract_occupancies()
  # print "dir(atoms): ", dir(atoms)
  umax = cctbx.adptbx.b_as_u(max_b_iso)
  uiso = atoms.extract_u_iso_or_u_equiv(model.get_xray_structure().unit_cell())
  sol_sel = model.solvent_selection()
  hd_sel = model.get_hd_selection()
  sel = sol_sel & hd_sel & ((occ < min_occupancy) | (uiso>umax))
  sel = ~sel
  return model.select( selection = sel)

# TODO code duplicate in model.add_hydrogens.insert_atoms
def insert_atom_into_model(xs, atom, atom_name, site_frac, occupancy, uiso, element):
  i_seq = atom.i_seq
  cart = xs.unit_cell().orthogonalize
  xyz = cart(site_frac)
  h = atom.detached_copy()
  h.name = atom_name
  h.xyz = xyz
  h.sigxyz = (0,0,0)
  h.occ = occupancy
  h.sigocc = 0
  h.b = cctbx.adptbx.u_as_b(uiso)
  h.sigb = 0
  h.uij = (-1,-1,-1,-1,-1,-1)
  if (iotbx.pdb.hierarchy.atom.has_siguij()):
    h.siguij = (-1,-1,-1,-1,-1,-1)
  h.element = "%2s" % element.strip()
  ag = atom.parent() # atom group
  rg = ag.parent()
  na =  ag.atoms_size()
  assert na==2 or na==1
  j_seq = i_seq+1
  if na==2:
    j_seq = j_seq +1
  #h.set_serial("") # 0) # j_seq+1)
  # h.i_seq = j_seq
  ag.append_atom(atom=h)
  appat = ag.atoms()[-1]
  scatterer = cctbx.xray.scatterer(
    label           = h.name,
    scattering_type = h.element.strip(),
    site            = site_frac,
    u               = uiso,
    occupancy       = h.occ)
  xs.add_scatterer(
    scatterer = scatterer,
    insert_at_index = j_seq)

def distances_to_peaks(xray_structure, sites_frac, peak_heights,
    distance_cutoff, use_selection=None):
  asu_mappings = xray_structure.asu_mappings(buffer_thickness =
    distance_cutoff)
  asu_mappings.process_sites_frac(sites_frac, min_distance_sym_equiv =
    xray_structure.min_distance_sym_equiv())
  pair_generator = cctbx.crystal.neighbors_fast_pair_generator(asu_mappings =
    asu_mappings, distance_cutoff = distance_cutoff)
  n_xray = xray_structure.scatterers().size()
  result = {}
  for pair in pair_generator:
    if(pair.i_seq < n_xray):
      if (pair.j_seq < n_xray): continue
      # i_seq = molecule
      # j_seq = site
      rt_mx_i = asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      i_seq_new_site_frac = pair.j_seq - n_xray
      new_site_frac = rt_mx_ji * sites_frac[i_seq_new_site_frac]
      jn = pair.i_seq
    else:
      if(pair.j_seq >= n_xray): continue
      # i_seq = site
      # j_seq = molecule
      rt_mx_i = asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = asu_mappings.get_rt_mx_j(pair)
      rt_mx_ij = rt_mx_j.inverse().multiply(rt_mx_i)
      i_seq_new_site_frac = pair.i_seq - n_xray
      new_site_frac = rt_mx_ij * sites_frac[i_seq_new_site_frac]
      jn = pair.j_seq
    if use_selection is not None:
      height = peak_heights[i_seq_new_site_frac]
      if(use_selection[jn]):
        if jn in result:
          result[jn].extend( [(height, new_site_frac)] )
        else:
          result[jn] = [(height, new_site_frac)]
  return result

def choose_h_for_water(unit_cell, o_site, xyz_h, h1_site=None, min_dod_angle=80,
    max_dod_angle=170):
  hh1_site = h1_site
  assert False
  if hh1_site is None:
    mx = max(xyz_h)
    xyz_h.remove(mx)
    hh1_site = mx[1]
  h_max = -1.e300
  h2_site = None
  d1 = unit_cell.distance(o_site, hh1_site)
  assert d1 > 0.5 and d1 < 1.5, "D-O distance is out of range: %f"%d1
  for h,s in xyz_h:
    d = unit_cell.distance(o_site, s)
    if d<0.5:
      continue
    a = unit_cell.angle(s, o_site, hh1_site)
    if a>min_dod_angle and a<max_dod_angle :
      if h>h_max:
        h_max = h
        h2_site = s
  result = []
  if h1_site is None:
    result = [hh1_site]
  if not h2_site is None :
    result.extend( [h2_site] )
  return result

def match_dod(unit_cell, xyz_h, min_od=0.5, max_od=1.35, min_dod_angle=85.,
    max_dod_angle=170.):
  n=len(xyz_h)
  besta=360.
  r = None
  for i in range(n):
    x1h=xyz_h[i]
    x1 = x1h[1]
    for j in range(i+1,n):
      x2h=xyz_h[j]
      x2 = x2h[1]
      d1 = unit_cell.distance(x1, x2)
      if d1<min_od:
        continue
      for k in range(j+1,n):
        x3h=xyz_h[k]
        x3 = x3h[1]
        d2 = unit_cell.distance(x2, x3)
        d3 = unit_cell.distance(x1, x3)
        if d2<min_od:
          continue
        if d3<min_od:
          continue
        a = [(x2,x1,x3), (x1,x2,x3), (x1,x3,x2) ]
        da = []
        for t in a:
          av = unit_cell.angle(t[0],t[1],t[2])
          td1 = unit_cell.distance(t[0],t[1])
          td2 = unit_cell.distance(t[1],t[2])
          if( td1>max_od or td2>max_od or av>max_dod_angle or av<min_dod_angle):
            da.append(360.)
          else:
            da.append(abs(av - 105.0))
        ibest = min(range(len(da)), key=da.__getitem__)
        if da[ibest] < besta and da[ibest]<90.:
          r = a[ibest]
          besta = da[ibest]
  return r

def match_od(unit_cell, xyz_h, min_od=0.5, max_od=1.35):
  n=len(xyz_h)
  besta=1.e10
  r = None
  for i in range(n):
    x1h=xyz_h[i]
    x1 = x1h[1]
    for j in range(i+1,n):
      x2h=xyz_h[j]
      x2 = x2h[1]
      d1 = unit_cell.distance(x1, x2)
      if d1<min_od or d1>max_od:
        continue
      dd1 =abs(d1-1.)
      if dd1 < besta:
        r = (x1,x2)
        besta = dd1
  return r

def len_peaks(set_of_peaks):
  r = 0
  for peaks in set_of_peaks.values():
    r += len(peaks)
  return r

def print_scats(xray_structure, fn):
  import cctbx
  ftmp = open(fn, "w")
  print(" scats ", file=ftmp)
  itmp=0
  for scat in xray_structure.scatterers():
    itmp=itmp+1
    print(itmp, ' ', scat.label.strip(), \
      ' ', scat.scattering_type.strip(), ' ', \
      xray_structure.unit_cell().orthogonalize(scat.site), \
      ' ', scat.occupancy, ' ', cctbx.adptbx.u_as_b(scat.u_iso), file=ftmp)
  ftmp.close()

def print_atom(out,atom):
  print('     ', atom.format_atom_record(), file=out)
  print("        atom.xyz:  ", atom.xyz, file=out)
  print("        atom.occ:  ", atom.occ, file=out)
  print("        atom.b:    ", atom.b, file=out)
  print('        atom.segid: "%s"' % atom.segid, file=out)
  print("        atom.i_seq: ", atom.i_seq, file=out)
  print("        atom.name: ", atom.name, file=out)
  print("        atom.element: ", atom.element.strip(), file=out)

def print_scat(out, s):
  print("scatterer  label : ", s.label.strip(), file=out)
  print("            site : ", s.site, file=out)
  print("         element : ", s.scattering_type.strip(), file=out)
  print("            biso : ", cctbx.adptbx.u_as_b(s.u_iso), file=out)
  print("             occ : ", s.occupancy, file=out)

def atom_scat(atom,scat):
  import StringIO
  s=StringIO.StringIO()
  print_atom(s,atom)
  print_scat(s,scat)
  # print>>s, "dir(s):\n",dir(scat)
  return s.getvalue()

def atom_as_str(atom):
  import StringIO
  s=StringIO.StringIO()
  print_atom(s,atom)
  return s.getvalue()

def assert_water_is_consistent(model):
  xs = model.get_xray_structure()
  unit_cell = xs.unit_cell()
  scatterers = xs.scatterers()
  hier = model.get_hierarchy()
  water_rgs = model.extract_water_residue_groups()
  for rg in water_rgs:
    if (rg.atom_groups_size() != 1):
      raise RuntimeError(
        "Not implemented: cannot handle water with alt. conf.")
    ag = rg.only_atom_group()
    atoms = ag.atoms()
    h_atoms = []
    o_atom=None
    if atoms.size()>0:
      for atom in atoms:
        if (atom.element.strip() == "O"):
          o_atom = atom
        else:
          h_atoms.append(atom)
    else:
      assert False
    o_i = o_atom.i_seq
    o_site = scatterers[o_i].site
    for hatom in h_atoms:
      hsite = scatterers[hatom.i_seq].site
      doh = unit_cell.distance(hsite, o_site)
      assert doh >0.35 and doh < 1.45, "%f\n%s"%(doh,atom_as_str(hatom))


def build_water_hydrogens_from_map(model, fmodel, params=None, log=None):
  self = model
  xs = self.xray_structure
  unit_cell = xs.unit_cell()
  scatterers = xs.scatterers()
  hier = model.get_hierarchy()
  mmtbx.utils.assert_model_is_consistent(model)
  if log is None:
    log = model.log
  if params is None:
    params = all_master_params().extract()
  max_od_dist = params.max_od_dist
  min_od_dist = params.min_od_dist
  assert min_od_dist >0.5
  assert max_od_dist >min_od_dist and max_od_dist<1.5
  keep_max = params.map_next_to_model.max_model_peak_dist
  keep_min = params.map_next_to_model.min_model_peak_dist
  params.map_next_to_model.max_model_peak_dist = max_od_dist * 1.8
  params.map_next_to_model.min_model_peak_dist = min_od_dist
  params.peak_search.min_cross_distance = 0.9 #params.min_od_dist # 0.5
  params.map_next_to_model.use_hydrogens = True
  params.map_next_to_model.min_peak_peak_dist = min_od_dist #0.8
  max_dod_angle = params.max_dod_angle
  min_dod_angle = params.min_dod_angle
  assert max_dod_angle<180 and min_dod_angle>30 and max_dod_angle>min_dod_angle
  peaks = find_hydrogen_peaks(
    fmodel = fmodel,
    pdb_atoms = model.pdb_atoms,
    params = params,
    log = log)
  if peaks is None:
    return model
  params.map_next_to_model.use_hydrogens = False
  hs = peaks.heights
  params.map_next_to_model.max_model_peak_dist = keep_max
  params.map_next_to_model.min_model_peak_dist = keep_min
  sol_O = model.solvent_selection().set_selected(
    model.get_hd_selection(), False)
  print("Number of solvent molecules: ", sol_O.count(True), file=log)
  sol_sel = model.solvent_selection()
  hd_sel = sol_sel & model.get_hd_selection()
  print("Number of water hydrogens: ", hd_sel.count(True), file=log)
  pks = distances_to_peaks(xs, peaks.sites, hs, max_od_dist, use_selection=sol_O)
  pkss = distances_to_peaks(xs, peaks.sites, hs, max_od_dist*1.7, use_selection=sol_O)
  print("Peaks to consider: ", len(pks.keys()), " : ", len(pkss.keys()), file=log)
  water_rgs = self.extract_water_residue_groups()
  water_rgs.reverse()
  element='D'
  next_to_i_seqs = []
  for rg in water_rgs:
    if (rg.atom_groups_size() != 1):
      raise RuntimeError(
        "Not implemented: cannot handle water with alt. conf.")
    ag = rg.only_atom_group()
    atoms = ag.atoms()
    h1=0
    h1_atom=None
    h1_site=None
    h1_i_seq=None
    if atoms.size()==2:
      o_atom = None
      h_atom = None
      for atom in atoms:
        if (atom.element.strip() == "O"): o_atom = atom
        else:                             h_atom = atom
      assert [o_atom, h_atom].count(None) == 0
      h1_i_seq = h_atom.i_seq
      h1_site=scatterers[h1_i_seq].site
      h1_atom=h_atom
      if h_atom.name.strip().endswith('1'):
        h1=1
      elif h_atom.name.strip().endswith('2'):
        h1=2
    elif atoms.size()==1:
      o_atom = atoms[0]
      h1_site = None
    elif atoms.size()==3:
      continue
    else:
      assert False
    o_i = o_atom.i_seq
    if o_i in pks:
      o_site = scatterers[o_i].site
      o_u = scatterers[o_i].u_iso_or_equiv(unit_cell)
      h_sites = pks[o_i]
      hh = choose_h_for_water(unit_cell, o_site, h_sites, h1_site=h1_site,
          min_dod_angle=min_dod_angle, max_dod_angle=max_dod_angle)
      nats = atoms.size() + len(hh)
      if nats==2:
        if len(hh)==0:
          h_swap = h1_site
        else:
          h_swap = hh[0]
        assert o_i in pkss
        h_sitess = pkss[o_i]
        hhs = choose_h_for_water(unit_cell, h_swap, h_sitess, h1_site=o_site,
            min_dod_angle=min_dod_angle, max_dod_angle=max_dod_angle)
        if len(hhs)>0:
          assert len(hhs)==1
          scatterers[o_i].site = h_swap
          o_atom.xyz = unit_cell.orthogonalize(h_swap)
          if h1_site is not None:
            scatterers[h1_i_seq].site = o_site
            h1_atom.xyz = unit_cell.orthogonalize(o_site)
          else:
            hhs.append(o_site)
            assert len(hhs)==2
          hh = hhs
          assert (atoms.size() + len(hh)) == 3
      for i,site_frac in enumerate(hh):
        assert (h1>0 and i<1) or (h1==0 and i<2)
        if h1==1: j=2
        elif h1==2: j=1
        else: j=i+1
        name = element+str(j)
        i_seq = o_atom.i_seq
        # this breaks atom sequence: o_atom.i_seq
        insert_atom_into_model(xs, atom=o_atom, atom_name=name, site_frac=site_frac,
          occupancy=1, uiso=o_u, element=element)
        next_to_i_seqs.append(i_seq)
  if( model.refinement_flags is not None and len(next_to_i_seqs)!=0):
    # TODO: adp_group=True according to params.dod_and_od_group_adp
    model.refinement_flags.add(
      next_to_i_seqs=next_to_i_seqs, # [i_seq], # ,
      sites_individual = True,
      s_occupancies    = False,
      adp_individual_iso=True)
  print("Number of H added:", len(next_to_i_seqs), file=log)
  # print "DEBUG! ", dir(model.refinement_flags) #.size()
  # print "DEBUG! ", dir(model.refinement_flags.sites_individual) #.size()
  #print "DEBUG! ", model.refinement_flags.sites_individual.size()
  model.reprocess_pdb_hierarchy_inefficient()
  np =  model.refinement_flags.sites_individual.size()
  assert np == model.get_number_of_atoms()
  assert model.refinement_flags.sites_individual.count(True) == np
  mmtbx.utils.assert_model_is_consistent(model)
  sol_sel = model.solvent_selection()
  hd_sel = sol_sel & model.get_hd_selection()
  assert hd_sel.count(True) >= len(next_to_i_seqs)
  return model

def select_one_water(water_residue_group, n_atoms):
  rg = water_residue_group
  if (rg.atom_groups_size() != 1):
    raise RuntimeError(
      "Not implemented: cannot handle water with alt. conf.")
  ag = rg.only_atom_group()
  atoms = ag.atoms()
  ws = []
  for atom in atoms:
    i = atom.i_seq
    ws.append(i)
  result = flex.bool()
  for j in range(0,n_atoms):
    if j in ws:
      result.append(True)
    else:
      result.append(False)
  nw = result.count(True)
  assert nw>0 and nw <= 3
  return result

def build_water_hydrogens_from_map2(model, fmodel, params=None, log=None):
  # self = model
  xs = model.get_xray_structure()
  unit_cell = xs.unit_cell()
  scatterers = xs.scatterers()
  hier = model.get_hierarchy()
  mmtbx.utils.assert_model_is_consistent(model)
  #assert_water_is_consistent(model)
  model.reprocess_pdb_hierarchy_inefficient()
  assert_water_is_consistent(model)
  if log is None:
    log = model.log
  if params is None:
    params = all_master_params().extract()
  max_od_dist = params.max_od_dist
  min_od_dist = params.min_od_dist
  assert min_od_dist >0.5
  assert max_od_dist >min_od_dist and max_od_dist<1.5
  keep_max = params.map_next_to_model.max_model_peak_dist
  keep_min = params.map_next_to_model.min_model_peak_dist
  params.map_next_to_model.max_model_peak_dist = max_od_dist * 1.8
  params.map_next_to_model.min_model_peak_dist = min_od_dist
  params.peak_search.min_cross_distance = params.min_od_dist # 0.5
  if params.peak_search.min_cross_distance < 0.7:
    params.peak_search.min_cross_distance = 0.7
  #if params.map_next_to_model.min_model_peak_dist < 0.7:
  #  params.map_next_to_model.min_model_peak_dist = 0.7
  params.map_next_to_model.use_hydrogens = True
  params.map_next_to_model.min_peak_peak_dist = min_od_dist
  #if params.map_next_to_model.min_peak_peak_dist <0.7:
  #  params.map_next_to_model.min_peak_peak_dist = 0.7
  max_dod_angle = params.max_dod_angle
  min_dod_angle = params.min_dod_angle
  assert max_dod_angle<180 and min_dod_angle>30 and max_dod_angle>min_dod_angle
  peaks = find_hydrogen_peaks(
    fmodel = fmodel,
    pdb_atoms = model.get_hierarchy().atoms(),
    params = params,
    log = log)
  if peaks is None:
    return model
  params.map_next_to_model.use_hydrogens = False
  hs = peaks.heights
  params.map_next_to_model.max_model_peak_dist = keep_max
  params.map_next_to_model.min_model_peak_dist = keep_min
  sol_O = model.solvent_selection().set_selected(
    model.get_hd_selection(), False)
  print("Number of solvent molecules: ", sol_O.count(True), file=log)
  sol_sel = model.solvent_selection()
  hd_sel = sol_sel & model.get_hd_selection()
  print("Number of water hydrogens: ", hd_sel.count(True), file=log)
  # pks = distances_to_peaks(xs, peaks.sites, hs, max_od_dist, use_selection=sol_O)
  obsmap = obs_map(fmodel, map_type=params.secondary_map_type)
  # filtered_peaks = map_peak_filter(peaks.sites, obsmap)
  #print>>log, "Number of filtered peaks: ", len(filtered_peaks)
  # pkss = distances_to_peaks(xs, peaks.sites, hs, max_od_dist*1.7, use_selection=sol_O)
  #pkss = distances_to_peaks(xs, filtered_peaks, hs, max_od_dist*1.7,
  #    use_selection=sol_O)
  cutoff2 = params.secondary_map_cutoff
  assert cutoff2 > 0. and cutoff2 < 100.
  pkss = make_peak_dict(peaks, sol_O, obsmap, cutoff2)
  npeaks=0
  for pk in pkss.values():
    npeaks = npeaks + len(pk)
  print("Peaks to consider: ", npeaks, file=log)
  water_rgs = model.extract_water_residue_groups()
  water_rgs.reverse()
  element='D'
  next_to_i_seqs = []
  for rg in water_rgs:
    if (rg.atom_groups_size() != 1):
      raise RuntimeError(
        "Not implemented: cannot handle water with alt. conf.")
    ag = rg.only_atom_group()
    atoms = ag.atoms()
    h_atoms = []
    o_atom=None
    if atoms.size()>0:
      for atom in atoms:
        if (atom.element.strip() == "O"):
          o_atom = atom
        else:
          h_atoms.append(atom)
    else:
      assert False
    o_i = o_atom.i_seq
    if o_i in pkss:
      o_site = scatterers[o_i].site
      site_symmetry = sgtbx.site_symmetry(xs.unit_cell(), xs.space_group(),
        o_site, 0.5, True)
      if(site_symmetry.n_matrices() != 1):
        special = True
        continue # TODO: handle this situation

      o_u = scatterers[o_i].u_iso_or_equiv(unit_cell)
      h_sites = pkss[o_i]
      # print_atom(sys.stdout, o_atom)
      val = obsmap.eight_point_interpolation(o_site)
      if val>cutoff2:
        h_sites.append((0., o_site))
      for hatom in h_atoms:
        hsite = scatterers[hatom.i_seq].site
        doh = unit_cell.distance(hsite, o_site)
        assert doh >0.5 and doh < 1.45, "%f\n%s"%(doh,atom_as_str(hatom))
        val = obsmap.eight_point_interpolation(hsite)
        if val>cutoff2:
          h_sites.append((0., hsite))
      dod = match_dod(unit_cell, h_sites, min_od=params.min_od_dist,
          max_od=params.max_od_dist, min_dod_angle=params.min_dod_angle,
          max_dod_angle=params.max_dod_angle)
      if dod is None:
        od= match_od(unit_cell, h_sites, min_od=params.min_od_dist,
            max_od=params.max_od_dist)
        if od is not None:
          scatterers[o_i].site = od[0]
          o_atom.xyz = unit_cell.orthogonalize(od[0])
          h=od[1]
          if len(h_atoms)>0:
            hatom = h_atoms.pop()
            hatom.xyz = unit_cell.orthogonalize(h)
            hatom.occ = 1
            hatom.b = cctbx.adptbx.u_as_b(o_u)
            hatom.name = "D1"
            h_i = hatom.i_seq
            scatterers[h_i].label = "D1"
            scatterers[h_i].site = h
            scatterers[h_i].occupancy = 1
            scatterers[h_i].u_iso = o_u
            if len(h_atoms)>0:
              # mark for deletion
              hatom = h_atoms.pop()
              hatom.name = "D2"
              hatom.occ = 0.
              h_i = hatom.i_seq
              scatterers[h_i].label="D2"
              scatterers[h_i].occupancy=0.
          else:
            insert_atom_into_model(xs, atom=o_atom, atom_name="D1",
              site_frac=h, occupancy=1, uiso=o_u, element='D')
            next_to_i_seqs.append(o_atom.i_seq)
      else:
        scatterers[o_i].site = dod[1]
        o_atom.xyz = unit_cell.orthogonalize(dod[1])
        ii=1
        for h in (dod[0],dod[2]):
          name = "D"+str(ii)
          ii=ii+1
          if len(h_atoms)>0:
            hatom = h_atoms.pop()
            hatom.xyz = unit_cell.orthogonalize(h)
            hatom.occ = 1
            hatom.b = cctbx.adptbx.u_as_b(o_u)
            hatom.name = name
            h_i = hatom.i_seq
            scatterers[h_i].label = name
            scatterers[h_i].site = h
            scatterers[h_i].occupancy = 1
            scatterers[h_i].u_iso = o_u
          else:
            insert_atom_into_model(xs, atom=o_atom, atom_name=name,
              site_frac=h, occupancy=1, uiso=o_u, element='D')
            next_to_i_seqs.append(o_atom.i_seq)
  if( model.refinement_flags is not None and len(next_to_i_seqs)!=0):
    # TODO: adp_group=True according to params.dod_and_od_group_adp
    model.refinement_flags.add(
      next_to_i_seqs=next_to_i_seqs, # [i_seq], # ,
      sites_individual = True,
      s_occupancies    = False,
      adp_individual_iso=True)
  print("Number of H added:", len(next_to_i_seqs), file=log)
  # print "DEBUG! ", dir(model.refinement_flags) #.size()
  # print "DEBUG! ", dir(model.refinement_flags.sites_individual) #.size()
  #print "DEBUG! ", model.refinement_flags.sites_individual.size()
  model.reprocess_pdb_hierarchy_inefficient()
  if model.refinement_flags.sites_individual is not None:
    np =  model.refinement_flags.sites_individual.size()
    assert np == model.get_number_of_atoms()
    assert model.refinement_flags.sites_individual.count(True) == np
  mmtbx.utils.assert_model_is_consistent(model)
  sol_sel = model.solvent_selection()
  hd_sel = sol_sel & model.get_hd_selection()
  assert hd_sel.count(True) >= len(next_to_i_seqs)
  assert_water_is_consistent(model)
  if False:
    model.idealize_h_minimization()
    model.get_hierarchy()
    mmtbx.utils.assert_model_is_consistent(model)
    assert_water_is_consistent(model)
  return model


def get_pdb_oxygen(water_residue_group):
  rg = water_residue_group
  if (rg.atom_groups_size() != 1):
    raise RuntimeError(
      "Not implemented: cannot handle water with alt. conf.")
  ag = rg.only_atom_group()
  atoms = ag.atoms()
  for atom in atoms:
    if atom.element.strip() == 'O':
      return atom
  assert False

def one_water_correlation(model, fmodels, water):
  import mmtbx.solvent.ordered_solvent as ordered_solvent
  from mmtbx import real_space_correlation
  params = ordered_solvent.master_params().extract()
  par = params.secondary_map_and_map_cc_filter
  rcparams = real_space_correlation.master_params().extract()
  rcparams.detail = "residue"
  if fmodels.fmodel_n is not None:
    fmodel = fmodels.fmodel_neutron()
  else:
    fmodel = fmodels.fmodel_xray()
  title = "xray"
  if fmodel.xray_structure.guess_scattering_type_neutron():
    title="neutron"
  scatterers = model.get_xray_structure().scatterers()
  assert scatterers is fmodel.xray_structure.scatterers()
  results = real_space_correlation.map_statistics_for_atom_selection(
    atom_selection = water,
    fmodel = fmodel,
    map1_type="Fo",
    map2_type="Fmodel"
  )
  return results.cc

def scatterers_info(scatterers, selection):
  r = str()
  for s,u in zip(scatterers,selection):
    if u:
      r += s.element_symbol().strip()
      r += " : "
      r += ("%4.2f %6.3f"%(s.occupancy,s.u_iso))
      r += " = "
      r += s.label
      r += ";   "
  return r

def obs_map(
    fmodel,
    map_type,
    resolution_factor=0.25):
  map_coeffs = fmodel.electron_density_map().map_coefficients(map_type)
  return map_coeffs.fft_map(
    resolution_factor=resolution_factor).apply_sigma_scaling().real_map()


def map_peak_filter(sites_frac, obs_map, cutoff):
  result = flex.vec3_double()
  for site_frac in sites_frac:
    val = obs_map.eight_point_interpolation(site_frac)
    if val>cutoff:
      result.append(site_frac)
  return result


def water_map_correlations(model, fmodels, log=None):
  print_statistics.make_header("Water real space correlations", out=log)
  fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
    update_f_calc=True, update_f_mask=True)
  fmodels.show_short()
  scatterers = model.get_xray_structure().scatterers()
  waters = model.solvent_selection()
  water_rgs = model.extract_water_residue_groups()
  n_atoms = len(scatterers)
  for rg in water_rgs:
    o_atom = get_pdb_oxygen(rg)
    o_scat = scatterers[o_atom.i_seq]
    water = select_one_water(rg, n_atoms)
    if water.count(True) < 2:
      continue
    for s,u in zip(scatterers,water):
      if u:
        e = s.element_symbol().strip()
        if e=='D':
          if s.occupancy <= 0.02 or s.occupancy >=0.98:
            keep_occ = s.occupancy
            neutron_cc = one_water_correlation(model, fmodels, water)
            if s.occupancy <= 0.02:
              s.occupancy = 1.0
              s.u_iso = o_scat.u_iso
            else:
              s.occupancy = 0.0
            fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
              update_f_calc=True, update_f_mask=True)
            ncc = one_water_correlation(model, fmodels, water)
            if ncc < neutron_cc:
              s.occupancy = keep_occ
              fmodels.update_xray_structure(xray_structure = model.get_xray_structure(),
                update_f_calc=True, update_f_mask=True)
            else:
              neutron_cc = ncc
  return True
