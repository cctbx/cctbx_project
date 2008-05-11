import sys, math
from cctbx.array_family import flex
from mmtbx import find_peaks
from mmtbx import utils
import iotbx.phil
from scitbx import matrix
from libtbx import adopt_init_args
from mmtbx.refinement import print_statistics
import mmtbx.utils
from cctbx import maptbx

master_params_part1 = iotbx.phil.parse("""\
map_type = mFobs-DFmodel
  .type = str
  .help = Map type to be used to find hydrogens
map_cutoff = 2.0
  .type = float
  .help = Map cutoff
angular_step = 3.0
  .type = float
  .help = Step in degrees for 6D rigid body search for best fit
""")

master_params_part2 = find_peaks.master_params.fetch(iotbx.phil.parse("""\
use_sigma_scaled_maps = True
resolution_factor = 1./4.
map_next_to_model
{
  min_model_peak_dist = 0.7
  max_model_peak_dist = 1.05
  min_peak_peak_dist = 1.0
  use_hydrogens = False
}
peak_search
{
  peak_search_level = 1
  min_cross_distance = 1.0
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
  fp_manager = find_peaks.manager(fmodel     = fmodel,
                                  map_type   = params.map_type,
                                  map_cutoff = params.map_cutoff,
                                  params     = params,
                                  log        = log)
  result = fp_manager.peaks_mapped()
  fp_manager.show_mapped(pdb_atoms = pdb_atoms)
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
    scatterers[water_and_peaks.i_seq_o ].site = o
    scatterers[water_and_peaks.i_seq_h1].site = h1
    scatterers[water_and_peaks.i_seq_h2].site = h2
  print >> log, "%6.3f"%d_best

def run(fmodel, model, log, params = None):
  print_statistics.make_header("Fit water hydrogens into residual map",
    out = log)
  if(params is None):
    params = all_master_params().extract()
  print_statistics.make_sub_header("find peak-candidates", out = log)
  peaks = find_hydrogen_peaks(
    fmodel = fmodel,
    pdb_atoms = model.pdb_atoms,
    params = params,
    log = log)
  waters_and_peaks = extract_hoh_peaks(
    peaks = peaks,
    pdb_hierarchy = model.pdb_hierarchy,
    pdb_atoms = model.pdb_atoms,
    xray_structure = model.xray_structure)
  print_statistics.make_sub_header("6D rigid body fit of HOH", out = log)
  print >> log, "Fit quality:"
  for water_and_peaks in waters_and_peaks:
    fit_water(water_and_peaks = water_and_peaks,
              xray_structure  = model.xray_structure,
              params          = params,
              log             = log)
  # adjust ADP for H
  u_isos = model.xray_structure.extract_u_iso_or_u_equiv()
  u_iso_mean = flex.mean(u_isos)
  sel_big = u_isos > u_iso_mean*2
  hd_sel = model.xray_structure.hd_selection()
  sel_big.set_selected(~hd_sel, False)
  model.xray_structure.set_u_iso(value = u_iso_mean, selection = sel_big)

def run_flip_hd(fmodel, model, log, params = None):
  if(params is None):
    params = all_master_params().extract()
  params.map_next_to_model.min_model_peak_dist = 0 # XXX
  print_statistics.make_sub_header("find peak-candidates", out = log)
  xray_structure_dc = fmodel.xray_structure.deep_copy_scatterers()
  xh_connectivity_table = model.xh_connectivity_table()
  scatterers = fmodel.xray_structure.scatterers()
  scattering_types = scatterers.extract_scattering_types()
  fmodel_dc = fmodel.deep_copy()
  for h_in_residue in model.refinement_flags.group_h:
    #print list(h_in_residue)
    x_atoms = []
    for i_seq_h in h_in_residue:
      for xh_connectivity_table_i in xh_connectivity_table:
        if(i_seq_h in xh_connectivity_table_i[:2]):
          x_atoms.append(xh_connectivity_table_i[0])
          break
    #print x_atoms
    xray_structure_dc_dc = xray_structure_dc.deep_copy_scatterers()
    xray_structure_dc_dc.set_occupancies(value = 0, selection = h_in_residue)
    fmodel_dc.update_xray_structure(
      xray_structure = xray_structure_dc_dc,
      update_f_calc  = True,
      update_f_mask  = True)
    #
    fft_map = fmodel_dc.electron_density_map(
      map_type          = "mFo-DFc",
      resolution_factor = 1/4.,
      symmetry_flags    = maptbx.use_space_group_symmetry)
    fft_map.apply_sigma_scaling()
    fft_map_data = fft_map.real_map_unpadded()
    for i_seq_h in h_in_residue:
      ed_val = maptbx.eight_point_interpolation(fft_map_data,
          scatterers[i_seq_h].site)
      assert fmodel_dc.xray_structure.scatterers()[i_seq_h].occupancy == 0
      if 1:#((ed_val < 0 and scattering_types[i_seq_h]=='D') or
         #(ed_val > 0 and scattering_types[i_seq_h]=='H')):
        print >> log, "%3d %6.2f %s" % (i_seq_h, ed_val, scattering_types[i_seq_h]), \
          model.pdb_atoms[i_seq_h].quote()
      #if(abs(ed_val) < 1.0):
      #  print >> log, "*** %3d %6.2f %s" % (i_seq_h, ed_val, scattering_types[i_seq_h]), \
      #  model.pdb_atoms[i_seq_h].quote()
      #  for map_cutoff in [3, -3]:
      #    params.map_cutoff = map_cutoff
      #    peaks = find_hydrogen_peaks(
      #      fmodel               = fmodel,
      #      pdb_atoms = model.pdb_atoms,
      #      params               = params,
      #      log                  = log)
