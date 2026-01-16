"""Compute map correlation coefficient given input model and reflection data
"""

from __future__ import absolute_import, division, print_function
import sys
from six.moves import range
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx.file_reader import any_file
import iotbx.phil
from cctbx.array_family import flex
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import Sorry, null_out
from libtbx import group_args
from mmtbx.command_line.map_comparison import get_mtz_labels, get_d_min,\
  get_crystal_symmetry
from cctbx.sgtbx import space_group_info
from iotbx import extract_xtal_data

core_params_str = """\
atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation. Determined automatically if \
          if None is given
  .expert_level = 2
hydrogen_atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation for H or D.
  .expert_level = 2
resolution_factor = 1./4
  .type = float
use_hydrogens = None
  .type = bool
"""

map_files_params_str = """\
map_file_name = None
  .type = path
  .help = A CCP4-formatted map
  .style = file_type:ccp4_map input_file
d_min = None
  .type = float
  .short_caption = Resolution
  .help = Resolution of map
map_coefficients_file_name = None
  .type = path
  .help = MTZ file containing map
  .style = file_type:ccp4_map input_file process_hkl child:map_labels:map_coefficients_label
map_coefficients_label = None
  .type = str
  .short_caption = Data label
  .help = Data label for complex map coefficients in MTZ file
  .style = renderer:draw_map_arrays_widget
"""

master_params_str = """\
%s
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
detail = atom residue *automatic
  .type = choice(multi=False)
  .help = Level of details to show CC for
map_1
  .help = First map to use in map CC calculation
{
 type = Fc
   .type = str
   .help = Electron density map type. Example xmFobs-yDFcalc (for \
           maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
           unweighted map), x and y are any real numbers.
 fill_missing_reflections = False
   .type = bool
 isotropize = False
   .type = bool
}
map_2
  .help = Second map to use in map CC calculation
{
 type = 2mFo-DFc
   .type = str
   .help = Electron density map type. Example xmFobs-yDFcalc (for \
           maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
           unweighted map), x and y are any real numbers.
 fill_missing_reflections = True
   .type = bool
 isotropize = True
   .type = bool
}
pdb_file_name = None
  .type = str
  .help = PDB/mmCIF file name.
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
high_resolution = None
  .type=float
low_resolution = None
  .type=float
%s
"""%(core_params_str, map_files_params_str)

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def extract_data_and_flags(params, crystal_symmetry=None):
  data_and_flags = None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name)
    reflection_file_server = reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry   = True,
      reflection_files = [reflection_file])
    parameters = extract_xtal_data.data_and_flags_master_params().extract()
    parameters.force_anomalous_flag_to_be_equal_to = False
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    if(params.high_resolution is not None):
      parameters.high_resolution = params.high_resolution
    if(params.low_resolution is not None):
      parameters.low_resolution = params.low_resolution
    data_and_flags = extract_xtal_data.run(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      extract_r_free_flags   = False) #XXX
  return data_and_flags

def compute_map_from_model(high_resolution, low_resolution, xray_structure,
                           grid_resolution_factor=None,
                           crystal_gridding = None):
  f_calc = xray_structure.structure_factors(d_min = high_resolution).f_calc()
  if (low_resolution is not None):
    f_calc = f_calc.resolution_filter(d_max = low_resolution)
  if (crystal_gridding is None):
    return f_calc.fft_map(
      resolution_factor = min(0.5,grid_resolution_factor),
      symmetry_flags    = None)
  return miller.fft_map(
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_calc)

def check_map_file(map_file, params):
  if ( (params.map_coefficients_file_name is not None) and
       (params.map_file_name is not None) ):
    raise Sorry('Please use map coefficients or the map, not both.')
  elif (params.map_coefficients_file_name is not None):
    file_handle = any_file(params.map_coefficients_file_name)
    if (file_handle.file_type != 'hkl'):
      raise Sorry('%s is not in MTZ format.' % file_handle.file_name)
    labels = get_mtz_labels(file_handle)
    if (params.map_coefficients_label is None):
      raise Sorry('Data labels for map coefficients are not specified for %s' %
                  params.map_coefficients_file_name)
    elif (params.map_coefficients_label not in labels):
      raise Sorry('%s labels for map_coefficients were not found in %s' %
                  (params.map_coefficients_label,
                   params.map_coefficients_file_name))
  elif (params.map_file_name is not None):
    file_handle = any_file(params.map_file_name)
    if (file_handle.file_type != 'ccp4_map'):
      raise Sorry('%s is not a CCP4-formatted map file.' %
                  file_handle.file_name)
  elif (map_file is not None):
    params.map_file_name = map_file.file_name

def broadcast(m, log):
  print("-"*79, file=log)
  print(m, file=log)
  print("*"*len(m), file=log)

def simple(fmodel, pdb_hierarchy, params=None, log=None, show_results=False):
  if(params is None): params = master_params().extract()
  if(log is None): log = sys.stdout

  crystal_gridding = None
  unit_cell = None
  d_min = 1.0
  map_1 = None
  map_2 = None

  # compute map_1 and map_2 if given F_obs (fmodel exists)
  if ( (params.map_file_name is None) and
       (params.map_coefficients_file_name is None) and
       (fmodel is not None) ):
    e_map_obj = fmodel.electron_density_map()
    coeffs_2 = e_map_obj.map_coefficients(
      map_type     = params.map_2.type,
      fill_missing = params.map_2.fill_missing_reflections,
      isotropize   = params.map_2.isotropize)
    fft_map_2 = coeffs_2.fft_map(resolution_factor = params.resolution_factor)
    crystal_gridding = fft_map_2
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()

    coeffs_1 = e_map_obj.map_coefficients(
      map_type     = params.map_1.type,
      fill_missing = params.map_1.fill_missing_reflections,
      isotropize   = params.map_1.isotropize)
    fft_map_1 = miller.fft_map(crystal_gridding = crystal_gridding,
                               fourier_coefficients = coeffs_1)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()

    unit_cell = fmodel.xray_structure.unit_cell()
    d_min = fmodel.f_obs().d_min()

  # or read map coefficents
  elif (params.map_coefficients_file_name is not None):
    map_handle = any_file(params.map_coefficients_file_name)
    crystal_symmetry = get_crystal_symmetry(map_handle)
    unit_cell = crystal_symmetry.unit_cell()
    d_min = get_d_min(map_handle)
    crystal_gridding = maptbx.crystal_gridding(
      crystal_symmetry.unit_cell(), d_min=d_min,
      resolution_factor=params.resolution_factor,
      space_group_info=crystal_symmetry.space_group_info())
    coeffs_2 = map_handle.file_server.get_miller_array(
      params.map_coefficients_label)
    fft_map_2 = miller.fft_map(crystal_gridding=crystal_gridding,
                               fourier_coefficients=coeffs_2)
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()

  # or read CCP4 map
  else:
    map_handle = any_file(params.map_file_name)
    unit_cell = map_handle.file_object.unit_cell()
    sg_info = space_group_info(map_handle.file_object.space_group_number)
    n_real = map_handle.file_object.unit_cell_grid
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell, space_group_info=sg_info, pre_determined_n_real=n_real)
    map_2 = map_handle.file_object.map_data()

    # check for origin shift
    # modified from phenix.command_line.real_space_refine
    # plan to centralize functionality in another location
    # -------------------------------------------------------------------------
    shift_manager = mmtbx.utils.shift_origin(
      map_data=map_2, pdb_hierarchy=pdb_hierarchy,
      crystal_symmetry=map_handle.crystal_symmetry())
    if (shift_manager.shift_cart is not None):
      print("Map origin is not at (0,0,0): shifting the map and model.", file=log)
    pdb_hierarchy = shift_manager.pdb_hierarchy
    map_2 = shift_manager.map_data
    # -------------------------------------------------------------------------

  # compute map_1 (Fc) if given a map (fmodel does not exist)
  if (map_1 is None):
    xray_structure = pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=crystal_gridding.crystal_symmetry())
    fft_map_1 = compute_map_from_model(d_min, None, xray_structure,
                                       crystal_gridding=crystal_gridding)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()

  # compute cc
  assert ( (map_1 is not None) and (map_2 is not None) )
  broadcast(m="Map correlation and map values", log=log)
  overall_cc = flex.linear_correlation(x = map_1.as_1d(),
    y = map_2.as_1d()).coefficient()
  print("  Overall map cc(%s,%s): %6.4f"%(params.map_1.type,
    params.map_2.type, overall_cc), file=log)
  detail, atom_radius = params.detail, params.atom_radius
  detail, atom_radius = set_detail_level_and_radius(
    detail=detail, atom_radius=atom_radius, d_min=d_min)
  use_hydrogens = params.use_hydrogens
  if(use_hydrogens is None):
    if(params.scattering_table == "neutron" or d_min <= 1.2):
      use_hydrogens = True
    else:
      use_hydrogens = False
  hydrogen_atom_radius = params.hydrogen_atom_radius
  if(hydrogen_atom_radius is None):
    if(params.scattering_table == "neutron"):
      hydrogen_atom_radius = atom_radius
    else:
      hydrogen_atom_radius = 1
  results = compute(
    pdb_hierarchy        = pdb_hierarchy,
    unit_cell            = unit_cell,
    fft_n_real           = map_1.focus(),
    fft_m_real           = map_1.all(),
    map_1                = map_1,
    map_2                = map_2,
    detail               = detail,
    atom_radius          = atom_radius,
    use_hydrogens        = use_hydrogens,
    hydrogen_atom_radius = hydrogen_atom_radius)
  if(show_results):
    show(log=log, results=results, params=params, detail=detail)
  return overall_cc, results

def show(log, results, detail, params=None, map_1_name=None, map_2_name=None):
  assert params is not None or [map_1_name,map_2_name].count(None)==0
  if([map_1_name,map_2_name].count(None)==2):
    map_1_name,map_2_name = params.map_1.type, params.map_2.type
  print(file=log)
  print("Rho1 = %s, Rho2 = %s"%(map_1_name, map_2_name), file=log)
  print(file=log)
  if(detail == "atom"):
    print(" <----id string---->  occ     ADP      CC   Rho1   Rho2", file=log)
  else:
    print("  <id string>    occ     ADP      CC   Rho1   Rho2", file=log)
  fmt = "%s %4.2f %7.2f %7.4f %6.2f %6.2f"
  for r in results:
    print(fmt%(r.id_str, r.occupancy, r.b, r.cc, r.map_value_1,
      r.map_value_2), file=log)

def compute(pdb_hierarchy,
            unit_cell,
            fft_n_real,
            fft_m_real,
            map_1,
            map_2,
            detail,
            atom_radius,
            use_hydrogens,
            hydrogen_atom_radius):
  assert detail in ["atom", "residue"]
  results = []
  for chain in pdb_hierarchy.chains():
    for residue_group in chain.residue_groups():
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          r_id_str = "%2s %1s %3s %4s %1s"%(chain.id, conformer.altloc,
            residue.resname, residue.resseq, residue.icode)
          r_sites_cart = flex.vec3_double()
          r_b          = flex.double()
          r_occ        = flex.double()
          r_mv1        = flex.double()
          r_mv2        = flex.double()
          r_rad        = flex.double()
          for atom in residue.atoms():
            a_id_str = "%s %4s"%(r_id_str, atom.name)
            if(atom.element_is_hydrogen()): rad = hydrogen_atom_radius
            else: rad = atom_radius
            if(not (atom.element_is_hydrogen() and not use_hydrogens)):
              map_value_1 = map_1.eight_point_interpolation(
                unit_cell.fractionalize(atom.xyz))
              map_value_2 = map_2.eight_point_interpolation(
                unit_cell.fractionalize(atom.xyz))
              r_sites_cart.append(atom.xyz)
              r_b         .append(atom.b)
              r_occ       .append(atom.occ)
              r_mv1       .append(map_value_1)
              r_mv2       .append(map_value_2)
              r_rad       .append(rad)
              if(detail == "atom"):
                sel = maptbx.grid_indices_around_sites(
                  unit_cell  = unit_cell,
                  fft_n_real = fft_n_real,
                  fft_m_real = fft_m_real,
                  sites_cart = flex.vec3_double([atom.xyz]),
                  site_radii = flex.double([rad]))
                cc = flex.linear_correlation(x=map_1.select(sel),
                  y=map_2.select(sel)).coefficient()
                result = group_args(
                  chain_id    = chain.id,
                  atom        = atom,
                  id_str      = a_id_str,
                  cc          = cc,
                  map_value_1 = map_value_1,
                  map_value_2 = map_value_2,
                  b           = atom.b,
                  occupancy   = atom.occ,
                  n_atoms     = 1)
                results.append(result)
          if(detail == "residue") and (len(r_mv1) > 0):
            sel = maptbx.grid_indices_around_sites(
              unit_cell  = unit_cell,
              fft_n_real = fft_n_real,
              fft_m_real = fft_m_real,
              sites_cart = r_sites_cart,
              site_radii = r_rad)
            cc = flex.linear_correlation(x=map_1.select(sel),
              y=map_2.select(sel)).coefficient()
            result = group_args(
              residue     = residue,
              chain_id    = chain.id,
              id_str      = r_id_str,
              cc          = cc,
              map_value_1 = flex.mean(r_mv1),
              map_value_2 = flex.mean(r_mv2),
              b           = flex.mean(r_b),
              occupancy   = flex.mean(r_occ),
              n_atoms     = r_sites_cart.size())
            results.append(result)
  return results

def set_detail_level_and_radius(detail, atom_radius, d_min):
  assert detail in ["atom","residue","automatic"]
  if(detail == "automatic"):
    if(d_min < 2.0): detail = "atom"
    else:            detail = "residue"
  if(atom_radius is None):
    if(d_min < 1.0):                    atom_radius = 1.0
    elif(d_min >= 1.0 and d_min<2.0):   atom_radius = 1.5
    elif(d_min >= 2.0 and d_min < 4.0): atom_radius = 2.0
    else:                               atom_radius = 2.5
  return detail, atom_radius

class selection_map_statistics_manager(object):
  """
  Utility class for performing repeated calculations on multiple maps.  Useful
  in post-refinement validation, ligand fitting, etc. where we want to collect
  both CC and values for 2mFo-DFc and mFo-DFc maps.
  """
  __slots__ = ["fft_m_real", "fft_n_real", "atom_selection", "sites",
               "sites_frac", "atom_radii", "map_sel"]

  def __init__(self,
      atom_selection,
      xray_structure,
      fft_m_real,
      fft_n_real,
      atom_radius=1.5,
      exclude_hydrogens=False):
    self.fft_m_real = fft_m_real
    self.fft_n_real = fft_n_real
    if (isinstance(atom_selection, flex.bool)):
      atom_selection = atom_selection.iselection()
    assert (len(atom_selection) == 1) or (not atom_selection.all_eq(0))
    if (exclude_hydrogens):
      not_hd_selection = (~(xray_structure.hd_selection())).iselection()
      atom_selection = atom_selection.intersection(not_hd_selection)
    assert (len(atom_selection) != 0)
    self.atom_selection = atom_selection
    self.sites = xray_structure.sites_cart().select(atom_selection)
    self.sites_frac = xray_structure.sites_frac().select(atom_selection)
    scatterers = xray_structure.scatterers().select(atom_selection)
    self.atom_radii = flex.double(self.sites.size(), atom_radius)
    for i_seq, sc in enumerate(scatterers):
      if (sc.element_symbol().strip().lower() in ["h","d"]):
        assert (not exclude_hydrogens)
        self.atom_radii[i_seq] = 1.0
    self.map_sel = maptbx.grid_indices_around_sites(
      unit_cell  = xray_structure.unit_cell(),
      fft_n_real = fft_n_real,
      fft_m_real = fft_m_real,
      sites_cart = self.sites,
      site_radii = self.atom_radii)

  def analyze_map(self, map, model_map=None, min=None, max=None,
      compare_at_sites_only=False):
    """
    Extract statistics for the given map, plus statistics for the model map
    if given.  The CC can either be calculated across grid points within the
    given radius of the sites, or at the sites directly.
    """
    assert (map.focus() == self.fft_n_real) and (map.all() == self.fft_m_real)
    map_sel = map.select(self.map_sel)
    map_values = flex.double()
    model_map_sel = model_map_mean = model_map_values = None
    if (model_map is not None):
      assert ((model_map.focus() == self.fft_n_real) and
              (model_map.all() == self.fft_m_real))
      model_map_sel = model_map.select(self.map_sel)
      model_map_values = flex.double()
    for site_frac in self.sites_frac:
      map_values.append(map.eight_point_interpolation(site_frac))
      if (model_map is not None):
        model_map_values.append(model_map.eight_point_interpolation(site_frac))
    cc = None
    if (model_map is not None):
      if (compare_at_sites_only):
        cc = flex.linear_correlation(x=map_values,
          y=model_map_values).coefficient()
      else :
        cc = flex.linear_correlation(x=map_sel, y=model_map_sel).coefficient()
      model_map_mean = flex.mean(model_map_values)
    n_above_max = n_below_min = None
    if (min is not None):
      n_below_min = (map_values < min).count(True)
    if (max is not None):
      n_above_max = (map_values > max).count(True)
    return group_args(
      cc=cc,
      min=flex.min(map_values),
      max=flex.max(map_values),
      mean=flex.mean(map_values),
      n_below_min=n_below_min,
      n_above_max=n_above_max,
      model_mean=model_map_mean)

def map_statistics_for_atom_selection(
    atom_selection,
    fmodel=None,
    resolution_factor=0.25,
    map1=None,
    map2=None,
    xray_structure=None,
    map1_type="2mFo-DFc",
    map2_type="Fmodel",
    atom_radius=1.5,
    exclude_hydrogens=False):
  """
  Simple-but-flexible function to give the model-to-map CC and mean density
  values (sigma-scaled, unless pre-calculated maps are provided) for any
  arbitrary atom selection.
  """
  assert (atom_selection is not None) and (len(atom_selection) > 0)
  if (fmodel is not None):
    assert (map1 is None) and (map2 is None) and (xray_structure is None)
    edm = fmodel.electron_density_map()
    map1_coeffs = edm.map_coefficients(map1_type)
    map1 = map1_coeffs.fft_map(
      resolution_factor=resolution_factor).apply_sigma_scaling().real_map()
    map2_coeffs = edm.map_coefficients(map2_type)
    map2 = map2_coeffs.fft_map(
      resolution_factor=resolution_factor).apply_sigma_scaling().real_map()
    xray_structure = fmodel.xray_structure
  else :
    assert (not None in [map1, map2, xray_structure])
    assert isinstance(map1, flex.double) and isinstance(map2, flex.double)
  if (exclude_hydrogens):
    hd_selection = xray_structure.hd_selection()
    if (type(atom_selection).__name__ == "size_t"):
      atom_selection_new = flex.size_t()
      for i_seq in atom_selection :
        if (not hd_selection[i_seq]):
          atom_selection_new.append(i_seq)
      atom_selection = atom_selection_new
      assert (len(atom_selection) > 0)
    else :
      assert (type(atom_selection).__name__ == "bool")
      atom_selection &= ~hd_selection
  manager = selection_map_statistics_manager(
    atom_selection=atom_selection,
    xray_structure=xray_structure,
    fft_n_real = map1.focus(),
    fft_m_real = map1.all(),
    exclude_hydrogens=exclude_hydrogens)
  stats = manager.analyze_map(
    map=map1,
    model_map=map2)
  return group_args(
    cc=stats.cc,
    map1_mean=stats.mean,
    map2_mean=stats.model_mean)

def map_statistics_for_fragment(fragment, **kwds):
  """
  Shortcut to map_statistics_for_atom_selection using a PDB hierarchy object
  to define the atom selection.
  """
  atoms = fragment.atoms()
  i_seqs = atoms.extract_i_seq()
  assert (not i_seqs.all_eq(0))
  return map_statistics_for_atom_selection(i_seqs, **kwds)

def find_suspicious_residues(
    fmodel,
    pdb_hierarchy,
    hetatms_only=True,
    skip_single_atoms=True,
    skip_alt_confs=True,
    min_acceptable_cc=0.8,
    min_acceptable_2fofc=1.0,
    max_frac_atoms_below_min=0.5,
    ignore_resnames=(),
    log=None):
  if (log is None) : log = null_out()
  xray_structure = fmodel.xray_structure
  assert (len(pdb_hierarchy.atoms()) == xray_structure.scatterers().size())
  edm = fmodel.electron_density_map()
  map_coeffs1 = edm.map_coefficients(
    map_type="2mFo-DFc",
    fill_missing=False)
  map1 = map_coeffs1.fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  map_coeffs2 = edm.map_coefficients(
    map_type="Fc",
    fill_missing=False)
  map2 = map_coeffs2.fft_map(
    resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()
  unit_cell = xray_structure.unit_cell()
  hd_selection = xray_structure.hd_selection()
  outliers = []
  for chain in pdb_hierarchy.models()[0].chains():
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) and (skip_alt_confs):
        continue
      for atom_group in residue_group.atom_groups():
        if (atom_group.resname in ignore_resnames):
          continue
        atoms = atom_group.atoms()
        assert (len(atoms) > 0)
        if (len(atoms) == 1) and (skip_single_atoms):
          continue
        if (hetatms_only):
          if (not atoms[0].hetero):
            continue
        map_stats = map_statistics_for_fragment(
          fragment=atom_group,
          map1=map1,
          map2=map2,
          xray_structure=fmodel.xray_structure,
          exclude_hydrogens=True)
        n_below_min = n_heavy = sum = 0
        for atom in atoms :
          if (hd_selection[atom.i_seq]):
            continue
          n_heavy += 1
          site = atom.xyz
          site_frac = unit_cell.fractionalize(site)
          map_value = map1.tricubic_interpolation(site_frac)
          if (map_value < min_acceptable_2fofc):
            n_below_min += 1
          sum += map_value
        map_mean = sum / n_heavy
        frac_below_min = n_below_min / n_heavy
        if ((map_stats.cc < min_acceptable_cc) or
            (frac_below_min > max_frac_atoms_below_min) or
            (map_mean < min_acceptable_2fofc)):
          residue_info = "%1s%3s%2s%5s" % (atom_group.altloc,
            atom_group.resname, chain.id, residue_group.resid())
          xyz_mean = atoms.extract_xyz().mean()
          outliers.append((residue_info, xyz_mean))
          print("Suspicious residue: %s" % residue_info, file=log)
          print("  Overall CC to 2mFo-DFc map = %.2f" % map_stats.cc, file=log)
          print("  Fraction of atoms where 2mFo-DFc < %.2f = %.2f" % \
            (min_acceptable_2fofc, frac_below_min), file=log)
          print("  Mean 2mFo-DFc value = %.2f" % map_mean, file=log)
  return outliers

def extract_map_stats_for_single_atoms(xray_structure, pdb_atoms, fmodel,
                                        selection=None, fc_map=None,
                                        two_fofc_map=None):
  """
  Memory-efficient routine for harvesting map values for individual atoms
  (e.g. waters).  Only one FFT'd map at a time is in memory.

  :param selection: optional atom selection (flex array)
  :returns: group_args object with various simple metrics
  """
  if (selection is None):
    selection = ~(xray_structure.hd_selection())
  sites_cart = xray_structure.sites_cart()
  sites_frac = xray_structure.sites_frac()
  unit_cell = xray_structure.unit_cell()
  def collect_map_values(map, get_selections=False):
    values = []
    selections = []
    if (map is None):
      assert (not get_selections)
      return [ None ] * len(pdb_atoms)
    for i_seq, atom in enumerate(pdb_atoms):
      if (selection[i_seq]):
        site_frac = sites_frac[i_seq]
        values.append(map.eight_point_interpolation(site_frac))
        if (get_selections):
          sel = maptbx.grid_indices_around_sites(
            unit_cell  = unit_cell,
            fft_n_real = map.focus(),
            fft_m_real = map.all(),
            sites_cart = flex.vec3_double([sites_cart[i_seq]]),
            site_radii = flex.double([1.5]))
          selections.append(map.select(sel))
      else :
        values.append(None)
        selections.append(None)
    if (get_selections):
      return values, selections
    else :
      return values
  def get_map(map_type):
    map_coeffs = fmodel.map_coefficients(map_type=map_type)
    if (map_coeffs is not None):
      return map_coeffs.fft_map(
        resolution_factor=0.25).apply_sigma_scaling().real_map_unpadded()

  # use maps
  if ( (fmodel is None) and (fc_map is not None) and
       (two_fofc_map is not None) ):
    fofc = [ 0.0 for i in range(len(pdb_atoms)) ]
    anom = [ None for i in range(len(pdb_atoms)) ]
    two_fofc, two_fofc_sel = collect_map_values(
      two_fofc_map, get_selections=True)
    f_model_val, f_model_sel = collect_map_values(
      fc_map, get_selections=True)

  # otherwise, use data to calculate maps
  else:
    two_fofc_map = get_map("2mFo-DFc")
    two_fofc, two_fofc_sel = collect_map_values(
      two_fofc_map, get_selections=True)
    del two_fofc_map
    fofc_map = get_map("mFo-DFc")
    fofc = collect_map_values(fofc_map)
    del fofc_map
    anom_map = get_map("anomalous")
    anom = collect_map_values(anom_map)
    del anom_map
    fmodel_map = get_map("Fmodel")
    f_model_val, f_model_sel = collect_map_values(
      fmodel_map, get_selections=True)
    del fmodel_map

  two_fofc_ccs = []
  for i_seq, atom in enumerate(pdb_atoms):
    if (selection[i_seq]):
      cc = flex.linear_correlation(x=two_fofc_sel[i_seq],
        y=f_model_sel[i_seq]).coefficient()
      two_fofc_ccs.append(cc)
    else :
      two_fofc_ccs.append(None)
  return group_args(
    two_fofc_values=two_fofc,
    fofc_values=fofc,
    anom_values=anom,
    fmodel_values=f_model_val,
    two_fofc_ccs=two_fofc_ccs)
