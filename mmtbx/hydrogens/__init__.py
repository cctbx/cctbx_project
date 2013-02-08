from __future__ import division
from mmtbx.utils import rotatable_bonds
from scitbx.matrix import rotate_point_around_axis
from cctbx.array_family import flex
from cctbx import maptbx


hydrogens_master_params_str = """
    refine = individual riding *Auto
      .type = choice
      .help = Choice for refinement: riding model or full (H is refined as \
              other atoms, useful at very high resolutions only)
      .short_caption = Hydrogen refinement model
      .expert_level=1
    optimize_scattering_contribution = True
      .type = bool
    contribute_to_f_calc = True
      .type = bool
      .help = Add H contribution to Xray (Fcalc) calculations
      .short_caption=Include hydrogens in Fcalc
      .expert_level=1
    high_resolution_limit_to_include_scattering_from_h = 1.6
      .type = float
      .short_caption = High-resolution limit to include scattering from H
      .expert_level=2
    real_space_optimize_x_h_orientation = True
      .type = bool
      .short_caption = Optimize X-H orientation in real-space
      .expert_level = 1
    xh_bond_distance_deviation_limit = 0.0
      .type = float
      .help = Idealize XH bond distances if deviation from ideal is greater \
              than xh_bond_distance_deviation_limit
      .short_caption = X-H bond distance deviation limit
      .expert_level=2
"""

def rotatable(pdb_hierarchy, mon_lib_srv):
  """
  General tool to identify rotatable H, such as C-O-H, C-H3, in any molecule.
  """
  result = []
  def analyze_group(g, atoms):
    for gi in g:
      assert len(gi[0])==2 # because this is axis
      assert len(gi[1])>0  # because these are atoms rotating about this axis
      # condition 1: axis does not contain H or D
      a1, a2 = atoms[gi[0][0]], atoms[gi[0][1]]
      e1 = a1.element.strip().upper()
      e2 = a2.element.strip().upper()
      condition_1 = [e1,e2].count("H")==0 and [e1,e2].count("D")==0
      # condition 2: all atoms to rotate are H or D
      condition_2 = True
      rot_atoms = []
      for gi1i in gi[1]:
        if(not atoms[gi1i].element.strip().upper() in ["H","D"]):
          condition_2 = False
          break
      rot_atoms = []
      axis = None
      if(condition_1 and condition_2):
        axis = [a1.i_seq, a2.i_seq]
        for gi1i in gi[1]:
          rot_atoms.append(atoms[gi1i].i_seq)
    if(axis is not None): return axis, rot_atoms
    else: return None
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            atoms = residue.atoms()
            fr = rotatable_bonds.axes_and_atoms_aa_specific(
              residue=residue, mon_lib_srv=mon_lib_srv,
              remove_clusters_with_all_h=False, log=None)
            if(fr is not None):
              r = analyze_group(g=fr, atoms=atoms)
              if(r is not None): result.append(r)
  return result

def fit_rotatable(pdb_hierarchy, mon_lib_srv, xray_structure, map_data):
  unit_cell = xray_structure.unit_cell()
  sel = rotatable(pdb_hierarchy=pdb_hierarchy, mon_lib_srv=mon_lib_srv)
  sites_cart = xray_structure.sites_cart()
  scatterers = xray_structure.scatterers()
  for sel_ in sel:
    ed_val = -1
    angle = 0.
    angular_step = 1
    axis = sel_[0]
    points_i_seqs = sel_[1]
    sites_frac_best = flex.vec3_double(len(points_i_seqs))
    while angle <= 360:
      sites_frac_tmp  = flex.vec3_double(len(points_i_seqs))
      ed_val_ = 0
      for i_seq, point_i_seq in enumerate(points_i_seqs):
        site_cart_new = rotate_point_around_axis(
          axis_point_1 = sites_cart[axis[0]],
          axis_point_2 = sites_cart[axis[1]],
          point        = sites_cart[point_i_seq],
          angle        = angle,
          deg          = True)
        site_frac_new = unit_cell.fractionalize(site_cart_new)
        ed_val_ += abs(maptbx.eight_point_interpolation(map_data,site_frac_new))
        sites_frac_tmp[i_seq] = site_frac_new
      if(ed_val_ > ed_val):
        ed_val = ed_val_
        sites_frac_best = sites_frac_tmp.deep_copy()
      angle += angular_step
    for i_seq, point_i_seq in enumerate(points_i_seqs):
      scatterers[point_i_seq].site = sites_frac_best[i_seq]
  pdb_hierarchy.adopt_xray_structure(xray_structure)

def run_fit_rotatable(fmodel, ref_model, angular_step, log):
  fft_map = fmodel.electron_density_map().fft_map(
    resolution_factor = 1./4.,
    map_type          = "2mFo-DFc",
    symmetry_flags    = maptbx.use_space_group_symmetry)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  xrs = fmodel.xray_structure
  from mmtbx import monomer_library
  mon_lib_srv = monomer_library.server.server() # XXX see if it's available around
  fit_rotatable(
    pdb_hierarchy  = ref_model.pdb_hierarchy(sync_with_xray_structure=True),
    mon_lib_srv    = mon_lib_srv,
    xray_structure = xrs,
    map_data       = map_data)
  fmodel.update_xray_structure(xray_structure = xrs, update_f_calc=True)
  ref_model.xray_structure = xrs
