from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math,sys
from libtbx import adopt_init_args
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.phil
from mmtbx import find_peaks
import mmtbx.utils
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from six.moves import range
from libtbx.utils import user_plus_sys_time
from libtbx.utils import null_out
from mmtbx.solvent import map_to_water
from libtbx import group_args
import string
import libtbx.log
from libtbx import Auto
from cctbx import geometry_restraints
from iotbx.pdb.mmcif import format_pdb_atom_name

def get_unique_altloc(exclude):
  for l in string.ascii_uppercase:
    if not l in exclude:
      return l

def get_unique_altloc2(available, exclude):
  l = None
  for l in available:
    if l in [' ', '']: continue
    if not l in exclude:
      return l
  return l

output_params_str = """
  output_residue_name = Auto
    .type=str
    .input_size = 50
  output_chain_id = S
    .type=str
    .input_size = 50
  output_atom_name = O
    .type=str
    .input_size = 50
  scattering_type = O
    .type=str
    .help = Defines scattering factors for newly added waters
    .expert_level=2
    .input_size = 50
"""

h_bond_params_str = """
  dist_min = 1.8
    .type = float
    .short_caption = Min distance between water and any atom
    .expert_level = 1
  dist_max = 3.2
    .type = float
    .short_caption = Max distance between water and any atom
    .expert_level = 1
  dist_min_altloc = 0.5
    .type = float
"""

adp_occ_params_str = """
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
    .short_caption = New solvent ADP type
    .expert_level = 1
  b_iso_min = 1.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
    .short_caption = Minimum B-factor
    .expert_level = 1
  b_iso_max = 80.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
    .short_caption = Maximum B-factor
    .expert_level = 1
  anisotropy_min = 0.1
    .type = float
    .help = For solvent refined as anisotropic: remove is less than this value
    .short_caption = Minimum anisotropic B-factor
    .expert_level = 1
  occupancy_min = 0.2
    .type=float
    .help = Minimum occupancy value, waters with smaller value will be rejected
    .short_caption = Minimum occupancy
  occupancy_max = 1.0
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
    .short_caption = Maximum occupancy
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
    .short_caption = Initial occupancy value
"""

master_params_str = """\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
    .short_caption = Minimum resolution
  mode = *second_half filter_only every_macro_cycle every_macro_cycle_after_first
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
    .short_caption = Mode
  add_hydrogens
  {
    d_min = 1.0
      .type = float
      .help = None = don't add. Resolution must be better than d_min for H building
    r_work = 0.1
      .type = float
      .help = Rwork must be less than r_work for water to be added
    min_dist = 0.5
      .type = float
    max_dist = 1.5
      .type = float
    min_angle = 80.0
      .type = float
    max_angle = 170.0
      .type = float
    diff_map_peak_threshold = 2.0
      .type = float
      .help = mFo-DFc map threshold for H peak search
  }
  mask_atoms_selection = protein and (name CA or name CB or name N or name C or name O)
    .type = str
    .help = Mask macromolecule atoms in peak picking map
  include_altlocs = False
    .type = bool
    .help = Search for water with altlocs
  n_cycles = 1
    .type = int
    .short_caption = Number of cycles
  %s
  primary_map_type = mFobs-DFmodel
    .type=str
    .help = Map used to identify candidate water sites - by default this is \
      the standard difference map.
  primary_map_cutoff = 3.
    .type=float
    .short_caption = Primary map cutoff (sigma)
  secondary_map_and_map_cc_filter
    .short_caption = Secondary map filter
    .style = auto_align box
  {
    cc_map_1_type = "Fmodel"
      .type = str
      .short_caption = Model map type for CC calculation
    cc_map_2_type = 2mFobs-DFmodel
      .type = str
      .short_caption = Experimental map type for CC calculation
    poor_cc_threshold = 0.70
      .type = float
      .short_caption = Minimum map correlation
    poor_map_value_threshold = 1.0
      .type = float
      .short_caption = Minimum map value (sigma)
  }
  %s
  refine_oat = False
    .type = bool
    .help = Q & B coarse grid search.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_occupancies = True
    .type = bool
    .help = Refine solvent occupancies.
    .expert_level = 1
  %s
  filter_at_start = True
    .type = bool
    .expert_level = 1
  ignore_final_filtering_step = False
    .type = bool
    .expert_level=2
  correct_drifted_waters = True
    .type = bool
    .expert_level=2
""" % (output_params_str, h_bond_params_str, adp_occ_params_str)

def master_params():
  return iotbx.phil.parse(master_params_str)

def new_solvent_sites_as_hierarchy_chain(
      model,
      sites,
      params,
      conformer_indices=None):
  uc = model.crystal_symmetry().unit_cell()
  ort = uc.orthogonalize
  b_iso = flex.mean(model.get_xray_structure(
    ).extract_u_iso_or_u_equiv()) * adptbx.u_as_b(1) * 1.5
  u_star = adptbx.u_iso_as_u_star(uc, adptbx.b_as_u(b_iso))
  u_cart = adptbx.u_star_as_u_cart(uc, u_star)
  size = model.size()
  # Figure out resname
  resname = params.output_residue_name
  if resname is Auto or resname is None:
    resname = model.get_hierarchy().get_water_resname()
    if resname is None:
      table = model.get_xray_structure().get_scattering_table()
      if table=="neutron": resname = "DOD"
      else:                resname = "HOH"
  # Figure out first water resseq
  first_water_resseq = model.get_hierarchy().get_water_max_resseq()+1
  #
  # Figure out atom name
  element = "O"
  atom_name = params.output_atom_name
  if atom_name is Auto or atom_name is None:
    atom_name = element
  atom_name = format_pdb_atom_name(atom_name, element)
  #
  new_chain = iotbx.pdb.hierarchy.chain(id=params.output_chain_id)
  for i_seq, site_frac in enumerate(sites):
    water_i_seq = size+i_seq+1 # same as resseq below
    new_atom = (iotbx.pdb.hierarchy.atom()
      .set_serial(new_serial=iotbx.pdb.hy36encode(width=5, value=water_i_seq))
      .set_name(new_name=atom_name)
      .set_xyz(new_xyz=ort(site_frac))
      .set_occ(new_occ=1)
      .set_element(element)
      .set_charge("")
      .set_hetero(new_hetero=True))
    if(  params.new_solvent == "isotropic"):
      new_atom.set_b(new_b=b_iso)
    elif(params.new_solvent == "anisotropic"):
      new_atom.set_uij(new_uij=u_cart)
    else: raise RuntimeError
    altloc=""
    if conformer_indices is not None:
      ci = conformer_indices.conformer_indices[i_seq]
      cm = conformer_indices.index_altloc_mapping
      if not ci in cm.values() and ci==0: altloc = ""
      else:
        altloc = list(cm.keys())[list(cm.values()).index(ci)]
    new_atom_group = iotbx.pdb.hierarchy.atom_group(
      altloc  = altloc,
      resname = resname)
    new_atom_group.append_atom(atom=new_atom)
    new_residue_group = iotbx.pdb.hierarchy.residue_group(
      resseq=iotbx.pdb.resseq_encode(value=first_water_resseq), icode=" ")
    first_water_resseq += 1
    new_residue_group.append_atom_group(atom_group=new_atom_group)
    new_chain.append_residue_group(residue_group=new_residue_group)
  if new_chain.residue_groups_size()>0: return new_chain
  else:                                 return None

def add_solvent_to_model_inplace(
      model,
      params,
      sites=None,
      water_residue_groups=None, # come together
      chain_id=None,             # come together
      conformer_indices=None):
  """
  A series of very specialzied manipulations on the model to add new solvent
  without re-creating the new model.
  """
  assert sites is None or [water_residue_groups, chain_id].count(None)==2
  # Get new solvent as pdb hierarchy chain
  if sites is None:
    new_solvent_chain = create_water_chain(
      water_residue_groups=water_residue_groups,
      chain_id=chain_id,
      start_resseq=1)
  else:
    new_solvent_chain = new_solvent_sites_as_hierarchy_chain(
      sites             = sites,
      model             = model,
      params            = params,
      conformer_indices = conformer_indices)
  if new_solvent_chain is None: return
  # Get new solvent selection
  n_new_solvent      = new_solvent_chain.atoms().size()
  initial_model_size = model.size()
  new_solvent_selection = flex.bool(initial_model_size, False)
  new_solvent_selection.extend(flex.bool(n_new_solvent, True))
  # Append new chain to the model
  #model.get_hierarchy().only_model().append_chain(chain=new_solvent_chain)
  add_chain_to_hierarchy(
    hierarchy = model.get_hierarchy(), new_chain = new_solvent_chain)
  model._update_atom_selection_cache()
  model.get_hierarchy().atoms().reset_i_seq()
  model.unset_processed_pdb_file()
  # Force-update xray_structure
  # This is done THIS WAY to keep scattering_table
  solvent_xray_structure = new_solvent_chain.as_new_hierarchy(
    ).extract_xray_structure(crystal_symmetry=model.crystal_symmetry())
  model._xray_structure = model._xray_structure.concatenate(
    solvent_xray_structure)
  sst = model.get_xray_structure().select(
    new_solvent_selection).site_symmetry_table()
  #
  # Adjust refinement flags
  #
  rfs = model.refinement_flags
  if rfs is not None:
    scatterers = solvent_xray_structure.scatterers()
    #for i_sc, sc in enumerate(solvent_xray_structure.scatterers()):
    for i, atom in enumerate(list(new_solvent_chain.atoms())):
      # occupancy
      if(params.refine_occupancies):
        sel = [[flex.size_t([initial_model_size+i])], ]
        if rfs.s_occupancies is None:
          rfs.s_occupancies = sel
        else:
          rfs.s_occupancies.extend(sel)
          #rfs.s_occupancies.append(sel)
      # sites
      if(rfs.individual_sites):
        rfs.sites_individual.append(True)
      # ADP
      if(rfs.individual_adp): # H?
        iso = rfs.adp_individual_iso
        ani = rfs.adp_individual_aniso

        # if water happened to be anisotropic already (= not really new!)
        if scatterers[i].flags.use_u_aniso():
          if iso is not None: iso.append(False)
          if ani is not None: ani.append(True)
          continue

        if(params.new_solvent == "isotropic"):
          if iso is not None: iso.append(True)
          if ani is not None: ani.append(False)
        if(params.new_solvent == "anisotropic"):
          if iso is not None: iso.append(False)
          if ani is not None: ani.append(True)
  #
  nonbonded_types = flex.std_string()
  for atom in new_solvent_chain.atoms():
    if(atom.element_is_hydrogen()): nonbonded_types.append('HOH2')
    else:                           nonbonded_types.append('OH2')
  #
  # Finally, add new atoms to restraints manager
  #
  rm = model.get_restraints_manager()
  if(rm is not None):
    geometry = rm.geometry
    if (geometry.model_indices is None):
      model_indices = None
    else:
      model_indices = flex.size_t(n_new_solvent, 0)
    if(geometry.conformer_indices is None):
      conformer_indices = None
    else:
      if conformer_indices is not None:
        assert conformer_indices.conformer_indices.size() == n_new_solvent
        conformer_indices = conformer_indices.conformer_indices
      else:
        conformer_indices = flex.size_t(n_new_solvent, 0)
    if (geometry.sym_excl_indices is None):
      sym_excl_indices = None
    else:
      sym_excl_indices = flex.size_t(n_new_solvent, 0)
    if (geometry.donor_acceptor_excl_groups is None):
      donor_acceptor_excl_groups = None
    else:
      donor_acceptor_excl_groups = flex.size_t(n_new_solvent, 0)
    geometry = geometry.new_including_isolated_sites(
      n_additional_sites         = n_new_solvent,
      model_indices              = model_indices,
      conformer_indices          = conformer_indices,
      sym_excl_indices           = sym_excl_indices,
      donor_acceptor_excl_groups = donor_acceptor_excl_groups,
      site_symmetry_table        = sst,
      nonbonded_types            = nonbonded_types,
      nonbonded_charges          = flex.int(n_new_solvent, 0))
    model.restraints_manager = mmtbx.restraints.manager(
      geometry              = geometry,
      cartesian_ncs_manager = rm.cartesian_ncs_manager,
      normalization         = rm.normalization)
    c_ncs_m = model.get_cartesian_NCS_manager()
    if(c_ncs_m is not None):
      c_ncs_m.register_additional_isolated_sites(
        number=n_new_solvent)
    model.restraints_manager.geometry.update_plain_pair_sym_table(
      sites_frac = model.get_sites_frac())
  #
  if model.get_restraints_manager() is not None:
    aproxies, bproxies = [], []
    get_class = iotbx.pdb.common_residue_names_get_class
    sct = model.get_scattering_table()
    if sct in ["neutron", "electron"]: target_bond = 0.98
    else:                              target_bond = 0.85
    for m in model.get_hierarchy().models():
      for c in m.chains():
        for conf in c.conformers():
          for r in conf.residues():
            if(not get_class(name=r.resname) == "common_water"): continue
            i_seqs = r.atoms().extract_i_seq()
            if i_seqs.size() != 3: continue
            io = None
            ihs = []
            for atom in r.atoms():
              if(atom.element_is_hydrogen()): ihs.append(atom.i_seq)
              else: io = atom.i_seq
            #
            if io is not None and len(ihs)==2:
              for ih in ihs:
                bproxy = geometry_restraints.bond_simple_proxy(
                  i_seqs         = (ih, io),
                  distance_ideal = target_bond,
                  weight         = 1./(0.02**2),
                  origin_id      = 0)
                bproxies.append(bproxy)
              aproxy = geometry_restraints.angle_proxy(
                i_seqs      = (ihs[0], io, ihs[1]),
                angle_ideal = 103.91,
                weight      = 1./(3.0**2),
                origin_id   = 0)
              aproxies.append(aproxy)
  if model.get_restraints_manager() is not None:
    g = model.get_restraints_manager().geometry
    g.add_new_bond_restraints_in_place(bproxies, model.get_sites_cart())
    g.add_angles_in_place(aproxies)
  #
  # Update H riding manager
  #
  if model.riding_h_manager is not None:
    new_riding_h_manager = model.riding_h_manager.update(
      pdb_hierarchy       = model.get_hierarchy(),
      geometry_restraints = geometry,
      n_new_atoms         = n_new_solvent)
    model.riding_h_manager = new_riding_h_manager
  return new_solvent_selection

from scitbx import matrix

def find_water_hydrogens(
      o_site,
      peak_heights,
      peak_coords,
      min_od=0.5,
      max_od=1.5,
      min_angle=80.0,
      max_angle=170.0):
    """
    AI generated.

    Finds the two most likely hydrogen positions for a water oxygen
    using Cartesian coordinates (sites_cart).
    Return None if less than 2 H are found.

    Prompt:
    I need a function that takes the position of a water oxygen and two lists:
    (1) the values of map peaks and (2) their coordinates.
    It should return the coordinates of the two points that most likely
    correspond to the positions of the two water hydrogen atoms. These hydrogen
    positions should be selected from the provided coordinate and peak heights
    list such that:

    a) They both correspond to the highest map peaks.

    b) The H-O-H angle is close to the expected value for water (about 103deg).

    c) The O-H distances are close to 1 A (approximately).
    """
    xyz_h = list(zip(peak_heights, peak_coords))
    if len(xyz_h) < 2: return None
    o_vec = matrix.col(o_site)
    # Identify H1: Filter candidates in Cartesian space
    valid_h1_candidates = [
      p for p in xyz_h
      if min_od < (o_vec - matrix.col(p[1])).length() < max_od
    ]
    if not valid_h1_candidates: return None
    # Explicitly extract the peak height (x[0]) to avoid vec3 comparison crashes
    best_h1 = max(valid_h1_candidates, key=lambda x: x[0])
    h1_height, h1_site = best_h1
    h1_vec = matrix.col(h1_site)
    d1 = (o_vec - h1_vec).length()
    xyz_h.remove(best_h1)
    # Trackers for H2 parameters
    h2_site = None
    h2_height = None
    d2_final = None
    angle_final = None
    best_angle_score = 1e10
    h_max = -1e300
    # Iterate through remaining peaks to find the best H2
    for h, s in xyz_h:
      s_vec = matrix.col(s)
      d2 = (o_vec - s_vec).length()
      if d2 < min_od or d2 > max_od: continue
      # Cartesian angle calculation: v1.angle(v2, deg=True)
      a = (s_vec - o_vec).angle((h1_vec - o_vec), deg=True)
      if min_angle < a < max_angle:
        angle_deviation = abs(a - 103.91)
        if h > h_max:
          h_max = h
          best_angle_score = angle_deviation
          h2_site = s
          h2_height = h
          d2_final = d2
          angle_final = a
        elif h == h_max and angle_deviation < best_angle_score:
          best_angle_score = angle_deviation
          h2_site = s
          h2_height = h
          d2_final = d2
          angle_final = a
    if h2_site is None: return None
    # Return the packaged group_args object
    return group_args(
      h_sites   = (h1_site, h2_site),
      map_peaks = (h1_height, h2_height),
      distances = (d1, d2_final),
      angle     = angle_final
    )

def get_detached_water_residue_group(o_atom, h_sites_cart):
    # Detach at the residue_group level (parent of the atom_group)
    residue_group = o_atom.parent().parent().detached_copy()
    # Get the corresponding atom_group in the new detached residue_group
    # (Assuming no alternate conformations for water, so index 0)
    atom_group = residue_group.atom_groups()[0]
    detached_o_atom = atom_group.atoms()[0]
    # Create and append the hydrogens
    for i, h_xyz in enumerate(h_sites_cart):
      h_atom = detached_o_atom.detached_copy()
      h_atom.name = f" H{i+1} "
      h_atom.element = "H"
      h_atom.xyz = h_xyz
      h_atom.uij = (-1, -1, -1, -1, -1, -1)
      atom_group.append_atom(h_atom)
    # Return the single residue_group representing the whole water molecule
    return residue_group

def create_water_chain(water_residue_groups, chain_id, start_resseq=1):
  """
  Creates and returns a new iotbx.pdb.hierarchy.chain populated with
  the provided detached water residue groups.
  """
  new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
  current_resseq = start_resseq
  for rg in water_residue_groups:
    new_rg = rg.detached_copy()
    # Assign a sequential sequence number to the water
    new_rg.resseq = iotbx.pdb.resseq_encode(value=current_resseq)
    new_rg.icode = " "
    current_resseq += 1
    # Enforce water-specific attributes
    for ag in new_rg.atom_groups():
      for atom in ag.atoms():
        atom.hetero = True
    # Suppresses the insertion of 'TER' records
    new_rg.link_to_previous = True
    new_chain.append_residue_group(residue_group=new_rg)
  return new_chain

def add_chain_to_hierarchy(hierarchy, new_chain):
  """
  Appends a new chain to the hierarchy as a completely separate chain object,
  even if a chain with the same ID already exists.
  It safely renumbers the new chain to prevent numbering collisions without
  disturbing existing data structures.
  """
  # Find the highest existing resseq for this chain ID across the hierarchy
  highest_resseq = 0
  for model in hierarchy.models():
    for chain in model.chains():
      if chain.id == new_chain.id:
        for rg in chain.residue_groups():
          try:
            seq_val = rg.resseq_as_int()
            if seq_val > highest_resseq:
              highest_resseq = seq_val
          except ValueError:
            pass
  # Renumber the new chain's residues to prevent collisions, preserving altlocs
  current_resseq = highest_resseq + 1
  previous_incoming_resid = None
  for rg in new_chain.residue_groups():
    incoming_resid = rg.resid()
    if(previous_incoming_resid is not None and
       incoming_resid != previous_incoming_resid):
      current_resseq += 1
    rg.resseq = iotbx.pdb.resseq_encode(value=current_resseq)
    rg.icode = " "
    previous_incoming_resid = incoming_resid
  # Append the new chain as a separate entity (least intrusive method)
  if hierarchy.models_size() == 0:
    hierarchy.append_model(iotbx.pdb.hierarchy.model())
  hierarchy.models()[0].append_chain(new_chain.detached_copy())
  # Re-index the atoms in the hierarchy
  hierarchy.atoms().reset_i_seq()
  return hierarchy

class maps(object):
  def __init__(self,
               fmodel,
               model,
               mask_atoms_selection,
               difference_map_type,
               model_map_type,
               data_map_type,
               grid_step=0.6,
               radius=2.0):
    self.radius = radius
    self.fmodel = fmodel
    self.model  = model
    self.r_work = fmodel.r_work()
    self.e_map = fmodel.electron_density_map()
    self.crystal_symmetry = fmodel.xray_structure.crystal_symmetry()
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = self.crystal_symmetry.unit_cell(),
      space_group_info = self.crystal_symmetry.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = grid_step)
    self.difference_map = self._get_real_map(map_type = difference_map_type)
    self.model_map      = self._get_real_map(map_type = model_map_type)
    self.data_map       = self._get_real_map(map_type = data_map_type)
    #self._estimate_diff_map_cutoff()
    # Compute mask in P1 to mask out desired regions
    if mask_atoms_selection is not None:
      bb_sel = self.model.selection(string=mask_atoms_selection)
      if bb_sel.count(True)>0:
        xrs = fmodel.xray_structure.select(bb_sel)
        #
        # Both ways to compute mask should be the same, but they are slighly not,
        # expectedly.
        #
        mask_p1 = mmtbx.masks.mask_from_xray_structure(
          xray_structure           = xrs,
          p1                       = True,
          for_structure_factors    = True,
          solvent_radius           = 0,
          shrink_truncation_radius = 0,
          atom_radius              = 1.2,
          n_real                   = self.crystal_gridding.n_real(),
          in_asu                   = False).mask_data
        maptbx.unpad_in_place(map=mask_p1)
        sel0 = mask_p1 > 0.1
        mask_p1 = mask_p1.set_selected(sel0,  1)
        mask_p1 = mask_p1.set_selected(~sel0, 0)

        #xrs = xrs.expand_to_p1()
        #import boost_adaptbx.boost.python as bp
        #cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
        #radii = flex.double(xrs.scatterers().size(), 1.2)
        #mask_p1 = cctbx_maptbx_ext.mask(
        #  sites_frac                  = xrs.sites_frac(),
        #  unit_cell                   = xrs.unit_cell(),
        #  n_real                      = self.crystal_gridding.n_real(),
        #  mask_value_inside_molecule  = 0,
        #  mask_value_outside_molecule = 1,
        #  radii                       = radii,
        #  wrapping                    = True)

        self.difference_map = self.difference_map * mask_p1

  def _get_real_map(self, map_type):
    coeffs = self.e_map.map_coefficients(
      map_type     = map_type,
      fill_missing = False,
      isotropize   = True)
    fft_map = coeffs.fft_map(crystal_gridding = self.crystal_gridding)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

  def score_atom(self, atom, min_cc, min_value):
    site_cart = atom.xyz
    uc = self.crystal_symmetry.unit_cell()
    site_frac = uc.fractionalize(site_cart)
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = uc,
      fft_n_real = self.model_map.focus(),
      fft_m_real = self.model_map.all(),
      sites_cart = flex.vec3_double([site_cart]),
      site_radii = flex.double([self.radius]))
    cc = flex.linear_correlation(
      x=self.model_map.select(sel),
      y=self.data_map.select(sel)).coefficient()
    value_2 = self.data_map.eight_point_interpolation(site_frac)
    diff_map_val = self.difference_map.eight_point_interpolation(site_frac)
    cutoff = -3
    result = (cc > min_cc and value_2 > min_value*atom.occ) and \
             not (diff_map_val < cutoff)
    return group_args(result = result, cc=cc, value_2=value_2)

  def _estimate_diff_map_cutoff(self):
    scatterers = self.fmodel.xray_structure.scatterers()
    sel = flex.random_bool(scatterers.size(), 0.1)
    result = flex.double()
    cntr=0
    for s, sc in zip(sel, scatterers):
      if not s: continue
      occ = sc.occupancy
      sc.occupancy=0
      self.fmodel.update_xray_structure(update_f_calc=True)

      coeffs = self.fmodel.electron_density_map().map_coefficients(
        map_type     = "mFobs-DFmodel",
        fill_missing = False,
        isotropize   = True)
      fft_map = coeffs.fft_map(crystal_gridding = self.crystal_gridding)
      fft_map.apply_sigma_scaling()
      map_data = fft_map.real_map_unpadded()
      mv = map_data.tricubic_interpolation(sc.site)
      result.append(mv)

      print(self.fmodel.r_work(), mv)
      sc.occupancy=occ
      self.fmodel.update_xray_structure(update_f_calc=True)
      print(self.fmodel.r_work())
      print()
      if cntr>100: break
      cntr+=1
    mean = flex.mean(result)
    if mean/2 > 3: cutoff = 3

def fix_altlocs_and_filter(model, dist_min=1.8, fix_only=False):
  present_altlocs = list(
    model.get_hierarchy().get_conformer_indices().index_altloc_mapping.keys())
  sps = model.crystal_symmetry().special_position_settings()
  get_class = iotbx.pdb.common_residue_names_get_class
  only_model = model.get_hierarchy().only_model()
  eps = 1.e-3
  dist_min = dist_min-eps
  sites_cart = model.get_sites_cart()
  atoms = only_model.atoms()
  remove_sel = flex.size_t()
  for agi in only_model.atom_groups(): # loop over water
    if(not get_class(name=agi.resname) == "common_water"): continue
    for ai in agi.atoms():
      if ai.element_is_hydrogen(): continue # skip water H
      #
      selection_around_ai = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_min, i_seq=ai.i_seq)
      if len(selection_around_ai) == 0: continue
      #
      altlocs_inside = []
      for j in selection_around_ai:
        altlocs_inside.append( atoms[j].parent().altloc )
      #
      skip = False
      for j in selection_around_ai:
        aj = atoms[j]
        if aj.element_is_hydrogen(): continue
        agj = aj.parent()
        if(agj.altloc in [' ', '']):
          if(get_class(name=agj.resname) != "common_water"):
            remove_sel.extend(agi.atoms().extract_i_seq())
            skip=True
          else:
            new_altloc = get_unique_altloc2(
              available = present_altlocs,
              exclude   = altlocs_inside+[agi.altloc])
            if new_altloc is not None: agj.altloc = new_altloc
            else:
              remove_sel.extend(agi.atoms().extract_i_seq())
              skip = True
      if skip: continue
      #
      altlocs_inside = []
      for j in selection_around_ai:
        altlocs_inside.append( atoms[j].parent().altloc )
      #
      if agi.altloc in altlocs_inside or agi.altloc in [' ', '']:
        #print(ai.i_seq, selection_around_ai, altlocs_inside, [agi.altloc])
        #new_altloc = get_unique_altloc(exclude=altlocs_inside+[agi.altloc])
        new_altloc = get_unique_altloc2(
              available = present_altlocs,
              exclude   = altlocs_inside+[agi.altloc])
        if new_altloc is not None:
          agi.altloc = new_altloc
        else:
          remove_sel.extend(agi.atoms().extract_i_seq())
  #
  if remove_sel.size() > 0 and not fix_only:
    remove_sel = ~flex.bool(model.size(), remove_sel)
    model = model.select(selection = remove_sel)
  return model

def get_sphere_selection(
      sites_cart_all, special_position_settings, radius, i_seq):
  sel = flex.bool(sites_cart_all.size(), False)
  sel[i_seq] = True
  selection_around_i_seq = special_position_settings.pair_generator(
    sites_cart      = sites_cart_all,
    distance_cutoff = radius
      ).neighbors_of(primary_selection = sel).iselection()
  selection_around_i_seq = list(selection_around_i_seq)
  if len(selection_around_i_seq) > 0:
    selection_around_i_seq.remove(i_seq)
  return selection_around_i_seq

def filter_by_distance(model, fix_altlocs_and_filter_was_run, dist_min=1.8,
                       dist_max=3.2):
  interaction_selection = model.selection(
    map_to_water.selection_string_interaction)
  sps = model.crystal_symmetry().special_position_settings()
  get_class = iotbx.pdb.common_residue_names_get_class
  only_model = model.get_hierarchy().only_model()
  sites_cart = model.get_sites_cart()
  atoms = only_model.atoms()
  remove_sel = flex.size_t()
  for agi in only_model.atom_groups(): # loop over water
    if(not get_class(name=agi.resname) == "common_water"): continue
    for ai in agi.atoms():
      if ai.element_is_hydrogen(): continue # skip water H
      # Get selections
      selection_around_ai_min = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_min, i_seq=ai.i_seq)
      selection_around_ai_max = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_max, i_seq=ai.i_seq)
      selection_shell = list( set(selection_around_ai_max) -
                              set(selection_around_ai_min) )
      #
      altloc_i = ai.parent().altloc.strip()
      # Make sure anything inside smaller sphere are altlocs
      for j in selection_around_ai_min:
        aj = atoms[j]
        altloc_j = aj.parent().altloc
        if fix_altlocs_and_filter_was_run:
          if altloc_i =="" or altloc_i==" " or altloc_j =="" or altloc_j==" " or altloc_i==altloc_j:
            remove_sel.extend(agi.atoms().extract_i_seq())

      # Check water inside shell dist_min < dist < dist_max
      found = False
      for j in selection_shell:
        if not interaction_selection[j]: continue
        aj = atoms[j]
        if aj.element_is_hydrogen(): continue
        agj = aj.parent()
        altloc_j = agj.altloc.strip()
        if altloc_i=="" or altloc_j=="":        found = True
        if altloc_i==altloc_j and altloc_i!="": found = True
        if get_class(name=agj.resname) == "common_water":
          if altloc_i!="" and altloc_j!="": found = True
      if not found:
        remove_sel.extend(agi.atoms().extract_i_seq())
  #
  if remove_sel.size() > 0:
    remove_sel = ~flex.bool(model.size(), remove_sel)
    model = model.select(selection = remove_sel)
  return model

class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     is_neutron_scat_table,
                     params = master_params().extract(),
                     find_peaks_params = None,
                     log = sys.stdout):
    adopt_init_args(self, locals())
    # XXX Rationalize this:
    self.find_peaks_params.map_next_to_model.min_peak_peak_dist=self.params.dist_max
    if self.params.include_altlocs:
      self.find_peaks_params.peak_search.min_cross_distance=0.5
      self.find_peaks_params.map_next_to_model.min_model_peak_dist=0.5
      self.find_peaks_params.map_next_to_model.min_peak_peak_dist=0.5
    #
    self.ma         = libtbx.log.manager(log = self.log)
    self.total_time = 0
    self._maps      = None
    self._peaks     = None
    self.n_water    = None
    self.model_size_init = self.model.size()
    self.new_solvent_selection = None
    #
    self._call(msg="Start",             func=None)
    self._call(msg="Filter (dist)",     func=self._filter_dist_fix_altlocs)
    self._call(msg="Filter (q & B)",    func=self._filter_q_b)
    self._call(msg="Compute maps",      func=self._get_maps)
    self._call(msg="Filter (map)",      func=self._filter_map)
    self._call(msg="Find peaks",        func=self._find_peaks)
    self._call(msg="Add new water",     func=self._add_new_solvent)
    self._call(msg="Refine new water",  func=self._refine)
    self._call(msg="Filter (q & B)",    func=self._filter_q_b)
    self._call(msg="Compute maps",      func=self._get_maps)
    self._call(msg="Build H",           func=self._build_h)
    self._call(msg="Filter (dist only)",func=self._filter_dist)
    #self._call(msg="Correct drifted",   func=self._correct_drifted_waters)

  def _call(self, msg, func = None):
    timer = user_plus_sys_time()
    self.ma.add_and_show(msg)
    self._assert_same_model()
    if(func is not None): func()
    self._get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps()
    self._assert_same_model()
    t = timer.elapsed()
    self.total_time += t
    self._add_to_message(this_step_time=t)

  def _add_to_message(self, this_step_time):
    rs="r_work=%6.4f r_free=%6.4f"%(self.fmodel.r_work(), self.fmodel.r_free())
    nw="n_water=%3d"%(self.n_water)
    self.ma.add_and_show("  %s | %s | time (s): %s (total time: %s)"%(rs, nw,
      ("%8.3f"%this_step_time).strip(), ("%8.3f"%self.total_time).strip()))

  def _get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps(self):
    n_water = self.model.solvent_selection().count(True)
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = True)
    self._get_maps()
    self.n_water = n_water

  def _assert_same_model(self):
    mmtbx.utils.assert_xray_structures_equal( # XXX MAKE METHOD OF XRS
      x1 = self.model.get_xray_structure(),
      x2 = self.fmodel.xray_structure,
      eps=1.e-3)

  def _build_h(self):
    ahp = self.params.add_hydrogens
    if ahp.d_min is None: return
    if self.fmodel.f_obs().d_min() > ahp.d_min: return
    if self.fmodel.r_work() > ahp.r_work: return
    mpl = maptbx.MapPeakLocator(
      map_data    = self._maps.difference_map,
      unit_cell   = self.model.crystal_symmetry().unit_cell(),
      is_periodic = True)
    get_class = iotbx.pdb.common_residue_names_get_class
    self.model.get_hierarchy().atoms().reset_i_seq()
    self.model.get_hierarchy().atoms_reset_serial()
    hoh_residue_groups = []
    keep_sel = flex.bool(self.model.size(), True)
    chain_ids = []
    for m in self.model.get_hierarchy().models():
      for c in m.chains():
        for conf in c.conformers():
          for r in conf.residues():
            if(not get_class(name=r.resname) == "common_water"): continue
            i_seqs = r.atoms().extract_i_seq()
            # skip water with H
            has_h = False
            for atom in r.atoms():
              if(atom.element_is_hydrogen()): has_h = True
            if has_h: continue
            #
            for atom in r.atoms():
              if(atom.element.strip().upper()=="O"):
                sites_cart, peaks_heights = mpl.get_peaks_within_radius(
                  target_cart = atom.xyz,
                  R           = ahp.max_dist,
                  threshold   = ahp.diff_map_peak_threshold)
                wo = find_water_hydrogens(
                  o_site       = atom.xyz,
                  peak_heights = peaks_heights,
                  peak_coords  = sites_cart,
                  min_od       = ahp.min_dist,
                  max_od       = ahp.max_dist,
                  min_angle    = ahp.min_angle,
                  max_angle    = ahp.max_angle)
                if wo is not None:
                  hoh_rg = get_detached_water_residue_group(
                    o_atom       = atom,
                    h_sites_cart = wo.h_sites)
                  hoh_residue_groups.append(hoh_rg)
                  keep_sel[atom.i_seq] = False
                  chain_ids.append(atom.parent().parent().parent().id)
    if len(hoh_residue_groups) > 0:
      self.model = self.model.select(keep_sel)
      assert len(list(set(chain_ids)))==1, list(set(chain_ids))
      add_solvent_to_model_inplace(
        model                = self.model,
        params               = self.params,
        water_residue_groups = hoh_residue_groups, # come together
        chain_id             = chain_ids[0],       # come together
        conformer_indices    = None)

  def _get_maps(self):
    p = self.params
    self._maps = maps(
      fmodel               = self.fmodel,
      model                = self.model,
      mask_atoms_selection = p.mask_atoms_selection,
      difference_map_type  = p.primary_map_type,
      model_map_type       = p.secondary_map_and_map_cc_filter.cc_map_1_type,
      data_map_type        = p.secondary_map_and_map_cc_filter.cc_map_2_type)

  def _filter_dist_fix_altlocs(self):
    if self.params.include_altlocs:
      self.model = fix_altlocs_and_filter(
        model    = self.model,
        dist_min = self.params.dist_min)
    self.model = filter_by_distance(
      model                          = self.model,
      fix_altlocs_and_filter_was_run = self.params.include_altlocs,
      dist_min                       = self.params.dist_min,
      dist_max                       = self.params.dist_max)

  def _filter_dist(self):
    self.model = filter_by_distance(
      model                          = self.model,
      fix_altlocs_and_filter_was_run = self.params.include_altlocs,
      dist_min                       = self.params.dist_min,
      dist_max                       = self.params.dist_max)

  def _filter_q_b(self):
    self._filter(filter_occ=True, filter_adp=True)

  def _filter_map(self):
    self._filter(filter_map=True)

  def _filter(self,
              filter_map=False,
              filter_occ=False,
              filter_adp=False):
    mfp = self.params.secondary_map_and_map_cc_filter
    selection = flex.bool(self.model.size(), True)
    get_class = iotbx.pdb.common_residue_names_get_class
    scatterers = self.model.get_xray_structure().scatterers()
    occ = scatterers.extract_occupancies()
    b_isos = scatterers.extract_u_iso_or_u_equiv(
      self.model.get_xray_structure().unit_cell()) * math.pi**2*8
    anisotropy = scatterers.anisotropy(unit_cell =
      self.model.get_xray_structure().unit_cell())
    for m in self.model.get_hierarchy().models():
      for c in m.chains():
        for conf in c.conformers():
          for r in conf.residues():
            if(not get_class(name=r.resname) == "common_water"): continue
            i_seqs = r.atoms().extract_i_seq()
            keep = True
            has_oxygen = False # catch only H/D with no O water
            for atom in r.atoms():
              if(atom.element.strip().upper()=="O"): has_oxygen = True
              if(atom.element_is_hydrogen()): continue
              i_seq = atom.i_seq
              # Occupancy
              if filter_occ:
                if(atom.occ > self.params.occupancy_max or
                   atom.occ < self.params.occupancy_min): keep = False
                assert approx_equal(atom.occ, occ[i_seq], 1.e-3)
              # ADP
              if filter_adp:
                if(anisotropy[i_seq] < self.params.anisotropy_min): keep = False
                if(b_isos[i_seq] < self.params.b_iso_min or
                   b_isos[i_seq] > self.params.b_iso_max): keep = False
              #
              if filter_map:
                good_map = self._maps.score_atom(
                  atom      = atom,
                  min_cc    = mfp.poor_cc_threshold,
                  min_value = mfp.poor_map_value_threshold)
                if(not good_map.result):
                  keep = False
            if(not has_oxygen): keep=False
            if(not keep):
              selection = selection.set_selected(i_seqs, False)
    self.model = self.model.select(selection)

  def _refine(self):
    if(self.params.mode == "filter_only"): return
    if(self.model.size() == self.model_size_init or
       self.n_water==0): return
    for i in range(self.params.n_cycles):
      self.refine_oat()

  def _find_peaks(self):
    if(self.params.mode == "filter_only"): return
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    self.find_peaks_params.max_number_of_peaks=self.model.get_number_of_atoms()
    assert self.params.primary_map_type is not None
    self._peaks = find_peaks.manager(
      xray_structure = self.fmodel.xray_structure,
      map_data       = self._maps.difference_map, # diff-map
      map_cutoff     = self.params.primary_map_cutoff,
      params         = self.find_peaks_params,
      log            = null_out()).peaks_mapped()

  def _write_pdb_file(self, file_name="tmp.pdb", sites_frac=None):
    if sites_frac is not None:
      fmt = "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f  1.00  0.00           O"
      uc = self.fmodel.xray_structure.crystal_symmetry().unit_cell()
      with open(file_name, "w") as fo:
        for i, sf in enumerate(sites_frac):
          sc = uc.orthogonalize(sf)
          print(fmt%(i,i,sc[0],sc[1],sc[2]), file=fo)
    else:
      with open(file_name, "w") as fo:
        fo.write(self.model.model_as_pdb())

  def _add_new_solvent(self, conformer_indices=None):
    """
    A series of very specialzied manipulations on the model to add new solvent
    without re-creating the new model.
    """
    if(self._peaks is None): return
    if(self._peaks.sites.size()==0): return
    if(self.params.mode == "filter_only"): return

    self.new_solvent_selection = add_solvent_to_model_inplace(
      sites             = self._peaks.sites,
      model             = self.model,
      params            = self.params,
      conformer_indices = None)
    if self.params.include_altlocs:
      self.model = fix_altlocs_and_filter(model=self.model, fix_only=True)
      ss = self.model.solvent_selection()
      ms = self.model.select(ss)
      self.model = self.model.select(~ss)
      self.new_solvent_selection = add_solvent_to_model_inplace(
        sites             = ms.get_sites_frac(),
        model             = self.model,
        params            = self.params,
        conformer_indices = ms.get_hierarchy().get_conformer_indices())

  def refine_oat(self):
    if self.new_solvent_selection is None: return
    if(self.params.refine_oat and self.new_solvent_selection.count(True)>0):
      from phenix.programs import oat
      from cctbx import adptbx
      atoms = self.model.get_hierarchy().atoms()
      scatterers = self.fmodel.xray_structure.scatterers()
      for i, sel in enumerate(self.new_solvent_selection):
        if not sel: continue
        r_start = self.fmodel.r_work()
        scatterers[i].occupancy=0
        scatterers[i].u_iso=0
        self.fmodel.update_xray_structure(update_f_calc=True)
        r_omit = self.fmodel.r_work()
        fmodel_dc = self.fmodel.deep_copy()
        oo = oat.loop(
          fmodel  = fmodel_dc,
          site    = atoms[i].xyz,
          label   = "O",
          qs = flex.double([0.3, 0.6, 0.9]),
          bs = flex.double([10, 30, 60])
          )
        scatterers[i].occupancy = oo.o_best
        scatterers[i].u_iso     = adptbx.b_as_u(oo.b_best)
        self.fmodel.update_xray_structure(update_f_calc=True)
        r_final = self.fmodel.r_work()
        #print(atoms[i].quote(), "%8.6f %8.6f %8.6f %8.6f"%(
        #  r_start, r_omit, r_final, oo.rw_best), "|", oo.o_best, oo.b_best)
      self.model.adopt_xray_structure(
        xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), OAT, final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)
      #
    if([self.params.refine_adp , self.params.refine_occupancies].count(True)>0 and
       self.new_solvent_selection.count(True)>0):
      from mmtbx.refinement import wrappers
      if self.params.new_solvent == "anisotropic":
        raise RuntimeError("Not implemented: new_solvent=anisotropic")
      o = wrappers.unrestrained_qbr_fsr(
        fmodel     = self.fmodel,
        model      = self.model,
        refine_xyz = False,
        refine_q   = self.params.refine_occupancies,
        refine_b   = self.params.refine_adp,
        selection  = self.new_solvent_selection,
        q_min      = 0.004,
        b_max      = 60,
        b_min      = self.params.b_iso_min,
        log        = self.log)
      self.model.adopt_xray_structure(
          xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), MIN, final r_work=%6.4f r_free=%6.4f"%(
          self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def _correct_drifted_waters(self):
    if(self.params.mode != "filter_only"): return
    if(not self.params.correct_drifted_waters): return
    sol_sel = self.model.solvent_selection()
    hd_sel  = self.model.get_hd_selection()
    hd_sol = sol_sel & hd_sel
    if(hd_sol.count(True)>0): return
    map_cutoff = self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold/2
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    peaks = find_peaks.manager(
      xray_structure = self.fmodel.xray_structure,
      map_data       = self._maps.data_map,
      map_cutoff     = map_cutoff,
      params         = find_peaks_params_drifted,
      log            = null_out()).peaks_mapped()
    if(peaks is not None and self.fmodel.r_work() > 0.01):
      sites_frac, heights = peaks.sites, peaks.heights
      model_sites_frac = self.model.get_xray_structure().sites_frac()
      solvent_selection = self.model.solvent_selection()
      mmtbx.utils.correct_drifted_waters(
        sites_frac_all   = model_sites_frac,
        sites_frac_peaks = sites_frac,
        water_selection  = solvent_selection,
        unit_cell        = self.model.get_xray_structure().unit_cell())
      self.model.get_xray_structure().set_sites_frac(
        sites_frac = model_sites_frac)
