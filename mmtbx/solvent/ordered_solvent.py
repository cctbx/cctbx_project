from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math,sys
import copy
from libtbx import adopt_init_args
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.phil
from mmtbx import find_peaks
import mmtbx.utils
from cctbx import maptbx
from cctbx import crystal
from libtbx.test_utils import approx_equal
from six.moves import range
from libtbx.utils import user_plus_sys_time
from libtbx.utils import null_out
from mmtbx.solvent import map_to_water
from libtbx import group_args
import libtbx.log
from libtbx import Auto
from cctbx import geometry_restraints
from iotbx.pdb.mmcif import format_pdb_atom_name
from iotbx.pdb import water_altlocs
from cctbx import xray
from scitbx import matrix

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
    .short_caption = Min distance between alternative locations
    .help = Min distance between water and atom in another alternative \
            conformation. Two alternative water sites must also be at \
            least 0.7*d_min apart.
  dist_min_h_nonpolar = 2.0
    .type = float
    .short_caption = Min distance between water and H bonded to carbon
    .help = Closer contacts are clashes (O and H radii minus 0.4 A overlap)
    .expert_level = 1
  dist_min_h_polar = 1.5
    .type = float
    .short_caption = Min distance between water and polar H (on N, O, S)
    .help = Polar H can make H-bonds with water, so shorter than for non-polar H
    .expert_level = 1
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

# ensemble refinement includes output_params_str and adp_occ_params_str as
# they are; phenix.refine water building uses the two variants below: no
# output chain ID (water takes the ID of the nearest chain) and Auto choice
# for ADP of new water.
def _without_chain_id(s):
  lines = s.splitlines()
  result = []
  skip = 0
  for l in lines:
    if(l.strip().startswith("output_chain_id")):
      skip = 2 # its .type and .input_size lines
      continue
    if(skip > 0):
      skip -= 1
      continue
    result.append(l)
  return "\n".join(result)+"\n"

output_params_no_chain_str = _without_chain_id(output_params_str)

adp_occ_params_auto_str = adp_occ_params_str.replace(
  """  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
""",
  """  new_solvent = *Auto isotropic anisotropic
    .type = choice
    .help = ADP of added water. Auto: anisotropic if resolution is better \
            than anisotropic_d_min and the model is anisotropic, isotropic \
            otherwise.
  anisotropic_d_min = 1.2
    .type = float
    .help = Resolution limit for new_solvent=Auto: added water is \
            anisotropic if data resolution is better than this.
    .short_caption = Resolution for anisotropic new water
    .expert_level = 1
""")
assert adp_occ_params_auto_str != adp_occ_params_str

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
  keep_input_water = False
    .type = bool
    .help = Do not remove water that comes with input model, and do not \
            change its coordinates, occupancy, ADP and altloc in water \
            update. Its chain ID and residue number may change, and it may \
            be moved to a symmetry copy next to its chain.
    .short_caption = Keep input water
  coordination_restraints
    .short_caption = Water coordination restraints
    .style = auto_align box
  {
    enabled = True
      .type = bool
      .help = Restrain distance between each water and the atom that \
              coordinates it to the current value. This also excludes the \
              pair from nonbonded repulsion.
    sigma = 0.5
      .type = float
      .help = Sigma of water coordination restraint
  }
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
""" % (output_params_no_chain_str, h_bond_params_str,
       adp_occ_params_auto_str)

def master_params():
  return iotbx.phil.parse(master_params_str)

def new_solvent_sites_as_hierarchy_chain(
      model,
      sites,
      params,
      chain_id="S",
      altlocs=None,
      occupancies=None):
  """
  Chain of new water molecules (O atoms with isotropic B) at sites
  (fractional). altlocs, occupancies: per site, optional.
  """
  uc = model.crystal_symmetry().unit_cell()
  ort = uc.orthogonalize
  b_iso = flex.mean(model.get_xray_structure(
    ).extract_u_iso_or_u_equiv()) * adptbx.u_as_b(1) * 1.5
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
  new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
  for i_seq, site_frac in enumerate(sites):
    water_i_seq = size+i_seq+1 # same as resseq below
    occ = 1
    if(occupancies is not None): occ = occupancies[i_seq]
    new_atom = (iotbx.pdb.hierarchy.atom()
      .set_serial(new_serial=iotbx.pdb.hy36encode(width=5, value=water_i_seq))
      .set_name(new_name=atom_name)
      .set_xyz(new_xyz=ort(site_frac))
      .set_occ(new_occ=occ)
      .set_b(new_b=b_iso)
      .set_element(element)
      .set_charge("")
      .set_hetero(new_hetero=True))
    altloc=""
    if(altlocs is not None): altloc = altlocs[i_seq]
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
      water_residue_groups=None,
      chain_id=None,
      conformer_indices=None,
      nonbonded_types=None,
      nonbonded_charges=None,
      solvent_xray_structure=None,
      altlocs=None,
      occupancies=None,
      keep_resseq=False,
      update_plain_pair_sym_table=True):
  """
  A series of very specialzied manipulations on the model to add new solvent
  without re-creating the new model.

  New water is given either as sites (fractional coordinates; with altlocs
  and occupancies per site, optional) or as water_residue_groups. chain_id:
  ID of the chain it is added to ("S" if not given). keep_resseq: residue
  numbers of water_residue_groups are final. conformer_indices: object with
  attribute conformer_indices (for atoms being added); if not given they
  follow from altlocs of the atoms. update_plain_pair_sym_table=False: the
  caller does it (it is expensive, and once is enough if several chains are
  added one after another).
  """
  assert [sites, water_residue_groups].count(None)==1
  if chain_id is None: chain_id = "S"
  # altloc -> conformer index as used in restraints (before adding atoms)
  altloc_to_conformer_index = {}
  if(conformer_indices is None and model.get_restraints_manager() is not None
     and model.get_restraints_manager().geometry.conformer_indices is not None):
    cis = model.get_restraints_manager().geometry.conformer_indices
    for ag in model.get_hierarchy().atom_groups():
      if(ag.altloc in altloc_to_conformer_index or ag.atoms_size()==0):
        continue
      altloc_to_conformer_index[ag.altloc] = cis[ag.atoms()[0].i_seq]
  # Get new solvent as pdb hierarchy chain
  if sites is None:
    new_solvent_chain = create_water_chain(
      water_residue_groups=water_residue_groups,
      chain_id=chain_id,
      start_resseq=1,
      keep_resseq=keep_resseq)
  else:
    new_solvent_chain = new_solvent_sites_as_hierarchy_chain(
      sites       = sites,
      model       = model,
      params      = params,
      chain_id    = chain_id,
      altlocs     = altlocs,
      occupancies = occupancies)
  if new_solvent_chain is None: return
  new_altlocs = [a.parent().altloc for a in new_solvent_chain.atoms()]
  # Get new solvent selection
  n_new_solvent      = new_solvent_chain.atoms().size()
  initial_model_size = model.size()
  new_solvent_selection = flex.bool(initial_model_size, False)
  new_solvent_selection.extend(flex.bool(n_new_solvent, True))
  # Append new chain to the model
  #model.get_hierarchy().only_model().append_chain(chain=new_solvent_chain)
  add_chain_to_hierarchy(
    hierarchy = model.get_hierarchy(), new_chain = new_solvent_chain,
    keep_resseq = keep_resseq)
  model._update_atom_selection_cache()
  model.get_hierarchy().atoms().reset_i_seq()
  # Existing atoms keep their old serials across model.select() (the filter
  # steps), so serials derived from model.size() can collide with them.
  # Renumber everything to keep serials unique in PDB/mmCIF output.
  model.get_hierarchy().atoms_reset_serial()
  model.unset_processed_pdb_file()
  # Force-update xray_structure
  # This is done THIS WAY to keep scattering_table
  if solvent_xray_structure is None:
    solvent_xray_structure = new_solvent_chain.as_new_hierarchy(
      ).extract_xray_structure(crystal_symmetry=model.crystal_symmetry())
  assert solvent_xray_structure.scatterers().size() == n_new_solvent
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
      if(rfs.sites_torsion_angles is not None):
        rfs.sites_torsion_angles.append(True)
      # ADP
      if(rfs.individual_adp): # H?
        iso = rfs.adp_individual_iso
        ani = rfs.adp_individual_aniso

        # if water happened to be anisotropic already (= not really new!)
        if scatterers[i].flags.use_u_aniso():
          if iso is not None: iso.append(False)
          if ani is not None: ani.append(True)
          continue

        # isotropic unless the atom comes as anisotropic
        if iso is not None: iso.append(True)
        if ani is not None: ani.append(False)
  #
  if nonbonded_types is None:
    nonbonded_types = flex.std_string()
    for atom in new_solvent_chain.atoms():
      if(atom.element_is_hydrogen()): nonbonded_types.append('HOH2')
      else:                           nonbonded_types.append('OH2')
  if nonbonded_charges is None:
    nonbonded_charges = flex.int(n_new_solvent, 0)
  #
  # Finally, add new atoms to restraints manager
  #
  geometry = None
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
        conformer_indices = flex.size_t()
        for altloc in new_altlocs:
          conformer_indices.append(altloc_to_conformer_index.get(altloc, 0))
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
      nonbonded_charges          = nonbonded_charges)
    model.restraints_manager = mmtbx.restraints.manager(
      geometry              = geometry,
      cartesian_ncs_manager = rm.cartesian_ncs_manager,
      normalization         = rm.normalization)
    c_ncs_m = model.get_cartesian_NCS_manager()
    if(c_ncs_m is not None):
      c_ncs_m.register_additional_isolated_sites(
        number=n_new_solvent)
    if(update_plain_pair_sym_table):
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
      # new atoms are at the end: in the last chain
      for c in m.chains()[-1:]:
        for conf in c.conformers():
          for r in conf.residues():
            if(not get_class(name=r.resname) == "common_water"): continue
            i_seqs = r.atoms().extract_i_seq()
            if i_seqs.size() not in [2, 3]: continue
            # only new water: existing one already has restraints
            if flex.min(i_seqs) < initial_model_size: continue
            io = None
            ihs = []
            for atom in r.atoms():
              if(atom.element_is_hydrogen()): ihs.append(atom.i_seq)
              else: io = atom.i_seq
            #
            if io is not None and len(ihs) in [1, 2]:
              for ih in ihs:
                bproxy = geometry_restraints.bond_simple_proxy(
                  i_seqs         = (ih, io),
                  distance_ideal = target_bond,
                  weight         = 1./(0.02**2),
                  origin_id      = 0)
                bproxies.append(bproxy)
              if(len(ihs) == 2):
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
  if model.riding_h_manager is not None and geometry is not None:
    new_riding_h_manager = model.riding_h_manager.update(
      pdb_hierarchy       = model.get_hierarchy(),
      geometry_restraints = geometry,
      n_new_atoms         = n_new_solvent)
    model.riding_h_manager = new_riding_h_manager
  return new_solvent_selection

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

def ideal_water_h_sites(o_site, h_sites, distances, angle):
  """
  Sites of two H of water with given O-H distances and H-O-H angle, oriented
  as the input H sites are: in the plane of O and the two input sites,
  symmetrically about the bisector of their directions from O.

  Parameters
  ----------
  o_site : sequence of three floats
  h_sites : sequence of two sites
    Approximate H sites, difference map peaks for example. Their directions
    from O must be neither the same nor opposite.
  distances : sequence of two floats
    O-H distances.
  angle : float
    H-O-H angle, degrees.

  Returns
  -------
  tuple
    Two H sites, in the order of h_sites.
  """
  o = matrix.col(o_site)
  u1 = (matrix.col(h_sites[0]) - o).normalize()
  u2 = (matrix.col(h_sites[1]) - o).normalize()
  bisector = (u1 + u2).normalize()
  across   = (u1 - u2).normalize()
  half = math.radians(angle/2)
  c, s = math.cos(half), math.sin(half)
  return (o + distances[0]*(c*bisector + s*across)).elems, \
         (o + distances[1]*(c*bisector - s*across)).elems

def ideal_water_h_geometry(model, water_i_seqs):
  """
  Ideal O-H distances and H-O-H angles of water molecules, from restraints.

  Parameters
  ----------
  model : mmtbx.model.manager
  water_i_seqs : list
    (i_seq of O, i_seq of H, i_seq of H) for each water.

  Returns
  -------
  list
    For each water: ((distance, distance), angle), or None if the model has
    no restraints or they do not define the geometry of this water.
  """
  result = [None]*len(water_i_seqs)
  rm = model.get_restraints_manager()
  if(rm is None): return result
  table = rm.geometry.bond_params_table
  angles = {}
  o_i_seqs = set([i_seqs[0] for i_seqs in water_i_seqs])
  for proxy in rm.geometry.angle_proxies:
    if(proxy.i_seqs[1] in o_i_seqs):
      angles[(proxy.i_seqs[1], min(proxy.i_seqs[0], proxy.i_seqs[2]),
        max(proxy.i_seqs[0], proxy.i_seqs[2]))] = proxy.angle_ideal
  for i, (o, h1, h2) in enumerate(water_i_seqs):
    b1 = table.lookup(o, h1)
    b2 = table.lookup(o, h2)
    angle = angles.get((o, min(h1, h2), max(h1, h2)))
    if(b1 is None or b2 is None or angle is None): continue
    result[i] = ((b1.distance_ideal, b2.distance_ideal), angle)
  return result

def get_detached_water_residue_group(residue_group, h_sites_cart):
    """
    Copy of water residue_group with H atoms added.

    Parameters
    ----------
    residue_group : iotbx.pdb.hierarchy.residue_group
    h_sites_cart : dict
      Index of atom group of residue_group -> sites of H atoms to add to it.

    Returns
    -------
    tuple
      Detached residue group and, for each of its atoms, i_seq of the atom
      it is a copy of or -(i_seq+1) for a new H of water oxygen i_seq.
    """
    result = residue_group.detached_copy()
    sources = []
    for i_ag, ag in enumerate(residue_group.atom_groups()):
      new_ag = result.atom_groups()[i_ag]
      o_atom, new_o_atom = None, None
      for atom, new_atom in zip(ag.atoms(), new_ag.atoms()):
        sources.append(atom.i_seq)
        if(atom.element.strip().upper()=="O"):
          o_atom, new_o_atom = atom, new_atom
      for i, h_xyz in enumerate(h_sites_cart.get(i_ag, [])):
        h_atom = new_o_atom.detached_copy()
        h_atom.name = " H%d "%(i+1)
        h_atom.element = "H"
        h_atom.xyz = h_xyz
        h_atom.uij = (-1, -1, -1, -1, -1, -1)
        new_ag.append_atom(h_atom)
        sources.append(-(o_atom.i_seq+1))
    return result, sources

def create_water_chain(water_residue_groups, chain_id, start_resseq=1,
                       keep_resseq=False):
  """
  Creates and returns a new iotbx.pdb.hierarchy.chain populated with
  the provided detached water residue groups. Residue groups are numbered
  from start_resseq unless keep_resseq is True.
  """
  new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
  current_resseq = start_resseq
  for rg in water_residue_groups:
    new_rg = rg.detached_copy()
    # Assign a sequential sequence number to the water
    if(not keep_resseq):
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

def add_chain_to_hierarchy(hierarchy, new_chain, keep_resseq=False):
  """
  Adds the residues of new_chain to the hierarchy. If the last chain of the
  model has the same ID, the residues are appended to that chain object so
  that all waters of one chain stay in one chain object (one label_asym_id
  in mmCIF). Otherwise new_chain is appended as a separate chain. Either
  way the new atoms end up at the end of the atom sequence, which the
  caller relies on. The new residues are renumbered to prevent resseq
  collisions without disturbing existing data structures, unless keep_resseq
  is True.
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
    if(keep_resseq): break
    incoming_resid = rg.resid()
    if(previous_incoming_resid is not None and
       incoming_resid != previous_incoming_resid):
      current_resseq += 1
    rg.resseq = iotbx.pdb.resseq_encode(value=current_resseq)
    rg.icode = " "
    previous_incoming_resid = incoming_resid
  if hierarchy.models_size() == 0:
    hierarchy.append_model(iotbx.pdb.hierarchy.model())
  chains = hierarchy.models()[0].chains()
  # keep_resseq: water is being regrouped by chain ID, it stays a chain
  # object of its own (not part of the chain of the same ID before it)
  if not keep_resseq and len(chains) > 0 and chains[-1].id == new_chain.id:
    # Merge into the existing chain object. Only the last chain qualifies:
    # appending to an earlier chain would put the new atoms in the middle of
    # the atom sequence.
    target_chain = chains[-1]
    for rg in new_chain.residue_groups():
      new_rg = rg.detached_copy()
      new_rg.link_to_previous = True # no BREAK record between waters
      target_chain.append_residue_group(residue_group=new_rg)
  else:
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

get_class = iotbx.pdb.common_residue_names_get_class

def is_water(atom_group):
  return get_class(name=atom_group.resname) == "common_water"

def is_water_residue_group(residue_group):
  for ag in residue_group.atom_groups():
    if(not is_water(ag)): return False
  return True

def model_altlocs(model):
  """
  Altloc identifiers that water is allowed to have: those of non-water atoms
  of the model, sorted. Empty list: water is built in single conformation.
  """
  result = set()
  for ag in model.get_hierarchy().atom_groups():
    if(is_water(ag)): continue
    altloc = ag.altloc.strip()
    if(altloc != ""): result.add(altloc)
  return sorted(result)

def water_strength(occ, b_iso):
  """
  Number proportional to the peak density of an atom: used to decide which
  of two waters that cannot coexist stays.
  """
  return occ / (max(b_iso, 0.) + 10.)**1.5

def get_water_neighbors(model, distance_cutoff):
  """
  Neighbors (accounting for symmetry) of non-H atoms of water.
  Returns dictionary: i_seq of water atom -> list of (j_seq, distance) for
  atoms within distance_cutoff, excluding atoms of the same residue group
  (own H and alternative conformations of the same water). Symmetry copies
  of the water atom itself are included as (i_seq, distance): a water close
  to a symmetry element clashes with its own copy.
  One neighbor search for all atoms: much faster than one search per water.
  """
  atoms = model.get_hierarchy().atoms()
  water_i_seqs = []
  for ag in model.get_hierarchy().atom_groups():
    if(not is_water(ag)): continue
    for a in ag.atoms():
      if(not a.element_is_hydrogen()): water_i_seqs.append(a.i_seq)
  result = {}
  if(len(water_i_seqs) == 0): return result
  closest = {}
  for i in water_i_seqs: closest[i] = {}
  for i, j, rt_mx, d in crystal.pairs_within(model.get_sites_cart(),
      distance_cutoff, model.crystal_symmetry(), i_seqs=water_i_seqs):
    if(i != j and atoms[j].parent().parent().memory_id() ==
       atoms[i].parent().parent().memory_id()): continue
    if(j not in closest[i] or d < closest[i][j]): closest[i][j] = d
  for i in water_i_seqs:
    result[i] = sorted(closest[i].items())
  return result

def is_polar_h(atom):
  """
  True if H atom is bonded (closest non-H atom of the same residue group within
  1.3 A) to N, O or S, i.e. can make H-bonds.
  """
  best = None
  for a in atom.parent().parent().atoms():
    if(a.element_is_hydrogen()): continue
    d = a.distance(atom)
    if(d < 1.3 and (best is None or d < best[0])):
      best = (d, a.element.strip().upper())
  return best is not None and best[1] in ["N", "O", "S"]

class exclusion_distances(object):
  """
  Smallest distance at which a water can be present together with an atom.

  Water against water, against an atom in alternative conformation and
  against an atom with partial occupancy: iotbx.pdb.water_altlocs cutoff
  (such a pair is a pair of alternates, water_altlocs would label it).
  Metals: dist_min whatever their occupancy or altloc (coordination).
  H: dist_min_h_polar or dist_min_h_nonpolar. Everything else: dist_min.
  """
  def __init__(self, dist_min, dist_min_h_nonpolar, dist_min_h_polar):
    adopt_init_args(self, locals())
    self.alt = max(water_altlocs.DISTANCE_CUTOFF, dist_min)
    self._h_is_polar = {}

  def max(self):
    return max(self.alt, self.dist_min, self.dist_min_h_nonpolar,
      self.dist_min_h_polar)

  def __call__(self, atom):
    if(atom.element_is_hydrogen()):
      is_polar = self._h_is_polar.get(atom.i_seq)
      if(is_polar is None):
        is_polar = is_polar_h(atom)
        self._h_is_polar[atom.i_seq] = is_polar
      if(is_polar): return self.dist_min_h_polar
      return self.dist_min_h_nonpolar
    if(water_altlocs.is_metal(atom)): return self.dist_min
    ag = atom.parent()
    if(is_water(ag) or ag.altloc.strip() != "" or atom.occ < 1.):
      return self.alt
    return self.dist_min

def exclusion_distance(atom, params):
  """
  Smallest distance at which a water can be present together with atom, see
  exclusion_distances.
  """
  return exclusion_distances(
    dist_min            = params.dist_min,
    dist_min_h_nonpolar = params.dist_min_h_nonpolar,
    dist_min_h_polar    = params.dist_min_h_polar)(atom)

def can_coexist(altloc_i, altloc_j):
  """
  Atoms are present at the same time unless they are in different
  alternative conformations.
  """
  return altloc_i == "" or altloc_j == "" or altloc_i == altloc_j

def min_water_separation(params, d_min, anisotropic=False):
  """
  Two alternative water sites closer than this are one water: resolution
  does not support two sites closer than 0.7*d_min, and if water is
  anisotropic two sites closer than 1 A are one anisotropic water.
  """
  result = params.dist_min_altloc
  if(d_min is not None): result = max(result, 0.7*d_min)
  if(anisotropic): result = max(result, 1.0)
  return result

def assign_conformers(model, sites_cart, heights, params, d_min=None,
                      fixed_water_selection=None, anisotropic=False):
  """
  Decide which candidate water sites (map peaks) are added to the model and
  in which alternative conformation.

  Altloc identifiers come from non-water atoms of the model; without them
  water is built in single conformation. A site that cannot be present
  together with an atom (see exclusion_distances) is rejected if that atom
  has no altloc, otherwise it takes an altloc that none of the atoms it
  excludes has: preferably that of an atom in alternative conformation it
  can H-bond with, else the most populated one of the residue it excludes,
  else the first available. Two sites that exclude each other get different
  altlocs if at least min_water_separation apart; otherwise the weaker one
  is rejected. Sites are processed from strongest to weakest.

  Parameters
  ----------
  model : mmtbx.model.manager
  sites_cart : flex.vec3_double
    Candidate sites.
  heights : flex.double
    Peak heights of the sites.
  params : ordered_solvent parameters
  d_min : float, optional
    Data resolution.
  fixed_water_selection : flex.bool, optional
    Atoms of water that must not be relabelled.
  anisotropic : bool, optional
    New water is going to be anisotropic, see min_water_separation.

  Returns
  -------
  group_args
    keep (flex.bool over sites), altlocs (list of str over sites),
    occupancies (flex.double over sites), relabel (dict: i_seq of non-H atom
    of existing water -> its new altloc).
  """
  n_sites = sites_cart.size()
  assert heights.size() == n_sites
  keep = flex.bool(n_sites, True)
  altlocs = [""]*n_sites
  occupancies = flex.double(n_sites, params.occupancy)
  result = group_args(keep=keep, altlocs=altlocs, occupancies=occupancies,
    relabel={})
  if(n_sites == 0): return result
  letters = model_altlocs(model)
  atoms = model.get_hierarchy().atoms()
  n_atoms = atoms.size()
  if(fixed_water_selection is None):
    fixed_water_selection = flex.bool(n_atoms, False)
  interaction_selection = model.selection(
    map_to_water.selection_string_interaction)
  excl = exclusion_distances(
    dist_min            = params.dist_min,
    dist_min_h_nonpolar = params.dist_min_h_nonpolar,
    dist_min_h_polar    = params.dist_min_h_polar)
  separation = min_water_separation(params=params, d_min=d_min,
    anisotropic=anisotropic)
  all_sites = model.get_sites_cart().deep_copy()
  all_sites.extend(sites_cart)
  site_i_seqs = list(range(n_atoms, n_atoms+n_sites))
  # Closest distance for each (site, atom or site) pair
  contacts = [dict() for i in range(n_sites)]
  for i, j, rt_mx, d in crystal.pairs_within(all_sites,
      max(excl.max(), params.dist_max), model.crystal_symmetry(),
      i_seqs=site_i_seqs):
    c = contacts[i-n_atoms]
    if(j not in c or d < c[j]): c[j] = d
  #
  def reject(k):
    keep[k] = False
    altlocs[k] = ""
  def occupancy_of(atom_groups, letter):
    occ = flex.double()
    for ag in atom_groups:
      if(ag.altloc.strip() != letter): continue
      for a in ag.atoms():
        if(not a.element_is_hydrogen()): occ.append(a.occ)
    if(occ.size() == 0): return None
    return flex.mean(occ)
  # Model atoms each site excludes / can interact with
  site_model_excl  = [[] for i in range(n_sites)] # atoms with altloc
  site_water_excl  = [[] for i in range(n_sites)] # blank water, can relabel
  site_votes       = [[] for i in range(n_sites)] # altlocs of H-bond partners
  site_site_excl   = [[] for i in range(n_sites)] # (other site, distance)
  for k in range(n_sites):
    for j, d in contacts[k].items():
      if(j == k+n_atoms): # own symmetry copy
        if(d < excl.alt): reject(k)
        continue
      if(j >= n_atoms):
        if(d < excl.alt): site_site_excl[k].append((j-n_atoms, d))
        continue
      aj = atoms[j]
      agj = aj.parent()
      altloc_j = agj.altloc.strip()
      if(d >= excl(aj)):
        if(altloc_j != "" and d <= params.dist_max and
           interaction_selection[j] and not aj.element_is_hydrogen()):
          site_votes[k].append(altloc_j)
        continue
      if(d < params.dist_min_altloc and not aj.element_is_hydrogen()):
        reject(k)
      elif(altloc_j != ""):
        site_model_excl[k].append(j)
      elif(is_water(agj) and not aj.element_is_hydrogen() and
           not fixed_water_selection[j] and
           d >= separation and len(letters) > 0):
        site_water_excl[k].append(j)
      else:
        reject(k)
  # Strongest first
  relabel = result.relabel
  def letters_excluded_by_model(k):
    r = set()
    for j in site_model_excl[k]: r.add(atoms[j].parent().altloc.strip())
    for j in site_water_excl[k]:
      if(j in relabel): r.add(relabel[j])
    return r
  def choose(k, excluded):
    candidates = [l for l in letters if l not in excluded]
    if(len(candidates) == 0): return None, None
    rgs = {}
    for j in site_model_excl[k]:
      rg = atoms[j].parent().parent()
      rgs[rg.memory_id()] = rg
    def residue_occupancy(l):
      r = None
      for rg in rgs.values():
        occ = occupancy_of(rg.atom_groups(), l)
        if(occ is not None and (r is None or occ > r)): r = occ
      return r
    # altloc of H-bond partners
    best, best_count = None, 0
    for l in candidates:
      n = site_votes[k].count(l)
      if(n > best_count): best, best_count = l, n
    # most populated conformer of the residues the site excludes
    if(best is None):
      best_occ = None
      for l in candidates:
        occ = residue_occupancy(l)
        if(occ is not None and (best_occ is None or occ > best_occ)):
          best, best_occ = l, occ
    if(best is None): return candidates[0], None
    occ = residue_occupancy(best)
    if(occ is None):
      # occupancy of the H-bond partner
      for j, d in contacts[k].items():
        if(j >= n_atoms): continue
        aj = atoms[j]
        if(aj.parent().altloc.strip() == best and
           not aj.element_is_hydrogen() and d >= excl(aj)):
          occ = aj.occ
          break
    return best, occ
  order = flex.sort_permutation(heights, reverse=True)
  accepted = set()
  for k in order:
    if(not keep[k]): continue
    excluded = letters_excluded_by_model(k)
    blank_sites  = [] # accepted blank sites this one excludes
    too_close = False
    for q, d in site_site_excl[k]:
      if(not keep[q] or q not in accepted): continue
      if(d < separation or len(letters) == 0):
        too_close = True
        break
      if(altlocs[q] == ""): blank_sites.append(q)
      else: excluded.add(altlocs[q])
    if(too_close):
      reject(k)
      continue
    blank_waters = [j for j in site_water_excl[k] if j not in relabel]
    if(len(excluded) == 0 and len(blank_sites) == 0 and
       len(blank_waters) == 0):
      accepted.add(k)
      continue
    # Partners that are blank so far need altlocs too. If the site is
    # restricted in its choice, it chooses first and partners take any other
    # altloc; otherwise partners (they are stronger) choose first.
    trial = {}
    letter, occ = None, None
    site_first = len(excluded) > 0
    if(site_first):
      letter, occ = choose(k, excluded)
      if(letter is None):
        reject(k)
        continue
    for q in blank_sites:
      taken = set([letter]) if site_first else set(trial.values())
      l, occ_q = choose(q, taken)
      if(l is None): break
      trial[("site", q)] = l
    for j in blank_waters:
      taken = set([letter]) if site_first else set(trial.values())
      candidates = [l for l in letters if l not in taken]
      if(len(candidates) == 0): break
      trial[("water", j)] = candidates[0]
    if(len(trial) != len(blank_sites)+len(blank_waters)):
      reject(k)
      continue
    if(not site_first):
      letter, occ = choose(k, set(trial.values()))
      if(letter is None):
        reject(k)
        continue
    altlocs[k] = letter
    if(occ is None): occ = 0.5
    occupancies[k] = min(occ, params.occupancy)
    for key, l in trial.items():
      if(key[0] == "site"):
        altlocs[key[1]] = l
        occupancies[key[1]] = 0.5
      else:
        relabel[key[1]] = l
    accepted.add(k)
  return result

def conformer_index_by_altloc(model):
  """
  Altloc -> conformer index as used in geometry restraints of the model;
  None if there are no restraints or no conformer indices.
  """
  rm = model.get_restraints_manager()
  if(rm is None or rm.geometry.conformer_indices is None): return None
  conformer_indices = rm.geometry.conformer_indices
  result = {"": 0}
  for ag in model.get_hierarchy().atom_groups():
    altloc = ag.altloc.strip()
    if(altloc in result or ag.atoms_size() == 0): continue
    result[altloc] = conformer_indices[ag.atoms()[0].i_seq]
  return result

def set_altloc(model, atom_group, altloc, index_by_altloc=None):
  """
  Change altloc of atom_group of model in-place, keeping conformer indices
  of geometry restraints consistent. altloc must be blank or one that atoms
  of the model already have. index_by_altloc: result of
  conformer_index_by_altloc for the model as it was before any change of
  altlocs (to change many without looking at all atoms each time).
  """
  if(index_by_altloc is None):
    index_by_altloc = conformer_index_by_altloc(model)
  atom_group.altloc = altloc
  if(index_by_altloc is None): return
  g = model.get_restraints_manager().geometry
  index = index_by_altloc[altloc.strip()]
  for a in atom_group.atoms():
    g.conformer_indices[a.i_seq] = index
  g.reset_internals()

def distance_filter_selection(model, params=None, fixed_water_selection=None,
                       dist_min=None, dist_max=None, dist_min_h_nonpolar=None,
                       dist_min_h_polar=None):
  """
  Selection of atoms that stay after removing water that:
    - cannot be present together with another atom (see exclusion_distances)
      unless the two are in different alternative conformations; of two
      waters the weaker one is removed;
    - is closer than that to its own symmetry copy;
    - has no interacting atom (N, O, S, metals; not H) at distance up to
      dist_max that can be present at the same time.
  Atoms of the same residue group (own H) are ignored. Water selected by
  fixed_water_selection is never removed. Distances are taken from params
  unless given explicitly.
  """
  if(params is None): params = master_params().extract()
  if(dist_min is None): dist_min = params.dist_min
  if(dist_max is None): dist_max = params.dist_max
  if(dist_min_h_nonpolar is None):
    dist_min_h_nonpolar = params.dist_min_h_nonpolar
  if(dist_min_h_polar is None): dist_min_h_polar = params.dist_min_h_polar
  if(fixed_water_selection is None):
    fixed_water_selection = flex.bool(model.size(), False)
  interaction_selection = model.selection(
    map_to_water.selection_string_interaction)
  excl = exclusion_distances(
    dist_min            = dist_min,
    dist_min_h_nonpolar = dist_min_h_nonpolar,
    dist_min_h_polar    = dist_min_h_polar)
  only_model = model.get_hierarchy().only_model()
  atoms = only_model.atoms()
  neighbors = get_water_neighbors(model=model,
    distance_cutoff=max(dist_max, excl.max()))
  xrs = model.get_xray_structure()
  b_isos = xrs.extract_u_iso_or_u_equiv()*adptbx.u_as_b(1)
  occs = xrs.scatterers().extract_occupancies()
  def strength(i):
    if(fixed_water_selection[i]): return 1.e+9
    return water_strength(occ=occs[i], b_iso=b_isos[i])
  water_o = [] # (atom group, atom)
  for agi in only_model.atom_groups():
    if(not is_water(agi)): continue
    for ai in agi.atoms():
      if(not ai.element_is_hydrogen()): water_o.append((agi, ai))
  removed = set() # i_seqs of removed water non-H atoms
  # Cannot be present together with something that stays
  for agi, ai in sorted(water_o, key=lambda x: -strength(x[1].i_seq)):
    i = ai.i_seq
    if(fixed_water_selection[i]): continue
    altloc_i = agi.altloc.strip()
    for j, d in neighbors[i]:
      if(j == i):
        if(d < excl.alt):
          removed.add(i)
          break
        continue
      if(j in removed): continue
      aj = atoms[j]
      if(d >= excl(aj)): continue
      agj = aj.parent()
      if(not can_coexist(altloc_i, agj.altloc.strip())): continue
      if(is_water(agj) and not aj.element_is_hydrogen()):
        sj, si = strength(j), strength(i)
        if(sj < si or (sj == si and j > i)): continue # the other one goes
      removed.add(i)
      break
  # No interacting atom: repeat, water can lose its only partner
  while True:
    n_removed = len(removed)
    for agi, ai in water_o:
      i = ai.i_seq
      if(i in removed or fixed_water_selection[i]): continue
      altloc_i = agi.altloc.strip()
      found = False
      for j, d in neighbors[i]:
        if(d > dist_max or j in removed): continue
        if(not interaction_selection[j]): continue
        aj = atoms[j]
        if(aj.element_is_hydrogen()): continue
        if(d < excl(aj) and j != i): continue
        if(can_coexist(altloc_i, aj.parent().altloc.strip())):
          found = True
          break
      if(not found): removed.add(i)
    if(len(removed) == n_removed): break
  #
  remove_sel = flex.size_t()
  for agi, ai in water_o:
    if(ai.i_seq in removed):
      remove_sel.extend(agi.atoms().extract_i_seq())
  return ~flex.bool(model.size(), remove_sel)

def filter_by_distance(model, params=None, fixed_water_selection=None,
                       dist_min=None, dist_max=None, dist_min_h_nonpolar=None,
                       dist_min_h_polar=None):
  """
  Model without water removed by distance criteria, see
  distance_filter_selection.
  """
  selection = distance_filter_selection(
    model                 = model,
    params                = params,
    fixed_water_selection = fixed_water_selection,
    dist_min              = dist_min,
    dist_max              = dist_max,
    dist_min_h_nonpolar   = dist_min_h_nonpolar,
    dist_min_h_polar      = dist_min_h_polar)
  if(selection.count(False) > 0): model = model.select(selection)
  return model

def water_key(residue_group):
  """
  Identifier of water molecule: chain ID and residue ID.
  """
  return (residue_group.parent().id, residue_group.resid())

class water_state(object):
  """
  What needs to be remembered about water between water updates: which
  water came with input model and must be kept (keep_input_water=True).
  """
  def __init__(self, model=None, keep_input_water=False):
    self.input_water_keys = set()
    if(keep_input_water and model is not None):
      for rg in model.get_hierarchy().residue_groups():
        if(is_water_residue_group(rg)):
          self.input_water_keys.add(water_key(rg))

  def fixed_selection(self, model):
    """
    Selection of atoms of water that must be kept as is.
    """
    result = flex.bool(model.size(), False)
    if(len(self.input_water_keys) == 0): return result
    for rg in model.get_hierarchy().residue_groups():
      if(not is_water_residue_group(rg)): continue
      if(water_key(rg) in self.input_water_keys):
        result.set_selected(rg.atoms().extract_i_seq(), True)
    return result

  def update_keys(self, old_to_new):
    """
    old_to_new: dictionary old key -> new key, for water that changed chain
    ID or residue number.
    """
    result = set()
    for key in self.input_water_keys:
      result.add(old_to_new.get(key, key))
    self.input_water_keys = result

water_state_class = water_state # manager has argument water_state

def _map_flags(x, old_to_new, n, default=False):
  """
  Refinement flags x (flex.bool, flex.size_t or nested lists of flex.size_t)
  for atoms reordered as old_to_new (flex.size_t: new index of old atom) in
  a model of n atoms; default: value of flex.bool for atoms that are new.
  """
  if(x is None): return None
  if(isinstance(x, flex.bool)):
    result = flex.bool(n, default)
    result.set_selected(old_to_new, x)
    return result
  if(isinstance(x, flex.size_t)):
    return flex.size_t(sorted(old_to_new.select(x)))
  return [_map_flags(x_, old_to_new, n) for x_ in x]

def _copy_flags(x):
  if(x is None): return None
  if(isinstance(x, (flex.bool, flex.size_t))): return x.deep_copy()
  return [_copy_flags(x_) for x_ in x]

_flag_names = ["sites_individual", "sites_torsion_angles", "sites_rigid_body",
  "adp_individual_iso", "adp_individual_aniso", "adp_group", "group_h",
  "adp_tls", "s_occupancies"]

def readd_water_residue_groups(model, params, groups):
  """
  Remove water residue groups from model and add them (possibly changed)
  at the end of the model, chain by chain. Scatterers, nonbonded types and
  conformer indices of atoms that existed are carried over.

  Parameters
  ----------
  model : mmtbx.model.manager
  params : ordered_solvent parameters
  groups : list
    (chain_id, [(residue_group, sources), ...]) in output order.
    residue_group: detached residue group to add, residue number set.
    sources: for each of its atoms, i_seq of the atom of model it is a copy
    of, or -(i_seq+1) for a new H atom of water with non-H atom i_seq.

  Refinement flags and TLS groups of atoms that existed stay as they were
  (water does not lose its occupancy group, TLS group, ADP type); new H is
  refined as its O is (coordinates, isotropic ADP).

  Returns
  -------
  tuple
    New model, and dictionary: old i_seq -> new i_seq for carried atoms.
  """
  n_old = model.size()
  flags_old = None
  if(model.refinement_flags is not None):
    flags_old = dict([(name, _copy_flags(getattr(model.refinement_flags, name)))
      for name in _flag_names])
  tls_old = None
  if(model.tls_groups is not None and model.tls_groups.iselections is not None):
    tls_old = _copy_flags(model.tls_groups.iselections)
  xrs = model.get_xray_structure()
  scatterers = xrs.scatterers()
  u_isos = xrs.extract_u_iso_or_u_equiv()
  uc = xrs.unit_cell()
  rm = model.get_restraints_manager()
  nb_types, nb_charges, conf_indices = None, None, None
  if(rm is not None):
    nb_types     = rm.geometry.nonbonded_types
    nb_charges   = rm.geometry.nonbonded_charges
    conf_indices = rm.geometry.conformer_indices
  keep_sel = flex.bool(model.size(), True)
  prepared = []
  for chain_id, rgs_sources in groups:
    rgs = []
    srcs = []
    new_scatterers = flex.xray_scatterer()
    types, charges, cis = flex.std_string(), flex.int(), flex.size_t()
    for rg, sources in rgs_sources:
      rg_atoms = rg.atoms()
      assert rg_atoms.size() == len(sources)
      rgs.append(rg)
      for atom, src in zip(rg_atoms, sources):
        srcs.append(src)
        if(src >= 0):
          keep_sel[src] = False
          new_scatterers.append(scatterers[src])
          i_ref = src
          if(nb_types is not None): types.append(nb_types[src])
        else:
          i_ref = -src-1
          new_scatterers.append(xray.scatterer(
            label           = atom.name,
            site            = uc.fractionalize(atom.xyz),
            u               = u_isos[i_ref],
            occupancy       = scatterers[i_ref].occupancy,
            scattering_type = atom.element.strip().upper()))
          if(nb_types is not None): types.append('HOH2')
        if(nb_charges is not None): charges.append(nb_charges[i_ref])
        if(conf_indices is not None): cis.append(conf_indices[i_ref])
    sxrs = xray.structure(
      special_position_settings = xrs,
      scatterers                = new_scatterers)
    prepared.append((chain_id, rgs, srcs, sxrs, types, charges, cis))
  model = model.select(keep_sel)
  old_to_new = {}
  prepared = [p for p in prepared if len(p[1]) > 0]
  for i_chain, (chain_id, rgs, srcs, sxrs, types, charges, cis) in enumerate(
      prepared):
    n_before = model.size()
    add_solvent_to_model_inplace(
      model                  = model,
      params                 = params,
      water_residue_groups   = rgs,
      chain_id               = chain_id,
      conformer_indices      = (group_args(conformer_indices=cis)
                                if conf_indices is not None else None),
      nonbonded_types        = (types if nb_types is not None else None),
      nonbonded_charges      = (charges if nb_charges is not None else None),
      solvent_xray_structure = sxrs,
      keep_resseq            = True,
      update_plain_pair_sym_table = (i_chain == len(prepared)-1))
    for k, src in enumerate(srcs):
      if(src >= 0): old_to_new[src] = n_before+k
  # Atoms that existed are reordered: flags go with them
  n_new = model.size()
  permutation = flex.size_t(n_old, 0)
  k = 0
  for i in range(n_old):
    j = old_to_new.get(i)
    if(j is None):
      j = k
      k += 1
    permutation[i] = j
  if(flags_old is not None):
    rfs = model.refinement_flags
    for name in _flag_names:
      default = name in ["sites_individual", "sites_torsion_angles",
                         "adp_individual_iso"]
      setattr(rfs, name, _map_flags(flags_old[name], permutation, n_new,
        default))
  if(tls_old is not None and model.tls_groups is not None):
    model.tls_groups.iselections = _map_flags(tls_old, permutation, n_new)
  return model, old_to_new

def nearest_polymer_atoms(model, i_seqs, distance_cutoffs=(6., 12., 24.),
                          margin=0.25):
  """
  For each atom i_seq of i_seqs: the nearest non-H atom of a polymer (or, if
  model has none, of anything that is not water), symmetry copies included.
  Atoms farther than the first cutoff are assigned by the nearest CA or P
  atom (if there are any): neighbor search over all atoms at a large
  distance is expensive.

  Not to move water that is about as far from two chains (or from two
  copies of one) back and forth between updates: an atom of the chain the
  water is in is taken if it is less than margin farther than the nearest
  one, and then one that needs no symmetry operation, if any.

  Returns
  -------
  dict
    i_seq -> (distance, j_seq, rt_mx) such that atom i_seq is at distance
    from rt_mx applied to atom j_seq. Atoms with nothing within the largest
    cutoff are absent.
  """
  atoms = model.get_hierarchy().atoms()
  target = flex.bool(atoms.size(), False)
  other  = flex.bool(atoms.size(), False)
  for ag in model.get_hierarchy().atom_groups():
    if(is_water(ag)): continue
    sel = ag.atoms().extract_i_seq()
    other.set_selected(sel, True)
    if(get_class(name=ag.resname) in water_altlocs.POLYMER):
      target.set_selected(sel, True)
  if(target.count(True) == 0): target = other
  result = {}
  if(target.count(True) == 0): return result
  target = target & ~model.get_hd_selection()
  trace = flex.bool(atoms.size(), False)
  for i in target.iselection():
    if(atoms[i].name.strip() in ["CA", "P"]): trace[i] = True
  if(trace.count(True) == 0): trace = target
  todo = list(i_seqs)
  sites_cart = model.get_sites_cart()
  def chain_id(i_seq):
    return atoms[i_seq].parent().parent().parent().id
  def better(new, old):
    return old is None or new[0] < old[0]
  for i_cutoff, cutoff in enumerate(distance_cutoffs):
    if(len(todo) == 0): break
    if(i_cutoff == 0): reference = target.iselection()
    else:              reference = trace.iselection()
    # sites of interest first, then reference atoms
    indices = flex.size_t(todo)
    indices.extend(reference)
    n_todo = len(todo)
    # nearest: any; of the same chain; without symmetry operation, by chain
    best_any, best_same, best_unit = {}, {}, {}
    for i, j, rt_mx, d in crystal.pairs_within(sites_cart.select(indices),
        cutoff, model.crystal_symmetry(), i_seqs=range(n_todo)):
      if(j < n_todo): continue
      i_seq, j_seq = todo[i], indices[j]
      c = (d, j_seq, rt_mx)
      if(better(c, best_any.get(i_seq))): best_any[i_seq] = c
      cid = chain_id(j_seq)
      if(cid == chain_id(i_seq) and better(c, best_same.get(i_seq))):
        best_same[i_seq] = c
      if(rt_mx.is_unit_mx() and better(c, best_unit.get((i_seq, cid)))):
        best_unit[(i_seq, cid)] = c
    for i_seq, c in best_any.items():
      same = best_same.get(i_seq)
      if(same is not None and same[0] < c[0]+margin): c = same
      unit = best_unit.get((i_seq, chain_id(c[1])))
      if(unit is not None and unit[0] < c[0]+margin): c = unit
      result[i_seq] = c
    todo = [i for i in todo if i not in result]
  return result

def assign_water_chains(model, params=None, water_state=None, log=None):
  """
  Move each water molecule to the symmetry copy that is nearest to a polymer
  chain and give it the ID of that chain. Water forms one chain object per
  chain ID after all other chains. Water that stays in its chain keeps its
  residue number, water that changes chain gets the next number available
  in the new chain. Water bonded to anything else (a link to a metal) is not
  touched.

  Parameters
  ----------
  model : mmtbx.model.manager
  params : ordered_solvent parameters, optional
  water_state : water_state, optional
    Its keys are updated for water that changes chain ID or number.
  log : file, optional

  Returns
  -------
  mmtbx.model.manager
    The same object if nothing needed to change.
  """
  if(params is None): params = master_params().extract()
  hierarchy = model.get_hierarchy()
  if(hierarchy.models_size() != 1): return model
  atoms = hierarchy.atoms()
  rm = model.get_restraints_manager()
  bonds = None
  if(rm is not None): bonds = rm.geometry.shell_sym_tables[0]
  # Water molecules that can be handled
  waters = [] # (residue group, i_seq of reference atom)
  for rg in hierarchy.residue_groups():
    if(not is_water_residue_group(rg)): continue
    i_ref = None
    for a in rg.atoms():
      if(i_ref is None and not a.element_is_hydrogen()): i_ref = a.i_seq
    if(i_ref is None): continue
    waters.append((rg, i_ref))
  if(bonds is not None and len(waters) > 0):
    # bonds are stored once, with either atom: check the other direction
    water_atom = {}
    for rg, i_ref in waters:
      for i in rg.atoms().extract_i_seq(): water_atom[i] = rg.memory_id()
    bonded = set()
    for i in range(bonds.size()):
      for j in bonds[i].keys():
        for a, b in [(i, j), (j, i)]:
          if(a in water_atom and water_atom.get(b) != water_atom[a]):
            bonded.add(water_atom[a])
    waters = [w for w in waters if w[0].memory_id() not in bonded]
  if(len(waters) == 0): return model
  nearest = nearest_polymer_atoms(model=model, i_seqs=[w[1] for w in waters])
  if(len(nearest) == 0): return model
  # Move to symmetry copy next to the chain
  xrs = model.get_xray_structure()
  uc = xrs.unit_cell()
  scatterers = xrs.scatterers()
  n_moved = 0
  new_chain_id = {}
  for rg, i_ref in waters:
    best = nearest.get(i_ref)
    if(best is None):
      new_chain_id[rg.memory_id()] = rg.parent().id
      continue
    d, j, rt_mx = best
    new_chain_id[rg.memory_id()] = atoms[j].parent().parent().parent().id
    if(rt_mx.is_unit_mx()): continue
    op = rt_mx.inverse()
    r = matrix.sqr(op.r().as_double())
    for a in rg.atoms():
      sc = scatterers[a.i_seq]
      sc.site = op * sc.site
      if(sc.flags.use_u_aniso()):
        sc.u_star = (r * matrix.sym(sym_mat3=sc.u_star) *
          r.transpose()).as_sym_mat3()
        a.uij = adptbx.u_star_as_u_cart(uc, sc.u_star)
    n_moved += 1
  if(n_moved > 0):
    model.set_sites_cart(sites_cart = xrs.sites_cart())
    if(rm is not None): rm.geometry.reset_internals()
  # Is water already where it has to be?
  chain_order = []
  for chain in hierarchy.chains():
    is_water_chain = True
    for rg in chain.residue_groups():
      if(not is_water_residue_group(rg)): is_water_chain = False
    if(not is_water_chain and chain.id not in chain_order):
      chain_order.append(chain.id)
  for rg, i_ref in waters:
    cid = new_chain_id[rg.memory_id()]
    if(cid not in chain_order): chain_order.append(cid)
  handled = set([rg.memory_id() for rg, i_ref in waters])
  layout = []
  ok = True
  seen_water_chain = False
  for chain in hierarchy.chains():
    n_w, n_other = 0, 0
    for rg in chain.residue_groups():
      if(rg.memory_id() in handled):
        n_w += 1
        if(new_chain_id[rg.memory_id()] != chain.id): ok = False
      elif(not is_water_residue_group(rg)): n_other += 1
    if(n_w > 0 and n_other > 0): ok = False
    if(n_w > 0):
      seen_water_chain = True
      layout.append(chain.id)
    elif(seen_water_chain): ok = False
  if(ok):
    expected = [cid for cid in chain_order if cid in layout]
    if(layout != expected): ok = False
  if(ok): return model
  # Rebuild water chains
  highest = {}
  for chain in hierarchy.chains():
    for rg in chain.residue_groups():
      if(rg.memory_id() in handled and
         new_chain_id[rg.memory_id()] != chain.id): continue
      highest[chain.id] = max(highest.get(chain.id, 0), rg.resseq_as_int())
  stay, move = {}, {}
  for rg, i_ref in waters:
    cid = new_chain_id[rg.memory_id()]
    if(cid == rg.parent().id): stay.setdefault(cid, []).append(rg)
    else:                      move.setdefault(cid, []).append(rg)
  groups = []
  key_changes = {}
  for cid in chain_order:
    rgs_sources = []
    for rg in stay.get(cid, []):
      rgs_sources.append(
        (rg.detached_copy(), list(rg.atoms().extract_i_seq())))
    for rg in move.get(cid, []):
      new_rg = rg.detached_copy()
      highest[cid] = highest.get(cid, 0)+1
      new_rg.resseq = iotbx.pdb.resseq_encode(value=highest[cid])
      new_rg.icode = " "
      key_changes[water_key(rg)] = (cid, new_rg.resid())
      rgs_sources.append((new_rg, list(rg.atoms().extract_i_seq())))
    if(len(rgs_sources) > 0): groups.append((cid, rgs_sources))
  n_old = model.size()
  model, water_old_to_new = readd_water_residue_groups(
    model  = model,
    params = params,
    groups = groups)
  assert model.size() == n_old
  if(water_state is not None): water_state.update_keys(key_changes)
  if(log is not None):
    print("  water: moved to symmetry copy: %d, changed chain: %d"%(
      n_moved, len(key_changes)), file=log)
  return model

def coordination_origin_id():
  from cctbx.geometry_restraints.linking_class import linking_class
  return linking_class().get_origin_id('water coordination')

def remove_coordination_restraints(model):
  """
  Remove restraints built by update_coordination_restraints. Returns their
  number.
  """
  rm = model.get_restraints_manager()
  if(rm is None): return 0
  g = rm.geometry
  origin_id = coordination_origin_id()
  sites_cart = model.get_sites_cart()
  pair_proxies = g.pair_proxies(sites_cart = sites_cart)
  old = set()
  if(pair_proxies.bond_proxies is not None):
    for p in pair_proxies.bond_proxies.simple:
      if(p.origin_id == origin_id): old.add(tuple(sorted(p.i_seqs)))
    for p in pair_proxies.bond_proxies.asu:
      if(p.origin_id == origin_id):
        old.add(tuple(sorted((p.i_seq, p.j_seq))))
  if(len(old) > 0):
    g.remove_bond_restraints_in_place(
      bonded_pairs = sorted(old), sites_cart = sites_cart)
  return len(old)

def update_coordination_restraints(model, params=None):
  """
  Restrain the distance between each water and the atom that coordinates it
  to its current value (bond-like restraint with sigma
  params.coordination_restraints.sigma): water does not drift away from
  weak density, and the pair (only the pair) is not subject to nonbonded
  repulsion. Restraints built by previous call are removed first.

  The coordinating atom is the nearest N, O, S or metal (not H, symmetry
  copies included) within params.dist_max that can be present at the same
  time as the water; another water only if there is nothing else. Water on
  special position gets no restraint.

  Returns
  -------
  int
    Number of restraints.
  """
  if(params is None): params = master_params().extract()
  remove_coordination_restraints(model = model)
  rm = model.get_restraints_manager()
  if(rm is None): return 0
  if(not params.coordination_restraints.enabled): return 0
  g = rm.geometry
  origin_id = coordination_origin_id()
  sites_cart = model.get_sites_cart()
  # Build new ones
  atoms = model.get_hierarchy().atoms()
  interaction_selection = model.selection(
    map_to_water.selection_string_interaction)
  sst = g.site_symmetry_table
  # The restraint is made to the copy of the atom that is found by
  # add_new_bond_restraints_in_place: it has to be the only one in reach.
  reach = params.dist_max+0.1
  water_i_seqs = []
  for ag in model.get_hierarchy().atom_groups():
    if(not is_water(ag)): continue
    for a in ag.atoms():
      if(not a.element_is_hydrogen()): water_i_seqs.append(a.i_seq)
  copies = {} # (i, j) -> distances to copies of j in reach
  if(len(water_i_seqs) > 0):
    for i, j, rt_mx, d in crystal.pairs_within(sites_cart, reach,
        model.crystal_symmetry(), i_seqs=water_i_seqs):
      if(i == j or not interaction_selection[j]): continue
      copies.setdefault(i, {}).setdefault(j, []).append(d)
  pairs = {}
  for i in water_i_seqs:
    if(sst.is_special_position(i_seq=i) or i not in copies): continue
    ai = atoms[i]
    altloc_i = ai.parent().altloc.strip()
    rg_id = ai.parent().parent().memory_id()
    best, best_water = None, None
    for j in sorted(copies[i].keys()):
      ds = copies[i][j]
      d = min(ds)
      if(len(ds) > 1 or d > params.dist_max): continue
      aj = atoms[j]
      if(aj.element_is_hydrogen()): continue
      if(aj.parent().parent().memory_id() == rg_id): continue
      if(not can_coexist(altloc_i, aj.parent().altloc.strip())): continue
      if(is_water(aj.parent())):
        if(best_water is None or d < best_water[0]): best_water = (d, j)
      else:
        if(best is None or d < best[0]): best = (d, j)
    if(best is None): best = best_water
    if(best is None): continue
    d, j = best
    key = (min(i, j), max(i, j))
    if(key in pairs or g.is_bonded_atoms(i, j)): continue
    pairs[key] = d
  proxies = []
  weight = 1./params.coordination_restraints.sigma**2
  for key in sorted(pairs.keys()):
    proxies.append(geometry_restraints.bond_simple_proxy(
      i_seqs         = key,
      distance_ideal = pairs[key],
      weight         = weight,
      origin_id      = origin_id))
  if(len(proxies) > 0):
    # Only the two atoms of a restraint are excluded from nonbonded
    # interactions. Adding bonds makes atoms bonded to the coordinating atom
    # (and other water it coordinates) 1-3 and 1-4 neighbors of the water,
    # which they are not: keep 1-3 and 1-4 tables as they are.
    shells_1_3_and_1_4 = g.shell_sym_tables[1:]
    g.add_new_bond_restraints_in_place(
      proxies    = proxies,
      sites_cart = sites_cart,
      max_distance_between_connecting_atoms = params.dist_max+0.1)
    g.shell_sym_tables = [g.shell_sym_tables[0]] + shells_1_3_and_1_4
    # all-atom table of pairs depends on coordinates only: no need to
    # recompute it (expensive)
    plain_pair_sym_table = g.plain_pair_sym_table
    g.reset_internals()
    g.plain_pair_sym_table = plain_pair_sym_table
    g.pair_proxies(sites_cart = sites_cart)
  return len(proxies)

def resolve_new_solvent(params, d_min, model):
  """
  ADP type of new water: "isotropic" or "anisotropic". Auto: anisotropic if
  d_min is better than params.anisotropic_d_min and more than half of
  non-water non-H atoms of model are anisotropic.
  """
  choice = params.new_solvent
  if(choice in ["isotropic", "anisotropic"]): return choice
  if(d_min is None or params.anisotropic_d_min is None or
     d_min >= params.anisotropic_d_min): return "isotropic"
  sel = ~(model.solvent_selection() | model.get_hd_selection())
  n = sel.count(True)
  if(n == 0): return "isotropic"
  use_u_aniso = model.get_xray_structure().use_u_aniso()
  if((use_u_aniso & sel).count(True)*2 > n): return "anisotropic"
  return "isotropic"

class manager(object):
  """
  Update of ordered water: remove water that does not fit, add new water at
  difference map peaks (in alternative conformations if model has them),
  refine occupancy and B of new water, add H, assign chains, restrain water
  to its coordinating atom. Result: self.model.
  """
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     is_neutron_scat_table,
                     params = None,
                     find_peaks_params = None,
                     water_state = None,
                     log = sys.stdout):
    adopt_init_args(self, locals())
    if(self.params is None): self.params = master_params().extract()
    # Copy: do not change caller's parameters
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    else:
      self.find_peaks_params = copy.deepcopy(self.find_peaks_params)
    if(self.water_state is None): self.water_state = water_state_class()
    #
    self.ma         = libtbx.log.manager(log = self.log)
    self.total_time = 0
    self._maps      = None # computed when needed, see _get_maps()
    self._model_state = None # model state fmodel and maps correspond to
    self._peaks     = None
    self.n_water    = None
    self.model_size_init = self.model.size()
    self.new_solvent_selection = None
    self._labelled_selection = None # water that got altloc in this update
    self._saved_occ = None # of existing water that got altloc in this update
    self._saved_u_iso = None
    #
    self._call(msg="Start",             func=self._unrestrain)
    self._call(msg="Filter (dist)",     func=self._filter_dist)
    self._call(msg="Filter (q & B)",    func=self._filter_q_b)
    self._call(msg="Filter (map)",      func=self._filter_map)
    self._call(msg="Find peaks",        func=self._find_peaks)
    self._call(msg="Add new water",     func=self._add_new_solvent)
    self._call(msg="Refine new water",  func=self._refine)
    self._call(msg="Filter (q & B)",    func=self._filter_q_b_new)
    self._call(msg="Build H",           func=self._build_h)
    self._call(msg="Filter (dist only)",func=self._filter_dist)
    self._call(msg="Assign chains",     func=self._assign_chains)
    self._call(msg="Restrain water",    func=self._restrain)

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

  def _get_model_state(self):
    xrs = self.model.get_xray_structure()
    scatterers = xrs.scatterers()
    return group_args(
      size  = scatterers.size(),
      sites = xrs.sites_frac().as_double(),
      occ   = scatterers.extract_occupancies(),
      u     = scatterers.extract_u_iso_or_u_equiv(xrs.unit_cell()),
      n_aniso = xrs.use_u_aniso().count(True))

  def _model_changed(self):
    s1, s2 = self._model_state, self._get_model_state()
    if(s1 is None or s1.size != s2.size): return True
    if(s1.n_aniso != s2.n_aniso): return True
    for a, b in [(s1.sites, s2.sites), (s1.occ, s2.occ), (s1.u, s2.u)]:
      if(a.size()>0 and flex.max(flex.abs(a - b)) > 0): return True
    return False

  def _get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps(self):
    """
    Update fmodel (and invalidate maps) only if model changed: this is
    expensive (mask, structure factors).
    """
    self.n_water = self.model.solvent_selection().count(True)
    if(not self._model_changed()): return
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = True)
    self._maps = None
    self._model_state = self._get_model_state()

  def _assert_same_model(self):
    mmtbx.utils.assert_xray_structures_equal( # XXX MAKE METHOD OF XRS
      x1 = self.model.get_xray_structure(),
      x2 = self.fmodel.xray_structure,
      eps=1.e-3)

  def _select(self, selection):
    """
    Select atoms of the model and of everything that is per atom.
    """
    if(selection.count(False) == 0): return
    self.model = self.model.select(selection)
    if(self.new_solvent_selection is not None):
      self.new_solvent_selection = self.new_solvent_selection.select(selection)
    if(self._labelled_selection is not None):
      self._labelled_selection = self._labelled_selection.select(selection)
    if(self._saved_occ is not None):
      self._saved_occ   = self._saved_occ.select(selection)
      self._saved_u_iso = self._saved_u_iso.select(selection)

  def _fixed_selection(self):
    return self.water_state.fixed_selection(model = self.model)

  def _d_min(self):
    return self.fmodel.f_obs().d_min()

  def _anisotropic(self):
    return resolve_new_solvent(params=self.params, d_min=self._d_min(),
      model=self.model) == "anisotropic"

  def _build_h(self):
    """
    Add H to water (also to water in alternative conformations) at
    difference map peaks. Water on special position gets no H. H that water
    has already and that are not in density are moved to the peaks.

    The peaks are weak and distorted, so they are accepted with loose
    tolerances (add_hydrogens parameters). They define the orientation of
    water only: H are placed with ideal O-H distances and H-O-H angle (as
    restraints define them, see ideal_water_h_sites), so refinement that
    follows starts from sensible geometry. Nothing is done in filter_only
    mode (final clean up: no refinement follows).
    """
    if(self.params.mode == "filter_only"): return
    ahp = self.params.add_hydrogens
    if ahp.d_min is None: return
    if self._d_min() > ahp.d_min: return
    if self.fmodel.r_work() > ahp.r_work: return
    self._move_misplaced_h()
    mpl = maptbx.MapPeakLocator(
      map_data    = self._get_maps().difference_map,
      unit_cell   = self.model.crystal_symmetry().unit_cell(),
      is_periodic = True)
    self.model.get_hierarchy().atoms().reset_i_seq()
    self.model.get_hierarchy().atoms_reset_serial()
    fixed = self._fixed_selection()
    sst = self.model.get_xray_structure().site_symmetry_table()
    hoh_residue_groups = {} # chain id -> [(water residue group with H, sources)]
    chain_ids = []
    for m in self.model.get_hierarchy().models():
      for c in m.chains():
        for rg in c.residue_groups():
          if(not is_water_residue_group(rg)): continue
          used = [] # H peaks taken by other conformers of this water
          h_sites = {} # index of atom group -> H sites
          for i_ag, ag in enumerate(rg.atom_groups()):
            o_atom = None
            skip = False
            for atom in ag.atoms():
              if(atom.element_is_hydrogen()): skip = True # has H already
              elif(atom.element.strip().upper()=="O"): o_atom = atom
            if(skip or o_atom is None): continue
            if(fixed[o_atom.i_seq]): continue
            if(sst.is_special_position(i_seq=o_atom.i_seq)): continue
            sites_cart, peaks_heights = mpl.get_peaks_within_radius(
              target_cart = o_atom.xyz,
              R           = ahp.max_dist,
              threshold   = ahp.diff_map_peak_threshold)
            sites_, heights_ = [], []
            for s, h in zip(sites_cart, peaks_heights):
              taken = False
              for u in used:
                if(abs(matrix.col(s)-matrix.col(u)) < 1.e-3): taken = True
              if(not taken):
                sites_.append(s)
                heights_.append(h)
            wo = find_water_hydrogens(
              o_site       = o_atom.xyz,
              peak_heights = heights_,
              peak_coords  = sites_,
              min_od       = ahp.min_dist,
              max_od       = ahp.max_dist,
              min_angle    = ahp.min_angle,
              max_angle    = ahp.max_angle)
            if wo is not None:
              h_sites[i_ag] = wo.h_sites
              used.extend(wo.h_sites)
          if(len(h_sites) > 0):
            hoh_rg, sources = get_detached_water_residue_group(
              residue_group = rg,
              h_sites_cart  = h_sites)
            if(c.id not in chain_ids): chain_ids.append(c.id)
            hoh_residue_groups.setdefault(c.id, []).append((hoh_rg, sources))
    if len(hoh_residue_groups) > 0:
      self.model, old_to_new = readd_water_residue_groups(
        model  = self.model,
        params = self.params,
        groups = [(cid, hoh_residue_groups[cid]) for cid in chain_ids])
      # Water that got H is at the end of the model now
      n_changed = 0
      for cid in chain_ids:
        for hoh_rg, sources in hoh_residue_groups[cid]:
          n_changed += len(sources)
      self._idealize_new_h(
        i_seq_first    = self.model.size() - n_changed,
        i_seqs_existed = set(old_to_new.values()))
      self._model_state = None # fmodel needs the new xray structure
      self.new_solvent_selection = None
      self._labelled_selection = None
      self._saved_occ, self._saved_u_iso = None, None

  def _idealize_new_h(self, i_seq_first, i_seqs_existed):
    """
    Give ideal geometry to water that just got its two H, keeping the
    orientation. The new H are atoms from i_seq_first on that are not in
    i_seqs_existed. Restraints for the new H exist only now, in the model
    with H added.
    """
    self.model.get_hierarchy().atoms().reset_i_seq()
    water_i_seqs = []
    for rg in self.model.get_hierarchy().residue_groups():
      if(not is_water_residue_group(rg)): continue
      for ag in rg.atom_groups():
        o_atoms, h_atoms = [], []
        for atom in ag.atoms():
          if(atom.element_is_hydrogen()):
            if(atom.i_seq >= i_seq_first and
               atom.i_seq not in i_seqs_existed): h_atoms.append(atom)
          elif(atom.element.strip().upper()=="O"): o_atoms.append(atom)
        if(len(o_atoms) != 1 or len(h_atoms) != 2): continue
        water_i_seqs.append(
          (o_atoms[0].i_seq, h_atoms[0].i_seq, h_atoms[1].i_seq))
    sites_cart = self.model.get_sites_cart()
    ideal = ideal_water_h_geometry(
      model = self.model, water_i_seqs = water_i_seqs)
    for (o, h1, h2), geometry in zip(water_i_seqs, ideal):
      if(geometry is None): continue
      sites_cart[h1], sites_cart[h2] = ideal_water_h_sites(
        o_site    = sites_cart[o],
        h_sites   = (sites_cart[h1], sites_cart[h2]),
        distances = geometry[0],
        angle     = geometry[1])
    self.model.set_sites_cart(sites_cart = sites_cart)

  def _difference_map_without(self, i_seqs):
    """
    Difference map computed with atoms i_seqs left out of the model.
    """
    xrs = self.model.get_xray_structure()
    xrs_omit = xrs.deep_copy_scatterers()
    xrs_omit.set_occupancies(value=0, selection=flex.size_t(i_seqs))
    self.fmodel.update_xray_structure(
      xray_structure = xrs_omit,
      update_f_calc  = True,
      update_f_mask  = True)
    p = self.params
    result = maps(
      fmodel               = self.fmodel,
      model                = self.model,
      mask_atoms_selection = p.mask_atoms_selection,
      difference_map_type  = p.primary_map_type,
      model_map_type       = p.secondary_map_and_map_cc_filter.cc_map_1_type,
      data_map_type        = p.secondary_map_and_map_cc_filter.cc_map_2_type
      ).difference_map
    self.fmodel.update_xray_structure(
      xray_structure = xrs,
      update_f_calc  = True,
      update_f_mask  = True)
    return result

  def _move_misplaced_h(self):
    """
    Move H of water that are not at difference map peaks to the peaks.

    Water that comes with two H is not given new H. Its H are refined as
    individual atoms: if they point the wrong way they stay there, and their
    occupancy refines to zero. Here H peaks are searched for such water in the
    same way as for water without H, in a difference map computed without H
    of water. The two peaks found define where H should be: at the peaks, or,
    if restraints define ideal geometry of the water, at the sites that have
    the ideal geometry and the orientation of the peaks (see _build_h). If an
    H is further than add_hydrogens.min_dist from where it should be, both H
    are moved there and get occupancy of their O. H are left as they are if
    no two peaks are found, and in filter_only mode (final clean up: no
    refinement follows).
    """
    if(self.params.mode == "filter_only"): return
    ahp = self.params.add_hydrogens
    self.model.get_hierarchy().atoms().reset_i_seq()
    fixed = self._fixed_selection()
    sst = self.model.get_xray_structure().site_symmetry_table()
    waters = [] # (O atom, [H atom, H atom])
    for rg in self.model.get_hierarchy().residue_groups():
      if(not is_water_residue_group(rg)): continue
      for ag in rg.atom_groups():
        o_atoms, h_atoms = [], []
        for atom in ag.atoms():
          if(atom.element_is_hydrogen()): h_atoms.append(atom)
          elif(atom.element.strip().upper()=="O"): o_atoms.append(atom)
        if(len(o_atoms) != 1 or len(h_atoms) != 2): continue
        if(fixed[o_atoms[0].i_seq]): continue
        if(sst.is_special_position(i_seq=o_atoms[0].i_seq)): continue
        # H (not D) is a negative peak in nuclear map
        if(self.is_neutron_scat_table and
           [a.element.strip().upper() for a in h_atoms].count("H") > 0):
          continue
        waters.append((o_atoms[0], h_atoms))
    if(len(waters) == 0): return
    h_i_seqs = []
    for o_atom, h_atoms in waters:
      for h in h_atoms: h_i_seqs.append(h.i_seq)
    mpl = maptbx.MapPeakLocator(
      map_data    = self._difference_map_without(i_seqs = h_i_seqs),
      unit_cell   = self.model.crystal_symmetry().unit_cell(),
      is_periodic = True)
    sites_cart = self.model.get_sites_cart()
    occupancies = self.model.get_occ()
    ideal = ideal_water_h_geometry(
      model        = self.model,
      water_i_seqs = [(o_atom.i_seq, h_atoms[0].i_seq, h_atoms[1].i_seq)
        for o_atom, h_atoms in waters])
    used = {} # residue group -> H peaks taken by its other conformers
    n_moved = 0
    for (o_atom, h_atoms), geometry in zip(waters, ideal):
      peak_sites, peak_heights = mpl.get_peaks_within_radius(
        target_cart = o_atom.xyz,
        R           = ahp.max_dist,
        threshold   = ahp.diff_map_peak_threshold)
      rg_id = o_atom.parent().parent().memory_id()
      if(rg_id not in used): used[rg_id] = []
      used_rg = used[rg_id]
      sites_, heights_ = [], []
      for s, h in zip(peak_sites, peak_heights):
        taken = False
        for u in used_rg:
          if(abs(matrix.col(s)-matrix.col(u)) < 1.e-3): taken = True
        if(not taken):
          sites_.append(s)
          heights_.append(h)
      wo = find_water_hydrogens(
        o_site       = o_atom.xyz,
        peak_heights = heights_,
        peak_coords  = sites_,
        min_od       = ahp.min_dist,
        max_od       = ahp.max_dist,
        min_angle    = ahp.min_angle,
        max_angle    = ahp.max_angle)
      if(wo is None): continue
      used_rg.extend(wo.h_sites)
      # each H has the peak that makes the two H move least
      h_sites = wo.h_sites
      if(geometry is not None):
        h_sites = ideal_water_h_sites(
          o_site    = o_atom.xyz,
          h_sites   = h_sites,
          distances = geometry[0],
          angle     = geometry[1])
      a, b = [matrix.col(s) for s in h_sites]
      h1, h2 = [matrix.col(h.xyz) for h in h_atoms]
      if(abs(h1-a)+abs(h2-b) > abs(h1-b)+abs(h2-a)): a, b = b, a
      if(max(abs(h1-a), abs(h2-b)) < ahp.min_dist): continue # H are at peaks
      for h, site in zip(h_atoms, [a, b]):
        sites_cart[h.i_seq] = site.elems
        occupancies[h.i_seq] = o_atom.occ
      n_moved += 1
    if(n_moved == 0): return
    self.model.set_sites_cart(sites_cart = sites_cart)
    self.model.set_occupancies(values = occupancies)
    self.ma.add_and_show("  water with H moved to peaks: %d"%n_moved)

  def _get_maps(self):
    """
    Maps for current model; computed only if model changed since last time.
    """
    if(self._model_changed()):
      self._get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps()
    if(self._maps is None):
      p = self.params
      self._maps = maps(
        fmodel               = self.fmodel,
        model                = self.model,
        mask_atoms_selection = p.mask_atoms_selection,
        difference_map_type  = p.primary_map_type,
        model_map_type       = p.secondary_map_and_map_cc_filter.cc_map_1_type,
        data_map_type        = p.secondary_map_and_map_cc_filter.cc_map_2_type)
    return self._maps

  def _filter_dist(self):
    self._select(distance_filter_selection(
      model                 = self.model,
      params                = self.params,
      fixed_water_selection = self._fixed_selection()))

  def _filter_q_b(self):
    self._filter(filter_occ=True, filter_adp=True)

  def _filter_q_b_new(self):
    self._filter_q_b()
    self._revert_orphan_altlocs()
    self._set_adp_type()

  def _filter_map(self):
    self._filter(filter_map=True)

  def _filter(self,
              filter_map=False,
              filter_occ=False,
              filter_adp=False):
    r_work = self.fmodel.r_work()
    mfp = self.params.secondary_map_and_map_cc_filter
    maps_ = None
    if filter_map: maps_ = self._get_maps()
    selection = flex.bool(self.model.size(), True)
    fixed = self._fixed_selection()
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
            if(fixed.select(i_seqs).count(True) > 0): continue
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
                # B at a limit stays: a value refined with that limit as
                # the bound comes back as the limit give or take rounding
                # errors (-4.e-15 for b_iso_min=0)
                if(b_isos[i_seq] < self.params.b_iso_min - 1.e-6 or
                   b_isos[i_seq] > self.params.b_iso_max + 1.e-6): keep = False
              #
              if filter_map and r_work*100 > 5.:
                good_map = maps_.score_atom(
                  atom      = atom,
                  min_cc    = mfp.poor_cc_threshold,
                  min_value = mfp.poor_map_value_threshold)
                if(not good_map.result):
                  keep = False
            if(not has_oxygen): keep=False
            if(not keep):
              selection = selection.set_selected(i_seqs, False)
    self._select(selection)

  def _revert_orphan_altlocs(self):
    """
    Water that got altloc in this update because of another water that is
    gone by now: back to no altloc.
    """
    if(self._labelled_selection is None or
       self._labelled_selection.count(True) == 0): return
    atoms = self.model.get_hierarchy().atoms()
    excl = exclusion_distances(
      dist_min            = self.params.dist_min,
      dist_min_h_nonpolar = self.params.dist_min_h_nonpolar,
      dist_min_h_polar    = self.params.dist_min_h_polar)
    neighbors = get_water_neighbors(model=self.model,
      distance_cutoff=excl.max())
    scatterers = self.model.get_xray_structure().scatterers()
    index_by_altloc = conformer_index_by_altloc(self.model)
    changed = False
    for i in self._labelled_selection.iselection():
      ai = atoms[i]
      if(ai.element_is_hydrogen()): continue
      altloc_i = ai.parent().altloc.strip()
      if(altloc_i == ""): continue
      justified = False
      for j, d in neighbors.get(i, []):
        if(j == i): continue
        aj = atoms[j]
        altloc_j = aj.parent().altloc.strip()
        if(altloc_j != "" and altloc_j != altloc_i and d < excl(aj)):
          justified = True
          break
      if(not justified):
        set_altloc(model=self.model, atom_group=ai.parent(), altloc="",
          index_by_altloc=index_by_altloc)
        changed = True
        # existing water: as it was before this update
        if(self._saved_occ is not None and self._saved_occ[i] >= 0):
          for a in ai.parent().atoms():
            k = a.i_seq
            if(self._saved_occ[k] < 0): continue
            sc = scatterers[k]
            sc.occupancy = self._saved_occ[k]
            a.occ = sc.occupancy
            if(self._saved_u_iso[k] >= 0):
              sc.u_iso = self._saved_u_iso[k]
              a.b = adptbx.u_as_b(sc.u_iso)
    if(changed): self.model._update_atom_selection_cache()

  def _set_adp_type(self):
    """
    New water is refined with isotropic ADP in this update; make it
    anisotropic if this is what is asked for.
    """
    if(self.new_solvent_selection is None): return
    if(not self._anisotropic()): return
    xrs = self.model.get_xray_structure()
    sel = self.new_solvent_selection & ~self.model.get_hd_selection() & \
      ~xrs.use_u_aniso()
    if(self._saved_occ is not None):
      sel = sel & (self._saved_occ < 0) # not water that was there before
    if(sel.count(True) == 0): return
    xrs.convert_to_anisotropic(selection = sel)
    uc = xrs.unit_cell()
    scatterers = xrs.scatterers()
    atoms = self.model.get_hierarchy().atoms()
    for i in sel.iselection():
      atoms[i].uij = adptbx.u_star_as_u_cart(uc, scatterers[i].u_star)
    rfs = self.model.refinement_flags
    if(rfs is not None and rfs.individual_adp):
      if(rfs.adp_individual_iso is not None):
        rfs.adp_individual_iso.set_selected(sel, False)
      if(rfs.adp_individual_aniso is not None):
        rfs.adp_individual_aniso.set_selected(sel, True)

  def _refine(self):
    if(self.params.mode == "filter_only"): return
    if(self.new_solvent_selection is None or
       self.new_solvent_selection.count(True)==0 or
       self.n_water==0): return
    for i in range(self.params.n_cycles):
      self.refine_oat()

  def _find_peaks(self):
    if(self.params.mode == "filter_only"): return
    fpp = self.find_peaks_params
    fpp.max_number_of_peaks=self.model.get_number_of_atoms()
    if(len(model_altlocs(self.model)) > 0):
      # sites can be next to atoms in alternative conformations and to each
      # other: assign_conformers decides
      fpp.peak_search.min_cross_distance = min_water_separation(
        params=self.params, d_min=self._d_min(),
        anisotropic=self._anisotropic())
      fpp.map_next_to_model.min_model_peak_dist = self.params.dist_min_altloc
    assert self.params.primary_map_type is not None
    self._peaks = find_peaks.manager(
      xray_structure = self.fmodel.xray_structure,
      map_data       = self._get_maps().difference_map, # diff-map
      map_cutoff     = self.params.primary_map_cutoff,
      params         = fpp,
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

  def _add_new_solvent(self):
    """
    Add water at the peaks, in alternative conformations where needed.
    Existing atoms stay as they are, except water that gets altloc because a
    new water is its alternative.
    """
    if(self._peaks is None): return
    if(self._peaks.sites.size()==0): return
    if(self.params.mode == "filter_only"): return
    uc = self.model.crystal_symmetry().unit_cell()
    ac = assign_conformers(
      model                 = self.model,
      sites_cart            = uc.orthogonalize(self._peaks.sites),
      heights               = self._peaks.heights,
      params                = self.params,
      d_min                 = self._d_min(),
      fixed_water_selection = self._fixed_selection(),
      anisotropic           = self._anisotropic())
    if(ac.keep.count(True) == 0): return
    atoms = self.model.get_hierarchy().atoms()
    relabelled = flex.size_t()
    index_by_altloc = conformer_index_by_altloc(self.model)
    for i_seq in sorted(ac.relabel.keys()):
      ag = atoms[i_seq].parent()
      set_altloc(model=self.model, atom_group=ag, altloc=ac.relabel[i_seq],
        index_by_altloc=index_by_altloc)
      relabelled.extend(ag.atoms().extract_i_seq())
    n_old = self.model.size()
    xrs = self.model.get_xray_structure()
    saved_occ   = xrs.scatterers().extract_occupancies()
    saved_u_iso = xrs.extract_u_iso_or_u_equiv()
    saved_occ.set_selected(~flex.bool(n_old, relabelled), -1)
    saved_u_iso.set_selected(xrs.use_u_aniso(), -1) # not refined
    altlocs = [a for a, k in zip(ac.altlocs, ac.keep) if k]
    self.new_solvent_selection = add_solvent_to_model_inplace(
      sites       = self._peaks.sites.select(ac.keep),
      model       = self.model,
      params      = self.params,
      altlocs     = altlocs,
      occupancies = ac.occupancies.select(ac.keep))
    self._labelled_selection = flex.bool(self.model.size(), False)
    for k, altloc in enumerate(altlocs):
      if(altloc != ""): self._labelled_selection[n_old+k] = True
    # relabelled water: its occupancy and B need to be refined too
    self.new_solvent_selection.set_selected(relabelled, True)
    self._labelled_selection.set_selected(relabelled, True)
    n_new = self.model.size()-n_old
    saved_occ.extend(flex.double(n_new, -1))
    saved_u_iso.extend(flex.double(n_new, -1))
    self._saved_occ, self._saved_u_iso = saved_occ, saved_u_iso

  def refine_oat(self):
    if self.new_solvent_selection is None: return
    selection = self.new_solvent_selection & ~self.model.get_hd_selection()
    if(self.params.refine_oat and selection.count(True)>0):
      from phenix.programs import oat
      atoms = self.model.get_hierarchy().atoms()
      scatterers = self.fmodel.xray_structure.scatterers()
      for i, sel in enumerate(selection):
        if not sel: continue
        scatterers[i].occupancy=0
        scatterers[i].u_iso=0
        self.fmodel.update_xray_structure(update_f_calc=True)
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
      self.model.adopt_xray_structure(
        xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), OAT, final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)
      #
    if([self.params.refine_adp , self.params.refine_occupancies].count(True)>0 and
       selection.count(True)>0):
      from mmtbx.refinement import wrappers
      wrappers.unrestrained_qbr_fsr(
        fmodel     = self.fmodel,
        model      = self.model,
        refine_xyz = False,
        refine_q   = self.params.refine_occupancies,
        refine_b   = self.params.refine_adp,
        selection  = selection,
        selection_b = selection & ~self.fmodel.xray_structure.use_u_aniso(),
        q_min      = 0.004,
        q_max      = self.params.occupancy_max,
        # beyond the limit of the filter that follows: water that does not
        # belong there gets a B that removes it
        b_max      = self.params.b_iso_max+10,
        b_min      = self.params.b_iso_min,
        freeze_unselected = True,
        log        = self.log)
      self.model.adopt_xray_structure(
          xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), MIN, final r_work=%6.4f r_free=%6.4f"%(
          self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def _assign_chains(self):
    model = assign_water_chains(
      model       = self.model,
      params      = self.params,
      water_state = self.water_state,
      log         = self.log)
    if(model is not self.model):
      # atoms are reordered: fmodel needs the new xray structure even if
      # coordinates happen to be in the same order
      self.model = model
      self._model_state = None
    self.new_solvent_selection = None
    self._labelled_selection = None
    self._saved_occ, self._saved_u_iso = None, None

  def _unrestrain(self):
    remove_coordination_restraints(model = self.model)
    # atom on special position is exactly there in xray structure
    xrs = self.model.get_xray_structure()
    if(xrs.special_position_indices().size() > 0):
      self.model.set_sites_cart(sites_cart = xrs.sites_cart())

  def _restrain(self):
    n = update_coordination_restraints(
      model = self.model, params = self.params)
    if(n > 0):
      print("  water coordination restraints: %d"%n, file=self.log)
