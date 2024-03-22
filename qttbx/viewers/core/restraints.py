from cctbx.array_family import flex
import sys
import pandas as pd
import numpy as np
from cctbx import geometry_restraints
import mmtbx
from iotbx.data_manager import DataManager
from collections import defaultdict
import gemmi
from mmtbx.monomer_library.pdb_interpretation import (
geometry_restraints_proxy_registries,
ener_lib_as_nonbonded_params,
master_params,
)
from .experimental_interpret import (
  DummySourceInfo,
)
from cctbx import crystal
from cctbx.crystal import special_position_settings
from libtbx import group_args
from .views import ModelView

def make_proxies(model=None,
                bonds=None,
                angles=None,
                torsions=None,
                planes=None,
                chirals = None,
                params=None,
                ):
  """
  Make restraint proxies.
  This function assumes all the necessary information is already present
  and provided as Pandas Dataframes. This function's role is purely to translate
  that information into cctbx data structures
  """
  # Make sure have params
  if params is None:
    params = master_params.extract()

  # extract basic objects
  hierarchy = model.get_hierarchy()
  crystal_symmetry = model.crystal_symmetry()
  n_seq = len(hierarchy.atoms())




  # Init restraints proxy registry
  geometry_proxy_registries = geometry_restraints_proxy_registries(
    n_seq=n_seq,
    strict_conflict_handling=False)
  geometry_proxy_registries.initialize_tables()
  sites_cart = model.get_sites_cart()


  # Bonds
  if bonds is not None:
    # name resolution
    if "ideal" in bonds and not "value_dist" in bonds:
      bonds["value_dist"] = bonds["ideal"]

    bond_source_info = DummySourceInfo(n_expected_atoms=2)
    for t in bonds.itertuples():

      #print(t.i_seq_1,t.i_seq_2)
      proxy = geometry_restraints.bond_simple_proxy(
              i_seqs=[int(t.i_seq_1),int(t.i_seq_2)],
              distance_ideal=float(t.value_dist),
              weight=float(t.weight) if t.weight is not pd.NA else 0,
              origin_id=0,
              )
          #r = geometry_restraints.bond(sites_cart=sites_cart, proxy=proxy)

      registry_process_result = geometry_proxy_registries.bond_simple.process(
                    source_info=bond_source_info,
                    proxy=proxy)


  # Angles
  if angles is not None:
    # Convienence name resolution
    if "ideal" in angles and not "value_angle" in angles:
      angles["value_angle"] = angles["ideal"]

    angle_source_info = DummySourceInfo(n_expected_atoms=3)
    for t in angles.itertuples():
      proxy = geometry_restraints.angle_proxy(
      i_seqs=[int(t.i_seq_1),int(t.i_seq_2),int(t.i_seq_3)],
      angle_ideal=float(t.value_angle),
      weight=float(t.weight) if t.weight is not pd.NA else 0,
      origin_id=0)


      registry_process_result = geometry_proxy_registries.angle.process(
                        source_info=angle_source_info,
                        proxy=proxy)

      # evaluate_registry_process_result(
      #         proxy_label="angle", m_i=m_i, m_j=m_j, i_seqs=i_seqs,
      #         registry_process_result=registry_process_result)

  if params is None:
    params = master_params.extract()


  # torsions
  if torsions is not None:
    # Convienence name resolution
    if "ideal" in torsions and not "value_angle" in torsions:
      torsions["value_angle"] = torsions["ideal"]
    if "harmonic" in torsions and not "periodicity" in torsions:
      torsions["periodicity"] = torsions["harmonic"]
      torsions["periodicity"].replace({None:0},inplace=True)
    if "alt_angle_ideals" not in torsions:
      torsions["alt_angle_ideals"] = pd.NA

    torsion_source_info = DummySourceInfo(n_expected_atoms=4)
    for i,t in enumerate(torsions.itertuples()):
      i_seqs=[int(t.i_seq_1),int(t.i_seq_2),int(t.i_seq_3),int(t.i_seq_4)]

      proxy = geometry_restraints.dihedral_proxy(
      i_seqs=i_seqs,
      angle_ideal=float(t.value_angle),
      weight=float(t.weight),
      periodicity=int(t.periodicity),
      alt_angle_ideals=t.alt_angle_ideals if t.alt_angle_ideals  is not pd.NA else None,
      origin_id=0)
      #proxy = proxy.sort_i_seqs()
      #assert proxy.i_seqs == tuple(i_seqs)


      tor_registry_process_result = geometry_proxy_registries.dihedral.process(
                        source_info=torsion_source_info,
                        proxy=proxy)

      geometry_proxy_registries.dihedral.proxies[i] = proxy
    #assert prox_out.i_seqs == proxy.i_seqs, f"i_seqs out: {proxy.i_seqs}, i_seqs out: {prox_out.i_seqs}"
  # chirals
  if chirals is not None:
    # Convienence name resolution
    if "ideal" in chirals and not "volume_ideal" in chirals:
      chirals["volume_ideal"] = chirals["ideal"]

    chiral_source_info = DummySourceInfo(n_expected_atoms=4)
    for i,t in enumerate(chirals.itertuples()):
      i_seqs=[int(t.i_seq_1),int(t.i_seq_2),int(t.i_seq_3),int(t.i_seq_4)]

      proxy = geometry_restraints.chirality_proxy(
      i_seqs=i_seqs,
      volume_ideal=float(t.volume_ideal),
      weight=int(t.weight),
      both_signs = bool(t.both_signs))


      chir_registry_process_result = geometry_proxy_registries.chirality.process(
                        source_info=chiral_source_info,
                        proxy=proxy)

      geometry_proxy_registries.chirality.proxies[i] = proxy
  #planes
  if planes is not None:

    plane_source_info = DummySourceInfo(n_expected_atoms=4) # Need other than 4?
    for t in planes.itertuples():
      i_seqs = flex.size_t([int(getattr(t,f"i_seq_{i}")) for i in range(len(planes.columns)) if hasattr(t,f"i_seq_{i}") and (not pd.isna(getattr(t,f"i_seq_{i}")))])
      weights = flex.double([int(getattr(t,f"weight_{i}")) for i in range(len(planes.columns)) if hasattr(t,f"weight_{i}") and (not pd.isna(getattr(t,f"weight_{i}")))])

      proxy = geometry_restraints.planarity_proxy(
        i_seqs=i_seqs,
        weights=weights)


      registry_process_result = geometry_proxy_registries.planarity.process(
                        source_info=plane_source_info,
                        proxy=proxy)
  # Finish
  return geometry_proxy_registries


def prep_for_grm(sites,model,geometry_proxy_registries,params=None,nonbonded_params=None):
  """
  Prepare necessary inputs to initialize the geometry restraints manager.
  """
  # Make sure have params
  if params is None:
    params = master_params.extract()


  if nonbonded_params is None:
    ener_lib = mmtbx.monomer_library.server.ener_lib() # slow
    nonbonded_params = ener_lib_as_nonbonded_params(
          ener_lib=ener_lib,
          assume_hydrogens_all_missing=True,
          factor_1_4_interactions=params.vdw_1_4_factor,
          default_distance=params.default_vdw_distance,
          minimum_distance=params.min_vdw_distance,
          const_shrink_donor_acceptor=params.const_shrink_donor_acceptor)


  # Basic structures
  hierarchy = model.get_hierarchy()
  pdb_atoms = hierarchy.atoms()
  sites_cart = model.get_sites_cart()
  crystal_symmetry = model.crystal_symmetry()
  n_seq = len(hierarchy.atoms())
  modelv = ModelView(model)


  # sym excl indices
  sym_excl_residue_groups = []
  sym_excl_indices = flex.size_t(n_seq, 0)
  for rg in hierarchy.residue_groups():
    rg_atoms = rg.atoms()
    if (rg_atoms.extract_occ().all_lt(1.0)):
      sym_excl_residue_groups.append(rg)
      sym_excl_indices.set_selected(
        rg_atoms.extract_i_seq(), len(sym_excl_residue_groups))

  # donor acceptor excl groups
  donor_acceptor_excl_groups = flex.size_t(n_seq, 0)
  counter = 0
  for rg in hierarchy.residue_groups():
    rg_atoms = rg.atoms()
    donor_acceptor_excl_groups.set_selected(
      rg_atoms.extract_i_seq(), counter)
    counter += 1

  geometry_proxy_registries.initialize_tables()
  # site_symmetry_table
  sps = special_position_settings(crystal_symmetry)
  site_symmetry_table = \
        sps.site_symmetry_table(
          sites_cart=sites_cart,
          unconditional_general_position_flags=(
            pdb_atoms.extract_occ() != 1))

  # Init tables
  # bond params table
  bond_params_table = geometry_restraints.extract_bond_params(
        n_seq=sites_cart.size(),
        bond_simple_proxies=geometry_proxy_registries.bond_simple.proxies)

  geometry_proxy_registries.initialize_tables()

  max_bond_distance = params.max_reasonable_bond_distance


  asu_mappings = sps.asu_mappings(
        buffer_thickness=max_bond_distance*3,
        )
  asu_mappings.process_sites_cart(
        original_sites=sites_cart,
        site_symmetry_table=site_symmetry_table)

  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  geometry_restraints.add_pairs(
        bond_asu_table, geometry_proxy_registries.bond_simple.proxies)

  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
        pair_asu_table=bond_asu_table,
        max_shell=3)

  shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
    for shell_asu_table in shell_asu_tables]

  #ener_lib = mmtbx.monomer_library.server.ener_lib() # slow



  # ortho matrix
  #orthogonalization_matrix=crystal_symmetry.unit_cell().orthogonalization_matrix()

  # Finish
  return {
    "bond_params_table":bond_params_table,
    "shell_asu_tables":shell_asu_tables,
    "shell_sym_tables":shell_sym_tables,
    "site_symmetry_table":site_symmetry_table,
    "model_indices":modelv.model_indices,
    "conformer_indices":modelv.conformer_indices,
    "sym_excl_indices":sym_excl_indices,
    "donor_acceptor_excl_groups":donor_acceptor_excl_groups,
    "nonbonded_params":nonbonded_params,
    "nonbonded_types":flex.std_string(sites["nonbonded_types"].fillna("")),
    "nonbonded_charges":flex.int(sites["nonbonded_charges"].fillna(0).astype(int)),
    "nonbonded_function":geometry_restraints.prolsq_repulsion_function(c_rep=100),
    "nonbonded_distance_cutoff":params.nonbonded_distance_cutoff,
    "nonbonded_buffer":params.nonbonded_buffer,
    "orthogonalization_matrix":None,#orthogonalization_matrix,
    "params":params,
  }



def build_grm(model, # mmtbx.model.manager
              geometry_proxy_registries, # proxies composed as registry
              nonbonded_types,
              nonbonded_charges,
              prep, # dict of params from prep function
              ):
  """
  Build the geometry restraints manager.
  The maximum amount of pre-processing that can occur should occur before this
  function
  """
  params = prep["params"]
  grm = geometry_restraints.manager.manager(
      crystal_symmetry=model.crystal_symmetry(),
      model_indices=prep["model_indices"],
      conformer_indices=prep["conformer_indices"],
      sym_excl_indices=prep["sym_excl_indices"],
      donor_acceptor_excl_groups=prep["donor_acceptor_excl_groups"],
      site_symmetry_table=prep["site_symmetry_table"],
      bond_params_table=prep["bond_params_table"],
      shell_sym_tables=prep["shell_sym_tables"],
      nonbonded_params=prep["nonbonded_params"],
      nonbonded_types=prep["nonbonded_types"],
      nonbonded_charges=prep["nonbonded_charges"],
      nonbonded_function=prep["nonbonded_function"],
      nonbonded_distance_cutoff=prep["nonbonded_distance_cutoff"],
      nonbonded_buffer=prep["nonbonded_buffer"],

      angle_proxies=geometry_proxy_registries.angle.proxies,
      dihedral_proxies=geometry_proxy_registries.dihedral.proxies,
      chirality_proxies=None,#self.geometry_proxy_registries.chirality.proxies,
      planarity_proxies=geometry_proxy_registries.planarity.proxies,
      parallelity_proxies=geometry_proxy_registries.parallelity.proxies,
      ramachandran_manager=None,
      external_energy_function=None,
      max_reasonable_bond_distance=params.max_reasonable_bond_distance,
      plain_pairs_radius=5,
      log=sys.stdout)
  return grm


class RestraintsBuilder:
  def __init__(self):
    pass


  def build_grm_from_dataframe_mol(self,dataframe_mol):
    geometry_proxy_registries = make_proxies(model=dataframe_mol.model,
                                           bonds=dataframe_mol.bonds,
                                           angles=dataframe_mol.angles,
                                           torsions=dataframe_mol.torsions,
                                           planes=dataframe_mol.planes,
                                           )
    prep_dict = prep_for_grm(dataframe_mol.sites,
                                    dataframe_mol.model,
                                    geometry_proxy_registries)


    grm = build_grm(model=dataframe_mol.model,
                    geometry_proxy_registries=geometry_proxy_registries,
                    nonbonded_types=dataframe_mol.sites.nonbonded_types,
                    nonbonded_charges=dataframe_mol.sites.nonbonded_charges,
                    prep=prep_dict
                    )
    return grm

  def build_rm_from_grm(self,grm):
    rm = mmtbx.restraints.manager(
            geometry=grm,
            cartesian_ncs_manager=None,
            normalization=False,
            use_afitt=False,
            afitt_object=None,
        )
    return rm