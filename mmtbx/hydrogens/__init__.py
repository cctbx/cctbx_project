from __future__ import division
from mmtbx.utils import rotatable_bonds

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
            r = analyze_group(g=fr, atoms=atoms)
            if(r is not None): result.append(r)
  return result
