from __future__ import division
from cctbx.array_family import flex

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

def rotatable(pdb_hierarchy, xray_structure, log = None):
  import mmtbx.monomer_library
  mon_lib_srv = mmtbx.monomer_library.server.server()
  from mmtbx.utils import rotatable_bonds
  hd_selection = xray_structure.hd_selection()
  sel = flex.bool(hd_selection.size(), False)
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        if(len(conformers)>1): continue
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          fr = rotatable_bonds.axes_and_atoms_aa_specific(
            residue=residue, mon_lib_srv=mon_lib_srv,
            remove_clusters_with_all_h=False, log=log)
          atoms = residue.atoms()
          if(fr is not None):
            for fr_ in fr:
              fr1 = fr_[1]
              if(len(fr1)==1 and atoms[fr1[0]].element.strip().upper() == "H"):
                #print "    ", atoms[fr_[1][0]].i_seq, \
                #  hd_selection[atoms[fr_[1][0]].i_seq], atoms[fr_[1][0]].element
                sel[atoms[fr1[0]].i_seq]=True
              if(len(fr1)==3 and atoms[fr1[0]].element.strip().upper() == "H" and
                 atoms[fr1[1]].element.strip().upper() == "H" and
                 atoms[fr1[2]].element.strip().upper() == "H"):
                #print "    ", atoms[fr_[1][0]].i_seq, \
                #  hd_selection[atoms[fr_[1][0]].i_seq], atoms[fr_[1][0]].element
                sel[atoms[fr1[0]].i_seq]=True
                sel[atoms[fr1[1]].i_seq]=True
                sel[atoms[fr1[2]].i_seq]=True
  return sel
