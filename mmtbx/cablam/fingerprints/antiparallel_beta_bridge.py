from __future__ import absolute_import, division, print_function
from mmtbx.cablam import cablam_fingerprints

#Antiparallel Beta Bridge/loose definitions, uses 4 bonds
#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
# g (h) i (j) k
# r (q) p (o) n
antiparallel_beta_bridge_close = cablam_fingerprints.motif(
  motif_name = "antiparallel_beta_bridge_close",
  residue_names = {"i":"antiparallel_beta_bridge_close", "p":"antiparallel_beta_bridge_close"})

residue1 = antiparallel_beta_bridge_close.add_residue(
  bond_move = 'p',
  index = 'i')
bond1 = residue1.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'p')
bond1.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond2 = residue1.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'p')
bond2.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue2 = antiparallel_beta_bridge_close.add_residue(
  sequence_move = 2,
  index = 'p')
bond3 = residue2.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'i')
bond3.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond4 = residue2.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'i')
bond4.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue3 = antiparallel_beta_bridge_close.add_residue(
  bond_move = 'g',
  index = 'r')
bond5 = residue3.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'g')
bond5.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue4 = antiparallel_beta_bridge_close.add_residue(
  sequence_move = 4,
  index = 'g')
bond6 = residue4.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'r')
bond6.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue5 = antiparallel_beta_bridge_close.add_residue(
  bond_move = 'n',
  index = 'k')
bond7 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'n')
bond7.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = antiparallel_beta_bridge_close.add_residue(
  sequence_move = 2,
  index = 'n')
bond8 = residue6.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'k')
bond8.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue7 = antiparallel_beta_bridge_close.add_residue(
  end_of_motif = True,
  index = 'p')
bond9 = residue7.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'i')
bond9.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond10 = residue7.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'i')
bond10.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

#-------------------------------------------------------------------------------

#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
# (g) h i j (k)
# (r) q p o (n)

antiparallel_beta_bridge_wide = cablam_fingerprints.motif(
  motif_name = "antiparallel_beta_bridge_wide",
  residue_names = {"i":"antiparallel_beta_bridge_wide", "p":"antiparallel_beta_bridge_wide"})

residue1 = antiparallel_beta_bridge_wide.add_residue(
  sequence_move = 1,
  index = 'i')

residue2 = antiparallel_beta_bridge_wide.add_residue(
  bond_move = 'o',
  index = 'j')
bond1 = residue2.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'o')
bond1.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond2 = residue2.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'o')
bond2.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue3 = antiparallel_beta_bridge_wide.add_residue(
  sequence_move = 1,
  index = 'o')
bond3 = residue3.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'j')
bond3.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond4 = residue3.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'j')
bond4.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue4 = antiparallel_beta_bridge_wide.add_residue(
  sequence_move = 1,
  index = 'p')

residue5 = antiparallel_beta_bridge_wide.add_residue(
  bond_move = 'h',
  index = 'q')
bond5 = residue5.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'h')
bond5.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond6 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'h')
bond6.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = antiparallel_beta_bridge_wide.add_residue(
  sequence_move = 1,
  index = 'h')
bond7 = residue6.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'q')
bond7.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond8 = residue6.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'q')
bond8.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue7 = antiparallel_beta_bridge_wide.add_residue(
  end_of_motif = True,
  index = 'i')
if __name__ == "__main__":
  cablam_fingerprints.make_pickle(antiparallel_beta_bridge_close)
  cablam_fingerprints.make_pickle(antiparallel_beta_bridge_wide)
