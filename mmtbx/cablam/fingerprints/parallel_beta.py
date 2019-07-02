from __future__ import absolute_import, division, print_function
from mmtbx.cablam import cablam_fingerprints

#Parallel beta
#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
#   (g) h* i* j* (k)
# m (n) o* p* q* (r) s

parallel_beta = cablam_fingerprints.motif(
  motif_name = "parallel_beta",
  residue_names = {"i":"parallel_beta_close", "p":"parallel_beta_wide"})

residue1 = parallel_beta.add_residue(
  sequence_move = 1,
  index = 'i')
bond1 = residue1.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'q')
bond1.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond2 = residue1.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'o')
bond2.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue2 = parallel_beta.add_residue(
  sequence_move = 1,
  index = 'j')

residue3 = parallel_beta.add_residue(
  bond_move = 's',
  index = 'k')
bond3 = residue3.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'q')
bond3.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)
bond4 = residue3.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 's')
bond4.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue4 = parallel_beta.add_residue(
  sequence_move = -2,
  index = 's')
bond5 = residue4.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'k')
bond5.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue5 = parallel_beta.add_residue(
  sequence_move = -1,
  index = 'q')
bond6 = residue5.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'k')
bond6.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond7 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'i')
bond7.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = parallel_beta.add_residue(
  sequence_move = -1,
  index = 'p')

residue7 = parallel_beta.add_residue(
  sequence_move = -2,
  index = 'o')
bond8 = residue7.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'i')
bond8.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond9 = residue7.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'g')
bond9.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue8 = parallel_beta.add_residue(
  bond_move = 'g',
  index = 'm')
bond10 = residue8.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'g')
bond10.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue9 = parallel_beta.add_residue(
  sequence_move = 1,
  index = 'g')
bond11 = residue9.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'o')
bond11.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond12 = residue9.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'm')
bond12.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue10 = parallel_beta.add_residue(
  sequence_move = 1,
  index = 'h')

residue11 = parallel_beta.add_residue(
  end_of_motif = True,
  index = 'i')

#-------------------------------------------------------------------------------

#Parallel beta bridges
#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
# (g) h i j (k)
# (n) o p q (r)

parallel_beta_bridge = cablam_fingerprints.motif(
  motif_name = "parallel_beta_bridge",
  residue_names = {"i":"parallel_beta_bridge_close", "p":"parallel_beta_bridge_wide"})

residue1 = parallel_beta_bridge.add_residue(
  sequence_move = 2,
  index = 'i')
bond1 = residue1.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'q')
bond1.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond2 = residue1.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'o')
bond2.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue2 = parallel_beta_bridge.add_residue(
  bond_move = 'q',
  index = 'k')
bond3 = residue2.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'q')
bond3.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue3 = parallel_beta_bridge.add_residue(
  sequence_move = -1,
  index = 'q')
bond4 = residue3.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'k')
bond4.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond5 = residue3.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'i')
bond5.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue4 = parallel_beta_bridge.add_residue(
  sequence_move = -1,
  index = 'p')

residue5 = parallel_beta_bridge.add_residue(
  bond_move = 'g',
  index = 'o')
bond6 = residue5.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'i')
bond6.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond7 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'g')
bond7.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = parallel_beta_bridge.add_residue(
  sequence_move = 2,
  index = 'g')
bond8 = residue6.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'o')
bond8.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue7 = parallel_beta_bridge.add_residue(
  end_of_motif = True,
  index = 'i')

if __name__ == "__main__":
  cablam_fingerprints.make_pickle(parallel_beta)
  cablam_fingerprints.make_pickle(parallel_beta_bridge)
