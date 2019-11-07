from __future__ import absolute_import, division, print_function
from mmtbx.cablam import cablam_fingerprints

#Antiparallel beta, close
#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
# g (h)* i* (j)* k
# r (q)* p* (o)* n

antiparallel_beta_wcw = cablam_fingerprints.motif(
  motif_name = "antiparallel_beta_wcw",
  residue_names = {"i":"antiparallel_beta_close", "p":"antiparallel_beta_close"})

residue1 = antiparallel_beta_wcw.add_residue(
  bond_move = 'p',
  end_of_motif = False,
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

residue2 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
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

residue3 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'q')

residue4 = antiparallel_beta_wcw.add_residue(
  bond_move = 'g',
  end_of_motif = False,
  index = 'r')
bond5 = residue4.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'g')
bond5.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond6 = residue4.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'g')
bond6.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue5 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'g')
bond7 = residue5.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'r')
bond7.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond8 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'r')
bond8.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = antiparallel_beta_wcw.add_residue(
  sequence_move = 2,
  end_of_motif = False,
  index = 'h')

residue7 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'j')

residue8 = antiparallel_beta_wcw.add_residue(
  bond_move = 'n',
  end_of_motif = False,
  index = 'k')
bond9 = residue8.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'n')
bond9.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond10 = residue8.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'n')
bond10.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue9 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'n')
bond11 = residue9.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'k')
bond11.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond12 = residue9.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'k')
bond12.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue10 = antiparallel_beta_wcw.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'o')
residue11 = antiparallel_beta_wcw.add_residue(
  end_of_motif = True,
  index = 'p')
bond13 = residue11.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'i')
bond13.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond14 = residue11.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'i')
bond14.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

#-------------------------------------------------------------------------------

#Antiparallel beta wide
#Original by Christopher Williams, converted to new format by Danny Oh
#Two strands:
# e (f) g* (h)* i* (j) k
# t (s) r* (q)* p* (o) n

antiparallel_beta_cwc = cablam_fingerprints.motif(
  motif_name = "antiparallel_beta_cwc",
  residue_names = {"q":"antiparallel_beta_wide", "h":"antiparallel_beta_wide"})

residue1 = antiparallel_beta_cwc.add_residue(
  bond_move = 'p',
  end_of_motif = False,
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

residue2 = antiparallel_beta_cwc.add_residue(
  sequence_move = 1,
  end_of_motif = False,
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

residue3 = antiparallel_beta_cwc.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'q')


residue4 = antiparallel_beta_cwc.add_residue(
  bond_move = 'g',
  end_of_motif = False,
  index = 'r')
bond5 = residue4.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'g')
bond5.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)
bond6 = residue4.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'g')
bond6.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue5 = antiparallel_beta_cwc.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'g')
bond7 = residue5.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'r')
bond7.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond8 = residue5.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'r')
bond8.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue6 = antiparallel_beta_cwc.add_residue(
  sequence_move = 1,
  end_of_motif = False,
  index = 'h')

residue7 = antiparallel_beta_cwc.add_residue(
  sequence_move = 2,
  end_of_motif = False,
  index = 'i')
bond9 = residue7.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'p')
bond9.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)
bond10 = residue7.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'p')
bond10.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue8 = antiparallel_beta_cwc.add_residue(
  bond_move = 'n',
  end_of_motif = False,
  index = 'k')
bond11 = residue8.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'n')
bond11.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue9 = antiparallel_beta_cwc.add_residue(
  sequence_move = 6,
  end_of_motif = False,
  index = 'n')
bond12 = residue9.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 'k')
bond12.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue10 = antiparallel_beta_cwc.add_residue(
  bond_move = 'e',
  end_of_motif = False,
  index = 't')
bond13 = residue10.add_bond(
  required = True,
  src_atom = ' H  ',
  trg_index = 'e')
bond13.add_target_atom(
  atomname = ' O  ',
  anyseqdist = True)

residue11 = antiparallel_beta_cwc.add_residue(
  sequence_move = 2,
  end_of_motif = False,
  index = 'e')
bond14 = residue11.add_bond(
  required = True,
  src_atom = ' O  ',
  trg_index = 't')
bond14.add_target_atom(
  atomname = ' H  ',
  anyseqdist = True)

residue12 = antiparallel_beta_cwc.add_residue(
  end_of_motif = True,
  index = 'g')

if __name__ == "__main__":
  cablam_fingerprints.make_pickle(antiparallel_beta_wcw)
  cablam_fingerprints.make_pickle(antiparallel_beta_cwc)
