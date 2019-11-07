from __future__ import absolute_import, division, print_function
from mmtbx.cablam import cablam_fingerprints

#Loose helix definitions, n-terminal end

alpha_helix_3os = cablam_fingerprints.motif(
  motif_name = 'alpha_helix_3os',
  residue_names = {'b':'alpha_helix_3os'})
res1 = alpha_helix_3os.add_residue(
  sequence_move = 1)
bond1 = res1.add_bond(
  src_atom = ' O  ')
bond1.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)

res2 = alpha_helix_3os.add_residue(
  sequence_move = 1,
  index = 'b')
bond2 = res2.add_bond(
  src_atom = ' O  ')
bond2.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)

res3 = alpha_helix_3os.add_residue(
  end_of_motif = True)
bond3 = res3.add_bond(
  src_atom = ' O  ')
bond3.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)

#-------------------------------------------------------------------------------

#Loose helix definitions, c-terminal end

alpha_helix_3hs = cablam_fingerprints.motif(
  motif_name = 'alpha_helix_3hs',
  residue_names = {'b':'alpha_helix_3hs'})
res1 = alpha_helix_3hs.add_residue(
  sequence_move = 1)
bond1 = res1.add_bond(
  src_atom = ' H  ')
bond1.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

res2 = alpha_helix_3hs.add_residue(
  sequence_move = 1,
  index = 'b')
bond2 = res2.add_bond(
  src_atom = ' H  ')
bond2.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

res3 = alpha_helix_3hs.add_residue(
  end_of_motif = True)
bond3 = res3.add_bond(
  src_atom = ' H  ')
bond3.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

#-------------------------------------------------------------------------------

#Regular helix definition

alpha_helix_3_full = cablam_fingerprints.motif(
  motif_name='alpha_helix_3_full',
  residue_names = {'b':'alpha_helix_3_full'})
res1 = alpha_helix_3_full.add_residue(
  sequence_move = 1)
bond1 = res1.add_bond(
  src_atom = ' O  ')
bond1.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)
bond2 = res1.add_bond(
  src_atom = ' H  ')
bond2.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

res2 = alpha_helix_3_full.add_residue(
  sequence_move = 1,
  index = 'b')
bond3 = res2.add_bond(
  src_atom = ' O  ')
bond3.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)
bond4 = res2.add_bond(
  src_atom = ' H  ')
bond4.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

res3 = alpha_helix_3_full.add_residue(
  end_of_motif = True)
bond5 = res3.add_bond(
  src_atom = ' O  ')
bond5.add_target_atom(
  atomname = ' H  ',
  seqdist = 4)
bond6 = res3.add_bond(
  src_atom = ' H  ')
bond6.add_target_atom(
  atomname = ' O  ',
  seqdist = -4)

#-------------------------------------------------------------------------------

threeten_general = cablam_fingerprints.motif(
  motif_name = 'threeten_general',
  residue_names = {'a':'threeten_gen_n_end','b':'threeten_gen_n','c':'threeten_gen_center','d':'threeten_gen_c','e':'threeten_gen_c_end'})
res1 = threeten_general.add_residue(
  sequence_move = 1,
  index = 'a')
bond1 = res1.add_bond(
  src_atom = ' O  ')
bond1.add_target_atom(
  atomname = ' H  ',
  seqdist = 3)

res2 = threeten_general.add_residue(
  sequence_move = 1,
  index = 'b')
bond2 = res2.add_bond(
  src_atom = ' O  ')
bond2.add_target_atom(
  atomname = ' H  ',
  seqdist = 3)

res3 = threeten_general.add_residue(
  sequence_move = 1,
  index = 'c')
bond3a = res3.add_bond(
  allow_bifurcated = True,
  src_atom = ' O  ')
bond3a.add_target_atom(
  atomname = ' H  ',
  seqdist = 3)
bond3b = res3.add_bond(
  allow_bifurcated = True,
  src_atom = ' H  ')
bond3b.add_target_atom(
  atomname = ' O  ',
  seqdist = -3)

res4 = threeten_general.add_residue(
  sequence_move = 1,
  index = 'd')
bond4 = res4.add_bond(
  src_atom = ' H  ')
bond4.add_target_atom(
  atomname = ' O  ',
  seqdist = -3)

res5 = threeten_general.add_residue(
  sequence_move = 1,
  index = 'e',
  end_of_motif = True)
bond5 = res5.add_bond(
  src_atom = ' H  ')
bond5.add_target_atom(
  atomname = ' O  ',
  seqdist = -3)
if __name__ == "__main__":
  cablam_fingerprints.make_pickle(alpha_helix_3os)
  cablam_fingerprints.make_pickle(alpha_helix_3hs)
  cablam_fingerprints.make_pickle(alpha_helix_3_full)
  cablam_fingerprints.make_pickle(threeten_general)
