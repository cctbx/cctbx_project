from __future__ import absolute_import, division, print_function
from cctbx.sgtbx import space_group_info

#take care of R vs. H nomenclature, as well as ambiguous I-centered groups

#For every possible Patterson symmetry there is one key.
# The key is the conventional lookup symbol (CCP4 convention)
# representing the simplest chiral space group whose diffraction pattern
# has that Patterson symmetry.

# The values serve the purpose of defining those sets of space
#  groups which cannot be distinguished based on the symmetry of
#  their diffraction patterns.

sgtbx_chiral_adaptor = {
  'P1':['P1'],
  'P2':['P2'],
  'C2':['C2'],
  'P222':['P222'],
#  'C222':['C222','B222','A222'],
  'C222':['C222',],
  'I222':['I222','I212121'],
  'F222':['F222'],
  'P4':['P4'],
  'P422':['P422'],
  'I4':['I4'],
  'I422':['I422'],
  'H3':['R3'],
  'H32':['R32'],
  'P3':['P3'],
  'P321':['P321'],
  'P312':['P312'],
  'P6':['P6'],
  'P622':['P622'],
  'P23':['P23'],
  'P432':['P432'],
  'I23':['I23','I213'],
  'I432':['I432'],
  'F23':['F23'],
  'F432':['F432'],
}

sgtbx_redefine = {}
for key in sgtbx_chiral_adaptor:
  sgtbx_redefine[key] = sgtbx_chiral_adaptor[key][0]

patterson_lookup = {}
def get_patterson_group(conventional_symbol):
  global patterson_lookup
  SGI = space_group_info(sgtbx_redefine[conventional_symbol])
  if SGI in patterson_lookup:
    PG = patterson_lookup[SGI]
  else:
    SG = SGI.group()
    PG = SG.build_derived_patterson_group()
    patterson_lookup[SGI]=PG
  return PG
