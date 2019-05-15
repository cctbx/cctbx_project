from __future__ import absolute_import, division, print_function
from six.moves import range
def LATT_SYMM(s, space_group, decimal=False):
  Z = space_group.conventional_centring_type_symbol()
  Z_dict = {
    "P": 1,
    "I": 2,
    "R": 3,
    "F": 4,
    "A": 5,
    "B": 6,
    "C": 7,
  }
  try:
    LATT_N = Z_dict[Z]
  except Exception:
    raise RuntimeError("Error: Lattice type not supported by SHELX.")
  # N must be made negative if the structure is non-centrosymmetric.
  if (space_group.is_centric()):
    if (not space_group.is_origin_centric()):
      raise RuntimeError("Error: " \
        + " SHELX manual: If the structure is centrosymmetric, the" \
        + " origin MUST lie on a center of symmetry.")
  else:
    LATT_N = -LATT_N;
  print("LATT", LATT_N, file=s)
  # The operator x,y,z is always assumed, so MUST NOT be input.
  for i in range(1, space_group.n_smx()):
    print("SYMM", space_group(i).as_xyz(
      decimal=decimal,
      t_first=False,
      symbol_letters="XYZ",
      separator=","), file=s)

def shelxd(s,
      title,
      crystal_symmetry,
      n_sites,
      scattering_type,
      d_min,
      mind_mdis=-3.5,
      mind_mdeq=None):
  print("TITL", title, file=s)
  print("CELL 1.0 %.6g %.6g %.6g %.6g %.6g %.6g" \
                % crystal_symmetry.unit_cell().parameters(), file=s)
  LATT_SYMM(s, crystal_symmetry.space_group())
  print("SFAC", scattering_type, file=s)
  print("UNIT", n_sites * 4, file=s)
  print("SHEL 999 %.2f" % d_min, file=s)
  print("PSMF  pres -%.2f" % d_min, file=s)
  print("PATS  np  100   npt     99999   nf     5", file=s)
  print("FIND", n_sites, file=s)
  print("MIND", mind_mdis, end=' ', file=s)
  if (mind_mdeq is not None):
    print(mind_mdeq, end=' ', file=s)
  print(file=s)
  if (n_sites >= 10): # following advice in shelx-de.pdf manual
    print("WEED 0.3", file=s)
    print("SKIP 0.5", file=s)
  print("NTRY 100", file=s)
  print("HKLF 3", file=s)
  print("END", file=s)
