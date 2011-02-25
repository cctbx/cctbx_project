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
  except:
    raise RuntimeError, "Error: Lattice type not supported by SHELX."
  # N must be made negative if the structure is non-centrosymmetric.
  if (space_group.is_centric()):
    if (not space_group.is_origin_centric()):
      raise RuntimeError, "Error: " \
        + " SHELX manual: If the structure is centrosymmetric, the" \
        + " origin MUST lie on a center of symmetry."
  else:
    LATT_N = -LATT_N;
  print >> s, "LATT", LATT_N
  # The operator x,y,z is always assumed, so MUST NOT be input.
  for i in xrange(1, space_group.n_smx()):
    print >> s, "SYMM", space_group(i).as_xyz(
      decimal=decimal,
      t_first=False,
      symbol_letters="XYZ",
      separator=",")

def shelxd(s,
      title,
      crystal_symmetry,
      n_sites,
      scattering_type,
      d_min,
      mind_mdis=-3.5,
      mind_mdeq=None):
  print >> s, "TITL", title
  print >> s, "CELL 1.0 %.6g %.6g %.6g %.6g %.6g %.6g" \
                % crystal_symmetry.unit_cell().parameters()
  LATT_SYMM(s, crystal_symmetry.space_group())
  print >> s, "SFAC", scattering_type
  print >> s, "UNIT", n_sites * 4
  print >> s, "SHEL 999 %.2f" % d_min
  print >> s, "PSMF  pres -%.2f" % d_min
  print >> s, "PATS  np  100   npt     99999   nf     5"
  print >> s, "FIND", n_sites
  print >> s, "MIND", mind_mdis,
  if (mind_mdeq is not None):
    print >> s, mind_mdeq,
  print >> s
  if (n_sites >= 10): # following advice in shelx-de.pdf manual
    print >> s, "WEED 0.3"
    print >> s, "SKIP 0.5"
  print >> s, "NTRY 100"
  print >> s, "HKLF 3"
  print >> s, "END"
