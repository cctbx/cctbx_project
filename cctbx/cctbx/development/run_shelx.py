from cctbx.development.fmt_utils import *
from cctbx_boost import adptbx

def LATT_SYMM(SgOps):
  lines = []
  l = lines.append
  Z = SgOps.getConventionalCentringTypeSymbol()
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
  if (SgOps.isCentric()):
    if (not SgOps.isOriginCentric()):
      raise RuntimeError, "Error: " \
        + " SHELX manual: If the structure is centrosymmetric, the" \
        + " origin MUST lie on a center of symmetry."
  else:
    LATT_N = -LATT_N;
  l("LATT %d" % (LATT_N,))
  # The operator x,y,z is always assumed, so MUST NOT be input.
  for i in xrange(1, SgOps.nSMx()):
    l("SYMM %s" % (SgOps(i).as_xyz(0, 0, "XYZ", ","),))
  return lines

def calculate_cell_content(xtal):
  result = {}
  for site in xtal.Sites:
    lbl = site.CAASF().Label()
    if (not result.has_key(lbl)):
      result[lbl] = site.Occ() * site.M()
    else:
      result[lbl] += site.Occ() * site.M()
  return result

def SFAC_DISP_UNIT(xtal, short_sfac):
  import string
  lines = []
  l = lines.append
  UNIT = []
  if (short_sfac):
    celcon = calculate_cell_content(xtal)
    l("SFAC " + string.join(celcon.keys()))
    for sf in celcon.keys():
      l("DISP %s 0 0 0" % (sf,))
      UNIT.append(str(max(1, int(celcon[sf] + 0.5))))
  else:
    from cctbx_boost.eltbx.caasf_it1992 import CAASF_IT1992
    for site in xtal.Sites:
      sf = CAASF_IT1992(site.CAASF().Label())
      l("SFAC %s %.6g %.6g %.6g %.6g %.6g %.6g =" %
        (site.Label(),
         sf.a(0), sf.b(0),
         sf.a(1), sf.b(1),
         sf.a(2), sf.b(2)))
      l("     %.6g %.6g %.6g %.6g %.6g 0 1 1" %
        (sf.a(3), sf.b(3), sf.c(),
         site.fpfdp().real, site.fpfdp().imag))
      UNIT.append(str(max(1, int(site.Occ() * site.M() + 0.5))))
  l("UNIT " + string.join(UNIT))
  return lines

def NOFIX(x):
  return x

def FIX(x):
  if (x < 0.): return -10. + x
  return 10. + x

def atoms(xtal, short_sfac):
  if (short_sfac):
    celcon = calculate_cell_content(xtal).keys()
  lines = []
  l = lines.append
  i = 0
  for site in xtal.Sites:
    i += 1
    lbl = site.CAASF().Label() + str(i)
    if (short_sfac):
      sfac = celcon.index(site.CAASF().Label()) + 1
    else:
      sfac = i
    coor = []
    for x in site.Coordinates(): coor.append(NOFIX(x))
    coor = dot5fdot_list(coor)
    sof = NOFIX(site.w())
    if (not site.isAnisotropic()):
      l("%-4s %d %s %s %s" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot(NOFIX(site.Uiso()))))
    else:
      U = adptbx.Ustar_as_Uuvrs(xtal.UnitCell, site.Uaniso())
      Ufix = []
      for c in U: Ufix.append(NOFIX(c))
      U = Ufix
      l("%-4s %d %s %s %s =" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot_list(U[:2])))
      l("    %s" % dot6gdot_list((U[2], U[5], U[4], U[3])))
  return lines

def HKLF(fcalc):
  lines = []
  l = lines.append
  l("HKLF -3")
  for i in xrange(len(fcalc.H)):
    F = abs(fcalc.F[i])
    s = "%8.2f" % (F,)
    assert  len(s) == 8, "F does not fit f8.2 format."
    l("%4d%4d%4d%s%8.2f" % (fcalc.H[i] + (s, 0.01)))
  l("   0   0   0    0.00    0.00")
  return lines

def pre_check(xtal):
  if (len(xtal.Sites) > 99):
    # SHELX WPDB will mess up atom labels.
    raise RuntimeError, "Cannot handle more than 99 sites."
  for site in xtal.Sites:
    if (site.Occ() > 1.1):
      raise RuntimeError, "Error: occ too large: %s: %.6g" % (
        site.Label(), site.Occ())
    if (site.Uiso() > 1.0):
      raise RuntimeError, "Error: Uiso too large: %s: %.6g" % (
        site.Label(), site.Uiso())

def check_r1(SgInfo, miller_indices, shelx_lst):
  import string
  for l in shelx_lst:
    if (string.find(l, "R1 = ") >= 0):
      flds = string.split(l)
      R1 = string.atof(flds[9])
      nData = string.atof(flds[12])
      if (len(miller_indices) != nData):
        raise RuntimeError, "Shelx lost Miller indices."
      print "R1", R1, SgInfo.BuildLookupSymbol()
      if (R1 > 0.01):
        raise RuntimeError, "Error: " + l[:-1]
      return
  raise RuntimeError, "R1 not found in Shelx .lst file."

def check_anisou(shelx_titl, xtal, shelx_pdb):
  import string
  # SHELXL WPDB does not include H atoms. Therefore we
  # need a dictionary of labels to map back to the index
  # in the xtal.Sites list.
  lbl_dict = {}
  i = 0
  for site in xtal.Sites:
    i += 1
    lbl = string.upper(site.CAASF().Label() + str(i))
    lbl_dict[lbl] = i - 1
  TotalANISOU = 0
  TotalMismatches = 0
  for l in shelx_pdb[4:]:
    if (l[:6] == "ANISOU"):
      TotalANISOU += 1
      lbl = string.strip(l[11:16])
      i = lbl_dict[lbl]
      assert xtal.Sites[i].isAnisotropic()
      U = l[28:70]
      Ucart = adptbx.Ustar_as_Ucart(xtal.UnitCell, xtal.Sites[i].Uaniso())
      mismatch = 0
      s = ""
      for i in xrange(6):
        ushelx = string.atoi(U[i*7:(i+1)*7])
        uadptbx = int(round(Ucart[i] * 1.e+4,))
        s += "%7d" % uadptbx
        if (abs(ushelx - uadptbx) > 1): mismatch = 1
      if (mismatch != 0):
        print l[:-1]
        print U
        print s
        print "Error: ANISOU mismatch."
        TotalMismatches += 1
  print shelx_titl + (": ANISOU mismatches: %d of %d" % (
    TotalMismatches, TotalANISOU))
  assert TotalMismatches == 0

def run_shelx(shelx_titl, xtal, fcalc, short_sfac=0):
  import os, sys
  pre_check(xtal)
  lines = []
  l = lines.append
  l("TITL " + shelx_titl)
  l("CELL 1.0 " + dot6gdot_list(xtal.UnitCell.getParameters()))
  l("ZERR 1 0.01 0.01 0.01 0 0 0")
  lines += LATT_SYMM(xtal.SgOps)
  lines += SFAC_DISP_UNIT(xtal, short_sfac)
  l("FVAR 1.")
  l("L.S. 1")
  l("BLOC 0")
  l("SPEC -0.1")
  l("WPDB 2")
  lines += atoms(xtal, short_sfac)
  lines += HKLF(fcalc)
  f = open("tmp.ins", "w")
  for l in lines:
    print l
    f.write(l + "\n")
  f.close()
  sys.stdout.flush()
  sys.stderr.flush()
  try: os.unlink("tmp.lst")
  except: pass
  os.system("shelxl tmp")
  f = open("tmp.lst", "r")
  shelx_lst = f.readlines()
  f.close()
  sys.stderr.flush()
  f = open("tmp.pdb", "r")
  shelx_pdb = f.readlines()
  f.close()
  sys.stderr.flush()
  for l in shelx_lst:
    print l[:-1]
  sys.stdout.flush()
  check_r1(xtal.SgInfo, fcalc.H, shelx_lst)
  check_anisou(shelx_titl, xtal, shelx_pdb)
