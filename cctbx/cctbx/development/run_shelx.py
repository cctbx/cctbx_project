from cctbx.web.shelx import LATT_SYMM
from cctbx import adptbx
from cctbx.eltbx.caasf import it1992
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.development.fmt_utils import *
from scitbx.python_utils import dicts
import sys, os

def check_shelx_availability():
  shelxl_out = []
  try:
    f = os.popen("shelxl", "r")
    shelxl_out = f.readlines()
    f.close()
  except:
    pass
  if (len(shelxl_out) == 0):
    print "SHELX not available."
    sys.exit(1)

def calculate_cell_content(xray_structure):
  result = dicts.with_default_value(0)
  for sc in xray_structure.scatterers():
    result[sc.caasf.label()] += sc.occupancy * sc.multiplicity()
  return result

def SFAC_DISP_UNIT(xray_structure, short_sfac):
  lines = []
  l = lines.append
  UNIT = []
  if (short_sfac):
    celcon = calculate_cell_content(xray_structure)
    l("SFAC " + " ".join(celcon.keys()))
    for sf in celcon.keys():
      l("DISP %s 0 0 0" % (sf,))
      UNIT.append(str(max(1, int(celcon[sf] + 0.5))))
  else:
    for scatterer in xray_structure.scatterers():
      caasf = it1992(scatterer.caasf.label())
      a = caasf.a()
      b = caasf.b()
      l("SFAC %s %.6g %.6g %.6g %.6g %.6g %.6g =" %
        (scatterer.label,
         a[0], b[0],
         a[1], b[1],
         a[2], b[2]))
      l("     %.6g %.6g %.6g %.6g %.6g 0 1 1" %
        (a[3], b[3], caasf.c(),
         scatterer.fp_fdp.real, scatterer.fp_fdp.imag))
      UNIT.append(
        str(max(1, int(scatterer.occupancy * scatterer.multiplicity() + 0.5))))
  l("UNIT " + " ".join(UNIT))
  return lines

def NOFIX(x):
  return x

def FIX(x):
  if (x < 0.): return -10. + x
  return 10. + x

def atoms(xray_structure, short_sfac):
  if (short_sfac):
    celcon = calculate_cell_content(xray_structure).keys()
  lines = []
  l = lines.append
  i = 0
  for scatterer in xray_structure.scatterers():
    i += 1
    lbl = scatterer.caasf.label() + str(i)
    if (short_sfac):
      sfac = celcon.index(scatterer.caasf.label()) + 1
    else:
      sfac = i
    coor = []
    for x in scatterer.site: coor.append(NOFIX(x))
    coor = dot5fdot_list(coor)
    sof = NOFIX(scatterer.weight())
    if (not scatterer.anisotropic_flag):
      l("%-4s %d %s %s %s" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot(NOFIX(scatterer.u_iso))))
    else:
      u = adptbx.u_star_as_u_cif(xray_structure.unit_cell(), scatterer.u_star)
      u_fix = []
      for c in u: u_fix.append(NOFIX(c))
      u = u_fix
      l("%-4s %d %s %s %s =" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot_list(u[:2])))
      l("    %s" % dot6gdot_list((u[2], u[5], u[4], u[3])))
  return lines

def HKLF(fcalc_array):
  lines = []
  l = lines.append
  l("HKLF -3")
  for i,h in fcalc_array.indices().items():
    f = abs(fcalc_array.data()[i])
    s = "%8.2f" % (f,)
    assert  len(s) == 8, "structure factor does not fit f8.2 format."
    l("%4d%4d%4d%s%8.2f" % (h + (s, 0.01)))
  l("   0   0   0    0.00    0.00")
  return lines

def pre_check(xray_structure):
  if (len(xray_structure.scatterers()) > 99):
    # SHELX WPDB will mess up atom labels.
    raise RuntimeError, "Cannot handle more than 99 scatterer."
  for scatterer in xray_structure.scatterers():
    if (scatterer.occupancy > 1.1):
      raise RuntimeError, "Error: occupancy too large: %s: %.6g" % (
        scatterer.label, scatterer.occupancy)
    if (scatterer.u_iso > 1.0):
      raise RuntimeError, "Error: u_iso too large: %s: %.6g" % (
        scatterer.label, scatterer.u_iso)

def check_r1(miller_set, shelx_lst, verbose):
  for l in shelx_lst:
    if (l.find("R1 = ") >= 0):
      flds = l.split()
      R1 = float(flds[9])
      n_data = int(flds[12])
      if (len(miller_set.indices()) != n_data):
        raise RuntimeError, "Shelx lost Miller indices."
      if (0 or verbose):
        print "R1", R1, miller_set.space_group_info()
      if (R1 > 0.01):
        raise RuntimeError, "Error: " + l[:-1]
      return
  raise RuntimeError, "R1 not found in Shelx .lst file."

def check_anisou(shelx_titl, xray_structure, shelx_pdb, verbose):
  # SHELXL WPDB does not include H atoms. Therefore we
  # need a dictionary of labels to map back to the index
  # in the xray_structure.scatterers() list.
  lbl_dict = {}
  i = 0
  for scatterer in xray_structure.scatterers():
    i += 1
    lbl = (scatterer.caasf.label() + str(i)).upper()
    lbl_dict[lbl] = i - 1
  TotalANISOU = 0
  TotalMismatches = 0
  for l in shelx_pdb[4:]:
    if (l[:6] == "ANISOU"):
      TotalANISOU += 1
      lbl = l[11:16].strip()
      i = lbl_dict[lbl]
      assert xray_structure.scatterers()[i].anisotropic_flag
      u = l[28:70]
      u_cart = adptbx.u_star_as_u_cart(
        xray_structure.unit_cell(), xray_structure.scatterers()[i].u_star)
      mismatch = 0
      s = ""
      for i in xrange(6):
        u_shelx = int(u[i*7:(i+1)*7])
        u_adptbx = int(round(u_cart[i] * 1.e+4,))
        s += "%7d" % u_adptbx
        if (abs(u_shelx - u_adptbx) > 1): mismatch = 1
      if (mismatch != 0):
        print l[:-1]
        print u
        print s
        print "Error: ANISOU mismatch."
        TotalMismatches += 1
  if (0 or verbose or TotalMismatches > 0):
    print shelx_titl + (": ANISOU mismatches: %d of %d" % (
      TotalMismatches, TotalANISOU))
  assert TotalMismatches == 0

def run_shelx(shelx_titl, structure_factors, short_sfac=False, verbose=0):
  xray_structure = structure_factors.xray_structure()
  assert xray_structure.scatterers().size() > 0
  pre_check(xray_structure)
  f_calc_array = structure_factors.f_calc_array()
  lines = []
  l = lines.append
  l("TITL " + shelx_titl)
  l("CELL 1.0 " + dot6gdot_list(xray_structure.unit_cell().parameters()))
  l("ZERR 1 0.01 0.01 0.01 0 0 0")
  lines += LATT_SYMM(xray_structure.space_group())
  lines += SFAC_DISP_UNIT(xray_structure, short_sfac)
  l("FVAR 1.")
  l("L.S. 1")
  l("BLOC 0")
  l("SPEC -0.1")
  l("WPDB 2")
  lines += atoms(xray_structure, short_sfac)
  lines += HKLF(f_calc_array)
  f = open("tmp.ins", "w")
  for l in lines:
    if (0 or verbose): print l
    f.write(l + "\n")
  f.close()
  sys.stdout.flush()
  sys.stderr.flush()
  try: os.unlink("tmp.lst")
  except: pass
  f = os.popen("shelxl tmp", "r")
  shelx_out = f.readlines()
  f.close()
  if (0 or verbose):
    for l in shelx_out: print l[:-1]
  f = open("tmp.lst", "r")
  shelx_lst = f.readlines()
  f.close()
  sys.stderr.flush()
  f = open("tmp.pdb", "r")
  shelx_pdb = f.readlines()
  f.close()
  sys.stderr.flush()
  if (0 or verbose):
    for l in shelx_lst: print l[:-1]
  sys.stdout.flush()
  check_r1(f_calc_array, shelx_lst, verbose)
  check_anisou(shelx_titl, xray_structure, shelx_pdb, verbose)

def exercise(space_group_info,
             anomalous_flag=False,
             anisotropic_flag=False,
             d_min=2.,
             verbose=0):
  structure_factors = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    anisotropic_flag=anisotropic_flag,
    random_f_prime_d_min=1.0,
    random_f_double_prime=anomalous_flag
    ).structure_factors_direct(
        anomalous_flag=anomalous_flag, d_min=d_min)
  if (0 or verbose):
    structure_factors.xray_structure().show_summary()
  shelx_titl = str(space_group_info) \
             + ", anomalous=" + str(anomalous_flag) \
             + ", anisotropic_flag=" + str(anisotropic_flag)
  run_shelx(shelx_titl, structure_factors, verbose=verbose)

def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    for anisotropic_flag in (False, True):
      exercise(space_group_info, anomalous_flag, anisotropic_flag,
               verbose=flags.Verbose)

def run():
  check_shelx_availability()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
