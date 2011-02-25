from __future__ import division
from cctbx import adptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.development.fmt_utils import dot5fdot_list, dot6gdot, dot6gdot_list
from cctbx.development.run_shelx import FIX, NOFIX, HKLF
from iotbx.shelx.write_ins import LATT_SYMM
from libtbx.test_utils import is_below_limit
import libtbx.path
from libtbx import easy_run
from cStringIO import StringIO
import sys

def sfac_unit(lapp, xray_structure):
  reg = xray_structure.scattering_type_registry()
  unit_cell_occupancy_sums = reg.unit_cell_occupancy_sums(
    xray_structure.scatterers())
  st_sorted = sorted(reg.as_type_gaussian_dict().keys())
  lapp("     scattering types: %s" % " ".join(st_sorted))
  sfac_indices = {}
  unit = []
  for i_st,st in enumerate(st_sorted):
    sfac_indices[st] = i_st+1
    gaussian = reg.gaussian(st)
    assert gaussian.n_terms() == 4
    assert gaussian.use_c()
    a = gaussian.array_of_a()
    b = gaussian.array_of_b()
    lapp("SFAC %.6g %.6g %.6g %.6g %.6g %.6g =" % (
      a[0], b[0],
      a[1], b[1],
      a[2], b[2]))
    lapp("     %.6g %.6g %.6g 0 0 0 1" % (
      a[3], b[3], gaussian.c()))
    unit.append("%.0f" % unit_cell_occupancy_sums[reg.unique_index(st)])
  lapp("UNIT " + " ".join(unit))
  return sfac_indices

def atoms(lapp, sfac_indices, xray_structure):
  for scatterer in xray_structure.scatterers():
    st = scatterer.scattering_type
    lbl = st + str(sfac_indices[st])
    sfac = sfac_indices[st]
    coor = []
    for x in scatterer.site: coor.append(NOFIX(x))
    coor = dot5fdot_list(coor)
    sof = FIX(scatterer.weight())
    if (not scatterer.flags.use_u_aniso()):
      lapp("%-4s %d %s %s %s" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot(NOFIX(scatterer.u_iso))))
    else:
      u = adptbx.u_star_as_u_cif(xray_structure.unit_cell(), scatterer.u_star)
      u_fix = []
      for c in u: u_fix.append(NOFIX(c))
      u = u_fix
      lapp("%-4s %d %s %s %s =" % (lbl, sfac, coor, dot6gdot(sof),
        dot6gdot_list(u[:2])))
      lapp("    %s" % dot6gdot_list((u[2], u[5], u[4], u[3])))

def write_shelx76_ls(xray_structure, f_obs, titl=None, l_s_parameters="0"):
  assert xray_structure.scatterers().size() > 0
  lines = []
  lapp = lines.append
  if (titl is None):
    titl = str(xray_structure.space_group_info())
  lapp("TITL " + titl)
  lapp("CELL 0.7 " + dot6gdot_list(xray_structure.unit_cell().parameters()))
  s = StringIO()
  LATT_SYMM(s, xray_structure.space_group(), decimal=True)
  lapp(s.getvalue()[:-1])
  sfac_indices = sfac_unit(lapp, xray_structure)
  lapp("FVAR 1.")
  lapp("L.S. %s" % l_s_parameters)
  atoms(lapp, sfac_indices, xray_structure)
  HKLF(lapp, f_obs, skip_zeros=True)
  lapp("END")
  print >> open("tmp.ins", "w"), "\n".join(lines)

def run_shelx76(titl, xray_structure, f_obs):
  write_shelx76_ls(xray_structure, f_obs, titl)
  shelx_out = easy_run.fully_buffered(command="shelx76 < tmp.ins") \
    .raise_if_errors() \
    .stdout_lines
  reflections_key = "REFLEXIONS READ, OF WHICH"
  residuals_key = "RESIDUALS BEFORE CYCLE   1 FOR"
  r = None
  lines = iter(shelx_out)
  for line in lines:
    if (line.find(reflections_key) >= 0):
      flds = line.split()
      assert len(flds) == 7
      assert flds[6] == "REJECTED"
      assert flds[5] == "0"
    elif (line.find(residuals_key) >= 0):
      assert len(lines.next().strip()) == 0
      flds = lines.next().split()
      assert len(flds) == 12
      r = float(flds[2])
  if (r is None):
    raise RuntimeError("Not found in shelx76 output: %s" % residuals_key)
  assert is_below_limit(value=r, limit=0.005)

def exercise(space_group_info, d_min=1.):
  xray_structure = random_structure.xray_structure(
    space_group_info,
    elements=("N", "C", "C", "O"),
    random_u_iso=True)
  xray_structure.scattering_type_registry(table="it1992")
  f_obs = xray_structure.structure_factors(
    anomalous_flag=False,
    d_min=d_min,
    algorithm="direct",
    cos_sin_table=False).f_calc().amplitudes()
  titl = str(space_group_info)
  run_shelx76(titl, xray_structure, f_obs)

def run_call_back(flags, space_group_info):
  sg = space_group_info.group()
  if (sg.is_centric() and not sg.is_origin_centric()):
    print "Skipping space group: centre of inversion is not at origin."
    return
  exercise(space_group_info)

def run(args):
  if (libtbx.path.full_command_path(command="shelx76") is None):
    print "shelx76 not available."
    return
  debug_utils.parse_options_loop_space_groups(args, run_call_back)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
