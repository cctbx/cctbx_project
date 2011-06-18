from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from cmath import cos, sin, pi


def miller_export_as_shelx_fcf(self, f_calc, file_object=None):
  """ Export self and the miller array f_calc as ShelX would do with
  the instruction LIST 6
  """
  assert self.is_real_array()
  assert f_calc.is_complex_array()
  assert self.indices().all_eq(f_calc.indices())
  assert self.anomalous_flag() is f_calc.anomalous_flag()
  if (file_object is None): file_object = sys.stdout
  fo = self.data()
  fc = f_calc.data()
  f = file_object
  print >> f, """
loop
_symmetry_equiv_pos_as_xyz
"""
  for op in self.space_group():
    print >> f, "'%s'" % op.as_xyz()

  print >> f
  cell_labels = [ "_cell_length_%s" % s for s in ('a', 'b', 'c') ]\
              + [ "_cell_angle_%s" % s for s in ('alpha', 'beta', 'gamma') ]
  for lbl, param in zip(cell_labels, self.unit_cell().parameters()):
    print >> f, "%s\t%f" % param

  print >> f, """
loop_
 _refln_index_h
 _refln_index_k
 _refln_index_l
 _refln_F_squared_meas
 _refln_F_squared_sigma
 _refln_F_calc
 _refln_phase_calc
"""

miller.array.export_as_shelx_fcf = miller_export_as_shelx_fcf


def list_6_as_miller_arrays(file_name):
  """ Read the file of given name and return a pair of miller arrays
  (F_obs^2, F_cal) """
  # potentially iotbx.cif could be used here
  fcf = iter(open(file_name))
  space_group = sgtbx.space_group()
  unit_cell_params = {}
  indices = flex.miller_index()
  f_obs_squares = flex.double()
  sigma_f_obs_squares = flex.double()
  f_calc_amplitudes = flex.double()
  f_calc_phases = flex.double()
  for li in fcf:
    if li.startswith('loop_'):
      for li in fcf:
        li = li.strip()
        if li == '_symmetry_equiv_pos_as_xyz':
          for li in fcf:
            li = li.strip()
            if not li: break
            space_group.expand_smx(li[1:-1])
        else:
          for i in xrange(6): fcf.next()
          for li in fcf:
            items = li.split()
            if not items: break
            h,k,l, fo, sig_fo, fc, phase = items
            indices.append((int(h), int(k), int(l)))
            f_obs_squares.append(float(fo))
            sigma_f_obs_squares.append(float(sig_fo))
            f_calc_amplitudes.append(float(fc))
            f_calc_phases.append(float(phase))
        if not li: break
    elif li.startswith('_cell'):
      lbl, value = li.split()
      unit_cell_params[lbl] = float(value)

  unit_cell = uctbx.unit_cell(
    [ unit_cell_params[p]
      for p in ( "_cell_length_a","_cell_length_b","_cell_length_c",
                 "_cell_angle_alpha","_cell_angle_beta","_cell_angle_gamma" )
    ])
  crystal_symmetry = crystal.symmetry(
    unit_cell=unit_cell,
    space_group=space_group)
  f_calc_phases *= pi/180
  f_calc = flex.complex_double(
    reals=f_calc_amplitudes * flex.cos(f_calc_phases),
    imags=f_calc_amplitudes * flex.sin(f_calc_phases) )
  miller_set = miller.set(
      crystal_symmetry=crystal_symmetry,
      indices=indices).auto_anomalous()
  f_obs_squares = miller.array(
    miller_set=miller_set,
    data=f_obs_squares,
    sigmas=sigma_f_obs_squares)
  f_obs_squares.set_observation_type_xray_intensity()
  f_obs_squares.set_info(miller.array_info(
    source=file_name,
    labels=["F_squared_meas", "F_squared_sigma"]))
  f_calc = miller.array(
    miller_set=miller_set,
    data=f_calc)
  f_obs_squares.set_info(miller.array_info(
    source=file_name,
    labels=["F_calc"]))
  return f_obs_squares, f_calc
