from __future__ import absolute_import, division, print_function
from six.moves import range
def site_constraints_special_op_simplified(O, fvars, site, p_tolerance):
  assert len(fvars) > 0
  new_fvars = []
  def add_fvar(fv):
    new_fvars.append(fv)
    return 10 * (len(fvars) + len(new_fvars))
  coded_variables = [None]*3
  i_var_i_terms = [[], [], []]
  for i_term,term in enumerate(O.terms):
    if (len(term.i_vars) == 0):
      c = term.constant
      if (c < -4 or c > 4): c %= 1
      coded_variables[i_term] = 10 + float(c)
      continue
    if (len(term.i_vars) > 1):
      return None
    i_var_i_terms[term.i_vars[0]].append(i_term)
  def check_p(array):
    for p in array:
      if (p_tolerance > abs(p) > 5-p_tolerance):
        return None
  def set_cv(i_term, value):
    assert coded_variables[i_term] is None
    coded_variables[i_term] = value
  for i_var in range(3):
    i_terms = i_var_i_terms[i_var]
    if (len(i_terms) == 1):
      i_term = i_terms[0]
      set_cv(i_term, site[i_term])
    elif (len(i_terms) == 2):
      i0, i1 = i_terms
      t0, t1 = [O.terms[_] for _ in i_terms]
      assert t0.is_identity()
      m1, c1 = t1.multipliers[0], t1.constant
      if (c1 == 0):
        if (abs(m1) <= 1):
          p0 = 1
          p1 = float(m1)
        else:
          p0 = float(1 / m1)
          p1 = 1
        s1 = 1
      else:
        p0 = float(-c1 / m1)
        p1 = float(-c1)
        s1 = -1
      check_p([p0, p1])
      fv = site[i0] / p0
      m10 = add_fvar(fv)
      set_cv(i0,      m10 + p0)
      set_cv(i1, s1 * m10 + p1)
    elif (len(i_terms) == 3):
      assert i_terms == [0,1,2]
      t0, t1, t2 = O.terms
      t0.is_identity()
      m1, c1 = t1.multipliers[0], t1.constant
      m2, c2 = t2.multipliers[0], t2.constant
      if (c1 != 0):
        if (c1 != 0 and c2 != 0):
          if (c1 / m1 != c2 / m2):
            return None
          p0 = float(-c1 / m1)
          p1 = float(-c1)
          p2 = float(-c2)
          s1, s2 = -1, -1
        else:
          p0 = float(-c1 / m1)
          p1 = float(-c1)
          p2 = float(m2) * p0
          s1, s2 = -1, 1
      elif (c2 != 0):
        p0 = float(-c2 / m2)
        p1 = float(m1) * p0
        p2 = float(-c2)
        s1, s2 = 1, -1
      else:
        p0 = 1 / float(max(1, abs(m1), abs(m2)))
        p1 = float(m1) * p0
        p2 = float(m2) * p0
        s1, s2 = 1, 1
      check_p([p0, p1, p2])
      fv = site[0] / p0
      m10 = add_fvar(fv)
      set_cv(0,      m10 + p0)
      set_cv(1, s1 * m10 + p1)
      set_cv(2, s2 * m10 + p2)
  assert coded_variables.count(None) == 0
  fvars.extend(new_fvars)
  return coded_variables

def raise_special_position_constraints_cannot_be_encoded(sos):
  from iotbx.shelx.errors import error
  raise error(
    "The special position constraints %s cannot be encoded using the"
    " SHELX FVAR mechanism." % str(sos), line=None)

def site_constraints_site_symmetry_ops(O, fvars, site, p_tolerance):
  sos = O.special_op_simplified()
  result = sos.shelx_fvar_encoding(
    fvars=fvars, site=site, p_tolerance=p_tolerance)
  if (result is None):
    raise_special_position_constraints_cannot_be_encoded(sos)
  return result

unit_shifts_prioritized = (
  (0,0,0),
  (0,0,1),
  (0,1,0),
  (1,0,0),
  (0,0,-1),
  (0,-1,0),
  (-1,0,0),
  (1,1,0),
  (0,1,1),
  (1,0,1),
  (1,-1,0),
  (0,1,-1),
  (1,0,-1),
  (-1,1,0),
  (0,-1,1),
  (-1,0,1),
  (-1,-1,0),
  (0,-1,-1),
  (-1,0,-1),
  (1,1,1),
  (1,1,-1),
  (1,-1,1),
  (-1,1,1),
  (1,-1,-1),
  (-1,-1,1),
  (-1,1,-1),
  (-1,-1,-1))

def move_sites_if_necessary_for_shelx_fvar_encoding(xray_structure):
  from cctbx import xray
  from scitbx import matrix
  xs = xray_structure
  scs = xs.scatterers().deep_copy()
  sstab = xs.site_symmetry_table()
  for i_sc in xs.special_position_indices():
    site = scs[i_sc].site
    ss = sstab.get(i_sc)
    sos = sos_orig = ss.special_op_simplified()
    fvars = [None]
    coded_variables = sos.shelx_fvar_encoding(site=site, fvars=fvars)
    if (coded_variables is None):
      site_0 = matrix.col(site).each_mod_short()
      for u in unit_shifts_prioritized:
        site = site_0 + matrix.col(u)
        ss = xs.site_symmetry(site)
        sos = ss.special_op_simplified()
        coded_variables = sos.shelx_fvar_encoding(site=site, fvars=fvars)
        if (coded_variables is not None):
          scs[i_sc].site = site
          break
      else:
        raise_special_position_constraints_cannot_be_encoded(sos_orig)
  result = xray.structure(
    special_position_settings=xs,
    scatterers=scs,
    non_unit_occupancy_implies_min_distance_sym_equiv_zero
      =xs.non_unit_occupancy_implies_min_distance_sym_equiv_zero(),
    scattering_type_registry=xs.scattering_type_registry())
  assert result.special_position_indices().all_eq(xs.special_position_indices())
  return result

def dev_build_shelx76_fvars(xray_structure, validation_tolerance=1e-10):
  from iotbx.shelx.parsers import decode_variables
  fvars = [1]
  encoded_sites = []
  scs = xray_structure.scatterers()
  sstab = xray_structure.site_symmetry_table()
  for i_sc in xray_structure.special_position_indices():
    site = scs[i_sc].site
    ss = sstab.get(i_sc)
    coded_variables = ss.shelx_fvar_encoding(site=site, fvars=fvars)
    if (validation_tolerance is not None):
      values, behaviors = decode_variables(
        free_variable=fvars, coded_variables=coded_variables)
      mismatch = xray_structure.unit_cell().distance(site, values)
      assert mismatch < validation_tolerance
    encoded_sites.append(coded_variables)
  return fvars, encoded_sites
