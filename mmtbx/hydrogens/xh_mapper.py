"""
Map between the nuclear and the X-ray positions of H and D atoms.

With the X-ray (n_gaussian, wk1995, it1992, ...) scattering table the
electron density of H is centered on the electron cloud, which sits closer to
the parent atom than the nucleus. The mapper moves every H or D bonded to one
heavy atom along its bond to either the X-ray or the nuclear ideal length.
"""
import copy
from cctbx.array_family import flex
from libtbx.utils import Sorry, null_out

class mapper(object):
  """
  X-H bond length mapper.

  Parameters
  ----------
  model : mmtbx.model.manager
    Processed model with restraints. Not modified.
  scattering_table : str, optional
    Scattering table; the model's own one when None.

  Attributes
  ----------
  is_identity : bool
    True for the neutron and electron tables: nothing is mapped.
  i_h, i_x : flex.size_t
    Indices of the mapped H/D and of their parent atoms, i_h increasing.
  d_xray, d_nuclear : flex.double
    Ideal X-ray and nuclear X-H lengths, one per mapped H.
  n_atoms : int
    Number of atoms in the model.
  """

  def __init__(self, model, scattering_table=None):
    if scattering_table is None:
      scattering_table = model.get_scattering_table()
    if scattering_table is None:
      raise Sorry("scattering_table must be set")
    self.is_identity = scattering_table in ("neutron", "electron")
    self.n_atoms = model.size()
    self.i_h = flex.size_t()
    self.i_x = flex.size_t()
    self.d_xray = flex.double()
    self.d_nuclear = flex.double()
    if self.is_identity:
      return
    pip = model.get_current_pdb_interpretation_params()
    neutron_first = pip.pdb_interpretation.use_neutron_distances
    own = self._table(model)
    other_model = model.deep_copy()
    other_pip = copy.deepcopy(pip)
    other_pip.pdb_interpretation.use_neutron_distances = not neutron_first
    other_model.set_log(null_out())
    other_model.process(make_restraints=True,
      pdb_interpretation_params=other_pip)
    other = self._table(other_model)
    if neutron_first:
      nuclear, xray = own, other
    else:
      nuclear, xray = other, own
    assert sorted(nuclear.keys()) == sorted(xray.keys()), (
      "X-H tables differ: H i_seqs only in one of them: %s" % str(
        sorted(set(nuclear.keys()).symmetric_difference(xray.keys()))))
    for i_h in sorted(nuclear.keys()):
      assert nuclear[i_h][0] == xray[i_h][0], (
        "H i_seq %d: parent i_seq %d (nuclear) and %d (X-ray) differ" % (
          i_h, nuclear[i_h][0], xray[i_h][0]))
      self.i_h.append(i_h)
      self.i_x.append(nuclear[i_h][0])
      self.d_xray.append(xray[i_h][1])
      self.d_nuclear.append(nuclear[i_h][1])

  @staticmethod
  def _table(model):
    """{i_h: (i_x, ideal length)} for the H/D that are in exactly one bond."""
    import mmtbx.model
    table = mmtbx.model.xh_connectivity_table(
      geometry=model.get_restraints_manager(),
      xray_structure=model.get_xray_structure()).table
    count = {}
    for entry in table:
      count[entry[1]] = count.get(entry[1], 0) + 1
    return dict((e[1], (e[0], e[3])) for e in table if count[e[1]] == 1)

  def _set_lengths(self, sites_cart, lengths):
    assert sites_cart.size() == self.n_atoms
    result = sites_cart.deep_copy()
    for i_h, i_x, length in zip(self.i_h, self.i_x, lengths):
      x = sites_cart[i_x]
      v = tuple(sites_cart[i_h][k] - x[k] for k in range(3))
      r = (v[0]**2 + v[1]**2 + v[2]**2)**0.5
      if r < 1.e-6:
        continue
      result[i_h] = tuple(x[k] + v[k] * length / r for k in range(3))
    return result

  def shrink(self, sites_cart):
    """
    Move every mapped H/D along its bond to the X-ray ideal length.

    Parameters
    ----------
    sites_cart : flex.vec3_double
      Coordinates of all atoms; not modified.

    Returns
    -------
    flex.vec3_double
      New coordinates of all atoms.
    """
    return self._set_lengths(sites_cart, self.d_xray)

  def expand(self, sites_cart):
    """
    Move every mapped H/D along its bond to the nuclear ideal length.

    Parameters
    ----------
    sites_cart : flex.vec3_double
      Coordinates of all atoms; not modified.

    Returns
    -------
    flex.vec3_double
      New coordinates of all atoms.
    """
    return self._set_lengths(sites_cart, self.d_nuclear)

  def selection(self):
    """
    Returns
    -------
    flex.bool
      True for the mapped H/D atoms.
    """
    result = flex.bool(self.n_atoms, False)
    result.set_selected(self.i_h, True)
    return result

  def gradients(self, sites_cart, gradients):
    """
    Gradient with respect to the nuclear sites of a target evaluated at the
    shrunk sites.

    For a mapped H with parent P, bond vector v = H - P, r = |v|, u = v / r,
    X-ray length d, s = d / r and Q = I - u u^T, the target gradients g_H and
    g_P at the shrunk sites give dT/dH = s Q g_H and
    dT/dP = g_P + (I - s Q) g_H. Bonds with r < 1e-6 are left as given.

    Parameters
    ----------
    sites_cart : flex.vec3_double
      Nuclear coordinates of all atoms; not modified.
    gradients : flex.vec3_double
      Target gradients at shrink(sites_cart), one per atom; not modified.

    Returns
    -------
    flex.vec3_double
      Gradients with respect to sites_cart.
    """
    assert sites_cart.size() == self.n_atoms
    assert gradients.size() == sites_cart.size()
    result = gradients.deep_copy()
    for i_h, i_x, d in zip(self.i_h, self.i_x, self.d_xray):
      x = sites_cart[i_x]
      v = tuple(sites_cart[i_h][k] - x[k] for k in range(3))
      r = (v[0]**2 + v[1]**2 + v[2]**2)**0.5
      if r < 1.e-6:
        continue
      u = tuple(c / r for c in v)
      g = gradients[i_h]
      ug = sum(u[k] * g[k] for k in range(3))
      s = d / r
      qg = tuple(g[k] - u[k] * ug for k in range(3))
      result[i_h] = tuple(s * c for c in qg)
      p = result[i_x]
      result[i_x] = tuple(p[k] + g[k] - s * qg[k] for k in range(3))
    return result
