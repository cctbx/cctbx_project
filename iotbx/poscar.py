"http://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html"
from __future__ import absolute_import, division, print_function

from libtbx import slots_getstate_setstate
from six.moves import range
from six.moves import zip

class reader(slots_getstate_setstate):

  __slots__ = """
    title
    scale
    lattice_vectors
    types
    type_counts
    sites
  """.split()

  def __init__(O, lines, source_info=None):
    assert len(lines) >= 7
    O.title = lines[0]
    scale_str = lines[1].split()
    assert len(scale_str) == 1
    O.scale = float(scale_str[0])
    O.lattice_vectors = []
    for i in [2,3,4]:
      vec_str = lines[i].split()
      assert len(vec_str) == 3
      vec = [float(_) for _ in vec_str]
      O.lattice_vectors.append(vec)
    i_type_counts = 5
    type_counts_str = lines[i_type_counts].split()
    assert len(type_counts_str) > 0
    _ = type_counts_str[0]
    try:
      int(_)
    except ValueError:
      O.types = type_counts_str
      i_type_counts = 6
      type_counts_str = lines[i_type_counts].split()
      assert len(type_counts_str) == len(O.types)
    else:
      O.types = None
    O.type_counts = [int(_) for _ in type_counts_str]
    key = lines[i_type_counts+1].strip()
    if (key != "Direct"):
      from libtbx.utils import Sorry
      if (source_info is not None):
        from libtbx.str_utils import show_string
        s = " of %s" % show_string(source_info)
      else:
        s = ""
      raise Sorry('POSCAR "%s" is not supported (line %d%s).'
        % (key, i_type_counts+1+1, s))
    n_sites = sum(O.type_counts)
    assert len(lines) >= i_type_counts+2+n_sites
    O.sites = []
    for i in range(i_type_counts+2, i_type_counts+2+n_sites):
      site_str = lines[i].split()
      assert len(site_str) >= 3
      site = [float(_) for _ in site_str[:3]]
      O.sites.append(site)

  def unit_cell(O):
    from scitbx import matrix
    vecs = [matrix.col(_) for _ in O.lattice_vectors]
    params = [_.length() for _ in vecs]
    for i in range(3):
      params.append(vecs[(i+1)%3].angle(vecs[(i+2)%3], deg=True))
    from cctbx import uctbx
    return uctbx.unit_cell(params)

  def make_up_types_if_necessary(O, default="Si"):
    if (O.types is None):
      O.types = [default]*len(O.type_counts)
    return O

  def scatterers(O, u_iso=0):
    assert O.types is not None
    from cctbx import xray
    from cctbx.array_family import flex
    result = flex.xray_scatterer()
    sites = iter(O.sites)
    for type,count in zip(O.types, O.type_counts):
      for _ in range(count):
        result.append(xray.scatterer(
          label="%s%d"%(type, len(result)+1),
          scattering_type=type,
          site=next(sites),
          u=u_iso))
    assert len(result) == len(O.sites)
    return result

  def xray_structure(O, u_iso=0):
    from cctbx import xray
    from cctbx import crystal
    return xray.structure(
      crystal_symmetry=crystal.symmetry(
        unit_cell=O.unit_cell(), space_group_symbol="P1"),
      scatterers=O.scatterers(u_iso=u_iso))
