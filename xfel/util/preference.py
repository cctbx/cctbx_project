from __future__ import division

import abc
from collections import Counter, deque, UserDict
from dataclasses import dataclass
import glob
from numbers import Number
from typing import Any, Iterable, List, Sequence, Tuple, TypeVar, Union
import sys

from cctbx import sgtbx
from dxtbx.model import ExperimentList
from iotbx.phil import parse
from libtbx import Auto, table_utils
from libtbx.mpi4py import MPI
from xfel.util.drift import params_from_phil, read_experiments

import numpy as np
import scipy as sp


message = """
This utility tool aims to determinate, characterise, and quantify the degree
of preferential orientation in crystals. To this aim, it investigates the
the distribution of various crystallographic directions on a sphere in 3D.
The code assumes each set of vectors follows Wilson Distribution, and attempts
to model said distribution by fitting its parameter `mu` and `kappa`.

Wilson distribution describes a bimodal arrangement of points / unit vectors
on a sphere around a central axis called `mu`. The distribution is invariant
to any rotation around `mu` and inversion, and its exact type depends on the
concentration parameter `kappa`. For `kappa` > 0, the points are focused in
the polar region around +/- `mu`. In case of  `kappa` < 0, the points
concentrate mostly in a equatorial region far from `mu`. `kappa` close to 0
describes a distribution uniform on sphere: no preferential orientation.

This code has been prepared using the following books & papers as references:
- http://palaeo.spb.ru/pmlibrary/pmbooks/mardia&jupp_2000.pdf, section 10.3.2
- https://www.sciencedirect.com/science/article/pii/S0047259X12002084, sect. 2
- https://www.tandfonline.com/doi/abs/10.1080/03610919308813139
""".strip()


phil_scope_str = """
  input {
    glob = None
      .type = str
      .multiple = True
      .help = glob which matches all expt files to be investigated.
    exclude = None
      .type = str
      .multiple = True
      .help = glob which matches all expt files to be excluded from input.
    space_group = Auto
      .type = space_group
      .help = Reject all expts that are not described by this space group. \
              Investigate only directions that are symmetrically equivalent \
              according to the point symmetry of given group. \
              By default, look at the most common space group only.
    symmetrize = True
      .type = bool
      .help = Apply point group symmetry extracted from `space_group` to \
              extracted unit cell bases to avoid bias introduced by indexing
  }
  plot {
    style = none ascii *hammer hedgehog
      .type = choice
      .help = Which kind of plot should be produced: \
              Ascii writes a crude heatmap of distribution in x/y coords.\
              Hammer plots heatmap of distribution on Hammer projection.\
              Hedgehog draws all individual vectors (use for small data only).
  }
"""
phil_scope = parse(phil_scope_str)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ CONVENIENCE AND TYPING ~~~~~~~~~~~~~~~~~~~~~~~~~ #


COMM = MPI.COMM_WORLD
Int3 = Tuple[int, int, int]
Number3 = Sequence[Number]
sg_p1 = sgtbx.space_group('P 1')  # noqa
pg_i1 = sg_p1.build_derived_laue_group()
SgtbxPointGroup = Any
SgtbxSpaceGroup = Any
SgtbxSymmOp = Any
T = TypeVar('T')


def flatten(sequence: Sequence[Sequence[T]]) -> List[T]:
  """Flatten a sequence of sequences into a 1-dimensional list"""
  return [element for sub_sequence in sequence for element in sub_sequence]


def locate_paths(include_globs: Iterable[str],
                 exclude_globs: Iterable[str],
                 ) -> List[str]:
  """Return a list of expt paths in `include_globs`, outside `exclude_globs`"""
  include_paths = flatten([glob.glob(ig) for ig in include_globs])
  exclude_paths = flatten([glob.glob(eg) for eg in exclude_globs])
  return [path for path in include_paths if path not in exclude_paths]


def space_group_auto(expts: Iterable[ExperimentList],
                     comm: type(COMM) = None,
                     ) -> Tuple[SgtbxSpaceGroup, str]:
  """Return the most common space groups across comm world, and summary str"""
  counter = Counter([e.crystal.get_space_group().make_tidy() for e in expts])
  if comm is not None:
    counters = comm.allgather(counter)
    counter = sum(counters, Counter())
  message_ = ''
  for sg, sg_count in counter.items():
    message_ += f'Found {sg_count:6d} expts with space group {sg.info()}\n'
  most_common = counter.most_common(1)[0][0]
  message_ += f'Evaluating the most common space group {most_common.info()}'
  return most_common, message_


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SYMMETRY HANDLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def transform(vectors: Union[Number3, Iterable[Number3]],
              symm_op: SgtbxSymmOp,
              ) -> np.ndarray[Number3]:
  """Transform Number3 or Nx3 Iterable of Number3-s using symm_op """
  vectors = np.array(list(vectors))  # read generators, avoid tuples
  symm_op_m3 = np.array(symm_op.as_double_array()[:9]).reshape((3, 3))
  return vectors @ symm_op_m3.T


# ~~~~~~~~~~~~~~~~~~~~~~~~~~ ORIENTATION SCRAPPING ~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class DirectSpaceBases(np.ndarray):
  """
  This class is responsible for scraping and storing vectors a, b, c;
  Preferably, these come from expts, but the class can be initiated with
  a raw numpy array as well. It is 3-dimensional, with individual dimensions
  representing:

  - [n, :, :] - access n-th 3x3 array of a, b, c vectors;
  - [:, n, :] - access an Nx3 array of a (n=0), b (n=1), or c (n=2) vectors;
  - [:, :, n] - access an Nx3 array of x (n=0), y (n=1), or z (n=2) components.

  Consequently, the object is always (and must be init. using) a Nx3n3 array.
  """
  def __new__(cls, abcs: np.ndarray):
    if abcs.ndim != 3 or abcs.shape[1] != 3 or abcs.shape[2] != 3:
      msg = 'DirectSpaceVectors must be init with a Nx3x3 array of abc vectors'
      raise ValueError(msg)
    return super().__new__(cls, abcs.shape, dtype=float, buffer=abcs)

  @classmethod
  def from_expts(cls,
                 expts: ExperimentList,
                 space_group: SgtbxSpaceGroup = sg_p1,
                 ) -> 'DirectSpaceBases':
    """Extract N [a, b, c] arrays from N expts into a Nx3x3 ndarray, return"""
    abcs = []
    for expt in expts:
      expt_sg = expt.crystal.get_space_group()
      if expt_sg != space_group:
        continue
      abcs.append(expt.crystal.get_real_space_vectors().as_numpy_array())
    new = np.stack(abcs, axis=0) if abcs else np.empty(shape=(0, 3, 3))
    return cls(new)

  @property
  def a(self) -> np.ndarray:
    return np.array(self[:, 0, :])

  @property
  def b(self) -> np.ndarray:
    return np.array(self[:, 1, :])

  @property
  def c(self) -> np.ndarray:
    return np.array(self[:, 2, :])

  @property
  def x(self) -> np.ndarray:
    return np.array(self[:, :, 0])

  @property
  def y(self) -> np.ndarray:
    return np.array(self[:, :, 1])

  @property
  def z(self) -> np.ndarray:
    return np.array(self[:, :, 2])

  def transform(self, symm_op: SgtbxSymmOp) -> 'DirectSpaceBases':
    """Transform all vectors in self using sgtbx rt_mx-type object"""
    if self.size:
      x = transform(self.x, symm_op)
      y = transform(self.y, symm_op)
      z = transform(self.z, symm_op)
      new = np.stack([x, y, z], axis=2)
    else:
      new = np.empty(shape=(0, 3, 3))
    return self.__class__(new)

  def symmetrize(self, point_group: SgtbxPointGroup) -> 'DirectSpaceBases':
    """Transform all vectors in self using all symm. ops. in point group"""
    transformed = [self.transform(symm_op) for symm_op in point_group]
    return self.__class__(np.concatenate(transformed, axis=0))


# ~~~~~~~~~~~~~~~~~~~ PREFERENTIAL ORIENTATION CALCULATOR ~~~~~~~~~~~~~~~~~~~ #


def cart2sph(xyz: np.ndarray) -> np.ndarray:
  """Convert an array of xyz vectors into an array of r, polar, azim vectors"""
  r = np.linalg.norm(xyz, axis=1)
  polar = np.arccos(xyz[:, 2] / r)
  azim = np.arctan2(xyz[:, 1], xyz[:, 0])
  return np.vstack([r, polar, azim]).T


class SphericalDistribution:
  """
  General class for handling distribution of unit vectors in 3D. Operates in
  and provides methods to transform between different coordinate systems:
  - `cart` - Cartesian coordinates X, Y, Z in laboratory reference frame;
  - `sph` - Spherical ref. system with e1 = Z, polar e2 from Z to X, azim. XY;
  - `mu_sph` - Spherical reference system with e1 = mu and e2, e3 arbitrary;
  """
  E1 = np.array([1, 0, 0])
  E2 = np.array([0, 1, 0])
  E3 = np.array([0, 0, 1])

  def __init__(self):
    self.vectors: np.ndarray = np.empty((0, 3), dtype=float)
    self.mu: np.ndarray = np.array([1., 0., 0.], dtype=float)

  @staticmethod
  def normalized(vectors: np.ndarray, axis: int = -1) -> np.ndarray:
    """Return `vectors` normalized using standard l2 norm along `axis` """
    l2 = np.atleast_1d(np.linalg.norm(vectors, 2, axis))
    l2[l2 == 0] = 1
    return vectors / np.expand_dims(l2, axis)

  @staticmethod
  def are_parallel(v: np.ndarray, w: np.ndarray, eps: float = 1e-8) -> bool:
    return abs(np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w))) < 1 - eps

  @property
  def mu_basis_vectors(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Basis vector of cartesian system in which e1 = mu; e2 & e3 arbitrary"""
    e1 = self.mu / np.linalg.norm(self.mu)
    e0 = self.E1 if not self.are_parallel(e1, self.E1) else self.E2
    e2 = np.cross(e1, e0) / np.linalg.norm(np.cross(e1, e0))
    e3 = np.cross(e1, e2) / np.linalg.norm(np.cross(e1, e2))
    return e1, e2, e3

  @property
  def mu_ring(self) -> np.ndarray:
    """A closed 3D ring of points around vector self.mu used for plotting"""
    azim = np.linspace(0.0, 2*np.pi, num=100, endpoint=True)
    r = np.ones_like(azim)
    polar = np.full_like(azim, np.pi / 2)
    return self.mu_sph2cart(np.vstack([r, polar, azim]).T)

  def mu_sph2cart(self, r_polar_azim: np.ndarray) -> np.ndarray:
    """Convert spherical coordinates r, polar, azim to cartesian in mu basis"""
    r, polar, azim = np.hsplit(r_polar_azim, 3)
    e1, e2, e3 = self.mu_basis_vectors
    e1_component = e1 * np.cos(polar)
    e2_component = e2 * np.sin(polar) * np.cos(azim)
    e3_component = e3 * np.sin(polar) * np.sin(azim)
    return r * (e1_component + e2_component + e3_component)


class WatsonDistribution(SphericalDistribution):
  """Class for holding, fitting, and generating Watson distribution.
  Description names reflect those used in respective references:
  - https://arxiv.org/pdf/1104.4422.pdf, page 3
  - http://palaeo.spb.ru/pmlibrary/pmbooks/mardia&jupp_2000.pdf, section 10.3.2
  - https://www.tandfonline.com/doi/abs/10.1080/03610919308813139"""
  def __init__(self, mu: np.ndarray = None, kappa: float = None) -> None:
    super().__init__()
    self.kappa: float = kappa
    self.mu: np.ndarray = mu
    self.nll: float = np.Infinity

  def __str__(self) -> str:
    return f'Watson Distribution around mu={self.mu} with kappa={self.kappa}'

  @classmethod
  def from_vectors(cls, vectors: np.ndarray) -> 'WatsonDistribution':
    """Define the distribution by fitting it to a list of vectors"""
    new = WatsonDistribution()
    new.fit(vectors=cls.normalized(vectors))
    return new

  @property
  def r_avg(self) -> float:
    """Length of the not-normalized mean direction of vectors, bar{R}"""
    return np.linalg.norm(self.x_avg)

  @property
  def x_avg(self) -> np.ndarray:
    """Sum of vectors divided by their count, bar{x}"""
    return np.sum(self.vectors, axis=0) / self.vectors.shape[0]

  @staticmethod
  def kummer_function(a: float, b: float, kappa: float) -> float:
    """Confluent hypergeometric function 1F1, a.k.a. Kummer function"""
    return sp.special.hyp1f1(a, b, kappa)

  @property
  def scatter_matrix(self) -> np.ndarray:
    """Scatter matrix of `vectors` distribution (9.2.10)"""
    return np.matmul(self.vectors.T, self.vectors) / len(self.vectors)

  def log_likelihood(self, kappa: float, mu: np.ndarray) -> float:
    """Log likelihood of given mu, kappa given current vectors. (10.3.30)"""
    t = self.scatter_matrix
    m = self.kummer_function(1/2, 3/2, kappa)
    return len(self.vectors) * (kappa * mu.T @ t @ mu - np.log(m))

  def nll_of_kappa(self, kappa: float, mu: np.ndarray) -> float:
    """Negative log likelihood of this Watson Distribution as a function
    of kappa, with `mu` and `vectors` fixed and given in `params`"""
    return -self.log_likelihood(kappa=kappa, mu=mu)

  def fit(self, vectors: np.ndarray) -> None:
    """Fit distribution to `vectors`, update `self.mu` and `self.kappa`"""
    self.vectors = vectors
    eig_val, eig_vec = np.linalg.eig(self.scatter_matrix)
    fitted = {'mu': np.array([1., 0., 0.]), 'kappa': 0., 'nll': np.inf}
    for eig_val, eig_vec in zip(eig_val, eig_vec.T):
      result = sp.optimize.minimize(self.nll_of_kappa, x0=0., args=eig_vec)
      nll = result['fun']
      if nll < fitted['nll']:
        fitted = {'mu': eig_vec, 'kappa': result['x'][0], 'nll': nll}
    self.kappa = fitted['kappa']
    self.mu = fitted['mu']
    self.nll = fitted['nll']

  def sample(self, n: int, seed: int = 42) -> np.ndarray:
    """Sample `n` vectors from self, based on doi 10.1080/03610919308813139"""
    if n < 0:
      return
    k = self.kappa
    rho = (4 * k) / (2 * k + 3 + ((2 * k + 3) ** 2 - 16 * k) ** 0.5)
    r = ((3 * rho) / (2 * k)) ** 3 * np.exp(-3 + 2 * k / rho)
    rng = np.random.default_rng(seed=seed)

    def cos2_of_polar(_n: int) -> np.ndarray:
      u0 = rng.uniform(size=2*_n)
      u1 = rng.uniform(size=2*_n)
      s = u0 ** 2 / (1 - rho * (1 - u0 ** 2))
      v = (r * u1 ** 2) / (1 - rho * s) ** 3
      good_s = s[v <= np.exp(2 * k * s)]
      return good_s[:_n] if len(good_s) >= _n else \
          np.concatenate([good_s, cos2_of_polar(_n-len(good_s))], axis=None)
    u2 = rng.uniform(size=n)
    theta = np.arccos(cos2_of_polar(n) ** 0.5)
    phi = 4 * np.pi * u2
    theta[u2 < 0.5] = np.pi - theta[u2 < 0.5]
    phi[u2 >= 0.5] = 2 * np.pi * (2 * u2[u2 >= 0.5] - 1)
    self.vectors = self.mu_sph2cart(np.vstack([np.ones_like(theta), theta, phi]).T)


class ZoneAxisFamily(tuple):
  """Class for handling "crystal forms", i.e. sets of symm-equiv. zone axes."""

  def __str__(self) -> str:
    return f'{{{self[0]},{self[1]},{self[2]}}}'


class UniquePseudoNodeGenerator:
  """
  This class generates a list of unique pseudo-nodes; each pseudo-node
  represents a single pseudo-vector expressed using integer coordinates
  in cartesian space. They can be used to express all possible unique
  lattice directions or zone axes with indices up to `radius`.
  For example, pseudo-nodes [1, 1, 0], [-1, -1, 0], and [2, 2, 0] all express
  the same pseudo-vector [1, 1, 0], independent of symmetry
  """
  def __init__(self, laue_group: SgtbxPointGroup = pg_i1) -> None:
    self.point_group = laue_group
    self.nodes_to_yield = deque()
    self.nodes_considered = set()
    self.expand(around=(0, 0, 0), radius=3)

  def __iter__(self) -> 'UniquePseudoNodeGenerator':
    return self

  def __next__(self) -> Int3:
    if self.nodes_to_yield:
      return self.nodes_to_yield.popleft()
    raise StopIteration

  def add(self, nodes: Iterable[Int3]) -> None:
    """Add new pseudo-nodes, but only if they hadn't been yielded yet"""
    for node in nodes:
      node = tuple(node)
      symmetry_equivalents = {tuple(transform(node, symm_op))
                              for symm_op in self.point_group}
      if not symmetry_equivalents.intersection(self.nodes_considered):
        self.nodes_to_yield.append(node)
        self.nodes_considered.add(node)

  def expand(self, around: Int3, radius: int = 2) -> None:
    """Generate new direction pseudo-vectors in a `RADIUS` around `around`."""
    p_range = np.arange(around[0] - radius, around[0] + radius + 1)
    q_range = np.arange(around[1] - radius, around[1] + radius + 1)
    r_range = np.arange(around[2] - radius, around[2] + radius + 1)
    pqr_mesh = np.meshgrid(p_range, q_range, r_range)
    pqr = np.column_stack([mesh_comp.ravel() for mesh_comp in pqr_mesh])
    pqr = pqr[np.linalg.norm(pqr, axis=1) <= radius]
    p, q, r = pqr.T
    pqr = pqr[(p > 0) | ((p == 0) & (q > 0)) | ((p == 0) & (q == 0) & (r == 1))]
    pqr = pqr // np.gcd(np.gcd(pqr[:, 0], pqr[:, 1]), pqr[:, 2])[:, np.newaxis]
    self.add(np.unique(pqr, axis=0))


class PreferentialDistributionResults(UserDict[Int3, WatsonDistribution]):
  """UserDict holding information about fit results with convenient methods"""

  @property
  def best(self) -> Tuple[Int3, WatsonDistribution]:
    """Return a tuple w/ best fit (most offending) direction & distribution"""
    return list(self.sorted.items())[0]

  @property
  def sorted(self) -> 'PreferentialDistributionResults':
    return self.__class__(sorted(self.items(), key=lambda i: i[1].nll))

  def plot(self, kind: str = 'hedgehog'):
    """Plot all results in self as a hedgehog or hammer plot """
    artists = {'hedgehog': HedgehogArtist, 'hammer': HammerArtist}
    artist = artists[kind]()
    for direction, distribution in self.items():
      hh = Hedgehog(distribution=distribution, color='r', name=str(direction))
      artist.register_hedgehog(hh)
    artist.plot()

  @property
  def table(self) -> str:
    """Prepare a pretty string for logging"""
    table_data = [['Direction', 'kappa', 'mu', 'NLL']]
    for dir_, v in self.sorted.items():
      kappa = f'{v.kappa:+.3f}'
      mu = f'[{v.mu[0]:+.3f},{v.mu[1]:+.3f},{v.mu[2]:+.3f}]'
      nll = f'{v.nll:.2E}'
      table_data.append([str(dir_), kappa, mu, nll])
    return table_utils.format(table_data, has_header=1, delim='  ')


def find_preferential_distribution(
        dsv: DirectSpaceBases,
        space_group: SgtbxSpaceGroup
) -> PreferentialDistributionResults:
  """Look for a preferential orientation along any unique zone axis {hkl}"""
  laue_group = space_group.build_derived_laue_group()
  unique_pseudo_node_generator = UniquePseudoNodeGenerator(laue_group)
  results = PreferentialDistributionResults()
  for upn in unique_pseudo_node_generator:
    a_star = np.cross(dsv.b, dsv.c)  # not normalized by volume!
    b_star = np.cross(dsv.c, dsv.a)  # not normalized by volume!
    c_star = np.cross(dsv.a, dsv.b)  # not normalized by volume!
    vectors = a_star * upn[0] + b_star * upn[1] + c_star * upn[2]
    results[ZoneAxisFamily(upn)] = WatsonDistribution.from_vectors(vectors)
  return results


# ~~~~~~~~~~~~~~~~~~~~~~~~~ ORIENTATION VISUALIZING ~~~~~~~~~~~~~~~~~~~~~~~~~ #

@dataclass
class Hedgehog:
  """Class for holding any `SphericalDistribution` with its metadata"""
  distribution: SphericalDistribution
  color: str
  name: str


class BaseDistributionArtist(abc.ABC):
  """Base class for plotting hedgehogs of vector distributions"""
  PROJECTION: str = NotImplemented

  def __init__(self) -> None:
    self.hedgehogs = []
    from mpl_toolkits.mplot3d import Axes3D  # noqa: required to use 3D axes
    import matplotlib.pyplot as plt
    self.plt = plt
    self._init_figure()

  def _init_figure(self) -> None:
    self.fig = self.plt.figure(constrained_layout=True)
    self.axes = []

  def _generate_axes(self) -> None:
    from matplotlib.gridspec import GridSpec
    len_ = len(self.hedgehogs)
    axes_grid_width = np.ceil(np.sqrt(len_)).astype(int)
    axes_grid_height = np.ceil(len_ / axes_grid_width).astype(int)
    gs = GridSpec(axes_grid_height, axes_grid_width, figure=self.fig)
    for h in range(axes_grid_height):
      for w in range(axes_grid_width):
        ax = self.fig.add_subplot(gs[h, w], projection=self.PROJECTION)
        if len(self.axes) >= len_:
          ax.set_axis_off()
        self.axes.append(ax)

  def register_hedgehog(self, hedgehog: Hedgehog) -> None:
    self.hedgehogs.append(hedgehog)

  @abc.abstractmethod
  def plot(self) -> None:
    pass


class HedgehogArtist(BaseDistributionArtist):
  """Class responsible for drawing distribution of vectors as "hedgehogs"."""
  PROJECTION = '3d'

  def _plot_hedgehog(self, axes: 'plt.Axes', hedgehog: Hedgehog) -> None:
    origin = [0., 0., 0.]
    name = hedgehog.name
    v = hedgehog.distribution.vectors
    mu = hedgehog.distribution.mu
    mu_ring = hedgehog.distribution.mu_ring
    alpha = 1 / np.log2(len(v) + 1E-8)
    axes.quiver(*origin, v[:, 0], v[:, 1], v[:, 2], colors=hedgehog.color,
                alpha=alpha, arrow_length_ratio=0.0)
    axes.quiver(*-mu, *2*mu, colors='k', arrow_length_ratio=0.1)
    axes.plot(mu_ring[:, 0], mu_ring[:, 1], mu_ring[:, 2], color='k')
    axes.set_xlim([-1, 1])
    axes.set_ylim([-1, 1])
    axes.set_zlim([-1, 1])
    axes.set_label(axes.get_label() + ' ' + name if axes.get_label() else name)

  def plot(self):
    self._generate_axes()
    for ax, hedgehog in zip(self.axes, self.hedgehogs):
      self._plot_hedgehog(axes=ax, hedgehog=hedgehog)
    self.plt.show()


def calculate_geographic_heat(vectors: np.ndarray, n_bins: int = 10,
                              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
  """
  For an array of cartesian vectors, calculate and return:
    - edges of bins in azimuth coordinates,
    - edges of bins in polar coordinates,
    - heatmap on sphere in azimuth/polar coordinates,
  Assuming a geographic definition of polar angle (i.e. from pi/2 to -pi/2).
  """
  polar_edges = np.linspace(-np.pi / 2., np.pi / 2., n_bins + 1)
  azim_edges = np.linspace(-np.pi, np.pi, n_bins + 1)
  r_polar_azim = cart2sph(vectors)
  polar = np.pi / 2 - r_polar_azim[:, 1]
  azim = r_polar_azim[:, 2]
  polar_centers = (polar_edges[:-1] + polar_edges[1:]) / 2
  heat, azim_edges, polar_edges = np.histogram2d(
    x=azim, y=polar, bins=[azim_edges, polar_edges])
  heat = np.divide(heat, np.tile(np.cos(polar_centers), (n_bins, 1)))
  return azim_edges, polar_edges, heat.T  # heat transposed for x/y plotting


class HammerArtist(BaseDistributionArtist):
  """Class responsible for drawing distributions as hammer heatmaps"""
  PROJECTION = 'hammer'

  def _plot_hammer(self, ax: 'plt.Axes', hedgehog: Hedgehog) -> None:
    geo_heat = calculate_geographic_heat(vectors=hedgehog.distribution.vectors)
    ax.grid(False)
    ax.pcolor(*geo_heat, cmap=self.plt.get_cmap('viridis'))
    axes_params = {'ls': '', 'marker': 'o', 'mec': 'w'}  # lab x, y, and z-axes
    ax.plot(0., 0., c='r', **axes_params)
    ax.plot([-np.pi / 2, np.pi / 2], [0., 0.], c='g', **axes_params)
    ax.plot([0., 0.], [-np.pi / 2, np.pi / 2], c='b', **axes_params)
    ax.tick_params(labelbottom=False, labelleft=False)
    ax.set_title(hedgehog.name)

  def plot(self) -> None:
    self._generate_axes()
    for ax, hedgehog in zip(self.axes, self.hedgehogs):
      self._plot_hammer(ax=ax, hedgehog=hedgehog)
    self.plt.show()


def ascii_plot(vectors: np.ndarray, n_bins: int = 10) -> str:
  """A string with geographic heat plot on a simple xy cartesian coords"""
  px_width = 4
  _, _, heat = calculate_geographic_heat(vectors=vectors)
  minh, maxh = np.min(heat), np.max(heat)
  int_heat = np.rint(4.0 / (maxh - minh) * (heat.T - minh)).astype(int)
  colormap = ' ░▒▓█'
  plot_array = np.empty((px_width * n_bins + 2, n_bins + 2,), dtype=str)  # x/y
  plot_array[0, 0] = '┌'
  plot_array[-1, 0] = '┐'
  plot_array[0, -1] = '└'
  plot_array[-1, -1] = '┘'
  plot_array[1:-1, 0] = '─'
  plot_array[1:-1, -1] = '─'
  plot_array[0, 1:-1] = '│'
  plot_array[-1, 1:-1] = '│'
  for azim_i in range(n_bins):
    for polar_i in range(n_bins):
      azim_from = px_width * azim_i + 1
      azim_to = px_width * (azim_i + 1) + 1
      color = colormap[int_heat[azim_i, polar_i]]
      plot_array[azim_from:azim_to, polar_i + 1] = color
  mx = 1 + (px_width * n_bins) // 2
  my = 1 + n_bins // 2
  plot_array[0, my] = 'X'
  plot_array[mx, my] = 'X'
  plot_array[-1, my] = 'X'
  plot_array[mx // 2, my] = 'Y'
  plot_array[mx + mx // 2, my] = 'Y'
  plot_array[mx, 0] = 'Z'
  plot_array[mx, -1] = 'Z'
  return '\n'.join(''.join(c for c in line) for line in plot_array.T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ENTRY POINTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def run(params_) -> None:
  expt_paths = locate_paths(params_.input.glob, params_.input.exclude)
  expt_paths = expt_paths[COMM.rank::COMM.size]
  expts = read_experiments(*expt_paths)
  sgi = params_.input.space_group
  space_group = space_group_auto(expts, COMM)[0] if sgi is Auto else sgi.group()
  abc_stack = DirectSpaceBases.from_expts(expts, space_group)
  if params_.input.symmetrize:
    abc_stack = abc_stack.symmetrize(space_group.build_derived_point_group())
  abc_stacks = COMM.gather(abc_stack)
  if COMM.rank != 0:
    return
  abc_stack = DirectSpaceBases(np.concatenate(abc_stacks, axis=0))
  distributions = find_preferential_distribution(abc_stack, space_group)
  print(distributions.table)

  plot_style = params_.plot.style
  if plot_style != 'none':
    if plot_style == 'ascii':
      for direction, distribution in distributions.items():
        print(f'Ascii distribution heat plot for direction {direction}:')
        print(ascii_plot(distribution.vectors))
    else:
      distributions.plot(plot_style)


def exercise_watson_distribution() -> None:
  ha = HammerArtist()
  for kappa in [-10.0, -1.0, -0.1, -0.01, .000001, 0.01, 0.1, 1.0, 10.]:
    wd = WatsonDistribution(mu=np.array([0, 0, 1]), kappa=kappa)
    wd.sample(1_000_000)
    wd.fit(wd.vectors)
    print(wd)
    hh = Hedgehog(distribution=wd, color='r', name='kappa=' + str(kappa))
    ha.register_hedgehog(hh)
    print(ascii_plot(hh.distribution.vectors))
  ha.plot()


params = []
if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(phil_scope, sys.argv[1:])
  run(params)
