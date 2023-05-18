from __future__ import division
from dataclasses import dataclass
import glob
from typing import List
import sys

from dxtbx.model import ExperimentList
from xfel.util.drift import params_from_phil, read_experiments

from mpl_toolkits.mplot3d import Axes3D  # noqa: required to use 3D axes
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import scipy as sp


message = """
This utility tool aims to determinate, characterise, and quantify the degree
of preferential orientation in crystals. To this aim, it investigates
the distribution of vectors a, b, and c on a directions sphere in 3D.
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

This is a work in progress.
""".strip()


phil_scope_str = """
  scrap {
    input {
      glob = None
        .type = str
        .multiple = True
        .help = glob which matches all expt files to be investigated.
      exclude = None
        .type = str
        .multiple = True
        .help = glob which matches all expt files to be excluded from input.
    }
"""


############################ ORIENTATION SCRAPPING ############################


class DirectSpaceVectorScraper:
  """Class responsible for scraping vectors a, b, c from expt files"""
  def __init__(self, parameters) -> None:
    self.parameters = parameters

  @staticmethod
  def assemble_abc_stack(expts: ExperimentList) -> np.ndarray:
    """Extract N vectors a, b, c from N expts into a 3xNx3 numpy array.
    return_value[0] is a list of vectors "a"; [1] of "b", and [2] of "c"."""
    abc = [e.crystal.get_real_space_vectors().as_numpy_array() for e in expts]
    return np.stack(abc, axis=1)

  def locate_input_paths(self) -> List:
    """Return a list of expt paths in scrap.input.glob, but not in exclude"""
    input_paths, exclude_paths = [], []
    for ig in self.parameters.scrap.input.glob:
      input_paths.extend(glob.glob(ig))
    for ie in self.parameters.scrap.input.exclude:
      exclude_paths.extend(glob.glob(ie))
    return [it for it in input_paths if it not in exclude_paths]

  def scrap(self) -> np.ndarray:
    """Read and return a Nx3x3 orientation matrix based on scrap parameters"""
    expt_paths = self.locate_input_paths()
    expts = read_experiments(expt_paths)
    return self.assemble_abc_stack(expts)


##################### PREFERENTIAL ORIENTATION CALCULATOR #####################

class MisesFisherCalculator:
  EPSILON: float = 1e-6
  GUESS: np.ndarray = np.eye(3)
  MAX_ITER: int = 100

  @staticmethod
  def a3(kappa: float) -> float:
    """The normalization constant for the Bessel function of 1st kind in 3D."""
    return sp.special.iv(1.5, kappa) / sp.special.iv(0.5, kappa)

  def x_avg(self, vectors: np.ndarray) -> np.ndarray:
    """Sum of vectors divided by their count, bar{x}"""
    return np.sum(vectors, axis=0) / vectors.shape[0]

  def mu(self, vectors: np.ndarray) -> np.ndarray:
    """Mean direction of `vectors` normalized to 1, mu."""
    return self.x_avg(vectors) / self.r_avg(vectors)

  def r_avg(self, vectors: np.ndarray) -> float:
    """Length of the not-normalized mean direction of vectors, bar{R}"""
    return np.linalg.norm(self.x_avg(vectors))

  def kappa0(self, vectors: np.ndarray) -> float:
    """Simple approximation of kappa, following (Sra, 2011), hat{kappa}"""
    r = self.r_avg(vectors)
    return r * (3 - r ** 2) / (1 - r ** 2)

  def kappa(self, vectors):
    """Von Mises-Fisher concentration parameter of vectors on sphere, kappa"""
    k = self.kappa0(vectors)
    if k != 0:
      for i in range(10):
        a3 = self.a3(k)
        k = k - (a3 - self.r_avg(vectors)) / (1 - a3 ** 2 - (2 * a3 / k))
        print(k)
    return k


class SphericalDistribution:
  E1 = np.array([1, 0, 0])
  E2 = np.array([0, 1, 0])
  E3 = np.array([0, 0, 1])

  def __init__(self):
    self.vectors: np.ndarray = None
    self.mu: np.ndarray = np.array([1, 0, 0])

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
  def mu_basis_vectors(self):
    """Basis vector of cartesian system in which e1 = mu; e2 & e3 arbitrary"""
    e1 = self.mu / np.linalg.norm(self.mu)
    e0 = self.E1 if not self.are_parallel(e1, self.E1) else self.E2
    e2 = np.cross(e1, e0)
    e3 = np.cross(e1, e2)
    return e1, e2, e3

  def mu_sph2cart(self, vectors: np.ndarray):
    """Convert spherical coordinates r, polar, azim to cartesian in mu basis"""
    r, polar, azim = np.hsplit(vectors, 3)
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

  def __str__(self):
    return f'Watson Distribution around mu={self.mu} with kappa={self.mu}'

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

  def log_likelihood(self, mu: np.ndarray, kappa: float) -> float:
    """Log likelihood of given mu, kappa given current vectors. (10.3.30)"""
    t = self.scatter_matrix
    m = self.kummer_function(1/2, 3/2, kappa)
    return len(self.vectors) * (kappa * mu.T @ t @ mu - np.log(m))

  def nll_of_kappa(self, kappa: float, mu: np.ndarray) -> float:
    """Negative log likelihood of this Watson Distribution as a function
    of kappa, with `mu` and `vectors` fixed and given in `params`"""
    return -self.log_likelihood(mu=mu, kappa=kappa)

  def fit(self, vectors: np.ndarray) -> None:
    """Fit distribution to `vectors`, update `self.mu` and `self.kappa`"""
    self.vectors = vectors
    eig_val, eig_vec = np.linalg.eig(self.scatter_matrix)
    fitted = {'mu': np.array([1., 0., 0.]), 'kappa': 0., 'nll': np.inf}
    for eig_val, eig_vec in zip(eig_val, eig_vec.T):
        res = sp.optimize.minimize(self.nll_of_kappa, x0=0., args=eig_vec)
        if (nll := res['fun']) < fitted['nll']:
            fitted = {'mu': eig_vec, 'kappa': res['x'][0], 'nll': nll}
    self.mu = fitted['mu']
    self.kappa = fitted['kappa']

  def sample(self, n: int) -> np.ndarray:
    """Sample `n` vectors from self, based on doi 10.1080/03610919308813139"""
    if n < 0:
        return
    k = self.kappa
    rho = (4 * k) / (2 * k + 3 + ((2 * k + 3) ** 2 - 16 * k) ** 0.5)
    r = ((3 * rho) / (2 * k)) ** 3 * np.exp(-3 + 2 * k / rho)
    rng = np.random.default_rng()
    def cos2_of_polar_angle(_n: int) -> np.ndarray:
      u0 = rng.uniform(size=2*_n)
      u1 = rng.uniform(size=2*_n)
      s = u0 ** 2 / (1 - rho * (1 - u0 ** 2))
      v = (r * u1 ** 2) / (1 - rho * s) ** 3
      w = k * s
      good_s = s[v <= np.exp(2 * w)]
      return good_s[:_n] if (lgs := len(good_s)) > _n else \
          np.concatenate([good_s, cos2_of_polar_angle(_n-lgs)], axis=None)
    u2 = rng.uniform(size=n)
    theta = np.arccos(cos2_of_polar_angle(n) ** 0.5)
    phi = 4 * np.pi * u2
    theta[u2 < 0.5] = np.pi - theta[u2 < 0.5]
    phi[u2 >= 0.5] = 2 * np.pi * (2 * u2 - 1)
    return self.mu_sph2cart(np.vstack(np.ones_like(theta), theta, phi).T)


########################### ORIENTATION VISUALIZING ###########################

@dataclass
class Hedgehog:
  """Class for holding any `SphericalDistribution` with its metadata"""
  distribution: SphericalDistribution
  color: str
  name: str


class HedgehogArtist:
  """Class responsible for drawing distribution of vectors as "hedgehogs"."""
  def __init__(self, parameters) -> None:
    self.parameters = parameters
    self.hedgehogs = []
    self._init_figure()

  def __len__(self) -> int:
    return len(self.hedgehogs)

  def _init_figure(self) -> None:
    self.fig = plt.figure()
    self.axes = []

  def _generate_axes(self) -> None:
    gs_width = np.ceil(np.sqrt(len(self))).astype(int)
    gs_height = np.ceil(len(self) / gs_width).astype(int)
    gs = GridSpec(gs_height, gs_width, hspace=0, wspace=0)
    for h in range(gs_height):
      for w in range(gs_width):
        self.axes.append(self.fig.add_subplot(gs[h, w], projection='3d'))

  def _plot_hedgehog(self, axes: plt.Axes, hedgehog: Hedgehog) -> None:
    origin = [0., 0., 0.]
    name = hedgehog.name
    v = hedgehog.distribution.vectors
    axes.quiver(*origin, v[:, 0], v[:, 1], v[:, 2], colors=hedgehog.color)
    axes.set_xlim([-1, 1])
    axes.set_ylim([-1, 1])
    axes.set_zlim([-1, 1])
    axes.set_label(axes.get_label() + ' ' + name if axes.get_label() else name)

  def register_hedgehog(self, hedgehog: Hedgehog) -> None:
    self.hedgehogs.append(hedgehog)

  def plot(self):
    self._generate_axes()
    for axes, hedgehog in zip(self.axes, self.hedgehogs):
      self._plot_hedgehog(axes=axes, hedgehog=hedgehog)
    plt.show()


################################ ENTRY POINTS #################################


def run(params_):
    abc_stack = DirectSpaceVectorScraper(parameters=params_).scrap()
    hha = HedgehogArtist(parameters=params_)
    for vectors, color, name in zip(abc_stack, 'rgb', 'abc'):
      wd = WatsonDistribution.from_vectors(vectors)
      print(name + ': ' + str(wd))
      hh = Hedgehog(distribution=wd, color=color, name=name)
      hha.register_hedgehog(hh)
    hha.plot()


params = []
if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)
