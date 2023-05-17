from __future__ import division
import glob
import sys
from typing import Dict, List, Sequence

# from dxtbx.model import ExperimentList
# from xfel.util.drift import params_from_phil, read_experiments

import numpy as np
import scipy as sp
import scipy.spatial.transform
from scipy.linalg import expm, logm


message = """
This utility tool aims to determinate, characterise, and quantify the degree
of preferential orientation in crystals. To this aim, it first uses
the Iterative Reweighted Least Squares (IRLS) algorithm to estimate
the concentration matrix A of the Bingham distributio
from all experimental crystal matrices. Then it reports the type and degree
of anisotropy as expressed in the eigenvalues of said matrix.

This is work in progress.
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


# class OrientationScraper:
#   """Class responsible for scraping orientation matrices from expt files"""
#   def __init__(self, parameters) -> None:
#     self.parameters = parameters
#
#   @staticmethod
#   def assemble_orientation_stack(expts: ExperimentList) -> np.ndarray:
#     """Extract N orientation matrices from N expts into a Nx3x3 numpy array"""
#     ori_list = [np.array(expt.crystal.get_U()).reshape(3,3) for expt in expts]
#     return np.stack(arrays=ori_list, axis=0)
#
#   def locate_input_paths(self) -> List:
#     """Return a list of expt paths in scrap.input.glob, but not in exclude"""
#     input_paths, exclude_paths = [], []
#     for ig in self.parameters.scrap.input.glob:
#       input_paths.extend(glob.glob(ig))
#     for ie in self.parameters.scrap.input.exclude:
#       exclude_paths.extend(glob.glob(ie))
#     return [it for it in input_paths if it not in exclude_paths]
#
#   def scrap(self) -> np.ndarray:
#     """Read and return a Nx3x3 orientation matrix based on scrap parameters"""
#     expt_paths = self.locate_input_paths()
#     expts = read_experiments(expt_paths)
#     return self.assemble_orientation_stack(expts)


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


class WatsonDistribution:
  """The equations numbers given here refer to numbering in the book
  "Directional Statistics" by Kanti V. Mardia and Peter E. Jupp, Willey 2000"""
  def __init__(self):
    self.kappa: float = None
    self.mu: np.ndarray = None
    self.vectors: np.ndarray = None

  @property
  def r_avg(self) -> float:
    """Length of the not-normalized mean direction of vectors, bar{R}"""
    return np.linalg.norm(self.x_avg)

  @property
  def x_avg(self) -> np.ndarray:
    """Sum of vectors divided by their count, bar{x}"""
    return np.sum(self.vectors, axis=0) / self.vectors.shape[0]

  @staticmethod
  def normalized(vectors: np.ndarray, axis: int = -1,
                 order: int = 2) -> np.ndarray:
    """Return `vectors` normalized using `order` along `axis` """
    l2 = np.atleast_1d(np.linalg.norm(vectors, order, axis))
    l2[l2 == 0] = 1
    return vectors / np.expand_dims(l2, axis)

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

  def fit(self, vectors: np.ndarray):
    self.vectors = vectors
    eig_val, eig_vec = np.linalg.eig(self.scatter_matrix)
    fitted = {'mu': np.array([1., 0., 0.]), 'kappa': 0., 'nll': np.inf}
    for eig_val, eig_vec in zip(eig_val, eig_vec.T):
        res = sp.optimize.minimize(self.nll_of_kappa, x0=0., args=eig_vec)
        if (nll := res['fun']) < fitted['nll']:
            fitted = {'mu': eig_vec, 'kappa': res['x'][0], 'nll': nll}
    return fitted


########################### ORIENTATION VISUALIZING ###########################


################################ ENTRY POINTS #################################


# def run(params_):
#     os = OrientationScraper(parameters=params_)
#     poc = PreferentialOrientationCalculator()
#     ori = os.scrap()
#     avg_ori = poc.irls_orientations(orientations=ori)
#     print(avg_ori)


# params = []
# if __name__ == '__main__':
#   if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
#     print(message)
#     exit()
#   params = params_from_phil(sys.argv[1:])
#   run(params)


def normalized(a, axis=-1, order=2):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)


def main3():
  x = normalized(
      np.array([(0.01, 0.6, 0.6),
                (0, -0.5, 0.6),
                (0, 0.6, -0.6),
                (0, -0.6, -0.6),
                (0, 1, 0),
                (0, -1, 0),
                (0, 0, -1),
                (0, 0, 1),
                      ]))
  x = normalized(np.array([(0.01, 1, 0.01),
                (0.01, 1, -0.01),
                (-0.01, 1, 0.01),
                (-0.01, 1, -0.01),
                (0.01, -1, 0.01),
                (-0.01, -1, 0.01),
                (0.01, -1, -0.01)]))
  # Set initial guess for the parameters
  mu0 = np.mean(x, axis=0)
  kappa0 = 1.0


def main4():
    x = normalized(np.array([
        (0.01000124, 1, 0.01000234),
        (0.0100065, 1, -0.010005),
        (-0.010008, 1, 0.0100066),
        #(-0.0100098, 1, -0.01000423),
        #(0.01000234, -1, 0.010006),
        (-0.0100065, -1, 0.010004),
        (0.01000234, -1, -0.0100076),
        (-0.0100035, -1, -0.01000)]))
    # x = normalized(
    #     np.array([(0.01, 0.6, 0.6),
    #               (0, -0.5, 0.6),
    #               (0, 0.6, -0.6),
    #               (0, -0.6, -0.6),
    #               (0, 1, 0),
    #               (0, -1, 0),
    #               (0, 0, -1),
    #               (0, 0, 1),
    #               ]))
    wd = WatsonDistribution()
    res = wd.fit(x)
    print(res)
    assert 0
    print(x.T)
    print(x)
    scatter = np.matmul(x.T, x) / len(x)
    print(scatter)
    print(np.linalg.eig(scatter)[0])
    print(np.linalg.eig(scatter)[1].T)
    for val, vec in zip(np.linalg.eig(scatter)[0], np.linalg.eig(scatter)[1].T):
        print(f'val: {val:16f}, {vec=}')

"""
I can't quite get it to work. For references I used, see:
- https://arxiv.org/pdf/1104.4422.pdf, page 3
- http://palaeo.spb.ru/pmlibrary/pmbooks/mardia&jupp_2000.pdf, section 10.3.2
 
"""


if __name__ == '__main__':
  main4()