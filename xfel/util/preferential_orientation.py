from __future__ import division
import glob
import sys
from typing import List

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


def main2():
  mfc = MisesFisherCalculator()
  vectors = np.array([(1, 0, 0),
                      (0, 1, 0),
                      (0, 0, 1),
                      (-0.9, 0, 0),
                      (0, -1.1, 0),
                      (0, 0, -1)
                      ])
  print(f"{mfc.mu(vectors)}")
  print(f"{mfc.kappa(vectors)}")

if __name__ == '__main__':
  main2()