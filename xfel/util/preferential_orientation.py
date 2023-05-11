from __future__ import division
import glob
import sys
from typing import List

from dxtbx.model import ExperimentList
from xfel.util.drift import params_from_phil, read_experiments

import numpy as np
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


class OrientationScraper:
  """Class responsible for scraping orientation matrices from expt files"""
  def __init__(self, parameters) -> None:
    self.parameters = parameters

  @staticmethod
  def assemble_orientation_stack(expts: ExperimentList) -> np.ndarray:
    """Extract N orientation matrices from N expts into a Nx3x3 numpy array"""
    ori_list = [np.array(expt.crystal.get_U()).reshape(3,3) for expt in expts]
    return np.stack(arrays=ori_list, axis=0)

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
    return self.assemble_orientation_stack(expts)


##################### PREFERENTIAL ORIENTATION CALCULATOR #####################

class PreferentialOrientationCalculator:
  """Calculator of the preferential orientation degree using IRLS algorithm"""
  EPSILON: float = 1e-6
  GUESS: np.eye(3)
  MAX_ITER: int = 100

  def irls_orientations(self, orientations: np.ndarray) -> np.ndarray:
    """Estimates the average orientation matrix from a dataset of
    3D orientation matrices using IRLS algorithm."""

    # Initial guess for A (arithmetic mean)
    A = self.GUESS

    # Initialize weights to identity matrices
    W = np.tile(np.eye(3), (len(orientations), 1, 1))

    for i in range(self.MAX_ITER):
      # Compute residuals
      R = np.array([logm(A.T @ d) for d in orientations])
      R = np.array([r - r.T for r in R])

      # Compute Jacobian
      J = np.zeros((3, 3))
      for r, w in zip(R, W):
        J += w @ r @ r.T @ w

      # Update A
      A_new = expm((np.linalg.inv(J) @ np.sum(R * W[:, np.newaxis], axis=0)).T)

      # Check convergence
      if np.allclose(A, A_new, rtol=self.EPSILON):
        break

      A = A_new

      # Update weights
      for j in range(len(orientations)):
        r = logm(A.T @ orientations[j])
        w = np.eye(3) - (2 / 3) * (r + r.T - np.trace(r) * np.eye(3))
        W[j] = np.linalg.inv(w)

    return A


########################### ORIENTATION VISUALIZING ###########################


################################ ENTRY POINTS #################################


def run(params_):
    os = OrientationScraper(parameters=params_)
    poc = PreferentialOrientationCalculator()
    ori = os.scrap()
    avg_ori = poc.irls_orientations(orientations=ori)
    print(avg_ori)


params = []
if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print(message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)