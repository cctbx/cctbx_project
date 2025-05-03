import sys
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Optional
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm

class SpotPair:

    def __init__(self, q1: float, q2: float, theta: float, preserve_order=False):
        # Ensure q1 <= q2
        if not preserve_order and q1 > q2:
            q1, q2 = q2, q1
        self.q1 = q1
        self.q2 = q2
        self.theta = theta

    def area(self) -> float:
        """Calculate twice the area of the triangle formed by the vectors."""
        v1 = np.array([self.q1, 0.0])
        v2 = self.q2 * np.array([np.cos(self.theta), np.sin(self.theta)])
        return abs(np.cross(v1, v2))

@dataclass
class GeneratedSpotPair():
    hkl1: np.ndarray
    hkl2: np.ndarray
    q1: np.float64
    q2: np.float64
    theta: np.float64

class IndexedSpotPair2d(SpotPair):
    def initialize_orientation(self, basis_vectors: np.ndarray) -> None:
        """Initialize orientation by computing angle between calculated and observed vectors."""
        # Get calculated vectors
        q1_calc = self.indices1[0]*basis_vectors[0] + self.indices1[1]*basis_vectors[1]
        q2_calc = self.indices2[0]*basis_vectors[0] + self.indices2[1]*basis_vectors[1]

        # Compute angle of q1_calc from x-axis
        phi_calc = np.arctan2(q1_calc[1], q1_calc[0])

        # Our observed q1 is along x-axis, so this is our basic rotation
        self.init_phi = phi_calc

        # Check if we need to flip orientation by comparing second vector
        q2_obs = self.q2 * np.array([np.cos(self.theta), np.sin(self.theta)])
        R = np.array([[np.cos(self.init_phi), -np.sin(self.init_phi)],
                      [np.sin(self.init_phi), np.cos(self.init_phi)]])
        q2_obs_rot = R @ q2_obs

        # If distance is large, try flipping
        if np.linalg.norm(q2_obs_rot - q2_calc) > np.linalg.norm(q2_obs_rot + q2_calc):
            self.init_phi += np.pi

        # Optional: small local optimization to refine this initial guess
        from scipy.optimize import minimize_scalar
        result = minimize_scalar(
            lambda dphi: self.compute_cost(basis_vectors, dphi),
            bounds=(-0.1, 0.1),  # small range around initial guess
            method='bounded'
        )
        self.init_phi += result.x

    def compute_cost(self, basis_vectors: np.ndarray, delta_phi: float) -> float:
        """Compute distance between rotated observed and calculated positions."""
        # Compute vectors from indices (fixed)
        q1_calc = self.indices1[0]*basis_vectors[0] + self.indices1[1]*basis_vectors[1]
        q2_calc = self.indices2[0]*basis_vectors[0] + self.indices2[1]*basis_vectors[1]

        # Create observed vectors in standard orientation (q1 along x)
        q1_obs = np.array([self.q1, 0.0])
        q2_obs = self.q2 * np.array([np.cos(self.theta), np.sin(self.theta)])

        # Total rotation to apply to observed vectors
        total_phi = self.init_phi + delta_phi
        R = np.array([[np.cos(total_phi), -np.sin(total_phi)],
                     [np.sin(total_phi), np.cos(total_phi)]])

        # Rotate observed vectors to match calculated
        q1_obs_rot = R @ q1_obs
        q2_obs_rot = R @ q2_obs

        return (np.linalg.norm(q1_obs_rot - q1_calc) + 
                np.linalg.norm(q2_obs_rot - q2_calc))


class VectorPairMatch:
  @classmethod
  def from_pairs(cls, gen_pair, obs_pair):
    return

class PairMatch2d(VectorPairMatch):
    def __init__(self, hkl1, q1, hkl2, q2, theta_rad):
        self.hkl1 = hkl1
        self.hkl2 = hkl2
        self.q1 = q1
        self.q2 = q2
        self.theta_rad = theta_rad
        self.is_outlier = False

class OneVectorMatch(VectorPairMatch):
    """
    q1 is the indexed vector and hkl1 is the corresponding indices in the sublattice.
    """
    def __init__(self, hkl1, q1, q2, theta_rad):
        self.hkl1 = hkl1
        self.q1 = q1
        self.q2 = q2
        self.theta_rad = theta_rad


# A couple helper functions

def angle_between(v1, v2):
    v1_u = v1/np.linalg.norm(v1)
    v2_u = v2/np.linalg.norm(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

def third_vector(q, v1, v2, theta1, theta2):
    """
    Compute a vector v3 with length q that forms angles theta1 and theta2
    with vectors v1 and v2 respectively.

    Parameters:
    q (float): Desired length of the output vector
    v1 (array-like): First reference vector
    v2 (array-like): Second reference vector
    theta1 (float): Desired angle with v1 (in radians)
    theta2 (float): Desired angle with v2 (in radians)

    Returns:
    numpy.ndarray: The computed vector v3
    """
    # Convert inputs to numpy arrays and normalize vectors
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    v1_norm = np.linalg.norm(v1)
    v2_norm = np.linalg.norm(v2)
    v1 = v1 / v1_norm
    v2 = v2 / v2_norm

    # The vector we're looking for can be written as a linear combination
    # of v1, v2, and their cross product: v3 = a*v1 + b*v2 + c*(v1 × v2)

    # First, get the cross product and normalize it
    v1xv2 = np.cross(v1, v2)
    if np.allclose(v1xv2, 0):
        raise ValueError("Input vectors are parallel, solution is not unique")
    v1xv2 = v1xv2 / np.linalg.norm(v1xv2)

    # The conditions are:
    # q*cos(theta1) = a + b*cos(gamma)
    # q*cos(theta2) = a*cos(gamma) + b
    # a^2 + b^2 + c^2 = q^2
    # where gamma is the angle between v1 and v2

    cos_gamma = np.dot(v1, v2)
    sin_gamma = np.sqrt(1 - cos_gamma**2)

    # Solve for a and b
    A = np.array([[1, cos_gamma],
                  [cos_gamma, 1]])
    b = q * np.array([np.cos(theta1), np.cos(theta2)])

    try:
        a, b = np.linalg.solve(A, b)

        # Now solve for c using the Pythagorean theorem
        c_sq = q**2 - (a**2 + b**2 + 2*a*b*cos_gamma)
        if c_sq < 0:
            raise ValueError("No solution exists for these angles")
        c = np.sqrt(c_sq)

        # There are two possible solutions (±c)
        # We'll return the positive c solution
        v3 = a*v1 + b*v2 + c*v1xv2

        return v3

    except np.linalg.LinAlgError:
        raise ValueError("No solution exists for these angles")

def find_third_vector(v1: np.ndarray, v2: np.ndarray,
                     q: float, theta1: float, theta2: float) -> Tuple[np.ndarray, np.ndarray]:
    """Find a vector v3 given its length and angles with v1 and v2.

    Args:
        v1, v2: Two known vectors in xy-plane (2D arrays)
        q: Length of vector to find
        theta1: Angle between v1 and v3
        theta2: Angle between v2 and v3

    Returns:
        Two possible 3D positions for v3 (above/below v1-v2 plane)
    """
    # Convert 2D vectors to 3D
    v1_3d = np.array([v1[0], v1[1], 0.0])
    v2_3d = np.array([v2[0], v2[1], 0.0])

    # Normalize
    v1_unit = v1_3d / np.linalg.norm(v1_3d)
    v2_unit = v2_3d / np.linalg.norm(v2_3d)

    # Get normal to v1-v2 plane (will be along z-axis)
    n = np.cross(v1_unit, v2_unit)
    n = n / np.linalg.norm(n)  # should be [0, 0, ±1]

    # Solve for components
    v1v2 = np.dot(v1_unit, v2_unit)
    A = np.array([[1, v1v2], [v1v2, 1]])
    b = q * np.array([np.cos(theta1), np.cos(theta2)])
    a, b = np.linalg.solve(A, b)

    # Find c from length condition
    c_sq = q*q - (a*a + b*b + 2*a*b*v1v2)
    if c_sq < 0:
        raise ValueError("No solution exists for these constraints")
    c = np.sqrt(c_sq)

    # Return both possible 3D positions
    v3_plus = a*v1_unit + b*v2_unit + c*n
    v3_minus = a*v1_unit + b*v2_unit - c*n

    return v3_plus, v3_minus

def find_best_third_vector(v3_candidates: Tuple[np.ndarray, np.ndarray],
                          sub_basis: np.ndarray) -> np.ndarray:
    """From two possible v3 positions, generate all lattice-equivalent points
    and choose the closest to origin.

    Args:
        v3_candidates: Two possible positions for v3
        sub_basis: Current 2x2 sublattice basis

    Returns:
        The best choice for the third basis vector
    """
    a, b = (np.hstack((x,0)) for x in sub_basis)
    best_v3 = None
    min_length = float('inf')

    for v3 in v3_candidates:
        # Generate all equivalent points v3 + ha + kb
        # Check a generous range of h,k values
        for h in range(-3, 4):
            for k in range(-3, 4):
                v3_equiv = v3 + h*a + k*b
                length = np.linalg.norm(v3_equiv)
                #print(h, k, round(length,5))

                if length < min_length:
                    min_length = length
                    best_v3 = v3_equiv

    return best_v3


class Basis:
    def __init__(self, vectors: np.ndarray, qmax: float, 
                 q_tolerance: float = 0.001, theta_tol_degrees: float = 1.0):
        self.vectors = vectors
        self.qmax = qmax
        self.q_tolerance = q_tolerance
        self.theta_tolerance = np.radians(theta_tol_degrees)
        self.points = None
        self.point_indices = None
        self.pairs = None
        self.indexed_pairs = []

    @classmethod
    def from_vectors(cls, vectors: np.ndarray, reduce=True, qmax:float = 0.5, q_tol: float = 0.001,
                     theta_tol_deg: float = 1.0):
        assert vectors.shape in ((2,2), (3,3))
        if vectors.shape == (2,2):
            temp_basis = Basis2d(vectors, qmax, q_tol, theta_tol_deg)

            if not reduce: # Short circuit
                return temp_basis

            temp_basis.generate_points_and_pairs_fast()
            points = temp_basis.points

            reduced_vectors = []
            distances = np.linalg.norm(points, axis=1)
            sort_idx = np.argsort(distances)
            for i in sort_idx:
                p = points[i]
                if len(reduced_vectors) == 0 and np.linalg.norm(p) > 1e-6:
                    reduced_vectors.append(p)
                elif len(reduced_vectors) == 1:
                    cross_prod = np.cross(reduced_vectors[0], p)
                    if abs(cross_prod) > 1e-6:
                        reduced_vectors.append(p)
                        break
            reduced_vectors = np.vstack(reduced_vectors)
            return Basis2d(reduced_vectors, qmax, q_tol, theta_tol_deg)

        elif vectors.shape == (3,3):
            temp_basis = Basis3d(vectors, qmax, q_tol, theta_tol_deg)

            if not reduce: # Short circuit
                return temp_basis
            temp_basis.generate_points_and_pairs_fast()
            points = temp_basis.points
            reduced_vectors = []
            distances = np.linalg.norm(points, axis=1)
            sort_idx = np.argsort(distances)

            for i in sort_idx:
                p = points[i]
                if len(reduced_vectors) == 0 and np.linalg.norm(p) > 1e-6:
                    reduced_vectors.append(p)
                elif len(reduced_vectors) == 1:
                    cross_prod = np.cross(reduced_vectors[0], p)
                    if np.linalg.norm(cross_prod) > 1e-6:
                        reduced_vectors.append(p)
                elif len(reduced_vectors) == 2:
                    v1, v2 = reduced_vectors
                    det = np.dot(np.cross(v1, v2), p)
                    if abs(det) > 1e-6:
                        reduced_vectors.append(p)
                        break
            reduced_vectors = np.vstack(reduced_vectors)
            return Basis3d(reduced_vectors, qmax, q_tol, theta_tol_deg)

    @classmethod
    def from_params(cls, *params, reduce=True, **init_kwargs):
        """
        params: q1 (A-1), q2, ga* (degrees) or q1, q2, q3, al*, be*, ga*.
        reduce: if True, return the setting from the shortest suitable (non-collinear
          etc) vectors. If False, return the setting as given.
        """
        assert len(params) in (3,6)

        if len(params) == 3:
            # 2D case: q1, q2, gamma*
            q1, q2, gamma_star_rad = params

            # Create basis vectors
            v1 = np.array([q1, 0.0])
            v2 = q2 * np.array([np.cos(gamma_star_rad), np.sin(gamma_star_rad)])
            vectors = np.vstack([v1, v2])

        else:
            # 3D case: q1, q2, q3, alpha*, beta*, gamma*
            q1, q2, q3, alpha_star, beta_star, gamma_star = params
            # Convert angles to radians
            alpha_star_rad = np.radians(alpha_star)
            beta_star_rad = np.radians(beta_star)
            gamma_star_rad = np.radians(gamma_star)

            # Create basis vectors using crystallographic conventions
            # First vector along x
            v1 = np.array([q1, 0.0, 0.0])

            # Second vector in xy plane
            v2 = q2 * np.array([np.cos(gamma_star_rad),
                               np.sin(gamma_star_rad),
                               0.0])

            # Third vector using all angles
            cx = np.cos(beta_star_rad)
            cy = (np.cos(alpha_star_rad) -
                 np.cos(beta_star_rad)*np.cos(gamma_star_rad))/np.sin(gamma_star_rad)
            cz = np.sqrt(1.0 - cx*cx - cy*cy)
            v3 = q3 * np.array([cx, cy, cz])

            vectors = np.vstack([v1, v2, v3])

        if reduce:
            # Create temporary basis to generate points
            temp_basis = cls.from_vectors(vectors, **init_kwargs)
            temp_basis.generate_points_and_pairs_fast()
            points = temp_basis.points

            if len(params) == 3:
                # Find two shortest non-collinear vectors
                basis_vectors = []
                distances = np.linalg.norm(points, axis=1)
                sort_idx = np.argsort(distances)

                for i in sort_idx:
                    p = points[i]
                    if len(basis_vectors) == 0 and np.linalg.norm(p) > 1e-6:
                        basis_vectors.append(p)
                    elif len(basis_vectors) == 1:
                        cross_prod = np.cross(basis_vectors[0], p)
                        if abs(cross_prod) > 1e-6:
                            basis_vectors.append(p)
                            break
                vectors = np.vstack(basis_vectors)

            else:
                # 3D reduction (might want to implement a more sophisticated method)
                # For now, just find three shortest non-coplanar vectors
                basis_vectors = []
                distances = np.linalg.norm(points, axis=1)
                sort_idx = np.argsort(distances)

                for i in sort_idx:
                    p = points[i]
                    if len(basis_vectors) == 0 and np.linalg.norm(p) > 1e-6:
                        basis_vectors.append(p)
                    elif len(basis_vectors) == 1:
                        cross_prod = np.cross(basis_vectors[0], p)
                        if np.linalg.norm(cross_prod) > 1e-6:
                            basis_vectors.append(p)
                    elif len(basis_vectors) == 2:
                        v1, v2 = basis_vectors
                        det = np.dot(np.cross(v1, v2), p)
                        if abs(det) > 1e-6:
                            basis_vectors.append(p)
                            break
                vectors = np.vstack(basis_vectors)

        return cls.from_vectors(vectors, **init_kwargs)

    def match(self, pair: SpotPair) -> Tuple[Optional[dict], str]:
        raise NotImplementedError

    def fom_1d(self, pairs, delta=.001):
        """Calculate figure of merit as percentage of observed q-values within threshold of lattice points.

    Args:
        pairs: List of SpotPair objects with observed q-values
        delta: Tolerance threshold in Å-1

    Returns:
        Percentage (0-100) of observed q-values that match lattice points
    """
        if not hasattr(self, 'points') or self.points is None:
            self.generate_points_and_pairs_fast()

        # Extract all observed q-values
        q_obs = []
        for p in pairs: 
            q_obs.extend([p.q1, p.q2])
        q_obs = np.array(q_obs)

        # Calculate all q-values from lattice points
        q_calc = np.linalg.norm(self.points, axis=1)

        # Reshape for broadcasting
        q_obs_col = q_obs.reshape(-1, 1)  # shape: (n_obs, 1)
        q_calc_row = q_calc.reshape(1, -1)  # shape: (1, n_calc)

        # Compute differences between all pairs
        diffs = np.abs(q_obs_col - q_calc_row)  # shape: (n_obs, n_calc)

        has_match = np.any(diffs < delta, axis=1)

        # Compute percentage
        pct_matched = 100.0 * np.sum(has_match) / len(q_obs)

        return pct_matched

    def generate_points_and_pairs_fast_old(self):
        """Generate unique lattice point pairs using matrix operations for efficiency."""
        # Generate points
        max_indices = np.ceil(self.qmax / np.linalg.norm(self.vectors, axis=1))
        ranges = [np.arange(-n, n+1) for n in max_indices.astype(int)]

        # For 2D case
        if len(ranges) == 2:
            # Only use positive h values (for non-zero h)
            h_range = np.concatenate([[0], np.arange(1, max_indices[0] + 1)])
            k_range = np.arange(-max_indices[1], max_indices[1] + 1)
        else:  # 3D case
            h_range = np.concatenate([[0], np.arange(1, max_indices[0] + 1)])
            k_range = np.arange(-max_indices[1], max_indices[1] + 1)
            l_range = np.arange(-max_indices[2], max_indices[2] + 1)

        # Mesh grid for indices
        if len(ranges) == 2:
            H, K = np.meshgrid(h_range, k_range, indexing='ij')
            indices = np.column_stack((H.flatten(), K.flatten()))
        else:  # 3D
            H, K, L = np.meshgrid(h_range, k_range, l_range, indexing='ij')
            indices = np.column_stack((H.flatten(), K.flatten(), L.flatten()))

        # Generate points
        points = indices @ self.vectors
        self.points = points

        # Filter by magnitude
        magnitudes = np.linalg.norm(points, axis=1)
        mask = (magnitudes <= self.qmax) & (magnitudes > 0)  # Exclude origin

        # Keep only points within qmax
        filtered_points = points[mask]
        filtered_indices = indices[mask]
        filtered_magnitudes = magnitudes[mask]

        self.points = filtered_points
        self.point_indices = filtered_indices

        # Sort by magnitude
        sort_idx = np.argsort(filtered_magnitudes)
        sorted_points = filtered_points[sort_idx]
        sorted_indices = filtered_indices[sort_idx]
        sorted_magnitudes = filtered_magnitudes[sort_idx]

        # Create normalized vectors for dot product
        normalized_points = sorted_points / sorted_magnitudes[:, np.newaxis]

        # Calculate all pairwise dot products (for normal points)
        n_points = len(sorted_points)
        dot_products = np.dot(normalized_points, normalized_points.T)

        # Create inverted indices (-h,-k,...)
        inverted_indices = -sorted_indices
        inverted_points = inverted_indices @ self.vectors

        # Normalize inverted points
        inverted_magnitudes = np.linalg.norm(inverted_points, axis=1)
        normalized_inverted = inverted_points / inverted_magnitudes[:, np.newaxis]

        # Calculate dot products between normal and inverted points
        mixed_dot_products = np.dot(normalized_points, normalized_inverted.T)

        # Get upper triangle indices (excluding diagonal)
        rows, cols = np.triu_indices(n_points, k=1)

        # Normal pairs
        q1_normal = sorted_magnitudes[rows]
        q2_normal = sorted_magnitudes[cols]
        hkl1_normal = sorted_indices[rows]
        hkl2_normal = sorted_indices[cols]
        theta_normal = np.arccos(np.clip(dot_products[rows, cols], -1, 1))

        # Inverted pairs
        q1_inverted = sorted_magnitudes[rows]
        q2_inverted = sorted_magnitudes[cols]  # Magnitude is the same for ±hkl
        hkl1_inverted = sorted_indices[rows]
        hkl2_inverted = inverted_indices[cols]
        theta_inverted = np.arccos(np.clip(mixed_dot_products[rows, cols], -1, 1))

        # Stack normal and inverted pairs
        self.q1_values = np.concatenate([q1_normal, q1_inverted])
        self.q2_values = np.concatenate([q2_normal, q2_inverted])
        self.theta_values = np.concatenate([theta_normal, theta_inverted])
        self.hkl1_values = np.vstack([hkl1_normal, hkl1_inverted])
        self.hkl2_values = np.vstack([hkl2_normal, hkl2_inverted])

        # Ensure q1 <= q2
        swap_mask = self.q1_values > self.q2_values
        if np.any(swap_mask):
            raise RuntimeError
            # Swap q values
            self.q1_values[swap_mask], self.q2_values[swap_mask] = self.q2_values[swap_mask], self.q1_values[swap_mask].copy()

            # Swap hkl values
            self.hkl1_values[swap_mask], self.hkl2_values[swap_mask] = self.hkl2_values[swap_mask], self.hkl1_values[swap_mask].copy()

    def generate_points_and_pairs_fast(self):
        """Generate unique lattice point pairs preserving q1 ordering."""
        # Generate points
        max_indices = np.ceil(self.qmax / np.linalg.norm(self.vectors, axis=1))
        ranges = [np.arange(-n, n+1) for n in max_indices.astype(int)]

        # Generate mesh grid based on dimension
        if len(ranges) == 2:  # 2D
            h_range = np.concatenate([[0], np.arange(1, max_indices[0] + 1)])
            k_range = np.arange(-max_indices[1], max_indices[1] + 1)
            H, K = np.meshgrid(h_range, k_range, indexing='ij')
            indices = np.column_stack((H.flatten(), K.flatten()))
        else:  # 3D
            h_range = np.concatenate([[0], np.arange(1, max_indices[0] + 1)])
            k_range = np.arange(-max_indices[1], max_indices[1] + 1)
            l_range = np.arange(-max_indices[2], max_indices[2] + 1)
            H, K, L = np.meshgrid(h_range, k_range, l_range, indexing='ij')
            indices = np.column_stack((H.flatten(), K.flatten(), L.flatten()))

        # Generate points
        points = indices @ self.vectors

        # Filter by magnitude
        magnitudes = np.linalg.norm(points, axis=1)
        mask = (magnitudes <= self.qmax) & (magnitudes > 0)  # Exclude origin

        # Keep only points within qmax
        filtered_points = points[mask]
        filtered_indices = indices[mask]
        filtered_magnitudes = magnitudes[mask]
        self.points = filtered_points
        self.point_indices = filtered_indices

        # Sort by magnitude (this will make q1_values naturally sorted)
        sort_idx = np.argsort(filtered_magnitudes)
        sorted_points = filtered_points[sort_idx]
        sorted_indices = filtered_indices[sort_idx]
        sorted_magnitudes = filtered_magnitudes[sort_idx]

        # Create inverted versions
        inverted_indices = -sorted_indices
        inverted_points = inverted_indices @ self.vectors

        # Create normalized vectors for dot products
        n_points = len(sorted_points)
        normalized_points = sorted_points / sorted_magnitudes[:, np.newaxis]
        inverted_normalized = inverted_points / np.linalg.norm(inverted_points, axis=1)[:, np.newaxis]

        # Create full matrices with NaN padding outside upper triangle

        # 1. First create mask for upper triangle (excluding diagonal)
        triu_mask = np.triu(np.ones((n_points, n_points)), k=1)
        nan_mask = ~triu_mask.astype(bool)

        # 2. Create q matrices (sorted by construction)
        q1_matrix = np.broadcast_to(sorted_magnitudes[:, np.newaxis], (n_points, n_points))
        q2_matrix = np.broadcast_to(sorted_magnitudes[np.newaxis, :], (n_points, n_points))

        # 3. Calculate theta matrices
        dot_products_normal = np.dot(normalized_points, normalized_points.T)
        theta_matrix_normal = np.arccos(np.clip(dot_products_normal, -1, 1))

        dot_products_inverted = np.dot(normalized_points, inverted_normalized.T)
        theta_matrix_inverted = np.arccos(np.clip(dot_products_inverted, -1, 1))

        # 4. Create hkl matrices
        # For normal pairs
        hkl1_shape = (n_points, n_points, sorted_indices.shape[1])
        hkl1_matrix_normal = np.broadcast_to(sorted_indices[:, np.newaxis, :], hkl1_shape)
        hkl2_matrix_normal = np.broadcast_to(sorted_indices[np.newaxis, :, :], hkl1_shape)

        # For inverted pairs
        hkl1_matrix_inverted = hkl1_matrix_normal.copy()
        hkl2_matrix_inverted = np.broadcast_to(inverted_indices[np.newaxis, :, :], hkl1_shape)

        # 5. Apply NaN mask to matrices
        q1_matrix_normal = q1_matrix.copy()
        q1_matrix_normal[nan_mask] = np.nan

        q2_matrix_normal = q2_matrix.copy()
        q2_matrix_normal[nan_mask] = np.nan

        theta_matrix_normal[nan_mask] = np.nan

        q1_matrix_inverted = q1_matrix.copy()
        q1_matrix_inverted[nan_mask] = np.nan

        q2_matrix_inverted = q2_matrix.copy()
        q2_matrix_inverted[nan_mask] = np.nan

        theta_matrix_inverted[nan_mask] = np.nan

        # 6. Stack normal and inverted matrices horizontally
        q1_stacked = np.hstack([q1_matrix_normal, q1_matrix_inverted])
        q2_stacked = np.hstack([q2_matrix_normal, q2_matrix_inverted])
        theta_stacked = np.hstack([theta_matrix_normal, theta_matrix_inverted])

        # 7. Flatten stacked matrices
        q1_flat = q1_stacked.flatten()
        q2_flat = q2_stacked.flatten()
        theta_flat = theta_stacked.flatten()

        # 8. Remove NaN values
        valid_mask = ~np.isnan(q1_flat)
        self.q1_values = q1_flat[valid_mask]
        self.q2_values = q2_flat[valid_mask]
        self.theta_values = theta_flat[valid_mask]

        # Handle hkl values (more complex due to extra dimension)
        hkl_dim = sorted_indices.shape[1]

        # Apply NaN mask to hkl matrices via boolean indexing
        # Create "is NaN" array matching hkl shape
        nan_mask_expanded = np.broadcast_to(nan_mask[:, :, np.newaxis],
                                            (n_points, n_points, hkl_dim))

        # Set invalid hkls to a dummy value (will be filtered out later)
        hkl1_matrix_normal = np.where(nan_mask_expanded, -999, hkl1_matrix_normal)
        hkl2_matrix_normal = np.where(nan_mask_expanded, -999, hkl2_matrix_normal)
        hkl1_matrix_inverted = np.where(nan_mask_expanded, -999, hkl1_matrix_inverted)
        hkl2_matrix_inverted = np.where(nan_mask_expanded, -999, hkl2_matrix_inverted)

        # Stack hkl matrices horizontally
        hkl1_stacked = np.hstack([hkl1_matrix_normal.reshape(n_points, -1),
                                 hkl1_matrix_inverted.reshape(n_points, -1)])
        hkl2_stacked = np.hstack([hkl2_matrix_normal.reshape(n_points, -1),
                                 hkl2_matrix_inverted.reshape(n_points, -1)])

        # Reshape to get original structure back
        hkl1_flat = hkl1_stacked.reshape(-1, hkl_dim)
        hkl2_flat = hkl2_stacked.reshape(-1, hkl_dim)

        # Apply same valid mask as for q values
        self.hkl1_values = hkl1_flat[valid_mask]
        self.hkl2_values = hkl2_flat[valid_mask]

        # Ensure q1 <= q2 (swap if needed)
        swap_mask = self.q1_values > self.q2_values
        if np.any(swap_mask):
            print('swap')
            self.q1_values[swap_mask], self.q2_values[swap_mask] = self.q2_values[swap_mask], self.q1_values[swap_mask].copy()
            self.hkl1_values[swap_mask], self.hkl2_values[swap_mask] = self.hkl2_values[swap_mask], self.hkl1_values[swap_mask].copy()




class Basis2d(Basis):

    def fom_2d(self, pairs):
        hits = 0
        for p in pairs:
            result = self.match(p)
            if result[1]=='indexed_2d':
                hits += 1
        return 100 * hits/len(pairs)
    def match(self, pair: SpotPair) -> Tuple[Optional[dict], str]:

        if not hasattr(self, 'q1_values'):
            self.generate_points_and_pairs_fast()

        # Compute differences with input pair
        dq1 = np.abs(self.q1_values - pair.q1)
        dq2 = np.abs(self.q2_values - pair.q2)
        dtheta = np.abs(self.theta_values - pair.theta)

        # Create mask for matches within tolerance
        matches = (dq1 < self.q_tolerance) & (dq2 < self.q_tolerance) & (dtheta < self.theta_tolerance)

        if np.any(matches):
            # Get the first match (or could find best one later)
            idx = np.where(matches)[0][0]
            hkl1 = self.hkl1_values[idx]
            hkl2 = self.hkl2_values[idx]

            result = PairMatch2d(hkl1, pair.q1, hkl2, pair.q2, pair.theta)
            return result, 'indexed_2d'

        qmags = np.linalg.norm(self.points, axis=1)
        dq1 = np.abs(pair.q1 - qmags)
        dq2 = np.abs(pair.q2 - qmags)
        i_best = np.argmin(np.minimum(dq1, dq2))
        one_vec_match = False
        if dq1[i_best] < self.q_tolerance:
            one_vec_match = True
            q1 = pair.q1
            q2 = pair.q2
        elif dq2[i_best] < self.q_tolerance:
            one_vec_match = True
            q1 = pair.q2
            q2 = pair.q1
        if one_vec_match:
            return OneVectorMatch(self.point_indices[i_best], q1, q2, pair.theta), 'one_vector'
        return None, 'unindexed'



    def astar(self):
        return(np.linalg.norm(self.vectors[0]))
    def bstar(self):
        return(np.linalg.norm(self.vectors[1]))
    def gammastar(self):
        result = np.degrees(np.arccos(
            np.dot(self.vectors[0], self.vectors[1]) / (a * b)))
        return result

    def __str__(self):
        """Format 2D sublattice parameters."""
        # Calculate a, b, gamma from basis vectors
        a = np.linalg.norm(self.vectors[0])
        b = np.linalg.norm(self.vectors[1])
        gamma = np.degrees(np.arccos(
            np.dot(self.vectors[0], self.vectors[1]) / (a * b)))

        return f"a={a:.5f}, b={b:.5f}, gamma={gamma:.2f}°"

    def area(self):
        return abs(np.cross(*self.vectors))

    def compute_pair_cost(self, basis_vectors: np.ndarray, pair: PairMatch2d,
                         q_weight: float = 1.0, theta_weight: float = 1.0, use_outliers=False) -> float:
        """Compute cost for a single indexed pair using Miller indices."""
        if pair.is_outlier and not use_outliers: return 0
        # Compute q-vectors from Miller indices
        q1_calc = pair.hkl1[0]*basis_vectors[0] + pair.hkl1[1]*basis_vectors[1]
        q2_calc = pair.hkl2[0]*basis_vectors[0] + pair.hkl2[1]*basis_vectors[1]

        # Compare magnitudes
        q1_calc_mag = np.linalg.norm(q1_calc)
        q2_calc_mag = np.linalg.norm(q2_calc)
        q_error = (abs(q1_calc_mag - pair.q1) + 
                  abs(q2_calc_mag - pair.q2))

        # Compare angle
        cos_theta_calc = np.dot(q1_calc, q2_calc) / (q1_calc_mag * q2_calc_mag)
        theta_calc = np.arccos(np.clip(cos_theta_calc, -1, 1))
        theta_error = abs(theta_calc - pair.theta_rad)

        return q_weight * q_error + theta_weight * theta_error

    def cost(self, params, pairs, q_weight=1000, theta_weight=1.0):
        """Total cost for given lattice parameters."""
        a_star, b_star, gamma_star = params

        # Basis vectors in standard orientation
        v1 = np.array([a_star, 0.0])
        v2 = b_star * np.array([np.cos(gamma_star), np.sin(gamma_star)])
        basis_vectors = np.vstack([v1, v2])

        # Sum costs from all pairs
        return sum(self.compute_pair_cost(basis_vectors, pair, q_weight, theta_weight)
                  for pair in pairs)

    def plot_cost_histogram(self, pairs, q_weight=1000, theta_weight=1.0):

        #copied from refine, sue me
        a_star = np.linalg.norm(self.vectors[0])
        b_star = np.linalg.norm(self.vectors[1])
        gamma_star = np.arccos(np.dot(self.vectors[0], self.vectors[1]) / 
                              (a_star * b_star))
        v1 = np.array([a_star, 0.0])
        v2 = b_star * np.array([np.cos(gamma_star), np.sin(gamma_star)])
        basis_vectors = np.vstack([v1, v2])

        costs = [self.compute_pair_cost(basis_vectors, p, q_weight, theta_weight) for p in pairs]
        plt.hist(costs, bins=100)
        plt.show()
    def refine(self, pairs, q_weight: float = 1000.0, theta_weight: float = 1.0) -> None:
        """Refine lattice parameters in place."""
        from scipy.optimize import minimize

        # Initial parameter vector
        a_star = np.linalg.norm(self.vectors[0])
        b_star = np.linalg.norm(self.vectors[1])
        gamma_star = np.arccos(np.dot(self.vectors[0], self.vectors[1]) / 
                              (a_star * b_star))
        initial_params = np.array([a_star, b_star, gamma_star])
        ip = copy.deepcopy(initial_params)

        # Run optimization
        result = minimize(
            lambda p: self.cost(p, pairs, q_weight, theta_weight),
            initial_params,
            method='BFGS',
            #options={'gtol': 1e-8}
        )

        if not result.success:
            print(result.message)

        # Extract refined parameters and update basis
        a_star, b_star, gamma_star = result.x
        if gamma_star < np.pi/2: gamma_star = np.pi - gamma_star
        self.vectors = np.array([[a_star, 0.0],
                               [b_star * np.cos(gamma_star), 
                                b_star * np.sin(gamma_star)]])

        # Regenerate points and pairs with new basis
        self.generate_points_and_pairs_fast()


class Basis3d(Basis):

    def volume(self) -> float:
        """Calculate the volume of the reciprocal space unit cell.

        Returns:
            Volume in Å⁻³
        """
        return np.abs(np.linalg.det(self.vectors))

    def __str__(self):
        # Reciprocal cell parameters
        a_star = np.linalg.norm(self.vectors[0])
        b_star = np.linalg.norm(self.vectors[1])
        c_star = np.linalg.norm(self.vectors[2])

        alpha_star = np.degrees(np.arccos(
            np.dot(self.vectors[1], self.vectors[2]) / (b_star * c_star)))
        beta_star = np.degrees(np.arccos(
            np.dot(self.vectors[0], self.vectors[2]) / (a_star * c_star)))
        gamma_star = np.degrees(np.arccos(
            np.dot(self.vectors[0], self.vectors[1]) / (a_star * b_star)))

        # Direct cell parameters
        direct = self.compute_direct_cell_params(self.vectors)

        return \
            f"Reciprocal cell:\n" + \
            f"  a*={a_star:.3f}, b*={b_star:.3f}, c*={c_star:.3f} Å⁻¹\n" + \
            f"  α*={alpha_star:.2f}°, β*={beta_star:.2f}°, γ*={gamma_star:.2f}°\n" + \
            f"Direct cell:\n" + \
            f"  a={direct['a']:.3f}, b={direct['b']:.3f}, c={direct['c']:.3f} Å\n" + \
            f"  α={direct['alpha']:.2f}°, β={direct['beta']:.2f}°, γ={direct['gamma']:.2f}°\n" + \
            ','.join( [str(round(direct[x], 3)) for x in ['a','b','c','alpha','beta','gamma']] )


    def compute_direct_cell_params(self, recip_basis: np.ndarray) -> dict:
        """Compute direct cell parameters from reciprocal space basis vectors.

        Args:
            recip_basis: 3x3 matrix of reciprocal space basis vectors

        Returns:
            dict with a,b,c (in Å) and alpha,beta,gamma (in degrees)
        """
        # Compute reciprocal metric tensor G* = B B^T
        G_star = recip_basis @ recip_basis.T

        # Invert to get direct metric tensor G = (G*)^-1
        G = np.linalg.inv(G_star)

        # Extract cell parameters
        a = np.sqrt(G[0,0])
        b = np.sqrt(G[1,1])
        c = np.sqrt(G[2,2])

        alpha = np.degrees(np.arccos(G[1,2] / (b*c)))
        beta = np.degrees(np.arccos(G[0,2] / (a*c)))
        gamma = np.degrees(np.arccos(G[0,1] / (a*b)))

        return {
            'a': a,
            'b': b,
            'c': c,
            'alpha': alpha,
            'beta': beta,
            'gamma': gamma
        }

    def match(self, pair: SpotPair) -> Tuple[Optional[dict], str]:
        """Find a matching pair in the lattice using binary search on sorted values."""

        if not hasattr(self, 'q1_values'):
            self.generate_points_and_pairs_fast()

        # Use binary search to find range of indices where q1 is within tolerance
        q1_min = pair.q1 - self.q_tolerance
        q1_max = pair.q1 + self.q_tolerance

        # Find indices where q1_values are in range [q1_min, q1_max]
        left_idx = np.searchsorted(self.q1_values, q1_min, side='left')
        right_idx = np.searchsorted(self.q1_values, q1_max, side='right')

        # If no values in range, return unindexed
        if left_idx >= right_idx:
            return None, 'unindexed'

        # Get the subset of potential matches
        subset_slice = slice(left_idx, right_idx)
        q2_subset = self.q2_values[subset_slice]
        theta_subset = self.theta_values[subset_slice]

        # Check q2 and theta matches in the smaller subset
        match_dq2 = (q2_subset > pair.q2 - self.q_tolerance) & \
                    (q2_subset < pair.q2 + self.q_tolerance)
        match_theta = (theta_subset > pair.theta - self.theta_tolerance) & \
                     (theta_subset < pair.theta + self.theta_tolerance)

        matches = match_dq2 & match_theta

        if np.any(matches):
            # Get the first match
            match_idx = np.where(matches)[0][0]
            full_idx = left_idx + match_idx  # Adjust back to original array index

            hkl1 = self.hkl1_values[full_idx]
            hkl2 = self.hkl2_values[full_idx]

            result = PairMatch2d(hkl1, pair.q1, hkl2, pair.q2, pair.theta)
            return result, 'indexed_3d'

        return None, 'unindexed'

    def doubled_cells(self) -> list:
        """Generate 6 cell-doubled variants of the current basis.

        Returns:
            List of 6 Basis3d objects with doubled cells:
            [0] Double a: a'=2a, b'=b, c'=c
            [1] Double b: a'=a, b'=2b, c'=c
            [2] Double c: a'=a, b'=b, c'=2c
            [3] Double ab plane: a'=a+b, b'=a-b, c'=c
            [4] Double ac plane: a'=a+c, b'=b, c'=a-c
            [5] Double bc plane: a'=a, b'=b+c, c'=b-c
        """
        # Get current reciprocal basis vectors
        a_star = self.vectors[0].copy()
        b_star = self.vectors[1].copy()
        c_star = self.vectors[2].copy()

        # Create transformation matrices in reciprocal space
        # When we double a direct space vector, the corresponding reciprocal vector is halved

        # 1. Double a: a*' = a*/2, b*' = b*, c*' = c*
        basis1 = self.vectors.copy()
        basis1[0] = a_star / 2

        # 2. Double b: a*' = a*, b*' = b*/2, c*' = c*
        basis2 = self.vectors.copy()
        basis2[1] = b_star / 2

        # 3. Double c: a*' = a*, b*' = b*, c*' = c*/2
        basis3 = self.vectors.copy()
        basis3[2] = c_star / 2

        # 4. Double ab plane: a'=a+b, b'=a-b in direct space
        # In reciprocal space: a*' = (a*+b*)/2, b*' = (a*-b*)/2, c*' = c*
        basis4 = self.vectors.copy()
        basis4[0] = (a_star + b_star) / 2
        basis4[1] = (a_star - b_star) / 2

        # 5. Double ac plane: a'=a+c, b'=b, c'=a-c in direct space
        # In reciprocal space: a*' = (a*+c*)/2, b*' = b*, c*' = (a*-c*)/2
        basis5 = self.vectors.copy()
        basis5[0] = (a_star + c_star) / 2
        basis5[2] = (a_star - c_star) / 2

        # 6. Double bc plane: a'=a, b'=b+c, c'=b-c in direct space
        # In reciprocal space: a*' = a*, b*' = (b*+c*)/2, c*' = (b*-c*)/2
        basis6 = self.vectors.copy()
        basis6[1] = (b_star + c_star) / 2
        basis6[2] = (b_star - c_star) / 2

        # 7. Body-centered (I-centered)
        # In reciprocal space: a*'=a*, b*'=b*, c*'=a*/2+b*/2+c*/2
        basis7 = self.vectors.copy()
        basis7[2] = (a_star + b_star + c_star) / 2

        # Create new Basis3d objects
        result = []
        for basis in [basis1, basis2, basis3, basis4, basis5, basis6, basis7]:
            new_basis = type(self).from_vectors(
                basis, 
                qmax=self.qmax,
                q_tol=self.q_tolerance,
                theta_tol_deg=np.degrees(self.theta_tolerance)
            )
            result.append(new_basis)

        return result

    def match_old(self, pair: SpotPair) -> Tuple[Optional[dict], str]:
        """Find a matching pair in the lattice using vectorized operations."""

        if not hasattr(self, 'q1_values'):
            self.generate_points_and_pairs_fast()

        # Compute differences with input pair
        dq1 = np.abs(self.q1_values - pair.q1)
        dq2 = np.abs(self.q2_values - pair.q2)
        dtheta = np.abs(self.theta_values - pair.theta)

        # Create mask for matches within tolerance
        match_dq1 = (self.q1_values > pair.q1-self.q_tolerance) & \
            (self.q1_values < pair.q1+self.q_tolerance)
        match_dq2 = (self.q2_values > pair.q2-self.q_tolerance) & \
            (self.q2_values < pair.q2+self.q_tolerance)
        match_theta = (self.theta_values > pair.theta-self.theta_tolerance) & \
            (self.theta_values < pair.theta+self.theta_tolerance)
        matches = match_dq1 & match_dq2 & match_theta
        #matches = (dq1 < self.q_tolerance) & (dq2 < self.q_tolerance) & (dtheta < self.theta_tolerance)

        if np.any(matches):
            # Get the first match (or could find best one later)
            idx = np.where(matches)[0][0]
            hkl1 = self.hkl1_values[idx]
            hkl2 = self.hkl2_values[idx]

            result = PairMatch2d(hkl1, pair.q1, hkl2, pair.q2, pair.theta)
            return result, 'indexed_3d'

        return None, 'unindexed'

    def match_multiple(self, pairs: List[SpotPair]) -> List[Tuple[Optional[dict], str]]:
        """Find matching pairs for a batch of input pairs.

        Args:
            pairs: List of SpotPair objects to match

        Returns:
            List of (result, status) tuples for each input pair
        """
        n_pairs = len(pairs)
        n_lattice_pairs = len(self.q1_values)

        # Extract arrays from input pairs
        q1_inputs = np.array([p.q1 for p in pairs])
        q2_inputs = np.array([p.q2 for p in pairs])
        theta_inputs = np.array([p.theta for p in pairs])

        # Reshape for broadcasting
        q1_inputs = q1_inputs.reshape(-1, 1)  # shape: (n_pairs, 1)
        q2_inputs = q2_inputs.reshape(-1, 1)
        theta_inputs = theta_inputs.reshape(-1, 1)

        # Compute differences with all lattice pairs
        # Each row is one input pair, each column is one lattice pair
        dq1 = np.abs(q1_inputs - self.q1_values)  # shape: (n_pairs, n_lattice_pairs)
        print(dq1.shape)
        dq2 = np.abs(q2_inputs - self.q2_values)
        dtheta = np.abs(theta_inputs - self.theta_values)

        # Find matches within tolerance
        matches = (dq1 < self.q_tolerance) & (dq2 < self.q_tolerance) & (dtheta < self.theta_tolerance)

        # Prepare results
        results = []

        # For each input pair
        for i in range(n_pairs):
            pair_matches = matches[i]

            if np.any(pair_matches):
                # Get the first match for this pair
                idx = np.where(pair_matches)[0][0]
                hkl1 = self.hkl1_values[idx]
                hkl2 = self.hkl2_values[idx]

                result = PairMatch2d(hkl1, pairs[i].q1, hkl2, pairs[i].q2, pairs[i].theta)
                results.append((result, 'indexed_3d'))
            else:
                results.append((None, 'unindexed'))

        return results

    def compute_pair_cost(self, basis_vectors: np.ndarray, pair: PairMatch2d,
                         q_weight: float = 1.0, theta_weight: float = 1.0, use_outliers=False) -> float:
        """Compute cost for a single indexed pair using Miller indices."""
        if not use_outliers and pair.is_outlier: return 0
        # Compute q-vectors from Miller indices
        q1_calc = pair.hkl1[0]*basis_vectors[0] + pair.hkl1[1]*basis_vectors[1] + pair.hkl1[2]*basis_vectors[2]
        q2_calc = pair.hkl2[0]*basis_vectors[0] + pair.hkl2[1]*basis_vectors[1] + pair.hkl2[2]*basis_vectors[2]

        # Compare magnitudes
        q1_calc_mag = np.linalg.norm(q1_calc)
        q2_calc_mag = np.linalg.norm(q2_calc)
        q_error = (abs(q1_calc_mag - pair.q1) +
                  abs(q2_calc_mag - pair.q2))

        # Compare angle
        cos_theta_calc = np.dot(q1_calc, q2_calc) / (q1_calc_mag * q2_calc_mag)
        theta_calc = np.arccos(np.clip(cos_theta_calc, -1, 1))
        theta_error = abs(theta_calc - pair.theta_rad)

        return q_weight * q_error + theta_weight * theta_error

    def vectors_from_params(self, params):
        """Create basis vectors from refinement parameters."""
        a_star, b_star, c_star, alpha_star, beta_star, gamma_star = params

        # First vector along x
        v1 = np.array([a_star, 0.0, 0.0])

        # Second vector in xy plane
        v2 = b_star * np.array([np.cos(gamma_star),
                               np.sin(gamma_star),
                               0.0])

        # Third vector using all angles
        cx = np.cos(beta_star)
        cy = (np.cos(alpha_star) -
             np.cos(beta_star)*np.cos(gamma_star))/np.sin(gamma_star)
        cz_sq = 1.0 - cx*cx - cy*cy
        if cz_sq < 0:
            # This happens with invalid angle combinations
            cz = 0.0
            # Add large penalty to cost function
            return None
        else:
            cz = np.sqrt(cz_sq)
        v3 = c_star * np.array([cx, cy, cz])

        return np.vstack([v1, v2, v3])

    def cost(self, params, pairs, q_weight=1000, theta_weight=1.0):
        """Total cost for given lattice parameters."""
        # Generate basis vectors from parameters
        basis_vectors = self.vectors_from_params(params)

        # If invalid parameters, return large penalty
        if basis_vectors is None:
            return 1.0e6

        # Sum costs from all pairs
        return sum(self.compute_pair_cost(basis_vectors, pair, q_weight, theta_weight)
                  for pair in pairs)

    def plot_cost_histogram(self, pairs, q_weight=1000, theta_weight=1.0):

        #copied from refine, sue me
        a_star = np.linalg.norm(self.vectors[0])
        b_star = np.linalg.norm(self.vectors[1])
        c_star = np.linalg.norm(self.vectors[2])

        alpha_star = np.arccos(np.dot(self.vectors[1], self.vectors[2]) /
                              (b_star * c_star))
        beta_star = np.arccos(np.dot(self.vectors[0], self.vectors[2]) /
                             (a_star * c_star))
        gamma_star = np.arccos(np.dot(self.vectors[0], self.vectors[1]) /
                              (a_star * b_star))

        current_params = np.array([a_star, b_star, c_star,
                                  alpha_star, beta_star, gamma_star])
        current_vectors = self.vectors_from_params(current_params)
        costs = [self.compute_pair_cost(current_params, p, q_weight, theta_weight) for p in pairs]
        plt.hist(costs, bins=100)
        plt.show()

    def refine(self, pairs, q_weight: float = 1000.0, theta_weight: float = 1.0) -> None:
        """Refine lattice parameters in place."""
        from scipy.optimize import minimize
        import copy

        # Initial parameter vector from current basis
        a_star = np.linalg.norm(self.vectors[0])
        b_star = np.linalg.norm(self.vectors[1])
        c_star = np.linalg.norm(self.vectors[2])

        alpha_star = np.arccos(np.dot(self.vectors[1], self.vectors[2]) /
                              (b_star * c_star))
        beta_star = np.arccos(np.dot(self.vectors[0], self.vectors[2]) /
                             (a_star * c_star))
        gamma_star = np.arccos(np.dot(self.vectors[0], self.vectors[1]) /
                              (a_star * b_star))

        initial_params = np.array([a_star, b_star, c_star,
                                  alpha_star, beta_star, gamma_star])

        # Define bounds to prevent unphysical values
        bounds = [
            (0.01, None),    # a_star > 0
            (0.01, None),    # b_star > 0
            (0.01, None),    # c_star > 0
            (0.1, np.pi-0.1),  # alpha_star in (0, pi)
            (0.1, np.pi-0.1),  # beta_star in (0, pi)
            (0.1, np.pi-0.1)   # gamma_star in (0, pi)
        ]

        # Run optimization
        result = minimize(
            lambda p: self.cost(p, pairs, q_weight, theta_weight),
            initial_params,
            method='L-BFGS-B',  # Use bounded optimization
            bounds=bounds,
            options={'ftol': 1e-10}
        )

        if not result.success:
            print(f"Refinement warning: {result.message}")

        # Extract refined parameters and update basis
        refined_params = result.x
        refined_vectors = self.vectors_from_params(refined_params)

        if refined_vectors is not None:
            self.vectors = refined_vectors

            # Regenerate points and pairs with new basis
            self.generate_points_and_pairs_fast()

            # Print refinement results
            refined_vols = self.volume()
            print(f"Refinement complete. Final cost: {result.fun:.6f}")
            print(f"Reciprocal volume: {refined_vols:.6f} Å⁻³")
            print(f"Direct cell volume: {1/refined_vols:.1f} Å³")


class LatticeReconstruction:
    def __init__(self, qmax: float, q_tolerance: float = 0.001, 
                 theta_tol_degrees: float = 1):
        self.qmax = qmax
        self.q_tolerance = q_tolerance
        self.theta_tol_degrees = theta_tol_degrees

        self.all_pairs = []
        self.sublattice_indexed = []
        self.sublattice_1vector = []
        self.unindexed = []
        self.current_basis = None
        self.sub_basis = None

    def calculate_sublattice_area(self) -> float:
        """Calculate area of current sublattice"""
        if self.sub_basis is None:
            return float('inf')
        return np.abs(np.linalg.det(self.sub_basis))

    def read_triplet(self, q1: float, q2: float, theta: float) -> None:
        """Read a new q1,q2,theta triplet and update reconstruction."""
        pair = SpotPair(q1, q2, np.radians(theta))  # assume input theta in degrees
        self.all_pairs.append(pair)
        self.update()

    def store_triplet(self, q1, q2, theta):
        pair = SpotPair(q1, q2, np.radians(theta))
        self.all_pairs.append(pair)

    def generate_2d_bases(self, scan_pts=21, scan_range=1, max_axis=25):
        scan_range=np.radians(scan_range)
        self.basis_candidates_2d = []
        for pair in self.all_pairs:
            deltas = np.linspace(-scan_range, scan_range, scan_pts)
            for d in deltas:
                b = Basis.from_params(pair.q1, pair.q2, pair.theta+d)
                if 1/b.astar() < max_axis and 1/b.bstar() < max_axis:
                    self.basis_candidates_2d.append(b)


    def update(self) -> None:
        """Main update method implementing the algorithm."""
        pair = self.all_pairs[-1] # Most recently read

        # Initialize first sublattice
        if len(self.all_pairs) == 1:
            self._initialize_first_sublattice(pair)
            self.print_status()
            return

        # Next, test for smaller sublattice
        elif pair.area() < self.sub_basis.area():
            trial_basis = Basis.from_params(pair.q1, pair.q2, pair.theta)
            print(f'Found smaller 2d basis: {trial_basis}')
            print(f'From pair {pair.q1}, {pair.q2}, {np.degrees(pair.theta)}')
            reproc = input('Reprocess? y/[n] ')=='y'
            if reproc:
                self.current_basis = None # Wipe out any 3d basis
                self._initialize_first_sublattice(pair)
                self.reprocess_all_pairs()
            else:
                pass

        # Try 3d indexing if possible
        elif self.current_basis is not None:
            result, status = self.current_basis.match(pair)
            if status=='indexed':
                self._store_result(result, status)
                self.print_status()
                return


        # Finally, try matching in the current sublattice
        else:

            # Process new pair
            result, status = self.sub_basis.match(pair)
            self._store_result(result, status)

        self.print_status()

    def _initialize_first_sublattice(self, pair: SpotPair) -> None:
        """Set up initial 2D sublattice from first pair using least oblique cell."""
        self.sub_basis = Basis.from_params(pair.q1, pair.q2, pair.theta)
        self.sl_from_pair = pair

    def set_sub_basis(self, basis):
        self.sub_basis = basis
        self.reprocess_all_pairs()

    def summarize_half_indexed_pairs(self) -> str:
        """Create summary table of half-indexed pairs."""
        if not self.sublattice_1vector:
          return "No half-indexed pairs found."

        lines = ["Half-indexed pairs for sublattice:",
                 "i | q[uni] | q[idx] | theta | err [Å⁻¹]"]
        entries = []

        for i, pair in enumerate(self.sublattice_1vector):
          # For each pair, determine which q is indexed
          q_idx, q_uni = pair.q1, pair.q2

          # Find nearest sublattice point to the indexed vector
          diffs = np.abs(np.linalg.norm(self.sub_basis.points, axis=1) - q_idx)
          err = np.min(diffs)

          #lines.append(f"{i:2d}   {q_uni:.4f}    {q_idx:.4f}    {np.degrees(pair.theta_rad):.1f}   {err:.4f}")
          entries.append((i, q_uni, q_idx, pair.theta_rad, err))

        for i, q_uni, q_idx, theta_rad, err in sorted(
            entries, 
            #key=lambda x:(round(x[1], 3), x[4])
            key=lambda x:(round(x[4], 4), round(x[1], 3))
            )[:50]:
            lines.append(f"{i:2d}\t{q_uni:.4f}\t{q_idx:.4f}\t{np.degrees(theta_rad):.1f}\t{err:.4f}")
        result = "\n".join(lines)

        return result

    def find_matching_pairs(self) -> List[Tuple[int, int, float]]:
        """Find pairs of half-indexed vectors with matching unindexed q values.

        Returns:
          List of (i1, i2, q) tuples where:
            i1, i2: indices into sublattice_1vector
            q: the matching q value
        """
        matches = []
        n = len(self.sublattice_1vector)

        for i in range(n):
          pair_i = self.sublattice_1vector[i]
          # Get unindexed q value
          q_i_uni = pair_i.q2 # The unindexed vector
          q_i_idx = pair_i.q1 # The indexed vector
          #q_i = pair_i.q2 if np.array_equal(pair_i.indexed_vector, pair_i.v1) else pair_i.q1

          for j in range(i+1, n):
            pair_j = self.sublattice_1vector[j]
            # Get unindexed q value
            q_j_uni = pair_j.q2
            q_j_idx = pair_j.q1

            if abs(q_i_uni - q_j_uni) < self.q_tolerance and abs(q_i_idx - q_j_idx) > self.q_tolerance:
              # Use average q value for the match
              q_match = (q_i_uni + q_j_uni) / 2
              matches.append((i, j, q_match))

        return matches

    def flag_outliers_2d(self, multiplier=2, q_weight=1000, theta_weight=1):
        costs = []
        for p in self.sublattice_indexed:
            costs.append(self.sub_basis.compute_pair_cost(
              self.sub_basis.vectors, p, q_weight, theta_weight))
        for p, c in zip(self.sublattice_indexed, costs):
            p.is_outlier = c > multiplier*np.median(costs)
    def plot_costs_2d(self, q_weight=1000, theta_weight=1):
        costs = []
        for p in self.sublattice_indexed:
            costs.append(self.sub_basis.compute_pair_cost(
              self.sub_basis.vectors, p, q_weight, theta_weight, use_outliers=True))
        costs_inlier = [c for c,p in zip(costs, self.sublattice_indexed) if not p.is_outlier]
        costs_outlier = [c for c,p in zip(costs, self.sublattice_indexed) if p.is_outlier]
        bins = np.linspace(0, max(costs), 50)
        plt.hist(costs_inlier, bins=bins)
        plt.hist(costs_outlier, bins=bins, color='red')
        plt.show()
        
        
    def reprocess_all_pairs(self) -> None:
        """Clear all categorizations and reprocess all pairs except the last."""
        stored_pairs = self.all_pairs # [:-1]
        self.sublattice_indexed = []
        self.sublattice_1vector = []
        self.unindexed = []
        for old_pair in stored_pairs:
            result, status = self.sub_basis.match(old_pair)
            self._store_result(result, status)

    def _store_result(self, result: SpotPair, status: str) -> None:
        """Store a processed pair in appropriate category."""
        if status == 'indexed_2d':
            self.sublattice_indexed.append(result)
        elif status == 'indexed_3d':
            self.indexed.append(result)

        elif status == 'one_vector':
            self.sublattice_1vector.append(result)
        else:  # unindexed
            self.unindexed.append(result)

    def generate_3d_basis_indices(self, pair1_idx, pair2_idx, delta_theta_1=None, delta_theta_2=None, inv=False):
        pair1 = self.sublattice_1vector[pair1_idx]
        pair2 = self.sublattice_1vector[pair2_idx]
        return generate_3d_basis_pairs(pair1, pair2, delta_theta_1, delta_theta_2, inv)

    def generate_3d_basis_pairs(self, pair1, pair2, delta_theta_1=None, delta_theta_2=None, inv=False):
        """Generate a 3D basis from two matching one-vector pairs."""

        # Get indexed vectors from sublattice (2D)
        v1 = np.hstack((pair1.hkl1 @ self.sub_basis.vectors, 0))
        v2 = np.hstack((pair2.hkl1 @ self.sub_basis.vectors, 0))

        # Get common q value and both angles
        q = pair1.q2
        theta1 = pair1.theta_rad
        theta2 = pair2.theta_rad

        parity = -1 if inv else 1
        if delta_theta_1 is not None:
            theta1 += np.radians(delta_theta_1)
        if delta_theta_2 is not None:
            theta2 += np.radians(delta_theta_2)
        # Find possible positions for third vector (returns 3D vectors)
        try:
            v3_candidates = [third_vector(q, v1, parity*v2, theta1, theta2)]
        except ValueError:
            return None

        # Choose best third vector
        v3 = find_best_third_vector(v3_candidates, self.sub_basis.vectors)

        # Create full 3D basis by extending 2D vectors
        v1v2 = np.hstack((self.sub_basis.vectors, np.zeros((2,1))))
        basis_vectors = np.vstack((v1v2, v3))

        result = Basis3d.from_vectors(
            vectors=basis_vectors,
            qmax=self.sub_basis.qmax,
            q_tol=self.sub_basis.q_tolerance,
            theta_tol_deg=np.degrees(self.sub_basis.theta_tolerance)
        )
        return result

    def print_status(self, verbose=False) -> None:
        """Print current status of reconstruction."""
        print(f"\nCurrent minimum sublattice: {self.sub_basis}")
        if hasattr(self, 'full_basis'):
            print(f"Current full lattice:\n{self.full_basis}")
        else:
            print("No full lattice determined yet.")

        print("\n---")
        if not verbose: return
        print(self.summarize_half_indexed_pairs())
        return

        # If we have matching pairs, show possible 3D cells
        matches = self.find_matching_pairs()
        if matches:
            print("\nLattice completion possibilities:")
            print("entry | i1 | i2 | cell | vol (Å³) | % idx")
            for i, (i1, i2, q) in enumerate(matches):
                try:
                    basis3d = self.generate_3d_basis_indices(i1, i2)
                    # Calculate metrics
                    volume = np.abs(np.linalg.det(basis3d.vectors))
                    n_indexed = len([p for p in self.all_pairs
                                   if basis3d.match(p)[1] == 'indexed_3d'])
                    pct_indexed = 100 * n_indexed / len(self.all_pairs)

                    print(f"{i}) {i1:2d}  {i2:2d}   {basis3d}  {volume:.1f}  {pct_indexed:.1f}")
                except ValueError as e:
                    # Skip if no solution exists
                    #print(f"{i}) {i1:2d}  {i2:2d}   No valid solution")
                    pass

# Grid search stuff, temporary

# Define grid search parameters
def grid_search_3d_basis(recon, pair1_idx, pair2_idx,
                         delta_range=(-2.0, 2.0), steps=21, inv=False):
    # Create grid of delta theta values
    delta_values = np.linspace(delta_range[0], delta_range[1], steps)
    grid_shape = (steps, steps)
    fom = np.zeros(grid_shape)

    # Total number of pairs
    total_pairs = len(recon.all_pairs)

    # Grid search
    for i, delta1 in tqdm(enumerate(delta_values)):
        for j, delta2 in tqdm(enumerate(delta_values), leave=False):
            try:
                # Generate 3D basis with current deltas
                basis3d = recon.generate_3d_basis_pairs(
                    pair1_idx, pair2_idx,
                    delta_theta_1=delta1,
                    delta_theta_2=delta2,
                    inv=inv
                )

                # Count indexed pairs
                hits = 0
                vol = basis3d.volume()
                if vol > .0002:
                    for p in recon.all_pairs:
                        if basis3d.match(p)[1] != 'unindexed':
                            hits += 1

                # Compute percentage
                fom[i, j] = 100*hits / total_pairs * vol

            except Exception as e:
                # Failed to generate basis (e.g., no valid solution)
                fom[i, j] = 0

    return delta_values, fom

# Run grid search for both inv=False and inv=True
def run_both_grid_searches(recon, pair1_idx, pair2_idx,
                          delta_range=(-2.0, 2.0), steps=21):
    print(f"Running grid search for pair indices {pair1_idx} and {pair2_idx}...")
    print(f"Delta range: {delta_range}, steps: {steps}")

    delta_values, results_normal = grid_search_3d_basis(
        recon, pair1_idx, pair2_idx, delta_range, steps, inv=False
    )

    delta_values, results_inv = grid_search_3d_basis(
        recon, pair1_idx, pair2_idx, delta_range, steps, inv=True
    )

    return delta_values, results_normal, results_inv

# Plot the results as heatmaps
def plot_grid_search_results(delta_values, results_normal, results_inv):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Custom colormap (white to blue)
    colors = [(1, 1, 1), (0, 0, 1)]
    cmap = LinearSegmentedColormap.from_list('WhiteToBlue', colors)

    # Determine shared vmax for consistent coloring
    vmax = max(np.max(results_normal), np.max(results_inv))

    # Plot normal results
    im1 = ax1.imshow(results_normal, extent=[delta_values[0], delta_values[-1],
                                           delta_values[0], delta_values[-1]],
                    origin='lower', cmap=cmap, vmin=0, vmax=vmax)
    ax1.set_title('Normal (inv=False)')
    ax1.set_xlabel('Delta theta 2 (degrees)')
    ax1.set_ylabel('Delta theta 1 (degrees)')

    # Plot inverted results
    im2 = ax2.imshow(results_inv, extent=[delta_values[0], delta_values[-1],
                                        delta_values[0], delta_values[-1]],
                   origin='lower', cmap=cmap, vmin=0, vmax=vmax)
    ax2.set_title('Inverted (inv=True)')
    ax2.set_xlabel('Delta theta 2 (degrees)')

    # Add colorbar
    cbar = fig.colorbar(im1, ax=[ax1, ax2], orientation='vertical', shrink=0.8)
    cbar.set_label('Indexing percentage (%)')

    # Add max value annotations
    max_normal = np.max(results_normal)
    max_normal_idx = np.unravel_index(np.argmax(results_normal), results_normal.shape)
    delta1_normal = delta_values[max_normal_idx[0]]
    delta2_normal = delta_values[max_normal_idx[1]]

    max_inv = np.max(results_inv)
    max_inv_idx = np.unravel_index(np.argmax(results_inv), results_inv.shape)
    delta1_inv = delta_values[max_inv_idx[0]]
    delta2_inv = delta_values[max_inv_idx[1]]

    ax1.plot(delta2_normal, delta1_normal, 'r+', markersize=10)
    ax1.text(delta2_normal, delta1_normal, f' {max_normal:.1f}%', color='red')

    ax2.plot(delta2_inv, delta1_inv, 'r+', markersize=10)
    ax2.text(delta2_inv, delta1_inv, f' {max_inv:.1f}%', color='red')

    plt.tight_layout()
    return fig

# Function to run everything and return the best parameters
def find_best_3d_basis(recon, pair1_idx, pair2_idx,
                       delta_range=(-2.0, 2.0), steps=21):
    # Run grid searches
    delta_values, results_normal, results_inv = run_both_grid_searches(
        recon, pair1_idx, pair2_idx, delta_range, steps
    )

    # Plot results
    fig = plot_grid_search_results(delta_values, results_normal, results_inv)

    # Find best parameters
    max_normal = np.max(results_normal)
    max_normal_idx = np.unravel_index(np.argmax(results_normal), results_normal.shape)
    delta1_normal = delta_values[max_normal_idx[0]]
    delta2_normal = delta_values[max_normal_idx[1]]

    max_inv = np.max(results_inv)
    max_inv_idx = np.unravel_index(np.argmax(results_inv), results_inv.shape)
    delta1_inv = delta_values[max_inv_idx[0]]
    delta2_inv = delta_values[max_inv_idx[1]]

    # Choose best overall parameters
    if max_normal >= max_inv:
        best_params = {
            'delta_theta_1': delta1_normal,
            'delta_theta_2': delta2_normal,
            'inv': False,
            'fom': max_normal
        }
    else:
        best_params = {
            'delta_theta_1': delta1_inv,
            'delta_theta_2': delta2_inv,
            'inv': True,
            'fom': max_inv
        }

    # Generate the best basis
    best_basis = recon.generate_3d_basis_pairs(
        pair1_idx, pair2_idx,
        delta_theta_1=best_params['delta_theta_1'],
        delta_theta_2=best_params['delta_theta_2'],
        inv=best_params['inv']
    )

    idx_pct = best_params['fom']/best_basis.volume()
    print("\nBest parameters:")
    print(f"delta_theta_1 = {best_params['delta_theta_1']:.2f} degrees")
    print(f"delta_theta_2 = {best_params['delta_theta_2']:.2f} degrees")
    print(f"inv = {best_params['inv']}")
    print(f"Indexing percentage = {idx_pct:.1f}%")


    return best_basis, best_params, fig




def run():
    recon = LatticeReconstruction(qmax=.5)
    for l in open(sys.argv[1]):
        recon.store_triplet(*[float(x) for x in l.split()])
    print('generate')
    recon.generate_2d_bases()
    print('fom1')
    basis_fom = [
        (b, b.fom_1d(recon.all_pairs) )
        for b in recon.basis_candidates_2d
    ]
    basis_fom.sort(key=lambda x:x[1], reverse=True)
    top_1d = basis_fom[:50]
    print('fom2')
    basis_fom1_fom2 = [
        (b, f, b.fom_2d(recon.all_pairs))
        for b, f in top_1d
    ]
    basis_fom1_fom2.sort(key=lambda x:x[2], reverse=True)
    for i, b in enumerate(basis_fom1_fom2[:20]):
        print(i, b[0], round(b[1], 1), round(b[2], 1))
    i_sb = int(input('use what sub-basis?'))
    best_sb = basis_fom1_fom2[i_sb]
    print('using 2d basis: ', best_sb[0], round(best_sb[1],2), round(best_sb[2],2))
    recon.set_sub_basis(best_sb[0])
    import IPython;IPython.embed()
    for _ in range(2):
        recon.sub_basis.refine(recon.sublattice_indexed)
        recon.reprocess_all_pairs()

    recon.print_status(verbose=True)

    from cluster2 import ManualClusterer
    fn2 = sys.argv[2]
    data = np.load(fn2)['triplets'][:,1:4]
    data[:,0] = 1/data[:,0]
    data[:,1] = 1/data[:,1]
    data2 = np.vstack((data[:,1], data[:,0], data[:,2])).transpose()
    data = np.vstack((data, data2))
    sublattice_qvals = np.linalg.norm(recon.sub_basis.points, axis=1)
    accept = False
    while not accept:
        cl = ManualClusterer(data, n_maxima=0, mark_qvals=sublattice_qvals)
        triplets = cl.select_triplets()
        assert len(triplets)==2
        pairs = [SpotPair(q1,q2,np.radians(th)) for q1,q2,th in triplets]
        match_results = [recon.sub_basis.match(p) for p in pairs]
        matches_1v = []
        for m in match_results:
            assert m[1] == 'one_vector'
            matches_1v.append(m[0])
        best_basis, best_params, fig = find_best_3d_basis(recon, *matches_1v, delta_range=(-2.0, 2.0), steps=11)
        matches = []
        hits = 0
        for p in recon.all_pairs:
            result = best_basis.match(p)
            if result[1]=='indexed_3d':
                hits += 1
                matches.append(result[0])


        import IPython;IPython.embed()




    exit()
        #recon.read_triplet(*[float(x) for x in l.split()])
    for _ in range(2):
        recon.sub_basis.refine(recon.sublattice_indexed)
        recon.reprocess_all_pairs()

if __name__=="__main__":
    run()
