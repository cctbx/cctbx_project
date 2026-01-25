from __future__ import division
import numpy as np
import copy

from dials.array_family import flex
from cctbx import uctbx, miller
from scitbx import matrix
from scipy.optimize import minimize


class PowderGeometryRefiner:
    """
    Refines detector geometry using powder diffraction d-spacings.

    Uses a 5-parameter model:
    - dist: distance from origin to detector plane along normal (mm)
    - shift1: shift along detector fast axis (mm)
    - shift2: shift along detector slow axis (mm)
    - tau2: rotation around fast axis (mrad) - tilt
    - tau3: rotation around slow axis (mrad) - tilt

    NOT refined: tau1 (rotation around normal) - indeterminate from powder
    """

    def __init__(self, experiments, reflections, params):
        self.experiments = experiments
        self.reflections = reflections
        self.params = params

        # Compute reference d-spacings from unit_cell and space_group if provided
        if params.unit_cell is not None and params.space_group is not None:
            print(f"Computing d-spacings from unit_cell={params.unit_cell} "
                  f"and space_group={params.space_group.info()}")

            # Get beam and detector for resolution calculation
            beam = experiments[0].beam
            detector = experiments[0].detector
            d_min = params.d_min
            d_max = params.d_max

            # Average unit cell over space group
            unit_cell = params.space_group.average_unit_cell(params.unit_cell)

            # Generate Miller indices within resolution range
            generator = miller.index_generator(unit_cell, params.space_group.type(), False, d_min)
            indices = generator.to_array()

            # Compute d-spacings and filter by resolution range
            all_spacings = unit_cell.d(indices)
            spacings_in_range = flex.sorted([d for d in all_spacings if d_min <= d <= d_max])

            self.reference_d = np.array(spacings_in_range)
            print(f"Generated {len(self.reference_d)} reference d-spacings in range "
                  f"[{d_min:.3f}, {d_max:.3f}] A")
        else:
            self.reference_d = np.array(params.reference_d_spacings)

        # Store initial detector state
        self.detector = experiments[0].detector
        self.initial_state = self._get_detector_state()

        # Track which parameters to refine
        self.refine_distance = params.refine.distance
        self.refine_shift = params.refine.shift
        self.refine_tilt = params.refine.tilt

        # Build parameter vector and bounds
        self._setup_parameters()

        # Prepare reflections
        self._prepare_reflections()

    def _get_detector_state(self):
        """Extract fast, slow, origin, and center from detector panel."""
        panel = self.detector[0]
        fast = matrix.col(panel.get_fast_axis())
        slow = matrix.col(panel.get_slow_axis())
        origin = matrix.col(panel.get_origin())
        normal = fast.cross(slow)

        # Compute panel center
        size = panel.get_image_size()
        pixel_size = panel.get_pixel_size()
        center = origin + (size[0]/2 * pixel_size[0]) * fast + (size[1]/2 * pixel_size[1]) * slow

        return {
            'fast': fast,
            'slow': slow,
            'origin': origin,
            'normal': normal,
            'center': center,
        }

    def _setup_parameters(self):
        """Set up parameter vector based on what is being refined."""
        # Initial values: [dist, shift1, shift2, tau2, tau3]
        # All start at 0 (representing no change from initial geometry)
        self.param_names = []
        self.initial_params = []
        self.bounds = []

        if self.refine_distance:
            self.param_names.append('dist')
            self.initial_params.append(0.0)
            self.bounds.append((-50.0, 50.0))  # +/- 50 mm

        if self.refine_shift:
            self.param_names.append('shift1')
            self.initial_params.append(0.0)
            self.bounds.append((-10.0, 10.0))  # +/- 10 mm

            self.param_names.append('shift2')
            self.initial_params.append(0.0)
            self.bounds.append((-10.0, 10.0))  # +/- 10 mm

        if self.refine_tilt:
            self.param_names.append('tau2')
            self.initial_params.append(0.0)
            self.bounds.append((-50.0, 50.0))  # +/- 50 mrad

            self.param_names.append('tau3')
            self.initial_params.append(0.0)
            self.bounds.append((-50.0, 50.0))  # +/- 50 mrad

        self.initial_params = np.array(self.initial_params)

    def _prepare_reflections(self):
        """Filter reflections by d-range and proximity to reference d-spacings."""
        refls = self.reflections

        # Compute initial d-spacings to filter
        refls.centroid_px_to_mm(self.experiments)
        refls.map_centroids_to_reciprocal_space(self.experiments)
        d_star_sq = flex.pow2(refls['rlp'].norms())
        refls['d'] = uctbx.d_star_sq_as_d(d_star_sq)

        # Filter by d-range
        d_min = self.params.d_min
        d_max = self.params.d_max
        sel = (refls['d'] >= d_min) & (refls['d'] <= d_max)
        refls_in_range = refls.select(sel)

        print(f"Found {len(refls_in_range)} reflections in d-range "
              f"[{d_min:.3f}, {d_max:.3f}] A")

        # Filter by proximity to reference d-spacings (in inverse angstroms)
        max_dist = self.params.max_distance_inv_ang
        dvals = refls_in_range['d'].as_numpy_array()
        dvals_inv = 1.0 / dvals
        ref_d_inv = 1.0 / self.reference_d

        # For each reflection, find minimum distance to any reference
        min_distances = np.min(
            np.abs(dvals_inv[:, np.newaxis] - ref_d_inv[np.newaxis, :]),
            axis=1
        )
        close_to_ref = min_distances <= max_dist
        sel_close = flex.bool(close_to_ref.tolist())
        self.refls_filtered = refls_in_range.select(sel_close)

        print(f"Using {len(self.refls_filtered)} reflections within "
              f"{max_dist:.4f} inv. A of reference d-spacings")

        # Store pixel positions for later use
        self.xyzobs_px = self.refls_filtered['xyzobs.px.value']
        self.panels = self.refls_filtered['panel']
        self.ids = self.refls_filtered['id']

    def apply_params(self, x):
        """Apply parameter vector to detector geometry."""
        # Parse parameter vector
        params_dict = {}
        for i, name in enumerate(self.param_names):
            params_dict[name] = x[i]

        # Convert to Python floats for scitbx matrix compatibility
        dist = float(params_dict.get('dist', 0.0))
        shift1 = float(params_dict.get('shift1', 0.0))
        shift2 = float(params_dict.get('shift2', 0.0))
        tau2 = float(params_dict.get('tau2', 0.0)) / 1000.0  # mrad to rad
        tau3 = float(params_dict.get('tau3', 0.0)) / 1000.0  # mrad to rad

        # Get initial state
        d1 = self.initial_state['fast']
        d2 = self.initial_state['slow']
        dn = self.initial_state['normal']
        origin = self.initial_state['origin']
        center = self.initial_state['center']

        # Apply rotations around the panel CENTER (not origin)
        # tau2: rotation around fast axis (d1)
        # tau3: rotation around slow axis (d2)
        if abs(tau2) > 1e-10 or abs(tau3) > 1e-10:
            # Use axis-angle rotation around the actual detector axes
            # This correctly handles the left-handed detector frame
            R2 = d1.axis_and_angle_as_r3_rotation_matrix(tau2, deg=False)
            R3 = d2.axis_and_angle_as_r3_rotation_matrix(tau3, deg=False)

            # Combined rotation: first tau2 around fast, then tau3 around slow
            R = R3 * R2

            # Apply rotation to axes
            d1_new = R * d1
            d2_new = R * d2

            # Rotate origin around center to preserve panel center position
            origin_rotated = center + R * (origin - center)
        else:
            d1_new = d1
            d2_new = d2
            origin_rotated = origin

        # Apply shifts and distance change (relative to initial axes)
        new_origin = (origin_rotated +
                      shift1 * d1 +
                      shift2 * d2 +
                      dist * dn)

        # Update detector panel
        panel = self.detector[0]
        panel.set_frame(
            d1_new.elems,
            d2_new.elems,
            new_origin.elems
        )

    def compute_dvals(self):
        """Compute d-spacings for all reflections with current geometry."""
        # Reset mm coordinates to force recalculation
        if 'xyzobs.mm.value' in self.refls_filtered:
            del self.refls_filtered['xyzobs.mm.value']
        if 'rlp' in self.refls_filtered:
            del self.refls_filtered['rlp']

        self.refls_filtered.centroid_px_to_mm(self.experiments)
        self.refls_filtered.map_centroids_to_reciprocal_space(self.experiments)
        d_star_sq = flex.pow2(self.refls_filtered['rlp'].norms())
        self.dvals = uctbx.d_star_sq_as_d(d_star_sq)
        return self.dvals

    def find_nearest_reference(self, d_obs):
        """Find closest reference d-spacing."""
        distances = np.abs(self.reference_d - d_obs)
        return self.reference_d[np.argmin(distances)]

    def objective(self, x):
        """Compute sum of squared residuals."""
        self.apply_params(x)
        dvals = self.compute_dvals()

        residuals = []
        for d_obs in dvals:
            d_ref = self.find_nearest_reference(d_obs)
            residuals.append(d_obs - d_ref)

        return np.sum(np.array(residuals) ** 2)

    def run(self):
        """Run refinement and return results."""
        if len(self.param_names) == 0:
            print("No parameters selected for refinement!")
            return None

        # Compute initial objective
        initial_obj = self.objective(self.initial_params)
        print(f"\nInitial objective: {initial_obj:.6f}")
        print(f"Initial RMS residual: {np.sqrt(initial_obj / len(self.refls_filtered)):.6f} A")

        # Run minimization using Powell method (derivative-free, more robust)
        print(f"\nRefining {len(self.param_names)} parameters: {self.param_names}")
        result = minimize(
            self.objective,
            x0=self.initial_params,
            method='Powell',
            options={'maxiter': 100, 'disp': True, 'xtol': 0.001, 'ftol': 0.0001}
        )

        # Apply final parameters
        self.apply_params(result.x)
        final_dvals = self.compute_dvals()

        # Compute final statistics
        final_obj = result.fun
        print(f"\nFinal objective: {final_obj:.6f}")
        print(f"Final RMS residual: {np.sqrt(final_obj / len(self.refls_filtered)):.6f} A")

        # Report parameter changes
        print("\nRefined parameters:")
        for i, name in enumerate(self.param_names):
            value = result.x[i]
            unit = "mm" if name in ['dist', 'shift1', 'shift2'] else "mrad"
            print(f"  {name}: {value:.4f} {unit}")

        return result

    def get_refined_experiments(self):
        """Return experiments with refined detector."""
        return self.experiments

    def report_geometry_changes(self):
        """Print summary of geometry changes."""
        panel = self.detector[0]
        new_origin = matrix.col(panel.get_origin())
        new_fast = matrix.col(panel.get_fast_axis())
        new_slow = matrix.col(panel.get_slow_axis())

        old_origin = self.initial_state['origin']
        old_fast = self.initial_state['fast']
        old_slow = self.initial_state['slow']

        origin_delta = new_origin - old_origin

        # Compute tilt changes from axis differences
        # tau2 (tilt around fast) shows in slow axis z-component
        # tau3 (tilt around slow) shows in fast axis z-component
        tau2_change = np.arcsin(new_slow[2]) - np.arcsin(old_slow[2])
        tau3_change = np.arcsin(-new_fast[2]) - np.arcsin(-old_fast[2])

        print("\nGeometry changes:")
        print(f"  Origin shift: ({origin_delta[0]:.4f}, {origin_delta[1]:.4f}, {origin_delta[2]:.4f}) mm")
        print(f"  Tilt changes: tau2={tau2_change*1000:.3f} mrad, tau3={tau3_change*1000:.3f} mrad")
