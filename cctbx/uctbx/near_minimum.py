from __future__ import absolute_import, division, print_function
import numpy as np


def cell_to_metric_tensor(cell):
    """
    Convert cell parameters to metric tensor G.

    Parameters
    ----------
    cell : tuple or unit_cell
        (a, b, c, alpha, beta, gamma) with angles in degrees,
        or an object with .parameters() method

    Returns
    -------
    G : ndarray, shape (3, 3)
        The metric tensor
    """
    # Handle cctbx unit_cell objects
    if hasattr(cell, 'parameters'):
        cell = cell.parameters()

    a, b, c, alpha, beta, gamma = cell
    alpha_rad = np.radians(alpha)
    beta_rad = np.radians(beta)
    gamma_rad = np.radians(gamma)

    cos_alpha = np.cos(alpha_rad)
    cos_beta = np.cos(beta_rad)
    cos_gamma = np.cos(gamma_rad)

    G = np.array([
        [a**2, a*b*cos_gamma, a*c*cos_beta],
        [a*b*cos_gamma, b**2, b*c*cos_alpha],
        [a*c*cos_beta, b*c*cos_alpha, c**2]
    ])
    return G


def metric_tensor_to_cell(G):
    """
    Extract cell parameters from metric tensor.

    Parameters
    ----------
    G : ndarray, shape (3, 3)
        The metric tensor

    Returns
    -------
    cell : tuple
        (a, b, c, alpha, beta, gamma) with angles in degrees
    """
    a = np.sqrt(G[0, 0])
    b = np.sqrt(G[1, 1])
    c = np.sqrt(G[2, 2])

    cos_alpha = np.clip(G[1, 2] / (b * c), -1, 1)
    cos_beta = np.clip(G[0, 2] / (a * c), -1, 1)
    cos_gamma = np.clip(G[0, 1] / (a * b), -1, 1)

    alpha = np.degrees(np.arccos(cos_alpha))
    beta = np.degrees(np.arccos(cos_beta))
    gamma = np.degrees(np.arccos(cos_gamma))

    return (a, b, c, alpha, beta, gamma)


def cell_distance(cell1, cell2):
    """
    Simple L1 distance between cell parameters.

    Parameters
    ----------
    cell1, cell2 : array-like, shape (6,)
        Cell parameters (a, b, c, alpha, beta, gamma)

    Returns
    -------
    float
        Sum of absolute differences
    """
    s = 0.0
    for x1, x2 in zip(cell1, cell2):
        d = x2 - x1
        s += d if d >= 0 else -d
    return s


def lean_min_cell_params(crystal_symmetry):
    """
    Compute a crystal.symmetry's minimum-cell parameters via the C++
    reduction result directly (``red.as_unit_cell()``), skipping the
    ``sgtbx.change_of_basis_op`` construction that
    ``crystal.symmetry.change_of_basis_op_to_minimum_cell()`` performs.

    Bit-identical (to floating-point noise, ~1e-13) to::

        cb = crystal_symmetry.change_of_basis_op_to_minimum_cell()
        crystal_symmetry.unit_cell().change_basis(cb).parameters()

    but ~3x faster over a large database, since it never builds the
    rot_mx/rt_mx/change_of_basis_op/inverse/new_denominators/multiply chain
    of boost.python objects -- only the reduced cell's parameters are
    needed, not a reusable operator. Intended for bulk callers (e.g.
    bulk_lean_min_cell_params below); for a single cell where the
    change-of-basis operator itself is also needed, use
    change_of_basis_op_to_minimum_cell() instead.

    Parameters
    ----------
    crystal_symmetry : crystal.symmetry
        Must expose .space_group() and .unit_cell(), as any crystal.symmetry
        does.

    Returns
    -------
    tuple of 6 floats
        (a, b, c, alpha, beta, gamma) of the minimum cell.
    """
    z2p_op = crystal_symmetry.space_group().z2p_op()
    r_inv = z2p_op.c_inv().r()
    p_cell = crystal_symmetry.unit_cell().change_basis(r_inv.num(), r_inv.den())
    red = p_cell.minimum_reduction()
    return red.as_unit_cell().parameters()


def bulk_lean_min_cell_params(crystal_symmetries):
    """
    Apply lean_min_cell_params to a sequence of crystal.symmetry objects.

    A row whose reduction raises (e.g. a space group genuinely incompatible
    with its cell) is excluded rather than propagating the exception, matching
    the old per-row try/except loop this replaces: such a row simply never
    becomes a bulk-prefilter candidate. The caller is responsible for mapping
    the returned indices back to whatever per-row bookkeeping (pdb codes,
    etc.) it maintains alongside crystal_symmetries.

    Parameters
    ----------
    crystal_symmetries : sequence of crystal.symmetry, length N

    Returns
    -------
    (params, valid_indices) : (ndarray, ndarray)
        params has shape (M, 6): minimum-cell parameters for each row that
        didn't raise, in the same relative order as the input. valid_indices
        has shape (M,), dtype int: params[i] corresponds to
        crystal_symmetries[valid_indices[i]]. M <= N; both are shape (0, 6)
        / (0,) if crystal_symmetries is empty or every row fails.
    """
    params = []
    valid_indices = []
    for i, cs in enumerate(crystal_symmetries):
        try:
            p = lean_min_cell_params(cs)
        except Exception:
            continue
        params.append(p)
        valid_indices.append(i)
    if not params:
        # np.array([]) on an empty list collapses to shape (0,), not
        # (0, 6), which then breaks cells[:, i] downstream -- return the
        # correctly-shaped empty arrays explicitly instead.
        return np.empty((0, 6)), np.empty((0,), dtype=int)
    return np.array(params), np.array(valid_indices, dtype=int)


def cell_to_metric_tensor_vec(cells):
    """
    Vectorized version of cell_to_metric_tensor.

    Parameters
    ----------
    cells : ndarray, shape (N, 6)
        Cell parameters (a, b, c, alpha, beta, gamma), angles in degrees.

    Returns
    -------
    ndarray, shape (N, 3, 3)
        Metric tensor for each input cell.
    """
    cells = np.asarray(cells)
    if cells.shape[0] == 0:
        # An (0,) or (0, 6) empty input both collapse here; cells[:, i]
        # below would otherwise raise IndexError on a bare shape-(0,) array
        # (e.g. from np.array([]) on an empty list of cells).
        return np.empty((0, 3, 3))
    a, b, c, al, be, ga = [cells[:, i] for i in range(6)]
    al, be, ga = np.radians(al), np.radians(be), np.radians(ga)
    G = np.zeros((cells.shape[0], 3, 3))
    G[:, 0, 0] = a * a
    G[:, 1, 1] = b * b
    G[:, 2, 2] = c * c
    G[:, 0, 1] = G[:, 1, 0] = a * b * np.cos(ga)
    G[:, 0, 2] = G[:, 2, 0] = a * c * np.cos(be)
    G[:, 1, 2] = G[:, 2, 1] = b * c * np.cos(al)
    return G


def metric_tensor_to_cell_vec(G):
    """
    Vectorized version of metric_tensor_to_cell.

    Parameters
    ----------
    G : ndarray, shape (..., 3, 3)
        Metric tensor(s); any number of leading batch dimensions.

    Returns
    -------
    ndarray, shape (..., 6)
        Cell parameters (a, b, c, alpha, beta, gamma), angles in degrees.
    """
    a = np.sqrt(G[..., 0, 0])
    b = np.sqrt(G[..., 1, 1])
    c = np.sqrt(G[..., 2, 2])
    cos_al = np.clip(G[..., 1, 2] / (b * c), -1, 1)
    cos_be = np.clip(G[..., 0, 2] / (a * c), -1, 1)
    cos_ga = np.clip(G[..., 0, 1] / (a * b), -1, 1)
    return np.stack([
        a, b, c,
        np.degrees(np.arccos(cos_al)),
        np.degrees(np.arccos(cos_be)),
        np.degrees(np.arccos(cos_ga)),
    ], axis=-1)


def bulk_query_frame_distances(self_params, cbi_near_ops, min_cell_params):
    """
    Vectorized query-frame L1 distance between a query cell and N candidate
    cells (given as minimum-cell parameters), minimized over the query's n
    cached nearly-reduced settings.

    This is a fast lower bound for the distance that
    crystal.symmetry.change_of_basis_op_to_nearest_setting eventually
    reports for each candidate: for each cached setting i (with composed
    change-of-basis operator f_i = cbi_near_ops[i]), the candidate's minimum
    cell is transformed by f_i via the metric-tensor congruence
    G' = M^-T G M^-1 (M = f_i's rotation matrix; this convention, not the
    naive M G M^T or M^T G M, is the one that reproduces
    uctbx.unit_cell.change_basis, verified against boost.python ground
    truth on synthetic test cells), and the L1 distance to self_params is
    taken; the minimum over settings is returned per candidate.

    It is a lower bound, not always bit-identical to the value
    change_of_basis_op_to_nearest_setting's caller would compute, because
    that method picks its "best" setting via a different (minimum-cell-frame)
    criterion first and only then reports whatever query-frame distance that
    already-chosen setting happens to produce. The two selection criteria
    provably coincide when self is already its own minimum cell, and
    empirically coincide for close matches in general; they can diverge for
    poor/far matches, but only with this function's value being <= the true
    value -- so a top-M prefilter ranked by this function's output cannot
    silently drop a true top-N match. See cctbx/command_line/
    search_nearest_cell.py for the intended bulk-prefilter-then-exact-verify
    usage.

    Parameters
    ----------
    self_params : array-like, shape (6,)
        The query's own (untransformed) unit cell parameters.
    cbi_near_ops : sequence of sgtbx.change_of_basis_op, length n
        The composed change-of-basis operators for the query's cached
        nearly-reduced settings, e.g. from
        crystal.symmetry.near_minimum_settings_and_cb_ops().
    min_cell_params : ndarray, shape (N, 6)
        Minimum-cell parameters for each of the N candidate rows, e.g. from
        bulk_lean_min_cell_params().

    Returns
    -------
    ndarray, shape (N,)
        The bulk query-frame distance for each candidate row.
    """
    G_row = cell_to_metric_tensor_vec(min_cell_params)  # (N,3,3)
    Minvs = []
    for f_i in cbi_near_ops:
        r = f_i.c().r()
        M = np.array(r.num(), dtype=float).reshape(3, 3, order='C') / r.den()
        Minvs.append(np.linalg.inv(M))
    Minvs = np.array(Minvs)  # (n,3,3)
    MinvsT = Minvs.transpose(0, 2, 1)  # (n,3,3)

    # G'_{s,r} = Minv_s^T @ G_row[r] @ Minv_s, for all settings s, rows r.
    G_new = np.einsum('sab,rbc,scd->srad', MinvsT, G_row, Minvs)  # (n,N,3,3)
    cells_new = metric_tensor_to_cell_vec(G_new)  # (n,N,6)
    self_params = np.asarray(self_params, dtype=float)
    diffs = np.abs(cells_new - self_params[None, None, :])
    dists_query_frame = diffs.sum(axis=-1)  # (n,N)
    return dists_query_frame.min(axis=0)  # (N,)


def _quick_crossing_check(cell, vec1, vec2_array, delta_length_frac, delta_angle_deg):
    """
    Vectorized check if vec2's can become shorter than vec1 under perturbations.

    Uses analytical sensitivity formulas to determine if a candidate vector
    could become shorter than the reference vector under cell parameter
    perturbations within the given tolerances.

    Parameters
    ----------
    cell : tuple
        (a, b, c, alpha, beta, gamma)
    vec1 : ndarray, shape (3,)
        Reference vector as integer indices
    vec2_array : ndarray, shape (N, 3)
        Candidate vectors as integer indices
    delta_length_frac : float
        Fractional length perturbation (e.g., 0.03 means +/-3%)
    delta_angle_deg : float
        Angle perturbation in degrees

    Returns
    -------
    can_cross : ndarray, shape (N,)
        Boolean array indicating which vec2's can become shorter than vec1
    """
    a, b, c, alpha, beta, gamma = cell

    cos_alpha = np.cos(np.radians(alpha))
    cos_beta = np.cos(np.radians(beta))
    cos_gamma = np.cos(np.radians(gamma))
    sin_alpha = np.sin(np.radians(alpha))
    sin_beta = np.sin(np.radians(beta))
    sin_gamma = np.sin(np.radians(gamma))

    G = cell_to_metric_tensor(cell)

    # Current squared lengths
    len1_sq = vec1 @ G @ vec1
    len2_sq = np.sum((vec2_array @ G) * vec2_array, axis=1)
    current_diff = len1_sq - len2_sq

    # Extract components
    n1_1, n1_2, n1_3 = vec1[0], vec1[1], vec1[2]
    n2_1 = vec2_array[:, 0]
    n2_2 = vec2_array[:, 1]
    n2_3 = vec2_array[:, 2]

    # Sensitivity of |v|^2 to fractional axis changes
    # d|v|^2/d(delta_a) = 2a * n1 * (n1*a + n2*b*cos_gamma + n3*c*cos_beta)
    sens1_a = 2 * a * n1_1 * (n1_1 * a + n1_2 * b * cos_gamma + n1_3 * c * cos_beta)
    sens1_b = 2 * b * n1_2 * (n1_2 * b + n1_1 * a * cos_gamma + n1_3 * c * cos_alpha)
    sens1_c = 2 * c * n1_3 * (n1_3 * c + n1_1 * a * cos_beta + n1_2 * b * cos_alpha)

    sens2_a = 2 * a * n2_1 * (n2_1 * a + n2_2 * b * cos_gamma + n2_3 * c * cos_beta)
    sens2_b = 2 * b * n2_2 * (n2_2 * b + n2_1 * a * cos_gamma + n2_3 * c * cos_alpha)
    sens2_c = 2 * c * n2_3 * (n2_3 * c + n2_1 * a * cos_beta + n2_2 * b * cos_alpha)

    axis_contrib = (np.abs(sens1_a - sens2_a) +
                    np.abs(sens1_b - sens2_b) +
                    np.abs(sens1_c - sens2_c)) * delta_length_frac

    # Angle sensitivities
    delta_angle_rad = np.radians(delta_angle_deg)

    sens1_gamma = 2 * n1_1 * n1_2 * a * b * sin_gamma
    sens1_beta = 2 * n1_1 * n1_3 * a * c * sin_beta
    sens1_alpha = 2 * n1_2 * n1_3 * b * c * sin_alpha

    sens2_gamma = 2 * n2_1 * n2_2 * a * b * sin_gamma
    sens2_beta = 2 * n2_1 * n2_3 * a * c * sin_beta
    sens2_alpha = 2 * n2_2 * n2_3 * b * c * sin_alpha

    angle_contrib = (np.abs(sens1_gamma - sens2_gamma) +
                     np.abs(sens1_beta - sens2_beta) +
                     np.abs(sens1_alpha - sens2_alpha)) * delta_angle_rad

    max_change = axis_contrib + angle_contrib
    can_cross = np.abs(current_diff) <= max_change

    return can_cross


def find_near_minimum_settings_detail(
    cell,
    length_tolerance=0.03,
    angle_tolerance=3.0,
    max_search_index=10,
):
    """
    Find nearly-reduced settings using efficient vectorized filtering.

    This function finds all basis choices that are "nearly reduced" - i.e.,
    bases where the vectors could become the shortest under small perturbations
    to the cell parameters. This is useful for comparing cells near reduction
    boundaries.

    Algorithm:
    1. Generate all lattice vectors up to +/-max_search_index
    2. Find shortest vector s1
    3. Use analytical crossing check to find vectors that could become
       shorter than s1 under the given perturbations -> a-candidates
    4. For each a-candidate, repeat for b (non-parallel) and c (non-coplanar)
    5. Construct right-handed bases and return transformation matrices

    Parameters
    ----------
    cell : tuple or unit_cell
        (a, b, c, alpha, beta, gamma) - axes in Angstrom, angles in degrees
    length_tolerance : float
        Fractional tolerance for length perturbations (e.g., 0.03 = 3%)
    angle_tolerance : float
        Tolerance for angle perturbations in degrees
    max_search_index : int
        Search lattice vectors up to +/-max_search_index

    Returns
    -------
    list of dict
        Each dict contains:
        - 'P': transformation matrix (3x3 integer array)
        - 'cell': transformed cell parameters (tuple of 6 floats)
        - 'G': transformed metric tensor (3x3 array)
    """
    # Handle cctbx unit_cell objects
    if hasattr(cell, 'parameters'):
        cell = cell.parameters()

    G = cell_to_metric_tensor(cell)

    # Generate all lattice vectors
    grid = np.mgrid[-max_search_index:max_search_index+1,
                    -max_search_index:max_search_index+1,
                    -max_search_index:max_search_index+1]
    indices = grid.reshape(3, -1).T
    mask = ~np.all(indices == 0, axis=1)
    indices = indices[mask]

    lengths = np.sqrt(np.sum((indices @ G) * indices, axis=1))
    sort_order = np.argsort(lengths)
    indices = indices[sort_order]
    lengths = lengths[sort_order]

    # Find a-candidates
    s1_vec = indices[0]
    can_cross = _quick_crossing_check(
        cell, s1_vec, indices, length_tolerance, angle_tolerance
    )
    a_candidates = indices[can_cross]

    results = []

    for a_vec in a_candidates:
        # Filter to non-parallel vectors
        crosses = np.cross(a_vec, indices)
        cross_norms = np.linalg.norm(crosses, axis=1)
        non_parallel_indices = indices[cross_norms > 1e-10]
        if len(non_parallel_indices) == 0:
            continue

        s2_vec = non_parallel_indices[0]
        can_cross_b = _quick_crossing_check(
            cell, s2_vec, non_parallel_indices, length_tolerance, angle_tolerance
        )
        b_candidates = non_parallel_indices[can_cross_b]

        for b_vec in b_candidates:
            # Filter to non-coplanar vectors
            test_matrices = np.stack([
                np.tile(a_vec, (len(indices), 1)),
                np.tile(b_vec, (len(indices), 1)),
                indices
            ], axis=2)
            dets = np.linalg.det(test_matrices)
            is_non_coplanar = np.abs(dets) > 1e-10
            non_coplanar_indices = indices[is_non_coplanar]

            if len(non_coplanar_indices) == 0:
                continue

            s3_vec = non_coplanar_indices[0]
            can_cross_c = _quick_crossing_check(
                cell, s3_vec, non_coplanar_indices, length_tolerance, angle_tolerance
            )
            c_candidates = non_coplanar_indices[can_cross_c]

            for c_vec in c_candidates:
                P = np.column_stack([a_vec, b_vec, c_vec]).astype(int)
                det = int(np.round(np.linalg.det(P)))

                # Only right-handed bases
                if det != 1:
                    continue

                G_new = P.T @ G @ P
                cell_new = metric_tensor_to_cell(G_new)

                results.append({
                    'P': P,
                    'cell': cell_new,
                    'G': G_new,
                })

    # Sort by sum of basis vector lengths
    results.sort(key=lambda x: sum(x['cell'][:3]))

    return results


def find_near_minimum_settings(
    cell,
    length_tolerance=0.03,
    angle_tolerance=3.0,
    max_search_index=10,
    test_multiples=False
):
    """
    Find nearly-reduced settings across order-2 super- and subcells.

    This function extends find_near_minimum_settings by searching not only the
    input cell, but also 7 order-2 supercells and 7 order-2 subcells (14 total
    additional cells). This helps find alternative settings that might be missed
    due to reduction boundaries when the true cell might be a multiple or
    submultiple of the measured cell.

    Parameters
    ----------
    cell : tuple or unit_cell
        (a, b, c, alpha, beta, gamma) - axes in Angstrom, angles in degrees
    length_tolerance : float
        Fractional tolerance for length perturbations (e.g., 0.03 = 3%)
    angle_tolerance : float
        Tolerance for angle perturbations in degrees
    max_search_index : int
        Search lattice vectors up to +/-max_search_index

    Returns
    -------
    list of dict
        Each dict contains:
        - 'P': transformation matrix (3x3 integer or float array) including
               the super/subcell transformation
        - 'cell': transformed cell parameters (tuple of 6 floats)
        - 'G': transformed metric tensor (3x3 array)
    """
    # Handle cctbx unit_cell objects
    if hasattr(cell, 'parameters'):
        cell = cell.parameters()

    # Define 7 standard order-2 supercell transformations (det = 2)
    supercell_transforms = [
        np.array([[2, 0, 0], [0, 1, 0], [0, 0, 1]]),  # double a
        np.array([[1, 0, 0], [0, 2, 0], [0, 0, 1]]),  # double b
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]]),  # double c
        np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]]),  # C-face diagonal
        np.array([[1, 0, 1], [0, 1, 0], [-1, 0, 1]]),  # B-face diagonal
        np.array([[1, 0, 0], [0, 1, 1], [0, -1, 1]]),  # A-face diagonal
        np.array([[1, 1, 1], [-1, 1, 0], [-1, 0, 1]]),  # body diagonal
    ]

    # Subcell transforms are the inverses of supercell transforms
    # (det = 1/2, map from original to smaller cell)
    subcell_transforms = [np.linalg.inv(S) for S in supercell_transforms]

    all_results = []

    # Process original cell
    original_results = find_near_minimum_settings_detail(
        cell, length_tolerance, angle_tolerance, max_search_index
    )
    all_results.extend(original_results)

    G_original = cell_to_metric_tensor(cell)

    if test_multiples:
        # Process supercells
        for S in supercell_transforms:
            # Transform to supercell: G_super = S^T @ G @ S
            G_super = S.T @ G_original @ S
            cell_super = metric_tensor_to_cell(G_super)

            # Find near-minimum settings in the supercell
            super_results = find_near_minimum_settings_detail(
                cell_super, length_tolerance, angle_tolerance, max_search_index
            )

            # Compose transformations: P_total = S @ P_super
            for result in super_results:
                P_total = S @ result['P']
                G_total = P_total.T @ G_original @ P_total
                cell_total = metric_tensor_to_cell(G_total)

                all_results.append({
                    'P': P_total,
                    'cell': cell_total,
                    'G': G_total,
                })

        # Process subcells
        for T in subcell_transforms:
            # Transform to subcell: G_sub = T^T @ G @ T
            G_sub = T.T @ G_original @ T
            cell_sub = metric_tensor_to_cell(G_sub)

            # Find near-minimum settings in the subcell
            sub_results = find_near_minimum_settings_detail(
                cell_sub, length_tolerance, angle_tolerance, max_search_index
            )

            # Compose transformations: P_total = T @ P_sub
            for result in sub_results:
                P_total = T @ result['P']
                G_total = P_total.T @ G_original @ P_total
                cell_total = metric_tensor_to_cell(G_total)

                all_results.append({
                    'P': P_total,
                    'cell': cell_total,
                    'G': G_total,
                })

    # Sort by sum of basis vector lengths
    all_results.sort(key=lambda x: sum(x['cell'][:3]))

    return all_results
