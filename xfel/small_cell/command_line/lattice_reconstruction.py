from dials.util import Sorry
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util import show_mail_handle_errors
from dials.array_family import flex
from iotbx.phil import parse
import numpy as np
from dxtbx.model import ExperimentList
from sklearn.cluster import DBSCAN, HDBSCAN
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import RectangleSelector
from sklearn.mixture import GaussianMixture
from scipy.optimize import minimize
import hashlib
import pickle
from cctbx.crystal import symmetry
from cctbx import sgtbx
import math
from dxtbx import flumpy


# Top-level phil scope
phil_scope = parse("""
mode = *initialize continue
  .type = choice
max_cycles = 1
  .type = int
  .help = Number of recruitment/refinement cycles to perform
triplets {
  d_min = 2.0
    .type = float
    .help = Compute d1,d2,theta triplets to this resolution
  load_npz = None
    .type = str
    .help = Load precomputed triplet array from file
  save_npz = None
    .type = str
    .help = Save computed triplet array to file
}
initialization {
  triangle_finding {
    mahalanobis_threshold = 2.5
      .type = float
      .help = Distance cutoff for triangle spec matching
    method = interactive clustering *box
      .type = choice
    save_clusters = None
      .type = str
      .help = Save clusters to pickle file
    load_clusters = None
      .type = str
      .help = Load clusters from pickle file
    interactive {
      min_points = 100
        .type = int
        .help = "Do not attempt clustering on ranges with fewer than this many points."
      cluster_safety_margin = 3.0
        .type = float
        .help = Multiplier for selected ranges when filtering data for GMM fit.
      gmm_params {
        use = False
          .type = bool
        q1_centers = None
          .type = floats
        q2_centers = None
          .type = floats
        radii = None
          .type = floats
        theta_centers = None
          .type = floats
        theta_widths = None
          .type = floats
      }
    }
    box {
      d1_range = None
        .type = floats(size=2)
        .help = "Resolution range (d_max, d_min) for reference spot"
      d2_range = None
        .type = floats(size=2)
        .help = "Resolution range for second spot"
      theta_range = None
        .type = floats(size=2)
        .help = "Angular range between spots (degrees)"
    }
  }
}
input {
  aligned_experiments = None
    .type = path
    .help = "Partially reconstructed experiment list (for continue mode)"
  aligned_reflections = None
    .type = path
    .help = "Partially reconstructed reflection list (for continue mode)"
}
output {
  experiments = aligned_experiments.expt
    .type = path
  reflections = aligned_reflections.refl
    .type = path
}
debug {
  unit_cell = None
    .type = unit_cell
    .help = sldkfj
}
""")

help_message = """
Nobody knows what this does or why it was created
"""

class ReconstructionManager:
  def __init__(self, experiments, reflections, triplets, params):
    self.reconstructions = []
    self.relationships = {}
    self.experiments = experiments
    self.reflections = reflections
    self.triplets = triplets
    self.params = params

  def initialize(self, cluster_results):
    """Interactive selection of initial reconstruction."""
    candidates = []
    for cluster in cluster_results:
      # Create left-handed reconstruction
      spec = TriangleSpec.from_cluster(cluster, self.params, hand=-1)
      lr_left = spec.reconstruct(self.triplets, self.experiments, self.reflections)
      candidates.append(lr_left)

      # Create right-handed reconstruction
      spec = TriangleSpec.from_cluster(cluster, self.params, hand=1)
      lr_right = spec.reconstruct(self.triplets, self.experiments, self.reflections)
      candidates.append(lr_right)

    current_id = 0
    while current_id < len(candidates):
      print(f"\nCandidate {current_id}: {len(candidates[current_id].reflections)} reflections")

      result = self.plot_2d_radial(
        reconstruction=candidates[current_id],
        mode='screen',
        navigation_info={
          'current_id': current_id,
          'total': len(candidates),
          'candidates': candidates
        }
      )

      if result['accepted']:
        self.reconstructions.append(candidates[current_id])
        return

      if result['next_id'] is not None:
        current_id = result['next_id']


  def plot_2d_radial(self, reconstruction=None, reconstruction_id=None, mode='view',
                     navigation_info=None):
    """Plot 2D radial view of a reconstruction with mode-specific interactions.

    Args:
      reconstruction: LatticeReconstruction to plot (for viewing candidates)
      reconstruction_id: Index in self.reconstructions (for viewing stored reconstructions)
      mode: Interaction mode
      navigation_info: Optional dict with:
        'current_id': Current position in sequence
        'total': Total number of items
        'candidates': List of reconstructions (if navigating through candidates)
    """
    if (reconstruction is None) == (reconstruction_id is None):
      raise Sorry("Must provide exactly one of reconstruction or reconstruction_id")

    if reconstruction is None:
      reconstruction = self.reconstructions[reconstruction_id]

    # Initialize result dictionary
    result = {
      'next_id': None,
      'accepted': None,    # for screen mode
      'boxes': [],         # for augment mode
      'points': []         # for select mode
    }
    # Set up single plot
    fig, ax = plt.subplots(figsize=(10,8))

    # Plot points
    for i_exp in range(len(reconstruction.experiments)):
      exp_sel = reconstruction.reflections['id'] == i_exp
      exp_refs = reconstruction.reflections.select(exp_sel)
      q_vecs = exp_refs['q.aligned']
      x, y = q_vecs.parts()[0], q_vecs.parts()[1]

      # Calculate theta using arctan2(y,x) to match theta_refine
      theta = np.rad2deg(np.arctan2(y, x))
      q_mags = np.sqrt(q_vecs.parts()[0]**2 +
                      q_vecs.parts()[1]**2 +
                      q_vecs.parts()[2]**2)

      # Plot spots
      ax.scatter(theta, q_mags, alpha=0.2, s=3, c='blue')

      # Plot reference spots
      i1, i2 = reconstruction.reference_spots[i_exp]
      ref_theta = np.rad2deg(np.arctan2([y[i1], y[i2]], [x[i1], x[i2]]))
      ref_q = q_mags[[i1, i2]]
      ax.scatter(ref_theta, ref_q, c=['red', 'green'], s=50)

    # Configure axis
    ax.set_xlim(-180, 180)
    ax.yaxis.set_major_formatter(tick.FuncFormatter(
      lambda q, _: f"{1/q:.1f}" if q > 0 else ""
    ))
    ax.grid(True, linestyle=':')
    ax.set_xlabel('Azimuthal angle θ (degrees)')
    ax.set_ylabel('Resolution d (Å)')

    # Mode-specific title
    title = "Left-handed" if reconstruction.hand == -1 else "Right-handed"
    ax.set_title(title)

    # Basic plotting code...

    def on_key(event):
      if event.key == 'j':  # Previous reconstruction
        if navigation_info is not None:
          if navigation_info['current_id'] > 0:
            result['next_id'] = navigation_info['current_id'] - 1
            plt.close(fig)
        elif reconstruction_id is not None and reconstruction_id > 0:
          result['next_id'] = reconstruction_id - 1
          plt.close(fig)

      elif event.key == 'k':  # Next reconstruction
        if navigation_info is not None:
          if navigation_info['current_id'] < navigation_info['total'] - 1:
            result['next_id'] = navigation_info['current_id'] + 1
            plt.close(fig)
        elif reconstruction_id is not None and reconstruction_id < len(self.reconstructions) - 1:
          result['next_id'] = reconstruction_id + 1
          plt.close(fig)

      else:
        # Mode-specific key handling
        if mode == 'screen':
          if event.key in ['y', 'n']:
            result['accepted'] = (event.key == 'y')
            plt.close(fig)

        elif mode == 'augment':
          if event.key == 'd':  # Done selecting
            plt.close(fig)

        elif mode == 'select':
          if event.key == 'd':  # Done selecting
            plt.close(fig)

    # Mode-specific mouse handling
    if mode == 'augment':
      def on_select(eclick, erelease):
        # Handle box selection...
        pass

    elif mode == 'select':
      def on_click(event):
        # Handle point selection...
        pass

    # Connect event handlers
    fig.canvas.mpl_connect('key_press_event', on_key)
    if mode == 'augment':
      rect_selector = RectangleSelector(...)
    elif mode == 'select':
      fig.canvas.mpl_connect('button_press_event', on_click)

    plt.show()
    return result


  def augment(self):
    """Allow selection of boxes to create new reconstruction."""
    # Show menu
    print("\nCurrent reconstructions:")
    for i, recon in enumerate(self.reconstructions):
      print(f"{i}: {len(recon.reflections)} reflections")

    cmd = input("\nEnter reconstruction number (or q to finish): ")
    if cmd == 'q':
      return True

    try:
      recon_id = int(cmd)
      if self.select_boxes(recon_id):
        return False  # Continue augmenting
    except ValueError:
      print("Invalid input")

    return False

  def select_boxes(self, recon_id):
    """Allow selection of two boxes and create new reconstruction."""
    recon = self.reconstructions[recon_id]
    # Modified plot_2d_radial to allow box selection
    boxes = recon.plot_2d_radial(select_boxes=True)
    if len(boxes) != 2:
      return False

    # Determine hand of selected points
    hand = self.determine_hand(boxes[0], boxes[1])
    # Create TriangleSpec with opposite hand
    spec = TriangleSpec.from_boxes(boxes[0], boxes[1], hand=-hand)

    # Create new reconstruction
    new_recon = spec.reconstruct(self.triplets, self.experiments, self.reflections)
    self.reconstructions.append(new_recon)
    self.relationships[(recon_id, len(self.reconstructions)-1)] = {
      'parent_boxes': boxes,
      'hand': -hand
    }
    return True


class LatticeReconstruction:
  def __init__(self, params):
    self.experiments = ExperimentList()
    self.reflections = flex.reflection_table()
    self.orientations = {}
    self.used_triangles = []  # List of (d1,d2,theta) tuples
    self.clusters = None      # Will store cluster analysis results
    self.reference_spots = {}
    self.params = params
    self.hand = None # Will be set to -1 or 1
    
  def add(self, expt, refl, i1, i2, xyz1=None, xyz2=None):
    """Add transformed reflections from one experiment

    Parameters:
      expt: single experiment
      refl: reflection table for the experiment
      i1, i2: indices of reference reflections in refl
      xyz1, xyz2: target coordinates for reference reflections after transformation
    """
    # First get reciprocal space vectors for all reflections
    s0 = expt.beam.get_s0()
    s1_vectors = refl['s1'].as_numpy_array()
    q_vectors = s1_vectors - s0

    # Get reference vectors
    q1 = q_vectors[i1]
    q2 = q_vectors[i2]
    k1 = np.linalg.norm(q1)
    k2 = np.linalg.norm(q2)

    if xyz1 is None or xyz2 is None:
      # Get magnitudes
      cos_theta = np.dot(q1, q2)/(k1*k2)
      theta = np.arccos(min(1.0, max(-1.0, cos_theta)))

      # Set xyz2 along x-axis
      xyz2 = np.array([k2, 0, 0])
      # Set xyz1 in x-y plane at angle theta from xyz2
      xyz1 = np.array([k1*np.cos(theta), k1*np.sin(theta), 0])

    # Find rotation matrix
    # First align q1 to xyz1
    v1 = q1/k1  # unit vector
    t1 = xyz1/np.linalg.norm(xyz1)  # unit vector
    rot1_axis = np.cross(v1, t1)
    if np.linalg.norm(rot1_axis) < 1e-8:
      # Vectors already aligned or anti-parallel
      rot1 = np.eye(3) if np.dot(v1, t1) > 0 else -np.eye(3)
    else:
      rot1_axis = rot1_axis/np.linalg.norm(rot1_axis)
      cos_alpha = np.dot(v1, t1)
      alpha = np.arccos(min(1.0, max(-1.0, cos_alpha)))
      # Rodrigues rotation formula
      K = np.array([[0, -rot1_axis[2], rot1_axis[1]],
                    [rot1_axis[2], 0, -rot1_axis[0]],
                    [-rot1_axis[1], rot1_axis[0], 0]])
      rot1 = np.eye(3) + np.sin(alpha)*K + (1-np.cos(alpha))*K.dot(K)

    # Now rotate q2 around xyz1 to get as close as possible to xyz2
    q2_rot = rot1.dot(q2)
    # Project q2_rot onto plane perpendicular to xyz1
    n1 = xyz1/np.linalg.norm(xyz1)
    q2_proj = q2_rot - np.dot(q2_rot, n1)*n1
    t2_proj = xyz2 - np.dot(xyz2, n1)*n1

    # Normalize the projections
    q2_proj = q2_proj/np.linalg.norm(q2_proj)
    t2_proj = t2_proj/np.linalg.norm(t2_proj)

    # Use cross product to determine rotation direction
    cross_prod = np.cross(q2_proj, t2_proj)
    rotation_direction = np.sign(np.dot(cross_prod, n1))

    # Find angle between projections
    cos_beta = np.dot(q2_proj, t2_proj)
    beta = rotation_direction * np.arccos(min(1.0, max(-1.0, cos_beta)))

    # Create rotation matrix around xyz1 (n1)
    if np.linalg.norm(q2_proj) < 1e-8 or np.linalg.norm(t2_proj) < 1e-8:
      rot2 = np.eye(3)
    else:
      rot2_axis = n1
      K = np.array([[0, -rot2_axis[2], rot2_axis[1]],
                    [rot2_axis[2], 0, -rot2_axis[0]],
                    [-rot2_axis[1], rot2_axis[0], 0]])
      rot2 = np.eye(3) + np.sin(beta)*K + (1-np.cos(beta))*K.dot(K)

    # Combined rotation
    rot = rot2.dot(rot1)

    # Apply to all vectors
    q_vectors_transformed = np.array([rot.dot(q) for q in q_vectors])

    # Determine handedness from transformed beam direction
    s0_transformed = rot.dot(s0)
    hand = int(s0_transformed[2] > 0)  # 1 if beam points in +z, 0 if -z

    # Add experiment and reflections to reconstruction
    new_id = len(self.experiments)
    self.experiments.append(expt)

    # Copy reflection table and add transformed coordinates
    new_refl = refl.copy()
    new_refl['transformed_xyz'] = flex.vec3_double(q_vectors_transformed)
    new_refl['q.aligned'] = flex.vec3_double(q_vectors_transformed)
    new_refl['hand'] = flex.int(len(refl), hand)

#    # Create detector-aligned coordinates (xyzobs.aligned)
#    # Project q-vectors onto plane perpendicular to beam
#    beam_dir = s0_transformed/np.linalg.norm(s0_transformed)
#    detector_coords = []
#    for q in q_vectors_transformed:
#      # Project out beam direction
#      q_proj = q - np.dot(q, beam_dir)*beam_dir
#      detector_coords.append([q_proj[0], q_proj[1], q_proj[2]])
#    new_refl['xyzobs.aligned'] = flex.vec3_double(detector_coords)

    new_refl['hand'] = flex.int(len(refl), hand)

    new_refl['id'] = flex.int(len(refl), new_id)

    if len(self.reflections) == 0:
      self.reflections = new_refl
    else:
      self.reflections.extend(new_refl)

    # Store rotation matrix
    self.orientations[new_id] = rot
    self.reference_spots[new_id] = (i1, i2)

  def add_2d(self, expt, refl, i1, i2):
    """Add experiment with initial alignment based on detector coordinates

    Args:
      expt: experiment to add
      refl: reflection table
      i1, i2: indices of reference reflections (i2 should be longer vector)
    """
    # Get detector coordinates
    xyz = refl['xyzobs.mm.value']
    beam_x, beam_y = expt.detector[0].get_beam_centre(expt.beam.get_s0())

    # Get coordinates relative to beam center
    x = xyz.parts()[0] - beam_x
    y = xyz.parts()[1] - beam_y

    # Get reference spot coordinates
    x1, y1 = x[i1], y[i1]
    x2, y2 = x[i2], y[i2]

    # Align second spot (i2) to 3:00 position
    phi = -np.arctan2(y2, x2)

    # Apply rotation to all spots
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    x_aligned = x * cos_phi - y * sin_phi
    y_aligned = x * sin_phi + y * cos_phi

    # Store aligned coordinates
    aligned_xyz = flex.vec3_double(
      x_aligned, y_aligned, xyz.parts()[2]
    )
    refl['xyzobs.aligned'] = aligned_xyz

    # Apply same rotation to s1 vectors
    s1_vectors = refl['s1'].as_numpy_array()
    # Create rotation matrix around beam direction
    s0 = expt.beam.get_s0()
    beam_axis = s0/np.linalg.norm(s0)
    K = np.array([[0, -beam_axis[2], beam_axis[1]],
                  [beam_axis[2], 0, -beam_axis[0]],
                  [-beam_axis[1], beam_axis[0], 0]])
    rot = np.eye(3) + np.sin(phi)*K + (1-np.cos(phi))*K.dot(K)

    # Apply rotation to s1 vectors
    s1_aligned = flex.vec3_double([rot.dot(s1) for s1 in s1_vectors])
    refl['s1.aligned'] = s1_aligned

    # Determine and store handedness
    hand = flex.int(len(refl), int(y_aligned[i1] > 0))
    refl['hand'] = hand

    # Store experiment and reflections
    new_id = len(self.experiments)
    self.experiments.append(expt)
    refl['id'] = flex.int(len(refl), new_id)

    if len(self.reflections) == 0:
      self.reflections = refl
    else:
      self.reflections.extend(refl)

    # Store reference spots
    self.reference_spots[new_id] = (i1, i2)

  def plot_2d_radial_disabled(self):
    """Plot aligned detector coordinates in radial format.

    Transforms (x,y) detector coordinates into (θ,q) coordinates where:
    - θ is the azimuthal angle from the y-axis (horizontal axis)
    - q is the distance from beam center (vertical axis)
    - q axis is labeled in d-spacing units
    """
    fig, axes = plt.subplots(1,2, figsize=(16,8))

    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_refs = self.reflections.select(exp_sel)
      xyz = exp_refs['xyzobs.aligned']
      x, y = xyz.parts()[0], xyz.parts()[1]
      hand = exp_refs['hand'][0]  # All spots in experiment have same hand
      ax = axes[hand]

      # Calculate q vectors
      #q_vectors = exp_refs['s1'].as_numpy_array() - self.experiments[i_exp].beam.get_s0()
      q_vectors = exp_refs['q.aligned'].as_numpy_array()
      q_mags = np.sqrt(np.sum(q_vectors**2, axis=1))

      # Calculate theta with 0 at y-axis, increasing clockwise
      theta = np.rad2deg(np.arctan2(x, y))  # arctan2(x,y) gives angle from y-axis

      # Plot spots with different colors for different hands
      color = 'blue' if hand else 'red'
      ax.scatter(theta, q_mags, alpha=0.5, s=1, c=color)

      # Plot reference spots
      i1, i2 = self.reference_spots[i_exp]
      ref_theta = np.rad2deg(np.arctan2([x[i1], x[i2]], [y[i1], y[i2]]))
      ref_q = q_mags[[i1, i2]]
      ax.scatter(ref_theta, ref_q, c=['red', 'green'], s=50)

    # Configure axes
    for ax in axes:
      # Set theta limits from -180 to +180
      ax.set_xlim(-180, 180)

      # Format q-axis as d-spacing
      ax.yaxis.set_major_formatter(tick.FuncFormatter(
        lambda q, _: f"{1/q:.3f}" if q > 0 else ""
      ))

      # Add grid
      ax.grid(True, linestyle=':')

      # Labels
      ax.set_xlabel('Azimuthal angle θ (degrees)')
      ax.set_ylabel('Resolution d (Å)')

    # Titles
    axes[0].set_title('Left-handed')
    axes[1].set_title('Right-handed')

    plt.tight_layout()
    plt.show()

  def update_reflection_indices(self):
    """Update the i_refl column to match current table positions."""
    self.reflections['i_refl'] = flex.size_t(range(len(self.reflections)))

  def remove_poorly_sampled_experiments(self):
    """Remove experiments that have reflections in fewer than 3 selected groups."""
    # Count groups per experiment
    group_counts = {}
    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_refs = self.reflections.select(exp_sel)
      groups = set(exp_refs['refinement_group']) - {0}  # Exclude unselected
      group_counts[i_exp] = len(groups)

    # Identify experiments to remove
    to_remove = [i_exp for i_exp, count in group_counts.items() if count < 3]
    if not to_remove:
      return

    print(f"Removing {len(to_remove)} experiments with < 3 groups")

    # Create new experiment list excluding removed experiments
    new_experiments = ExperimentList([
      expt for i, expt in enumerate(self.experiments)
      if i not in to_remove
    ])

    # Create lookup for new experiment indices
    old_to_new = {}
    new_id = 0
    for old_id in range(len(self.experiments)):
      if old_id not in to_remove:
        old_to_new[old_id] = new_id
        new_id += 1

    # Select reflections to keep and update their experiment IDs
    keep_sel = flex.bool(len(self.reflections), False)
    new_ids = flex.int(len(self.reflections), 0)
    for i in range(len(self.reflections)):
      old_id = self.reflections['id'][i]
      if old_id not in to_remove:
        keep_sel[i] = True
        new_ids[i] = old_to_new[old_id]

    # Update reflection table
    self.reflections = self.reflections.select(keep_sel)
    self.reflections['id'] = new_ids.select(keep_sel)

    # Update experiment list
    self.experiments = new_experiments

    # Update reference_spots dictionary
    new_reference_spots = {}
    for old_id, (i1, i2) in self.reference_spots.items():
      if old_id not in to_remove:
        new_reference_spots[old_to_new[old_id]] = (i1, i2)
    self.reference_spots = new_reference_spots

    # Update i_refl column to match new table positions
    self.update_reflection_indices()

  def theta_refine(self):
    """Refine experiment orientations by minimizing theta spread of selected groups.

    For each hand (left/right), simultaneously optimize rotations of all experiments
    to minimize the spread of theta values within reflection groups.
    """
    for hand in [0, 1]:
      # Get experiments for this hand
      expts_this_hand = []
      for i_exp in range(len(self.experiments)):
        exp_sel = self.reflections['id'] == i_exp
        if len(self.reflections.select(exp_sel)) > 0:  # Make sure experiment has reflections
          if self.reflections.select(exp_sel)['hand'][0] == hand:
            expts_this_hand.append(i_exp)

      if not expts_this_hand:
        continue

      print(f"\nRefining {'left' if hand==0 else 'right'}-handed experiments")

      def rotate_expt(i_exp, dphi):
        """Apply rotation to an experiment's aligned q-vectors"""
        exp_sel = self.reflections['id'] == i_exp
        q_vecs = self.reflections['q.aligned'].select(exp_sel)

        # Debug: print a vector before and after rotation
        if len(q_vecs) > 0:
          q0 = q_vecs[0]
          theta_before = np.rad2deg(np.arctan2(q0[1], q0[0]))

          # Create rotation matrix around z
          phi_rad = np.deg2rad(dphi)
          rot = np.array([[np.cos(phi_rad), -np.sin(phi_rad), 0],
                         [np.sin(phi_rad), np.cos(phi_rad), 0],
                         [0, 0, 1]])
          # Apply rotation
          q_new = flex.vec3_double([rot.dot(q) for q in q_vecs])

          q0_new = q_new[0]
          theta_after = np.rad2deg(np.arctan2(q0_new[1], q0_new[0]))

          self.reflections['q.aligned'].set_selected(exp_sel, q_new)

      def get_experiment_stats(i_exp):
        """Compute alignment statistics for an experiment"""
        exp_sel = self.reflections['id'] == i_exp
        exp_refs = self.reflections.select(exp_sel)
        q_vecs = exp_refs['q.aligned']
        # Get angles in xy plane
        exp_thetas = np.rad2deg(np.arctan2(q_vecs.parts()[1], q_vecs.parts()[0]))

        stats = {'n_refs': 0, 'rmsd': 0}
        # For each group this experiment contributes to
        for group in set(exp_refs['refinement_group']) - {0}:
          # Get all thetas for this group (from all experiments)
          group_sel = self.reflections['refinement_group'] == group
          group_refs = self.reflections.select(group_sel)
          group_q = group_refs['q.aligned']
          group_thetas = np.rad2deg(np.arctan2(
            group_q.parts()[1], group_q.parts()[0]))

          # Map all angles to [-180,180] before computing mean
          group_thetas = (group_thetas + 180) % 360 - 180
          group_mean = np.mean(group_thetas)

          # Get this experiment's contributions
          exp_group_sel = exp_refs['refinement_group'] == group
          exp_group_thetas = exp_thetas[exp_group_sel]
          exp_group_thetas = (exp_group_thetas + 180) % 360 - 180

          # Add to statistics
          diffs = exp_group_thetas - group_mean
          diffs = (diffs + 180) % 360 - 180  # Map differences to [-180,180]
          stats['n_refs'] += len(exp_group_thetas)
          stats['rmsd'] += np.sum(diffs**2)

        if stats['n_refs'] > 0:
          stats['rmsd'] = np.sqrt(stats['rmsd'] / stats['n_refs'])
        return stats

      def target_and_gradient(params):
        """Compute target function and gradient for current state"""
        # Apply rotations for this set of parameters
        for i_exp, dphi in zip(expts_this_hand, params):
          rotate_expt(i_exp, dphi)

        # Get groups for this hand
        sign = -1 if hand == 0 else 1
        groups = [g for g in set(self.reflections['refinement_group'])
                 if g * sign > 0]  # Filter for this hand's groups

        residual = 0
        gradient = np.zeros_like(params)

        # For each group
        for group in groups:
          group_sel = self.reflections['refinement_group'] == group
          group_refs = self.reflections.select(group_sel)

          # Get theta values and print some diagnostics
          q_vecs = group_refs['q.aligned']
          thetas_raw = np.rad2deg(np.arctan2(q_vecs.parts()[1], q_vecs.parts()[0]))
          #print(f"\nGroup {group}:")
          #print(f"Raw theta range: {min(thetas_raw):.1f}° to {max(thetas_raw):.1f}°")

          thetas = (thetas_raw + 180) % 360 - 180  # Map to [-180,180]
          #print(f"Mapped theta range: {min(thetas):.1f}° to {max(thetas):.1f}°")

          # Also look at some example vectors
          #print("Example q-vectors and their angles:")
          for i in range(min(5, len(q_vecs))):
            q = q_vecs[i]
            theta_raw = np.rad2deg(np.arctan2(q[1], q[0]))
            theta_mapped = (theta_raw + 180) % 360 - 180
            #print(f"q = [{q[0]:.3f}, {q[1]:.3f}, {q[2]:.3f}], "
                  #f"theta = {theta_raw:.1f}° -> {theta_mapped:.1f}°")


          # For each reflection in group
          for i, theta in enumerate(thetas):
            # Get experiment index in params array
            i_exp = group_refs['id'][i]
            param_idx = expts_this_hand.index(i_exp)

            # Compute differences with all other reflections in group
            diffs = thetas - theta
            diffs = (diffs + 180) % 360 - 180  # Map differences to [-180,180]

            # Add to residual (sum of absolute differences)
            residual += np.sum(np.abs(diffs))

            # For gradient, sum the differences (larger angles contribute more)
            gradient[param_idx] += np.sum(diffs) / len(thetas)

        # Undo the rotations for the next iteration
        for i_exp, dphi in zip(expts_this_hand, params):
          rotate_expt(i_exp, -dphi)

        return residual, gradient

      # Initial parameters (zero rotation for each experiment)
      x0 = np.zeros(len(expts_this_hand))

      # Optimize
      from scipy.optimize import minimize
      result = minimize(
        lambda x: target_and_gradient(x)[0],
        x0,
        jac=lambda x: target_and_gradient(x)[1],
        method='L-BFGS-B',
        options={'maxiter': 100, 'disp': True}  # Add disp=True
      )

      if result.success:
        print(f"Optimization succeeded after {result.nit} iterations")
        print(f"Final residual: {result.fun:.1f}")
        # Apply final rotations and print statistics
        print("\nExperiment statistics:")
        for i_exp, phi in zip(expts_this_hand, result.x):
          rotate_expt(i_exp, phi)
          stats = get_experiment_stats(i_exp)
          print(f"  Experiment {i_exp:3d}: rotation = {phi:6.2f}°, "
                f"RMSD = {stats['rmsd']:5.2f}° from {stats['n_refs']:3d} reflections")
      else:
        print(f"Optimization failed: {result.message}")




  def plot_2d_radial(self):
    """Plot aligned detector coordinates in radial format with selection capability."""

    # Initialize or reset refinement_group column
    self.reflections['refinement_group'] = flex.int(len(self.reflections), 0)  # Reset all to zero

    # Ensure we have reflection indices
    if 'i_refl' not in self.reflections:
      self.reflections['i_refl'] = flex.size_t(range(len(self.reflections)))

    fig, axes = plt.subplots(1,2, figsize=(16,8))

    # Store additional experiment info in plot_data
    plot_data = {
      0: {'theta': [], 'q': [], 'i_refl': [], 'exp_id': [], 'scatter_plots': {}},  # left-handed
      1: {'theta': [], 'q': [], 'i_refl': [], 'exp_id': [], 'scatter_plots': {}}   # right-handed
    }
    print(f"Total number of experiments: {len(self.experiments)}")

    # Plot points and collect data
    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_refs = self.reflections.select(exp_sel)
      q_vecs = exp_refs['q.aligned']
      x, y = q_vecs.parts()[0], q_vecs.parts()[1]  # Get x,y components
      hand = exp_refs['hand'][0]
      ax = axes[hand]

      # Calculate theta using arctan2(y,x) to match theta_refine
      theta = np.rad2deg(np.arctan2(y, x))
      q_mags = np.sqrt(q_vecs.parts()[0]**2 + q_vecs.parts()[1]**2 + q_vecs.parts()[2]**2)

      # Plot spots
      color = 'blue' if hand else 'red'
      scatter = ax.scatter(theta, q_mags, alpha=0.2, s=3, c=color)
      plot_data[hand]['scatter_plots'][i_exp] = scatter  # This line is crucial!

      # Plot reference spots
      i1, i2 = self.reference_spots[i_exp]
      ref_theta = np.rad2deg(np.arctan2([y[i1], y[i2]], [x[i1], x[i2]]))
      ref_q = q_mags[[i1, i2]]
      ax.scatter(ref_theta, ref_q, c=['red', 'green'], s=50)

      # Store data for selections
      plot_data[hand]['theta'].extend(theta)
      plot_data[hand]['q'].extend(q_mags)
      plot_data[hand]['i_refl'].extend(exp_refs['i_refl'])
      plot_data[hand]['exp_id'].extend([i_exp] * len(exp_refs))

    # Convert lists to arrays for easier indexing
    for hand in [0, 1]:
      plot_data[hand]['theta'] = np.array(plot_data[hand]['theta'])
      plot_data[hand]['q'] = np.array(plot_data[hand]['q'])
      plot_data[hand]['i_refl'] = np.array(plot_data[hand]['i_refl'])

    # Configure axes
    for ax in axes:
      ax.set_xlim(-180, 180)
      ax.set_ylim(.1, .5)
      ax.yaxis.set_major_formatter(tick.FuncFormatter(
        lambda q, _: f"{1/q:.3f}" if q > 0 else ""
      ))
      ax.grid(True, linestyle=':')
      ax.set_xlabel('Azimuthal angle θ (degrees)')
      ax.set_ylabel('Resolution d (Å)')

    axes[0].set_title('Left-handed')
    axes[1].set_title('Right-handed')

    # Selection handling
    next_group = [1, 1]  # [left_hand, right_hand] group counters
    selections = []

    def on_select(eclick, erelease):
      ax = eclick.inaxes
      if ax not in axes:
        return

      hand = 0 if ax == axes[0] else 1
      sign = -1 if hand == 0 else 1

      # Get selection bounds
      x1, y1 = eclick.xdata, eclick.ydata
      x2, y2 = erelease.xdata, erelease.ydata
      xmin, xmax = min(x1, x2), max(x1, x2)
      ymin, ymax = min(y1, y2), max(y1, y2)

      # Find points in selection
      data = plot_data[hand]
      mask = ((data['theta'] >= xmin) & (data['theta'] <= xmax) &
              (data['q'] >= ymin) & (data['q'] <= ymax))

      if np.any(mask):
        # Assign group number to selected reflections
        group_num = sign * next_group[hand]
        selected_indices = np.array(data['i_refl'])[mask]
        self.reflections['refinement_group'].set_selected(
          flex.size_t(selected_indices),
          flex.int(len(selected_indices), group_num)
        )

#        # Color all points from experiments that have selected reflections
        selected_exp_ids = set(np.array(data['exp_id'])[mask])
#        for exp_id in selected_exp_ids:
#          plot_data[hand]['scatter_plots'][exp_id].set_color('red')

        # Highlight selection
        rect = plt.Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,
                            fill=False, color='green')
        ax.add_patch(rect)
        selections.append((rect, selected_exp_ids, hand))  # Store exp_ids with selection

        next_group[hand] += 1
        fig.canvas.draw_idle()

    def on_key(event):
      if event.key == 'u' and selections:  # Undo last selection
        rect, exp_ids, hand = selections.pop()
        rect.remove()
        # Reset colors for experiments if they're not in any remaining selections
        remaining_exp_ids = set()
        for _, exp_ids_sel, hand_sel in selections:
          if hand_sel == hand:  # Only consider selections in same hand
            remaining_exp_ids.update(exp_ids_sel)
        for exp_id in exp_ids:
          if exp_id not in remaining_exp_ids:
            plot_data[hand]['scatter_plots'][exp_id].set_color('black')
        fig.canvas.draw_idle()
      elif event.key == 'r':  # Reset all selections
        self.reflections['refinement_group'] = flex.int(len(self.reflections), 0)
        for rect, _, _ in selections:
          rect.remove()
        selections.clear()
        # Reset all colors to black
        for hand in [0, 1]:
          for scatter in plot_data[hand]['scatter_plots'].values():
            scatter.set_color('black')
        next_group[0] = next_group[1] = 1
        fig.canvas.draw_idle()

    # Connect event handlers
    rect_selector = RectangleSelector(
      axes[0], on_select, interactive=False,
      useblit=True,
      props=dict(facecolor='none', edgecolor='green', alpha=0.5)
    )
    rect_selector2 = RectangleSelector(
      axes[1], on_select, interactive=False,
      useblit=True,
      props=dict(facecolor='none', edgecolor='green', alpha=0.5)
    )
    fig.canvas.mpl_connect('key_press_event', on_key)

    plt.tight_layout()
    print("\nSelection controls:")
    print("- Click and drag to select reflections")
    print("- Press 'u' to undo last selection")
    print("- Press 'r' to reset all selections")
    plt.show()

  def select_group_dspacings(self):
    """For each refinement group, plot d-spacing histogram and collect user selections."""
    group_nums = sorted(set(self.reflections['refinement_group']) - {0})
    group_info = {}

    for group_num in group_nums:
      # Get reflections for this group
      group_sel = self.reflections['refinement_group'] == group_num
      group_refs = self.reflections.select(group_sel)

      # Get q-vectors and compute d-spacings
      q_vectors = group_refs['q.aligned'].as_numpy_array()
      q_mags = np.sqrt(np.sum(q_vectors**2, axis=1))
      d_spacings = 1/q_mags

      # Print diagnostic information
      print(f"\nGroup {group_num}:")
      print(f"Number of reflections: {len(group_refs)}")
      print(f"d-spacing range: {min(d_spacings):.2f} - {max(d_spacings):.2f} Å")
      print(f"Contributing experiments: {sorted(set(group_refs['id']))}")

      # Check for outliers (more than 3 sigma from mean)
      mean_d = np.mean(d_spacings)
      std_d = np.std(d_spacings)
      std_q = np.std(q_mags)
      outliers = np.abs(d_spacings - mean_d) > 3*std_d
      if np.any(outliers):
        print(f"WARNING: {np.sum(outliers)} outlier d-spacings:")
        print(f"Outlier values: {d_spacings[outliers]}")

      # Create histogram
      fig, ax = plt.subplots()
      hist, edges, _ = ax.hist(d_spacings, bins=100)
      ax.set_xlabel('d-spacing (Å)')
      ax.set_ylabel('Count')
      ax.set_title(f'Group {group_num}: Click on peak or press q to skip')

      # Add mean and std lines
      ax.axvline(mean_d, color='r', linestyle='--', label=f'Mean: {mean_d:.2f} Å')
      ax.axvline(mean_d + std_d, color='g', linestyle=':', label='+1σ')
      ax.axvline(mean_d - std_d, color='g', linestyle=':', label='-1σ')
      ax.legend()

      # Store data for callback
      click_data = {'d_selected': None}

      def on_click(event):
        if event.inaxes == ax:
          click_data['d_selected'] = event.xdata
          plt.close(fig)

      def on_key(event):
        if event.key == 'q':
          plt.close(fig)

      # Connect event handlers
      fig.canvas.mpl_connect('button_press_event', on_click)
      fig.canvas.mpl_connect('key_press_event', on_key)

      plt.show()

      if click_data['d_selected'] is not None:
        # Store group information
        q_mean = np.mean(q_vectors, axis=0)
        group_info[group_num] = {
          'd_spacing': click_data['d_selected'],
          'q_vector': q_mean,  # Store mean q-vector for later use
          'hand': -1 if group_num < 0 else 1,
          'std_q': std_q
        }

    self.group_info = group_info
    print(f"\nCollected d-spacings for {len(group_info)} groups")

  def analyze_group_angles(self, triplets):
    """Analyze angles between pairs of groups using full dataset."""
    if not hasattr(self, 'group_info'):
      raise Sorry("No group information available. Run select_group_dspacings first.")

    groups = sorted(self.group_info.keys())
    n_groups = len(groups)
    dot_products = np.zeros((n_groups, n_groups))

    # Process all pairs
    for i, group1 in enumerate(groups):
      for j, group2 in enumerate(groups[:i+1]):
        if i == j:
          continue

        info1 = self.group_info[group1]
        info2 = self.group_info[group2]

        # Get q magnitudes for search and ensure d1 > d2 order
        q1 = 1/info1['d_spacing']
        q2 = 1/info2['d_spacing']
        if q1 > q2:  # d1 < d2, need to swap
          q1, q2 = q2, q1
          info1, info2 = info2, info1

        # Get observed angle between groups
        v1 = info1['q_vector']
        v2 = info2['q_vector']
        obs_angle = np.rad2deg(np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))

        print(f"\nAnalyzing groups {group1} and {group2}")
        print(f"d-spacings: {info1['d_spacing']:.2f}Å, {info2['d_spacing']:.2f}Å")
        print(f"Observed angle: {obs_angle:.1f}°")

        # Find matching pairs in full dataset using stored standard deviations
        q1_matches = np.abs(1/triplets.triplets[:,1] - q1) < info1['std_q']
        q2_matches = np.abs(1/triplets.triplets[:,2] - q2) < info2['std_q']
        pair_matches = q1_matches & q2_matches

        if not np.any(pair_matches):
          print("No matching pairs found in dataset")
          continue

        # Get angles for matching pairs and create mirrored set
        angles = triplets.triplets[pair_matches, 3]
        mirrored_angles = 180 - angles  # Mirror across 90°
        all_angles = np.concatenate([angles, mirrored_angles])
        print(f"Found {len(angles)} matching pairs ({len(all_angles)} with mirrors)")

        # Plot histogram with observed angle marked
        fig, ax = plt.subplots()
        hist, edges, _ = ax.hist(all_angles, bins=180, range=(0, 180))
        ax.axvline(obs_angle, color='r', linestyle='--',
                  label=f'Observed: {obs_angle:.1f}°')

        # If unit cell is provided, compute theoretical angles
        if self.params.debug.unit_cell is not None:
          uc = self.params.debug.unit_cell
          # Get d-spacing ranges
          d1_range = (1/(q1 + 1*info1['std_q']), 1/(q1 - 1*info1['std_q']))
          d2_range = (1/(q2 + 1*info2['std_q']), 1/(q2 - 1*info2['std_q']))

          # Generate Miller indices within these ranges
          from cctbx import miller
          cs = symmetry(unit_cell=uc, space_group=sgtbx.space_group('P1'))
          ms = miller.build_set(
            crystal_symmetry=cs,
            anomalous_flag=True,
            d_min=min(d1_range[0], d2_range[0])
          )
          indices = ms.indices()
          d_spacings = uc.d(indices)

          # Find pairs within ranges
          mask1 = (d_spacings >= d1_range[0]) & (d_spacings <= d1_range[1])
          mask2 = (d_spacings >= d2_range[0]) & (d_spacings <= d2_range[1])
          indices1 = indices.select(mask1)
          indices2 = indices.select(mask2)

          # Compute angles between all pairs
          theory_angles = []
          for h1 in indices1:
            for h2 in indices2:
              # Skip if same reflection
              if h1 == h2:
                continue
              v1 = uc.reciprocal_space_vector(h1)
              v2 = uc.reciprocal_space_vector(h2)
              try:
                angle = math.degrees(math.acos(np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)))
                #angle = uc.angle(h1, h2)
                theory_angles.append(float(angle))
                theory_angles.append(180 - float(angle))  # Add mirror
              except Exception:
                print('math domain error for ', h1, h2)
                continue


          # Plot theoretical angles
          ymax = ax.get_ylim()[1]
          for angle in theory_angles:
            ax.plot([angle, angle], [-ymax*.1, 0], 'r-', linewidth=1)

        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('Count')
        ax.set_title(f'Groups {group1} and {group2}\nClick on peak or press q to skip')
        ax.legend()

        # Store data for callback
        click_data = {'angle_selected': None}

        def on_click(event):
          if event.inaxes == ax:
            click_data['angle_selected'] = event.xdata
            plt.close(fig)

        def on_key(event):
          if event.key == 'q':
            plt.close(fig)

        # Connect event handlers
        fig.canvas.mpl_connect('button_press_event', on_click)
        fig.canvas.mpl_connect('key_press_event', on_key)

        plt.show()


  def refine_group_vectors(self):
    """Refine positions of group vectors using dot product data."""
    if not hasattr(self, 'dot_products'):
      raise Sorry("No dot product data. Run analyze_group_angles first.")

    groups = self.dot_product_groups
    n_groups = len(groups)

    # Initialize parameters: 2 angles per vector
    x0 = []
    q_mags = []  # Store magnitudes for conversion
    for group in groups:
      v = self.group_info[group]['q_vector']
      q_mag = np.linalg.norm(v)
      q_mags.append(q_mag)
      # Convert to spherical coordinates
      theta = np.arccos(v[2]/q_mag)
      phi = np.arctan2(v[1], v[0])
      x0.extend([theta, phi])

    def target(params):
      """Compute residual from current parameters"""
      # Convert parameters to vectors
      vectors = []
      for i in range(n_groups):
        theta, phi = params[2*i:2*i+2]
        sin_theta = np.sin(theta)
        v = q_mags[i] * np.array([
          sin_theta * np.cos(phi),
          sin_theta * np.sin(phi),
          np.cos(theta)
        ])
        vectors.append(v)

      # Compute residual from dot products
      residual = 0
      for i in range(n_groups):
        for j in range(i+1):
          target_dot = self.dot_products[i,j]
          if target_dot != 0:  # Skip unobserved pairs
            calc_dot = np.dot(vectors[i], vectors[j])
            residual += (calc_dot - target_dot)**2

      return residual

    # Minimize
    result = minimize(
      target, x0,
      method='L-BFGS-B',
      options={'maxiter': 100, 'disp': True}
    )

    if result.success:
      print(f"Refinement succeeded after {result.nit} iterations")
      print(f"Final residual: {result.fun:.6f}")

      # Store refined vectors
      refined_vectors = {}
      for i, group in enumerate(groups):
        theta, phi = result.x[2*i:2*i+2]
        sin_theta = np.sin(theta)
        v = q_mags[i] * np.array([
          sin_theta * np.cos(phi),
          sin_theta * np.sin(phi),
          np.cos(theta)
        ])
        refined_vectors[group] = v

      self.refined_vectors = refined_vectors

      # Print angle matrices
      print("\nInitial angles:")
      for i, group1 in enumerate(groups):
        angles = []
        v1 = self.group_info[group1]['q_vector']
        for group2 in groups:
          v2 = self.group_info[group2]['q_vector']
          angle = np.rad2deg(np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
          angles.append(f"{angle:6.1f}")
        print(f"Group {group1:3d}: {' '.join(angles)}")

      print("\nTarget angles (from dot products):")
      for i, group1 in enumerate(groups):
        angles = []
        for j, group2 in enumerate(groups):
          if i == j:
            angle = 0
          else:
            dot = self.dot_products[i,j]
            q1 = 1/self.group_info[group1]['d_spacing']
            q2 = 1/self.group_info[group2]['d_spacing']
            angle = np.rad2deg(np.arccos(dot/(q1*q2)))
          angles.append(f"{angle:6.1f}")
        print(f"Group {group1:3d}: {' '.join(angles)}")

      print("\nFinal angles:")
      for i, group1 in enumerate(groups):
        angles = []
        v1 = refined_vectors[group1]
        for group2 in groups:
          v2 = refined_vectors[group2]
          angle = np.rad2deg(np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))))
          angles.append(f"{angle:6.1f}")
        print(f"Group {group1:3d}: {' '.join(angles)}")
    else:
      print(f"Refinement failed: {result.message}")

  def plot_vectors(self):
    """3D scatter plot comparing original and refined q-vectors with full dataset."""

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot all aligned vectors from dataset
#    s1_vectors = self.reflections['s1.aligned'].as_numpy_array()
#    s0 = np.array([expt.beam.get_s0() for expt in self.experiments])
#    exp_ids = self.reflections['id']
#    q_vectors = s1_vectors - s0[exp_ids]
    q_vectors = self.reflections['q.aligned'].as_numpy_array()

    ax.scatter(q_vectors[:,0], q_vectors[:,1], q_vectors[:,2],
              c='gray', alpha=0.1, s=1, label='Dataset')

    # Plot original group vectors
    if hasattr(self, 'dot_product_groups'):
      orig_x, orig_y, orig_z = [], [], []
      for group in self.dot_product_groups:
        v = self.group_info[group]['q_vector']
        orig_x.append(v[0])
        orig_y.append(v[1])
        orig_z.append(v[2])
      ax.scatter(orig_x, orig_y, orig_z, c='blue', s=100, label='Original')

    if hasattr(self, 'refined_vectors'):
      # Plot refined vectors
      ref_x, ref_y, ref_z = [], [], []
      for group in self.dot_product_groups:
        v = self.refined_vectors[group]
        ref_x.append(v[0])
        ref_y.append(v[1])
        ref_z.append(v[2])
      ax.scatter(ref_x, ref_y, ref_z, c='red', s=100, label='Refined')

      # Draw lines connecting original to refined positions
      for i in range(len(self.dot_product_groups)):
        ax.plot([orig_x[i], ref_x[i]],
                [orig_y[i], ref_y[i]],
                [orig_z[i], ref_z[i]], 'k-', alpha=0.3)

    ax.set_xlabel('qx')
    ax.set_ylabel('qy')
    ax.set_zlabel('qz')
    ax.legend()
    plt.show()



  def plot_2d(self):
    """Plot all aligned detector coordinates"""
    fig, axes = plt.subplots(1,2)
    fig.set_figwidth(16)
    fig.set_figheight(8)

    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_refs = self.reflections.select(exp_sel)
      xyz = exp_refs['q.aligned']
      x, y = xyz.parts()[0], xyz.parts()[1]
      hand = exp_refs['hand'][0]  # All spots in experiment have same hand
      ax = axes[hand]

      # Plot spots with different colors for different hands
      color = 'blue' if hand else 'red'
      ax.scatter(x, y, alpha=0.5, s=1, c=color)

      # Plot reference spots
      i1, i2 = self.reference_spots[i_exp]
      ax.scatter([x[i1], x[i2]], [y[i1], y[i2]],
                  c=['red', 'green'], s=50)

    plt.axhline(y=0, color='k', linestyle=':')
    plt.axvline(x=0, color='k', linestyle=':')
    plt.axis('equal')
    plt.title('Aligned detector coordinates')
#    for ax in axes:
#      ax.set_xlim((-75,75))
#      ax.set_ylim((-75,75))
    plt.show()
  def refine_orientations(self, max_angle_rad=np.pi/36):
    """Refine orientation matrices to improve cluster agreement"""
    if self.clusters is None:
      self.update_clusters()

    xyz = self.reflections['transformed_xyz'].as_numpy_array()
    labels = self.clusters['labels']

    from scipy.stats import chi2

    # Compute covariance matrices and store with centers
    cluster_stats = {}
    mahalanobis_cutoff = chi2.ppf(0.99, df=3)  # 99th percentile for 3D

    for k in set(labels) - {-1}:
      mask = labels == k
      cluster_points = xyz[mask]
      center = np.mean(cluster_points, axis=0)
      cov = np.cov(cluster_points, rowvar=False)
      cluster_stats[k] = {'center': center, 'cov': cov, 'cov_inv': np.linalg.inv(cov)}

    def small_rotations_to_matrix(rotx, roty, rotz):
      """Convert small rotation angles to matrix"""
      Rx = np.array([[1, 0, 0],
                     [0, np.cos(rotx), -np.sin(rotx)],
                     [0, np.sin(rotx), np.cos(rotx)]])
      Ry = np.array([[np.cos(roty), 0, np.sin(roty)],
                     [0, 1, 0],
                     [-np.sin(roty), 0, np.cos(roty)]])
      Rz = np.array([[np.cos(rotz), -np.sin(rotz), 0],
                     [np.sin(rotz), np.cos(rotz), 0],
                     [0, 0, 1]])
      return Rz @ Ry @ Rx

    def residual_for_experiment(i_exp, rot_params):
      """Compute Euclidean residuals for points within Mahalanobis cutoff"""
      exp_sel = self.reflections['id'] == i_exp
      exp_xyz = xyz[exp_sel]

      # First find which points and centers to use (without any rotation)
      selected_points = []
      target_centers = []
      for pt in exp_xyz:
        best_dist = np.inf
        best_center = None
        for k, stats in cluster_stats.items():
          diff = pt - stats['center']
          maha_dist = np.sqrt(diff @ stats['cov_inv'] @ diff)
          if maha_dist < best_dist:
            best_dist = maha_dist
            best_center = stats['center']

        if best_dist < mahalanobis_cutoff:
          selected_points.append(pt)
          target_centers.append(best_center)

      if not selected_points:
        return np.array([0])  # Return dummy residual if no points selected

      # Now compute residuals for selected points
      selected_points = np.array(selected_points)
      target_centers = np.array(target_centers)

      # Apply rotation adjustment
      rot_adjust = small_rotations_to_matrix(*rot_params)
      adjusted_xyz = np.array([rot_adjust @ v for v in selected_points])

      # Return cartesian residuals
      return (adjusted_xyz - target_centers).ravel()

    from scipy.optimize import least_squares

    # Refine each experiment
    refined_orientations = {}
    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_xyz = xyz[exp_sel]

#      # Create diagnostic plot
#      fig = plt.figure()
#      ax = fig.add_subplot(111, projection='3d')
#      ax.scatter(0, 0, 0, c='black', s=100, label='Origin')
#
#      # Plot cluster centers
#      for k, stats in cluster_stats.items():
#        ax.scatter(*stats['center'], c='red', s=100, marker='*', label=f'Center {k}')
#
#      # Plot initial positions
#      ax.scatter(*exp_xyz.T, alpha=0.5, label='Initial')
#
#      ax.set_title(f'Experiment {i_exp} - Before refinement')
#      ax.legend()

      result = least_squares(
        lambda x: residual_for_experiment(i_exp, x),
        x0=[0, 0, 0],
        method='lm'
      )

      if result.success:
        adjustment = small_rotations_to_matrix(*result.x)
        final_xyz = np.array([adjustment @ v for v in exp_xyz])

#        # Show final positions
#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
#        ax.scatter(0, 0, 0, c='black', s=100, label='Origin')
#
#        # Plot centers again
#        for k, stats in cluster_stats.items():
#          ax.scatter(*stats['center'], c='red', s=100, marker='*', label=f'Center {k}')
#
#        # Plot final positions
#        ax.scatter(*final_xyz.T, alpha=0.5, label='Final')
#
#        ax.set_title(f'Experiment {i_exp} - After refinement')
#        ax.legend()

        #print(f"  Final rms = {np.sqrt(np.mean(result.fun**2)):.6f}")
        refined_orientations[i_exp] = adjustment @ self.orientations[i_exp]
      else:
        print(f"  Refinement failed: {result.message}")
        refined_orientations[i_exp] = self.orientations[i_exp]

    # Update orientations and recompute transformed coordinates
    self.orientations = refined_orientations

    # Update transformed coordinates
    new_xyz = []
    for i_exp in range(len(self.experiments)):
      exp_sel = self.reflections['id'] == i_exp
      exp_refl = self.reflections.select(exp_sel)
      q_vectors = exp_refl['s1'].as_numpy_array() - self.experiments[i_exp].beam.get_s0()
      transformed = np.array([self.orientations[i_exp] @ q for q in q_vectors])
      new_xyz.append(transformed)

    self.reflections['transformed_xyz'] = flex.vec3_double(np.vstack(new_xyz))
    self.clusters = None  # Force recomputation of clusters

  def plot_polar(self):
    """Plot reconstructed spots in polar coordinates"""
    # Get reciprocal space vectors
    xyz = self.reflections['transformed_xyz'].as_numpy_array()

    # Convert to q, phi, theta
    q = np.sqrt(np.sum(xyz**2, axis=1))
    phi = np.arctan2(xyz[:,1], xyz[:,0])  # -pi to pi
    theta = np.arccos(xyz[:,2]/q)         # 0 to pi

    # First make powder pattern
    fig1 = plt.figure(figsize=(12,4))
    ax1 = plt.gca()

    hist, edges = np.histogram(q, bins=2000)
    ax1.plot(edges[:-1], hist)

    # Format x-axis as d-spacing
    ax1.get_xaxis().set_major_formatter(tick.FuncFormatter(
      lambda x, _: "{:.3f}".format(1/x)))
    ax1.set_xlabel('d-spacing (Å)')
    ax1.set_ylabel('Counts')
    ax1.set_title('Select up to 9 ranges, press q when done')

    # Allow interactive range selection
    spans = []
    q_ranges = []
    selection_done = False

    def on_select(qmin, qmax):
      if len(q_ranges) >= 9:
        return
      q_ranges.append(tuple(sorted([qmin, qmax])))
      span = ax1.axvspan(qmin, qmax, color='red', alpha=0.2)
      spans.append(span)
      fig1.canvas.draw_idle()

    def on_key(event):
      nonlocal selection_done
      if event.key == 'q':
        selection_done = True
        plt.close(fig1)

    selector = SpanSelector(
      ax1, on_select, 'horizontal',
      props=dict(alpha=0.3, facecolor='red'),
      interactive=False,
      drag_from_anywhere=True
    )

    fig1.canvas.mpl_connect('key_press_event', on_key)
    plt.show()

    if not q_ranges:
      return

    # Create grid of angular plots
    fig2 = plt.figure(figsize=(15,15))
    for i, (qmin, qmax) in enumerate(q_ranges[:9]):
      ax = fig2.add_subplot(3, 3, i+1)

      mask = (q >= qmin) & (q <= qmax)
      phi_deg = np.rad2deg(phi[mask])
      theta_deg = np.rad2deg(theta[mask])

      ax.scatter(phi_deg, theta_deg, alpha=0.5, s=1)
      ax.set_xlabel('φ (degrees)')
      ax.set_ylabel('θ (degrees)')
      ax.set_title(f'd: {1/qmax:.2f}-{1/qmin:.2f}Å')

    plt.tight_layout()
    plt.show()

  def plot_debug(self, title=None):
    """Simple 3D scatter plot of all transformed reflections"""
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xyz = self.reflections['transformed_xyz'].as_numpy_array()

    if self.clusters is None:
      # Plot all reflections in blue
      ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2],
                c='blue', alpha=0.5, s=1, label='All spots')

      # Plot reference reflections in red and green
      ref_spots = []
      for i_exp in range(len(self.experiments)):
        exp_sel = self.reflections['id'] == i_exp
        exp_refs = self.reflections.select(exp_sel)
        i1, i2 = self.reference_spots[i_exp]
        ref_spots.append(exp_refs['transformed_xyz'][i1])
        ref_spots.append(exp_refs['transformed_xyz'][i2])

      ref_spots = np.array(ref_spots)
      ax.scatter(ref_spots[::2,0], ref_spots[::2,1], ref_spots[::2,2],
                c='red', s=10, label='Reference A')
      ax.scatter(ref_spots[1::2,0], ref_spots[1::2,1], ref_spots[1::2,2],
                c='green', s=10, label='Reference B')

    else:
      # Print cluster information
      labels = self.clusters['labels']
      unique_labels = set(labels)
      print("\nCluster information:")

      # Handle noise points
      noise = labels == -1
      if np.any(noise):
        noise_points = xyz[noise]
        print(f"Noise points: {len(noise_points)} reflections")

      # Print info for each cluster
      for k in sorted(unique_labels):
        if k == -1:
          continue
        mask = labels == k
        cluster_xyz = xyz[mask]
        center = np.mean(cluster_xyz, axis=0)
#        print(f"\nCluster {k}:")
#        print(f"  Points: {len(cluster_xyz)}")
#        print(f"  Center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")

      # Analyze experiment contributions to clusters
      exp_cluster_counts = {}  # Dict mapping experiment id to number of clusters it hits
      for i_exp in range(len(self.experiments)):
        exp_sel = self.reflections['id'] == i_exp
        exp_labels = set(labels[exp_sel])  # unique clusters this experiment contributes to
        if -1 in exp_labels:  # don't count noise as a cluster
          exp_labels.remove(-1)
        exp_cluster_counts[i_exp] = len(exp_labels)

      # Summarize the results
      contribution_summary = {}  # Dict mapping number of clusters to count of experiments
      for n_clusters in exp_cluster_counts.values():
        contribution_summary[n_clusters] = contribution_summary.get(n_clusters, 0) + 1

      print("\nExperiments contributing to clusters:")
      for n_clusters in sorted(contribution_summary.keys()):
        print(f"{n_clusters} clusters: {contribution_summary[n_clusters]} experiments")



      # Color by cluster
      labels = self.clusters['labels']
      # Get a color for each cluster (-1 is black for noise)
      unique_labels = set(labels)
      colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_labels)-1))

      # Plot noise points first in black
      noise = labels == -1
      if np.any(noise):
        ax.scatter(xyz[noise,0], xyz[noise,1], xyz[noise,2],
                  c='black', alpha=0.1, s=1, label='Noise')

      # Plot each cluster
      color_idx = 0
      for k in unique_labels:
        if k == -1:
          continue
        cluster_mask = labels == k
        ax.scatter(xyz[cluster_mask,0], xyz[cluster_mask,1], xyz[cluster_mask,2],
                  c=[colors[color_idx]], alpha=0.5, s=1,
                  label=f'Cluster {k}')
        color_idx += 1

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    if title: 
      ax.set_title(title)
    return fig
    
  def update_clusters(self, eps=.01, min_samples=5):
    """Perform cluster analysis on current reflection set

    Args:
      eps: Maximum distance between points in a cluster
      min_samples: Minimum points to form a cluster
    """

    if self.clusters is not None:
      return

    xyz = self.reflections['transformed_xyz'].as_numpy_array()
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(xyz)

    self.clusters = {
      'labels': clustering.labels_,  # -1 means noise
      'n_clusters': len(set(clustering.labels_)) - (1 if -1 in clustering.labels_ else 0)
    }
    
  def is_similar_triangle(self, d1, d2, theta, tolerance=0.1):
    """Check if a triangle is too similar to one already used
    
    Args:
      d1, d2: d-spacings
      theta: angle between spots
      tolerance: minimum distance in d1,d2,theta space
      
    Returns:
      bool: True if this triangle is too close to a used one
    """
    for used_d1, used_d2, used_theta in self.used_triangles:
      dist = np.sqrt((d1-used_d1)**2 + (d2-used_d2)**2 + (theta-used_theta)**2)
      if dist < tolerance:
        return True
    return False
    
  def find_candidate_triangles(self):
    """Find potential new triangles from current reconstruction"""
    self.update_clusters()
    # Use cluster analysis to find potential triangles
    # Filter using is_similar_triangle

  def match(self, other):
    """Check if two reconstructions represent the same lattice

    Args:
      other: Another LatticeReconstruction

    Returns:
      (bool, transform): Tuple of match success and transformation matrix
                        to align other to self
    """
    # Find common triangles in transformed space
    # Determine if alignment is possible
    # Return transformation matrix if successful
    return False, None

  def extend(self, other, transform):
    """Combine another reconstruction with this one

    Args:
      other: Another LatticeReconstruction
      transform: Transformation matrix from match()
    """
    # Apply transformation to other's reflections
    # Merge experiments, reflections, used_triangles
    # Invalidate clusters
    self.clusters = None

class SelectionDone(Exception):
  """Raised when user requests early termination of selection process"""
  pass

# Convenience function for axis labels
def q_to_d(x, _):
  return f"{1/x:.3f}"

class TripletData:
  """Generate and store triplets of spots with their geometric relationships
  """
  def __init__(self, experiments, reflections, params):
    """
    Args:
      experiments: dxtbx experiment list
      reflections: reflection table
      params: phil parameters
    """
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    self.triplets = None

    if params.triplets.load_npz is not None:
      self.load_triplets()
    else:
      self.compute()

  def load_triplets(self):
    """Load pre-computed triplets from npz file"""
    print(f"Loading triplets from {self.params.triplets.load_npz}")
    loaded = np.load(self.params.triplets.load_npz)
    self.triplets = loaded['triplets']
    print(f"Loaded {len(self.triplets)} triplets")

  def save_triplets(self):
    """Save triplets to npz file if requested"""
    if self.params.triplets.save_npz is not None:
      print(f"Saving triplets to {self.params.triplets.save_npz}")
      np.savez(self.params.triplets.save_npz, triplets=self.triplets)

  def compute(self):
    """Compute all d1,d2,theta triplets from input data

    Fills triplets array with columns:
      (frame_id, d1, d2, theta, spot1_idx, spot2_idx, hand)
    """
    if self.triplets is not None:
      return

    print(f"Computing triplets with d_min = {self.params.triplets.d_min}")

    triplets = []

    for i_expt, expt in enumerate(self.experiments):
      if i_expt % 1000 == 0:
        print(f"Processing frame {i_expt}")

      # Get spots for this frame
      sel = flex.bool(self.reflections['id'] == i_expt)
      frame_refls = self.reflections.select(sel)

      # Get detector coordinates relative to beam center
      xyz = frame_refls['xyzobs.mm.value']
      beam_x, beam_y = expt.detector[0].get_beam_centre(expt.beam.get_s0())
      x = xyz.parts()[0] - beam_x
      y = xyz.parts()[1] - beam_y

      # Compute d-spacings and filter
      s0 = expt.beam.get_s0()
      s1_vectors = frame_refls['s1'].as_numpy_array()
      d_spacings = np.array([1/np.linalg.norm(s1 - s0) for s1 in s1_vectors])

      valid_spots = d_spacings >= self.params.triplets.d_min
      if not np.any(valid_spots):
        continue

      d_spacings = d_spacings[valid_spots]
      s1_vectors = s1_vectors[valid_spots]
      spot_indices = np.arange(len(frame_refls))[valid_spots]
      flex_valid_spots = flumpy.from_numpy(valid_spots)
      x = x.select(flex_valid_spots)
      y = y.select(flex_valid_spots)

      # Compute all pairs
      n_spots = len(d_spacings)
      for i in range(n_spots):
        for j in range(i+1, n_spots):
          # Compute angle
          vec_a = s1_vectors[i] - s0
          vec_b = s1_vectors[j] - s0
          cos_angle = np.dot(vec_a, vec_b) / (np.linalg.norm(vec_a) * np.linalg.norm(vec_b))
          angle = np.rad2deg(np.arccos(min(1.0, max(-1.0, cos_angle))))

          # Make sure d1 > d2 and track coordinates
          if d_spacings[i] > d_spacings[j]: # correct order
            d1 = d_spacings[i]
            d2 = d_spacings[j]
            i1 = spot_indices[i]
            i2 = spot_indices[j]
            x1, y1 = x[i], y[i]
            x2, y2 = x[j], y[j]
          else: # Switch them
            d1 = d_spacings[j]
            d2 = d_spacings[i]
            i1 = spot_indices[j]
            i2 = spot_indices[i]
            x1, y1 = x[j], y[j]
            x2, y2 = x[i], y[i]

          # Determine handedness from cross product z-component
          # For d1 > d2:
          # cross_z > 0 means counterclockwise rotation from 1 to 2 (left-handed)
          cross_z = x1*y2 - y1*x2
          hand = -1 if cross_z > 0 else 1

          triplets.append((i_expt, d1, d2, angle, i1, i2, hand))

    self.triplets = np.array(triplets)

    print(f"Found {len(self.triplets)} triplets from {len(self.experiments)} frames")
    print(f"Average triplets per frame: {len(self.triplets)/len(self.experiments):.1f}")

    self.save_triplets()

    
  def find_matches(self, spec):
    """Find frames containing triplets that match the specification
    
    Args:
      spec: TriangleSpec object
    
    Returns:
      dict mapping frame_ids to lists of (spot1_idx, spot2_idx) pairs
    """
    print(f"\nSearching triplets in ranges:")
    print(f"  d1: {1/spec.q1:.2f}Å")
    print(f"  d2: {1/spec.q2:.2f}Å")
    print(f"  theta: {spec.theta:.1f}°")
    if self.triplets is None:
      self.compute()
      
    # Test all triplets against specification
    matches = spec.test(self.triplets[:,1],  # d1
                               self.triplets[:,2],  # d2
                               self.triplets[:,3])  # theta
    # Filter by hand if specified
    if spec.hand is not None:
      matches &= self.triplets[:,6] == spec.hand
    print(f"Found {np.sum(matches)} matching triplets")
    
    # Organize matching triplets by frame
    result = {}
    for triplet in self.triplets[matches]:
      frame_id = int(triplet[0])
      spot_pair = (int(triplet[4]), int(triplet[5]))
      if frame_id not in result:
        result[frame_id] = []
      result[frame_id].append(spot_pair)
      
    return result

  def interactive_clustering(self):
    """Interactive clustering guided by user selection"""

    if self.params.initialization.triangle_finding.interactive.gmm_params.use:
      # Use provided parameters
      gmm_params = self.params.initialization.triangle_finding.interactive.gmm_params
      cluster_results = []
      for q1, q2, r, th_c, th_w in zip(
        gmm_params.q1_centers,
        gmm_params.q2_centers,
        gmm_params.radii,
        gmm_params.theta_centers,
        gmm_params.theta_widths
      ):
        cluster_results.append({
          'q1_center': q1,
          'q2_center': q2,
          'radius': r,
          'theta_center': th_c,
          'theta_width': th_w
        })
    else:

      # Get q-magnitudes from both d1 and d2 columns
      q_mags = np.concatenate([
        1/self.triplets[:,1],  # from d1
        1/self.triplets[:,2]   # from d2
      ])

      # Create histogram
      hist, edges = np.histogram(q_mags, bins=2000)

      # Plot with d-spacing x-axis
      fig = plt.figure(figsize=(12,4))
      ax = plt.gca()
      ax.plot(edges[:-1], hist)

      # Format x-axis as d-spacing
      ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))
      ax.set_xlabel('d-spacing (Å)')
      ax.set_ylabel('Counts')
      ax.set_title('Select d-spacing ranges (press d when done, u to undo last selection)')

      # Store selected ranges
      q_ranges = []
      spans = []  # Store span patches for visualization

      def on_select(qmin, qmax):
        # Add range (sorted to ensure qmin < qmax)
        q_ranges.append(tuple(sorted([qmin, qmax])))
        # Add span patch to show selection
        span = ax.axvspan(qmin, qmax, color='red', alpha=0.2)
        spans.append(span)
        fig.canvas.draw_idle()

      selector = SpanSelector(
        ax, on_select, 'horizontal',
        props=dict(alpha=0.3, facecolor='red'),
        interactive=False,
        drag_from_anywhere=True
      )

      # Handle keyboard input
      selection_done = False
      def on_key(event):
        nonlocal selection_done
        if event.key == 'd':
          selection_done = True
          plt.close(fig)
        elif event.key == 'u' and q_ranges:  # Only if there are ranges to delete
          # Remove last range and span
          q_ranges.pop()
          spans[-1].remove()
          spans.pop()
          fig.canvas.draw_idle()

      fig.canvas.mpl_connect('key_press_event', on_key)

      plt.show()

      if not q_ranges:
        raise Sorry("No ranges selected")

      print(f"Selected {len(q_ranges)} ranges:")
      for qmin, qmax in q_ranges:
        print(f"  {1/qmax:.3f} - {1/qmin:.3f} Å")

      # Generate all pairs of ranges
      range_pairs = []
      for i in range(len(q_ranges)):
        for j in range(i+1, len(q_ranges)):
          range_pairs.append((q_ranges[i], q_ranges[j]))

      # Compute max range for consistent plotting
      max_range = max(q_max - q_min for q_min, q_max in q_ranges)
      print(f"\nGenerating {len(range_pairs)} q1-q2 plots...")

      cluster_results = []
      # Process each range pair
      for (q1min, q1max), (q2min, q2max) in range_pairs:
        clusters = self._select_clusters_for_ranges(
          (q1min, q1max),
          (q2min, q2max),
          max_range
        )
        cluster_results.extend(clusters)

      # Create labels array (all noise to start)
      labels = np.full(len(self.triplets), -1)

      # Assign cluster labels based on points_mask
      for i, cluster in enumerate(cluster_results):
        labels[cluster['points_mask']] = i

      # Store in expected format
      self.cluster_data = {
        'labels': labels,
        'points': np.column_stack((
          1/self.triplets[:,1],  # q1
          1/self.triplets[:,2],  # q2
          self.triplets[:,3]     # theta
        ))
      }

      # Store centers in expected format
      self.cluster_centers = {}
      for i, cluster in enumerate(cluster_results):
        self.cluster_centers[i] = np.array([
          cluster['q1_center'],
          cluster['q2_center'],
          cluster['theta_center']
        ])

      # Use all points (no sampling)
      self.sampled_indices = np.arange(len(self.triplets))

      self.cluster_results = cluster_results


  def _select_clusters_for_ranges(self, q1_range, q2_range, max_range):
    """Interactive cluster selection for a given q1,q2 range pair"""
    fig = plt.figure(figsize=(8,8))
    ax = plt.gca()
    
    # Select points in these ranges
    q1min, q1max = q1_range
    q2min, q2max = q2_range
    mask = ((1/self.triplets[:,1] >= q1min) & (1/self.triplets[:,1] <= q1max) &
            (1/self.triplets[:,2] >= q2min) & (1/self.triplets[:,2] <= q2max))
    
    # Get points
    q1 = 1/self.triplets[mask,1]
    q2 = 1/self.triplets[mask,2]
    
    # Set axis center points
    q1_center = (q1min + q1max) / 2
    q2_center = (q2min + q2max) / 2
    
    while True:  # Main loop for 2D fitting
      # Circle selection phase
      circles = []
      drawing = False
      center = None
      current_circle = None
      selection_done = False
      
      # Clear and reset plot
      ax.clear()
      ax.scatter(q1, q2, alpha=0.1, s=1)
      ax.set_xlim(q1_center - max_range/2, q1_center + max_range/2)
      ax.set_ylim(q2_center - max_range/2, q2_center + max_range/2)
      ax.get_xaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))
      ax.get_yaxis().set_major_formatter(tick.FuncFormatter(
        lambda x, _: "{:.3f}".format(1/x)))
      ax.set_xlabel('d1 (Å)')
      ax.set_ylabel('d2 (Å)')
      ax.set_title(f"Draw circles. u to undo, d when done")
      plt.draw()
      
      def on_press(event):
        nonlocal drawing, center, current_circle
        if event.inaxes != ax:
          return
        drawing = True
        center = (event.xdata, event.ydata)
        current_circle = plt.Circle(center, 0, fill=False, color='red')
        ax.add_patch(current_circle)
      
      def on_motion(event):
        nonlocal current_circle
        if not drawing or event.inaxes != ax:
          return
        dx = event.xdata - center[0]
        dy = event.ydata - center[1]
        radius = np.sqrt(dx*dx + dy*dy)
        current_circle.set_radius(radius)
        fig.canvas.draw_idle()
      
      def on_release(event):
        nonlocal drawing, current_circle
        if not drawing:
          return
        drawing = False
        if current_circle is not None:
          circles.append(current_circle)
          current_circle = None
      
      def on_key(event):
        nonlocal selection_done
        if event.key == 'u' and circles:  # Delete last circle
          circles[-1].remove()
          circles.pop()
          fig.canvas.draw_idle()
        elif event.key == 'd':  # Done selecting
          selection_done = True
      
      # Connect event handlers
      fig.canvas.mpl_connect('button_press_event', on_press)
      fig.canvas.mpl_connect('button_release_event', on_release)
      fig.canvas.mpl_connect('motion_notify_event', on_motion)
      fig.canvas.mpl_connect('key_press_event', on_key)
      
      while not selection_done:
        plt.pause(0.1)
      
      if not circles:
        plt.close(fig)
        return []
      
      fitting_response = None
      def on_fit_key(event):
        nonlocal fitting_response
        if event.key == 'y':
          fitting_response = 'accept'
        elif event.key == 'n':
          fitting_response = 'reject'
        elif event.key == 'r':
          fitting_response = 'retry'
      
      # Collect points from circles and fit Gaussians directly
      means = []
      sigmas = []
      for circle in circles:
        center = circle.center
        radius = circle.radius
        dist_sq = (q1 - center[0])**2 + (q2 - center[1])**2
        circle_mask = dist_sq <= radius**2
        points_2d = np.column_stack((q1[circle_mask], q2[circle_mask]))

        # Simple mean and standard deviation for each dimension
        mean = np.mean(points_2d, axis=0)
        sigma = np.std(points_2d, axis=0)
        means.append(mean)
        sigmas.append(sigma)

      # Remove old circles
      for circle in circles:
        circle.remove()

      # Draw confidence ellipses
      threshold = self.params.initialization.triangle_finding.mahalanobis_threshold
      from matplotlib.patches import Ellipse
      for mean, sigma in zip(means, sigmas):
        ellipse = Ellipse(
          xy=mean,
          width=2*sigma[0]*threshold,
          height=2*sigma[1]*threshold,
          fill=False, color='red'
        )
        ax.add_artist(ellipse)

      
      ax.set_title("Press: y to accept, n to reject, r to retry fit")
      plt.draw()
      
      # Connect key handler and wait for response
      cid = fig.canvas.mpl_connect('key_press_event', on_fit_key)
      while fitting_response is None:
        plt.pause(0.1)
      fig.canvas.mpl_disconnect(cid)
      
      if fitting_response == 'accept':
        break
      elif fitting_response == 'reject':
        plt.close(fig)
        return []
      else:  # retry
        continue
    
    # Now proceed with theta fitting
    # After 2D fits are accepted...
    theta_results = []

    for i, (mean, sigma) in enumerate(zip(means, sigmas)):
      # Get points near this cluster center
      dist_sq = (q1 - mean[0])**2 + (q2 - mean[1])**2
      cluster_mask = dist_sq <= (2*max(sigma))**2
      theta_values = self.triplets[mask][cluster_mask,3]

      # Create theta histogram figure
      fig_theta = plt.figure(figsize=(8,4))
      ax_theta = plt.gca()

      hist, edges = np.histogram(theta_values, bins=180, range=(0,180))

      while True:  # Main loop for theta selection
        # Clear everything by clearing axis
        ax_theta.clear()
        spans = []  # Reset spans list
        theta_ranges = []  # Reset ranges

        # Create/recreate histogram
        ax_theta.bar(edges[:-1], hist, width=np.diff(edges), align='edge')
        ax_theta.set_xlabel('Angle (degrees)')
        ax_theta.set_ylabel('Counts')
        ax_theta.set_title(f'Select angle range for cluster {i+1}\n' +
                          f'd1≈{1/mean[0]:.2f}Å, d2≈{1/mean[1]:.2f}Å')
        plt.draw()

        # Get initial ranges from user
        spans = []
        selection_done = False

        def on_select(tmin, tmax):
          theta_ranges.append(tuple(sorted([tmin, tmax])))
          span = ax_theta.axvspan(tmin, tmax, color='red', alpha=0.2)
          spans.append(span)
          fig_theta.canvas.draw_idle()

        selector = SpanSelector(
          ax_theta, on_select, 'horizontal',
          props=dict(alpha=0.3, facecolor='red'),
          interactive=False,
          drag_from_anywhere=True
        )

        # Handle keyboard input
        def on_key(event):
          nonlocal selection_done
          if event.key == 'u' and theta_ranges:  # Delete last range
            theta_ranges.pop()
            spans[-1].remove()
            spans.pop()
            fig_theta.canvas.draw_idle()
          elif event.key == 'd':  # Done selecting
            selection_done = True

        fig_theta.canvas.mpl_connect('key_press_event', on_key)

        while not selection_done:
          plt.pause(0.1)

        if not theta_ranges:
          plt.close(fig_theta)
          break

        # Fit sum of Gaussians with bounds
        n_gaussians = len(theta_ranges)
        init_params = []
        bounds = []
        for tmin, tmax in theta_ranges:
          mean_ = (tmin+tmax)/2
          std_ = (tmax-tmin)/4
          init_params.extend([
            (tmin + tmax) / 2,  # mean
            (tmax - tmin) / 4,  # std
            max(hist)           # scale
          ])
          bounds.extend([
            (None, None),      # mean can be any value
            (std_-1e-5,std_+1e-5),       # std must be positive
            (0.1, None)        # scale must be positive
          ])

        def multi_gaussian(x, params):
          """Sum of n Gaussians"""
          result = np.zeros_like(x)
          for i in range(n_gaussians):
            mean = params[i*3]
            std = params[i*3 + 1]
            scale = params[i*3 + 2]
            result += scale * np.exp(-(x - mean)**2 / (2*std**2))
          return result

        def objective(params):
          """Least squares between histogram and sum of Gaussians"""
          predicted = multi_gaussian(edges[:-1], params)
          return np.sum((hist - predicted)**2)

        result = minimize(
          objective,
          init_params,
          method='L-BFGS-B',  # Method that handles bounds
          bounds=bounds
        )
        fitted_params = result.x

        # Plot fit
        ax_theta.clear()
        ax_theta.bar(edges[:-1], hist, width=np.diff(edges), align='edge')
        x = np.linspace(0, 180, 1000)
        y_total = multi_gaussian(x, fitted_params)
        ax_theta.plot(x, y_total, 'r-', label='Total fit')

        # Plot individual components
        colors = ['g--', 'b--', 'c--', 'm--']  # For up to 4 components
        for j in range(n_gaussians):
          component_params = np.zeros(len(fitted_params))
          component_params[j*3:(j+1)*3] = fitted_params[j*3:(j+1)*3]
          y_component = multi_gaussian(x, component_params)
          ax_theta.plot(x, y_component, colors[j % len(colors)],
                       label=f'Component {j+1}')

        ax_theta.legend()
        ax_theta.set_title(f'Fitted angle distribution for cluster {i+1}\n' +
                          f'd1≈{1/mean[0]:.2f}Å, d2≈{1/mean[1]:.2f}Å\n' +
                          'Press: y to accept, n to reject, r to retry selection')

        # Handle accept/reject/retry
        fitting_response = None
        def on_fit_key(event):
          nonlocal fitting_response
          if event.key == 'y':
            fitting_response = 'accept'
          elif event.key == 'n':
            fitting_response = 'reject'
          elif event.key == 'r':
            fitting_response = 'retry'

        cid = fig_theta.canvas.mpl_connect('key_press_event', on_fit_key)
        plt.draw()

        while fitting_response is None:
          plt.pause(0.1)

        fig_theta.canvas.mpl_disconnect(cid)

        if fitting_response == 'accept':
          # Store parameters for all Gaussian components
          components = []
          for j in range(n_gaussians):
            components.append({
              'mean': fitted_params[j*3],
              'std': fitted_params[j*3 + 1],
              'scale': fitted_params[j*3 + 2]
            })
          theta_results.append(components)
          break
        elif fitting_response == 'reject':
          break
        # else retry - loop will continue

      plt.close(fig_theta)

    plt.close(fig) # close the q1-q2 plot
    if not theta_results:
      return []

    # Combine 2D and theta results
    cluster_results = []
    for (mean_2d, sigma_2d), theta_components in zip(zip(means, sigmas), theta_results):
      # For each Gaussian component in theta fit
      for theta_params in theta_components:
        theta_mean = theta_params['mean']
        theta_sigma = theta_params['std']
        theta_scale = theta_params['scale']

        # Create 3D covariance matrix
        cov_3d = np.zeros((3,3))
        cov_3d[0,0] = sigma_2d[0]**2  # q1 variance
        cov_3d[1,1] = sigma_2d[1]**2  # q2 variance
        cov_3d[2,2] = theta_sigma**2   # theta variance

        # Find points within 1 sigma using Mahalanobis distance
        center_3d = np.array([mean_2d[0], mean_2d[1], theta_mean])
        inv_cov = np.linalg.inv(cov_3d)

        # Get all points near this cluster
        points_3d = np.column_stack((
          1/self.triplets[:,1],  # q1
          1/self.triplets[:,2],  # q2
          self.triplets[:,3]     # theta
        ))

        # Compute Mahalanobis distances
        diff = points_3d - center_3d
        dist_sq = np.sum(diff @ inv_cov * diff, axis=1)
        cluster_mask = dist_sq <= 1.0  # One sigma

        cluster_results.append({
          'q1_center': mean_2d[0],
          'q2_center': mean_2d[1],
          'q1_sigma': sigma_2d[0],
          'q2_sigma': sigma_2d[1],
          'theta_center': theta_mean,
          'theta_sigma': theta_sigma,
          'theta_scale': theta_scale,
          'covariance': cov_3d,
          'points_mask': cluster_mask
        })

    return cluster_results

  def cluster_triplets(self, max_points=50000, dmin=None, dmax=None, 
                      d_eps=0.0003, min_d_samples=100,  # for d-spacing clustering
                      angle_eps=.5, min_angle_samples=25):  # for angle subclustering
    """Find clusters in reciprocal space using two-stage clustering
    
    Args:
      max_points: Maximum points to sample for clustering
      dmin, dmax: Optional d-spacing range to analyze
      d_eps: DBSCAN epsilon for d1,d2 clustering (in q-space)
      min_d_samples: Minimum points for d1,d2 clusters
      angle_eps: Angular separation (degrees) for subclustering
      min_angle_samples: Minimum points for angle subclusters
    """

    DEBUG_FIRST_STAGE_ONLY = False  # Add this flag

    # Data selection as before
    q1 = 1/self.triplets[:,1]
    q2 = 1/self.triplets[:,2]
    theta = self.triplets[:,3]
    
    # Apply d-spacing filter and sampling as before
    selected = np.arange(len(q1))
    if dmin is not None or dmax is not None:
      dmin = dmin or 0
      dmax = dmax or float('inf')
      mask = ((1/q1 >= dmin) & (1/q1 <= dmax) &
              (1/q2 >= dmin) & (1/q2 <= dmax))
      selected = np.where(mask)[0]
      print(f"Selected {len(selected)} points in d-spacing range {dmin:.2f}-{dmax:.2f}Å")
    
    if len(selected) > max_points:
      self.sampled_indices = np.random.choice(selected, size=max_points, replace=False)
      print(f"Randomly sampled {max_points} points for clustering")
    else:
      self.sampled_indices = selected
      print(f"Using all {len(selected)} points for clustering")
    
    # Get selected points
    q1 = q1[self.sampled_indices]
    q2 = q2[self.sampled_indices]
    theta = theta[self.sampled_indices]
    
    # First stage: cluster in d1,d2 space
    d_points = np.column_stack((q1, q2))
    d_clustering = DBSCAN(eps=d_eps, min_samples=min_d_samples).fit(d_points)
    #d_clustering = HDBSCAN(min_cluster_size=min_d_samples, min_samples=25).fit(d_points)
    
    if DEBUG_FIRST_STAGE_ONLY:
      # Just use d-clustering labels directly
      labels = d_clustering.labels_
    else:
      # Second stage: subcluster by angle within each d1,d2 cluster
      labels = np.full(len(q1), -1)  # Start all as noise
      next_label = 0
      
      for d_label in set(d_clustering.labels_) - {-1}:
        d_mask = d_clustering.labels_ == d_label
        cluster_thetas = theta[d_mask]
        
        # Cluster angles within this d1,d2 cluster
        angle_clustering = DBSCAN(
          eps=angle_eps, 
          min_samples=min_angle_samples
        ).fit(cluster_thetas.reshape(-1, 1))
        
        # Assign new labels to subclusters
        for a_label in set(angle_clustering.labels_) - {-1}:
          subcluster_mask = d_mask.copy()
          subcluster_mask[d_mask] = angle_clustering.labels_ == a_label
          labels[subcluster_mask] = next_label
          next_label += 1
    
    # Store results as before
    self.cluster_data = {
      'labels': labels,
      'points': np.column_stack((q1, q2, theta))
    }
    
    # Compute and store cluster statistics
    self.cluster_centers = {}
    self.cluster_covariances = {}
    
    for k in set(labels) - {-1}:
      mask = labels == k
      cluster_points = self.cluster_data['points'][mask]
      self.cluster_centers[k] = np.mean(cluster_points, axis=0)
      self.cluster_covariances[k] = np.cov(cluster_points, rowvar=False)



  def plot(self, plot_clusters=False, n_max=None):
    """Interactive plot of triplet data
    
    Args:
      plot_clusters: If True, color by cluster assignment
      n_max: If set, randomly subsample to this many points
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as tick
    import numpy as np
    
    # Get data
    q1 = 1/self.triplets[:,1]
    q2 = 1/self.triplets[:,2]
    theta = self.triplets[:,3]
    
#    # Subsample if requested
#    if n_max and len(q1) > n_max:
#      idx = np.random.choice(len(q1), size=n_max, replace=False)
#      q1 = q1[idx]
#      q2 = q2[idx]
#      theta = theta[idx]
    
    if plot_clusters:
      if not hasattr(self, 'cluster_data'):
        raise Sorry("No clustering results available. Run cluster_triplets first.")

      # Print cluster statistics
      labels = self.cluster_data['labels']
      print("\nCluster populations:")

      # Count points in each cluster
      cluster_sizes = {}
      for k in set(labels) - {-1}:
        cluster_sizes[k] = np.sum(labels == k)

      # Sort by size and print top 10
      print("\nTop 10 clusters by population:")
      print("Cluster  Points  Center (d1 Å, d2 Å, theta°)")
      print("-" * 45)
      for k, size in sorted(cluster_sizes.items(), key=lambda x: x[1], reverse=True)[:10]:
        center = self.cluster_centers[k]
        print(f"{k:7d}  {size:6d}  ({1/center[0]:.2f}, {1/center[1]:.2f}, {center[2]:.1f})")

      # Number of noise points
      n_noise = np.sum(labels == -1)
      print(f"\nNoise points: {n_noise} ({100*n_noise/len(labels):.1f}%)")

      colors = plt.cm.rainbow(np.linspace(0, 1, len(cluster_sizes)))

    # Apply sampled_indices mask
    if hasattr(self, 'sampled_indices'):
      q1 = q1[self.sampled_indices]
      q2 = q2[self.sampled_indices]
      theta = theta[self.sampled_indices]
      if plot_clusters:
        labels = labels[self.sampled_indices]

    # Create figure without tight_layout
    fig = plt.figure(figsize=(12, 8))

    # Create axes in triangle layout
    ax1 = fig.add_axes([0.1, 0.1, 0.35, 0.35])     # bottom left (d1,d2)
    ax2 = fig.add_axes([0.6, 0.1, 0.35, 0.35])     # bottom right (theta,d2)
    ax3 = fig.add_axes([0.1, 0.55, 0.35, 0.35])    # top left (d1,theta)

    # Create plots
    if plot_clusters:
      # Plot noise points
      noise_mask = labels == -1
      scatter1_points = []
      scatter2_points = []
      scatter3_points = []

      scatter1_points.append(ax1.scatter(q1[noise_mask], q2[noise_mask],
                                       c='gray', alpha=0.5, s=1))
      scatter2_points.append(ax2.scatter(theta[noise_mask], q2[noise_mask],
                                       c='gray', alpha=0.5, s=1))
      scatter3_points.append(ax3.scatter(q1[noise_mask], theta[noise_mask],
                                       c='gray', alpha=0.5, s=1))

      # Plot each cluster
      for i, k in enumerate(sorted(cluster_sizes.keys())):
        mask = labels == k
        scatter1_points.append(ax1.scatter(q1[mask], q2[mask],
                                         c=[colors[i]], alpha=0.5, s=1))
        scatter2_points.append(ax2.scatter(theta[mask], q2[mask],
                                         c=[colors[i]], alpha=0.5, s=1))
        scatter3_points.append(ax3.scatter(q1[mask], theta[mask],
                                         c=[colors[i]], alpha=0.5, s=1))
    else:
      scatter1_points = [ax1.scatter(q1, q2, c='blue', alpha=0.5, s=1)]
      scatter2_points = [ax2.scatter(theta, q2, c='blue', alpha=0.5, s=1)]
      scatter3_points = [ax3.scatter(q1, theta, c='blue', alpha=0.5, s=1)]


    # Set initial axis limits
    ax1.set_xlim(min(q1), max(q1))
    ax1.set_ylim(min(q2), max(q2))
    ax2.set_xlim(min(theta), max(theta))
    ax2.set_ylim(min(q2), max(q2))
    ax3.set_xlim(min(q1), max(q1))
    ax3.set_ylim(min(theta), max(theta))

    # Store text annotations
    text_boxes = {'ax1': None, 'ax2': None, 'ax3': None}

    # Set up axes
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(q_to_d))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(q_to_d))
    ax1.set_xlabel('d1 (Å)')
    ax1.set_ylabel('d2 (Å)')

    # Bottom right plot (theta,d2)
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(q_to_d))
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()
    ax2.set_xlabel('Angle (degrees)')
    ax2.set_ylabel('d2 (Å)')

    # Top left plot (d1,theta)
    ax3.xaxis.set_major_formatter(tick.FuncFormatter(q_to_d))
    ax3.xaxis.set_label_position('top')
    ax3.xaxis.tick_top()
    ax3.set_xlabel('d1 (Å)')
    ax3.set_ylabel('Angle (degrees)')

    def update_text_boxes(event_ax):
      # Clear previous text boxes
      for txt in text_boxes.values():
        if txt is not None:
          txt.remove()

      # Get current visible data based on axis limits
      xlim1, ylim1 = ax1.get_xlim(), ax1.get_ylim()
      xlim2, ylim2 = ax2.get_xlim(), ax2.get_ylim()
      xlim3, ylim3 = ax3.get_xlim(), ax3.get_ylim()

      # Use same masking logic as update_plots
      if event_ax == ax1:  # d1,d2 plot
        d1_lim = xlim1
        d2_lim = ylim1
        theta_lim = xlim2
      elif event_ax == ax2:  # theta,d2 plot
        theta_lim = xlim2
        d2_lim = ylim2
        d1_lim = xlim1
      elif event_ax == ax3:  # d1,theta plot
        d1_lim = xlim3
        theta_lim = ylim3
        d2_lim = ylim1

      # Apply all three filters
      mask = ((q1 >= d1_lim[0]) & (q1 <= d1_lim[1]) &
              (q2 >= d2_lim[0]) & (q2 <= d2_lim[1]) &
              (theta >= theta_lim[0]) & (theta <= theta_lim[1]))

      # Update text boxes with ranges from filtered data
      visible_theta = theta[mask]
      visible_q1 = q1[mask]
      visible_q2 = q2[mask]

      if len(visible_theta) > 0:
        text_boxes['ax1'] = ax1.text(0.02, 0.98,
          f'θ: {min(visible_theta):.1f}°-{max(visible_theta):.1f}°',
          transform=ax1.transAxes, va='top')

      if len(visible_q1) > 0:
        text_boxes['ax2'] = ax2.text(0.02, 0.98,
          f'd1: {1/max(visible_q1):.1f}Å-{1/min(visible_q1):.1f}Å',
          transform=ax2.transAxes, va='top')

      if len(visible_q2) > 0:
        text_boxes['ax3'] = ax3.text(0.02, 0.98,
          f'd2: {1/max(visible_q2):.1f}Å-{1/min(visible_q2):.1f}Å',
          transform=ax3.transAxes, va='top')

    updating = False
    def update_plots(event_ax):
      """Update all plots based on zoom in event_ax"""
      nonlocal updating
      if updating: return
      updating = True

      try:
        xlim = event_ax.get_xlim()
        ylim = event_ax.get_ylim()

        if event_ax == ax1:  # d1,d2 plot
          d1_lim = xlim
          d2_lim = ylim
          theta_lim = ax2.get_xlim()
        elif event_ax == ax2:  # theta,d2 plot
          theta_lim = xlim
          d2_lim = ylim
          d1_lim = ax1.get_xlim()
        elif event_ax == ax3:  # d1,theta plot
          d1_lim = xlim
          theta_lim = ylim
          d2_lim = ax1.get_ylim()

        # Apply all three filters to every update
        mask = ((q1 >= d1_lim[0]) & (q1 <= d1_lim[1]) &
                (q2 >= d2_lim[0]) & (q2 <= d2_lim[1]) &
                (theta >= theta_lim[0]) & (theta <= theta_lim[1]))

        # Update all scatter plots
        if plot_clusters:
          # Update each scatter collection
          for i, k in enumerate([-1] + sorted(cluster_sizes.keys())):
            submask = mask & (labels == k)
            scatter1_points[i].set_offsets(np.c_[q1[submask], q2[submask]])
            scatter2_points[i].set_offsets(np.c_[theta[submask], q2[submask]])
            scatter3_points[i].set_offsets(np.c_[q1[submask], theta[submask]])
        else:
          scatter1_points[0].set_offsets(np.c_[q1[mask], q2[mask]])
          scatter2_points[0].set_offsets(np.c_[theta[mask], q2[mask]])
          scatter3_points[0].set_offsets(np.c_[q1[mask], theta[mask]])

        # Update all axis limits
        if event_ax == ax1:
          ax3.set_xlim(xlim)        # d1 limits
          ax2.set_ylim(ylim)        # d2 limits
        elif event_ax == ax2:
          ax1.set_ylim(ylim)        # d2 limits
          ax3.set_ylim(xlim)        # theta limits
        elif event_ax == ax3:
          ax1.set_xlim(xlim)        # d1 limits
          ax2.set_xlim(ylim)        # theta limits

        fig.canvas.draw_idle()
        update_text_boxes(event_ax)
      finally:
        updating = False

    
    # Connect event handlers
    for ax in [ax1, ax2, ax3]:
      ax.callbacks.connect('xlim_changed', lambda evt: update_plots(evt.axes))
      ax.callbacks.connect('ylim_changed', lambda evt: update_plots(evt.axes))
    
    print('show...')
    plt.show()
    print('showed')

  def cluster(self):
    """Run clustering based on method in params"""
    method = self.params.initialization.triangle_finding.method
    if method == 'interactive':
      return self.interactive_clustering()
    elif method == 'clustering':
      return self.automatic_clustering()
    elif method == 'box':
      return self.box_clustering()
    else:
      raise Sorry(f"Unknown clustering method: {method}")


  # Some methods to allow checkpointing after interactive run...
  def compute_checksum(self):
    """Compute checksum of triplets array for validation"""
    return hashlib.sha256(self.triplets.tobytes()).hexdigest()

  def save_clusters(self, filename):
    """Save clustering results and metadata"""
    data = {
      'cluster_data': self.cluster_data,
      'cluster_centers': self.cluster_centers,
      'cluster_results': self.cluster_results,
      'sampled_indices': self.sampled_indices,
      'triplets_checksum': self.compute_checksum()
    }
    with open(filename, 'wb') as f:
      pickle.dump(data, f)

  def load_clusters(self, filename):
    """Load and validate clustering results"""
    with open(filename, 'rb') as f:
      data = pickle.load(f)
    if data['triplets_checksum'] != self.compute_checksum():
      raise Sorry("Triplets data does not match saved results")
    self.cluster_data = data['cluster_data']
    self.cluster_centers = data['cluster_centers']
    self.cluster_results = data['cluster_results']
    self.sampled_indices = data['sampled_indices']



class TriangleSpec:
  """Base class for triangle specification and matching
  """
  def __init__(self, q1, q2, theta, hand=None):
    self.q1 = q1
    self.q2 = q2
    self.theta = theta
    self.hand = hand # -1 for left-handed, 1 for right-handed, None for either

  @classmethod
  def from_cluster(cls, cluster, params, hand=None):
    """Factory method to create appropriate TriangleSpec from cluster data"""
    if 'covariance' in cluster:
      # For error propagation: σ_d = σ_q * (d^2)
      q1, q2 = cluster['q1_center'], cluster['q2_center']
      d1, d2 = 1/q1, 1/q2
      sigma_d1 = cluster['q1_sigma'] * (d1**2)
      sigma_d2 = cluster['q2_sigma'] * (d2**2)

      print(f"Creating GaussianTriangleSpec with:")
      print(f"  center: ({d1:.2f}Å, {d2:.2f}Å, {cluster['theta_center']:.1f}°)")
      print(f"  sigmas (q-space): ({cluster['q1_sigma']:.5f}, {cluster['q2_sigma']:.5f}, {cluster['theta_sigma']:.1f})")
      print(f"  sigmas (d-space): ({sigma_d1:.3f}Å, {sigma_d2:.3f}Å, {cluster['theta_sigma']:.1f}°)")

      return GaussianTriangleSpec(
        center=[q1, q2, cluster['theta_center']],
        covariance=cluster['covariance'],
        params=params,
        hand=hand if hand is not None else cluster.get('hand')
      )
    # ... other specializations based on cluster properties

  def test(self, d1, d2, theta):
    """
    Test if a given triangle matches this specification

    Args:
      d1, d2: d-spacings of two spots
      theta: angle between spots (degrees)

    Returns:
      bool: Whether the triangle matches the specification
    """
    raise NotImplementedError()

  def reconstruct(self, triplet_data, experiments, reflections):
    """Create LatticeReconstruction from spots matching this triangle spec

    Args:
      triplet_data: TripletData object containing matching information
      experiments: dxtbx experiment list
      reflections: reflection table

    Returns:
      LatticeReconstruction object
    """
    # Find matching frames using this spec
    matches = triplet_data.find_matches(self)
    print(f"\nFound {len(matches)} frames with matching triangles:")
    print("Frame IDs:", sorted(matches.keys()))

    # Create reconstruction
    reconstruction = LatticeReconstruction(self.params)
    reconstruction.hand = self.hand

    # Add each matching frame
    for frame_id, spot_pairs in matches.items():
      expt = experiments[frame_id]
      sel = flex.bool(reflections['id'] == frame_id)
      refl = reflections.select(sel)

      # For now, just use first matching triangle from each frame
      i1, i2 = spot_pairs[0]
      #reconstruction.add_2d(expt, refl, i1, i2)
      for i1, i2 in spot_pairs:
        reconstruction.add(expt, refl, i1, i2)

    return reconstruction



class BoxTriangleSpec(TriangleSpec):
  """Triangle specification using box ranges
  """
  def __init__(self, d1_range, d2_range, theta_range, params):
    self.d1_min, self.d1_max = sorted(d1_range)
    self.d2_min, self.d2_max = sorted(d2_range)
    self.theta_min, self.theta_max = sorted(theta_range)
    # Use midpoints for base class center values
    super().__init__(
      (self.d1_min + self.d1_max)/2,
      (self.d2_min + self.d2_max)/2,
      (self.theta_min + self.theta_max)/2
    )

  def test(self, d1, d2, theta):
    """Test if given triangles match this specification

    Args:
      d1, d2: arrays of d-spacings
      theta: array of angles (degrees)

    Returns:
      array of bool: Whether each triangle matches the specification
    """
    return ((self.d1_min <= d1) & (d1 <= self.d1_max) &
            (self.d2_min <= d2) & (d2 <= self.d2_max) &
            (self.theta_min <= theta) & (theta <= self.theta_max))


class GaussianTriangleSpec(TriangleSpec):
  """Triangle specification using Mahalanobis distance
  """
  def __init__(self, center, covariance, params, hand=None):
    super().__init__(*center)
    self.covariance = covariance
    self._inv_cov = np.linalg.inv(covariance)
    self.params = params
    self.threshold = params.initialization.triangle_finding.mahalanobis_threshold

  def test(self, d1, d2, theta):
    """
    Returns array of bool for points within Mahalanobis distance threshold
    """
    q1, q2 = 1/d1, 1/d2
    x = np.vstack([q1, q2, theta]).T - np.array([self.q1, self.q2, self.theta])
    dist = np.sqrt(np.sum(x.dot(self._inv_cov) * x, axis=1))
    min_dist = np.min(dist)
    print(f"Minimum Mahalanobis distance: {min_dist:.3f}")
    matches = dist < self.threshold
    if not np.any(matches):
      print("WARNING: No points within threshold!")
    return matches



def find_triangle_spec(experiments, reflections, params):
  """Find a triangle specification from the data.
  Strategy determined by params.

  Args:
    experiments: dxtbx experiment list
    reflections: reflection table
    params: phil scope parameters

  Returns:
    TriangleSpec object
  """
  if params.initialization.triangle_finding.method == "box":
    box_params = params.initialization.triangle_finding.box
    if None in (box_params.d1_range, box_params.d2_range, box_params.theta_range):
      raise Sorry("Box method requires d1_range, d2_range, and theta_range to be specified")

    return BoxTriangleSpec(
      box_params.d1_range,
      box_params.d2_range,
      box_params.theta_range
    )

  elif params.initialization.triangle_finding.method == "interactive":
    raise NotImplementedError("Interactive triangle finding not yet implemented")

  elif params.initialization.triangle_finding.method == "clustering":
    raise NotImplementedError("Clustering-based triangle finding not yet implemented")

  else:
    raise Sorry(f"Unknown triangle finding method: {params.initialization.triangle_finding.method}")

def run_initialize(experiments, reflections, triplets, params):
  """Initialize lattice reconstruction from unindexed spots

  Args:
    experiments: dxtbx experiment list
    reflections: reflection table
    params: phil parameters
  """
  # Get triangle specification from parameters
  if params.initialization.triangle_finding.method == "box":
    box_params = params.initialization.triangle_finding.box
    if None in (box_params.d1_range, box_params.d2_range, box_params.theta_range):
      raise Sorry("Box method requires d1_range, d2_range, and theta_range to be specified")

    spec = BoxTriangleSpec(
      box_params.d1_range,
      box_params.d2_range,
      box_params.theta_range
    )

  else:
    raise Sorry("Only box method is currently implemented")


  # Find matching frames
  matches = triplets.find_matches(spec)
  print(f"\nFound {len(matches)} frames with matching triangles:")
  print("Frame IDs:", sorted(matches.keys()))

  # Create reconstruction
  reconstruction = LatticeReconstruction()

  # Add each matching frame
  for frame_id, spot_pairs in matches.items():
    expt = experiments[frame_id]
    sel = flex.bool(reflections['id'] == frame_id)
    refl = reflections.select(sel)

    # For now, just use the first matching triangle from each frame
    i1, i2 = spot_pairs[0]
    reconstruction.add(expt, refl, i1, i2)

  return reconstruction

#  for frame_id, spot_pairs in matches.items():
#    print(f"Frame {frame_id}: {len(spot_pairs)} matching triangles")


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
  """Reconstruct lattice from unindexed spots."""

  usage = "dials.reconstruct_lattice [options] experiments.expt reflections.refl"
  parser = ArgumentParser(
    usage=usage,
    phil=phil,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)
    
  params, options = parser.parse_args(args, show_diff_phil=True)
  assert len(params.input.experiments) == len(params.input.reflections) == 1
  experiments = params.input.experiments[0].data
  reflections = params.input.reflections[0].data

    
  if len(experiments) == 0 or len(reflections) == 0:
    parser.print_help()
    raise Sorry("No experiments or reflections found in input.")
      
  triplets = TripletData(experiments, reflections, params)
  plc = params.initialization.triangle_finding.load_clusters 
  psc = params.initialization.triangle_finding.save_clusters 
  if plc is not None:
    triplets.load_clusters(plc)
  else:
    assert psc is not None
    triplets.cluster()
    triplets.save_clusters(psc)

  manager = ReconstructionManager(experiments, reflections, triplets, params)
  manager.initialize(triplets.cluster_results)
  augment_done = False
  while not augment_done:
    done = manager.augment()
  manager.select_points()
  manager.select_angles()
  manager.refine()
  manager.plot_3d()
#  partial_lattices = []
#  reconstruction = LatticeReconstruction(params)
#  for cluster in triplets.cluster_results:
#    spec = TriangleSpec.from_cluster(cluster, params)
#    lr = spec.reconstruct(triplets, experiments, reflections)
#    partial_lattices.append(lr)
#  import IPython;IPython.embed()
#  exit()

  triplets.cluster_triplets(dmax=5.03,dmin=4.88)
  triplets.plot(plot_clusters=True)
  exit()

  if params.mode == "initialize":
    # Call reconstruction methods
    reconstruction = run_initialize(experiments, reflections, triplets, params)
    reconstruction.update_clusters()
    fig1 = reconstruction.plot_debug(title='before')
    for i in range(10):
      print(f'refine cycle {i}')
      reconstruction.refine_orientations()
      reconstruction.update_clusters()
    fig2 = reconstruction.plot_debug(title='after')
    plt.show()
  elif params.mode == "continue":
    # Handle continuation mode
    pass


if __name__ == "__main__":
  run()
