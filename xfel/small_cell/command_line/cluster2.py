import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde
from scipy.ndimage import maximum_filter
from scipy.ndimage import generate_binary_structure
from sklearn.neighbors import KernelDensity
import sys
from scipy.optimize import minimize

import matplotlib
matplotlib.use('TkAgg')


class ManualClusterer:
    def __init__(self, data, bandwidths=[0.001, 0.001, 1.0], n_maxima=500,
                 qvals_1=None, qvals_2=None, sb1_callback=None, recon=None,
                 qmin=.1, qmax=.5, points=None):
        self.data = data
        self.bandwidths = np.array(bandwidths)
        self.n_maxima = n_maxima
        self.qvals_1 = qvals_1
        self.qvals_2 = qvals_2
        self.sb1_qvals = None
        self.qmin = qmin
        self.qmax = qmax
        self.points = points

        self.current_selection = None
        self.final_selection = None
        self.means = None
        self.span_patch = None
        self.rect_patch = None
        self.selected_triplets = []

        self.sb1_callback = sb1_callback
        self.recon = recon
        
        # Compute KDE maxima
        if n_maxima > 0:
            self.compute_slice_maxima()
        

    def select_qvals(self, title="Select q-value ranges"):
        """
        Simple q-value selection mode using 1D histogram on dimension 0.
        
        Parameters:
        -----------
        title : str
            Title for the histogram window
        
        Returns:
        --------
        list : Selected q-values (medians of selected ranges)
        """
        # Create figure for histogram
        self.fig_qval, self.ax_qval = plt.subplots(num='Q-value Selection', figsize=(16, 4))
        
        # Initialize list to store selected q-values
        self.selected_qvals = []
        self.qval_spans = []  # Store span patches for visualization
        
        # Plot histogram
        self._plot_qval_histogram(title)
        
        # Set up span selector
        self.qval_span = SpanSelector(
            self.ax_qval,
            self._on_qval_span_select,
            'horizontal',
            useblit=True,
            props=dict(alpha=0.5, facecolor='blue')
        )
        
        # Add done button
        done_ax = self.fig_qval.add_axes([0.9, 0.01, 0.09, 0.05])
        self.qval_done_button = plt.Button(done_ax, 'Done')
        self.qval_done_button.on_clicked(self._on_qval_done)
        
        # Add clear button
        clear_ax = self.fig_qval.add_axes([0.8, 0.01, 0.09, 0.05])
        self.qval_clear_button = plt.Button(clear_ax, 'Clear Last')
        self.qval_clear_button.on_clicked(self._on_qval_clear_last)
        
        # Connect key press events
        self.fig_qval.canvas.mpl_connect('key_press_event', self._on_qval_key_press)
        
        # Show plot and wait for user interaction
        plt.show(block=True)
        
        # Close figure
        plt.close(self.fig_qval)
        
        return self.selected_qvals

    def _plot_qval_histogram(self, title):
        """Plot histogram for q-value selection."""
        self.ax_qval.clear()
        
        # Create histogram
        plot_range = .1, self.qmax
        hist, edges = np.histogram(self.data[:, 0], bins=2000, range=plotrange)
        self.ax_qval.plot(edges[:-1], hist)
        self.ax_qval.set_title('abcd')
        self.ax_qval.set_xlabel('q-value')
        self.ax_qval.set_ylabel('Count')
        
        # Add tick marks for previously selected q-values
        if self.qvals_1 is not None:
            self.ax_qval.plot(self.qvals_1, np.zeros_like(self.qvals_1),
                             '|', color='blue', markersize=25, markeredgewidth=2, label='qvals_1')
        if self.qvals_2 is not None:
            self.ax_qval.plot(self.qvals_2, np.zeros_like(self.qvals_2),
                             '|', color='red', markersize=25, markeredgewidth=2, label='qvals_2')
        
        # Plot already selected q-values from this session
        if self.selected_qvals:
            self.ax_qval.plot(self.selected_qvals, np.zeros_like(self.selected_qvals),
                             'o', color='purple', markersize=10, label='Selected')
        
        # Redraw span patches
        for patch in self.qval_spans:
            self.ax_qval.add_patch(patch)
        
        if any([self.qvals_1 is not None, self.qvals_2 is not None, self.selected_qvals]):
            self.ax_qval.legend()
        
        self.fig_qval.canvas.draw_idle()

    def _on_qval_span_select(self, xmin, xmax):
        """Handle span selection for q-values."""
        # Select data within the span
        mask = (self.data[:, 0] >= xmin) & (self.data[:, 0] <= xmax)
        selected_data = self.data[mask, 0]
        
        if len(selected_data) > 0:
            # Calculate median of selected range
            median_qval = np.median(selected_data)
            
            # Add to selected q-values
            self.selected_qvals.append(median_qval)
            
            # Create span patch for visualization
            ylims = self.ax_qval.get_ylim()
            span_patch = Rectangle((xmin, ylims[0]), xmax-xmin, ylims[1]-ylims[0],
                                  alpha=0.2, color='purple')
            self.qval_spans.append(span_patch)
            
            # Update plot
            self._plot_qval_histogram(self.ax_qval.get_title())
            
            print(f"Selected q-value: {median_qval:.6f} (from {len(selected_data)} points in range [{xmin:.6f}, {xmax:.6f}])")
            
            # Save to file
            with open('selected_qvals.txt', 'a') as f:
                f.write(f"{median_qval:.6f}\n")

    def _on_qval_clear_last(self, event):
        """Remove the last selected q-value."""
        if self.selected_qvals:
            removed = self.selected_qvals.pop()
            if self.qval_spans:
                self.qval_spans.pop()
            print(f"Removed q-value: {removed:.6f}")
            
            # Update plot
            self._plot_qval_histogram(self.ax_qval.get_title())

    def _on_qval_key_press(self, event):
        """Handle key press events for q-value selection."""
        if event.key == 'c':
            # Clear last selection (same as button)
            self._on_qval_clear_last(None)
        elif event.key == 'd' or event.key == 'enter':
            # Done selecting (same as button)
            self._on_qval_done(None)
        elif event.key == 'escape':
            # Cancel and close without saving
            self.selected_qvals = []
            plt.close(self.fig_qval)

    def _on_qval_done(self, event):
        """Called when Done button is clicked for q-value selection."""
        print(f"Selected {len(self.selected_qvals)} q-values: {self.selected_qvals}")
        plt.close(self.fig_qval)

    def select_triplets(self, title=None):
        """Run the interactive selection process and return the selected triplets."""
        # Create the three windows with specific sizes
        self.fig1, self.ax1 = plt.subplots(num='Step 1: Histogram Selection', figsize=(16,3))
        self.fig2, self.ax2 = plt.subplots(num='Step 2: 2D Selection', figsize=(16,3))
        self.fig3, (self.ax3a, self.ax3b, self.ax3c) = plt.subplots(1, 3, num='Step 3: Final Selection',
                                                                    figsize=(16,3))

        # Set window positions to stack them vertically
        backend = plt.get_backend()
        assert 'Tk' in backend
        manager1 = self.fig1.canvas.manager
        manager2 = self.fig2.canvas.manager
        manager3 = self.fig3.canvas.manager
        dpi = self.fig1.dpi
        height1 = int(3 * dpi)  # 3 inches * dpi

        # Position windows with some spacing
        self.fig1.canvas.manager.window.wm_geometry("+100+50")
        self.fig2.canvas.manager.window.wm_geometry(f"+100+{50 + height1 + 40}")
        self.fig3.canvas.manager.window.wm_geometry(f"+100+{50 + 2*height1 + 80}")

        # Set fixed subplot sizes
        #self.fig3.set_tight_layout(False)
        self.fig1.subplots_adjust(bottom=.2)
        self.fig2.subplots_adjust(bottom=.2)
        self.fig3.subplots_adjust(bottom=.2)

        # Initialize the first window
        self.show_histogram(title=title)

        # Set up the done button
        done_ax = self.fig3.add_axes([0.9, 0.01, 0.09, 0.05])
        self.done_button = plt.Button(done_ax, 'Done')
        self.done_button.on_clicked(self.on_done)

        # Flag to track when selection is complete
        self.selection_done = False

        # Show plots and wait for user interaction
        plt.show(block=True)

        # Close all figures
        plt.close(self.fig1)
        plt.close(self.fig2)
        plt.close(self.fig3)

        # Return the selected triplets
        return self.selected_triplets

    def compute_kde_maxima(self):
        # Normalize data by bandwidths for anisotropic KDE
        normalized_data = self.data / self.bandwidths[np.newaxis, :]

        # Create KernelDensity object
        kde = KernelDensity(bandwidth=1, kernel='cosine')
        print('fit')
        kde.fit(normalized_data)
        print('done fit')

        # Take a random subsample (5%) for evaluation
        n_sample = max(int(len(normalized_data) * 0.05), 100000)  # At least 1000 points
        n_sample = 600000
        n_sample = min(n_sample, len(normalized_data))  # Can't sample more than we have
        print(f'{n_sample=}')

        # Random sampling without replacement
        sample_indices = np.random.choice(len(normalized_data), size=n_sample, replace=False)
        sample_data = normalized_data[sample_indices]

        # Evaluate KDE at sampled points
        print('score')
        sample_densities = np.exp(kde.score_samples(sample_data))
        print('done score')

        # Sort sampled points by density value
        sorted_indices = np.argsort(sample_densities)[::-1]

        # Filter to avoid maxima that are too close together
        safety_factor = 2.0  # How many bandwidths apart maxima should be
        min_distances = self.bandwidths * safety_factor

        # Initialize list to store maxima
        maxima_indices = []
        maxima_points = []
        maxima_densities = []

        # Function to be minimized (negative density)
        def negative_density(x):
            return -np.exp(kde.score_samples([x])[0])

        # Process points in order of decreasing density
        for i, idx in enumerate(sorted_indices):
            if i%1000==0:
                print('processed', i, ', kept', len(maxima_points))
            if len(maxima_points) >= self.n_maxima:
                break

            initial_point = sample_data[idx]
            q1, q2, th = initial_point
            if np.abs(q1-q2) < 2: continue
            if th<8 or th>160: continue

            optimized_point = initial_point
            optimized_density = 1
#            # Run optimization to find true maximum
#            result = minimize(
#                negative_density,
#                initial_point,
#                method='BFGS',
#                #options={'gtol': 1e-5}  # Gradient tolerance for convergence
#            )
#
#            if result.success:
#                optimized_point = result.x
#                optimized_density = -result.fun  # Negate back to get positive density


            # Check if this point is far enough from all accepted maxima
            too_close = False
            for accepted_point in maxima_points:
                dist = np.linalg.norm(optimized_point - accepted_point)
                if dist < 5:
                    too_close = True
                    break

            sus = False
            q1,q2,th = optimized_point
            if np.abs(q1-q2) < 2: sus = True
            if th < 5 or th > 160: sus = True

            if not too_close and not sus:
              maxima_points.append(optimized_point)
              maxima_densities.append(optimized_density)

        # Store the maxima and their density values
        self.kde_maxima = np.array(maxima_points) * self.bandwidths[np.newaxis, :]
        self.kde_values = np.array(maxima_densities)
        with open('autopeaks.txt', 'w') as f:
            for x in self.kde_maxima:
                print(round(x[0], 4), round(x[1], 4), round(x[2], 2), file=f)


        print(f"Found {len(self.kde_maxima)} KDE maxima from {n_sample} sampled points")

    def compute_slice_maxima(self, n_q1_peaks=30, min_q1=0.10, max_q1=0.35,
                             slice_width=0.002, max_peaks_per_slice=100,
                             max_peaks=500, verbose=True):
        """
        Find cluster centers using slice-based 2D KDE approach.

        This method finds q1 peaks using 1D KDE, then performs 2D KDE
        peak finding within each q1 slice, followed by 3D refinement.

        Final peaks are selected by score / sqrt(q1 * q2) to prioritize
        low-q peaks which are most valuable for unit cell determination.

        Parameters
        ----------
        n_q1_peaks : int
            Maximum number of q1 slices to process
        min_q1, max_q1 : float
            Range for q1 peak detection
        slice_width : float
            Width of each q1 slice
        max_peaks_per_slice : int
            Maximum peaks to find per slice
        verbose : bool
            Print progress
        """
        from scipy.ndimage import maximum_filter
        from scipy.optimize import minimize_scalar
        from scipy.spatial import cKDTree

        data = self.data
        bandwidths_3d = self.bandwidths
        bandwidths_2d = np.array([bandwidths_3d[1], bandwidths_3d[2]])

        # Step 1: Find q1 peaks using 1D Gaussian KDE
        if verbose:
            print("Finding q1 peaks...")

        q1_data = data[:, 0]
        q1_mask = (q1_data >= min_q1) & (q1_data <= max_q1)
        q1_filtered = q1_data[q1_mask]

        # Subsample for KDE
        np.random.seed(42)
        n_sample = min(1000000, len(q1_filtered))
        sample_idx = np.random.choice(len(q1_filtered), n_sample, replace=False)
        q1_sample = q1_filtered[sample_idx].reshape(-1, 1)

        # Fit 1D KDE
        kde_1d = KernelDensity(bandwidth=0.0002, kernel='cosine')
        kde_1d.fit(q1_sample)

        # Find peaks on fine grid
        n_eval = 10000
        q1_grid = np.linspace(min_q1, max_q1, n_eval)
        density_1d = np.exp(kde_1d.score_samples(q1_grid.reshape(-1, 1)))

        footprint = np.ones(5)
        local_max = maximum_filter(density_1d, footprint=footprint)
        peaks_mask = (density_1d == local_max) & (density_1d > np.percentile(density_1d, 50))
        peak_q1s_raw = q1_grid[np.where(peaks_mask)[0]]

        # Refine peak positions
        def neg_density(q1):
            return -kde_1d.score_samples([[q1]])[0]

        refined_q1s = []
        for q1 in peak_q1s_raw:
            result = minimize_scalar(neg_density, bounds=(q1-0.0005, q1+0.0005), method='bounded')
            refined_q1s.append(result.x)
        refined_q1s = np.array(refined_q1s)

        # Score and select with spacing
        peak_densities = np.exp(kde_1d.score_samples(refined_q1s.reshape(-1, 1)))
        scores = peak_densities / (refined_q1s ** 0.5)
        sort_idx = np.argsort(scores)[::-1]
        q1s_sorted = refined_q1s[sort_idx]

        min_spacing = 0.002
        q1_centers = []
        for q1 in q1s_sorted:
            if len(q1_centers) >= n_q1_peaks:
                break
            if all(abs(q1 - s) >= min_spacing for s in q1_centers):
                q1_centers.append(q1)
        q1_centers = np.sort(q1_centers)

        if verbose:
            print(f"Found {len(q1_centers)} q1 peaks")

        # Step 2: Process each q1 slice
        all_peaks = []
        all_scores = []

        for i, q1_center in enumerate(q1_centers):
            q1_min = q1_center - slice_width / 2
            q1_max = q1_center + slice_width / 2

            slice_mask = (data[:, 0] >= q1_min) & (data[:, 0] <= q1_max)
            slice_data = data[slice_mask]

            if len(slice_data) < 50:
                continue

            if verbose:
                print(f"Slice {i+1}/{len(q1_centers)}: q1=[{q1_min:.4f}, {q1_max:.4f}], {len(slice_data)} points", end="")

            # Filter diagonal and theta
            diag_threshold = 2.0
            theta_range = (8, 160)
            q_diff_norm = np.abs(slice_data[:, 0] - slice_data[:, 1]) / bandwidths_3d[0]
            off_diag_mask = q_diff_norm >= diag_threshold
            theta_mask = (slice_data[:, 2] >= theta_range[0]) & (slice_data[:, 2] <= theta_range[1])
            slice_data_filtered = slice_data[off_diag_mask & theta_mask]

            if len(slice_data_filtered) < 50:
                if verbose:
                    print(" -> 0 peaks (filtered)")
                continue

            # 3D KDE and tree on full slice
            normalized_3d = slice_data / bandwidths_3d
            kde_3d = KernelDensity(bandwidth=1, kernel='cosine')
            kde_3d.fit(normalized_3d)
            tree_3d = cKDTree(normalized_3d)

            # 2D KDE on filtered subsample
            max_sample_2d = 30000
            if len(slice_data_filtered) > max_sample_2d:
                sample_idx_2d = np.random.choice(len(slice_data_filtered), max_sample_2d, replace=False)
                slice_2d_sample = slice_data_filtered[sample_idx_2d, 1:3]
            else:
                slice_2d_sample = slice_data_filtered[:, 1:3]

            normalized_2d_sample = slice_2d_sample / bandwidths_2d
            kde_2d = KernelDensity(bandwidth=1, kernel='cosine')
            kde_2d.fit(normalized_2d_sample)
            scores_2d = kde_2d.score_samples(normalized_2d_sample)
            tree_2d = cKDTree(normalized_2d_sample)

            # Find 2D peaks
            def refine_2d(point_norm):
                current = point_norm.copy()
                for _ in range(5):
                    idx = tree_2d.query_ball_point(current, r=1.5)
                    if len(idx) < 5:
                        break
                    nearby = normalized_2d_sample[idx]
                    nearby_scores = scores_2d[idx]
                    weights = np.exp(nearby_scores - nearby_scores.max())
                    new_pos = np.average(nearby, axis=0, weights=weights)
                    if np.linalg.norm(new_pos - current) < 0.01:
                        break
                    current = new_pos
                return current

            top_k = min(5000, len(slice_2d_sample))
            top_idx = np.argsort(scores_2d)[-top_k:]

            min_distance_2d = 3.0
            peaks_2d_norm = []
            for idx in np.argsort(scores_2d[top_idx])[::-1]:
                point_idx = top_idx[idx]
                refined = refine_2d(normalized_2d_sample[point_idx])
                is_close = any(np.linalg.norm(refined - ex) < min_distance_2d for ex in peaks_2d_norm)
                if not is_close:
                    peaks_2d_norm.append(refined)
                if len(peaks_2d_norm) >= max_peaks_per_slice:
                    break

            # Refine in 3D
            def refine_3d(point_3d, max_neighbors=300):
                current = point_3d.copy() / bandwidths_3d
                for _ in range(10):
                    idx = tree_3d.query_ball_point(current, r=2.0)
                    if len(idx) < 5:
                        break
                    idx = np.array(idx)
                    if len(idx) > max_neighbors:
                        idx = idx[np.random.choice(len(idx), max_neighbors, replace=False)]
                    nearby = normalized_3d[idx]
                    nearby_scores = kde_3d.score_samples(nearby)
                    weights = np.exp(nearby_scores - nearby_scores.max())
                    new_pos = np.average(nearby, axis=0, weights=weights)
                    if np.linalg.norm(new_pos - current) < 0.01:
                        break
                    current = new_pos
                return current * bandwidths_3d

            min_distance_3d = 3.0
            slice_peaks = []
            slice_scores = []

            for p2d_norm in peaks_2d_norm:
                p2d = p2d_norm * bandwidths_2d
                dq2 = np.abs(slice_data[:, 1] - p2d[0])
                dtheta = np.abs(slice_data[:, 2] - p2d[1])
                nearby_mask = (dq2 < 0.003) & (dtheta < 4.0)

                if nearby_mask.sum() < 5:
                    continue

                nearby_scores_3d = kde_3d.score_samples(normalized_3d[nearby_mask])
                best_nearby = slice_data[nearby_mask][np.argmax(nearby_scores_3d)]

                refined_3d = refine_3d(best_nearby)
                refined_norm = refined_3d / bandwidths_3d

                # Check theta range
                if refined_3d[2] < theta_range[0] or refined_3d[2] > theta_range[1]:
                    continue

                # Check distance to existing peaks
                is_close = any(np.linalg.norm(refined_norm - (ex / bandwidths_3d)) < min_distance_3d
                              for ex in slice_peaks)
                if is_close:
                    continue

                refined_score = kde_3d.score_samples(refined_norm.reshape(1, -1))[0]
                slice_peaks.append(refined_3d)
                slice_scores.append(refined_score)

            all_peaks.extend(slice_peaks)
            all_scores.extend(slice_scores)

            if verbose:
                print(f" -> {len(slice_peaks)} peaks")

        # Deduplicate across slices
        if len(all_peaks) > 0:
            all_peaks = np.array(all_peaks)
            all_scores = np.array(all_scores)

            sort_idx = np.argsort(all_scores)[::-1]
            all_peaks = all_peaks[sort_idx]
            all_scores = all_scores[sort_idx]

            unique_peaks = []
            unique_scores = []
            q_tol, theta_tol = 0.0005, 2.0

            for p, s in zip(all_peaks, all_scores):
                is_dup = any(
                    abs(p[0] - u[0]) < q_tol and
                    abs(p[1] - u[1]) < q_tol and
                    abs(p[2] - u[2]) < theta_tol
                    for u in unique_peaks
                )
                if not is_dup:
                    unique_peaks.append(p)
                    unique_scores.append(s)

            unique_peaks = np.array(unique_peaks)
            unique_scores = np.array(unique_scores)

            # Select top peaks prioritizing low q1*q2 (log-space adjustment)
            if max_peaks is not None and len(unique_peaks) > max_peaks:
                q1_vals = unique_peaks[:, 0]
                q2_vals = unique_peaks[:, 1]
                # score - 0.5*log(q1*q2) is the log-space equivalent of score/sqrt(q1*q2)
                selection_scores = unique_scores - 0.5 * np.log(q1_vals * q2_vals)
                top_idx = np.argsort(selection_scores)[-max_peaks:]
                unique_peaks = unique_peaks[top_idx]
                unique_scores = unique_scores[top_idx]

                if verbose:
                    print(f"Selected top {max_peaks} peaks (prioritizing low-q)")

            self.kde_maxima = unique_peaks
            self.kde_values = unique_scores
        else:
            self.kde_maxima = np.array([]).reshape(0, 3)
            self.kde_values = np.array([])

        if verbose:
            print(f"Final: {len(self.kde_maxima)} maxima")

    def show_histogram(self, title=None):
        self.ax1.clear()
        plotrange=self.qmin, self.qmax
        self.ax1.set_xlim(plotrange)
        hist, edges = np.histogram(self.data[:,0], bins=2000, range=plotrange)
        self.ax1.plot(edges[:-1], hist)
        if title is None:
            title = "Select range in histogram"
        self.ax1.set_title(title)
        self.ax1.set_xlabel('q1')
        self.ax1.set_ylabel('counts')
        
        # Add tick marks for KDE maxima if available
        if hasattr(self, 'kde_maxima'):
            self.ax1.plot(self.kde_maxima[:, 0], np.zeros_like(self.kde_maxima[:, 0]),
                         '|', color='green', markersize=20)

        # Add tick marks for specified q-values
        if self.qvals_1 is not None:
            self.ax1.plot(self.qvals_1, np.zeros_like(self.qvals_1),
                         '|', color='blue', markersize=25, markeredgewidth=2)
        if self.qvals_2 is not None:
            self.ax1.plot(self.qvals_2, np.zeros_like(self.qvals_2),
                         '|', color='red', markersize=25, markeredgewidth=2)

        if self.sb1_qvals is not None:
            self.ax1.plot(self.sb1_qvals, np.zeros_like(self.sb1_qvals),
                         '|', color='orange', markersize=25, markeredgewidth=2)


        self.span = SpanSelector(
            self.ax1,
            self.on_span_select,
            'horizontal',
            useblit=True,
            props=dict(alpha=0.5, facecolor='red')
        )
        
        # Redraw the figure
        self.fig1.canvas.draw_idle()

    def on_span_select(self, xmin, xmax):
        # Remove previous span patch if it exists
        if self.span_patch is not None:
            try:
                self.span_patch.remove()
            except NotImplementedError:
                pass
            self.span_patch = None
        
        # Create new span patch
        ylims = self.ax1.get_ylim()
        self.span_patch = Rectangle((xmin, ylims[0]), xmax-xmin, ylims[1]-ylims[0],
                                  alpha=0.2, color='red')
        self.ax1.add_patch(self.span_patch)
        
        # Select data within the span
        mask = (self.data[:, 0] >= xmin) & (self.data[:, 0] <= xmax)
        self.current_selection = self.data[mask]
        
        # Select KDE maxima within the span
        if hasattr(self, 'kde_maxima'):
            kde_mask = (self.kde_maxima[:, 0] >= xmin) & (self.kde_maxima[:, 0] <= xmax)
            self.current_kde_selection = self.kde_maxima[kde_mask]

        # Select predicted points in the span
        if self.points is not None:
            points_mask = (self.points[:,0] >= xmin) & (self.points[:,0] <= xmax)
            self.current_points_selection = self.points[points_mask]
        else:
            self.current_points_selection = None
        
        # Show the second window
        self.show_scatter_2d()
        
        # Redraw both figures
        self.fig1.canvas.draw_idle()
        self.fig2.canvas.draw_idle()

    def show_scatter_2d(self):
        self.ax2.clear()
        plotrange = self.qmin, self.qmax
        self.ax2.set_xlim(plotrange)
        self.ax2.set_xlabel('q2')
        self.ax2.set_ylabel('theta')
        if self.current_selection is not None and len(self.current_selection) > 0:
            self.ax2.scatter(
                self.current_selection[:, 1],
                self.current_selection[:, 2],
                s=2, alpha=.2)
            self.ax2.set_title(f'Draw box to select points (n={len(self.current_selection)})')

            # Plot KDE maxima if available
            if hasattr(self, 'current_kde_selection') and len(self.current_kde_selection) > 0:
                self.ax2.scatter(self.current_kde_selection[:, 1], self.current_kde_selection[:, 2],
                               color='red', marker='x', s=50, label='KDE maxima')

            
            # Add vertical lines for q-values (second dimension)
            if self.qvals_2 is not None:
                ylim = self.ax2.get_ylim()
                for q in self.qvals_2:
                    if q >= self.current_selection[:, 1].min() and q <= self.current_selection[:, 1].max():
                        self.ax2.axvline(q, color='blue', linestyle='--', alpha=0.5)
            if self.qvals_1 is not None and self.points is None:
                for q in self.qvals_1:
                    if q >= self.current_selection[:, 1].min() and q <= self.current_selection[:, 1].max():
                        self.ax2.axvline(q, color='red', linestyle='--', alpha=0.5)
            if self.current_points_selection is not None:
                self.ax2.scatter(self.current_points_selection[:,1], self.current_points_selection[:,2],
                               color='orange', marker='o', s=20, label='points')

            # Create RectangleSelector only if it doesn't exist already
            if not hasattr(self, 'rect') or self.rect is None:
                self.rect = RectangleSelector(
                    self.ax2,
                    self.on_rect_select,
                    useblit=True,
                    props=dict(facecolor='red', alpha=0.2)
                )

            # If there was a previous rectangle, redraw it
            if self.rect_patch is not None:
                self.ax2.add_patch(self.rect_patch)
        else:
            self.ax2.set_title('No points selected')

        self.fig2.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig2.canvas.draw_idle()

    def on_rect_select(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Remove previous rectangle if it exists
        if self.rect_patch is not None:
            self.rect_patch.remove()

        # Create new rectangle patch
        self.rect_patch = Rectangle((min(x1, x2), min(y1, y2)),
                                  abs(x2-x1), abs(y2-y1),
                                  alpha=0.2, color='red')
        self.ax2.add_patch(self.rect_patch)
        self.fig2.canvas.draw_idle()

        # Select data within the rectangle
        mask = (
            (self.current_selection[:, 1] >= min(x1, x2)) &
            (self.current_selection[:, 1] <= max(x1, x2)) &
            (self.current_selection[:, 2] >= min(y1, y2)) &
            (self.current_selection[:, 2] <= max(y1, y2))
        )
        self.final_selection = self.current_selection[mask]

        # Select KDE maxima within the rectangle
        if hasattr(self, 'current_kde_selection'):
            kde_mask = (
                (self.current_kde_selection[:, 1] >= min(x1, x2)) &
                (self.current_kde_selection[:, 1] <= max(x1, x2)) &
                (self.current_kde_selection[:, 2] >= min(y1, y2)) &
                (self.current_kde_selection[:, 2] <= max(y1, y2))
            )
            self.final_kde_selection = self.current_kde_selection[kde_mask]

        # Show the third window
        self.show_final_plots()

    def show_final_plots(self):
        self.update_means()

        # Clear all axes
        for ax in [self.ax3a, self.ax3b, self.ax3c]:
            ax.clear()
            #ax.set_aspect('equal', adjustable='datalim')

        # Find the limits that contain all selected points
        x_min, x_max = self.final_selection[:, 0].min(), self.final_selection[:, 0].max()
        y_min, y_max = self.final_selection[:, 1].min(), self.final_selection[:, 1].max()
        z_min, z_max = self.final_selection[:, 2].min(), self.final_selection[:, 2].max()

        # Add a small margin
        margin = 0.1
        x_range = x_max - x_min
        y_range = y_max - y_min
        z_range = z_max - z_min

#        x_min -= margin * x_range
#        x_max += margin * x_range
#        y_min -= margin * y_range
#        y_max += margin * y_range
#        z_min -= margin * z_range
#        z_max += margin * z_range

        # Set the limits for each plot
        self.ax3a.set_xlim(x_min, x_max)
        self.ax3a.set_ylim(y_min, y_max)
        self.ax3a.set_xlabel('q1')
        self.ax3a.set_ylabel('q2')
        self.ax3b.set_xlabel('q1')
        self.ax3b.set_ylabel('theta')
        self.ax3c.set_xlabel('q2')
        self.ax3c.set_ylabel('theta')

        self.ax3b.set_xlim(x_min, x_max)
        self.ax3b.set_ylim(z_min, z_max)

        self.ax3c.set_xlim(y_min, y_max)
        self.ax3c.set_ylim(z_min, z_max)

        # Plot all points within the final selection
        self.ax3a.scatter(self.final_selection[:, 0], self.final_selection[:, 1], color='tab:blue',
                          s=5, alpha=.3)
        self.ax3b.scatter(self.final_selection[:, 0], self.final_selection[:, 2], color='tab:blue',
                          s=5, alpha=.3)
        self.ax3c.scatter(self.final_selection[:, 1], self.final_selection[:, 2], color='tab:blue',
                          s=5, alpha=.3)

        # Plot KDE maxima
        if hasattr(self, 'final_kde_selection') and len(self.final_kde_selection) > 0:
            self.ax3a.scatter(self.final_kde_selection[:, 0], self.final_kde_selection[:, 1], 
                            color='red', marker='x', s=50)
            self.ax3b.scatter(self.final_kde_selection[:, 0], self.final_kde_selection[:, 2], 
                            color='red', marker='x', s=50)
            self.ax3c.scatter(self.final_kde_selection[:, 1], self.final_kde_selection[:, 2], 
                            color='red', marker='x', s=50)
        
        # Plot means
        self.ax3a.plot(self.means[0], self.means[1], 'ro', markersize=6)
        self.ax3b.plot(self.means[0], self.means[2], 'ro', markersize=6)
        self.ax3c.plot(self.means[1], self.means[2], 'ro', markersize=6)

        # Add rectangle selectors to all plots
        self.rect_final = [
            RectangleSelector(ax, self.on_final_rect_select, useblit=True,
                            props=dict(facecolor='red', alpha=0.2))
            for ax in [self.ax3a, self.ax3b, self.ax3c]
        ]

        self.fig3.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig3.canvas.draw_idle()



    def on_final_rect_select(self, eclick, erelease):
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        # Get the current axis
        ax = eclick.inaxes

        # Create mask based on which plot was clicked
        if ax == self.ax3a:
            mask = (
                (self.final_selection[:, 0] >= min(x1, x2)) &
                (self.final_selection[:, 0] <= max(x1, x2)) &
                (self.final_selection[:, 1] >= min(y1, y2)) &
                (self.final_selection[:, 1] <= max(y1, y2))
            )
        elif ax == self.ax3b:
            mask = (
                (self.final_selection[:, 0] >= min(x1, x2)) &
                (self.final_selection[:, 0] <= max(x1, x2)) &
                (self.final_selection[:, 2] >= min(y1, y2)) &
                (self.final_selection[:, 2] <= max(y1, y2))
            )
        elif ax == self.ax3c:
            mask = (
                (self.final_selection[:, 1] >= min(x1, x2)) &
                (self.final_selection[:, 1] <= max(x1, x2)) &
                (self.final_selection[:, 2] >= min(y1, y2)) &
                (self.final_selection[:, 2] <= max(y1, y2))
            )

        # Save the limits
        axl = self.ax3a.get_xlim()
        ayl = self.ax3a.get_ylim()
        bxl = self.ax3b.get_xlim()
        byl = self.ax3b.get_ylim()
        cxl = self.ax3c.get_xlim()
        cyl = self.ax3c.get_ylim()
        # Clear all plots including mean markers
        for ax in [self.ax3a, self.ax3b, self.ax3c]:
            ax.clear()
            #ax.set_aspect('equal', adjustable='datalim')
            # Restore the original limits
            if ax == self.ax3a:
                ax.set_xlim(axl)
                ax.set_ylim(ayl)
            elif ax == self.ax3b:
                ax.set_xlim(bxl)
                ax.set_ylim(byl)
            else:
                ax.set_xlim(cxl)
                ax.set_ylim(cyl)

        # Update means based on the new selection
        selected_points = self.final_selection[mask]
        self.means = np.mean(selected_points, axis=0)

        # Replot everything
        self.ax3a.scatter(self.final_selection[~mask, 0], self.final_selection[~mask, 1],
                          color='gray', s=5, alpha=0.3)
        self.ax3a.scatter(self.final_selection[mask, 0], self.final_selection[mask, 1],
                          s=5, alpha=.3, color='tab:blue')
        self.ax3a.plot(self.means[0], self.means[1], 'ro', markersize=6)

        self.ax3b.scatter(self.final_selection[~mask, 0], self.final_selection[~mask, 2],
                          color='gray', s=5, alpha=0.3)
        self.ax3b.scatter(self.final_selection[mask, 0], self.final_selection[mask, 2],
                          s=5, alpha=.3, color='tab:blue')
        self.ax3b.plot(self.means[0], self.means[2], 'ro', markersize=6)

        self.ax3c.scatter(self.final_selection[~mask, 1], self.final_selection[~mask, 2],
                          color='gray', s=5, alpha=0.3)
        self.ax3c.scatter(self.final_selection[mask, 1], self.final_selection[mask, 2],
                          s=5, alpha=.3, color='tab:blue')
        self.ax3c.plot(self.means[1], self.means[2], 'ro', markersize=6)

        # Plot KDE maxima
        if hasattr(self, 'final_kde_selection') and len(self.final_kde_selection) > 0:
            self.ax3a.scatter(self.final_kde_selection[:, 0], self.final_kde_selection[:, 1], 
                            color='red', marker='x', s=50)
            self.ax3b.scatter(self.final_kde_selection[:, 0], self.final_kde_selection[:, 2], 
                            color='red', marker='x', s=50)
            self.ax3c.scatter(self.final_kde_selection[:, 1], self.final_kde_selection[:, 2], 
                            color='red', marker='x', s=50)
        
        # Plot means
        self.ax3a.plot(self.means[0], self.means[1], 'ro', markersize=6)
        self.ax3b.plot(self.means[0], self.means[2], 'ro', markersize=6)
        self.ax3c.plot(self.means[1], self.means[2], 'ro', markersize=6)

        self.ax3a.set_xlabel('q1')
        self.ax3a.set_ylabel('q2')
        self.ax3b.set_xlabel('q1')
        self.ax3b.set_ylabel('theta')
        self.ax3c.set_xlabel('q2')
        self.ax3c.set_ylabel('theta')

        self.fig3.canvas.draw_idle()

    def update_means(self):
        self.means = np.mean(self.final_selection, axis=0)

    def on_key_press(self, event):
        if event.key == 'a':
            # Add mean to selected triplets
            triplet = [self.means[0], self.means[1], self.means[2]]
            self.selected_triplets.append(triplet)

            print(f"Selected point {self.means[0]:.6f} {self.means[1]:.6f} {self.means[2]:.6f}. "
                  f"Total triplets: {len(self.selected_triplets)}")

            # Also append to file if desired
            with open('cluster_means.txt', 'a') as f:
                np.savetxt(f, [self.means], fmt='%.6f')

            if self.sb1_callback is not None:
                self.sb1_qvals, self.points = self.sb1_callback(triplet, self.recon)
                self.show_histogram(title=self.ax1.get_title())

    def on_done(self, event):
        """Called when the Done button is clicked."""
        self.selection_done = True
        plt.close('all')  # Close all open figures
        import gc;gc.collect()


# Example usage:
if __name__ == "__main__":
    # Janky format
    alldata = []
    for f in sys.argv[1:]:
      data = np.load(f)['triplets'][:,1:4]
      data[:,0] = 1/data[:,0]
      data[:,1] = 1/data[:,1]
      data2 = np.vstack((data[:,1], data[:,0], data[:,2])).transpose()
      data = np.vstack((data, data2))
      alldata.append(data)
    alldata = np.vstack(alldata)
    #alldata = alldata[:100000]

    bandwidths = [0.001, 0.001, 1.0]
    
    clusterer = ManualClusterer(alldata, bandwidths=bandwidths, n_maxima=500)
    clusterer.select_triplets()

