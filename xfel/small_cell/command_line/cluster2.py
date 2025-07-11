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


counter = 0

class ManualClusterer:
    def __init__(self, data, bandwidths=[0.001, 0.001, 1.0], n_maxima=500,
                 qvals_1=None, qvals_2=None, sb1_callback=None, recon=None):
        self.data = data
        self.bandwidths = np.array(bandwidths)
        self.n_maxima = n_maxima
        self.qvals_1 = qvals_1
        self.qvals_2 = qvals_2
        self.sb1_qvals = None

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
            self.compute_kde_maxima()
        
#        # Create the three windows with specific sizes
#        self.fig1, self.ax1 = plt.subplots(num='Step 1: Histogram Selection', figsize=(16, 3))
#        self.fig2, self.ax2 = plt.subplots(num='Step 2: 2D Selection', figsize=(16, 3))
#        self.fig3, (self.ax3a, self.ax3b, self.ax3c) = plt.subplots(1, 3, num='Step 3: Final Selection', 
#                                                                    figsize=(16, 3))
#
#        # Set fixed subplot sizes
#        self.fig3.set_tight_layout(False)
##        for i, ax in enumerate([self.ax3a, self.ax3b, self.ax3c]):
##            ax.set_aspect('equal', adjustable='datalim')
##            # Set a fixed position and size for each subplot
##            box = ax.get_position()
##            ax.set_position([0.1 + i*0.3, 0.1, 0.25, 0.25])  # [left, bottom, width, height]
#
#        # Initialize the first window
#        self.show_histogram()
#
##        # Arrange windows vertically
##        self.fig1.canvas.manager.window.move(0, 0)
##        self.fig2.canvas.manager.window.move(0, 300)
##        self.fig3.canvas.manager.window.move(0, 600)
#
#        plt.show()

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
        hist, edges = np.histogram(self.data[:, 0], bins=2000, range=(.1, .5))
        self.ax_qval.plot(edges[:-1], hist)
        self.ax_qval.set_title(title)
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
        self.fig1, self.ax1 = plt.subplots(num='Step 1: Histogram Selection', figsize=(16,2.5))
        self.fig2, self.ax2 = plt.subplots(num='Step 2: 2D Selection', figsize=(16,2.5))
        self.fig3, (self.ax3a, self.ax3b, self.ax3c) = plt.subplots(1, 3, num='Step 3: Final Selection',
                                                                    figsize=(16,2.5))

        # Set window positions to stack them vertically
        backend = plt.get_backend()
        assert 'Tk' in backend
        manager1 = self.fig1.canvas.manager
        manager2 = self.fig2.canvas.manager
        manager3 = self.fig3.canvas.manager
        dpi = self.fig1.dpi
        height1 = int(2.5 * dpi)  # 3 inches * dpi

        # Position windows with some spacing
        self.fig1.canvas.manager.window.wm_geometry("+100+50")
        self.fig2.canvas.manager.window.wm_geometry(f"+100+{50 + height1 + 40}")
        self.fig3.canvas.manager.window.wm_geometry(f"+100+{50 + 2*height1 + 80}")

        # Set fixed subplot sizes
        self.fig3.set_tight_layout(False)

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
        n_sample = min(n_sample, len(normalized_data))  # Can't sample more than we have

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

    def show_histogram(self, title=None):
        self.ax1.clear()
        hist, edges = np.histogram(self.data[:,0], bins=2000, range=(.1,.5))
        self.ax1.plot(edges[:-1], hist)
        if title is None:
            title = "Select range in histogram"
        self.ax1.set_title(title)
        
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
                         '|', color='green', markersize=25, markeredgewidth=2)


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
        
        # Show the second window
        self.show_scatter_2d()
        
        # Redraw both figures
        self.fig1.canvas.draw_idle()
        self.fig2.canvas.draw_idle()

    def show_scatter_2d(self):
        self.ax2.clear()
        self.ax2.set_xlim((.1, .5))
        if self.current_selection is not None and len(self.current_selection) > 0:
            self.ax2.scatter(
                self.current_selection[:, 1],
                self.current_selection[:, 2],
                s=2, alpha=.2)
            self.ax2.set_title(f'Draw box to select points (n={len(self.current_selection)})')

            # Plot KDE maxima if available
            if hasattr(self, 'current_kde_selection') and len(self.current_kde_selection) > 0:
                self.ax2.scatter(self.current_kde_selection[:, 1], self.current_kde_selection[:, 2],
                               color='green', marker='x', s=50, label='KDE maxima')

            
            # Add vertical lines for q-values (second dimension)
            if self.qvals_2 is not None:
                ylim = self.ax2.get_ylim()
                for q in self.qvals_2:
                    if q >= self.current_selection[:, 1].min() and q <= self.current_selection[:, 1].max():
                        self.ax2.axvline(q, color='blue', linestyle='--', alpha=0.5)
            if self.qvals_1 is not None:
                for q in self.qvals_1:
                    if q >= self.current_selection[:, 1].min() and q <= self.current_selection[:, 1].max():
                        self.ax2.axvline(q, color='red', linestyle='--', alpha=0.5)

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
        self.ax3a.plot(self.means[0], self.means[1], 'r*', markersize=15)
        self.ax3b.plot(self.means[0], self.means[2], 'r*', markersize=15)
        self.ax3c.plot(self.means[1], self.means[2], 'r*', markersize=15)
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
                self.sb1_qvals = self.sb1_callback(triplet, self.recon)
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

