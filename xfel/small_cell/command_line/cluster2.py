import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde
from scipy.ndimage import maximum_filter
from scipy.ndimage import generate_binary_structure
from sklearn.neighbors import KernelDensity
import sys


counter = 0

class ManualClusterer:
    def __init__(self, data, bandwidths=[0.001, 0.001, 1.0], n_maxima=500):
        self.data = data
        self.current_selection = None
        self.final_selection = None
        self.means = None
        self.span_patch = None
        self.rect_patch = None
        self.bandwidths = np.array(bandwidths)
        self.n_maxima = n_maxima
        
        # Compute KDE maxima
        self.compute_kde_maxima()
        
        # Create the three windows with specific sizes
        self.fig1, self.ax1 = plt.subplots(num='Step 1: Histogram Selection', figsize=(16, 3))
        self.fig2, self.ax2 = plt.subplots(num='Step 2: 2D Selection', figsize=(16, 3))
        self.fig3, (self.ax3a, self.ax3b, self.ax3c) = plt.subplots(1, 3, num='Step 3: Final Selection', 
                                                                    figsize=(16, 3))
        
        # Set fixed subplot sizes
        self.fig3.set_tight_layout(False)
#        for i, ax in enumerate([self.ax3a, self.ax3b, self.ax3c]):
#            ax.set_aspect('equal', adjustable='datalim')
#            # Set a fixed position and size for each subplot
#            box = ax.get_position()
#            ax.set_position([0.1 + i*0.3, 0.1, 0.25, 0.25])  # [left, bottom, width, height]
        
        # Initialize the first window
        self.show_histogram()
        
#        # Arrange windows vertically
#        self.fig1.canvas.manager.window.move(0, 0)
#        self.fig2.canvas.manager.window.move(0, 300)
#        self.fig3.canvas.manager.window.move(0, 600)
        
        plt.show()


    def compute_kde_maxima(self):
        # Normalize data by bandwidths for anisotropic KDE
        normalized_data = self.data / self.bandwidths[np.newaxis, :]

        # Create KernelDensity object
        kde = KernelDensity(bandwidth=1.0, kernel='epanechnikov')
        print('fit')
        kde.fit(normalized_data)
        print('done fit')

        # Take a random subsample (5%) for evaluation
        n_sample = max(int(len(normalized_data) * 0.05), 1000)  # At least 1000 points
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

        # Process points in order of decreasing density
        for idx in sorted_indices:
            if len(maxima_indices) >= self.n_maxima:
                break

            point = sample_data[idx]

            # Check if this point is far enough from all accepted maxima
            too_close = False
            for accepted_idx in maxima_indices:
                accepted_point = sample_data[accepted_idx]
                dist = np.linalg.norm(point - accepted_point)
                if dist < 5:
                    too_close = True
                    break

            sus = False
            q1,q2,th = point
            if np.abs(q1-q2) < 2: sus = True
            if th < 5 or th > 160: sus = True

            if not too_close and not sus:
                maxima_indices.append(idx)

        # Store the maxima and their density values
        self.kde_maxima = sample_data[maxima_indices] * self.bandwidths[np.newaxis, :]
        self.kde_values = sample_densities[maxima_indices]
        for x in self.kde_maxima[:10]:
            print(round(x[0], 4), round(x[1], 4), round(x[2], 2))


        print(f"Found {len(self.kde_maxima)} KDE maxima from {n_sample} sampled points")
#    def compute_kde_maxima(self):
#        # Normalize data by bandwidths for anisotropic KDE
#        normalized_data = self.data / self.bandwidths[np.newaxis, :]
#
#        # Create KernelDensity object
#        kde = KernelDensity(bandwidth=1.0, kernel='epanechnikov')
#        print('fit')
#        kde.fit(normalized_data)
#        print('done fit')
#
#        # Create grid for searching maxima
#        grid_points = 50  # balance between accuracy and performance
#        grid_ranges = []
#        for dim in range(3):
#            min_val, max_val = normalized_data[:, dim].min(), normalized_data[:, dim].max()
#            grid_ranges.append(np.linspace(min_val, max_val, grid_points))
#
#        # Create meshgrid
#        X, Y, Z = np.meshgrid(*grid_ranges, indexing='ij')
#        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
#        positions = self.data
#
#        #import IPython;IPython.embed()
#        # Evaluate KDE
#        print('score')
#        densities = np.exp(kde.score_samples(positions)).reshape(X.shape)
#        print('done score')
#
#        # Find local maxima
#        neighborhood = generate_binary_structure(3, 2)  # 3D, connectivity 2
#        local_max = maximum_filter(densities, footprint=neighborhood) == densities
#
#        # Get coordinates and values of local maxima
#        local_max_indices = np.where(local_max)
#        local_max_values = densities[local_max]
#
#        # Sort by density value
#        sorted_indices = np.argsort(local_max_values)[::-1][:self.n_maxima]
#
#        # Extract coordinates
#        maxima_coords = np.zeros((len(sorted_indices), 3))
#        for i, idx in enumerate(sorted_indices):
#            maxima_coords[i, 0] = grid_ranges[0][local_max_indices[0][idx]]
#            maxima_coords[i, 1] = grid_ranges[1][local_max_indices[1][idx]]
#            maxima_coords[i, 2] = grid_ranges[2][local_max_indices[2][idx]]
#
#        # Scale back to original coordinates
#        self.kde_maxima = maxima_coords * self.bandwidths[np.newaxis, :]
#        self.kde_values = local_max_values[sorted_indices]
#
#        print(f"Found {len(self.kde_maxima)} KDE maxima")

    def show_histogram(self):
        self.ax1.clear()
        hist, edges = np.histogram(self.data[:,0], bins=2000, range=(.1,.5))
        self.ax1.plot(edges[:-1], hist)
        self.ax1.set_title('Select range in histogram')
        
        # Add tick marks for KDE maxima
        self.ax1.plot(self.kde_maxima[:, 0], np.zeros_like(self.kde_maxima[:, 0]), 
                     '|', color='red', markersize=20)

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
            self.span_patch.remove()
        
        # Create new span patch
        ylims = self.ax1.get_ylim()
        self.span_patch = Rectangle((xmin, ylims[0]), xmax-xmin, ylims[1]-ylims[0],
                                  alpha=0.2, color='red')
        self.ax1.add_patch(self.span_patch)
        
        # Select data within the span
        mask = (self.data[:, 0] >= xmin) & (self.data[:, 0] <= xmax)
        self.current_selection = self.data[mask]
        
        # Select KDE maxima within the span
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

            # Plot KDE maxima
            if hasattr(self, 'current_kde_selection') and len(self.current_kde_selection) > 0:
                self.ax2.scatter(self.current_kde_selection[:, 1], self.current_kde_selection[:, 2], 
                               color='red', marker='x', s=50, label='KDE maxima')
            
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
        global counter
        # note counter is a global variable
        if event.key == 'a':
            # Append means to file
            counter += 1
            with open('cluster_means.txt', 'a') as f:
                np.savetxt(f, [self.means], fmt='%.6f')
            
            print(f"Selected point {self.means[0]:.6f} {self.means[1]:.6f} {self.means[2]:.6f}. "
                  f"{counter} lines have been written.")

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

