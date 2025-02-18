import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, RectangleSelector
from matplotlib.patches import Rectangle
import sys


counter = 0

class ManualClusterer:
    def __init__(self, data):
        self.data = data
        self.current_selection = None
        self.final_selection = None
        self.means = None
        self.span_patch = None
        self.rect_patch = None
        
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

    def show_histogram(self):
        self.ax1.clear()
        hist, edges = np.histogram(self.data[:,0], bins=2000, range=(.1,.5))
        self.ax1.plot(edges[:-1], hist)
        self.ax1.set_title('Select range in histogram')
        
        self.span = SpanSelector(
            self.ax1,
            self.on_span_select,
            'horizontal',
            useblit=True,
            rectprops=dict(alpha=0.5, facecolor='red')
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
            
            self.rect = RectangleSelector(
                self.ax2,
                self.on_rect_select,
                useblit=True,
                rectprops=dict(facecolor='red', alpha=0.2)
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

        # Select data within the rectangle
        mask = (
            (self.current_selection[:, 1] >= min(x1, x2)) &
            (self.current_selection[:, 1] <= max(x1, x2)) &
            (self.current_selection[:, 2] >= min(y1, y2)) &
            (self.current_selection[:, 2] <= max(y1, y2))
        )
        self.final_selection = self.current_selection[mask]

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

        x_min -= margin * x_range
        x_max += margin * x_range
        y_min -= margin * y_range
        y_max += margin * y_range
        z_min -= margin * z_range
        z_max += margin * z_range

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

        # Plot means
        self.ax3a.plot(self.means[0], self.means[1], 'ro', markersize=6)
        self.ax3b.plot(self.means[0], self.means[2], 'ro', markersize=6)
        self.ax3c.plot(self.means[1], self.means[2], 'ro', markersize=6)

        # Add rectangle selectors to all plots
        self.rect_final = [
            RectangleSelector(ax, self.on_final_rect_select, useblit=True,
                            rectprops=dict(facecolor='red', alpha=0.2))
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
    # Generate some sample data
    # Janky format
    data = np.load(sys.argv[1])['triplets'][:,1:4]
    data[:,0] = 1/data[:,0]
    data[:,1] = 1/data[:,1]
    data2 = np.vstack((data[:,1], data[:,0], data[:,2])).transpose()
    data = np.vstack((data, data2))

    
    clusterer = ManualClusterer(data)
