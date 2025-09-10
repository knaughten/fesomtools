from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

# Calculate and plot the number of layers in the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = filename for figure
def plot_num_layers (mesh_path, fig_name):

    # Plotting parameters
    circumpolar = True
    lat_max = -30 + 90
    font_sizes = [240, 192, 160]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Calculate the number of layers for each element
    num_layers = []
    for elm in elements:
        num_layers_elm = 0
        # Count the number of layers for each node
        for i in range(3):
            node = elm.nodes[i]
            num_layers_node = 1
            # Iterate until we reach the bottom
            while node.below is not None:
                num_layers_node += 1
                node = node.below
            num_layers_elm = max(num_layers_elm, num_layers_node)
        # Save the maximum number of layers across the 3 nodes
        num_layers.append(num_layers_elm)

    # Set up figure
    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(num_layers))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    title('Ice shelf draft (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=1, vmax=47)
    axis('off')

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to mesh directory: ")
    fig_name = input("Filename for figure: ")
    plot_num_layers(mesh_path, fig_name)
