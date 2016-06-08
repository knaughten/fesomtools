from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon

# Calculate and plot the ice shelf draft for the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = filename for figure
def plot_shelf (mesh_path, fig_name):

    # Plotting parameters
    circumpolar = True
    lat_max = -30 + 90
    font_sizes = [240, 192, 160]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Read the ice shelf draft for each node
    file = open(mesh_path + 'shelf.out', 'r')
    node_shelf = []
    for line in file:
        node_shelf.append(-float(line))
    file.close()

    # For each element, calculate the average shelf of the 3 component nodes
    elm_shelf = []
    for elm in elements:
        shelf1 = node_shelf[(elm.nodes[0]).id-1]
        shelf2 = node_shelf[(elm.nodes[1]).id-1]
        shelf3 = node_shelf[(elm.nodes[2]).id-1]
        elm_shelf.append(mean(array([shelf1, shelf2, shelf3])))

    # Set up figure
    fig = figure(figsize=(128, 96))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(elm_shelf))
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
    img.set_clim(vmin=min(elm_shelf), vmax=max(elm_shelf))
    axis('off')

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    fig_name = raw_input("Filename for figure: ")
    plot_shelf(mesh_path, fig_name)
