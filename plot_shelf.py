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
    font_sizes = [24, 20, 16]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Find the ice shelf draft at each element
    # (i.e. depth of surface nodes in ice shelf cavities)
    elm_shelf = []
    for elm in elements:
        if elm.cavity:
            shelf1 = (elm.nodes[0]).depth
            shelf2 = (elm.nodes[1]).depth
            shelf3 = (elm.nodes[2]).depth
            elm_shelf.append(mean(array([shelf1, shelf2, shelf3])))
        else:
            elm_shelf.append(0.0)

    # Set up figure
    fig = figure(figsize=(12, 9))
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
    axis('off')

    show()
    #savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = input("Path to mesh directory: ")
    fig_name = input("Filename for figure: ")
    plot_shelf(mesh_path, fig_name)
