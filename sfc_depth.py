from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon


# Calculate and plot the surface depth (i.e. to show if ice shelf
# cavities are working) for the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = path to desired output file for figure
def sfc_depth (mesh_path, fig_name):

    # Plotting parameters
    circumpolar = True
    lat_max = -60+90
    font_sizes = [18, 16,12]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Calculate the surface depth for each element
    values = []
    for elm in elements:
        values.append(mean(array([elm.nodes[0].depth, elm.nodes[1].depth, elm.nodes[2].depth])))

    # Set up figure
    fig = figure(figsize=(16, 12))
    ax = fig.add_subplot(1,1,1, aspect='equal')
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    xlim([-lat_max, lat_max])
    ylim([-lat_max, lat_max])
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    title('Depth of surface elements (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    fig_name = raw_input("File name for figure: ")
    sfc_depth(mesh_path, fig_name)
