from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon


# Calculate and plot the horizontal grid resolution (square root of the area 
# of each triangular 2D element) for the given FESOM grid, on a regular
# Cartesian grid (not circumpolar).
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = path to desired output file for figure
def grid_res_cartesian (mesh_path, fig_name):

    # Plotting parameters
    circumpolar = False
    lat_max = -60+90
    font_sizes = [30, 24, 20]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Calculate the grid resolution for each element
    values = []
    for elm in elements:
        values.append(sqrt(elm.area())*1e-3)

    # Set up figure
    fig = figure(figsize=(16, 8))
    ax = fig.add_subplot(1,1,1)
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    xlim([-180, 180])
    ylim([-90, 90])
    ax.get_xaxis().set_ticks(arange(-120,120+1,60))
    ax.get_yaxis().set_ticks(arange(-60,60+1,30))
    title('Horizontal grid resolution (km)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=0, vmax=max(values))

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    fig_name = raw_input("File name for figure: ")
    grid_res_cartesian(mesh_path, fig_name)
