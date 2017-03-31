from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon


# Calculate and plot the horizontal grid resolution (square root of the area 
# of each triangular 2D element) for the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# circumpolar = optional boolean (default True) indicating whether to plot a 
#               circumpolar Antarctic view or the global Cartesian view
# set_bound = optional boolean (default False) indicating whether to place a 
#             maximum on the colour scale
# bound = if set_bound is True, maximum bound for the colour scale
# save = optional boolean (default False) indicating whether to save the figure
#        to a file, or display it on the screen
# fig_name = if save is True, filename for the figure
def grid_res (mesh_path, circumpolar=True, set_bound=False, bound=None, save=False, fig_name=None):

    # Plotting parameters
    if circumpolar:
        lat_max = -30 + 90
    font_sizes = [30, 24, 20]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Calculate the grid resolution for each element
    values = []
    for elm in elements:
        values.append(sqrt(elm.area())*1e-3)

    # Set up figure    
    if circumpolar:
        fig = figure(figsize=(16,12))
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        fig = figure(figsize=(16, 8))
        ax = fig.add_subplot(1,1,1)
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(values))
    img.set_edgecolor('face')
    ax.add_collection(img)

    # Configure plot
    if circumpolar:
        xlim([-lat_max, lat_max])
        ylim([-lat_max, lat_max])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        axis('off')
    else:
        xlim([-180, 180])
        ylim([-90, 90])
        ax.get_xaxis().set_ticks(arange(-120,120+1,60))
        ax.get_yaxis().set_ticks(arange(-60,60+1,30))
    title('Horizontal grid resolution (km)', fontsize=font_sizes[0])
    if set_bound:
        cbar = colorbar(img, extend='max')
    else:
        cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    if set_bound:
        img.set_clim(vmin=0, vmax=bound)

    if save:
        fig.savefig(fig_name)
    else:
        fig.show()


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    domain = raw_input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    get_bound = raw_input("Set maximum bound for colour scale? (y/n) ")
    if get_bound == 'y':
        set_bound = True
        bound = float(raw_input("Maximum bound (km): "))
    else:
        set_bound = False
        bound = None
    action = raw_input("Save figure (s) or display on screen (d)? ")
    if action == 's':
        save = True
        fig_name = raw_input("Filename for figure: ")
    else:
        save = False
        fig_name = None
    grid_res(mesh_path, circumpolar, set_bound, bound, save, fig_name)
