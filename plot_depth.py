from patches import *
from netCDF4 import Dataset
from numpy import *
from matplotlib.collections import PatchCollection
from matplotlib.pyplot import *
from matplotlib.cm import *
from matplotlib.patches import Polygon


# Calculate and plot the seafloor depth for the given FESOM grid.
# Input:
# mesh_path = path to FESOM mesh directory
# fig_name = filename for figure
# circumpolar = boolean flag indicating whether to plot a circumpolar Antarctic
#               view or the global Cartesian view
def plot_depth (mesh_path, fig_name, circumpolar=True):

    # Plotting parameters
    if circumpolar:
        lat_max = -30 + 90
        font_sizes = [240, 192, 160]
    else:
        font_sizes = [30, 24, 20]

    # Build triangular patches for each element
    elements, patches = make_patches(mesh_path, circumpolar)

    # Find the depth of each element
    elm_depth = []
    for elm in elements:
        depth1 = (elm.nodes[0].find_bottom()).depth
        depth2 = (elm.nodes[1].find_bottom()).depth
        depth3 = (elm.nodes[2].find_bottom()).depth
        elm_depth.append(mean(array([depth1, depth2, depth3])))

    # Set up figure
    if circumpolar:
        fig = figure(figsize=(128, 96))
        ax = fig.add_subplot(1,1,1, aspect='equal')
    else:
        fig = figure(figsize=(16, 8))
        ax = fig.add_subplot(1,1,1)    
    # Set colours for patches and add them to plot
    img = PatchCollection(patches, cmap=jet)
    img.set_array(array(elm_depth))
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
    title('Seafloor depth (m)', fontsize=font_sizes[0])
    cbar = colorbar(img)
    cbar.ax.tick_params(labelsize=font_sizes[2])
    img.set_clim(vmin=0, vmax=max(elm_depth))

    savefig(fig_name)


# Command-line interface
if __name__ == "__main__":

    mesh_path = raw_input("Path to mesh directory: ")
    fig_name = raw_input("Filename for figure: ")
    domain = raw_input("Global (g) or circumpolar (c)? ")
    if domain == 'c':
        circumpolar = True
    elif domain == 'g':
        circumpolar = False
    plot_depth(mesh_path, fig_name, circumpolar)
